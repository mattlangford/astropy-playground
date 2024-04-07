# This is all copied over from the ipython notebook

import numpy as np
import scipy.optimize
import astropy.units as u
import astropy.time
import astropy.constants
import astropy.coordinates
import pytz
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib

# Makes queries against the National Weather Service API to get cloud cover
from datetime import datetime

import time
import googlemaps
import os
import urllib.request
import urllib.error
import json 
import pandas as pd
import bisect
import pytz

astropy.coordinates.solar_system_ephemeris.set("de430")

# Compute angular diistance between sun and moon
def distance_contact(
    location: astropy.coordinates.EarthLocation,
    time: astropy.time.Time
) -> u.Quantity:
    coordinate_sun = astropy.coordinates.get_sun(time)
    coordinate_moon = astropy.coordinates.get_body("moon", time)

    frame_local = astropy.coordinates.AltAz(obstime=time, location=location)

    alt_az_sun = coordinate_sun.transform_to(frame_local)
    alt_az_moon = coordinate_moon.transform_to(frame_local)
    return alt_az_moon.separation(alt_az_sun).degree

# Compute overlap percentage between sun and moon
def overlap_percent(
    location: astropy.coordinates.EarthLocation,
    time: astropy.time.Time
) -> float:    
    radius_sun = astropy.constants.R_sun
    radius_moon = 1737.4 * u.km
    
    coordinate_sun = astropy.coordinates.get_sun(time)
    coordinate_moon = astropy.coordinates.get_body("moon", time)
    
    frame_local = astropy.coordinates.AltAz(obstime=time, location=location)
    
    alt_az_sun = coordinate_sun.transform_to(frame_local)
    alt_az_moon = coordinate_moon.transform_to(frame_local)

    sun_angle = np.arctan(radius_sun / alt_az_sun.distance)
    moon_angle = np.arctan(radius_moon / alt_az_moon.distance)
    separation = alt_az_moon.separation(alt_az_sun)
    
    sun_min = -sun_angle.value
    sun_max = sun_angle.value
    moon_min = (separation - moon_angle).rad
    moon_max = (separation + moon_angle).rad
    
    # The apparent size of the moon is slightly large than the sun
    return np.maximum(0.0, np.minimum(sun_max, moon_max) - np.maximum(sun_min, moon_min)) / (2 * sun_angle.value)

def _query_with_retry(query, retries=5):
    count = 0
    sleep = 1
    while count < retries:
        count += 1
        try:
            with urllib.request.urlopen(query) as url:
                return json.load(url)
        except urllib.error.HTTPError as e:
            if e.code == 404:
                print(f"Unable to make query {query}")
                raise e
            if e.code == 500:
                print(f"{e}, retrying in {sleep}s")
                time.sleep(sleep)
                sleep *= 1.5
                continue
    raise Exception("Timeout, unable to make query")
                
def query_sky_cover(lat, lon, query_time):
    print (f"Querying at {lat}, {lon} t={query_time}")

    data = _query_with_retry(f"https://api.weather.gov/points/{lat},{lon}")
    tz = pytz.timezone(data["properties"]["timeZone"])
    forecast_url = data["properties"]["forecastGridData"]
    data = _query_with_retry(forecast_url)

    query_time = query_time.to_datetime(timezone=tz)
    dates = [datetime.fromisoformat(v["validTime"].split("/")[0]) for v in data["properties"]["skyCover"]["values"]]
    i = bisect.bisect_left(dates, query_time)
    if i > 0:
        return data["properties"]["skyCover"]["values"][i - 1]["value"]
    
    # Unable to find anything
    raise "Unable to find forecast"

def query_sky_cover_safe(lat, lon, t):
    try:
        return query_sky_cover(lat, lon, t)
    except:
        print(f"No data for {lat} {lon}")
        return 100


gmaps = googlemaps.Client(key=os.environ['GOOGLE_API_KEY'])

def get_travel_times(orig_lat, orig_lon, destinations, query_time):
    results = []
    s = ""
    for i, dest in enumerate(destinations):
        s += f"{dest[0]},{dest[1]}|"

        # Can only make a handful of queries per call
        if i + 2 < len(destinations) and (i == 0 or i % 20 != 0):
            continue
        s = s[:-1]
        directions_result = gmaps.distance_matrix(f"{orig_lat},{orig_lon}", s, mode="driving", departure_time=query_time)
        for element in directions_result['rows'][0]['elements']:
            if 'duration_in_traffic' in element:
                results.append((float(element['duration_in_traffic']['value']) * u.s).to(u.hr).value)
            else:
                results.append(0.0)
        s = ""
    return results

# Generate points within totality along the path of the eclipse
center_time = astropy.time.Time("2024-04-08 19:14:30")

lat_min, lon_min = 36.0, -90.7
lat_max, lon_max = 45.6, -70.0
lat_center = lat_min + (lat_max - lat_min) / 2.0
lon_center = lon_min + (lon_max - lon_min) / 2.0

points = 50
lats = np.linspace(lat_min, lat_max, points)
lons = np.linspace(lon_min, lon_max, points)
lat, lon = np.meshgrid(lats, lons)

dt = 16 * u.min
count = 20
centers = []
center_ts = []
ts = np.linspace(center_time - dt, center_time + dt, count)
for t in ts:
    locations = astropy.coordinates.EarthLocation(lat=lat * u.deg, lon=lon * u.deg, height=0 * u.m)
    overlap = overlap_percent(locations, t)
    indicies = overlap > 0.999
    for x, y in zip(lat[indicies], lon[indicies]):
        centers.append([x, y])
        center_ts.append(t)
    
centers = np.unique(np.array(centers), axis=0)
print(f"{len(centers)} totality points loaded")

# Query weather across the map
points = 12
cloud_lats = np.linspace(lat_min, lat_max, points)
cloud_lons = np.linspace(lon_min, lon_max, points)
cloud_lat, cloud_lon = np.meshgrid(cloud_lats, cloud_lons)

# Assume times are nearby
avg_t = ts[len(ts) // 2]

cover = np.zeros_like(cloud_lat)

t0 = time.time()
for i in range(cloud_lat.shape[0]):
    print(f"Querying i={i} Last iteration was { time.time() - t0}s")
    t0 = time.time()
    for j in range(cloud_lat.shape[1]):
        cover[i, j] = query_sky_cover_safe(cloud_lat[i, j], cloud_lon[i, j], avg_t)

tz = pytz.timezone('America/New_York')
weather_time = datetime.now(tz) 

# Query weather at each point within totality
center_covers = []
t0 = time.time()
for i, c in enumerate(centers):
    center_covers.append(query_sky_cover_safe(c[0], c[1], center_ts[i]))

    sample_rate = 20
    if i % sample_rate == 0:
        print(f"Making {sample_rate / (time.time() - t0)} queries/s")
        t0 = time.time()

center_covers = np.array(center_covers)

# Requery anything that timed out
missing = np.sum(cover == 100)
coords = np.where(cover == 100)
for (i, j) in zip(coords[0], coords[1]):
    cover[i, j] = query_sky_cover_safe(cloud_lat[i, j], cloud_lon[i, j], avg_t)

print(f"Found {np.sum(cover == 100) - missing} previous missing entries (out of {missing})")

for i, c in enumerate(centers):
    if center_covers[i] == 100:
        center_covers[i] = query_sky_cover_safe(c[0], c[1], center_ts[i])

# Query weather at each point within totality
leave_time = astropy.time.Time("2024-04-08 10:00:00") # 6am the morning of
center_traffic = []
pit_lat, pit_lon = 40.4406, -79.9959
center_traffic = get_travel_times(pit_lat, pit_lon, centers, leave_time.to_datetime())
center_traffic = np.array(center_traffic)
print(center_traffic)

# Generate best targets for a few different hour drive
best_targets = {}
for time_hour in [3.0, 5.0, 8.0, 10.0, 12.0, 15.0]:
    time_mask = np.logical_and(center_traffic < time_hour, center_traffic > 0.0)
    arr = np.ma.array(center_covers, mask=~time_mask)
    if arr.any():
        best_index = arr.argmin()
        best_targets[time_hour] = {'index': best_index, 'center': centers[best_index], 'cover': center_covers[best_index]}

print(best_targets)

# Plot the map, eclipse points, and weather

from mpl_toolkits.basemap import Basemap
import matplotlib

plt.figure(figsize=(15, 15))

# Lambert Conformal map of lower 48 states.
m = Basemap(llcrnrlon=lon_min, llcrnrlat=lat_min, urcrnrlon=lon_max, urcrnrlat=lat_max,
            lat_0=lat_center, lon_0=lon_center, projection='merc')

m.readshapefile('assets/cb_2018_us_state_500k/cb_2018_us_state_500k', 'states', drawbounds=True)
m.readshapefile('assets/tl_2017_us_primaryroads/tl_2017_us_primaryroads', 'roads', drawbounds=True, linewidth=0.1)

y, x = m(cloud_lon, cloud_lat)
plt.contourf(y, x, cover, 15, cmap='Blues_r', vmin=0, vmax=100)

y, x = m(centers[:, 1], centers[:, 0])
mask = np.logical_and(center_covers < 100, center_traffic > 0)
plt.scatter(y[mask], x[mask], c=center_covers[mask], s=70, cmap='Blues_r', vmin=0, vmax=100, edgecolors='k', linewidth=0.4, zorder=2)

cbar = plt.colorbar(shrink=0.4)
cbar.set_label('Sky Cover', rotation=270)

y, x = m(pit_lon, pit_lat)
plt.plot(y, x, 'x', c='r')

readme = ""

readme += "The last sky cover plot generated:\n"
readme += "![cover](cover.png)\n\n"
readme += "The best targets by hour:\n"

root_y, root_x = m(-81.0, 43.5)
for duration in sorted(best_targets.keys(), key=lambda i: (best_targets[i]['center'][1], i)):
    target = best_targets[duration]
    readme += f" - Less than {duration:.0f} hours: ({target['center'][0]:.5f}, {target['center'][1]:.5f}) cover: {target['cover']}\n"
    center = target['center']
    y, x = m(center[1], center[0])
    text = plt.text(root_y, root_x, f'<{duration:.0f}hr\nsk={target["cover"]}', horizontalalignment='center', verticalalignment='center', 
                    fontsize=7, color='black', bbox=dict(facecolor='white', edgecolor='black', boxstyle='round'))
    plt.plot([y, root_y], [x, root_x], c='r', linewidth=0.5)
    root_y += 100000
    root_x += 50000

print(readme)
with open('README.md', 'w') as file:
    file.write(readme)

plt.title(f"Sky Cover along Totality (queried_at={weather_time})")
plt.savefig('cover.png', bbox_inches="tight")