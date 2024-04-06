#!/bin/bash

# Path to the file in the repository you want to copy
FILE_PATH="cover.png"

# Start and end commits of the range
START_COMMIT=ef37ee811a4d94ae82a2e0803ee5f597ae7ddd47

# Directory where you want to copy the files
DEST_DIR="/tmp/cover/"
rm -r $DEST_DIR

# Create the destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

# Counter for file naming
COUNTER=0

# Loop through each commit in the range
for COMMIT_HASH in $(git rev-list --reverse $START_COMMIT^.. --author mattlangford@users.noreply.github.com); do
  # Check out the file at the current commit
  git checkout $COMMIT_HASH -- $FILE_PATH
  
  # Copy the file to the destination directory with a counter suffix
  printf -v FORMATTED_COUNTER "%02d" $COUNTER
  cp "$FILE_PATH" "$DEST_DIR/img_$FORMATTED_COUNTER.png"
  
  # Increment the counter
  let COUNTER=COUNTER+1
done

echo "Files copied successfully."

ffmpeg -framerate 10 -i $DEST_DIR/img_%02d.png -final_delay 300 -loop 0 -y cover.gif