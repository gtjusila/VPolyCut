#!/bin/bash

# Define variables
URL="https://miplib.zib.de/downloads/benchmark.zip"
ZIP_FILE="benchmark.zip"
DOWNLOAD_DIR="./benchmark"
DEST_DIR="./experiment_data/miplibbench"
LIST_FILE="$DEST_DIR/instances_list.txt"

# Clean the destination directory
echo "Cleaning $DEST_DIR..."
rm -rf "$DEST_DIR"
mkdir -p "$DEST_DIR"

# Remove previous instance list file if exists
rm -f "$LIST_FILE"

# Create a temporary download directory
mkdir -p "$DOWNLOAD_DIR"

# Download the zip file with a progress bar and wait for completion
echo "Downloading benchmark.zip..."
wget --progress=bar:force -O "$ZIP_FILE" "$URL" && echo "Download complete." || { echo "Download failed! Exiting."; exit 1; }

# Verify the zip file exists
if [ ! -f "$ZIP_FILE" ]; then
    echo "Error: benchmark.zip not found!"
    exit 1
fi

# Extract the zip file
echo "Extracting benchmark.zip..."
unzip -o "$ZIP_FILE" -d "$DOWNLOAD_DIR" || { echo "Extraction failed! Exiting."; exit 1; }

# Extract .mps.gz files, move them, and create the list
echo "Processing .mps.gz files..."
for file in "$DOWNLOAD_DIR"/*.mps.gz; do
    if [ -f "$file" ]; then
        output_file="$DEST_DIR/$(basename "${file%.gz}")"
        gunzip -c "$file" > "$output_file"
        echo "Extracted: $(basename "$output_file")"

        # Add file name to instances_list.txt without .mps extension
        echo "$(basename "$output_file" .mps)" >> "$LIST_FILE"
    fi
done

# Cleanup temporary files
echo "Cleaning up..."
rm -f "$ZIP_FILE"
rm -rf "$DOWNLOAD_DIR"

echo "Extraction and copying complete."
echo "Instance list saved in $LIST_FILE."

