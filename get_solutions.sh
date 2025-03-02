#!/bin/bash

# Define the URL and target directory
URL="https://miplib.zib.de/downloads/solutions.zip"
TARGET_DIR="./experiment_data/miplibbenchsolutions"
TEMP_DIR="./temp_solutions"

# Create necessary directories
mkdir -p "$TARGET_DIR"
mkdir -p "$TEMP_DIR"

# Download the zip file
ZIP_FILE="$TEMP_DIR/solutions.zip"
echo "Downloading solutions.zip..."
wget -O "$ZIP_FILE" "$URL"

# Extract the downloaded zip file
echo "Extracting solutions.zip..."
unzip -q "$ZIP_FILE" -d "$TEMP_DIR"

# Process each subfolder inside the extracted directory
echo "Processing extracted folders..."
found_solutions=0
for parent_dir in "$TEMP_DIR"/*/; do
    for dir in "$parent_dir"*/; do
        if [ -d "$dir/1" ]; then
            sol_gz_file=$(find "$dir/1" -name "*.sol.gz" | head -n 1)
            if [ -n "$sol_gz_file" ]; then
                sol_file="$(basename "$sol_gz_file" .gz)"
                gunzip -c "$sol_gz_file" > "$TARGET_DIR/$sol_file"
                echo "Transferred: $sol_file"
                ((found_solutions++))
            else
                echo "No .sol.gz file found in $dir/1"
            fi
        else
            echo "No '1' folder found in $dir"
        fi
    done
done

# Clean up temporary files
echo "Cleaning up temporary files..."
rm -rf "$TEMP_DIR"

echo "All .sol files have been extracted and moved to $TARGET_DIR. Total files transferred: $found_solutions."