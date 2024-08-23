#!/bin/bash

# Define the directories
not_sorted="not_sorted"
gomory_dir="Data/gomory"
vpc_dir="Data/vpc"

# Create the destination directories if they do not exist
mkdir -p "$gomory_dir"
mkdir -p "$vpc_dir"

# Loop through each folder in the not_sorted directory
for folder in "$not_sorted"/*; do
    # Check if it's a directory
    if [ -d "$folder" ]; then
        # Get the folder name
        folder_name=$(basename "$folder")
        
        # Move folders containing "gomory" to the gomory directory
        if [[ "$folder_name" == *gomory* ]]; then
            mv "$folder" "$gomory_dir/"
        
        # Move folders containing "vpc" to the vpc directory
        elif [[ "$folder_name" == *vpc* ]]; then
            mv "$folder" "$vpc_dir/"
        fi
    fi
done

echo "Folders moved successfully."

