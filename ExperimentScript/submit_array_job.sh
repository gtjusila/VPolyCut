#!/bin/bash

# Read the lines from the file into an array
mapfile -t lines < /home/htc/gtjusila/vpolycut/implementation_notes/instances_set.txt

git pull
# Number of lines in the file
num_lines=${#lines[@]}

# Check if num_lines is greater than 0
if [ $num_lines -eq 0 ]; then
    echo "No lines found in /home/htc/gtjusila/vpolycut/implementation_notes"
    exit 1
fi

echo "Submitting job arrays with $num_lines tasks"

# Submit the array job for the first Julia command (gomory)
sbatch --array=0-$(($num_lines-1)) --export=MODE=gomory /home/htc/gtjusila/vpolycut/array_job_script.sh

# Submit the array job for the second Julia command (vpc)
sbatch --array=0-$(($num_lines-1)) --export=MODE=vpc /home/htc/gtjusila/vpolycut/array_job_script.sh

