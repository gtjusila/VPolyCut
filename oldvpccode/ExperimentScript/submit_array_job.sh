#!/bin/bash

# Check if a command-line argument is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <constraint>"
    exit 1
fi

# Read the constraint from the command-line argument
constraint=$1

# Read the lines from the file into an array
mapfile -t lines < ../implementation_notes/instances_set.txt

git pull
# Number of lines in the file
num_lines=${#lines[@]}

# Check if num_lines is greater than 0
if [ $num_lines -eq 0 ]; then
    echo "No lines found in ../implementation_notes/instances_set.txt"
    exit 1
fi

echo "Submitting job arrays with $num_lines tasks"

# Submit the array job for the first Julia command (gomory) with constraint
#sbatch --array=0-$(($num_lines-1)) --export=MODE=gomory --constraint=$constraint --output=../temp/%x_%j.out --error=../temp/%x_%j.err array_job_script.sh
# Submit the array job for the second Julia command (vpc) with constraint
sbatch --array=0-$(($num_lines-1)) --export=MODE=vpc --constraint=$constraint --output=../temp/%x_%j.out --error=../temp/%x_%j.err array_job_script.sh
