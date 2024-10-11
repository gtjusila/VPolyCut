#!/bin/bash

# Default values for optional arguments
constraint=""
instances_file=""
result_path=""
modes=()

# Function to print the usage of the script
usage() {
    echo "Usage: $0 --constraint <constraint> --instances_file <instances_file> --result_path <result_path> --modes <mode1,mode2,...>"
    exit 1
}

# Parse named arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --constraint)
            constraint="$2"
            shift 2
            ;;
        --instances_file)
            instances_file="$2"
            shift 2
            ;;
        --result_path)
            result_path="$2"
            shift 2
            ;;
        --modes)
            IFS=',' read -r -a modes <<< "$2"
            # Trim whitespace from each mode
            for i in "${!modes[@]}"; do
                modes[$i]=$(echo "${modes[$i]}" | xargs)
            done
            shift 2
            ;;
        *)
            echo "Unknown parameter passed: $1"
            usage
            ;;
    esac
done

# Check if all required arguments are provided
if [[ -z "$constraint" || -z "$instances_file" || -z "$result_path" || ${#modes[@]} -eq 0 ]]; then
    echo "Error: Missing required arguments."
    usage
fi

# Check if the instances file exists
if [ ! -f "$instances_file" ]; then
    echo "Error: Instances file '$instances_file' not found."
    exit 1
fi

# Read the lines from the file into an array
mapfile -t lines < "$instances_file"

git pull
# Number of lines in the file
num_lines=${#lines[@]}

# Check if num_lines is greater than 0
if [ $num_lines -eq 0 ]; then
    echo "No lines found in $instances_file"
    exit 1
fi

# Create the results directory if it doesn't exist
mkdir -p "$result_path"

echo "Submitting job arrays with $num_lines tasks"

for mode in "${modes[@]}"; do
    for i in $(seq 0 $(($num_lines-1))); do
        # Generate a timestamp
        timestamp=$(date +"%Y%m%dT%H%M%S")

        # Get the instance (line) corresponding to the current task
        instance="${lines[$i]}"

        # Create an experiment label
        experiment_label="${timestamp}_${instance}_${mode}"
        
        # Define the experiment result path
        experiment_result_path="${result_path}/${experiment_label}"

        # Check if the experiment directory already exists
        if [ -d "$experiment_result_path" ]; then
            echo "Error: Experiment directory '$experiment_result_path' already exists."
            exit 1
        else
            # Create the experiment directory
            mkdir -p "$experiment_result_path"
        fi

        # Submit the SLURM job with the experiment folder and instance file as arguments
        sbatch --array=$i --export=MODE=$mode,EXPERIMENT_FOLDER=$experiment_result_path --constraint=$constraint --output=$experiment_result_path/%x_%j.out --error=$experiment_result_path/%x_%j.err array_job_script.sh "$instances_file"
    done
done

