#!/bin/bash
#SBATCH --job-name=vpolycut_experiment
#SBATCH --time=02:15:00
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G

# Load necessary modules (if any)
# module load julia
export JULIA_DEPOT_PATH="/scratch/htc/gtjusila/julia"

# Check if the necessary environment variables are set
if [ -z "$MODE" ] || [ -z "$EXPERIMENT_FOLDER" ]; then
    echo "Error: MODE or EXPERIMENT_FOLDER environment variables are not set."
    exit 1
fi

# Read the lines from the instances file, which is passed as the second argument
instances_file=$1

# Check if the instances file exists
if [ ! -f "$instances_file" ]; then
    echo "Error: Instances file '$instances_file' not found."
    exit 1
fi

# Read the lines from the instances file into an array
mapfile -t lines < "$instances_file"

# Get the line corresponding to this array task ID
line="${lines[$SLURM_ARRAY_TASK_ID]}"

# Run the appropriate Julia command based on the MODE environment variable
if [ "$MODE" == "gomory" ]; then
    julia --project run_experiment.jl -m=gomory -i="$line" -f="$EXPERIMENT_FOLDER"
elif [ "$MODE" == "intersection" ]; then
    julia --project run_experiment.jl -m=intersection -i="$line" -f="$EXPERIMENT_FOLDER"
elif [ "$MODE" == "vpc" ]; then
    julia --project run_experiment.jl -m=vpc -i="$line" -f="$EXPERIMENT_FOLDER"
else
    echo "Unknown mode: $MODE"
    exit 1
fi

