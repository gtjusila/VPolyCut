#!/bin/bash
#SBATCH --job-name=vpolycut_experiment
#SBATCH --time=01:15:00
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G

# Load necessary modules (if any)
# module load julia
export JULIA_DEPOT_PATH="/scratch/htc/gtjusila/julia"

mapfile -t lines < ../implementation_notes/instances_set.txt

# Get the line corresponding to this array task ID
line="${lines[$SLURM_ARRAY_TASK_ID]}"

# Run the appropriate julia command based on the MODE environment variable
if [ "$MODE" == "gomory" ]; then
    julia --project run_experiment.jl -m gomory -i "$line"
elif [ "$MODE" == "vpc" ]; then
    julia --project run_experiment.jl -m vpc -i "$line"
else
    echo "Unknown mode: $MODE"
    exit 1
fi
