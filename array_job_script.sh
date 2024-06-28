#!/bin/bash
#SBATCH --job-name=vpolycut_experiment
#SBATCH --output=/home/htc/gtjusila/vpolycut/slurm_output/%x_%A_%a.out
#SBATCH --error=/home/htc/gtjusila/vpolycut/slurm_output/%x_%A_%a.err
#SBATCH --time=01:15:00
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --constraint=E7-4830v3

# Load necessary modules (if any)
# module load julia
export JULIA_DEPOT_PATH="/scratch/htc/gtjusila/julia"

mapfile -t lines < /home/htc/gtjusila/vpolycut/implementation_notes/instances_set.txt

# Get the line corresponding to this array task ID
line="${lines[$SLURM_ARRAY_TASK_ID]}"

# Change to the working directory
cd /home/htc/gtjusila/vpolycut/vpolycut/

# Run the appropriate julia command based on the MODE environment variable
if [ "$MODE" == "gomory" ]; then
    julia --project=/home/htc/gtjusila/vpcenv /home/htc/gtjusila/vpolycut/vpolycut/experiment/run_experiment.jl -m gomory -i "$line"
elif [ "$MODE" == "vpc" ]; then
    julia --project=/home/htc/gtjusila/vpcenv /home/htc/gtjusila/vpolycut/vpolycut/experiment/run_experiment.jl -m vpc -i "$line"
else
    echo "Unknown mode: $MODE"
    exit 1
fi
