#!/bin/bash
#SBATCH --job-name=vpolyhedral_experiment
#SBATCH --output={{{EXPERIMENT_PATH}}}/slurm_output/output_%A_%a.out
#SBATCH --error={{{EXPERIMENT_PATH}}}/slurm_output/error_%A_%a.err
#SBATCH --array=1-{{{N}}}
#SBATCH --time=02:15:00
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

# Load necessary modules (if any)
# module load julia
export JULIA_DEPOT_PATH="/scratch/htc/gtjusila/julia"

id=$SLURM_ARRAY_TASK_ID

# Use awk to search for the line where the id matches
line=$(awk -F'\t' -v id="$id" '$1 == id {print; exit}' {{{EXPERIMENT_PATH}}}/experiment_list.tsv)

# Check if the line was found
if [ -z "$line" ]; then
    echo "Error: No matching line found in config.tsv for id: $id"
    exit 1
fi

# Extract fields from the line (assuming tab-delimited)
IFS=$'\t' read -r id instance_path output_path <<< "$line"

# Run the Julia script with the extracted parameters
julia --project /home/htc/gtjusila/Project/VPolyCut/scripts/vpc.jl -i="$instance_path" -o="$output_path" 
