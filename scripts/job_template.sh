#!/bin/bash
#SBATCH --job-name=vpolyhedral_experiment
#SBATCH --output={{{EXPERIMENT_PATH}}}/slurm_output/output_%A_%a.out
#SBATCH --error={{{EXPERIMENT_PATH}}}/slurm_output/error_%A_%a.err
#SBATCH --array=1-{{{N}}}
#SBATCH --time=02:15:00
#SBATCH --partition=small
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G

# Load necessary modules (if any)
# module load julia
export JULIA_DEPOT_PATH="{{{JULIA_DEPOT_PATH}}}"

id=$SLURM_ARRAY_TASK_ID

# Use awk to search for the line where the id matches
line=$(awk -F'\t' -v id="$id" '$1 == id {print; exit}' {{{EXPERIMENT_PATH}}}}/experiment_list.tsv)

# Check if the line was found
if [ -z "$line" ]; then
    echo "Error: No matching line found in config.tsv for id: $id"
    exit 1
fi

# Extract fields from the line (assuming tab-delimited)
IFS=$'\t' read -r id label mode instance_path output_path <<< "$line"

# Run the Julia script with the extracted parameters
julia --project "{{{PATH_TO_SCRIPT}}}" -m="$mode" -i="$instance_path" -o="$output_path"
