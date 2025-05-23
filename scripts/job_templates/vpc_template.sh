#!/bin/bash
#SBATCH --job-name=vpolyhedral_experiment
#SBATCH --array=1-{{{N}}}
#SBATCH --constraint=Gold6338
#SBATCH --time=02:15:00
#SBATCH --partition=big
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

set -euo pipefail

# -------- 1. environment --------
export JULIA_DEPOT_PATH="/scratch/htc/gtjusila/julia"
export JULIA_CPU_TARGET="generic;icelake-server"
id=$SLURM_ARRAY_TASK_ID

cd /home/htc/gtjusila/Project/VPolyCut/
# -------- 2. look up the row in the TSV --------
line=$(awk -F'\t' -v id="$id" '$1 == id {print; exit}' \
        {{{EXPERIMENT_PATH}}}/experiment_list.tsv)

if [ -z "$line" ]; then
    echo "Error: No matching line found in experiment_list.tsv for id: $id"
    exit 1
fi

IFS=$'\t' read -r _ instance_path output_path solution_path <<< "$line"

# -------- 3. redirect stdout & stderr to output_path --------
mkdir -p "$output_path"
echo "Parsed output_path: '$output_path'"

log_base="${output_path}/slurm_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
exec >"${log_base}.out" 2>"${log_base}.err"

echo "Logs for this task are now in ${log_base}.{out,err}"

# -------- 4. run the Julia experiment --------
julia --project /home/htc/gtjusila/Project/VPolyCut/scripts/vpc.jl \
      -i="$instance_path" \
      -o="$output_path" \
      -s="$solution_path" \
      -c="{{{CONFIG_PATH}}}"