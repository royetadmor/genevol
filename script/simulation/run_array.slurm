#!/bin/bash
#SBATCH --job-name=genevol
#SBATCH --array=1-10              # adjust to number of simulations
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --partition=short         # adjust to your cluster's partition

# Usage:
#   sbatch --array=1-20 run_array.slurm <mode> [name]
#   sbatch --array=1-20 run_array.slurm one_wgd
#   sbatch --array=1-20 run_array.slurm four_wgds run1

MODE=${1:?Error: mode argument required. Usage: sbatch run_array.slurm <mode> [name]}
NAME=${2:-}

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Load modules / activate environment as needed for your cluster:
# module load python/3.10
# source /path/to/venv/bin/activate

mkdir -p "${SCRIPT_DIR}/logs"

python3 "${SCRIPT_DIR}/run_single.py" \
    --mode   "${MODE}"                \
    --sim_id "${SLURM_ARRAY_TASK_ID}" \
    ${NAME:+--name "${NAME}"}
