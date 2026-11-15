#!/bin/bash
#SBATCH --job-name=genevol_bootstrap
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH -A itaym-users_v2
#SBATCH -p itaym-pool
#SBATCH --qos=owner

OUT_DIR="${1}"
PARAM_FILE="${2}"
BOOTSTRAP_SCRIPT="${3}"

python3 "${BOOTSTRAP_SCRIPT}" \
    --run-rep \
    --rep-id "$SLURM_ARRAY_TASK_ID" \
    --out-dir "$OUT_DIR" \
    --param-file "$PARAM_FILE"
