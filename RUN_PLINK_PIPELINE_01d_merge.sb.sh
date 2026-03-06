#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=ind-shared
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH -J plink_merge
#SBATCH -o logs/01d_merge_%A_%a.out
#SBATCH -e logs/01d_merge_%A_%a.err

set -euo pipefail

mkdir -p logs

MEMORY_MB=$(( SLURM_MEM_PER_NODE - 1000 ))
CHR=${SLURM_ARRAY_TASK_ID:?Array index required}

python PLINK_PIPELINE_01d.merge.py \
  --chr "${CHR}" \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --memory "${MEMORY_MB}" \
  "$@"
