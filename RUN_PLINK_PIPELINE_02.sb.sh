#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=ind-shared
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH -J plink_concat
#SBATCH -o logs/02_concat_%j.out
#SBATCH -e logs/02_concat_%j.err

set -euo pipefail

mkdir -p logs

MEMORY_MB=$(( SLURM_MEM_PER_NODE - 2000 ))

python PLINK_PIPELINE_02.concat.py \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --memory "${MEMORY_MB}" \
  "$@"
