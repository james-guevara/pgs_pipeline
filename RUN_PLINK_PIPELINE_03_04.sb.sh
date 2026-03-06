#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=ind-shared
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH -J plink_03_04
#SBATCH -o logs/03_04_%j.out
#SBATCH -e logs/03_04_%j.err

set -euo pipefail

mkdir -p logs

MEMORY_MB=$(( SLURM_MEM_PER_NODE - 2000 ))

echo "=== Step 03: Missingness QC ==="
python PLINK_PIPELINE_03.missingness.py \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --memory "${MEMORY_MB}" \
  "$@"

echo ""
echo "=== Step 04: Summary Stats ==="
python PLINK_PIPELINE_04.summary.py \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --memory "${MEMORY_MB}"
