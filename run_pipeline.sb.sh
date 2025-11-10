#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=ind-shared
#SBATCH --cpus-per-task=16
#SBATCH --time=1:00:00
#SBATCH -J pgs_pipeline_part1 
#SBATCH -o logs/out_%A_%a.log
#SBATCH -e logs/err_%A_%a.log

set -euo pipefail

CHR=${SLURM_ARRAY_TASK_ID:?Array index required}
VCF_DIR="vcfs"   # <-- edit this if needed
VCF="${VCF_DIR}/chr${CHR}_jointcall_VQSR_combined.vcf.gz"

mkdir -p logs

if [[ -f "${VCF}" ]]; then
  bash driver.sh "${CHR}" "${VCF}" 12000 16
else
  echo "[$(date)] Missing VCF for chr${CHR}: ${VCF}"
fi

