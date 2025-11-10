#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=ind-shared
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH -J plink_py
#SBATCH -o logs_py/py_%A_%a.out
#SBATCH -e logs_py/py_%A_%a.err

set -euo pipefail

CHR=${SLURM_ARRAY_TASK_ID:?Array index required}
VCF="vcfs/chr${CHR}_jointcall_VQSR_combined.vcf.gz"

mkdir -p logs_py

python PLINK_PIPELINE_01.py \
  --chr "${CHR}" \
  --vcf "${VCF}" \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --memory 16000 \
  --maf 0.01
