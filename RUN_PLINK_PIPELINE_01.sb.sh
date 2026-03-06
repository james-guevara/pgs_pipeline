#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=ind-shared
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH -J plink_01
#SBATCH -o logs/01_%A_%a.out
#SBATCH -e logs/01_%A_%a.err

set -euo pipefail

CHR=${SLURM_ARRAY_TASK_ID:?Array index required}

# Parse arguments
BATCH=""
VCF_PATTERN=""
FORMAT="pgen"
R2=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --batch) BATCH="$2"; shift 2;;
        --vcf-pattern) VCF_PATTERN="$2"; shift 2;;
        --format) FORMAT="$2"; shift 2;;
        --r2) R2="$2"; shift 2;;
        *) echo "Unknown option: $1"; exit 1;;
    esac
done

# Default VCF pattern if not provided
if [[ -z "$VCF_PATTERN" ]]; then
    VCF_PATTERN="vcfs/chr${CHR}.dose.vcf.gz"
else
    # Replace {CHR} placeholder in pattern
    VCF_PATTERN="${VCF_PATTERN//\{CHR\}/$CHR}"
fi

# Derive plink2 memory from SLURM allocation (leave 1G headroom)
MEMORY_MB=$(( SLURM_MEM_PER_NODE - 1000 ))

mkdir -p logs

CMD="python PLINK_PIPELINE_01.py \
  --chr ${CHR} \
  --vcf ${VCF_PATTERN} \
  --threads ${SLURM_CPUS_PER_TASK} \
  --memory ${MEMORY_MB} \
  --format ${FORMAT}"

# Add optional arguments
[[ -n "$BATCH" ]] && CMD="$CMD --batch $BATCH"
[[ -n "$R2" ]] && CMD="$CMD --r2 $R2"

echo "$CMD"
eval "$CMD"
