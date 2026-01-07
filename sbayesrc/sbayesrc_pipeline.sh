#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=0-06:00
#SBATCH --array=1-19
#SBATCH -o logs/out_%A_%a.log
#SBATCH -e logs/err_%A_%a.log
#SBATCH --job-name=sbayesrc

# SBayesRC pipeline for polygenic score weights
# Usage: sbatch sbayesrc_pipeline.sh
#
# Reads GWAS filenames from gwas_list.txt (one per line)
# Outputs SNP effect weights to output/main/*.txt

module load singularitypro

# Paths - adjust these for your setup
BASE="/expanse/projects/sebat1/s3/data/sebat/resources"
GWAS_DIR="${BASE}/sbayesrc_resources/cojo_sumstats/cojo"
OUT_TIDY="${BASE}/sbayesrc_output/tidy"
OUT_MAIN="${BASE}/sbayesrc_output/main"
LDM="${BASE}/sbayesrc_resources/ukbEUR_HM3"
ANNOT="${BASE}/sbayesrc_resources/annot_baseline2.2.txt"
SIF="${BASE}/sbayesrc_v0.2.6.sif"

# Get GWAS filename from list
GWAS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" gwas_list.txt)
BASENAME=${GWAS%.ma}

mkdir -p logs $OUT_TIDY $OUT_MAIN

echo "[$(date)] Starting $GWAS on $(hostname)"
echo "[$(date)] Task ID: ${SLURM_ARRAY_TASK_ID}"

# Step 1: Impute summary statistics using LD reference
# This fills in missing SNPs that are in the LD reference panel
singularity exec --bind /expanse/projects/sebat1 $SIF \
  sbayesrc \
  --ldm-eigen $LDM \
  --gwas-summary $GWAS_DIR/$GWAS \
  --impute-summary \
  --out $OUT_TIDY/${BASENAME}

# Check imputation succeeded
if [[ ! -f $OUT_TIDY/${BASENAME}.imputed.ma ]]; then
  echo "[$(date)] ERROR: Imputation failed for $GWAS"
  exit 1
fi

# Step 2: Run SBayesRC with functional annotations
# This produces SNP effect weights for PRS calculation
singularity exec --bind /expanse/projects/sebat1 $SIF \
  sbayesrc \
  --ldm-eigen $LDM \
  --gwas-summary $OUT_TIDY/${BASENAME}.imputed.ma \
  --sbayes RC \
  --annot $ANNOT \
  --out $OUT_MAIN/${BASENAME} \
  --threads ${SLURM_CPUS_PER_TASK}

echo "[$(date)] Finished $GWAS"
echo "[$(date)] Output: $OUT_MAIN/${BASENAME}.txt"
