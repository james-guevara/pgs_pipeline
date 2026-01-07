#!/bin/bash
#SBATCH --account=ddp195
#SBATCH --partition=ind-shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --array=1-19
#SBATCH -o logs/%A_%a.out
#SBATCH -e logs/%A_%a.err
#SBATCH --job-name=sbayesrc

#===============================================================================
# SBayesRC Pipeline
#
# Computes polygenic score weights from GWAS summary statistics.
#
# Usage:
#   1. Create gwas_list.txt with one .ma filename per line
#   2. Update --array above to match: wc -l gwas_list.txt
#   3. sbatch run_sbayesrc.sh
#
# Input:  input/*.ma (COJO format)
# Output: output/weights/*.txt (SNP effect weights for PRS)
#===============================================================================

set -euo pipefail

#-------------------------------------------------------------------------------
# CONFIG - edit these paths
#-------------------------------------------------------------------------------
WORK_DIR="/expanse/projects/sebat1/s3/data/sebat/pgs_pipeline"
INPUT_DIR="${WORK_DIR}/input"
OUTPUT_DIR="${WORK_DIR}/output"

# Resources (shared, don't need to change)
RESOURCES="/expanse/projects/sebat1/s3/data/sebat/resources"
LDM="${RESOURCES}/sbayesrc_resources/ukbEUR_HM3"
ANNOT="${RESOURCES}/sbayesrc_resources/annot_baseline2.2.txt"
SIF="${RESOURCES}/sbayesrc_v0.2.6.sif"

#-------------------------------------------------------------------------------
# SETUP
#-------------------------------------------------------------------------------
module load singularitypro

GWAS_LIST="${WORK_DIR}/gwas_list.txt"
GWAS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$GWAS_LIST")
NAME=${GWAS%.ma}

TIDY_DIR="${OUTPUT_DIR}/tidy"
WEIGHTS_DIR="${OUTPUT_DIR}/weights"
mkdir -p "$TIDY_DIR" "$WEIGHTS_DIR" "${WORK_DIR}/logs"

#-------------------------------------------------------------------------------
# CHECK IF ALREADY DONE
#-------------------------------------------------------------------------------
if [[ -f "${WEIGHTS_DIR}/${NAME}.txt" ]]; then
    echo "[$(date)] SKIP: ${NAME} already complete"
    exit 0
fi

echo "========================================"
echo "SBayesRC: ${NAME}"
echo "========================================"
echo "Start:    $(date)"
echo "Host:     $(hostname)"
echo "Task ID:  ${SLURM_ARRAY_TASK_ID}"
echo "Input:    ${INPUT_DIR}/${GWAS}"
echo "========================================"

#-------------------------------------------------------------------------------
# STEP 1: Impute missing SNPs
#-------------------------------------------------------------------------------
echo ""
echo "[$(date)] Step 1: Imputing summary statistics..."

if [[ -f "${TIDY_DIR}/${NAME}.imputed.ma" ]]; then
    echo "[$(date)] Using existing imputed file"
    GWAS_FOR_SBAYESRC="${TIDY_DIR}/${NAME}.imputed.ma"
else
    IMPUTE_LOG="${TIDY_DIR}/${NAME}.impute.log"
    singularity exec --bind /expanse/projects/sebat1 "$SIF" \
        gctb --ldm-eigen "$LDM" \
             --gwas-summary "${INPUT_DIR}/${GWAS}" \
             --impute-summary \
             --out "${TIDY_DIR}/${NAME}" 2>&1 | tee "$IMPUTE_LOG"

    # Check if imputation created output, or if no imputation was needed
    if [[ -f "${TIDY_DIR}/${NAME}.imputed.ma" ]]; then
        GWAS_FOR_SBAYESRC="${TIDY_DIR}/${NAME}.imputed.ma"
    elif grep -q "No SNP needs to be imputed" "$IMPUTE_LOG"; then
        echo "[$(date)] No imputation needed - using original file"
        GWAS_FOR_SBAYESRC="${INPUT_DIR}/${GWAS}"
    else
        echo "[$(date)] ERROR: Imputation failed"
        exit 1
    fi
fi

echo "[$(date)] Step 1 complete: $(wc -l < "$GWAS_FOR_SBAYESRC") SNPs"

#-------------------------------------------------------------------------------
# STEP 2: SBayesRC
#-------------------------------------------------------------------------------
echo ""
echo "[$(date)] Step 2: Running SBayesRC..."

singularity exec --bind /expanse/projects/sebat1 "$SIF" \
    gctb --ldm-eigen "$LDM" \
         --gwas-summary "$GWAS_FOR_SBAYESRC" \
         --sbayes RC \
         --annot "$ANNOT" \
         --out "${WEIGHTS_DIR}/${NAME}" \
         --thread "${SLURM_CPUS_PER_TASK}"

if [[ ! -f "${WEIGHTS_DIR}/${NAME}.txt" ]]; then
    echo "[$(date)] ERROR: SBayesRC failed"
    exit 1
fi

#-------------------------------------------------------------------------------
# DONE
#-------------------------------------------------------------------------------
echo ""
echo "========================================"
echo "COMPLETE: ${NAME}"
echo "End:      $(date)"
echo "Output:   ${WEIGHTS_DIR}/${NAME}.txt"
echo "========================================"
