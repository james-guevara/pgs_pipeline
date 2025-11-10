#!/bin/bash
set -euo pipefail

IN_DIR="06_maf"
OUT_DIR="07_prune"
THREADS=${1:-8}
MEMORY=${2:-16000}

mkdir -p "$OUT_DIR"

# Step 1: LD pruning per chromosome
for CHR in {1..22}; do
  IN_PREFIX="${IN_DIR}/chr${CHR}"
  OUT_PREFIX="${OUT_DIR}/chr${CHR}"
  if [[ -f "${IN_PREFIX}.pgen" ]]; then
    echo "[Chr${CHR}] Performing LD pruning..."
    plink2 \
      --pfile "${IN_PREFIX}" \
      --indep-pairwise 200 50 0.1 \
      --out "${OUT_PREFIX}" \
      --threads "$THREADS" \
      --memory "$MEMORY"
  else
    echo "⚠️  Skipping chr${CHR}: no file found."
  fi
done

echo "✅ LD pruning complete."

# Step 2: Extract pruned SNPs to new per-chromosome files
for CHR in {1..22}; do
  IN_PREFIX="${IN_DIR}/chr${CHR}"
  OUT_PREFIX="${OUT_DIR}/chr${CHR}_pruned"
  PRUNE_FILE="${OUT_DIR}/chr${CHR}.prune.in"
  if [[ -f "${PRUNE_FILE}" ]]; then
    echo "[Chr${CHR}] Extracting pruned SNPs..."
    plink2 \
      --pfile "${IN_PREFIX}" \
      --extract "${PRUNE_FILE}" \
      --make-pgen \
      --out "${OUT_PREFIX}" \
      --threads "$THREADS" \
      --memory "$MEMORY"
  fi
done

echo "✅ Extraction of pruned SNPs complete."
