#!/bin/bash
set -euo pipefail

IN_DIR="05_summary"
OUT_DIR="06_maf"
THREADS=${1:-8}
MEMORY=${2:-16000}

mkdir -p "$OUT_DIR"

for CHR in {1..22}; do
  IN_PREFIX="${IN_DIR}/chr${CHR}"
  OUT_PREFIX="${OUT_DIR}/chr${CHR}"
  if [[ -f "${IN_PREFIX}.pgen" ]]; then
    echo "[Chr${CHR}] Applying MAF ≥ 0.05 filter..."
    plink2 \
      --pfile "${IN_PREFIX}" \
      --maf 0.05 \
      --make-pgen \
      --out "${OUT_PREFIX}" \
      --threads "$THREADS" \
      --memory "$MEMORY"
  else
    echo "⚠️  Skipping chr${CHR}: no file found."
  fi
done

echo "✅ MAF 5% filtering complete. Results in ${OUT_DIR}/"
