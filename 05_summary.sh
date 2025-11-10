#!/bin/bash
set -euo pipefail

OUT_DIR="05_summary"
COHORT=${1:?Usage: $0 <cohort_name> [threads] [memory]}
THREADS=${2:-8}
MEMORY=${3:-16000}

mkdir -p "$OUT_DIR"

# Step 1 — Per-chromosome filtering (MAF ≥ 1%)
echo "Applying per-chromosome MAF filter..."
for CHR in {1..22}; do
  IN_PREFIX="04_missingness/chr${CHR}"
  OUT_PREFIX="${OUT_DIR}/chr${CHR}"
  if [[ -f "${IN_PREFIX}.pgen" ]]; then
    echo "  [Chr${CHR}] Filtering..."
    plink2 \
      --pfile "${IN_PREFIX}" \
      --maf 0.01 \
      --make-pgen \
      --out "${OUT_PREFIX}" \
      --threads "$THREADS" \
      --memory "$MEMORY"
  else
    echo "⚠️  Skipping chr${CHR}: no file found."
  fi
done

# Step 2 — Build merge list from filtered files
MERGE_LIST="${OUT_DIR}/merge_list.txt"
> "${MERGE_LIST}"

echo "Building merge list..."
MISSING_CHR=()
for CHR in {1..22}; do
  PREFIX="${OUT_DIR}/chr${CHR}"
  if [[ -f "${PREFIX}.pgen" ]]; then
    echo "${PREFIX}" >> "${MERGE_LIST}"
  else
    MISSING_CHR+=("${CHR}")
  fi
done

if ((${#MISSING_CHR[@]})); then
  echo "⚠️  Warning: missing chromosomes — ${MISSING_CHR[*]}"
fi

# Step 3 — Merge for summary QC
echo "Merging filtered data for final QC..."
plink2 \
  --pmerge-list "$MERGE_LIST" \
  --missing \
  --hardy \
  --freq \
  --out "${OUT_DIR}/${COHORT}" \
  --threads "$THREADS" \
  --memory "$MEMORY"

# Step 4 — Quick summary counts
{
  echo -e "Metric\tCount"
  echo -e "Samples\t$(($(wc -l < ${OUT_DIR}/${COHORT}.smiss) - 1))"
  echo -e "Variants\t$(($(wc -l < ${OUT_DIR}/${COHORT}.vmiss) - 1))"
} > "${OUT_DIR}/summary_counts.txt"

# Step 5 — Clean up merge list
rm -f "${MERGE_LIST}"

echo "✅ Summary and per-chromosome filtered data written to ${OUT_DIR}/"
