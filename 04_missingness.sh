#!/bin/bash
set -euo pipefail

MEMORY=${1:-4000}
THREADS=${2:-$(nproc)}
OUT_DIR="04_missingness"
mkdir -p "${OUT_DIR}"

MERGE_LIST="${OUT_DIR}/merge_list.txt"
> "${MERGE_LIST}"

echo "Building merge list..."
MISSING_CHR=()
for CHR in {1..22}; do
  PREFIX="03_maf/chr${CHR}"
  if [[ -f "${PREFIX}.pgen" ]]; then
    echo "${PREFIX}" >> "${MERGE_LIST}"
  else
    MISSING_CHR+=("${CHR}")
  fi
done

if ((${#MISSING_CHR[@]})); then
  echo "⚠️  Skipping missing chromosomes: ${MISSING_CHR[*]}"
fi
echo "Merge list written to ${MERGE_LIST}"
echo

echo "Running PLINK2 missingness..."
plink2 \
  --pmerge-list "${MERGE_LIST}" \
  --missing \
  --out "${OUT_DIR}/missingness" \
  --threads "${THREADS}" \
  --memory "${MEMORY}"

echo "Generating fail lists..."
awk 'NR>1 && $4>0.05 {print $1}' "${OUT_DIR}/missingness.smiss" > "${OUT_DIR}/fail_samples.txt"
awk 'NR>1 && $5>0.05 {print $2}' "${OUT_DIR}/missingness.vmiss" > "${OUT_DIR}/fail_variants.txt"

echo "Fail lists:"
echo "  $(wc -l < "${OUT_DIR}/fail_samples.txt") samples"
echo "  $(wc -l < "${OUT_DIR}/fail_variants.txt") variants"
echo

echo "Applying missingness filters..."
for CHR in {1..22}; do
  PREFIX="03_maf/chr${CHR}"
  if [[ -f "${PREFIX}.pgen" ]]; then
    echo "[Chr${CHR}] filtering..."
    plink2 \
      --pfile "${PREFIX}" \
      --remove "${OUT_DIR}/fail_samples.txt" \
      --exclude "${OUT_DIR}/fail_variants.txt" \
      --make-pgen \
      --out "${OUT_DIR}/chr${CHR}" \
      --threads "${THREADS}" \
      --memory "${MEMORY}"
  fi
done

echo "✅ Missingness QC + filtering complete. Filtered files in ${OUT_DIR}/"
