#!/bin/bash
set -euo pipefail

IN_DIR="08_relatedness"
OUT_DIR="09_pca"
COHORT=${1}
THREADS=${3:-8}
MEMORY=${4:-16000}
NUM_PCS=${5:-10}

mkdir -p "${OUT_DIR}"

# Build merge list of pruned per-chr prefixes
MERGE_LIST="${IN_DIR}/merge_list.txt"

echo "Running PCA on unrelateds..."

plink2 \
  --pmerge-list "${IN_DIR}/merge_list.txt" \
  --keep "${IN_DIR}/unrelateds.king.cutoff.in.id" \
  --freq counts \
  --pca allele-wts vcols=chrom,ref,alt \
  --out "${OUT_DIR}/ref" \
  --threads "${THREADS}" \
  --memory "${MEMORY}"

echo "✅ PCA complete: results in ${OUT_DIR}/unrelateds.eigenvec and .eigenval"

plink2 \
  --pmerge-list "${IN_DIR}/merge_list.txt" \
  --read-freq ${OUT_DIR}/ref.acount \
  --score ${OUT_DIR}/ref.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
  --score-col-nums 6-15 \
  --out "${OUT_DIR}/${COHORT}" \
  --threads "${THREADS}" \
  --memory "${MEMORY}"

echo "✅ PCA projection complete: results in ${OUT_DIR}/${COHORT}.eigenvec and .eigenval"
