#!/bin/bash
set -euo pipefail

CHR=${1:?Chromosome number required}
MEMORY=${2:-2000}
THREADS=${3:-$(nproc)}
OUT_DIR="03_maf"
mkdir -p "${OUT_DIR}"

plink2 \
  --pfile "02_annotate/chr${CHR}" \
  --maf 0.01 \
  --make-pgen \
  --out "${OUT_DIR}/chr${CHR}" \
  --threads "${THREADS}" \
  --memory "${MEMORY}"

