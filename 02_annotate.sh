#!/bin/bash
set -euo pipefail

CHR=${1:?Chromosome number required}
MEMORY=${2:-2000}
THREADS=${3:-$(nproc)}
OUT_DIR="02_annotate"
MAP_DIR="rsid_maps"
mkdir -p "${OUT_DIR}"

plink2 \
  --pfile "01_filter/chr${CHR}" \
  --update-name "${MAP_DIR}/chr${CHR}.map" 2 1 \
  --extract <(cut -f2 "${MAP_DIR}/chr${CHR}.map") \
  --write-snplist \
  --make-pgen \
  --out "${OUT_DIR}/chr${CHR}" \
  --threads "${THREADS}" \
  --memory "${MEMORY}"

