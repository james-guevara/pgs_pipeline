#!/bin/bash
set -euo pipefail

CHR=${1:?Chromosome number required}
VCF=${2:?VCF path required}
MEMORY=${3:-2000}
THREADS=${4:-$(nproc)}
OUT_DIR="01_filter"
mkdir -p "${OUT_DIR}"

plink2 \
  --vcf "${VCF}" \
  --snps-only just-acgt \
  --max-alleles 2 \
  --var-filter \
  --mac 10 \
  --geno 0.05 \
  --set-all-var-ids @:#:\$r:\$a \
  --make-pgen \
  --out "${OUT_DIR}/chr${CHR}" \
  --threads "${THREADS}" \
  --memory "${MEMORY}"
