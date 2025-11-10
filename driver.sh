#!/bin/bash
set -euo pipefail

CHR=${1:?Chromosome number required}
VCF=${2:?Path to VCF required}
MEMORY=${3:-16000}
THREADS=${4:-16}

echo "[$(date)] Starting chromosome ${CHR} on $(hostname)"

bash 01_filter.sh   "${CHR}" "${VCF}" "${MEMORY}" "${THREADS}"
bash 02_annotate.sh "${CHR}" "${MEMORY}" "${THREADS}"
bash 03_maf.sh      "${CHR}" "${MEMORY}" "${THREADS}"

echo "[$(date)] Done chromosome ${CHR}"
