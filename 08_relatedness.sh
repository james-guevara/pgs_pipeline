#!/bin/bash
set -euo pipefail

IN_DIR="07_prune"
OUT_DIR="08_relatedness"
COHORT=${1}
THREADS=${2:-8}
MEMORY=${3:-16000}

mkdir -p "${OUT_DIR}"

# Build merge list of pruned per-chr prefixes
MERGE_LIST="${OUT_DIR}/merge_list.txt"
> "$MERGE_LIST"
for CHR in {1..22}; do
  PREFIX="${IN_DIR}/chr${CHR}_pruned"
  [[ -f "${PREFIX}.pgen" ]] && echo "$PREFIX" >> "$MERGE_LIST"
done

plink2 \
  --pmerge-list "$MERGE_LIST" \
  --make-king-table \
  --out "$OUT_DIR"/${COHORT} \
  --threads "$THREADS" \
  --memory "$MEMORY"

plink2 \
  --pmerge-list "$MERGE_LIST" \
  --king-cutoff 0.0884 \
  --make-just-fam \
  --out $OUT_DIR/unrelateds
