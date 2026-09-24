#!/usr/bin/env bash
# Run inside a compute allocation with Nextflow, PLINK2 and Python3 on PATH.
# Usage: prepare_pgs.sh REPOSITORY TEST_ROOT CONFIG
set -euo pipefail
repo=$(realpath "$1"); root=$(realpath "$2"); config=$(realpath "$3")
mkdir -p "$root/fixture"
cd "$root/fixture"
python3 - <<'PY'
from pathlib import Path
samples=['s%d'%i for i in range(120)]
with open('input.vcf','w') as f:
 f.write('##fileformat=VCFv4.2\n##contig=<ID=22,length=100000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
 f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join(samples)+'\n')
 for i in range(20):
  g=['0/0','0/1','1/1']*40 if i<18 else ['0/0']*120
  f.write('22\t%d\t22:%d:A:G\tA\tG\t.\tPASS\t.\tGT\t%s\n'%(i+1,i+1,'\t'.join(g)))
Path('rsid.map').write_text(''.join('22:%d:A:G\trs%d\n'%(i+1,i+1) for i in range(20)))
Path('weights.tsv').write_text('ID\tA1\tBETA\n'+''.join('rs%d\tG\t0.1\n'%(i+1) for i in range(18)))
Path('scores.tsv').write_text('trait\tweights\tid_col\tallele_col\teffect_col\ntest\t%s\t1\t2\t3\n'%(Path('weights.tsv').resolve()))
PY
plink2 --vcf input.vcf --make-pgen vzs --out base --threads 1 --memory 1500
run_case() {
 local name=$1; shift
 mkdir -p "$root/$name"; cd "$root/$name"
 nextflow run "$repo/main.nf" -c "$config" -ansi-log false -work-dir "$PWD/work" \
  --outdir "$PWD/results" --report_dir "$PWD/info" --run_summary_qc false --cohort fixture "$@"
}
base="$root/fixture/base"; map="$root/fixture/rsid.map"; scores="$root/fixture/scores.tsv"
run_case prepare --input_pfile "$base" --prepare_pgs true --score_rsid_map "$map"
run_case prepare_direct --input_pfile "$base" --prepare_pgs true --score_rsid_map "$map" --direct_inputs true --plink_publish_mode link
run_case both --input_pfile "$base" --run_scores true --score_rsid_map "$map" --score_sheet "$scores"
prepared="$root/prepare/results/04_score_input/score_input"
run_case reuse --input_pgs_pfile "$prepared" --run_scores true --score_sheet "$scores"
run_case reuse_direct --input_pgs_pfile "$prepared" --run_scores true --score_sheet "$scores" --direct_inputs true
run_case neither --input_pfile "$base"
run_case no_map --input_pfile "$base" --prepare_pgs true
python3 "$repo/tests/integration/validate_prepare_pgs.py" "$root"
