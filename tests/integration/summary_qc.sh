#!/usr/bin/env bash
# Run in a compute allocation. Also exercises the existing preparation suite.
set -euo pipefail
repo=$(realpath "$1"); root=$(realpath "$2"); config=$(realpath "$3")
bash "$repo/tests/integration/prepare_pgs.sh" "$repo" "$root" "$config"
base="$root/fixture/base"; map="$root/fixture/rsid.map"; scores="$root/fixture/scores.tsv"
prepared="$root/prepare/results/04_score_input/score_input"
run_case() {
 local name=$1; shift
 mkdir -p "$root/$name/info"; cd "$root/$name"
 nextflow run "${PGS_ENTRYPOINT:-$repo/main.nf}" -c "$config" -ansi-log false -work-dir "$PWD/work" \
  -with-trace "$PWD/info/trace.txt" --outdir "$PWD/results" --report_dir "$PWD/info" \
  --run_pca false --run_summary_qc true --cohort fixture "$@"
}
run_case qc_prepare --input_pfile "$base" --prepare_pgs true --run_scores false --score_rsid_map "$map"
run_case qc_prepare_direct --input_pfile "$base" --prepare_pgs true --run_scores false --score_rsid_map "$map" --direct_inputs true
run_case qc_both --input_pfile "$base" --run_scores true --score_rsid_map "$map" --score_sheet "$scores"
run_case qc_both_direct --input_pfile "$base" --run_scores true --score_rsid_map "$map" --score_sheet "$scores" --direct_inputs true
run_case qc_reuse --input_pgs_pfile "$prepared" --run_scores false
run_case qc_reuse_direct --input_pgs_pfile "$prepared" --run_scores false --direct_inputs true
# Stale but valid PGS files must not redirect a base-only invocation.
for name in qc_base qc_base_direct; do
 mkdir -p "$root/$name/results/04_score_input"
 cp "$prepared".* "$root/$name/results/04_score_input/"
done
run_case qc_base --input_pfile "$base" --prepare_pgs false --run_scores false
run_case qc_base_direct --input_pfile "$base" --prepare_pgs false --run_scores false --direct_inputs true
python3 "$repo/tests/integration/validate_summary_qc.py" "$root"
