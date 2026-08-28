#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 11 ]]; then
    echo "usage: $0 PFILE ANCESTRY_TSV MIN_SAMPLES NUM_PCS KING_CUTOFF LD_WINDOW LD_STEP LD_R2 CPUS MEMORY_MB OUTPUT_DIR" >&2
    exit 2
fi

pfile=$1
ancestry_tsv=$2
min_samples=$3
num_pcs=$4
king_cutoff=$5
ld_window=$6
ld_step=$7
ld_r2=$8
cpus=$9
memory_mb=${10}
output_dir=${11}

mkdir -p "$output_dir"
status_file="$output_dir/status.tsv"
printf 'ancestry\tassigned_samples\tunrelated_training_samples\tpruned_variants\tpcs\tstatus\treliability\treason\n' > "$status_file"

for ancestry in AFR AMR EAS EUR SAS; do
    group_dir="$output_dir/$ancestry"
    mkdir -p "$group_dir"
    keep_file="$group_dir/keep.tsv"
    printf '#IID\n' > "$keep_file"
    awk -F '\t' -v ancestry="$ancestry" 'NR > 1 && $2 == ancestry {gsub(/\r/, "", $1); print $1}' "$ancestry_tsv" >> "$keep_file"
    assigned=$(( $(wc -l < "$keep_file") - 1 ))

    reliability=reliable
    if (( assigned < min_samples )); then
        reliability=unreliable_small_sample
    fi
    printf 'ancestry\tassigned_samples\treliability_threshold\treliability\n%s\t%d\t%d\t%s\n' \
      "$ancestry" "$assigned" "$min_samples" "$reliability" > "$group_dir/reliability.tsv"

    if (( assigned < 3 )); then
        printf 'ancestry\tassigned_samples\tminimum_computable_samples\treason\n%s\t%d\t3\tinsufficient_samples_for_pca\n' \
          "$ancestry" "$assigned" > "$group_dir/skipped.tsv"
        printf '%s\t%d\t0\t0\t0\tskipped\tnot_computable\tinsufficient_samples_for_pca\n' "$ancestry" "$assigned" >> "$status_file"
        continue
    fi

    plink2 \
      --pfile "$pfile" \
      --keep "$keep_file" \
      --indep-pairwise "$ld_window" "$ld_step" "$ld_r2" \
      --out "$group_dir/prune" \
      --threads "$cpus" \
      --memory "$memory_mb"

    pruned_variants=$(wc -l < "$group_dir/prune.prune.in")
    if (( pruned_variants == 0 )); then
        printf 'ancestry\tassigned_samples\tminimum_samples\treason\n%s\t%d\t%d\tno_variants_after_ld_pruning\n' \
          "$ancestry" "$assigned" "$min_samples" > "$group_dir/skipped.tsv"
        printf '%s\t%d\t0\t0\t0\tskipped\tnot_computable\tno_variants_after_ld_pruning\n' "$ancestry" "$assigned" >> "$status_file"
        continue
    fi

    plink2 \
      --pfile "$pfile" \
      --keep "$keep_file" \
      --extract "$group_dir/prune.prune.in" \
      --king-cutoff "$king_cutoff" \
      --out "$group_dir/king" \
      --threads "$cpus" \
      --memory "$memory_mb"

    unrelated=$(wc -l < "$group_dir/king.king.cutoff.in.id")
    if (( unrelated < 3 )); then
        printf 'ancestry\tassigned_samples\tunrelated_training_samples\treason\n%s\t%d\t%d\tinsufficient_unrelated_training_samples\n' \
          "$ancestry" "$assigned" "$unrelated" > "$group_dir/skipped.tsv"
        printf '%s\t%d\t%d\t%d\t0\tskipped\tnot_computable\tinsufficient_unrelated_training_samples\n' \
          "$ancestry" "$assigned" "$unrelated" "$pruned_variants" >> "$status_file"
        continue
    fi

    group_pcs=$num_pcs
    if (( group_pcs >= unrelated )); then
        group_pcs=$(( unrelated - 1 ))
    fi
    last_pc=$(( 5 + group_pcs ))

    plink2 \
      --pfile "$pfile" \
      --keep "$group_dir/king.king.cutoff.in.id" \
      --extract "$group_dir/prune.prune.in" \
      --freq counts \
      --out "$group_dir/training" \
      --threads "$cpus" \
      --memory "$memory_mb"

    awk 'NR > 1 && $5 > 0 && $5 < $6 {print $2}' "$group_dir/training.acount" > "$group_dir/training_polymorphic.ids"
    pca_variants=$(wc -l < "$group_dir/training_polymorphic.ids")
    if (( pca_variants == 0 )); then
        printf 'ancestry\tassigned_samples\tunrelated_training_samples\treason\n%s\t%d\t%d\tno_polymorphic_training_variants\n' \
          "$ancestry" "$assigned" "$unrelated" > "$group_dir/skipped.tsv"
        printf '%s\t%d\t%d\t0\t0\tskipped\tnot_computable\tno_polymorphic_training_variants\n' \
          "$ancestry" "$assigned" "$unrelated" >> "$status_file"
        continue
    fi

    plink2 \
      --pfile "$pfile" \
      --keep "$group_dir/king.king.cutoff.in.id" \
      --extract "$group_dir/training_polymorphic.ids" \
      --pca allele-wts "$group_pcs" \
      --out "$group_dir/training" \
      --threads "$cpus" \
      --memory "$memory_mb"

    plink2 \
      --pfile "$pfile" \
      --keep "$keep_file" \
      --extract "$group_dir/training_polymorphic.ids" \
      --read-freq "$group_dir/training.acount" \
      --score "$group_dir/training.eigenvec.allele" 2 5 header-read no-mean-imputation variance-standardize \
      --score-col-nums "6-$last_pc" \
      --out "$group_dir/projected" \
      --threads "$cpus" \
      --memory "$memory_mb"

    cp "$group_dir/projected.sscore" "$group_dir/pcs.tsv"
    printf '%s\t%d\t%d\t%d\t%d\tcompleted\t%s\t.\n' \
      "$ancestry" "$assigned" "$unrelated" "$pca_variants" "$group_pcs" "$reliability" >> "$status_file"
done
