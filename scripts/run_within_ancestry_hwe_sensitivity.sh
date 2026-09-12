#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
    echo "usage: $0 PFILE EXISTING_WITHIN_ANCESTRY_DIR PLINK2_CONTAINER THREADS OUTPUT_DIR" >&2
    exit 2
fi

pfile=$1
existing_dir=$2
plink2_container=$3
threads=$4
output_dir=$5
plink2=(singularity exec -B /expanse/projects/sebat1 "$plink2_container" plink2)

mkdir -p "$output_dir"
printf 'ancestry\tcondition\thwe_threshold\tassigned_samples\ttraining_samples\tvariants_after_hwe\tvariants_after_ld_pruning\tpcs\n' > "$output_dir/summary.tsv"

run_group() {
    local ancestry=$1
    local group_dir="$output_dir/$ancestry"
    local assigned_keep="$existing_dir/$ancestry/keep.tsv"
    local training_keep="$existing_dir/$ancestry/king.king.cutoff.in.id"
    local assigned training bad_ld bad_freq
    mkdir -p "$group_dir"

    assigned=$(( $(wc -l < "$assigned_keep") - 1 ))
    training=$(( $(wc -l < "$training_keep") - 1 ))
    bad_ld=()
    bad_freq=()
    if (( training < 50 )); then
        bad_ld=(--bad-ld)
        bad_freq=(--bad-freqs)
    fi

    for condition in none hwe_1e-6 hwe_1e-10 hwe_1e-15; do
        local condition_dir="$group_dir/$condition"
        local threshold=none
        local extract_args=()
        local variants_after_hwe
        mkdir -p "$condition_dir"

        if [[ "$condition" != none ]]; then
            threshold=${condition#hwe_}
            "${plink2[@]}" \
              --pfile "$pfile" \
              --keep "$training_keep" \
              --hwe "$threshold" 0 midp \
              --write-snplist \
              --out "$condition_dir/hwe" \
              --threads "$threads"
            extract_args=(--extract "$condition_dir/hwe.snplist")
            variants_after_hwe=$(wc -l < "$condition_dir/hwe.snplist")
        else
            variants_after_hwe=$(awk 'BEGIN{n=0} !/^#/ {n++} END{print n}' "$pfile.pvar")
        fi

        "${plink2[@]}" \
          --pfile "$pfile" \
          --keep "$training_keep" \
          "${extract_args[@]}" \
          "${bad_ld[@]}" \
          --indep-pairwise 200 50 0.2 \
          --out "$condition_dir/prune" \
          --threads "$threads"

        "${plink2[@]}" \
          --pfile "$pfile" \
          --keep "$training_keep" \
          --extract "$condition_dir/prune.prune.in" \
          --freq counts \
          --out "$condition_dir/training" \
          --threads "$threads"

        awk 'BEGIN { FS=OFS="\t" }
             NR == 1 {
                 for (i=1; i<=NF; i++) {
                     if ($i == "ID") id_col=i
                     else if ($i == "ALT_CTS") alt_cts_col=i
                     else if ($i == "OBS_CT") obs_ct_col=i
                 }
                 if (!id_col || !alt_cts_col || !obs_ct_col) exit 2
                 next
             }
             $alt_cts_col > 0 && $alt_cts_col < $obs_ct_col { print $id_col }' \
          "$condition_dir/training.acount" > "$condition_dir/training_polymorphic.ids"

        "${plink2[@]}" \
          --pfile "$pfile" \
          --keep "$training_keep" \
          --extract "$condition_dir/training_polymorphic.ids" \
          "${bad_freq[@]}" \
          --pca allele-wts 10 \
          --out "$condition_dir/training" \
          --threads "$threads"

        read -r score_id_col score_allele_col first_score_col last_score_col < <(
          awk '
            NR == 1 {
              for (i = 1; i <= NF; i++) {
                name = $i; sub(/^#/, "", name)
                if (name == "ID") id_col = i
                if (name == "A1") allele_col = i
                if (name == "PC1") first_col = i
                if (name == "PC10") last_col = i
              }
              if (!id_col || !allele_col || !first_col || !last_col) exit 2
              print id_col, allele_col, first_col, last_col
            }
          ' "$condition_dir/training.eigenvec.allele"
        )

        "${plink2[@]}" \
          --pfile "$pfile" \
          --keep "$assigned_keep" \
          --extract "$condition_dir/training_polymorphic.ids" \
          --read-freq "$condition_dir/training.acount" \
          --score "$condition_dir/training.eigenvec.allele" "$score_id_col" "$score_allele_col" header-read no-mean-imputation variance-standardize \
          --score-col-nums "$first_score_col-$last_score_col" \
          --out "$condition_dir/projected" \
          --threads "$threads"

        cp "$condition_dir/projected.sscore" "$condition_dir/pcs.tsv"
        printf '%s\t%s\t%s\t%d\t%d\t%d\t%d\t10\n' \
          "$ancestry" "$condition" "$threshold" "$assigned" "$training" \
          "$variants_after_hwe" "$(wc -l < "$condition_dir/training_polymorphic.ids")" \
          >> "$output_dir/summary.${ancestry}.tsv"
    done
}

for ancestry in AFR AMR EAS EUR SAS; do
    run_group "$ancestry" &
done
wait

for ancestry in AFR AMR EAS EUR SAS; do
    cat "$output_dir/summary.${ancestry}.tsv" >> "$output_dir/summary.tsv"
    rm "$output_dir/summary.${ancestry}.tsv"
done
