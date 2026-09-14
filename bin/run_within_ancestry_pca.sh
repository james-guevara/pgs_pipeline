#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 12 ]]; then
    echo "usage: $0 PFILE ANCESTRY_TSV MIN_SAMPLES NUM_PCS KING_CUTOFF LD_WINDOW LD_STEP LD_R2 MAF CPUS MEMORY_MB OUTPUT_DIR" >&2
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
maf=$9
cpus=${10}
memory_mb=${11}
output_dir=${12}

mkdir -p "$output_dir"
status_file="$output_dir/status.tsv"
printf 'ancestry\tassigned_samples\tunrelated_training_samples\tpruned_variants\tpcs\tstatus\treliability\treason\n' > "$status_file"

for ancestry in AFR AMR EAS EUR SAS; do
    group_dir="$output_dir/$ancestry"
    mkdir -p "$group_dir"
    keep_file="$group_dir/keep.tsv"
    # Preserve the confidence-qualified ANCESTRY field for reporting, but use
    # each sample's most likely group for within-ancestry PCA when available.
    # Build a two-column keep file from the source PSAM so nonzero FIDs are
    # handled correctly instead of silently dropping those samples.
    awk -F '\t' -v ancestry="$ancestry" '
      NR == FNR {
        if (FNR == 1) {
          for (i = 1; i <= NF; i++) {
            name = $i; sub(/^#/, "", name)
            if (name == "IID") iid_col = i
            if (name == "MOST_LIKELY_ANCESTRY") group_col = i
            if (name == "ANCESTRY") ancestry_col = i
          }
          if (!group_col) group_col = ancestry_col
          next
        }
        gsub(/\r/, "", $iid_col)
        if ($group_col == ancestry) wanted[$iid_col] = 1
        next
      }
      FNR == 1 {
        for (i = 1; i <= NF; i++) {
          name = $i; sub(/^#/, "", name)
          if (name == "FID") fid_col = i
          if (name == "IID") psam_iid_col = i
        }
        print "#FID\tIID"
        next
      }
      ($psam_iid_col in wanted) { print $fid_col "\t" $psam_iid_col }
    ' "$ancestry_tsv" "${pfile}.psam" > "$keep_file"
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

    ld_sample_override=()
    frequency_sample_override=()
    if (( assigned < 50 )); then
        # PLINK refuses LD estimation below 50 samples by default. The user has
        # explicitly requested best-effort PCs for these groups; their outputs
        # remain marked unreliable_small_sample.
        ld_sample_override=(--bad-ld)
        frequency_sample_override=(--bad-freqs)
    fi

    plink2 \
      --pfile "$pfile" \
      --keep "$keep_file" \
      --maf "$maf" \
      "${ld_sample_override[@]}" \
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

    unrelated=$(awk 'NR > 1 { n++ } END { print n + 0 }' "$group_dir/king.king.cutoff.in.id")
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
    plink2 \
      --pfile "$pfile" \
      --keep "$group_dir/king.king.cutoff.in.id" \
      --extract "$group_dir/prune.prune.in" \
      --freq counts \
      --out "$group_dir/training" \
      --threads "$cpus" \
      --memory "$memory_mb"

    # Resolve columns by name: PLINK releases may add fields (for example,
    # PROVISIONAL_REF?) before ALT_CTS. Positional parsing would then silently
    # classify every variant as monomorphic.
    awk '
      NR == 1 {
        for (i = 1; i <= NF; i++) {
          name = $i; sub(/^#/, "", name)
          if (name == "ID") id_col = i
          if (name == "ALT_CTS") alt_col = i
          if (name == "OBS_CT") obs_col = i
        }
        if (!id_col || !alt_col || !obs_col) exit 2
        next
      }
      $alt_col > 0 && $alt_col < $obs_col { print $id_col }
    ' "$group_dir/training.acount" > "$group_dir/training_polymorphic.ids"
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
      "${frequency_sample_override[@]}" \
      --pca allele-wts "$group_pcs" \
      --out "$group_dir/training" \
      --threads "$cpus" \
      --memory "$memory_mb"

    read -r score_id_col score_allele_col first_score_col last_score_col < <(
      awk -v final_pc="PC${group_pcs}" '
        NR == 1 {
          for (i = 1; i <= NF; i++) {
            name = $i; sub(/^#/, "", name)
            if (name == "ID") id_col = i
            if (name == "A1") allele_col = i
            if (name == "PC1") first_col = i
            if (name == final_pc) last_col = i
          }
          if (!id_col || !allele_col || !first_col || !last_col) exit 2
          print id_col, allele_col, first_col, last_col
        }
      ' "$group_dir/training.eigenvec.allele"
    )

    plink2 \
      --pfile "$pfile" \
      --keep "$keep_file" \
      --extract "$group_dir/training_polymorphic.ids" \
      --read-freq "$group_dir/training.acount" \
      --score "$group_dir/training.eigenvec.allele" "$score_id_col" "$score_allele_col" header-read no-mean-imputation variance-standardize \
      --score-col-nums "$first_score_col-$last_score_col" \
      --out "$group_dir/projected" \
      --threads "$cpus" \
      --memory "$memory_mb"

    cp "$group_dir/projected.sscore" "$group_dir/pcs.tsv"
    # Keep the exact PCA variant IDs, training sample IDs, eigenvalues, and
    # projected PCs. These larger files are regenerable task intermediates and
    # are not needed for downstream residualization or provenance.
    rm -f \
      "$group_dir/prune.prune.out" \
      "$group_dir/training.acount" \
      "$group_dir/training.eigenvec.allele" \
      "$group_dir/projected.sscore"
    printf '%s\t%d\t%d\t%d\t%d\tcompleted\t%s\t.\n' \
      "$ancestry" "$assigned" "$unrelated" "$pca_variants" "$group_pcs" "$reliability" >> "$status_file"
done
