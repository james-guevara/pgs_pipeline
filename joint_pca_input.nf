#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.dataset_manifest = null
params.common_markers = null
params.harmonization_dir = null
params.outdir = 'joint-pca-input'

process MAKE_MATCHED_BED {
    tag dataset_id
    label 'large'

    input:
    tuple val(dataset_id), val(pgen_prefix), val(priority)
    path common_markers
    val harmonization_dir

    output:
    tuple val(dataset_id), val(priority), path("${dataset_id}.bed"),
      path("${dataset_id}.bim"), path("${dataset_id}.fam")

    script:
    def memMb = Math.max(1000, task.memory.toMega() - 2000)
    """
    test -r '${pgen_prefix}.pgen' && test -r '${pgen_prefix}.pvar' && test -r '${pgen_prefix}.psam'
    remove_args=()
    removal='${harmonization_dir}/datasets/${dataset_id}/remove_samples.tsv'
    [[ -s "\$removal" ]] && remove_args=(--remove "\$removal")
    plink2 --pfile '${pgen_prefix}' --extract ${common_markers} "\${remove_args[@]}" \
      --make-bed --indiv-sort none --out '${dataset_id}' \
      --threads ${task.cpus} --memory ${memMb}
    """
}

process MERGE_PCA_BEDS {
    tag 'joint-pca-input'
    label 'large'
    publishDir params.outdir, mode: 'copy'

    input:
    path bed_files
    path dataset_manifest

    output:
    path 'joint_pca.bed'
    path 'joint_pca.bim'
    path 'joint_pca.fam'
    path 'participant_source.tsv'
    path 'overlapping_participants.tsv'

    script:
    """
    tail -n +2 ${dataset_manifest} | sort -t $'\t' -k3,3n -k1,1 > ordered.tsv
    : > seen.ids
    : > merge_prefixes.txt
    printf 'FID\tIID\tdataset_id\tpriority\n' > participant_source.tsv
    printf 'FID\tIID\tretained_dataset\texcluded_dataset\n' > overlapping_participants.tsv
    while IFS=$'\t' read -r dataset_id pgen_prefix priority rest; do
      awk 'NR==FNR {seen[\$1 SUBSEP \$2]=1; next} !seen[\$1 SUBSEP \$2] {print \$1, \$2}' \
        seen.ids "\${dataset_id}.fam" > "\${dataset_id}.keep"
      awk 'NR==FNR {owner[\$1 SUBSEP \$2]=\$3; next} owner[\$1 SUBSEP \$2] {print \$1"\t"\$2"\t"owner[\$1 SUBSEP \$2]"\t"dataset}' \
        participant_source.tsv dataset="\$dataset_id" "\${dataset_id}.fam" >> overlapping_participants.tsv
      if [[ ! -s "\${dataset_id}.keep" ]]; then
        continue
      fi
      plink --bfile "\$dataset_id" --keep "\${dataset_id}.keep" --make-bed --out "\${dataset_id}.dedup"
      awk -v dataset="\$dataset_id" -v priority="\$priority" '{print \$1"\t"\$2"\t"dataset"\t"priority}' \
        "\${dataset_id}.dedup.fam" >> participant_source.tsv
      awk '{print \$1, \$2}' "\${dataset_id}.dedup.fam" >> seen.ids
      printf '%s\n' "\${dataset_id}.dedup" >> merge_prefixes.txt
    done < ordered.tsv
    [[ -s merge_prefixes.txt ]] || { echo 'No samples remain for joint PCA input' >&2; exit 1; }
    first=\$(head -n 1 merge_prefixes.txt)
    tail -n +2 merge_prefixes.txt > remaining.txt
    if [[ -s remaining.txt ]]; then
      plink --bfile "\$first" --merge-list remaining.txt --make-bed --out joint_pca
    else
      cp "\${first}.bed" joint_pca.bed
      cp "\${first}.bim" joint_pca.bim
      cp "\${first}.fam" joint_pca.fam
    fi
    """
}

workflow {
    if (!params.dataset_manifest || !params.common_markers || !params.harmonization_dir) {
        error '--dataset_manifest, --common_markers, and --harmonization_dir are required'
    }
    datasets = Channel.fromPath(params.dataset_manifest, checkIfExists: true)
        .splitCsv(header: true, sep: '\t', strip: true)
        .map { row -> tuple(row.dataset_id.toString(), row.pgen_prefix.toString(), row.priority.toString().toInteger()) }
    common = Channel.value(file(params.common_markers, checkIfExists: true))
    MAKE_MATCHED_BED(datasets, common, params.harmonization_dir.toString())
    allBeds = MAKE_MATCHED_BED.out.flatten().filter { it instanceof Path }.collect()
    MERGE_PCA_BEDS(allBeds, file(params.dataset_manifest, checkIfExists: true))
}
