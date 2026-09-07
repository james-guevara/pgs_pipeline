#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.dataset_manifest = null
params.outdir = 'harmonized-release'
params.maf = 0.01
params.sample_miss = 0.05
params.python_container = 'python:3.11.16-bookworm'

process HARMONIZE_COHORTS {
    tag 'common-marker-release'
    label 'large'
    publishDir params.outdir, mode: 'copy'

    input:
    path dataset_manifest
    path harmonizer

    output:
    path 'release/common_markers.txt'
    path 'release/common_variants.tsv'
    path 'release/harmonization_manifest.json'
    path 'release/datasets'

    script:
    """
    python ${harmonizer} \
      --datasets ${dataset_manifest} \
      --output-dir release \
      --maf ${params.maf} \
      --sample-miss ${params.sample_miss}
    """
}

workflow {
    if (!params.dataset_manifest) {
        error '--dataset_manifest is required'
    }
    HARMONIZE_COHORTS(
        file(params.dataset_manifest, checkIfExists: true),
        file("${projectDir}/bin/harmonize_cohorts.py", checkIfExists: true)
    )
}
