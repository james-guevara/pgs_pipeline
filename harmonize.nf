#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_manifest = null
params.cohort = 'cohort'
params.outdir = 'results/harmonization'
params.harmonization_container = null

include { PGEN_HARMONIZATION_WORKFLOW } from './workflows/harmonization'

workflow {
    if (!params.input_manifest) {
        error '--input_manifest is required'
    }
    if (!params.harmonization_container) {
        error '--harmonization_container is required and must contain Python, numpy, pgenlib, pyarrow, and procps'
    }
    PGEN_HARMONIZATION_WORKFLOW(
        Channel.value(file(params.input_manifest, checkIfExists: true)),
        Channel.value(file("${projectDir}/bin/harmonize_pgen_samples.py", checkIfExists: true)),
        Channel.value(file("${projectDir}/bin/merge_pgen_samples.py", checkIfExists: true))
    )
}
