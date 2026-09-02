#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_pfile = null
params.pca_reference_sheet = null
params.cohort = 'cohort'
params.genome_build = 'GRCh38'
params.outdir = 'results'
params.num_pcs = 10
params.num_within_ancestry_pcs = 10
params.min_pca_variant_overlap = 0.90
params.min_ancestry_samples = 50
params.within_ancestry_king_cutoff = 0.0884
params.within_ancestry_ld_window = 200
params.within_ancestry_ld_step = 50
params.within_ancestry_ld_r2 = 0.2
params.python_container = 'python:3.11.16-bookworm'
params.fsx_direct = false
params.direct_inputs = false

include { ANCESTRY_WORKFLOW } from './workflows/pgs'

workflow {
    if (!params.input_pfile) {
        error '--input_pfile is required and must name a QCed PLINK 2 prefix'
    }

    def directPfileEnabled = (params.direct_inputs instanceof Boolean ? params.direct_inputs : params.direct_inputs?.toString()?.toBoolean()) ||
        (params.fsx_direct instanceof Boolean ? params.fsx_direct : params.fsx_direct?.toString()?.toBoolean())
    def pfilePrefix = params.input_pfile.toString()
    if (directPfileEnabled) {
        qcPfile = Channel.value(tuple(
            "${pfilePrefix}.pgen",
            "${pfilePrefix}.pvar",
            "${pfilePrefix}.psam"
        ))
    } else {
        qcPfile = Channel.value(tuple(
            file("${pfilePrefix}.pgen", checkIfExists: true),
            file("${pfilePrefix}.pvar", checkIfExists: true),
            file("${pfilePrefix}.psam", checkIfExists: true)
        ))
    }

    ANCESTRY_WORKFLOW(qcPfile, directPfileEnabled)
}
