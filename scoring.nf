#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_pfile = null
params.rsid_map = null
params.score_sheet = null
params.common_markers = null
params.cohort = 'cohort'
params.outdir = 'results'
params.maf = 0.01
params.sample_miss = 0.05
params.min_score_variant_match = 0.50
params.warn_score_variant_match = 0.80
params.direct_inputs = false
params.fsx_direct = false

include { SCORING_WORKFLOW } from './workflows/pgs'

workflow {
    if (!params.input_pfile || !params.rsid_map || !params.score_sheet) {
        error '--input_pfile, --rsid_map, and --score_sheet are required'
    }
    direct = (params.direct_inputs instanceof Boolean ? params.direct_inputs : params.direct_inputs.toString().toBoolean()) ||
        (params.fsx_direct instanceof Boolean ? params.fsx_direct : params.fsx_direct.toString().toBoolean())
    prefix = params.input_pfile.toString()
    if (direct) {
        pfile = Channel.value(tuple("${prefix}.pgen", "${prefix}.pvar", "${prefix}.psam"))
    } else {
        pfile = Channel.value(tuple(file("${prefix}.pgen", checkIfExists: true),
            file("${prefix}.pvar", checkIfExists: true), file("${prefix}.psam", checkIfExists: true)))
    }
    rsidMap = Channel.value(file(params.rsid_map, checkIfExists: true))
    commonMarkers = Channel.value(params.common_markers ?
        file(params.common_markers, checkIfExists: true) : file("${projectDir}/resources/empty_markers.txt"))
    SCORING_WORKFLOW(pfile, direct, rsidMap, commonMarkers)
}
