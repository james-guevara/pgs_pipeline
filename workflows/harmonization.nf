nextflow.enable.dsl=2

process HARMONIZE_PGEN_SAMPLES {
    tag params.cohort
    label 'large'
    container params.harmonization_container
    publishDir params.outdir, mode: 'copy'

    input:
    path input_manifest
    path harmonizer
    path merger

    output:
    tuple path('harmonized/harmonized.pgen'), path('harmonized/harmonized.pvar'),
        path('harmonized/harmonized.psam'), emit: pfile
    path 'harmonized/harmonization_summary.json', emit: summary
    path 'harmonized/marker_qc.parquet', emit: marker_qc
    path 'harmonized/sample_qc.tsv', emit: sample_qc

    script:
    """
    python '${harmonizer}' --manifest '${input_manifest}' --out-dir harmonized
    """
}

workflow PGEN_HARMONIZATION_WORKFLOW {
    take:
    input_manifest
    harmonizer
    merger

    main:
    HARMONIZE_PGEN_SAMPLES(input_manifest, harmonizer, merger)

    emit:
    pfile = HARMONIZE_PGEN_SAMPLES.out.pfile
    summary = HARMONIZE_PGEN_SAMPLES.out.summary
    marker_qc = HARMONIZE_PGEN_SAMPLES.out.marker_qc
    sample_qc = HARMONIZE_PGEN_SAMPLES.out.sample_qc
}
