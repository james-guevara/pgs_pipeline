#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.vcfs = null
params.rsid_maps = null
params.input_pfile = null
params.chromosomes = '1..22'
params.cohort = 'cohort'
params.genome_build = 'GRCh38'
params.outdir = 'results'
params.score_sheet = null
params.run_scores = false
params.run_summary_qc = true
params.min_score_variant_match = 0.50
params.warn_score_variant_match = 0.80
params.run_pca = false
params.pca_reference_sheet = null
params.r2 = null
params.aq = null
params.mac = 10
params.geno = 0.05
params.maf = 0.01
params.sample_miss = 0.05
params.variant_miss = 0.05
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
params.skip_rsid_annotation = false

include { PGS_WORKFLOW } from './workflows/pgs'

workflow {
    PGS_WORKFLOW()
}
