# Multi-fileset harmonization smoke test (2026-09-07)

Commit `15577e8` was tested on the persistent AWS host with the repository's
combined Python/PLINK container. Two derivative PGEN filesets were made from
real G2MH chr22 data: 500 variants in fileset A and 499 in fileset B.

The harmonizer completed sample-missingness measurement, removed samples above
the configured 0.05 threshold for frequency estimation, applied MAF >= 0.01 in
each fileset, and produced a frozen intersection of 400 markers. The detailed
table contained 401 lines including its header, and the JSON receipt reported
two datasets and 400 common markers. No genotype matrix was rewritten by the
harmonization operation.

Local validation also included seven Python/Nextflow structure tests and
successful Nextflow 26.04.6 previews for `main.nf`, `scoring.nf`, `ancestry.nf`,
`harmonize_cohorts.nf`, and `joint_pca_input.nf`.

The smoke test establishes the common-marker algorithm and runtime contents. A
full G2MH WGS/imputed release remains pending identification and staging of the
imputed PGEN input and the finalized genome-wide rsID map.
