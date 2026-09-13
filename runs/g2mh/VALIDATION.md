# G2MH release validation

The clean G2MH WGS-plus-imputed PGS and ancestry release completed on
2026-09-12 and is validated for reuse.

- Canonical release: `/expanse/projects/sebat1/s3/data/sebat/G2MH/combined/pgs-ancestry`
- Workflow commit: `e9473f98eb75e9d00b862bca50b5d0bd37f2fc89`
- Nextflow: `26.04.6`
- Participants: 1,784 (1,065 WGS and 719 imputed)
- Harmonized markers: 7,736,319
- Imputed dosage: preserved
- PGS: all 20 technical QC checks passed
- PCA panel overlap: 65,782/111,276 (59.12%; required minimum 50%)
- Global and within-ancestry PCA: completed

The release was compared with the independently completed R2=0.80 run from
2026-09-11. The genotype substrate, QC tables, all PGS outputs, all PCA and
ancestry outputs, final participant table, and data dictionary were
byte-for-byte identical. The machine-readable record and output checksums are
in `release-validation-2026-09-12.json`.
