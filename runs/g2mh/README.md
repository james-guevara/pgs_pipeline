# G2MH PGS and ancestry release

This directory records the cohort-specific inputs and parameters for the clean
G2MH WGS-plus-imputed rerun. The workflow code remains cohort-neutral.

## Fixed decisions

- Workflow revision: recorded from Git automatically when the release is run
- Nextflow: `26.04.6`
- Genome build: `GRCh38`
- Imputed-data inclusion threshold: `R2 >= 0.80`, applied while creating the
  reusable imputed base PGEN
- Imputed dosage field: `DS`, retained through sample-axis harmonization
- PGS scoring MAF: `0.01`
- Within-ancestry HWE filtering: disabled
- HWE statistics: reported by summary QC, never interpreted across pooled
  ancestry groups as a scientific exclusion criterion
- Global PCA minimum reference-panel overlap: `0.50` for this harmonized G2MH
  dataset; the validated overlap was 65,782/111,276 (59.12%)
- Within-ancestry PCA: 10 PCs, KING cutoff 0.0884, LD pruning 200/50/0.2
- Small ancestry groups: calculate when possible but flag groups with fewer
  than 50 assigned samples as unreliable

## Canonical products

The clean release belongs under:

`/expanse/projects/sebat1/s3/data/sebat/G2MH/combined/pgs-ancestry/`

The reusable inputs belong under `G2MH/WGS/plink-base/` and
`G2MH/imputed/plink-base/`. Do not delete the currently validated source
filesets until those canonical copies and the clean combined release have been
validated.

`harmonization-inputs.tsv` records the currently validated physical inputs.
`params.json` is passed to the downstream PGS/PCA workflow after harmonization.
`scores.tsv` is the provisional set of 20 currently scoring-ready weights. New
scores can be added later without rebuilding the base PGEN filesets.

The harmonizer accepts differing PSAM schemas at its input boundary and emits
one canonical `FID`, `IID`, `SEX`, then-extra-columns schema. Participant
overlap is checked by IID, not by an `(FID, IID)` pair.

## Validation gates

1. WGS, imputed, and harmonized sample counts are 1,065, 719, and 1,784.
2. WGS and imputed IID sets are disjoint.
3. The harmonized marker set is identical across all participants.
4. Fractional imputed dosages survive harmonization within PGEN precision.
5. The PCA reference overlap is at least 0.50.
6. Score matching QC is emitted for every requested score.
7. Global and within-ancestry PC outputs contain exactly the expected IIDs.
8. The participant analysis table has one row per IID and no missing-data-to-zero coercion.
9. Outputs are compared with the validated 2026-09-11 run before release.

The 2026-09-12 clean release passed all gates. See `VALIDATION.md` and
`release-validation-2026-09-12.json` for the regression result, job IDs,
timings, and checksums.
