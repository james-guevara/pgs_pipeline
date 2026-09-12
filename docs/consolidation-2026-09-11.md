# Workflow consolidation (2026-09-11)

This repository now uses the workflow line validated by the combined G2MH WGS
and imputed-data run as its authoritative implementation. That line ended at
commit `510c9b0` before consolidation.

## Canonical entry points

- `main.nf`: PLINK preparation, PGS, ancestry/PCA, or their combined workflow.
  Use `--run_scores` and `--run_pca` to select the desired branches.
- `harmonize.nf`: optional dosage-preserving harmonization across genotype
  sources.

Running `main.nf` with only `--run_pca true` is the ancestry-only workflow;
running it with only `--run_scores true` is the scoring-only workflow. Separate
top-level copies of those same implementations are intentionally not retained.

## Consolidated history

The former `refactor/reusable-subworkflows-main` branch contained an independent
rewrite of `workflows/pgs.nf`. Its resource registry was retained, while its
older duplicate workflow, older cohort harmonizer, and convenience wrappers
were superseded by the G2MH-validated implementation and its flags. The newer
dosage-preserving harmonizer, PLINK schema-safe parsers, sample-ID normalization,
score matching diagnostics, and G2MH-tested execution behavior are retained.

The expanded Expanse score catalog and the within-ancestry HWE sensitivity
driver were added during consolidation. The latter is an experiment; no HWE
threshold becomes a production default until its results are reviewed.

## Checkout policy

`/Users/jamesguevara/pgs_pipeline` is the canonical local checkout. Do not use
or recreate `/Users/jamesguevara/projects/pgs_pipeline`. New work starts from
GitHub `main` and is merged back through a short-lived branch.
