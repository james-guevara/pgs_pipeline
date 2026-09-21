# Plain/compressed PVAR input parity — 2026-09-21

Full PGS/PCA runs used identical 275-sample, 65,782-marker inputs with 100
synthetic scoring weights and the existing 1kg_grch38_v1 PCA reference. The
input fixture was reused from the earlier PGS reconciliation smoke test.
No participant-level data is included here.

Slurm job 54386156 ran on Expanse (one node, 8 CPUs, 24 GB allocation), with
Nextflow 26.04.6 and PLINK 2.0.0a.6.9 in the existing pinned Singularity image.
Each Nextflow process used at most 4 CPUs / 8 GB, local executor inside the
allocation, queue size 1. None of these three runs used resume:

1. Plain `.pvar`, standard staging (baseline).
2. Native `.pvar.zst`, standard staging.
3. Native `.pvar.zst`, direct shared-filesystem reads.

Compressed input was generated with `plink2 --pfile input --make-pgen vzs`.
Both compressed cases matched all 70 baseline scientific artifacts byte for
byte: summary QC, filtered score fileset, scores and explanations, global PCs,
ancestry probabilities, panel matching, all five within-ancestry PCA outputs,
and the final participant table. Runtime logs were excluded. See the aggregate
JSON alongside this file. The comparator is
`tests/integration/compare_pvar_runs.py` and fails on missing/different outputs.

Commands, all artifacts, Nextflow traces/reports, source snapshot, comparison
script, and aggregate validation are retained at:

`/expanse/projects/sebat1/j3guevar/integrated_genomics_pipeline/experiments/pgs-vzs-20260921`

## Reproduction

Use the same QCed input to create a compressed counterpart with pinned PLINK.
Run `main.nf` three times with separate work/output directories, the same score
sheet, rsID map, PCA reference, and site configuration, enabling
`--run_summary_qc true --run_scores true --run_pca true`. Use plain input with
`--direct_inputs false`, then compressed input with `false` and `true`.
The default `--input_pvar_format auto` should be used for this test. Place
results under `plain-staged/results`, `compressed-staged/results`, and
`compressed-direct/results`, then run:

```bash
python3 tests/integration/compare_pvar_runs.py /shared/experiment-root
```

This verifies input-storage equivalence; it does not independently validate
ancestry classifier accuracy or the scientific appropriateness of QC policy.
The original base fileset remains untouched; one work-directory metadata view
is produced for Python readers, while all genotype reads use native PLINK.
