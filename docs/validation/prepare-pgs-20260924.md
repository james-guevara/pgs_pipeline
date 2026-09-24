# Independent PGS preparation validation — 2026-09-24

Base: pgs_pipeline main 90e4b23. Implementation: feat/prepare-pgs-fileset.

Local regression suite: 34 passed, 1 skipped. The skipped test requires Bash 4
mapfile on macOS; score collation was exercised by the Linux integration test.

Real-container integration: PASS, Expanse Slurm job 54435170, Nextflow 26.04.6,
PLINK 2.0.0-a.6.9. Synthetic fixture: 120 samples, 20 biallelic SNPs, of which
18 survive the existing MAF policy. Ran seven independent Nextflow invocations:

- Preparation only, with rsID map and copy publication.
- Preparation only, direct input and hardlink publication (inode verified).
- Base input → preparation → scoring.
- Prepared input → scoring, staged input.
- Prepared input → scoring, direct input.
- Neither preparation nor scoring enabled, with summary QC disabled (no tasks).
- Preparation without an rsID map.

Validated native PVAR.ZST artifacts, expected counts, preparation provenance,
identical scoring tables for fresh preparation versus prepared-input reuse,
and identical retained genotype calls/sample order compared with the base.
Trace assertions verify that reuse runs do not execute preparation and that
preparation-only runs do not execute scoring or PCA.

Reports and traces:
`/expanse/lustre/scratch/j3guevar/temp_project/pgs-prepare-test-20260924/node-local-run-54435170/`

The test used node-local temporary storage for execution and copied its results
back to scratch. An earlier 12-sample test hit PLINK's small-sample frequency
estimation safeguard; increasing the fixture size resolved that without
changing production settings. An initial Lustre test was stopped to avoid slow
filesystem setup. Neither production WES nor SSC jobs was modified.

This test does not rerun PCA/ancestry or validate real-cohort rsID-map coverage.
Prepared-input match explanations are relative to the prepared fileset, not the
original base; preparation reports retain the earlier filtering information.
