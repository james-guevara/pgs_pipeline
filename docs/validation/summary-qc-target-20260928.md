# Summary-QC provenance follow-up (2026-09-28)

Review against the expanded requirements found that routing was already correct,
but the QC JSON lacked cohort, sample/variant counts, and preparation settings.
The follow-up adds `qc_target` (retaining `target`), `cohort`, numeric
`sample_count`/`variant_count`, and `preparation` containing the source input,
MAF threshold, and rsID-map path when preparation runs in this invocation.
Supplied PGS inputs have null preparation metadata: current defaults do not
establish their historical settings. No filesystem discovery was added.

Validation used Nextflow 26.04.6 and native macOS PLINK v2.0.0-a.7.9 M1
(27 Sep 2026), downloaded from the official PLINK site. Eight native cases pass:
base, prepared with mapping, prepared without mapping, and supplied PGS, each in
staged and direct mode. They check exact native report equivalence, 120 samples,
20 base versus 18 prepared variants, mapped IDs, source paths, all JSON fields,
and unknown supplied-input preparation settings despite a different current MAF.
The existing component suite also passes: 34 passed, one skipped for macOS's
Bash 3 lacking `mapfile`.

The existing 15-case container harness now asserts the additional JSON fields,
including preparation implied by scoring. It was not rerun for this follow-up;
the earlier container results apply to commit
`321b6904fd33f08671a1ccd7c068b7e4b0917ad9`, as recorded in
`summary-qc-target-20260925.md`. QC commands, routing, scoring, and PCA bindings
are unchanged by this follow-up. No real cohort validation was performed.

Reproduce the native checks with `pytest -q tests/test_summary_qc.py` when
Nextflow and PLINK2 are on PATH. Set `PGS_ENTRYPOINT` to a parent `pgs.nf` to run
the same checks through that entrypoint.
