# Summary QC target routing

On 2026-09-25, Expanse Slurm job **54454885** passed the expanded synthetic
integration suite using Nextflow 26.04.6 and the pinned PLINK 2.0.0-a.6.9 container.
The fixture contains 120 samples and 20 base variants; preparation retains 18
variants after MAF filtering and maps their IDs to rsIDs.

`tests/integration/summary_qc.sh` runs the seven existing preparation/scoring
cases and eight new summary-QC cases: preparation only, preparation implied by
scoring, supplied prepared input without scoring, and base-only input, each with
staged and direct input modes. `validate_summary_qc.py` verifies:

- Missingness, HWE, and frequency reports match independent native PLINK reports
  byte-for-byte for the selected fileset.
- PGS reports contain 18 mapped variants; base reports contain 20 coordinate-ID
  variants. Samples remain 120 in both cases.
- Provenance identifies the base/PGS target, actual source prefix/paths, and build.
- Freshly prepared outputs use staged QC even when preparation reads direct input.
- Existing valid PGS outputs cannot redirect a base-only invocation.
- Scoring outputs remain unchanged and disabled QC produces no reports.

Local checks: **34 passed, 1 skipped**. The skipped shell-collation test requires
Bash 4 `mapfile`, unavailable in the Mac's default Bash. The container integration
suite exercises collation through actual scoring. Existing PCA input guards pass;
PCA code and its base channel were not modified. QC and scoring consume sibling
channels after preparation; neither is gated on the other.

Evidence (synthetic only):

```text
/expanse/lustre/scratch/j3guevar/temp_project/pgs-summary-test-20260925/run-54454885.tar.gz
```

The archive contains `run/summary-qc-validation.json`, the preparation validation,
Nextflow traces and task logs. No production cohort jobs or artifacts were changed.
