# PCA panel matching: streaming PVAR reader

Date: 2026-10-12. Base: `c471fed2d49c14731b867167b97a1d0b346ec355`.
Local branch: `fix/stream-pca-panel`.

## Status

Prepared and locally tested; not pushed, merged, or deployed on Expanse.
PGS scoring and ancestry jobs have not been launched for new SPARK WGS.
Repository adoption is deferred for review, as requested by James.

## Motivation and implementation

The new SPARK WGS base has 57,002,418 variants and 45,178 samples. The Python
helper `bin/harmonize_pca_panel.py` previously materialized every PVAR row as a
dictionary, then built an index of the entire cohort by locus. This memory cost
scales with variant count, not sample count; this helper does not read genotypes.

The implementation uses Python's standard library, not Polars. `read_pvar`
now yields rows incrementally, optionally restricted to a set of panel loci.
The caller builds that set from the supplied reference panel and indexes only
cohort rows at those loci. All records at a panel locus remain available,
including incompatible and duplicate alleles. It still scans the full PVAR.

No cohort-specific rules, thresholds, allele matching, strand handling, duplicate
rules, output formats, or Nextflow process/channel contracts change. The existing
READ_PVAR_METADATA task still provides plain metadata from compressed PVAR input.
The helper has no new dependencies. Memory scales with the panel and cohort
records at panel loci; it is not strictly constant for unlimited duplicates.

## Validation

- Four regression tests pass (`tests/test_harmonize_pca_panel.py`).
- Covers exact, REF/ALT-swapped, duplicate-compatible, strand-complement-rejected,
  incompatible, absent, and duplicate-panel-ID cases; overlap gate and missing header.
- A 200,000-row synthetic PVAR retains two requested rows with peak traced Python
  allocation below 2 MiB. This is a regression bound, not real-cohort RSS measurement.
- All five output files are byte-identical to the baseline helper on the matching fixture.
- No real SPARK timing or peak-memory benchmark has been run yet.

Command:

```bash
python3 -m pytest -q tests/test_harmonize_pca_panel.py
```

## Follow-up before production / repository adoption

Review the diff, pin the resulting commit in the run snapshot, and measure actual
panel overlap, elapsed time, and peak RSS. Preserve the base PGEN and reference
checksums and all matching reports. Do not silently lower the panel overlap gate.

The separate `EXPLAIN_SCORE_MATCHING` helper is now addressed by the follow-up
recorded below. These two changes do not establish memory bounds for other workflow
processes.

## Other deferred change: FastVEP

FastVEP draft PR: https://github.com/james-guevara/fastVEP/pull/1
Commit: `b1b79335af455a4f300c9bd6fa68274f2bfb94ae`.
This repairs a CDS-coordinate/sequence mismatch for nonzero starting coding phase
in frameshift HGVSp construction. It was tested and used in the isolated SPARK
chr11 recovery container; adopting it as the standard pipeline pin remains deferred.
It is an annotation correction, separate from this matching-memory optimization.

Detailed evidence and scientific caveats are recorded in the rare-variant repo's
`docs/operations/zarr/SPARK_AUTOSOMES_20261010.md` on `feat/zarr-policy-backend`.
Other chromosomes retain their original annotation pin; this note does not claim
that all annotations have been regenerated with the patch.

## Follow-up: scoring-match report streaming

Implemented in `bin/explain_score_matching.py` on the same local branch;
not deployed or pushed. Nextflow arguments, channels, and stage boundaries are
unchanged. This helper produces reports; it does not modify genotypes or scores.

Previously it loaded the full rsID map, both PVAR tables, a second normalized
source dictionary, all weight rows, and all excluded-detail rows into memory.
Now it makes two streaming passes over weights (counts first, output second),
retains mappings only to requested weight IDs, scans source PVAR into a requested-ID
presence set, and stores scoring alleles only for requested IDs. Detail rows are
written directly. Memory scales with unique requested weight IDs and their relevant
map keys, not the total cohort variant count. Each trait still scans source metadata;
shared scans across traits and genotype I/O are not optimized by this change.

Semantics preserved: matched-ID priority, allele mismatch classes, MAF versus
mapping exclusion, no-map versus irrelevant-map distinction, duplicate weight-row
counts/order, last-key-wins mappings, and last-row-wins scoring alleles.

One explicit parser correction: map rows accept whitespace separators, including
spaces used by the canonical GRCh38 map. Previously the report parser recognized
only tabs, despite PLINK accepting the supplied space-delimited map. Existing TSV
behavior is preserved; reports for space-delimited maps are intentionally corrected.
This does not change actual PLINK scoring or filtering.

Validation:

- Eight targeted tests pass across both helpers and the existing scoring-report test.
- Both report files are byte-identical to baseline `c471fed` for the TSV fixture,
  covering all report categories, repeated weights, and overwritten map keys.
- Space-delimited and tab-delimited maps produce identical reports.
- Streaming 200,000 unrelated PVAR/map records with one requested ID peaks below
  2 MiB traced Python allocation. This is a synthetic regression bound, not a
  real SPARK RSS measurement or a promise for million-variant weight files.

```bash
python3 -m pytest -q tests/test_explain_score_matching.py \
  tests/test_score_matching_streaming.py tests/test_harmonize_pca_panel.py
```

Before adoption, review both helper changes, use a pinned run snapshot, and collect
real-data resource usage. Keep the FastVEP adoption decision separate.

## Deployment and user-approved overlap warning (2026-10-12)

Snapshot `a34e286` was deployed in the isolated SPARK scoring/ancestry run.
Controller 54862271 stopped because overlap was 89,739/111,276 (80.645422%)
against a 90% gate: 89,739 exact matches, 21,537 absent, no allele mismatches.
The user requested warning-only behavior. Low nonzero overlap now emits stderr
warning and report fields `overlap_policy=warn` and `low_overlap_warning`;
the legacy minimum_usable_fraction report key and parameter names remain for compatibility.
Zero usable markers remains an error. No matching rules or genotypes change.
The additional report fields intentionally extend the previous output format.
This supersedes the earlier instruction to preserve the fatal overlap gate.

Scoring continued independently with controller 54862278 using explicit session
0dfc629e-5148-4d8a-a586-49381880a566. Ancestry will use a separate pinned run snapshot.
Shared-repository adoption remains deferred.

## Follow-up: IID-only PSAM in within-ancestry PCA

SPARK's base PSAM has `#IID SEX`, a valid PLINK representation. The within-ancestry
keep-list AWK script assumed FID existed and evaluated `$0` when its FID column
index was unset, emitting the whole PSAM row as the first keep-list field. PLINK
therefore retained zero samples. Global projection and classification had completed.

The helper now emits FID=0 when the input has no FID, matching PLINK's implicit
family ID, and preserves explicit FIDs otherwise. It does not impute biological
family IDs or modify the base. Three regression cases cover IID-only, explicit
nonzero FID, and reordered columns using the actual shell helper's AWK program.
All twelve targeted tests pass. This additional generic fix remains on the local
pending-review branch; the ancestry run will use its pinned snapshot.
