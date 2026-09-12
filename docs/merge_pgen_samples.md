# Merge multiple cohorts with disjoint samples

`bin/merge_pgen_samples.py` uses PLINK's Python `pgenlib` reader/writer to
merge two or more `.pgen` + uncompressed `.pvar` + `.psam` filesets. This is a standalone
preparation step before passing the resulting prefix to `--input_pfile`.

Install in a virtual environment:

```sh
python3 -m venv .venv-merge
.venv-merge/bin/pip install pgenlib==0.94.1 numpy
.venv-merge/bin/python bin/merge_pgen_samples.py \
  --inputs /path/to/cohort_a /path/to/cohort_b /path/to/cohort_c \
  --dosage-mode preserve \
  --out-dir /path/to/new_merged_directory
```

The output prefix is `/path/to/new_merged_directory/merged`.
`merge_summary.json` records input paths, sample and variant counts, exclusions,
and phase handling. The `inputs` array contains counts and unmatched variants for
each fileset, relative to the final intersection. The output directory must not exist; failed runs do not
publish partial filesets or change inputs.

The original `--first A --second B` command remains supported for two inputs;
it cannot be combined with `--inputs`.

The merge keeps the intersection across **all inputs** of exact `(CHROM, POS, REF, ALT)` keys, in
first-input order, and generates `CHROM:POS:REF:ALT` IDs. Original variant IDs
need not agree. Chromosome strings must match exactly; no strand complementation,
allele swapping, chromosome renaming, or genome-build conversion is performed.
Use inputs on the same genome build; differing `##reference` metadata is rejected
among inputs that provide it, but missing metadata cannot establish the build.

All inputs must contain only biallelic A/C/G/T SNPs and unique variant keys.
Sample identity uses `(FID, IID)`, with omitted FID treated as `0`; duplicate or
overlapping identities are rejected even if SID differs. PSAM column names and
order must match across every input. All sample fields are copied in command-line
input order, preserving sample order within each fileset.

By default, missing hardcalls and hardcall phase are preserved, and retained
dosages which differ from hardcalls cause an error. With `--dosage-mode
preserve`, dosage vectors are concatenated and written to the output PGEN;
hardcall phase is not preserved in that mode. This is sufficient for the
biallelic SNP dosage inputs used by PGS and PCA, but is not a general merger for
all PGEN representations. Original
ploidy encoding is copied without chromosome-specific genotype corrections.
PVAR INFO/QUAL/FILTER and other cohort-specific annotations are omitted; the
output contains the five core columns and available reference metadata. REF
status is conservatively marked provisional in the PGEN.

Genotype buffers scale with sample count, not the whole genotype matrix.
Variant indexes and sample metadata are held in memory, so memory also scales
with input variant count. One PGEN reader remains open per input, so very large
input lists may hit operating-system file-descriptor or memory limits. The output
is written once, without intermediate pairwise merges.

## G2MH chr22 validation

The dosage-preserving path was validated on 2026-09-08 with 1,065 WGS samples
and 719 imputed samples. The imputed VCF was imported with `dosage=DS` after an
INFO `R2 >= 0.90` filter. Direct comparison of 72,619 nonmissing DS values at
101 chromosome-spanning variants found 33,951 fractional values and no
differences beyond PGEN precision (`1/16384`).

The sample-axis merge retained 83,965 exact shared variants and 1,784 samples
and completed in 3 seconds on Expanse. A second validation compared 180,183
source and merged dosages across 101 variants: sample order, variant
intersection/order, and every tested dosage were exact; 26,670 tested imputed
values were fractional.

Applying MAF >=1% after merging retained 79,206 variants. Requiring MAF >=1%
in both sources retained 76,955. All 2,251 combined-only variants had MAF
between 0.5% and 1% in the failing source; none was monomorphic. The median
absolute source-frequency difference was 0.00378, 95 variants differed by at
least 0.01, and 20 differed in missingness by at least 0.01. This supports a
combined-sample MAF policy for the unified G2MH analysis while retaining
source-specific frequencies as QC evidence.

The reusable `harmonize.nf` implementation was then validated on the same
inputs with Nextflow 26.04.6 and container digest
`sha256:67420bc37e066425bd85f3a26315a7213f233144409b54169104648edb155bfc`.
The task completed in 9.487 seconds and published 1,784 samples, 83,965 exact
shared markers, a 202,718-row marker-QC Parquet, and 1,784 sample-QC rows. Its
PGEN, PVAR, and PSAM were byte-identical to the independently validated
standalone merge. All retained marker-QC rows were present in both sources and
had combined MAF; all 118,753 excluded union rows had an exclusion reason.

The harmonized PGEN was also passed through the ordinary `main.nf` scoring
branch with direct PGEN input. The scoring-specific MAF >=1% view reproduced
the independently calculated reduction from 83,965 to 79,206 markers. A smoke
score matched all 5 requested variants and produced 1,784 participant scores;
the per-score and collated QC outputs both passed. This confirms that supplied
and harmonized PGENs now receive the same explicit scoring-view MAF policy,
without applying that PGS-specific filter to PCA.

Run synthetic round-trip and failure tests:

```sh
.venv-merge/bin/python -m unittest discover -s tests -p test_merge_pgen_samples.py -v
```
