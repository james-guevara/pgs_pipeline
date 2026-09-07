# Optional multi-fileset harmonization

This preparation workflow is needed only when participants from multiple
genotype filesets must be analyzed with an identical eligible marker set. It is
not part of ordinary single-fileset scoring.

`harmonize_cohorts.nf` accepts a TSV with `dataset_id`, `pgen_prefix`, and
`priority`. It applies sample missingness QC independently, calculates MAF after
those removals, and retains markers which are represented by the same
chromosome, position, REF, and ALT and have MAF >= the configured threshold in
every dataset. The default thresholds are sample missingness <= 0.05 and MAF >=
0.01.

The durable release contains `common_markers.txt`, an auditable
`common_variants.tsv` with per-dataset IDs and MAFs, per-dataset sample removal
files, and `harmonization_manifest.json`. PGS weight matching by rsID and effect
allele is deliberately downstream and is not performed here.

The authoritative base PGEN files are never modified. Scoring and ancestry can
optionally consume the common marker list. In harmonized mode, downstream MAF
filtering is repeated as a cheap guard and must not silently change the marker
set.

Harmonization reads each PVAR and runs only PLINK missingness/frequency reports;
it does not copy or rewrite the genotype matrices. The repository Dockerfile
pins Python, PLINK 2, and PLINK 1.9 in one runtime; PLINK 1.9 is used only for
the optional participant-wise PCA BED merge.

For PCA across multiple filesets, the matched common SNPs may additionally be
hard-called to BED and merged across participants. This PCA-only derivative is
appropriate for fitting shared within-ancestry axes; it must not replace the
dosage-preserving PGEN inputs used for PGS.
