# Fixed-reference PCA and ancestry contract

PCA is a required final workflow stage. It projects one already coherent cohort
onto a versioned 1000 Genomes reference and then calculates within-ancestry PCs.
Dataset-specific preparation such as G2MH WGS/GSA harmonization is upstream.

## Reference manifest

`examples/pca_reference.tsv` defines one immutable reference release. Its fields
are:

- `reference_id`: versioned release identifier.
- `genome_build`: build shared by the cohort and reference artifacts.
- `panel_variants`: canonical chromosome, position, ID, REF, ALT, and loading
  allele for every fixed-panel SNP.
- `allele_frequencies`: 1000 Genomes frequencies used for standardization.
- `loadings`: PLINK2 allele-weight file for the global PCs.
- `reference_scores`: 1000 Genomes coordinates in the same projected space.
- `ancestry_labels`: labeled 1000 Genomes training samples.
- `classifier`: versioned ancestry model and decision thresholds.
- `checksums`: SHA-256 manifest covering every reference artifact.

Paths may be staged object-store URIs or shared POSIX paths. The reference files
are inputs; the workflow must never modify them.

## Required process boundaries

1. **Validate reference**: verify the manifest, build, checksums recorded by the
   reference release, PC count, and matching sample IDs across scores and labels.
2. **Harmonize panel**: intersect cohort variants with `panel_variants` by build,
   chromosome, position, and alleles. Align dosage to the loading allele; record
   swaps and reject incompatible or unresolved strand-ambiguous variants.
3. **Check overlap**: report total panel SNPs, usable SNPs, and exclusions by
   reason. Fail below `min_pca_variant_overlap`.
4. **Project global PCs**: use reference allele frequencies and loadings, never
   cohort frequencies, and emit the exact used-variant list.
5. **Assign ancestry**: apply the versioned classifier to projected cohort PCs;
   retain per-class probabilities as well as the final label.
6. **Within-ancestry PCA**: for each assigned group with at least
   `min_ancestry_samples`, perform LD pruning, relatedness handling, PCA on an
   unrelated training subset, and projection of all samples in that group.
   Groups below the threshold receive a machine-readable skipped result.

If panel overlap is materially below the level on which the fixed classifier was
validated, the pipeline must stop. A separate reference-building workflow should
reproject 1000 Genomes and refit/recalibrate the classifier on the shared subset;
the cohort workflow should not silently create a cohort-specific reference.

## Harmonization outputs

The workflow must publish:

- `panel_overlap.tsv`: one row per reference SNP with status and reason.
- `panel_overlap_summary.tsv`: counts and fractions for matched, swapped,
  ambiguous, incompatible, duplicate, absent, and usable SNPs.
- `projection_variants.txt`: exact variants used for projection.
- `global_pcs.tsv`: sample IDs and projected PCs.
- `ancestry_probabilities.tsv`: sample IDs, final labels, and all class
  probabilities.
- `within_ancestry/<group>/pcs.tsv`, or `skipped.tsv` with sample count and
  threshold.
- Reference ID, file checksums, software/container versions, and effective
  parameters in the run provenance.

Raw PGS values are separate outputs. Scaling and residualization using
within-ancestry PCs belong to downstream analysis rather than the scoring step.
