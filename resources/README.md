# Pipeline resource registry

This directory records the versioned, reusable inputs required by the PGS and
PCA branches. Large genomic files are kept in shared storage rather than Git.
The manifests under `expanse/` and `aws/` describe the same logical resources
using environment-specific paths.

## Current bundle

`g2mh_provisional_v1` contains:

- the fixed `1kg_grch38_v1` PCA projection and ancestry-classification release;
- 20 provisional PGS weight files and their PLINK2 column definitions; and
- the pinned PLINK2 and Python container references used by the workflow.

The PGS weights are provisional and must not be interpreted as a finalized
scientific release. Cohort genotype files, rsID maps, and run outputs are not
part of this bundle: they are run-specific inputs or outputs.

Copies of this registry are published at:

- `/expanse/projects/sebat1/s3/data/sebat/pgs_pipeline/resources/g2mh_provisional_v1/`
- `s3://sebat-genomics-work/resources/pgs_pipeline/g2mh_provisional_v1/`

## Use

On Expanse:

```bash
nextflow run main.nf -profile slurm -c conf/expanse.config \
  --run_pca true \
  --pca_reference_sheet resources/expanse/pca_reference.tsv \
  --run_scores true \
  --score_sheet resources/expanse/scores.tsv
```

On AWS:

```bash
nextflow run main.nf -profile awsbatch \
  --run_pca true \
  --pca_reference_sheet resources/aws/pca_reference.tsv \
  --run_scores true \
  --score_sheet resources/aws/scores.tsv
```

`resource_catalog.tsv` is the human-readable inventory. The PCA release has its
own `SHA256SUMS`. A weight release should receive checksums and a new versioned
name when the weights are finalized; do not replace files in a named release.
