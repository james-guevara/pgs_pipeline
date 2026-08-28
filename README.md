# PGS pipeline

Containerized Nextflow pipeline for VCF filtering, rsID annotation, genome-wide
QC, polygenic scoring, and a specified fixed-reference PCA interface. It can
run unchanged with local Docker or AWS Batch.

## Workflow

1. Run PLINK2 filtering, rsID annotation, and MAF filtering independently for each chromosome.
2. Concatenate chromosome pfiles into one cohort pfile.
3. Apply two-pass variant and sample missingness QC.
4. Produce missingness, Hardy-Weinberg, and allele-frequency summaries.
5. Optionally calculate PGS values from externally generated weight files.
6. Project onto a fixed 1000 Genomes PCA reference and assign ancestry with
   probabilities, then calculate within-ancestry PCs using unrelated training
   subsets and project every confidently assigned sample in each eligible group.

The original Python and Slurm entrypoints remain for reference. `main.nf` is the
portable orchestration entrypoint; cloud queues and filesystem mounts are
deployment profiles, not part of the scientific workflow.

## Required inputs

- One bgzipped VCF per chromosome.
- One two-column rsID map per chromosome. Column 1 is the rsID and column 2 is the VCF variant ID.
- For scoring, a tab-separated manifest with `trait`, `weights`, `id_col`,
  `allele_col`, and `effect_col` columns. See `examples/scores.tsv`.

Input patterns must contain the literal `{chr}` placeholder. Local paths and S3
URIs are supported. The existing rsID maps are available from
[Google Drive](https://drive.google.com/drive/folders/199d80SdlSaYum8GNCQkFnh8zvhArmxnj?usp=sharing).

## Build and run locally

Requirements: Java 17+, Nextflow 26.04+, and a running Docker daemon.

```bash
nextflow run main.nf -profile local \
  --vcfs 'vcfs/chr{chr}.dose.vcf.gz' \
  --rsid_maps 'rsid_maps/chr{chr}.map' \
  --outdir results
```

Use `--chromosomes '1,2'` for a small test and add `-resume` to resume an
interrupted run. Reports are written locally under `pipeline_info` (override
with `--report_dir`).

The local profile pulls the pinned PLINK2 Biocontainer directly. Enable scoring
with `--run_scores true --score_sheet scores.tsv`. Enable fixed-reference PCA
with `--run_pca true --pca_reference_sheet pca_reference.tsv`. The PCA contract
is described in `docs/pca_ancestry_design.md`; the previous cohort-derived PCA
implementation was removed because it did not implement 1000 Genomes projection.

## Run on AWS Batch

AWS uses the same PLINK2 Biocontainer. The Batch host AMI must expose an AWS CLI
launcher at `/opt/nxf-aws-cli/bin/aws`; Nextflow mounts it into each task for S3
staging. Mirror
the public image to ECR if required by site policy, then launch with an IAM role:

```bash
AWS_ACCOUNT_ID=123456789012
AWS_REGION=us-west-2
ECR_REPOSITORY=pgs-pipeline

aws ecr get-login-password --region "$AWS_REGION" | \
  docker login --username AWS --password-stdin \
  "$AWS_ACCOUNT_ID.dkr.ecr.$AWS_REGION.amazonaws.com"
docker pull quay.io/biocontainers/plink2:2.0.0a.6.9--h9948957_0
docker tag quay.io/biocontainers/plink2:2.0.0a.6.9--h9948957_0 \
  "$AWS_ACCOUNT_ID.dkr.ecr.$AWS_REGION.amazonaws.com/$ECR_REPOSITORY:latest"
docker push \
  "$AWS_ACCOUNT_ID.dkr.ecr.$AWS_REGION.amazonaws.com/$ECR_REPOSITORY:latest"

nextflow run main.nf -profile awsbatch \
  --container "$AWS_ACCOUNT_ID.dkr.ecr.$AWS_REGION.amazonaws.com/$ECR_REPOSITORY:latest" \
  --aws_queue 'my-batch-queue' \
  --aws_region "$AWS_REGION" \
  --aws_work_dir 's3://my-pgs-bucket/nextflow-work' \
  --outdir 's3://my-pgs-bucket/results/run-001' \
  --vcfs 's3://my-pgs-bucket/vcfs/chr{chr}.dose.vcf.gz' \
  --rsid_maps 's3://my-pgs-bucket/rsid_maps/chr{chr}.map'
```

For data already resident on an FSx for Lustre mount, use the FSx profile and
pass absolute FSx paths. This avoids copying the VCF through the S3 work area:

```bash
nextflow run main.nf -profile awsbatch_fsx \
  --fsx_direct true \
  --fsx_mount /fsx \
  --vcfs '/fsx/path/to/chr{chr}.vcf.gz' \
  --rsid_maps '/fsx/path/to/chr{chr}.map' \
  --aws_queue 'my-fsx-batch-queue' \
  --aws_region us-east-1 \
  --aws_work_dir 's3://my-bucket/nextflow-work' \
  --outdir 's3://my-bucket/results/run-001'
```

The Batch compute environment launch template must mount FSx at `--fsx_mount`.
For an infrastructure-only smoke test, `--skip_rsid_annotation true` permits
conversion and QC without an rsID map. Do not use that option for PGS scoring,
which requires variants to be aligned to the weight files by rsID.

On Slurm with a shared filesystem and Apptainer:

```bash
nextflow run main.nf -profile slurm \
  --container /shared/containers/pgs-pipeline.sif \
  --shared_work_dir /shared/nextflow-work \
  --direct_inputs true \
  --vcfs '/shared/vcfs/chr{chr}.dose.vcf.gz' \
  --rsid_maps '/shared/rsid_maps/chr{chr}.map'
```

To run scoring and/or PCA from an existing QCed PLINK 2 dataset, provide its
prefix (the path before `.pgen`, `.pvar`, and `.psam`). This bypasses VCF
conversion, chromosome concatenation, and missingness filtering:

```bash
nextflow run main.nf -profile slurm \
  --input_pfile /shared/cohorts/cohort \
  --run_pca true \
  --pca_reference_sheet examples/pca_reference.tsv \
  --shared_work_dir /shared/nextflow-work \
  --outdir /shared/results/run-001
```

The input is treated as already QCed; use the VCF path when the pipeline should
apply its own preprocessing and missingness thresholds.

Add `--direct_inputs true` (or `--fsx_direct true` with the AWS FSx profile)
when that PGEN prefix is on a filesystem mounted at the same absolute path on
the controller and every worker. The workflow then reads the PGEN in place
instead of staging it through the Nextflow work directory.

`--direct_inputs` is executor-neutral: use it whenever controller and workers
share the same absolute POSIX paths. Without it, Nextflow stages declared input
files normally. Other schedulers can be added in a site-specific config passed
with `-c`, without changing `main.nf`. An Expanse adapter is provided in
`conf/expanse.config`; its account and queue are site policy, not workflow
requirements.

The Batch compute environment must be able to pull the ECR image. Its job role
must be able to read inputs and read/write the work and output S3 prefixes. Do
not put access keys in this repository; use IAM roles or the AWS credential
chain. The AWS profile follows Nextflow's documented model of a Batch queue, S3
work directory, and AWS CLI in the task image ([Nextflow documentation](https://training.nextflow.io/2.5.0/archive/basic_training/executors/)).

## Important parameters

| Parameter | Default | Meaning |
|---|---:|---|
| `chromosomes` | `1..22` | Inclusive range or comma-separated list |
| `input_pfile` | unset | Existing QCed PLINK 2 prefix; bypasses VCF preprocessing |
| `cohort` | `cohort` | Output prefix |
| `genome_build` | `GRCh38` | Cohort genome build; must match the PCA reference |
| `mac` | `10` | Minimum allele count |
| `geno` | `0.05` | Initial genotype missingness threshold |
| `maf` | `0.01` | Initial minor allele frequency threshold |
| `variant_miss` | `0.05` | Variant missingness threshold |
| `sample_miss` | `0.05` | Sample missingness threshold |
| `run_scores` | `false` | Run PLINK2 scoring branch |
| `run_summary_qc` | `true` | Generate cohort-wide missingness, HWE, and frequency summaries; set false for repeated scoring-only runs |
| `min_score_variant_match` | `0.50` | Fail a trait when fewer than this fraction of weight variants are scored |
| `warn_score_variant_match` | `0.80` | Flag a trait QC row below this match fraction |
| `run_pca` | `false` | Run reference validation, harmonization, global projection, and ancestry assignment |
| `pca_reference_sheet` | unset | One-row reference artifact manifest; see `examples/pca_reference.tsv` |
| `num_pcs` | `10` | Number of global reference PCs to project |
| `num_within_ancestry_pcs` | `10` | Maximum PCs calculated within each ancestry group |
| `min_pca_variant_overlap` | `0.90` | Minimum usable fraction of the fixed PCA panel |
| `min_ancestry_samples` | `50` | Minimum group size for within-ancestry PCA |
| `within_ancestry_king_cutoff` | `0.0884` | KING relatedness cutoff for the unrelated PCA training subset |
| `within_ancestry_ld_window` | `200` | LD-pruning window size in variants |
| `within_ancestry_ld_step` | `50` | LD-pruning step size in variants |
| `within_ancestry_ld_r2` | `0.2` | LD-pruning r-squared threshold |
| `r2` / `aq` | unset | Optional VCF INFO filter; R2 takes precedence |

Resource defaults live in `nextflow.config` and can be overridden with `-c`.
Scoring uses a dedicated portable default of 4 CPUs and 8 GB RAM; site adapters
can override the `scoring` process label without changing the workflow.
Each scoring run publishes per-trait QC, `combined_scores.tsv`, and
`score_qc_summary.tsv` alongside the PLINK score files.
`config.toml.example` is only for the legacy Python entrypoints.

## Weight-generation boundary

Weight generation is intentionally upstream of this workflow. SBayesRC,
PRS-CS, published PGS Catalog weights, or another method can run elsewhere. The
score manifest tells PLINK2 which one-based columns contain variant ID, effect
allele, and effect size, so the scoring branch is not tied to GCTB `.snpRes`.

The pipeline does not currently normalize allele conventions across arbitrary
weight formats; weight files must already be harmonized to the cohort build and
variant identifiers. The legacy `sbayesrc/` scripts remain as reference only.

The legacy `01d` multi-batch merge is also not represented in `main.nf`; inputs
are currently one VCF per chromosome.

Dataset-specific preparation, including G2MH WGS/GSA harmonization, is upstream
of this workflow. A pipeline run receives one coherent cohort dataset.
