# PGS pipeline

Containerized Nextflow pipeline for VCF filtering, rsID annotation, genome-wide
QC, polygenic scoring, relatedness, and PCA. It can run unchanged with local
Docker or AWS Batch.

## Workflow

1. Run PLINK2 filtering, rsID annotation, and MAF filtering independently for each chromosome.
2. Concatenate chromosome pfiles into one cohort pfile.
3. Apply two-pass variant and sample missingness QC.
4. Produce missingness, Hardy-Weinberg, and allele-frequency summaries.
5. Optionally calculate PGS values from externally generated weight files.
6. Run ancestry/PCA; this required final stage is temporarily gated while its
   statistical design is reconsidered.

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
with `--run_scores true --score_sheet scores.tsv`. PCA is temporarily off by
default while it is redesigned; `--run_ancestry true` enables the current
implementation for comparison testing.

## Run on AWS Batch

`Dockerfile.aws` is a thin AWS derivative of the same PLINK2 Biocontainer. It
adds AWS CLI for Nextflow's S3 staging and no scientific software. Build it in
AWS and push it to ECR, then launch Nextflow with an IAM role:

```bash
AWS_ACCOUNT_ID=123456789012
AWS_REGION=us-west-2
ECR_REPOSITORY=pgs-pipeline

aws ecr get-login-password --region "$AWS_REGION" | \
  docker login --username AWS --password-stdin \
  "$AWS_ACCOUNT_ID.dkr.ecr.$AWS_REGION.amazonaws.com"
docker build -f Dockerfile.aws -t pgs-pipeline:latest .
docker tag pgs-pipeline:latest \
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

`--direct_inputs` is executor-neutral: use it whenever controller and workers
share the same absolute POSIX paths. Without it, Nextflow stages declared input
files normally. Other schedulers can be added in a site-specific config passed
with `-c`, without changing `main.nf`.

The Batch compute environment must be able to pull the ECR image. Its job role
must be able to read inputs and read/write the work and output S3 prefixes. Do
not put access keys in this repository; use IAM roles or the AWS credential
chain. The AWS profile follows Nextflow's documented model of a Batch queue, S3
work directory, and AWS CLI in the task image ([Nextflow documentation](https://training.nextflow.io/2.5.0/archive/basic_training/executors/)).

## Important parameters

| Parameter | Default | Meaning |
|---|---:|---|
| `chromosomes` | `1..22` | Inclusive range or comma-separated list |
| `cohort` | `cohort` | Output prefix |
| `mac` | `10` | Minimum allele count |
| `geno` | `0.05` | Initial genotype missingness threshold |
| `maf` | `0.01` | Initial minor allele frequency threshold |
| `variant_miss` | `0.05` | Variant missingness threshold |
| `sample_miss` | `0.05` | Sample missingness threshold |
| `run_scores` | `false` | Run PLINK2 scoring branch |
| `run_ancestry` | `false` | Temporary gate for the required, pending-redesign KING/PCA stage |
| `r2` / `aq` | unset | Optional VCF INFO filter; R2 takes precedence |

Resource defaults live in `nextflow.config` and can be overridden with `-c`.
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
