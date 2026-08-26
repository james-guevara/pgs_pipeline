# PGS pipeline

Containerized Nextflow pipeline for VCF filtering, rsID annotation, genome-wide
QC, polygenic scoring, relatedness, and PCA. It can run unchanged with local
Docker or AWS Batch.

## Workflow

1. Run PLINK2 filtering, rsID annotation, and MAF filtering independently for each chromosome.
2. Concatenate chromosome pfiles into one cohort pfile.
3. Apply two-pass variant and sample missingness QC.
4. Produce missingness, Hardy-Weinberg, and allele-frequency summaries.
5. Optionally calculate PGS values from SBayesRC `.snpRes` weights.
6. Optionally run LD pruning, KING, PCA on unrelated samples, and projection of all samples.

The original Python and Slurm entrypoints remain for reference. `main.nf` is the
portable orchestration entrypoint.

## Required inputs

- One bgzipped VCF per chromosome.
- One two-column rsID map per chromosome. Column 1 is the rsID and column 2 is the VCF variant ID.
- For scoring, a tab-separated `trait<TAB>weight_file` sheet containing SBayesRC `.snpRes` files. See `examples/scores.tsv`.

Input patterns must contain the literal `{chr}` placeholder. Local paths and S3
URIs are supported. The existing rsID maps are available from
[Google Drive](https://drive.google.com/drive/folders/199d80SdlSaYum8GNCQkFnh8zvhArmxnj?usp=sharing).

## Build and run locally

Requirements: Java 17+, Nextflow 24.10+, and a running Docker daemon.

```bash
docker build -t pgs-pipeline:latest .

nextflow run main.nf -profile local \
  --vcfs 'vcfs/chr{chr}.dose.vcf.gz' \
  --rsid_maps 'rsid_maps/chr{chr}.map' \
  --outdir results
```

Use `--chromosomes '1,2'` for a small test and add `-resume` to resume an
interrupted run. Reports are written locally under `pipeline_info` (override
with `--report_dir`).

Enable scoring with `--run_scores true --score_sheet scores.tsv`. Disable the
ancestry branch with `--run_ancestry false`.

## Run on AWS Batch

Build the image, push it to ECR, then launch Nextflow from a machine with AWS
credentials or an IAM role:

```bash
AWS_ACCOUNT_ID=123456789012
AWS_REGION=us-west-2
ECR_REPOSITORY=pgs-pipeline

aws ecr get-login-password --region "$AWS_REGION" | \
  docker login --username AWS --password-stdin \
  "$AWS_ACCOUNT_ID.dkr.ecr.$AWS_REGION.amazonaws.com"
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
| `run_ancestry` | `true` | Run KING/PCA branch |
| `r2` / `aq` | unset | Optional VCF INFO filter; R2 takes precedence |

Resource defaults live in `nextflow.config` and can be overridden with `-c`.
`config.toml.example` is only for the legacy Python entrypoints.

## SBayesRC boundary and open decisions

The Nextflow workflow consumes completed SBayesRC `.snpRes` weights; it does not
yet generate them. The existing `sbayesrc/` job requires a GCTB image, UKB EUR
HM3 LD reference, and baseline annotation file that are not in this repository.

Settle these items before adding SBayesRC as an AWS process:

- Approved GCTB/SBayesRC image and redistribution terms.
- S3 locations for the LD matrix and annotations, and whether tasks stage them or use a shared filesystem.
- Whether the problematic 2024 PTSD input skips imputation, uses older weights, or is repaired around chromosome 8.
- Batch queue sizing and whether Spot instances are acceptable.
- Whether cohorts have multiple imputation batches. The legacy `01d` batch merge is not in `main.nf`; inputs are currently one VCF per chromosome.
