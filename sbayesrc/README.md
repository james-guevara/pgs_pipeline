# SBayesRC Pipeline

Computes polygenic score weights from GWAS summary statistics using SBayesRC.

## Files

- `run_sbayesrc.sh` - Main SLURM array job script (refactored)
- `gwas_list.txt` - List of .ma files to process (one per line)
- `sbayesrc_pipeline.sh` - Original script (reference only)

## Usage on Expanse

```bash
# 1. Copy files to Expanse
scp run_sbayesrc.sh gwas_list.txt expanse:/expanse/projects/sebat1/s3/data/sebat/pgs_pipeline/

# 2. Copy COJO files from sumstats/cojo/
scp ../sumstats/cojo/*.ma expanse:/expanse/projects/sebat1/s3/data/sebat/pgs_pipeline/input/

# 3. Update array size to match trait count
#    Edit --array=1-N where N = wc -l gwas_list.txt

# 4. Submit
cd /expanse/projects/sebat1/s3/data/sebat/pgs_pipeline
sbatch run_sbayesrc.sh
```

## Pipeline Steps

1. **Impute summary statistics** - Fill in missing SNPs using LD reference
2. **SBayesRC** - Compute SNP effect weights with 97 functional annotations

## Input/Output

- Input: `input/*.ma` (COJO format from sumstats wrangling)
- Intermediate: `output/tidy/*.imputed.ma`
- Output: `output/weights/*.snpRes` (SNP weights for PRS calculation)

## Known Issues

**Output format**: Uses `gctb` native output (`.snpRes`) which has more columns
than the `sbayesrc` R wrapper (`.txt`). Key column for PGS is `A1Effect`.

**PTSD (ptsd_2024.ma)**: Segfaults during imputation at LD block 300 (chr8)
regardless of memory (48G, 96G, 128G). Previous PTSD runs used different
input files (eur_ptsdcasecontrol_pcs_v4_aug3_2021) which completed successfully.
The 2024 sumstats may have data issues in the chr8 region. Options:
- Use existing weights from old PTSD run
- Try running without imputation step
- Debug chr8 SNPs in input file

## Resources (on Expanse)

```
/expanse/projects/sebat1/s3/data/sebat/resources/
├── sbayesrc_resources/
│   ├── ukbEUR_HM3/          # LD reference (sparse eigen decomposition)
│   └── annot_baseline2.2.txt # 97 functional annotations
└── sbayesrc_v0.2.6.sif       # Singularity container (GCTB 2.05beta)
```
