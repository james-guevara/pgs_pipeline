#!/usr/bin/env python3
import subprocess
from pathlib import Path
import argparse

from config import load_config, get


def run(cmd):
    """Run and print shell commands."""
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)


# -------------------------------------------------------------------------
# Individual steps
# -------------------------------------------------------------------------
def maf_filter(merged, out_dir, cohort, maf, threads, memory):
    """Apply MAF >= threshold filter."""
    out_prefix = Path(out_dir) / cohort
    print(f"\n[1/4] Applying MAF >= {maf} filter on merged dataset...")
    run([
        "plink2",
        "--pfile", str(merged),
        "--maf", str(maf),
        "--make-pgen",
        "--out", str(out_prefix),
        "--threads", str(threads),
        "--memory", str(memory),
    ])
    return out_prefix


def ld_prune(maf_prefix, out_dir, cohort, ld_window, ld_step, ld_r2, threads, memory):
    """Perform LD pruning and extract pruned SNPs."""
    prune_base = Path(out_dir) / cohort
    pruned_prefix = Path(out_dir) / f"{cohort}_pruned"

    print("\n[2/4] LD pruning on MAF-filtered dataset...")
    run([
        "plink2",
        "--pfile", str(maf_prefix),
        "--indep-pairwise", str(ld_window), str(ld_step), str(ld_r2),
        "--out", str(prune_base),
        "--threads", str(threads),
        "--memory", str(memory),
    ])

    run([
        "plink2",
        "--pfile", str(maf_prefix),
        "--extract", f"{prune_base}.prune.in",
        "--make-pgen",
        "--out", str(pruned_prefix),
        "--threads", str(threads),
        "--memory", str(memory),
    ])
    return pruned_prefix


def king_relatedness(pruned_prefix, out_dir, cohort, king_cutoff, threads, memory):
    """Compute KING relatedness and define unrelateds."""
    king_prefix = Path(out_dir) / cohort
    unrelated_prefix = Path(out_dir) / "unrelateds"

    print("\n[3/4] Running KING relatedness on pruned SNP set...")
    run([
        "plink2",
        "--pfile", str(pruned_prefix),
        "--make-king-table",
        "--out", str(king_prefix),
        "--threads", str(threads),
        "--memory", str(memory),
    ])

    run([
        "plink2",
        "--pfile", str(pruned_prefix),
        "--king-cutoff", str(king_cutoff),
        "--make-just-fam",
        "--out", str(unrelated_prefix),
        "--threads", str(threads),
        "--memory", str(memory),
    ])

    keep_file = unrelated_prefix.with_suffix(".king.cutoff.in.id")
    if not keep_file.exists():
        raise FileNotFoundError(f"Unrelated ID file not found: {keep_file}")
    return king_prefix, keep_file


def pca_analysis(pruned_prefix, keep_file, out_dir, cohort, threads, memory, num_pcs):
    """Run PCA on unrelateds and project all samples."""
    ref_prefix = Path(out_dir) / "ref"
    proj_prefix = Path(out_dir) / cohort

    print("\n[4/4] PCA on unrelateds, then projecting all samples...")
    run([
        "plink2",
        "--pfile", str(pruned_prefix),
        "--keep", str(keep_file),
        "--freq", "counts",
        "--pca", "allele-wts", "vcols=chrom,ref,alt",
        "--out", str(ref_prefix),
        "--threads", str(threads),
        "--memory", str(memory),
    ])

    first_pc_col = 6
    last_pc_col = 5 + num_pcs
    score_cols = f"{first_pc_col}-{last_pc_col}"

    run([
        "plink2",
        "--pfile", str(pruned_prefix),
        "--read-freq", f"{ref_prefix}.acount",
        "--score", f"{ref_prefix}.eigenvec.allele",
        "2", "5",
        "header-read",
        "no-mean-imputation",
        "variance-standardize",
        "--score-col-nums", score_cols,
        "--out", str(proj_prefix),
        "--threads", str(threads),
        "--memory", str(memory),
    ])

    return ref_prefix, proj_prefix


# -------------------------------------------------------------------------
# Main orchestration
# -------------------------------------------------------------------------
def ancestry_pipeline(cfg, merged_prefix, cohort, maf, king_cutoff, num_pcs,
                      ld_window, ld_step, ld_r2, threads, memory):
    """Full ancestry analysis pipeline from merged pfile."""
    merged = Path(merged_prefix)

    # Get directory names from config
    maf_dir = Path(get(cfg, "directories", "ancestry_maf"))
    prune_dir = Path(get(cfg, "directories", "ancestry_prune"))
    rel_dir = Path(get(cfg, "directories", "ancestry_relatedness"))
    pca_dir = Path(get(cfg, "directories", "ancestry_pca"))
    for d in (maf_dir, prune_dir, rel_dir, pca_dir):
        d.mkdir(exist_ok=True)

    maf_prefix = maf_filter(merged, maf_dir, cohort, maf, threads, memory)
    pruned_prefix = ld_prune(maf_prefix, prune_dir, cohort, ld_window, ld_step, ld_r2, threads, memory)
    king_prefix, keep_file = king_relatedness(pruned_prefix, rel_dir, cohort, king_cutoff, threads, memory)
    ref_prefix, proj_prefix = pca_analysis(pruned_prefix, keep_file, pca_dir, cohort, threads, memory, num_pcs)

    print("\nDone.")
    print(f"MAF-filtered dataset:      {maf_prefix}.*")
    print(f"Pruned dataset:            {pruned_prefix}.*")
    print(f"KING results:              {king_prefix}.king(.id/.kin0)")
    print(f"Unrelated IDs:             {keep_file}")
    print(f"PCA loadings (ref PCs):    {ref_prefix}.eigenvec.allele, {ref_prefix}.eigenval")
    print(f"PCA projection (all inds): {proj_prefix}.sscore (PCs in SCORE columns)")


def main():
    cfg = load_config()

    parser = argparse.ArgumentParser(description="Run ancestry pipeline (MAF -> prune -> KING -> PCA).")
    parser.add_argument("--config", type=Path, help="Path to config.toml")
    parser.add_argument("--merged-prefix", required=True, help="Prefix of merged pfile (e.g., 03_summary/cohort).")
    parser.add_argument("--cohort", default="cohort", help="Cohort name")
    parser.add_argument("--maf", type=float, default=get(cfg, "ancestry", "maf"), help="MAF threshold")
    parser.add_argument("--king-cutoff", type=float, default=get(cfg, "ancestry", "king_cutoff"),
                        help="KING cutoff")
    parser.add_argument("--num-pcs", type=int, default=get(cfg, "ancestry", "num_pcs"), help="Number of PCs")
    parser.add_argument("--ld-window", type=int, default=get(cfg, "ancestry", "ld_window"), help="LD window size")
    parser.add_argument("--ld-step", type=int, default=get(cfg, "ancestry", "ld_step"), help="LD step size")
    parser.add_argument("--ld-r2", type=float, default=get(cfg, "ancestry", "ld_r2"), help="LD r2 threshold")
    parser.add_argument("--threads", type=int, default=get(cfg, "defaults", "threads"), help="Threads")
    parser.add_argument("--memory", type=int, default=get(cfg, "defaults", "memory"), help="Memory (MB)")
    args = parser.parse_args()

    # Reload config if custom path provided
    if args.config:
        cfg = load_config(args.config)

    ancestry_pipeline(
        cfg=cfg,
        merged_prefix=args.merged_prefix,
        cohort=args.cohort,
        maf=args.maf,
        king_cutoff=args.king_cutoff,
        num_pcs=args.num_pcs,
        ld_window=args.ld_window,
        ld_step=args.ld_step,
        ld_r2=args.ld_r2,
        threads=args.threads,
        memory=args.memory,
    )


if __name__ == "__main__":
    main()
