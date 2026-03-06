#!/usr/bin/env python3
"""Step 04: Summary statistics on the missingness-filtered genome-wide pfile.

Computes missingness, Hardy-Weinberg, and allele frequency stats.
No concatenation needed — operates on the single filtered pfile from step 03.
"""
import subprocess
from pathlib import Path
import argparse

from config import load_config, get


def run(cmd):
    """Run and print shell commands."""
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)


def main():
    cfg = load_config()

    miss_dir = Path(get(cfg, "directories", "missingness"))
    default_input = str(miss_dir / "cohort")

    parser = argparse.ArgumentParser(description="Step 04: Summary QC stats.")
    parser.add_argument("--config", type=Path, help="Path to config.toml")
    parser.add_argument("--pfile", type=str, default=default_input,
                        help="Input pfile prefix (e.g., 03_missingness/cohort)")
    parser.add_argument("--cohort", default="cohort", help="Cohort name for output prefix")
    parser.add_argument("--out-dir", type=Path, default=Path(get(cfg, "directories", "summary")),
                        help="Output directory")
    parser.add_argument("--threads", type=int, default=get(cfg, "defaults", "threads"))
    parser.add_argument("--memory", type=int, default=get(cfg, "defaults", "memory"))
    args = parser.parse_args()

    if args.config:
        cfg = load_config(args.config)

    out_dir = args.out_dir
    out_dir.mkdir(exist_ok=True)

    out_prefix = str(out_dir / args.cohort)

    print("Computing summary statistics...")
    run([
        "plink2",
        "--pfile", args.pfile,
        "--missing",
        "--hardy",
        "--freq",
        "--out", out_prefix,
        "--threads", str(args.threads),
        "--memory", str(args.memory)
    ])

    # Write summary counts
    smiss = Path(f"{out_prefix}.smiss")
    vmiss = Path(f"{out_prefix}.vmiss")

    nsamp = max(sum(1 for i, line in enumerate(smiss.open()) if i > 0 and line.strip()), 0)
    nvar = max(sum(1 for i, line in enumerate(vmiss.open()) if i > 0 and line.strip()), 0)

    summary_file = out_dir / "summary_counts.txt"
    with summary_file.open("w") as f:
        f.write("Metric\tCount\n")
        f.write(f"Samples\t{nsamp}\n")
        f.write(f"Variants\t{nvar}\n")

    print(f"\nSummary QC complete: {nsamp} samples, {nvar} variants")
    print(f"Results in {out_dir}/")


if __name__ == "__main__":
    main()
