#!/usr/bin/env python3
"""Step 03: Missingness QC on the genome-wide concatenated pfile.

Two-pass approach on a single file (no re-concatenation):
  1. Compute missingness → identify high-missingness variants
  2. Exclude bad variants, recompute → identify high-missingness samples
  3. Apply both filters → filtered output
"""
import subprocess
import argparse
from pathlib import Path

from config import load_config, get


def run(cmd):
    """Run and print shell commands."""
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)


def count_lines(filepath):
    """Count non-empty lines in a file."""
    with open(filepath) as f:
        return sum(1 for line in f if line.strip())


def main():
    cfg = load_config()

    concat_dir = Path(get(cfg, "directories", "concat"))
    default_input = str(concat_dir / "cohort")

    parser = argparse.ArgumentParser(description="Step 03: Missingness QC on genome-wide pfile.")
    parser.add_argument("--config", type=Path, help="Path to config.toml")
    parser.add_argument("--pfile", type=str, default=default_input,
                        help="Input pfile prefix (e.g., 02_concat/cohort)")
    parser.add_argument("--out-dir", type=Path, default=Path(get(cfg, "directories", "missingness")),
                        help="Output directory")
    parser.add_argument("--cohort", default="cohort", help="Output prefix name")
    parser.add_argument("--threads", type=int, default=get(cfg, "defaults", "threads"))
    parser.add_argument("--memory", type=int, default=get(cfg, "defaults", "memory"))
    parser.add_argument("--sample-miss", type=float, default=get(cfg, "filters", "sample_miss"),
                        help="Sample missingness threshold")
    parser.add_argument("--variant-miss", type=float, default=get(cfg, "filters", "variant_miss"),
                        help="Variant missingness threshold")
    args = parser.parse_args()

    if args.config:
        cfg = load_config(args.config)

    out_dir = args.out_dir
    out_dir.mkdir(exist_ok=True)

    pfile = args.pfile
    fail_variants = out_dir / "fail_variants.txt"
    fail_samples = out_dir / "fail_samples.txt"

    # --- Pass 1: compute missingness, identify bad variants ---
    print("Pass 1: Computing missingness...")
    run([
        "plink2",
        "--pfile", pfile,
        "--missing",
        "--out", str(out_dir / "pass1"),
        "--threads", str(args.threads),
        "--memory", str(args.memory)
    ])

    vmiss = out_dir / "pass1.vmiss"
    smiss = out_dir / "pass1.smiss"

    # Count totals from plink2 output
    n_variants_total = sum(1 for i, line in enumerate(vmiss.open()) if i > 0 and line.strip())
    n_samples_total = sum(1 for i, line in enumerate(smiss.open()) if i > 0 and line.strip())

    # Generate variant fail list
    subprocess.run(
        f"awk 'NR>1 && $5>{args.variant_miss} {{print $2}}' {vmiss} > {fail_variants}",
        shell=True, check=True)
    n_fail_variants = count_lines(fail_variants)
    print(f"  Variants: {n_variants_total} total, {n_fail_variants} fail (>{args.variant_miss}), "
          f"{n_variants_total - n_fail_variants} remaining")

    # --- Pass 2: exclude bad variants, compute sample missingness ---
    print("\nPass 2: Recomputing sample missingness after variant filter...")
    run([
        "plink2",
        "--pfile", pfile,
        "--exclude", str(fail_variants),
        "--missing",
        "--out", str(out_dir / "pass2"),
        "--threads", str(args.threads),
        "--memory", str(args.memory)
    ])

    smiss2 = out_dir / "pass2.smiss"
    subprocess.run(
        f"awk 'NR>1 && $4>{args.sample_miss} {{print $1}}' {smiss2} > {fail_samples}",
        shell=True, check=True)
    n_fail_samples = count_lines(fail_samples)
    print(f"  Samples: {n_samples_total} total, {n_fail_samples} fail (>{args.sample_miss}), "
          f"{n_samples_total - n_fail_samples} remaining")

    # --- Apply both filters ---
    print("\nApplying filters...")
    out_prefix = str(out_dir / args.cohort)
    run([
        "plink2",
        "--pfile", pfile,
        "--exclude", str(fail_variants),
        "--remove", str(fail_samples),
        "--make-pgen",
        "--out", out_prefix,
        "--threads", str(args.threads),
        "--memory", str(args.memory)
    ])

    # Clean up intermediate files
    for name in ("pass1.smiss", "pass1.vmiss", "pass1.log",
                 "pass2.smiss", "pass2.vmiss", "pass2.log"):
        p = out_dir / name
        if p.exists():
            p.unlink()

    print(f"\nMissingness QC complete.")
    print(f"  Variants: {n_variants_total} -> {n_variants_total - n_fail_variants} "
          f"({n_fail_variants} removed)")
    print(f"  Samples:  {n_samples_total} -> {n_samples_total - n_fail_samples} "
          f"({n_fail_samples} removed)")
    print(f"Output: {out_prefix}.pgen/pvar/psam")


if __name__ == "__main__":
    main()
