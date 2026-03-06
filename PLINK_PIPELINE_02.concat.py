#!/usr/bin/env python3
"""Step 02: Concatenate per-chromosome pfiles into a single genome-wide pfile.

Reads per-chromosome pgen/pvar/psam files (same samples, different chromosomes)
and produces one merged pfile using plink2 --pmerge-list.
"""
import subprocess
import argparse
from pathlib import Path

from config import load_config, get


def run(cmd):
    """Run and print shell commands."""
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)


def main():
    cfg = load_config()

    # Default input: prefer merge dir if it has pgen files, else maf dir
    merge_dir = Path(get(cfg, "directories", "merge"))
    maf_dir = Path(get(cfg, "directories", "maf"))
    default_dir = merge_dir if (merge_dir / "chr1.pgen").exists() else maf_dir

    parser = argparse.ArgumentParser(description="Step 02: Concatenate per-chr pfiles into genome-wide pfile.")
    parser.add_argument("--config", type=Path, help="Path to config.toml")
    parser.add_argument("--dir", type=Path, default=default_dir,
                        help="Directory containing per-chromosome pgen files")
    parser.add_argument("--out-dir", type=Path, default=Path(get(cfg, "directories", "concat")),
                        help="Output directory")
    parser.add_argument("--cohort", default="cohort", help="Output prefix name")
    parser.add_argument("--threads", type=int, default=get(cfg, "defaults", "threads"))
    parser.add_argument("--memory", type=int, default=get(cfg, "defaults", "memory"))
    args = parser.parse_args()

    if args.config:
        cfg = load_config(args.config)

    input_dir = args.dir
    out_dir = args.out_dir
    out_dir.mkdir(exist_ok=True)

    # Build merge list
    merge_list = out_dir / "merge_list.txt"
    missing_chr = []
    with merge_list.open("w") as f:
        for chr_num in range(1, 23):
            prefix = input_dir / f"chr{chr_num}"
            if prefix.with_suffix(".pgen").exists():
                f.write(str(prefix) + "\n")
            else:
                missing_chr.append(chr_num)

    if missing_chr:
        print(f"Warning: Skipping missing chromosomes: {' '.join(map(str, missing_chr))}")

    n_chr = 22 - len(missing_chr)
    print(f"Concatenating {n_chr} chromosomes from {input_dir}...")

    out_prefix = str(out_dir / args.cohort)
    run([
        "plink2",
        "--pmerge-list", str(merge_list),
        "--make-pgen",
        "--out", out_prefix,
        "--threads", str(args.threads),
        "--memory", str(args.memory)
    ])

    # Clean up merge list
    if merge_list.exists():
        merge_list.unlink()

    # Report counts
    psam = Path(f"{out_prefix}.psam")
    pvar = Path(f"{out_prefix}.pvar")
    n_samples = sum(1 for line in psam.open() if not line.startswith("#")) - 1
    n_variants = sum(1 for line in pvar.open() if not line.startswith("#"))

    print(f"\nConcat complete: {n_samples} samples, {n_variants} variants")
    print(f"Output: {out_prefix}.pgen/pvar/psam")


if __name__ == "__main__":
    main()
