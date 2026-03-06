#!/usr/bin/env python3
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

    parser = argparse.ArgumentParser(description="Run PLINK2 missingness QC + filtering.")
    parser.add_argument("--config", type=Path, help="Path to config.toml")
    # Default input: prefer merge dir if it exists, else maf dir
    merge_dir = Path(get(cfg, "directories", "merge"))
    maf_dir = Path(get(cfg, "directories", "maf"))
    default_dir = merge_dir if (merge_dir / "chr1.pgen").exists() or (merge_dir / "chr1.bed").exists() else maf_dir
    parser.add_argument("--dir", type=Path, default=default_dir,
                        help="Directory containing per-chromosome .pgen/.psam/.pvar files")
    parser.add_argument("--out-dir", type=Path, default=Path(get(cfg, "directories", "missingness")),
                        help="Output directory for missingness results")
    parser.add_argument("--memory", type=int, default=get(cfg, "defaults", "memory"))
    parser.add_argument("--threads", type=int, default=get(cfg, "defaults", "threads"))
    parser.add_argument("--sample-miss", type=float, default=get(cfg, "filters", "sample_miss"),
                        help="Sample missingness threshold")
    parser.add_argument("--variant-miss", type=float, default=get(cfg, "filters", "variant_miss"),
                        help="Variant missingness threshold")
    args = parser.parse_args()

    input_dir = args.dir
    out_dir = args.out_dir
    memory = args.memory
    threads = args.threads

    out_dir.mkdir(exist_ok=True)
    merge_list = out_dir / "merge_list.txt"

    # Detect input format (pgen or bed)
    use_bed = (input_dir / "chr1.bed").exists() and not (input_dir / "chr1.pgen").exists()

    print(f"Building merge list from {input_dir}...")
    print(f"Detected format: {'bed' if use_bed else 'pgen'}")
    missing_chr = []
    with merge_list.open("w") as f:
        for chr_num in range(1, 23):
            prefix = input_dir / f"chr{chr_num}"
            ext = ".bed" if use_bed else ".pgen"
            if (prefix.with_suffix(ext)).exists():
                f.write(str(prefix) + "\n")
            else:
                missing_chr.append(chr_num)

    if missing_chr:
        print(f"Warning: Skipping missing chromosomes: {' '.join(map(str, missing_chr))}")
    print(f"Merge list written to {merge_list}\n")

    # Two-pass approach:
    #   Pass 1: compute missingness, filter variants first (removes batch-specific variants)
    #   Pass 2: recompute missingness on remaining variants, then filter samples

    merge_list_mode = "bfile" if use_bed else "pfile"

    # --- Pass 1: variant missingness ---
    print("Pass 1: Computing missingness to identify high-missingness variants...")
    run([
        "plink2",
        "--pmerge-list", str(merge_list), merge_list_mode,
        "--missing",
        "--out", str(out_dir / "pass1"),
        "--threads", str(threads),
        "--memory", str(memory)
    ])

    vmiss = out_dir / "pass1.vmiss"
    fail_variants = out_dir / "fail_variants.txt"
    subprocess.run(f"awk 'NR>1 && $5>{args.variant_miss} {{print $2}}' {vmiss} > {fail_variants}", shell=True, check=True)

    nvar_total_proc = subprocess.run(f"wc -l < {vmiss}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    nvar_total = max(int(nvar_total_proc.stdout.decode().strip() or 1) - 1, 0)

    nvar_proc = subprocess.run(f"wc -l < {fail_variants}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    nvar = int(nvar_proc.stdout.decode().strip() or 0)
    nvar_keep = nvar_total - nvar
    print(f"Pass 1: {nvar_total} total variants, {nvar} fail missingness > {args.variant_miss}, {nvar_keep} remaining")

    # Apply variant filter per chromosome into temp dir
    pass1_dir = out_dir / "pass1_filtered"
    pass1_dir.mkdir(exist_ok=True)

    for chr_num in range(1, 23):
        prefix = input_dir / f"chr{chr_num}"
        ext = ".bed" if use_bed else ".pgen"
        if (prefix.with_suffix(ext)).exists():
            print(f"[Pass1 Chr{chr_num}] removing high-missingness variants...")
            input_flag = "--bfile" if use_bed else "--pfile"
            run([
                "plink2",
                input_flag, str(prefix),
                "--exclude", str(fail_variants),
                "--make-pgen",
                "--out", str(pass1_dir / f"chr{chr_num}"),
                "--threads", str(threads),
                "--memory", str(memory)
            ])

    # --- Pass 2: sample missingness on variant-filtered data ---
    # Compute per-chromosome missingness, then aggregate in Python (avoids expensive concat)
    print("\nPass 2: Computing per-chromosome sample missingness on variant-filtered data...")
    from collections import defaultdict
    sample_missing = defaultdict(int)  # IID -> total missing count
    sample_obs = defaultdict(int)      # IID -> total observed count

    for chr_num in range(1, 23):
        prefix = pass1_dir / f"chr{chr_num}"
        if prefix.with_suffix(".pgen").exists():
            run([
                "plink2",
                "--pfile", str(prefix),
                "--missing",
                "--out", str(pass1_dir / f"chr{chr_num}_miss"),
                "--threads", str(threads),
                "--memory", str(memory)
            ])
            smiss_chr = pass1_dir / f"chr{chr_num}_miss.smiss"
            with smiss_chr.open() as f:
                header = f.readline()
                for line in f:
                    parts = line.strip().split()
                    if len(parts) >= 4:
                        iid = parts[0]
                        sample_missing[iid] += int(parts[1])
                        sample_obs[iid] += int(parts[2])

    # Write aggregated sample missingness and fail list
    smiss = out_dir / "sample_missingness.txt"
    fail_samples = out_dir / "fail_samples.txt"
    nsamp_total = len(sample_missing)
    nsamp = 0

    with smiss.open("w") as sf, fail_samples.open("w") as ff:
        sf.write("#IID\tMISSING_CT\tOBS_CT\tF_MISS\n")
        for iid in sorted(sample_missing):
            miss_ct = sample_missing[iid]
            obs_ct = sample_obs[iid]
            f_miss = miss_ct / obs_ct if obs_ct > 0 else 0
            sf.write(f"{iid}\t{miss_ct}\t{obs_ct}\t{f_miss:.6f}\n")
            if f_miss > args.sample_miss:
                ff.write(f"{iid}\n")
                nsamp += 1

    nsamp_keep = nsamp_total - nsamp
    print(f"Pass 2: {nsamp_total} total samples, {nsamp} fail missingness > {args.sample_miss}, {nsamp_keep} remaining")

    # Apply sample filter per chromosome (from pass1_filtered -> final output)
    print(f"\nApplying sample filter to produce final output in {out_dir}...")
    for chr_num in range(1, 23):
        prefix = pass1_dir / f"chr{chr_num}"
        if prefix.with_suffix(".pgen").exists():
            print(f"[Chr{chr_num}] removing high-missingness samples...")
            run([
                "plink2",
                "--pfile", str(prefix),
                "--remove", str(fail_samples),
                "--make-pgen",
                "--out", str(out_dir / f"chr{chr_num}"),
                "--threads", str(threads),
                "--memory", str(memory)
            ])

    # Clean up intermediate pass1 files
    import shutil
    shutil.rmtree(pass1_dir)
    for f in (out_dir / "pass1.smiss", out_dir / "pass1.vmiss",
              out_dir / "pass1.log", out_dir / "pass1.pgen", out_dir / "pass1.pvar",
              out_dir / "pass1.psam"):
        if f.exists():
            f.unlink()

    print(f"\nMissingness QC complete.")
    print(f"  Variants: {nvar_total} -> {nvar_keep} ({nvar} removed, missingness > {args.variant_miss})")
    print(f"  Samples:  {nsamp_total} -> {nsamp_keep} ({nsamp} removed, missingness > {args.sample_miss})")
    print(f"Filtered files in {out_dir}/")


if __name__ == "__main__":
    main()
