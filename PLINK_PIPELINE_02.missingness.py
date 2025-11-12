#!/usr/bin/env python3
import subprocess
from pathlib import Path
import argparse


def run(cmd):
    """Run and print shell commands."""
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser(description="Run PLINK2 missingness QC + filtering.")
    parser.add_argument("--dir", type=Path, default=Path("03_maf"),
                        help="Directory containing per-chromosome .pgen/.psam/.pvar files (default: 03_maf)")
    parser.add_argument("--out-dir", type=Path, default=Path("04_missingness"),
                        help="Output directory for missingness results (default: 04_missingness)")
    parser.add_argument("--memory", type=int, default=4000)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()

    input_dir = args.dir
    out_dir = args.out_dir
    memory = args.memory
    threads = args.threads

    out_dir.mkdir(exist_ok=True)
    merge_list = out_dir / "merge_list.txt"
    merge_list.write_text("")

    print(f"Building merge list from {input_dir}...")
    missing_chr = []
    with merge_list.open("a") as f:
        for chr_num in range(1, 23):
            prefix = input_dir / f"chr{chr_num}"
            if (prefix.with_suffix(".pgen")).exists():
                f.write(str(prefix) + "\n")
            else:
                missing_chr.append(chr_num)

    if missing_chr:
        print(f"⚠️  Skipping missing chromosomes: {' '.join(map(str, missing_chr))}")
    print(f"Merge list written to {merge_list}\n")

    # Run missingness
    print("Running PLINK2 missingness...")
    run([
        "plink2",
        "--pmerge-list", str(merge_list),
        "--missing",
        "--out", str(out_dir / "missingness"),
        "--threads", str(threads),
        "--memory", str(memory)
    ])

    # Generate fail lists
    print("Generating fail lists...")
    smiss = out_dir / "missingness.smiss"
    vmiss = out_dir / "missingness.vmiss"
    fail_samples = out_dir / "fail_samples.txt"
    fail_variants = out_dir / "fail_variants.txt"

    subprocess.run(f"awk 'NR>1 && $4>0.05 {{print $1}}' {smiss} > {fail_samples}", shell=True, check=True)
    subprocess.run(f"awk 'NR>1 && $5>0.05 {{print $2}}' {vmiss} > {fail_variants}", shell=True, check=True)

    nsamp_proc = subprocess.run(f"wc -l < {fail_samples}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    nsamp = int(nsamp_proc.stdout.decode().strip() or 0)
    
    nvar_proc = subprocess.run(f"wc -l < {fail_variants}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    nvar = int(nvar_proc.stdout.decode().strip() or 0)

    # Apply filters per chromosome
    print(f"Applying missingness filters to {input_dir}...")
    for chr_num in range(1, 23):
        prefix = input_dir / f"chr{chr_num}"
        if (prefix.with_suffix(".pgen")).exists():
            print(f"[Chr{chr_num}] filtering...")
            run([
                "plink2",
                "--pfile", str(prefix),
                "--remove", str(fail_samples),
                "--exclude", str(fail_variants),
                "--make-pgen",
                "--out", str(out_dir / f"chr{chr_num}"),
                "--threads", str(threads),
                "--memory", str(memory)
            ])

    print(f"✅ Missingness QC + filtering complete. Filtered files in {out_dir}/")


if __name__ == "__main__":
    main()
