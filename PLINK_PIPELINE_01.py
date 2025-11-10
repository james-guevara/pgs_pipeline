#!/usr/bin/env python3
import subprocess
import argparse
from pathlib import Path


def run(cmd):
    """Run a shell command and print it before executing."""
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)


def filter_vcf(chr_num, vcf, out_dir="01_filter", threads=16, memory=16000, r2=None, aq=None):
    Path(out_dir).mkdir(exist_ok=True)
    cmd = [
        "plink2",
        "--vcf", vcf,
        "--snps-only", "just-acgt",
        "--max-alleles", "2",
        "--var-filter",
        "--mac", "10",
        "--geno", "0.05",
        "--set-all-var-ids", "@:#:$r:$a",
        "--make-pgen",
        "--out", f"{out_dir}/chr{chr_num}",
        "--threads", str(threads),
        "--memory", str(memory)
    ]
    # Optional INFO filters
    if r2 is not None:
        cmd += ["--extract-if-info", "R2", ">=", str(r2)]
    elif aq is not None:
        cmd += ["--extract-if-info", "AQ", ">=", str(aq)]
    run(cmd)


def annotate(chr_num, out_dir="02_annotate", threads=16, memory=16000):
    Path(out_dir).mkdir(exist_ok=True)
    map_file = f"rsid_maps/chr{chr_num}.map"
    cmd = [
        "plink2",
        "--pfile", f"01_filter/chr{chr_num}",
        "--update-name", map_file, "2", "1",
        "--extract", f"<(cut -f2 {map_file})",
        "--write-snplist",
        "--make-pgen",
        "--out", f"{out_dir}/chr{chr_num}",
        "--threads", str(threads),
        "--memory", str(memory)
    ]
    run(["bash", "-c", " ".join(cmd)])  # allow process substitution


def maf_filter(chr_num, out_dir="03_maf", threads=16, memory=16000, maf=0.01):
    Path(out_dir).mkdir(exist_ok=True)
    cmd = [
        "plink2",
        "--pfile", f"02_annotate/chr{chr_num}",
        "--maf", str(maf),
        "--make-pgen",
        "--out", f"{out_dir}/chr{chr_num}",
        "--threads", str(threads),
        "--memory", str(memory)
    ]
    run(cmd)


def main():
    parser = argparse.ArgumentParser(description="Run PLINK pipeline in Python.")
    parser.add_argument("--chr", required=True, help="Chromosome number")
    parser.add_argument("--vcf", required=True, help="Path to VCF file")
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--memory", type=int, default=16000)
    parser.add_argument("--r2", type=float, default=None)
    parser.add_argument("--aq", type=float, default=None)
    parser.add_argument("--maf", type=float, default=0.01)
    args = parser.parse_args()

    filter_vcf(args.chr, args.vcf, threads=args.threads, memory=args.memory, r2=args.r2)
    annotate(args.chr, threads=args.threads, memory=args.memory)
    maf_filter(args.chr, threads=args.threads, memory=args.memory, maf=args.maf)

if __name__ == "__main__":
    main()
