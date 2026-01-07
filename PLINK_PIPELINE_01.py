#!/usr/bin/env python3
import subprocess
import argparse
from pathlib import Path

from config import load_config, get


def run(cmd):
    """Run a shell command and print it before executing."""
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)


def filter_vcf(chr_num, vcf, out_dir, threads, memory, mac, geno, vcf_half_call="missing", r2=None, aq=None):
    Path(out_dir).mkdir(exist_ok=True)
    cmd = [
        "plink2",
        "--vcf", vcf,
        "--vcf-half-call", vcf_half_call,
        "--snps-only", "just-acgt",
        "--max-alleles", "2",
        "--var-filter",
        "--mac", str(mac),
        "--geno", str(geno),
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


def annotate(chr_num, filter_dir, out_dir, rsid_maps_dir, threads, memory):
    Path(out_dir).mkdir(exist_ok=True)
    map_file = f"{rsid_maps_dir}/chr{chr_num}.map"
    cmd = [
        "plink2",
        "--pfile", f"{filter_dir}/chr{chr_num}",
        "--update-name", map_file, "2", "1",
        "--extract", f"<(cut -f2 {map_file})",
        "--write-snplist",
        "--make-pgen",
        "--out", f"{out_dir}/chr{chr_num}",
        "--threads", str(threads),
        "--memory", str(memory)
    ]
    run(["bash", "-c", " ".join(cmd)])  # allow process substitution


def maf_filter(chr_num, annotate_dir, out_dir, maf, threads, memory):
    Path(out_dir).mkdir(exist_ok=True)
    cmd = [
        "plink2",
        "--pfile", f"{annotate_dir}/chr{chr_num}",
        "--maf", str(maf),
        "--make-pgen",
        "--out", f"{out_dir}/chr{chr_num}",
        "--threads", str(threads),
        "--memory", str(memory)
    ]
    run(cmd)


def main():
    # Load config
    cfg = load_config()

    parser = argparse.ArgumentParser(description="Run PLINK pipeline in Python.")
    parser.add_argument("--chr", required=True, help="Chromosome number")
    parser.add_argument("--vcf", required=True, help="Path to VCF file")
    parser.add_argument("--config", type=Path, help="Path to config.toml")
    parser.add_argument("--threads", type=int, default=get(cfg, "defaults", "threads"))
    parser.add_argument("--memory", type=int, default=get(cfg, "defaults", "memory"))
    parser.add_argument("--r2", type=float, default=None, help="INFO R2 filter threshold")
    parser.add_argument("--aq", type=float, default=None, help="INFO AQ filter threshold")
    parser.add_argument("--maf", type=float, default=get(cfg, "filters", "maf"))
    parser.add_argument("--mac", type=int, default=get(cfg, "filters", "mac"))
    parser.add_argument("--geno", type=float, default=get(cfg, "filters", "geno"))
    parser.add_argument("--rsid-maps", default=get(cfg, "paths", "rsid_maps"),
                        help="Directory containing rsID map files")
    args = parser.parse_args()

    # Reload config if custom path provided
    if args.config:
        cfg = load_config(args.config)

    # Get directory names from config
    filter_dir = get(cfg, "directories", "filter")
    annotate_dir = get(cfg, "directories", "annotate")
    maf_dir = get(cfg, "directories", "maf")

    filter_vcf(args.chr, args.vcf, out_dir=filter_dir, threads=args.threads,
               memory=args.memory, mac=args.mac, geno=args.geno, r2=args.r2, aq=args.aq)
    annotate(args.chr, filter_dir=filter_dir, out_dir=annotate_dir,
             rsid_maps_dir=args.rsid_maps, threads=args.threads, memory=args.memory)
    maf_filter(args.chr, annotate_dir=annotate_dir, out_dir=maf_dir,
               maf=args.maf, threads=args.threads, memory=args.memory)

if __name__ == "__main__":
    main()
