#!/usr/bin/env python3
"""Step 01d: Merge multiple batches per chromosome using PLINK 1.9.

Expects bed/bim/fam files in batch subdirectories under the MAF output directory
(e.g., 01c_maf/batch1/chr1.bed, 01c_maf/batch2/chr1.bed, ...).

Produces merged bed/bim/fam files in the merge output directory
(e.g., 01d_merge/chr1.bed).
"""
import subprocess
import argparse
from pathlib import Path

from config import load_config, get


def run(cmd):
    """Run a shell command and print it before executing."""
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)


def detect_batches(maf_dir):
    """Auto-detect batch subdirectories in the MAF output directory."""
    maf_path = Path(maf_dir)
    batches = sorted([d.name for d in maf_path.iterdir() if d.is_dir()])
    if not batches:
        raise FileNotFoundError(f"No batch subdirectories found in {maf_dir}")
    return batches


def merge_chromosome(chr_num, maf_dir, out_dir, batches, threads, memory):
    """Merge a single chromosome across batches using PLINK 1.9."""
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # Verify bed files exist for all batches
    bed_files = []
    for batch in batches:
        bed = Path(maf_dir) / batch / f"chr{chr_num}.bed"
        if not bed.exists():
            print(f"Warning: {bed} not found, skipping batch {batch} for chr{chr_num}")
            continue
        bed_files.append(str(bed.with_suffix("")))  # prefix without extension

    if len(bed_files) < 2:
        if len(bed_files) == 1:
            print(f"Only one batch for chr{chr_num}, copying instead of merging.")
            prefix = bed_files[0]
            for ext in (".bed", ".bim", ".fam"):
                src = Path(f"{prefix}{ext}")
                dst = out_path / f"chr{chr_num}{ext}"
                if src.exists():
                    import shutil
                    shutil.copy2(src, dst)
            # Convert copied bed to pgen
            out_prefix = str(out_path / f"chr{chr_num}")
            pgen_cmd = [
                "plink2",
                "--bfile", out_prefix,
                "--make-pgen",
                "--out", out_prefix,
                "--threads", str(threads),
                "--memory", str(memory)
            ]
            run(pgen_cmd)
            return
        else:
            print(f"No bed files found for chr{chr_num}, skipping.")
            return

    # First file is the main input, rest go in merge list
    main_prefix = bed_files[0]
    merge_list = out_path / f"merge_list_chr{chr_num}.txt"
    with merge_list.open("w") as f:
        for prefix in bed_files[1:]:
            f.write(f"{prefix}.bed {prefix}.bim {prefix}.fam\n")

    out_prefix = str(out_path / f"chr{chr_num}")
    cmd = [
        "plink",
        "--bfile", main_prefix,
        "--merge-list", str(merge_list),
        "--make-bed",
        "--out", out_prefix,
        "--threads", str(threads),
        "--memory", str(memory)
    ]
    run(cmd)

    # Clean up merge list
    if merge_list.exists():
        merge_list.unlink()

    # Convert merged bed to pgen
    pgen_cmd = [
        "plink2",
        "--bfile", out_prefix,
        "--make-pgen",
        "--out", out_prefix,
        "--threads", str(threads),
        "--memory", str(memory)
    ]
    run(pgen_cmd)


def main():
    cfg = load_config()

    parser = argparse.ArgumentParser(
        description="Step 01d: Merge batches per chromosome using PLINK 1.9.")
    parser.add_argument("--config", type=Path, help="Path to config.toml")
    parser.add_argument("--maf-dir", type=Path,
                        default=Path(get(cfg, "directories", "maf")),
                        help="Directory containing batch subdirectories with bed/bim/fam files")
    parser.add_argument("--out-dir", type=Path,
                        default=Path(get(cfg, "directories", "merge")),
                        help="Output directory for merged files")
    parser.add_argument("--batches", nargs="+", default=None,
                        help="Batch names to merge (default: auto-detect)")
    parser.add_argument("--chr", type=int, default=None,
                        help="Single chromosome to merge (default: all 1-22)")
    parser.add_argument("--threads", type=int, default=get(cfg, "defaults", "threads"))
    parser.add_argument("--memory", type=int, default=get(cfg, "defaults", "memory"))
    args = parser.parse_args()

    if args.config:
        cfg = load_config(args.config)

    maf_dir = args.maf_dir
    out_dir = args.out_dir
    batches = args.batches or detect_batches(maf_dir)

    print(f"Merging batches: {', '.join(batches)}")
    print(f"Input directory: {maf_dir}")
    print(f"Output directory: {out_dir}\n")

    if args.chr is not None:
        chroms = [args.chr]
    else:
        chroms = list(range(1, 23))

    for chr_num in chroms:
        print(f"\n=== Chromosome {chr_num} ===")
        merge_chromosome(chr_num, maf_dir, out_dir, batches, args.threads, args.memory)

    print(f"\nMerge complete. Merged files in {out_dir}/")


if __name__ == "__main__":
    main()
