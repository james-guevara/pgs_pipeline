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

    parser = argparse.ArgumentParser(description="Run summary-level QC pipeline.")
    parser.add_argument("--config", type=Path, help="Path to config.toml")
    parser.add_argument("--cohort", default="cohort", help="Cohort name (used for output prefix)")
    parser.add_argument("--threads", type=int, default=get(cfg, "defaults", "threads"))
    parser.add_argument("--memory", type=int, default=get(cfg, "defaults", "memory"))
    parser.add_argument("--input-dir", type=Path, default=Path(get(cfg, "directories", "missingness")),
                        help="Directory containing per-chromosome filtered files")
    parser.add_argument("--out-dir", type=Path, default=Path(get(cfg, "directories", "summary")),
                        help="Output directory")
    args = parser.parse_args()

    cohort = args.cohort
    threads = args.threads
    memory = args.memory
    in_dir = args.input_dir
    out_dir = args.out_dir

    out_dir.mkdir(exist_ok=True)

    # Build merge list from input directory (MAF already filtered in step 01)
    merge_list = out_dir / "merge_list.txt"
    merge_list.write_text("")
    print("Building merge list...")

    missing_chr = []
    with merge_list.open("a") as f:
        for chr_num in range(1, 23):
            prefix = in_dir / f"chr{chr_num}"
            if (prefix.with_suffix(".pgen")).exists():
                f.write(str(prefix) + "\n")
            else:
                missing_chr.append(chr_num)

    if missing_chr:
        print(f"⚠️  Warning: missing chromosomes — {' '.join(map(str, missing_chr))}")

    # Merge for summary QC
    print("Merging filtered data for final QC...")
    run([
        "plink2",
        "--pmerge-list", str(merge_list),
        "--missing",
        "--hardy",
        "--freq",
        "--out", str(out_dir / cohort),
        "--threads", str(threads),
        "--memory", str(memory)
    ])

    # Quick summary counts
    smiss = out_dir / f"{cohort}.smiss"
    vmiss = out_dir / f"{cohort}.vmiss"
    summary_counts = out_dir / "summary_counts.txt"

    nsamp_proc = subprocess.run(f"wc -l < {smiss}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    nvar_proc = subprocess.run(f"wc -l < {vmiss}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    nsamp = max(int(nsamp_proc.stdout.decode().strip() or 1) - 1, 0)
    nvar = max(int(nvar_proc.stdout.decode().strip() or 1) - 1, 0)

    with summary_counts.open("w") as f:
        f.write("Metric\tCount\n")
        f.write(f"Samples\t{nsamp}\n")
        f.write(f"Variants\t{nvar}\n")

    # Clean up
    merge_list.unlink(missing_ok=True)
    print(f"✅ Summary QC complete. Results in {out_dir}/")


if __name__ == "__main__":
    main()
