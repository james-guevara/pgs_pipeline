#!/usr/bin/env python3
import subprocess
from pathlib import Path
import argparse


def run(cmd):
    """Run and print shell commands."""
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser(description="Run PLINK2 polygenic score calculations for multiple traits.")
    parser.add_argument("--pfile-prefix", type=Path, required=True,
                        help="Prefix for the merged PLINK2 fileset (e.g., 05_summary/cohort).")
    parser.add_argument("--score-list", type=Path, default=Path("sbayesrc_sumstats_filepaths.txt"),
                        help="Two-column file: <trait_name> <score_file_path>.")
    parser.add_argument("--out-dir", type=Path, default=Path("scores"),
                        help="Output directory for .sscore files.")
    parser.add_argument("--threads", type=int, default=32, help="Number of threads for PLINK2.")
    parser.add_argument("--memory", type=int, default=16000, help="Memory (MB) for PLINK2.")
    args = parser.parse_args()

    pfile_prefix = args.pfile_prefix
    score_list = args.score_list
    out_dir = args.out_dir
    threads = args.threads
    memory = args.memory

    out_dir.mkdir(exist_ok=True)

    # Read the two-column file: <trait_name> <score_file>
    score_entries = []
    with score_list.open() as f:
        for line in f:
            if line.strip():
                parts = line.strip().split()
                if len(parts) != 2:
                    raise ValueError(f"Invalid line in score list: {line.strip()}")
                trait, score_file = parts
                score_entries.append((trait, Path(score_file)))

    if not score_entries:
        raise ValueError("No valid entries found in the score list.")

    print(f"Found {len(score_entries)} score files to process.\n")

    for trait, score_file in score_entries:
        out_prefix = out_dir / trait

        print(f"[Running] {trait}")
        run([
            "plink2",
            "--pfile", str(pfile_prefix),
            "--score", str(score_file), "1", "2", "3", "header", "center", "list-variants", "no-mean-imputation",
            "--threads", str(threads),
            "--memory", str(memory),
            "--out", str(out_prefix)
        ])

    print(f"\n✅ All {len(score_entries)} polygenic scores computed. Results in {out_dir}/")


if __name__ == "__main__":
    main()
