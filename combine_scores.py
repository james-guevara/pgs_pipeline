#!/usr/bin/env python3
"""Combine multiple PLINK2 .sscore files into a single table."""

import argparse
from pathlib import Path
import polars as pl


def load_sscore(filepath: Path, trait_name: str | None = None) -> pl.DataFrame:
    """Load a PLINK2 .sscore file and extract IID + score."""
    if trait_name is None:
        trait_name = filepath.stem  # e.g., "height_2022" from "height_2022.sscore"

    df = pl.read_csv(filepath, separator="\t")

    # Rename columns: #IID -> IID, SCORE1_AVG -> trait_name
    df = df.select([
        pl.col("#IID").alias("IID"),
        pl.col("SCORE1_AVG").alias(trait_name)
    ])

    return df


def combine_scores(score_dir: Path, output: Path | None = None) -> pl.DataFrame:
    """Combine all .sscore files in a directory into one table."""
    score_files = sorted(score_dir.glob("*.sscore"))

    if not score_files:
        raise FileNotFoundError(f"No .sscore files found in {score_dir}")

    print(f"Found {len(score_files)} score files:")
    for f in score_files:
        print(f"  - {f.name}")

    # Load all score files
    dfs = [load_sscore(sf) for sf in score_files]

    # Concatenate horizontally by joining on IID
    combined = dfs[0]
    for df in dfs[1:]:
        combined = combined.join(df, on="IID", how="full", coalesce=True)

    print(f"\nCombined: {combined.shape[0]} samples × {combined.shape[1]} columns")

    if output:
        if output.suffix == ".parquet":
            combined.write_parquet(output)
        else:
            combined.write_csv(output, separator="\t")
        print(f"Saved to: {output}")

    return combined


def main():
    parser = argparse.ArgumentParser(description="Combine PLINK2 .sscore files into one table")
    parser.add_argument("score_dir", type=Path, help="Directory containing .sscore files")
    parser.add_argument("-o", "--output", type=Path, help="Output file (.tsv or .parquet)")
    parser.add_argument("--stats", action="store_true", help="Print summary statistics")
    args = parser.parse_args()

    combined = combine_scores(args.score_dir, args.output)

    if args.stats:
        print("\n=== Summary Statistics ===")
        # Get trait columns (all except IID)
        trait_cols = [c for c in combined.columns if c != "IID"]
        stats = combined.select(trait_cols).describe()
        print(stats)


if __name__ == "__main__":
    main()
