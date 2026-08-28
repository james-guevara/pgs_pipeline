#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def read_score(path: Path):
    trait = path.name.removesuffix(".sscore")
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "#IID" not in reader.fieldnames or "SCORE1_AVG" not in reader.fieldnames:
            raise ValueError(f"{path} must contain #IID and SCORE1_AVG columns")
        values = [(row["#IID"], row["SCORE1_AVG"]) for row in reader]
    return trait, values


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scores", nargs="+", type=Path, required=True)
    parser.add_argument("--qcs", nargs="+", type=Path, required=True)
    parser.add_argument("--scores-out", type=Path, required=True)
    parser.add_argument("--qc-out", type=Path, required=True)
    args = parser.parse_args()

    traits = []
    sample_order = None
    columns = {}
    for path in sorted(args.scores):
        trait, values = read_score(path)
        ids = [sample for sample, _ in values]
        if sample_order is None:
            sample_order = ids
        elif ids != sample_order:
            raise ValueError(f"Sample IDs or order in {path} differ from the other score files")
        if trait in columns:
            raise ValueError(f"Duplicate trait name: {trait}")
        traits.append(trait)
        columns[trait] = [value for _, value in values]

    with args.scores_out.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["IID", *traits])
        for index, sample in enumerate(sample_order or []):
            writer.writerow([sample, *(columns[trait][index] for trait in traits)])

    qc_rows = []
    fieldnames = None
    for path in sorted(args.qcs):
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if fieldnames is None:
                fieldnames = reader.fieldnames
            elif reader.fieldnames != fieldnames:
                raise ValueError(f"QC columns in {path} differ from the other QC files")
            qc_rows.extend(reader)
    with args.qc_out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames or [], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(sorted(qc_rows, key=lambda row: row["trait"]))


if __name__ == "__main__":
    main()
