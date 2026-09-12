#!/usr/bin/env python3
"""Join PGS, global ancestry, and within-ancestry PCA outputs by sample ID."""

import argparse
import csv
from pathlib import Path


ANCESTRIES = ("AFR", "AMR", "EAS", "EUR", "SAS")


def read_rows(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"{path} has no header")
        return reader.fieldnames, list(reader)


def sample_id(row):
    return row.get("IID") or row.get("#IID")


def index_rows(path):
    fields, rows = read_rows(path)
    id_candidates = ("IID", "#IID")
    resolved_id_column = next((candidate for candidate in id_candidates if candidate in fields), None)
    if resolved_id_column is None:
        raise ValueError(f"{path} is missing a sample identifier column ({', '.join(id_candidates)})")
    indexed = {}
    for row in rows:
        sample = row[resolved_id_column]
        if not sample:
            raise ValueError(f"{path} contains an empty sample ID")
        if sample in indexed:
            raise ValueError(f"Duplicate sample ID {sample} in {path}")
        indexed[sample] = row
    return fields, rows, indexed


def add_dictionary(rows, variable, data_type, nullable, description, source):
    rows.append({
        "variable": variable,
        "data_type": data_type,
        "nullable": "true" if nullable else "false",
        "description": description,
        "source": source,
    })


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--global-pcs", type=Path, required=True)
    parser.add_argument("--ancestry", type=Path, required=True)
    parser.add_argument("--within-dir", type=Path, required=True)
    parser.add_argument("--num-global-pcs", type=int, required=True)
    parser.add_argument("--num-within-pcs", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--dictionary", type=Path, required=True)
    args = parser.parse_args()

    score_fields, score_rows, scores = index_rows(args.scores)
    _, _, global_pcs = index_rows(args.global_pcs)
    ancestry_fields, _, ancestries = index_rows(args.ancestry)
    within_group_field = (
        "MOST_LIKELY_ANCESTRY"
        if "MOST_LIKELY_ANCESTRY" in ancestry_fields
        else "ANCESTRY"
    )
    score_traits = [field for field in score_fields if field not in {"IID", "#IID", "FID", "#FID"}]

    expected_ids = set(scores)
    for label, indexed in (("global PCs", global_pcs), ("ancestry", ancestries)):
        if set(indexed) != expected_ids:
            missing = sorted(expected_ids - set(indexed))[:5]
            extra = sorted(set(indexed) - expected_ids)[:5]
            raise ValueError(f"Sample mismatch for {label}; missing={missing}, extra={extra}")

    status_fields, status_rows = read_rows(args.within_dir / "status.tsv")
    required_status = {
        "ancestry", "assigned_samples", "unrelated_training_samples",
        "pruned_variants", "pcs", "status", "reliability", "reason",
    }
    if not required_status.issubset(status_fields):
        raise ValueError("Within-ancestry status file is missing required columns")
    status_by_group = {row["ancestry"]: row for row in status_rows}

    within_by_sample = {}
    for group in ANCESTRIES:
        pcs_path = args.within_dir / group / "pcs.tsv"
        status = status_by_group.get(group)
        if status and status["status"] == "completed" and not pcs_path.exists():
            raise ValueError(f"Completed ancestry group {group} has no pcs.tsv")
        if not pcs_path.exists():
            continue
        _, rows, indexed = index_rows(pcs_path)
        for sample, row in indexed.items():
            if sample not in expected_ids:
                raise ValueError(f"Within-ancestry sample {sample} is absent from scores")
            if ancestries[sample][within_group_field] != group:
                raise ValueError(
                    f"Within-ancestry group mismatch for {sample}: "
                    f"{group} vs {ancestries[sample][within_group_field]}"
                )
            if sample in within_by_sample:
                raise ValueError(f"Sample {sample} occurs in multiple within-ancestry files")
            within_by_sample[sample] = row

    probability_fields = [field for field in ancestry_fields if field.startswith("PROB_")]
    output_fields = ["IID", "ANCESTRY", "ANCESTRY_MOST_LIKELY", "ANCESTRY_MAX_PROBABILITY"]
    output_fields.extend(f"ANCESTRY_{field}" for field in probability_fields)
    output_fields.extend(f"GLOBAL_PC{index}" for index in range(1, args.num_global_pcs + 1))
    output_fields.extend([
        "WITHIN_PCA_STATUS", "WITHIN_PCA_RELIABILITY", "WITHIN_PCA_REASON",
        "WITHIN_PCA_N_ASSIGNED", "WITHIN_PCA_N_TRAINING",
        "WITHIN_PCA_N_VARIANTS", "WITHIN_PCA_N_PCS",
    ])
    output_fields.extend(f"WITHIN_PC{index}" for index in range(1, args.num_within_pcs + 1))
    output_fields.extend(f"PGS_{trait}" for trait in score_traits)

    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for score_row in score_rows:
            sample = sample_id(score_row)
            ancestry = ancestries[sample]
            global_row = global_pcs[sample]
            result = {
                "IID": sample,
                "ANCESTRY": ancestry["ANCESTRY"],
                "ANCESTRY_MOST_LIKELY": ancestry[within_group_field],
                "ANCESTRY_MAX_PROBABILITY": ancestry["MAX_PROBABILITY"],
            }
            for field in probability_fields:
                result[f"ANCESTRY_{field}"] = ancestry[field]
            for index in range(1, args.num_global_pcs + 1):
                result[f"GLOBAL_PC{index}"] = global_row[f"PC{index}_AVG"]

            status = status_by_group.get(ancestry[within_group_field])
            within_row = within_by_sample.get(sample)
            if status:
                result.update({
                    "WITHIN_PCA_STATUS": status["status"],
                    "WITHIN_PCA_RELIABILITY": status["reliability"],
                    "WITHIN_PCA_REASON": "" if status["reason"] == "." else status["reason"],
                    "WITHIN_PCA_N_ASSIGNED": status["assigned_samples"],
                    "WITHIN_PCA_N_TRAINING": status["unrelated_training_samples"],
                    "WITHIN_PCA_N_VARIANTS": status["pruned_variants"],
                    "WITHIN_PCA_N_PCS": status["pcs"],
                })
            if within_row:
                for index in range(1, args.num_within_pcs + 1):
                    result[f"WITHIN_PC{index}"] = within_row.get(f"PC{index}_AVG", "")
            for trait in score_traits:
                result[f"PGS_{trait}"] = score_row[trait]
            writer.writerow(result)

    dictionary = []
    add_dictionary(dictionary, "IID", "string", False, "Participant/sample identifier.", "All sample-level outputs")
    add_dictionary(dictionary, "ANCESTRY", "categorical", False, "Assigned global ancestry or uncertain label.", "ancestry_probabilities.tsv")
    add_dictionary(dictionary, "ANCESTRY_MOST_LIKELY", "categorical", False, "Highest-probability ancestry used to define the within-ancestry PCA group.", "ancestry_probabilities.tsv")
    add_dictionary(dictionary, "ANCESTRY_MAX_PROBABILITY", "float", False, "Highest ancestry assignment probability.", "ancestry_probabilities.tsv")
    for field in probability_fields:
        group = field.removeprefix("PROB_")
        add_dictionary(dictionary, f"ANCESTRY_{field}", "float", False, f"Probability of {group} ancestry.", "ancestry_probabilities.tsv")
    for index in range(1, args.num_global_pcs + 1):
        add_dictionary(dictionary, f"GLOBAL_PC{index}", "float", False, f"PC{index} projected from the fixed 1000 Genomes reference.", "global_pcs.tsv")
    add_dictionary(dictionary, "WITHIN_PCA_STATUS", "categorical", True, "Within-ancestry PCA completion status; missing for uncertain ancestry.", "within_ancestry/status.tsv")
    add_dictionary(dictionary, "WITHIN_PCA_RELIABILITY", "categorical", True, "Within-ancestry PCA reliability flag; missing for uncertain ancestry.", "within_ancestry/status.tsv")
    add_dictionary(dictionary, "WITHIN_PCA_REASON", "string", True, "Reason for a nonstandard within-ancestry PCA status.", "within_ancestry/status.tsv")
    add_dictionary(dictionary, "WITHIN_PCA_N_ASSIGNED", "integer", True, "Number assigned to the participant's ancestry group.", "within_ancestry/status.tsv")
    add_dictionary(dictionary, "WITHIN_PCA_N_TRAINING", "integer", True, "Unrelated samples used to train the within-ancestry PCA.", "within_ancestry/status.tsv")
    add_dictionary(dictionary, "WITHIN_PCA_N_VARIANTS", "integer", True, "Polymorphic LD-pruned variants used for within-ancestry PCA.", "within_ancestry/status.tsv")
    add_dictionary(dictionary, "WITHIN_PCA_N_PCS", "integer", True, "Number of within-ancestry PCs available for the group.", "within_ancestry/status.tsv")
    for index in range(1, args.num_within_pcs + 1):
        add_dictionary(dictionary, f"WITHIN_PC{index}", "float", True, f"Within-ancestry PC{index}; comparable only within the assigned ancestry.", "within_ancestry/<group>/pcs.tsv")
    for trait in score_traits:
        add_dictionary(dictionary, f"PGS_{trait}", "float", False, f"Centered PLINK SCORE1_AVG polygenic score for {trait}; not residualized or standardized.", "combined_scores.tsv")

    with args.dictionary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["variable", "data_type", "nullable", "description", "source"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(dictionary)


if __name__ == "__main__":
    main()
