#!/usr/bin/env python3
"""Explain where requested PGS variants were lost before or during scoring."""

import argparse
import csv
from collections import Counter
from pathlib import Path


def rows(path):
    with Path(path).open(newline="") as handle:
        yield from csv.reader(handle, delimiter="\t")


def pvar(path):
    variants = {}
    for row in rows(path):
        if not row or row[0].startswith("##"):
            continue
        if row[0].startswith("#"):
            continue
        variants[row[2]] = (row[3].upper(), row[4].upper())
    return variants


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trait", required=True)
    parser.add_argument("--weights", required=True)
    parser.add_argument("--id-col", required=True, type=int)
    parser.add_argument("--allele-col", required=True, type=int)
    parser.add_argument("--source-pvar", required=True)
    parser.add_argument("--score-pvar", required=True)
    parser.add_argument("--rsid-map", required=True)
    parser.add_argument("--matched-vars", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--details", required=True)
    args = parser.parse_args()

    old_to_new = {}
    mapped_rsids = set()
    for row in rows(args.rsid_map):
        if len(row) >= 2 and not row[0].startswith("#"):
            old_to_new[row[0]] = row[1]
            mapped_rsids.add(row[1])

    source = pvar(args.source_pvar)
    normalized_source = {old_to_new.get(variant_id, variant_id): alleles for variant_id, alleles in source.items()}
    scoring = pvar(args.score_pvar)
    matched = {row[0] for row in rows(args.matched_vars) if row}

    weight_rows = list(rows(args.weights))
    header = weight_rows[0]
    requested = weight_rows[1:]
    id_index = args.id_col - 1
    allele_index = args.allele_col - 1
    counts = Counter()
    details = []
    seen = Counter(row[id_index] for row in requested)

    for row in requested:
        variant_id = row[id_index]
        effect_allele = row[allele_index].upper()
        if variant_id in matched:
            category = "matched"
        elif variant_id in scoring:
            ref, alt = scoring[variant_id]
            if effect_allele not in {ref, alt}:
                category = "effect_allele_incompatible"
            else:
                category = "present_but_not_scored_other"
        elif variant_id in normalized_source:
            category = (
                "removed_by_score_view_maf"
                if not mapped_rsids or variant_id in mapped_rsids
                else "source_present_without_selected_mapping"
            )
        elif variant_id in mapped_rsids:
            category = "absent_from_harmonized_source"
        else:
            category = "no_selected_rsid_mapping_or_absent"
        counts[category] += 1
        if category != "matched":
            details.append((variant_id, effect_allele, category, seen[variant_id]))

    with Path(args.summary).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["trait", "category", "variants"])
        for category in sorted(counts):
            writer.writerow([args.trait, category, counts[category]])
        writer.writerow([args.trait, "requested_total", len(requested)])

    with Path(args.details).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["trait", "variant_id", "effect_allele", "category", "weight_rows_for_id"])
        for variant_id, effect_allele, category, duplicates in details:
            writer.writerow([args.trait, variant_id, effect_allele, category, duplicates])


if __name__ == "__main__":
    main()
