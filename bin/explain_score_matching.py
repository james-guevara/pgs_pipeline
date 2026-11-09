#!/usr/bin/env python3
"""Explain where requested PGS variants were lost before or during scoring."""

import argparse
import csv
from collections import Counter
from pathlib import Path


def rows(path):
    with Path(path).open(newline="") as handle:
        yield from csv.reader(handle, delimiter="\t")


def pvar(path, requested):
    variants = {}
    for row in rows(path):
        if not row or row[0].startswith("##"):
            continue
        if row[0].startswith("#"):
            continue
        if row[2] in requested:
            variants[row[2]] = (row[3].upper(), row[4].upper())
    return variants


def relevant_mapping(path, requested):
    """Keep only mappings affecting requested IDs, preserving last-key-wins."""
    old_to_new = {}
    mapped_rsids = set()
    remapped_requested = set()
    map_present = False
    with Path(path).open() as handle:
        for line in handle:
            row = line.split()
            if len(row) < 2 or row[0].startswith("#"):
                continue
            old, new = row[:2]
            map_present = True
            if old in requested:
                remapped_requested.add(old)
            if new in requested:
                old_to_new[old] = new
                mapped_rsids.add(new)
            else:
                old_to_new.pop(old, None)
    return old_to_new, mapped_rsids, remapped_requested, map_present


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

    id_index = args.id_col - 1
    allele_index = args.allele_col - 1
    # Count requested IDs, rather than storing the full weight rows.
    weight_rows = rows(args.weights)
    next(weight_rows)
    seen = Counter(row[id_index] for row in weight_rows)
    requested = set(seen)
    old_to_new, mapped_rsids, remapped_requested, map_present = relevant_mapping(
        args.rsid_map, requested
    )
    normalized_source = set()
    for row in rows(args.source_pvar):
        if not row or row[0].startswith("#"):
            continue
        variant_id = row[2]
        if variant_id in old_to_new:
            normalized_source.add(old_to_new[variant_id])
        elif variant_id in requested and variant_id not in remapped_requested:
            normalized_source.add(variant_id)
    scoring = pvar(args.score_pvar, requested)
    matched = {row[0] for row in rows(args.matched_vars) if row and row[0] in requested}
    counts = Counter()
    with Path(args.details).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["trait", "variant_id", "effect_allele", "category", "weight_rows_for_id"])
        weight_rows = rows(args.weights)
        next(weight_rows)
        for row in weight_rows:
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
                    if not map_present or variant_id in mapped_rsids
                    else "source_present_without_selected_mapping"
                )
            elif variant_id in mapped_rsids:
                category = "absent_from_harmonized_source"
            else:
                category = "no_selected_rsid_mapping_or_absent"
            counts[category] += 1
            if category != "matched":
                writer.writerow([args.trait, variant_id, effect_allele, category, seen[variant_id]])

    with Path(args.summary).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["trait", "category", "variants"])
        for category in sorted(counts):
            writer.writerow([args.trait, category, counts[category]])
        writer.writerow([args.trait, "requested_total", sum(seen.values())])



if __name__ == "__main__":
    main()
