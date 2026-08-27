#!/usr/bin/env python3
"""Match a cohort PVAR to a fixed PCA panel without strand guessing."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


COMPLEMENT = str.maketrans("ACGT", "TGCA")


def read_pvar(path):
    with path.open() as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                fields = line.lstrip("#").rstrip().split("\t")
                break
        else:
            raise ValueError(f"No #CHROM header in {path}")
        rows = []
        for line in handle:
            if not line.strip():
                continue
            values = line.rstrip().split("\t")
            rows.append(dict(zip(fields, values)))
    return rows


def classify(panel, cohort_rows):
    exact = [row for row in cohort_rows if row["REF"] == panel["REF"] and row["ALT"] == panel["ALT"]]
    swapped = [row for row in cohort_rows if row["REF"] == panel["ALT"] and row["ALT"] == panel["REF"]]
    compatible = exact + swapped
    if not cohort_rows:
        return "absent", None
    if len(compatible) > 1:
        return "duplicate_compatible", None
    if len(exact) == 1:
        return "exact", exact[0]
    if len(swapped) == 1:
        return "ref_alt_swapped", swapped[0]
    panel_ref_comp = panel["REF"].translate(COMPLEMENT)
    panel_alt_comp = panel["ALT"].translate(COMPLEMENT)
    if any(
        {row["REF"], row["ALT"]} == {panel_ref_comp, panel_alt_comp}
        for row in cohort_rows
    ):
        return "strand_complement_rejected", None
    return "incompatible_alleles", None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--cohort-pvar", type=Path, required=True)
    parser.add_argument("--min-overlap", type=float, required=True)
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    with args.panel.open(newline="") as handle:
        panel_rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"CHROM", "POS", "ID", "REF", "ALT"}
    if not panel_rows or not required.issubset(panel_rows[0]):
        raise ValueError("Panel must contain CHROM, POS, ID, REF, and ALT")

    cohort_by_locus = defaultdict(list)
    for row in read_pvar(args.cohort_pvar):
        cohort_by_locus[(row["CHROM"], row["POS"])].append(row)

    audit_rows = []
    usable = []
    counts = defaultdict(int)
    seen_panel_ids = set()
    for panel in panel_rows:
        if panel["ID"] in seen_panel_ids:
            status, cohort = "duplicate_panel_id", None
        else:
            seen_panel_ids.add(panel["ID"])
            status, cohort = classify(
                panel, cohort_by_locus[(panel["CHROM"], panel["POS"])]
            )
        counts[status] += 1
        is_usable = status in {"exact", "ref_alt_swapped"}
        audit_rows.append({
            "CHROM": panel["CHROM"],
            "POS": panel["POS"],
            "PANEL_ID": panel["ID"],
            "PANEL_REF": panel["REF"],
            "PANEL_ALT": panel["ALT"],
            "COHORT_ID": cohort["ID"] if cohort else "",
            "COHORT_REF": cohort["REF"] if cohort else "",
            "COHORT_ALT": cohort["ALT"] if cohort else "",
            "STATUS": status,
            "USABLE": str(is_usable).lower(),
        })
        if is_usable:
            usable.append((cohort["ID"], panel["ID"], panel["REF"]))

    audit_path = args.output_dir / "panel_overlap.tsv"
    with audit_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=audit_rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(audit_rows)

    with (args.output_dir / "usable_cohort_ids.txt").open("w") as ids, \
            (args.output_dir / "rename_to_panel_ids.tsv").open("w") as renames, \
            (args.output_dir / "reference_alleles.tsv").open("w") as alleles:
        for cohort_id, panel_id, panel_ref in usable:
            ids.write(f"{cohort_id}\n")
            renames.write(f"{panel_id}\t{cohort_id}\n")
            alleles.write(f"{panel_id}\t{panel_ref}\n")

    overlap = len(usable) / len(panel_rows)
    with (args.output_dir / "panel_overlap_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["panel_variants", len(panel_rows)])
        writer.writerow(["usable_variants", len(usable)])
        writer.writerow(["usable_fraction", f"{overlap:.8f}"])
        writer.writerow(["minimum_usable_fraction", args.min_overlap])
        for status in sorted(counts):
            writer.writerow([f"status_{status}", counts[status]])

    if overlap < args.min_overlap:
        raise SystemExit(
            f"PCA panel overlap {overlap:.4f} is below required {args.min_overlap:.4f}"
        )


if __name__ == "__main__":
    main()
