#!/usr/bin/env python3
"""Merge disjoint-sample PGENs and publish compact harmonization QC."""

import argparse
import csv
import json
import os
from pathlib import Path
import re
import tempfile
import time

import numpy as np
import pgenlib as pg
import pyarrow as pa
import pyarrow.parquet as pq

from merge_pgen_samples import merge_many, samples, variants


def read_manifest(path):
    with Path(path).open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) < 2 or not {"source", "pgen_prefix"}.issubset(rows[0]):
        raise ValueError("Manifest requires source and pgen_prefix columns and at least two rows")
    names = [row["source"].strip() for row in rows]
    if any(not name for name in names) or len(names) != len(set(names)):
        raise ValueError("Manifest source names must be nonempty and unique")
    safe = [re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_").lower() for name in names]
    if any(not name for name in safe) or len(safe) != len(set(safe)):
        raise ValueError("Manifest source names must remain unique after column-name normalization")
    prefixes = [Path(row["pgen_prefix"]).expanduser() for row in rows]
    for prefix in prefixes:
        for suffix in (".pgen", ".pvar", ".psam"):
            if not Path(str(prefix) + suffix).is_file():
                raise ValueError(f"Missing input: {prefix}{suffix}")
    return names, safe, prefixes


def marker_stats(reader, index, sample_count, buffer):
    reader.read_dosages(index, buffer)
    valid = buffer != -9
    observed = int(valid.sum())
    if not observed:
        return None, 1.0, valid
    af = float(buffer[valid].sum() / (2 * observed))
    return min(af, 1 - af), 1 - observed / sample_count, valid


def write_qc(names, safe_names, prefixes, output, batch_size=10000):
    indexes = [variants(prefix)[1] for prefix in prefixes]
    sample_data = [samples(prefix) for prefix in prefixes]
    sample_counts = [len(item[1]) for item in sample_data]
    union = list(indexes[0])
    seen = set(union)
    for index in indexes[1:]:
        union.extend(key for key in index if key not in seen)
        seen.update(index)
    retained = [key for key in indexes[0] if all(key in index for index in indexes[1:])]
    retained_set = set(retained)
    readers = [pg.PgenReader(os.fsencode(str(prefix) + ".pgen"), raw_sample_ct=n,
                             variant_ct=len(index))
               for prefix, n, index in zip(prefixes, sample_counts, indexes)]
    buffers = [np.empty(n, dtype=np.float64) for n in sample_counts]
    sample_missing = [np.zeros(n, dtype=np.int64) for n in sample_counts]
    fields = [pa.field("chrom", pa.string()), pa.field("pos", pa.int64()),
              pa.field("ref", pa.string()), pa.field("alt", pa.string()),
              pa.field("id", pa.string())]
    for safe in safe_names:
        fields.extend([pa.field(f"{safe}_present", pa.bool_()),
                       pa.field(f"{safe}_maf", pa.float64()),
                       pa.field(f"{safe}_missingness", pa.float64())])
    fields.extend([pa.field("combined_maf", pa.float64()),
                   pa.field("combined_missingness", pa.float64()),
                   pa.field("retained", pa.bool_()),
                   pa.field("exclusion_reason", pa.string())])
    schema = pa.schema(fields)
    writer = pq.ParquetWriter(output / "marker_qc.parquet", schema, compression="zstd")
    rows = []
    try:
        for key in union:
            row = {"chrom": key[0], "pos": int(key[1]), "ref": key[2], "alt": key[3],
                   "id": ":".join(key)}
            missing_sources = []
            combined_alt = combined_observed = 0.0
            for name, safe, reader, index, n, buffer, missing in zip(
                    names, safe_names, readers, indexes, sample_counts, buffers, sample_missing):
                present = key in index
                row[f"{safe}_present"] = present
                if not present:
                    row[f"{safe}_maf"] = None
                    row[f"{safe}_missingness"] = None
                    missing_sources.append(name)
                    continue
                maf, fmiss, valid = marker_stats(reader, index[key], n, buffer)
                row[f"{safe}_maf"] = maf
                row[f"{safe}_missingness"] = fmiss
                if key in retained_set:
                    missing += ~valid
                    combined_alt += float(buffer[valid].sum())
                    combined_observed += int(valid.sum())
            row["combined_maf"] = (min(combined_alt / (2 * combined_observed),
                                       1 - combined_alt / (2 * combined_observed))
                                   if key in retained_set and combined_observed else None)
            total_samples = sum(sample_counts)
            row["combined_missingness"] = (1 - combined_observed / total_samples
                                            if key in retained_set else None)
            row["retained"] = key in retained_set
            row["exclusion_reason"] = ("" if not missing_sources else
                                       "missing_from:" + ",".join(missing_sources))
            rows.append(row)
            if len(rows) >= batch_size:
                table = pa.Table.from_pylist(rows, schema=schema)
                writer.write_table(table)
                rows = []
        if rows:
            table = pa.Table.from_pylist(rows, schema=schema)
            writer.write_table(table)
    finally:
        writer.close()
        for reader in readers:
            reader.close()

    with (output / "sample_qc.tsv").open("w", newline="") as handle:
        columns = ["FID", "IID", "source", "missingness", "retained", "exclusion_reason"]
        out = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        out.writeheader()
        denominator = len(retained)
        for name, (_, records, _), missing in zip(names, sample_data, sample_missing):
            for record, count in zip(records, missing):
                out.writerow({"FID": record.get("FID", "0"), "IID": record["IID"],
                              "source": name,
                              "missingness": count / denominator if denominator else None,
                              "retained": "true", "exclusion_reason": ""})
    return len(union), len(retained)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()
    names, safe_names, prefixes = read_manifest(args.manifest)
    destination = Path(args.out_dir)
    if destination.exists():
        raise ValueError(f"Output directory already exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    with tempfile.TemporaryDirectory(prefix=".harmonize-", dir=destination.parent) as temp:
        staging = Path(temp) / "release"
        merge_summary = merge_many(prefixes, staging, dosage_mode="preserve")
        union_count, retained_count = write_qc(names, safe_names, prefixes, staging)
        for suffix in (".pgen", ".pvar", ".psam"):
            (staging / f"merged{suffix}").rename(staging / f"harmonized{suffix}")
        summary = {
            "schema_version": "1.0",
            "status": "validated",
            "sources": names,
            "matching": merge_summary["matching"],
            "dosage_policy": merge_summary["dosage_policy"],
            "input_count": len(prefixes),
            "output_samples": merge_summary["output_samples"],
            "union_markers": union_count,
            "output_markers": retained_count,
            "excluded_markers": union_count - retained_count,
            "elapsed_seconds": round(time.monotonic() - started, 3),
            "outputs": {"pgen_prefix": "harmonized", "marker_qc": "marker_qc.parquet",
                        "sample_qc": "sample_qc.tsv"},
            "inputs": merge_summary["inputs"],
        }
        (staging / "harmonization_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        (staging / "merge_summary.json").unlink()
        staging.rename(destination)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
