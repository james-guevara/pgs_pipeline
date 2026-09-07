#!/usr/bin/env python3
"""Build a common, allele-identical marker release across PLINK 2 filesets."""

import argparse
import csv
import json
import subprocess
from pathlib import Path


def read_datasets(path: Path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"dataset_id", "pgen_prefix", "priority"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("Dataset manifest requires dataset_id, pgen_prefix, and priority")
    ids = [row["dataset_id"] for row in rows]
    if len(ids) != len(set(ids)):
        raise ValueError("dataset_id values must be unique")
    return sorted(rows, key=lambda row: (int(row["priority"]), row["dataset_id"]))


def run_plink(prefix: Path, dataset_id: str, sample_miss: float, output_dir: Path):
    for suffix in (".pgen", ".pvar", ".psam"):
        if not Path(f"{prefix}{suffix}").is_file():
            raise FileNotFoundError(f"Missing {prefix}{suffix}")
    dataset_dir = output_dir / "datasets" / dataset_id
    dataset_dir.mkdir(parents=True, exist_ok=True)
    missing_prefix = dataset_dir / "missingness"
    subprocess.run(
        ["plink2", "--pfile", str(prefix), "--missing", "--out", str(missing_prefix)],
        check=True,
    )
    remove_path = dataset_dir / "remove_samples.tsv"
    with Path(f"{missing_prefix}.smiss").open() as source, remove_path.open("w") as target:
        header = source.readline().lstrip("#").split()
        fid_i = header.index("FID") if "FID" in header else None
        iid_i = header.index("IID")
        miss_i = header.index("F_MISS")
        for line in source:
            fields = line.split()
            if float(fields[miss_i]) > sample_miss:
                fid = fields[fid_i] if fid_i is not None else "0"
                target.write(f"{fid}\t{fields[iid_i]}\n")
    frequency_prefix = dataset_dir / "frequency"
    command = ["plink2", "--pfile", str(prefix)]
    if remove_path.stat().st_size:
        command += ["--remove", str(remove_path)]
    command += ["--freq", "--out", str(frequency_prefix)]
    subprocess.run(command, check=True)
    return remove_path, Path(f"{prefix}.pvar"), Path(f"{frequency_prefix}.afreq")


def read_eligible(pvar: Path, afreq: Path, maf: float):
    frequencies = {}
    with afreq.open() as handle:
        header = handle.readline().lstrip("#").split()
        id_i = header.index("ID")
        freq_i = header.index("ALT_FREQS")
        for line in handle:
            fields = line.split()
            try:
                alt_frequency = float(fields[freq_i])
            except ValueError:
                continue
            frequencies[fields[id_i]] = min(alt_frequency, 1.0 - alt_frequency)

    eligible = {}
    with pvar.open() as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.lstrip("#").rstrip().split("\t")
                break
        else:
            raise ValueError(f"No #CHROM header in {pvar}")
        for line in handle:
            fields = dict(zip(header, line.rstrip().split("\t")))
            marker_id = fields["ID"]
            marker_maf = frequencies.get(marker_id)
            if marker_maf is None or marker_maf < maf:
                continue
            key = (fields["CHROM"].removeprefix("chr"), int(fields["POS"]), fields["REF"].upper(), fields["ALT"].upper())
            if key in eligible:
                eligible[key] = None
            else:
                eligible[key] = (marker_id, marker_maf)
    return {key: value for key, value in eligible.items() if value is not None}


def chromosome_key(chromosome: str):
    order = {"X": 23, "Y": 24, "XY": 25, "MT": 26, "M": 26}
    return (order.get(chromosome, int(chromosome) if chromosome.isdigit() else 10_000), chromosome)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--datasets", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--maf", type=float, default=0.01)
    parser.add_argument("--sample-miss", type=float, default=0.05)
    args = parser.parse_args()
    if not 0 <= args.maf <= 0.5 or not 0 <= args.sample_miss <= 1:
        raise ValueError("MAF must be in [0, 0.5] and sample missingness in [0, 1]")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    datasets = read_datasets(args.datasets)
    eligible = {}
    removals = {}
    for row in datasets:
        dataset_id = row["dataset_id"]
        remove_path, pvar, afreq = run_plink(
            Path(row["pgen_prefix"]), dataset_id, args.sample_miss, args.output_dir
        )
        eligible[dataset_id] = read_eligible(pvar, afreq, args.maf)
        removals[dataset_id] = remove_path

    common = set.intersection(*(set(markers) for markers in eligible.values()))
    common = sorted(common, key=lambda key: (chromosome_key(key[0]), key[1], key[2], key[3]))
    common_ids = args.output_dir / "common_markers.txt"
    detailed = args.output_dir / "common_variants.tsv"
    with common_ids.open("w") as ids, detailed.open("w", newline="") as table:
        columns = ["canonical_id", "chromosome", "position", "ref", "alt"]
        for row in datasets:
            columns += [f"{row['dataset_id']}_id", f"{row['dataset_id']}_maf"]
        writer = csv.DictWriter(table, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for chrom, pos, ref, alt in common:
            canonical = f"{chrom}:{pos}:{ref}:{alt}"
            ids.write(f"{canonical}\n")
            record = {"canonical_id": canonical, "chromosome": chrom, "position": pos, "ref": ref, "alt": alt}
            for row in datasets:
                marker_id, marker_maf = eligible[row["dataset_id"]][(chrom, pos, ref, alt)]
                record[f"{row['dataset_id']}_id"] = marker_id
                record[f"{row['dataset_id']}_maf"] = f"{marker_maf:.10g}"
            writer.writerow(record)

    release = {
        "schema_version": 1,
        "maf_minimum_each_dataset": args.maf,
        "sample_missingness_maximum": args.sample_miss,
        "dataset_count": len(datasets),
        "common_marker_count": len(common),
        "datasets": [
            {
                "dataset_id": row["dataset_id"],
                "pgen_prefix": row["pgen_prefix"],
                "priority": int(row["priority"]),
                "eligible_marker_count": len(eligible[row["dataset_id"]]),
                "sample_remove_file": str(Path("datasets") / row["dataset_id"] / "remove_samples.tsv"),
            }
            for row in datasets
        ],
    }
    (args.output_dir / "harmonization_manifest.json").write_text(json.dumps(release, indent=2) + "\n")


if __name__ == "__main__":
    main()
