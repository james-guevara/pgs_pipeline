#!/usr/bin/env python3
"""Intersect exact SNP alleles and append samples from multiple PGEN filesets."""

import argparse
from contextlib import ExitStack, contextmanager
import gc
import json
import os
from pathlib import Path
import tempfile

import numpy as np
import pgenlib as pg


@contextmanager
def pgen_writer(*args, **kwargs):
    writer = pg.PgenWriter(*args, **kwargs)
    try:
        yield writer
    except BaseException:
        # pgenlib complains about incomplete output on close; retain the cause.
        output_path = Path(os.fsdecode(args[0]))
        try:
            writer.close()
        except RuntimeError:
            pass
        # Incomplete PGENs are never valid outputs. Explicit unlinking also
        # avoids delayed shared-filesystem cleanup failures when the enclosing
        # temporary directory is removed immediately after an expected error.
        del writer
        gc.collect()
        output_path.unlink(missing_ok=True)
        raise
    else:
        writer.close()


def read_table(path, marker):
    headers = []
    rows = []
    columns = None
    with path.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            if line.startswith("##"):
                headers.append(line.rstrip())
                continue
            if columns is None:
                columns = line.lstrip("#").split()
                if not line.startswith("#") or marker not in columns:
                    raise ValueError(f"{path}: expected a column header containing {marker}")
                if len(set(columns)) != len(columns):
                    raise ValueError(f"{path}: duplicate column names")
                continue
            values = line.split()
            if len(values) != len(columns):
                raise ValueError(f"{path}: row has wrong number of columns")
            rows.append(dict(zip(columns, values)))
    if not rows:
        raise ValueError(f"{path}: no records")
    return headers, columns, rows


def variants(prefix):
    headers, columns, rows = read_table(Path(str(prefix) + ".pvar"), "CHROM")
    if not {"CHROM", "POS", "ID", "REF", "ALT"}.issubset(columns):
        raise ValueError(f"{prefix}: incomplete PVAR columns")
    indexed = {}
    for i, row in enumerate(rows):
        ref, alt = row["REF"], row["ALT"]
        if ref not in ("A", "C", "G", "T") or alt not in ("A", "C", "G", "T") or ref == alt:
            raise ValueError(f"{prefix}: only biallelic A/C/G/T SNPs are supported (row {i + 1})")
        pos = int(row["POS"])
        if pos <= 0:
            raise ValueError(f"{prefix}: invalid position at row {i + 1}")
        key = (row["CHROM"], str(pos), ref, alt)
        if key in indexed:
            raise ValueError(f"{prefix}: duplicate variant key {key}")
        indexed[key] = i
    return headers, indexed


def samples(prefix):
    _, columns, rows = read_table(Path(str(prefix) + ".psam"), "IID")
    ids = set()
    for row in rows:
        key = (row.get("FID", "0"), row["IID"])
        if key in ids:
            raise ValueError(f"{prefix}: duplicate sample ID {key}")
        ids.add(key)
    return columns, rows, ids


def merge(first, second, output, dosage_mode="reject"):
    """Backward-compatible two-input API."""
    return merge_many([first, second], output, dosage_mode=dosage_mode)


def merge_many(prefixes, output, dosage_mode="reject"):
    """Intersect all inputs; publish a new output directory only after success."""
    if dosage_mode not in {"reject", "preserve"}:
        raise ValueError("dosage_mode must be reject or preserve")
    prefixes, output = [Path(prefix) for prefix in prefixes], Path(output)
    if len(prefixes) < 2:
        raise ValueError("At least two input filesets are required")
    if output.exists():
        raise ValueError(f"Output directory already exists: {output}")
    inputs, seen_ids, references = [], set(), set()
    columns = None
    offset = 0
    for prefix in prefixes:
        headers, index = variants(prefix)
        current_columns, rows, ids = samples(prefix)
        if columns is not None and current_columns != columns:
            raise ValueError(f"{prefix}: PSAM column names/order must match; harmonize the sample metadata first")
        columns = current_columns
        if seen_ids & ids:
            raise ValueError(f"{prefix}: Overlapping sample IDs: {len(seen_ids & ids)}")
        seen_ids.update(ids)
        current_refs = {h for h in headers if h.startswith("##reference=")}
        if references and current_refs and references != current_refs:
            raise ValueError(f"{prefix}: PVAR reference metadata differs")
        references.update(current_refs)
        inputs.append((prefix, index, rows, offset, offset + len(rows)))
        offset += len(rows)
    keys = list(inputs[0][1])
    for _, index, _, _, _ in inputs[1:]:
        keys = [key for key in keys if key in index]
    if not keys:
        raise ValueError("No matching variants across all inputs")
    input_summaries = [
        {"prefix": str(prefix.resolve()), "samples": len(rows),
         "variants": len(index), "unmatched": len(index) - len(keys)}
        for prefix, index, rows, _, _ in inputs
    ]
    summary = {
        "inputs": input_summaries,
        "input_count": len(inputs),
        "output_samples": offset, "output_variants": len(keys),
        "matching": "exact CHROM, POS, REF, ALT; intersection; first-input order",
        "dosage_policy": dosage_mode,
        "reference_status": "conservatively marked provisional",
    }
    # Retain the original summary fields for existing two-input callers.
    if len(inputs) == 2:
        for label, entry in zip(("first", "second"), input_summaries):
            summary.update({f"{label}_{field}": value for field, value in entry.items()})
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".pgen-merge-", dir=output.parent) as tmp:
        staging = Path(tmp) / "result"
        staging.mkdir()
        dest = staging / "merged"
        with ExitStack() as stack:
            readers = [stack.enter_context(pg.PgenReader(
                os.fsencode(str(prefix) + ".pgen"),
                raw_sample_ct=len(rows), variant_ct=len(index)))
                for prefix, index, rows, _, _ in inputs]
            phased = dosage_mode == "reject" and any(
                reader.hardcall_phase_present() for reader in readers
            )
            summary["phase_preserved"] = bool(phased)
            n = offset
            hardcalls = np.empty(n, dtype=np.int8)
            dosage_buffer = np.empty(max(len(item[2]) for item in inputs), dtype=np.float64)
            alleles = np.empty(2 * n, dtype=np.int32) if phased else None
            phase = np.empty(n, dtype=np.uint8) if phased else None
            merged_dosages = np.empty(n, dtype=np.float64) if dosage_mode == "preserve" else None
            with pgen_writer(
                os.fsencode(str(dest) + ".pgen"), n,
                variant_ct=len(keys), nonref_flags=True,
                hardcall_phase_present=phased,
                dosage_present=(dosage_mode == "preserve"),
            ) as writer:
                for key in keys:
                    for reader, (prefix, index, rows, start, end) in zip(readers, inputs):
                        idx = index[key]
                        dosage = dosage_buffer[:end - start]
                        reader.read(idx, hardcalls[start:end])
                        reader.read_dosages(idx, dosage)
                        if dosage_mode == "reject" and not np.array_equal(dosage, hardcalls[start:end]):
                            raise ValueError(f"{prefix}: Dosage differs from hardcall at {':'.join(key)}; hardcall-only merge refused")
                        if dosage_mode == "preserve":
                            merged_dosages[start:end] = dosage
                        if phased:
                            reader.read_alleles_and_phasepresent(
                                idx, alleles[2 * start:2 * end], phase[start:end])
                    if dosage_mode == "preserve":
                        writer.append_dosages(merged_dosages)
                    elif phased:
                        writer.append_partially_phased(alleles, phase)
                    else:
                        writer.append_biallelic(hardcalls)
        with Path(str(dest) + ".pvar").open("w") as handle:
            # Cohort-specific INFO/QUAL/FILTER values cannot describe all cohorts.
            for header in sorted(references):
                handle.write(header + "\n")
            handle.write("#CHROM\tPOS\tID\tREF\tALT\n")
            for chrom, pos, ref, alt in keys:
                handle.write(f"{chrom}\t{pos}\t{chrom}:{pos}:{ref}:{alt}\t{ref}\t{alt}\n")
        with Path(str(dest) + ".psam").open("w") as handle:
            handle.write("#" + "\t".join(columns) + "\n")
            for _, _, rows, _, _ in inputs:
                for row in rows:
                    handle.write("\t".join(row[col] for col in columns) + "\n")
        (staging / "merge_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
        # The existing-directory check above and rename avoid partial filesets.
        if output.exists():
            raise ValueError(f"Output directory appeared during merge: {output}")
        staging.rename(output)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", nargs="+", help="Two or more fileset prefixes, in sample output order")
    parser.add_argument("--first", help="First .pgen/.pvar/.psam prefix")
    parser.add_argument("--second", help="Second prefix; disjoint samples")
    parser.add_argument("--out-dir", required=True, help="New directory; contains merged.* and merge_summary.json")
    parser.add_argument(
        "--dosage-mode", choices=("reject", "preserve"), default="reject",
        help="Reject non-hardcall dosages (default), or preserve them in the merged PGEN",
    )
    args = parser.parse_args()
    if args.inputs is not None:
        if args.first is not None or args.second is not None:
            parser.error("Use --inputs or --first/--second, not both")
        if len(args.inputs) < 2:
            parser.error("--inputs requires at least two fileset prefixes")
        prefixes = args.inputs
    else:
        if args.first is None or args.second is None:
            parser.error("Provide --inputs with at least two prefixes, or both --first and --second")
        prefixes = [args.first, args.second]
    try:
        summary = merge_many(prefixes, args.out_dir, dosage_mode=args.dosage_mode)
    except (ValueError, RuntimeError, OSError) as exc:
        parser.exit(1, f"Merge failed: {exc}\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
