#!/usr/bin/env python3
"""Verify selected reference artifacts against a SHA256SUMS file."""

import argparse
import hashlib
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checksums", type=Path, required=True)
    parser.add_argument("--artifact", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    expected = {}
    for line in args.checksums.read_text().splitlines():
        digest, filename = line.split(None, 1)
        expected[Path(filename.strip()).name] = digest

    rows = []
    failed = False
    for artifact in args.artifact:
        observed = sha256(artifact)
        wanted = expected.get(artifact.name)
        status = "valid" if wanted == observed else ("not_listed" if wanted is None else "mismatch")
        failed = failed or status != "valid"
        rows.append((artifact.name, wanted or "", observed, status))

    with args.output.open("w") as handle:
        handle.write("artifact\texpected_sha256\tobserved_sha256\tstatus\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")
    if failed:
        raise SystemExit("PCA reference checksum validation failed")


if __name__ == "__main__":
    main()
