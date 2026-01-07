#!/usr/bin/env python3
"""
Download GWAS summary statistics from various sources.

This script helps manage downloads and documents how to obtain each file.
Some sources require manual download (clicking through portals), while others
can be fetched directly via wget/curl.

Usage:
    python download.py --list              # Show download status for all sources
    python download.py --id edu_2022       # Show download instructions for specific source
    python download.py --id edu_2022 --fetch   # Actually download (if wget-able)
    python download.py --wget-all          # Download all wget-able sources
"""

import argparse
import subprocess
from pathlib import Path

import yaml

SCRIPT_DIR = Path(__file__).parent
SOURCES_FILE = SCRIPT_DIR / "sources.yaml"
RAW_DIR = SCRIPT_DIR / "raw"


def load_sources() -> list[dict]:
    """Load sources from YAML config."""
    with open(SOURCES_FILE) as f:
        config = yaml.safe_load(f)
    return config["sources"]


def get_source(source_id: str) -> dict | None:
    """Get a specific source by ID."""
    sources = load_sources()
    for s in sources:
        if s["id"] == source_id:
            return s
    return None


def check_raw_file(source: dict) -> bool:
    """Check if raw file exists."""
    raw_file = source.get("raw_file")
    if not raw_file:
        return False
    return (RAW_DIR / raw_file).exists()


def show_download_info(source: dict):
    """Show download information for a source."""
    source_id = source["id"]
    phenotype = source["phenotype"]
    download = source.get("download", {})
    method = download.get("method", "unknown")
    url = download.get("url")
    command = download.get("command")
    notes = download.get("notes")
    raw_file = source.get("raw_file")
    has_file = check_raw_file(source)

    print(f"\n{'='*60}")
    print(f"Source: {source_id}")
    print(f"Phenotype: {phenotype}")
    print(f"Raw file: {raw_file} {'[EXISTS]' if has_file else '[MISSING]'}")
    print(f"Download method: {method}")
    print(f"{'='*60}")

    if url:
        print(f"\nURL:\n  {url}")

    if command:
        print(f"\nDownload command:\n  cd {RAW_DIR}\n  {command}")

    if notes:
        print(f"\nNotes:\n  {notes}")

    if method == "manual":
        print("\n** This source requires manual download. **")
        print("   1. Visit the URL above")
        print("   2. Follow instructions to download")
        print(f"   3. Save file as: {RAW_DIR}/{raw_file}")


def fetch_source(source: dict) -> bool:
    """Download a source if it has a wget command."""
    source_id = source["id"]
    download = source.get("download", {})
    method = download.get("method")
    command = download.get("command")
    raw_file = source.get("raw_file")

    if check_raw_file(source):
        print(f"[{source_id}] Raw file already exists: {raw_file}")
        return True

    if method != "wget" or not command:
        print(f"[{source_id}] Not wget-able (method={method})")
        return False

    print(f"[{source_id}] Downloading...")
    try:
        result = subprocess.run(
            command,
            shell=True,
            cwd=RAW_DIR,
            capture_output=True,
            text=True,
        )
        if result.returncode == 0:
            print(f"[{source_id}] Downloaded successfully")
            return True
        else:
            print(f"[{source_id}] Download failed: {result.stderr}")
            return False
    except Exception as e:
        print(f"[{source_id}] Error: {e}")
        return False


def list_downloads():
    """List all sources and their download status."""
    sources = load_sources()

    print(f"\n{'ID':<20} {'Phenotype':<20} {'Method':<10} {'Has File':<10} {'URL/Notes'}")
    print("-" * 100)

    for s in sources:
        source_id = s["id"]
        phenotype = s["phenotype"][:19]
        download = s.get("download", {})
        method = download.get("method", "?")
        has_file = "Y" if check_raw_file(s) else "N"
        url = download.get("url", "")[:40] if download.get("url") else "-"

        print(f"{source_id:<20} {phenotype:<20} {method:<10} {has_file:<10} {url}")

    # Summary
    total = len(sources)
    with_files = sum(1 for s in sources if check_raw_file(s))
    wgetable = sum(1 for s in sources if s.get("download", {}).get("method") == "wget")

    print(f"\nSummary: {with_files}/{total} have raw files, {wgetable} are wget-able")


def main():
    parser = argparse.ArgumentParser(description="Download GWAS summary statistics")
    parser.add_argument("--id", help="Show/fetch specific source")
    parser.add_argument("--list", action="store_true", help="List all sources and download status")
    parser.add_argument("--fetch", action="store_true", help="Actually download (with --id)")
    parser.add_argument("--wget-all", action="store_true", help="Download all wget-able sources")
    args = parser.parse_args()

    # Ensure raw directory exists
    RAW_DIR.mkdir(exist_ok=True)

    if args.list:
        list_downloads()
        return

    if args.id:
        source = get_source(args.id)
        if not source:
            print(f"Source not found: {args.id}")
            return
        if args.fetch:
            fetch_source(source)
        else:
            show_download_info(source)
        return

    if args.wget_all:
        sources = load_sources()
        success = 0
        for source in sources:
            if source.get("download", {}).get("method") == "wget":
                if fetch_source(source):
                    success += 1
        print(f"\nDownloaded {success} files")
        return

    parser.print_help()


if __name__ == "__main__":
    main()
