#!/usr/bin/env python3
"""
Wrangle GWAS summary statistics to COJO format for SBayesRC.

COJO format columns:
    SNP   - rsID
    A1    - effect allele
    A2    - other allele
    freq  - frequency of A1
    b     - beta (or log(OR) for case-control)
    se    - standard error
    p     - p-value
    N     - sample size

Usage:
    python wrangle.py                     # Process all sources with raw files
    python wrangle.py --id edu_2022       # Process specific source
    python wrangle.py --list              # List all sources and status
"""

import argparse
import gzip
import math
from pathlib import Path

import polars as pl
import yaml

# Paths
SCRIPT_DIR = Path(__file__).parent
SOURCES_FILE = SCRIPT_DIR / "sources.yaml"
RAW_DIR = SCRIPT_DIR / "raw"
COJO_DIR = SCRIPT_DIR / "cojo"


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


# =============================================================================
# FORMAT PARSERS
# =============================================================================


def parse_standard(raw_path: Path, source: dict) -> pl.DataFrame:
    """
    Parse standard tab/space-delimited sumstats.
    Works for: EA4, IQ, Height, BMI, most direct downloads.
    """
    mapping = source["column_mapping"]
    n_samples = source.get("n_override") or source.get("n_samples")

    # Read file as all strings to handle whitespace/scientific notation issues
    null_values = ["NA", "na", "N/A", ".", ""]
    df = pl.read_csv(
        raw_path,
        separator="\t",
        infer_schema_length=0,  # Don't infer, read as strings
        null_values=null_values,
    )

    # Handle space-separated files
    if len(df.columns) == 1:
        df = pl.read_csv(raw_path, separator=" ", infer_schema_length=0, null_values=null_values)

    # Helper to safely cast: strip whitespace then convert to float
    def to_float(col_name: str) -> pl.Expr:
        return pl.col(col_name).str.strip_chars().cast(pl.Float64)

    # Build COJO dataframe
    cojo = df.select(
        pl.col(mapping["snp"]).alias("SNP"),
        pl.col(mapping["a1"]).alias("A1"),
        pl.col(mapping["a2"]).alias("A2"),
        to_float(mapping["freq"]).alias("freq"),
        to_float(mapping["beta"]).alias("b"),
        to_float(mapping["se"]).alias("se"),
        to_float(mapping["p"]).alias("p"),
    )

    # Handle beta conversion (OR -> log(OR))
    if source.get("beta_is_or"):
        cojo = cojo.with_columns(pl.col("b").log().alias("b"))

    # Add N column (from original df or use fixed value)
    if mapping.get("n") and mapping["n"] in df.columns:
        n_col = pl.col(mapping["n"]).str.strip_chars().cast(pl.Float64).alias("N")
        cojo = cojo.with_columns(df.select(n_col))
    else:
        cojo = cojo.with_columns(pl.lit(n_samples).cast(pl.Float64).alias("N"))

    return cojo


def parse_daner(raw_path: Path, source: dict) -> pl.DataFrame:
    """
    Parse PGC daner format.
    These have OR (needs log transform), frequency columns like FRQ_U_XXXXX.
    """
    mapping = source["column_mapping"]
    n_samples = source.get("n_override") or source.get("n_samples")

    df = pl.read_csv(raw_path, separator="\t", infer_schema_length=10000)

    # Handle space-separated daner files
    if len(df.columns) == 1:
        df = pl.read_csv(raw_path, separator=" ", infer_schema_length=10000)

    # Find frequency column (may have variable suffix)
    freq_col = mapping["freq"]
    if freq_col not in df.columns:
        # Try to find column starting with FRQ_U_ or FRQ_A_
        for col in df.columns:
            if col.startswith("FRQ_U_") or col.startswith("FRQ_A_"):
                freq_col = col
                break

    # Build COJO dataframe
    cojo = df.select(
        pl.col(mapping["snp"]).alias("SNP"),
        pl.col(mapping["a1"]).alias("A1"),
        pl.col(mapping["a2"]).alias("A2"),
        pl.col(freq_col).cast(pl.Float64).alias("freq"),
        pl.col(mapping["beta"]).cast(pl.Float64).alias("b"),
        pl.col(mapping["se"]).cast(pl.Float64).alias("se"),
        pl.col(mapping["p"]).cast(pl.Float64).alias("p"),
    )

    # Convert OR to log(OR)
    if source.get("beta_is_or"):
        cojo = cojo.with_columns(pl.col("b").log().alias("b"))

    # Add N column
    if mapping.get("n") and mapping["n"] in df.columns:
        cojo = cojo.with_columns(df.select(pl.col(mapping["n"]).cast(pl.Float64).alias("N")))
    else:
        cojo = cojo.with_columns(pl.lit(n_samples).cast(pl.Float64).alias("N"))

    return cojo


def parse_pgc_vcf(raw_path: Path, source: dict) -> pl.DataFrame:
    """
    Parse PGC VCF-style TSV format (like PGC3 SCZ, PTSD, MDD).
    Has ## comment lines and header starting with #CHROM.
    """
    import io

    mapping = source["column_mapping"]
    n_samples = source.get("n_override") or source.get("n_samples")

    # Read file, skipping ## header lines (polars skip_rows doesn't work well with gzip)
    opener = gzip.open if str(raw_path).endswith(".gz") else open
    lines = []
    with opener(raw_path, "rt") as f:
        for line in f:
            if not line.startswith("##"):
                lines.append(line)

    df = pl.read_csv(
        io.StringIO("".join(lines)),
        separator="\t",
        infer_schema_length=0,  # Read as strings to handle scientific notation
        null_values=["NA", "na", "N/A", ".", ""],
    )

    # Clean column names (remove # prefix)
    df = df.rename({col: col.lstrip("#") for col in df.columns})

    # Build COJO dataframe
    cojo = df.select(
        pl.col(mapping["snp"]).alias("SNP"),
        pl.col(mapping["a1"]).alias("A1"),
        pl.col(mapping["a2"]).alias("A2"),
        pl.col(mapping["freq"]).cast(pl.Float64).alias("freq"),
        pl.col(mapping["beta"]).cast(pl.Float64).alias("b"),
        pl.col(mapping["se"]).cast(pl.Float64).alias("se"),
        pl.col(mapping["p"]).cast(pl.Float64).alias("p"),
    )

    # Handle beta conversion if needed
    if source.get("beta_is_or"):
        cojo = cojo.with_columns(pl.col("b").log().alias("b"))

    # Add N column
    if mapping.get("n") and mapping["n"] in df.columns:
        cojo = cojo.with_columns(df.select(pl.col(mapping["n"]).cast(pl.Float64).alias("N")))
    else:
        cojo = cojo.with_columns(pl.lit(n_samples).cast(pl.Float64).alias("N"))

    return cojo


def parse_opengwas_vcf(raw_path: Path, source: dict) -> pl.DataFrame:
    """
    Parse OpenGWAS VCF format.
    These have INFO field with AF, and FORMAT fields with ES, SE, LP, etc.
    """
    import io

    mapping = source["column_mapping"]
    n_samples = source.get("n_override") or source.get("n_samples")

    # Read VCF, skipping ## header lines (polars skip_rows doesn't work with gzip)
    opener = gzip.open if str(raw_path).endswith(".gz") else open
    lines = []
    with opener(raw_path, "rt") as f:
        for line in f:
            if not line.startswith("##"):
                lines.append(line)

    df = pl.read_csv(
        io.StringIO("".join(lines)),
        separator="\t",
        infer_schema_length=0,  # Read as strings
        null_values=["NA", "na", "N/A", ".", ""],
    )

    # Clean column names
    df = df.rename({col: col.lstrip("#") for col in df.columns})

    # Get sample column name BEFORE adding new columns
    sample_col = df.columns[-1]  # Last column is sample data

    # Extract AF from INFO field
    df = df.with_columns(
        pl.col("INFO").str.extract(r"AF=([0-9.]+)").cast(pl.Float64).alias("AF_extracted")
    )

    # Parse FORMAT field to get column order
    # FORMAT is like "ES:SE:LP:AF:ID" and the sample column has values "0.01:0.002:-2.5:0.3:rs123"

    # Get format fields from first row
    format_fields = df["FORMAT"][0].split(":")

    # Split sample column by :
    for i, field in enumerate(format_fields):
        df = df.with_columns(
            pl.col(sample_col).str.split(":").list.get(i).alias(field + "_val")
        )

    # Build COJO dataframe using select
    cojo = df.select(
        pl.col("ID").alias("SNP"),
        pl.col("ALT").alias("A1"),
        pl.col("REF").alias("A2"),
        pl.col("AF_extracted").alias("freq"),
        pl.col("ES_val").cast(pl.Float64).alias("b"),
        pl.col("SE_val").cast(pl.Float64).alias("se"),
        (10.0 ** (-pl.col("LP_val").cast(pl.Float64))).alias("p"),
        pl.lit(n_samples).cast(pl.Float64).alias("N"),
    )

    return cojo


def parse_gwas_catalog_harmonized(raw_path: Path, source: dict) -> pl.DataFrame:
    """
    Parse GWAS Catalog harmonized format.
    Standard columns: variant_id, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value
    """
    mapping = source["column_mapping"]
    n_samples = source.get("n_override") or source.get("n_samples")

    df = pl.read_csv(raw_path, separator="\t", infer_schema_length=10000)

    # Build COJO dataframe
    cojo = df.select(
        pl.col(mapping["snp"]).alias("SNP"),
        pl.col(mapping["a1"]).alias("A1"),
        pl.col(mapping["a2"]).alias("A2"),
        pl.col(mapping["freq"]).cast(pl.Float64).alias("freq"),
        pl.col(mapping["beta"]).cast(pl.Float64).alias("b"),
        pl.col(mapping["se"]).cast(pl.Float64).alias("se"),
        pl.col(mapping["p"]).cast(pl.Float64).alias("p"),
    )

    # Handle beta conversion if needed
    if source.get("beta_is_or"):
        cojo = cojo.with_columns(pl.col("b").log().alias("b"))

    # Add N column
    if mapping.get("n") and mapping["n"] in df.columns:
        cojo = cojo.with_columns(df.select(pl.col(mapping["n"]).cast(pl.Float64).alias("N")))
    else:
        cojo = cojo.with_columns(pl.lit(n_samples).cast(pl.Float64).alias("N"))

    return cojo


# Format dispatcher
PARSERS = {
    "standard": parse_standard,
    "daner": parse_daner,
    "pgc_vcf": parse_pgc_vcf,
    "opengwas_vcf": parse_opengwas_vcf,
    "gwas_catalog_harmonized": parse_gwas_catalog_harmonized,
}


# =============================================================================
# MAIN FUNCTIONS
# =============================================================================


def wrangle_source(source: dict, dry_run: bool = False) -> bool:
    """
    Wrangle a single source to COJO format.
    Returns True if successful.
    """
    source_id = source["id"]
    raw_file = source.get("raw_file")
    output_file = source.get("output_file")
    fmt = source.get("format")

    if not raw_file:
        print(f"  [{source_id}] No raw_file specified, skipping")
        return False

    raw_path = RAW_DIR / raw_file
    if not raw_path.exists():
        print(f"  [{source_id}] Raw file not found: {raw_path}")
        return False

    if not fmt or fmt not in PARSERS:
        print(f"  [{source_id}] Unknown format: {fmt}")
        return False

    if dry_run:
        print(f"  [{source_id}] Would process {raw_file} -> {output_file}")
        return True

    print(f"  [{source_id}] Processing {raw_file}...")

    try:
        parser = PARSERS[fmt]
        cojo_df = parser(raw_path, source)

        # Uppercase alleles (GCTB requires uppercase)
        cojo_df = cojo_df.with_columns(
            pl.col("A1").str.to_uppercase(),
            pl.col("A2").str.to_uppercase(),
        )

        # Filter out rows with missing values
        cojo_df = cojo_df.drop_nulls()

        # Write output
        output_path = COJO_DIR / output_file
        cojo_df.write_csv(output_path, separator="\t")

        print(f"  [{source_id}] Wrote {len(cojo_df):,} SNPs to {output_path.name}")
        return True

    except Exception as e:
        print(f"  [{source_id}] ERROR: {e}")
        return False


def list_sources():
    """List all sources and their status."""
    sources = load_sources()

    print(f"\n{'ID':<20} {'Phenotype':<25} {'Format':<15} {'Status':<10} {'Raw File'}")
    print("-" * 100)

    for s in sources:
        source_id = s["id"]
        phenotype = s["phenotype"][:24]
        fmt = s.get("format") or "?"
        status = s.get("sbayesrc_status", "?")

        raw_file = s.get("raw_file") or "-"
        raw_path = RAW_DIR / raw_file if raw_file != "-" else None
        has_raw = "Y" if raw_path and raw_path.exists() else "N"

        print(f"{source_id:<20} {phenotype:<25} {fmt:<15} {status:<10} [{has_raw}] {raw_file}")


def main():
    parser = argparse.ArgumentParser(description="Wrangle GWAS sumstats to COJO format")
    parser.add_argument("--id", help="Process specific source ID")
    parser.add_argument("--list", action="store_true", help="List all sources")
    parser.add_argument("--all", action="store_true", help="Process all sources with raw files")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be done")
    args = parser.parse_args()

    if args.list:
        list_sources()
        return

    # Ensure output directory exists
    COJO_DIR.mkdir(exist_ok=True)

    if args.id:
        source = get_source(args.id)
        if not source:
            print(f"Source not found: {args.id}")
            return
        wrangle_source(source, dry_run=args.dry_run)

    elif args.all:
        sources = load_sources()
        success = 0
        for source in sources:
            if wrangle_source(source, dry_run=args.dry_run):
                success += 1
        print(f"\nProcessed {success}/{len(sources)} sources")

    else:
        parser.print_help()


if __name__ == "__main__":
    main()
