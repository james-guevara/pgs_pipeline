#!/usr/bin/env python3
"""Build a provisional catalog from the score collection on Expanse."""

from __future__ import annotations

import csv
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCORES = ROOT / "resources" / "expanse" / "scores.tsv"
OUTPUT = ROOT / "resources" / "score_catalog.tsv"
EXPANSE_MANIFEST = ROOT / "resources" / "expanse" / "scores.tsv"
REMOTE_ROOT = "/expanse/projects/sebat1/a1sriniv/pgs_sumstats_organized"


def remote_files(subdir: str) -> list[str]:
    command = (
        f"find -L {REMOTE_ROOT}/{subdir} -maxdepth 1 -type f "
        "-printf '%f\\n' 2>/dev/null | sort"
    )
    result = subprocess.run(
        ["ssh", "expanse", command], check=True, capture_output=True, text=True
    )
    return [line for line in result.stdout.splitlines() if line]


def enabled_weights() -> dict[str, str]:
    with SCORES.open(newline="") as handle:
        return {
            Path(row["weights"]).name: row["trait"]
            for row in csv.DictReader(handle, delimiter="\t")
        }


def source_metadata() -> dict[str, dict[str, str]]:
    """Read the simple top-level fields from sources.yaml without a YAML dependency."""
    source_file = ROOT / "sumstats" / "sources.yaml"
    records: dict[str, dict[str, str]] = {}
    current: dict[str, str] | None = None
    for line in source_file.read_text().splitlines():
        id_match = re.match(r"^  - id:\s*(.+?)\s*$", line)
        if id_match:
            current = {"id": id_match.group(1).strip('"\'')}
            records[current["id"]] = current
            continue
        if current is None:
            continue
        field_match = re.match(
            r"^    (phenotype|publication|doi|year|ancestry):\s*(.*?)\s*(?:#.*)?$",
            line,
        )
        if field_match:
            value = field_match.group(2).strip().strip('"\'')
            current[field_match.group(1)] = "" if value in {"null", "None"} else value
    return records


def strip_extensions(name: str) -> str:
    return re.sub(r"(?:\.txt|\.tsv|\.tbl|\.ma|\.vcf|\.gz)+$", "", name)


def trait_id(weight: str) -> str:
    stem = strip_extensions(weight)
    stem = re.sub(r"_cojo(?:_freq)?$", "", stem, flags=re.I)
    stem = re.sub(r"\.bgen\.stats$", "", stem, flags=re.I)
    stem = re.sub(r"[^A-Za-z0-9]+", "_", stem).strip("_").lower()
    return stem


def category(identifier: str) -> str:
    rules = [
        ("growth_obesity_development", r"obes|bmi|height|birth|gest|preterm|postterm|pubertal|growth|headcircum"),
        ("psychiatric_neurologic", r"adhd|asd|bip|mdd|scz|ptsd|epilep|focal|neurodevelopmental|internalizing|compulsive|substanceuse|pfactor|anorexia"),
        ("sleep_circadian", r"sleep|insomnia|dozing|gettingup|morning|napping|snoring|nightmare|circadian|chronotype|alertness|noninsomnia"),
        ("substance_behavior", r"smok|cig|drink|cannabis|alcohol"),
        ("education_cognition", r"edu|iq|cog|intelligence"),
        ("environment", r"pm25|no2"),
        ("diet", r"carb|fat|protein|breakfast|macronutrient"),
    ]
    compact = identifier.replace("_", "")
    for label, pattern in rules:
        if re.search(pattern, compact, flags=re.I):
            return label
    return "other"


DISPLAY_REPLACEMENTS = {
    "adhd": "ADHD", "asd": "ASD", "bip": "Bipolar disorder",
    "bmi": "BMI", "cigperday": "Cigarettes per day", "edu": "Educational attainment",
    "iq": "Intelligence", "mdd": "Major depressive disorder", "no2": "NO2 exposure",
    "pm25": "PM2.5 exposure", "ptsd": "PTSD", "scz": "Schizophrenia",
    "smkces": "Smoking cessation", "afb": "Age at first birth",
}


STUDY_RULES = [
    {
        "pattern": r"^childhood_obesity_",
        "publication": "Bradfield et al. A trans-ancestral meta-analysis of genome-wide association studies reveals loci associated with childhood obesity",
        "year": "2019", "ancestry": "trans-ancestry",
        "doi": "10.1093/hmg/ddz161",
    },
    {
        "pattern": r"^(f[1-5]_|pfactor_2025$)",
        "publication": "Mapping the genetic landscape across 14 psychiatric disorders",
        "year": "2025", "ancestry": "EUR-like",
        "doi": "10.1038/s41586-025-09820-3",
    },
    {
        "pattern": r"^fetal_(early_preterm_birth|gest_duration|postterm_birth|preterm_birth)_ncomms2019$",
        "publication": "Variants in the fetal genome near pro-inflammatory cytokine genes on 2q13 associate with gestational duration",
        "year": "2019", "ancestry": "",
        "doi": "10.1038/s41467-019-11881-8",
    },
    {
        "pattern": r"^gcst9001187[45]_",
        "publication": "Demange et al. Investigating the genetic architecture of noncognitive skills using GWAS-by-subtraction",
        "year": "2021", "ancestry": "predominantly European",
        "doi": "10.1038/s41588-020-00754-2",
    },
    {
        "pattern": r"^(dozing_sumstats_jansenetal|gettingup_sumstats_jansenetal|morningness_sumstats_jansenetal|napping_sumstats_jansenetal|sleepdur_sumstats_jansenetal|snoring_sumstats_jansenetal|insomnia_sumstats_jansenetal_2019_main)$",
        "publication": "Jansen et al. Genome-wide analysis of insomnia and related sleep phenotypes",
        "year": "2019", "ancestry": "European",
        "doi": "10.1038/s41588-018-0333-3",
    },
    {
        "pattern": r"^(epilepsy_full|focal_epilepsy|genetic_generalized_epilepsy)$",
        "publication": "International League Against Epilepsy Consortium. GWAS meta-analysis of over 29,000 people with epilepsy identifies 26 risk loci and subtype-specific genetic architecture",
        "year": "2023", "ancestry": "multi-ancestry",
        "doi": "10.1038/s41588-023-01485-w",
    },
    {
        "pattern": r"^guan_2025_",
        "publication": "Guan et al. Family-based genome-wide association study designs for increased power and robustness",
        "year": "2025", "ancestry": "",
        "doi": "10.1038/s41588-025-02118-0",
    },
    {
        "pattern": r"^bolt_ss_.*_euro_hrc_1kg$",
        "publication": "Goodman et al. Genome-wide association analysis of composite sleep health scores in 413,904 individuals",
        "year": "2025", "ancestry": "European",
        "doi": "10.1038/s42003-025-07514-0",
    },
    {
        "pattern": r"^chronotype_raw_bolt_",
        "publication": "Jones et al. Genome-wide association analyses of chronotype in 697,828 individuals provides insights into circadian rhythms",
        "year": "2019", "ancestry": "European",
        "doi": "10.1038/s41467-018-08259-7",
    },
    {
        "pattern": r"^breakfast_skipping_2019_dashti$",
        "publication": "Dashti et al. Genome-wide association study of breakfast skipping links clock regulation with food timing",
        "year": "2019", "ancestry": "European",
        "doi": "10.1093/ajcn/nqz076",
    },
    {
        "pattern": r"^napping_2021_saxena$",
        "publication": "Dashti et al. Genetic determinants of daytime napping and effects on cardiometabolic health",
        "year": "2021", "ancestry": "European",
        "doi": "10.1038/s41467-020-20585-3",
    },
    {
        "pattern": r"^maternal_gestation_duration_",
        "publication": "Sole-Navais et al. Genetic effects on the timing of parturition and links to fetal birth weight",
        "year": "2023", "ancestry": "multi-ancestry",
        "doi": "10.1038/s41588-023-01343-9",
    },
    {
        "pattern": r"^(pubertal_growth_|sitar_longitudinal_)",
        "publication": "Bradfield et al. Trans-ancestral genome-wide association study of longitudinal pubertal height growth and shared heritability with adult health outcomes",
        "year": "2024", "ancestry": "trans-ancestry",
        "doi": "10.1186/s13059-023-03136-z",
    },
    {
        "pattern": r"^headcircum_",
        "publication": "Haworth et al. Genetics of early-life head circumference and shared genetic factors with adult brain traits",
        "year": "2022", "ancestry": "",
        "doi": "10.1186/s12920-022-01281-1",
    },
    {
        "pattern": r"^meta_(carb|fat|pro)_charge_ukbb",
        "publication": "Merino et al. Genetic analysis of dietary intake identifies new loci and functional links with metabolic traits",
        "year": "2022", "ancestry": "European",
        "doi": "10.1038/s41562-021-01182-w",
    },
    {
        "pattern": r"^birth_weight_maternal_2018$",
        "publication": "Beaumont et al. Genome-wide association study of offspring birth weight in 86,577 women identifies five novel loci and highlights maternal genetic effects",
        "year": "2018", "ancestry": "European",
        "doi": "10.1093/hmg/ddx429",
    },
    {
        "pattern": r"^birth_weight_fetal_2016$",
        "publication": "Horikoshi et al. Genome-wide associations for birth weight and correlations with adult disease",
        "year": "2016", "ancestry": "predominantly European",
        "doi": "10.1038/nature19806",
    },
    {
        "pattern": r"^pgcan2_2019_07_anorexia_vcf$",
        "publication": "Watson et al. Genome-wide association study identifies eight risk loci and implicates metabo-psychiatric origins for anorexia nervosa",
        "year": "2019", "ancestry": "European",
        "doi": "10.1038/s41588-019-0439-2",
    },
    {
        "pattern": r"^gcst90002409_childhood_bmi$",
        "publication": "Felix et al. Genome-wide association analysis identifies three new susceptibility loci for childhood body mass index",
        "year": "2016", "ancestry": "European",
        "doi": "10.1093/hmg/ddv472",
    },
    {
        "pattern": r"^birth_length_2014$",
        "publication": "van der Valk et al. A novel common variant in DCST2 is associated with length in early life and height in adulthood",
        "year": "2014", "ancestry": "European",
        "doi": "10.1093/hmg/ddu510",
    },
    {
        "pattern": r"^aowgwas_sumstats_gui2025$",
        "publication": "Gui et al. Genome-wide association meta-analysis of age at onset of walking in over 70,000 infants of European ancestry",
        "year": "2025", "ancestry": "European",
        "doi": "10.1038/s41562-025-02145-1",
    },
    {
        "pattern": r"^metaanalysis_nightmares_",
        "publication": "Ollila et al. Nightmares share genetic risk factors with sleep and psychiatric traits",
        "year": "2024", "ancestry": "predominantly European",
        "doi": "10.1038/s41398-023-02637-6",
    },
]


def family_metadata(identifier: str) -> dict[str, str]:
    for rule in STUDY_RULES:
        if re.search(rule["pattern"], identifier):
            return rule
    return {}


def display_name(identifier: str) -> str:
    if identifier in DISPLAY_REPLACEMENTS:
        return DISPLAY_REPLACEMENTS[identifier]
    text = identifier.replace("_", " ")
    text = re.sub(r"\bcojo\b", "", text, flags=re.I)
    text = re.sub(r"\s+", " ", text).strip()
    return text.title().replace("Bmi", "BMI").replace("Adhd", "ADHD").replace("Asd", "ASD")


def raw_match(weight: str, raw_files: list[str]) -> tuple[str, str]:
    key = strip_extensions(weight)
    key = re.sub(r"_cojo(?:_freq)?$", "", key, flags=re.I)
    key = strip_extensions(key)
    exact = [raw for raw in raw_files if strip_extensions(raw).lower() == key.lower()]
    if len(exact) == 1:
        return exact[0], "exact_filename_stem"
    prefix = [raw for raw in raw_files if strip_extensions(raw).lower().startswith(key.lower())]
    if len(prefix) == 1:
        return prefix[0], "filename_prefix"
    return "", "unresolved"


def main() -> None:
    weights = remote_files("processed/txt_weights")
    raw = remote_files("raw_sumstats")
    enabled = enabled_weights()
    metadata = source_metadata()
    rows = []
    seen: set[str] = set()
    for weight in weights:
        identifier = enabled.get(weight, trait_id(weight))
        if identifier in seen:
            raise ValueError(f"Duplicate generated trait_id: {identifier}")
        seen.add(identifier)
        raw_file, match_method = raw_match(weight, raw)
        metadata_id = {"asd_loSSC": "asd_2019"}.get(identifier, identifier)
        source = metadata.get(metadata_id, {})
        family = family_metadata(identifier)
        rows.append({
            "trait_id": identifier,
            "phenotype_label": source.get("phenotype", display_name(identifier)),
            "category": category(identifier),
            "weight_file": weight,
            "weight_path_expanse": f"{REMOTE_ROOT}/processed/txt_weights/{weight}",
            "raw_source_file": raw_file,
            "raw_match_method": match_method,
            "currently_enabled": "yes" if weight in enabled else "no",
            "processing_status": "final_weight_available",
            "catalog_status": "provisional",
            "ancestry": source.get("ancestry", family.get("ancestry", "")),
            "genome_build": "not_encoded_rsid_only",
            "publication": source.get("publication", family.get("publication", "")),
            "year": source.get("year", family.get("year", "")),
            "doi_or_source_url": source.get("doi", family.get("doi", "")),
            "notes": (
                "Imported from sumstats/sources.yaml; verify before release" if source
                else "Matched to a verified study family; remaining fields require verification" if family
                else "Metadata fields require source verification"
            ),
        })
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=rows[0].keys(), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    with EXPANSE_MANIFEST.open("w", newline="") as handle:
        fields = ["trait", "weights", "id_col", "allele_col", "effect_col"]
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({
                "trait": row["trait_id"],
                "weights": row["weight_path_expanse"],
                "id_col": 1,
                "allele_col": 2,
                "effect_col": 3,
            })
    print(f"Wrote {len(rows)} rows to {OUTPUT}")
    print(f"Wrote {len(rows)} rows to {EXPANSE_MANIFEST}")


if __name__ == "__main__":
    main()
