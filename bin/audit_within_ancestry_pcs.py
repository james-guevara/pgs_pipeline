#!/usr/bin/env python3
"""Audit within-ancestry PCs for cohort and technical covariates."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pcs", required=True, type=Path)
    parser.add_argument("--sample-qc", required=True, type=Path)
    parser.add_argument("--sample-manifest", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--manifest-id-columns", required=True)
    parser.add_argument("--categorical", default="source")
    parser.add_argument("--write-participant-table", action="store_true")
    return parser.parse_args()


def build_manifest_lookup(manifest: pd.DataFrame, id_columns: list[str]) -> pd.DataFrame:
    parts = []
    for column in id_columns:
        if column not in manifest:
            continue
        part = manifest.loc[manifest[column].notna()].copy()
        part["CURRENT_IID"] = part[column].astype(str)
        parts.append(part)
    if not parts:
        raise ValueError("none of the requested manifest identifier columns exist")
    lookup = pd.concat(parts, ignore_index=True)
    conflicts = lookup.groupby("CURRENT_IID").size()
    if (conflicts > 1).any():
        # Repeated aliases for the same manifest row are harmless; conflicting rows are not.
        check_columns = [column for column in manifest.columns if column not in id_columns]
        for identifier in conflicts[conflicts > 1].index:
            records = lookup.loc[lookup.CURRENT_IID == identifier, check_columns].drop_duplicates()
            if len(records) > 1:
                raise ValueError(f"manifest identifier maps to conflicting rows: {identifier}")
    return lookup.drop_duplicates("CURRENT_IID", keep="first")


def eta_squared(frame: pd.DataFrame, value: str, category: str) -> tuple[float, int, int]:
    subset = frame[[value, category]].dropna()
    groups = subset[category].nunique()
    if len(subset) < 2 or groups < 2:
        return np.nan, len(subset), groups
    grand_mean = subset[value].mean()
    between = sum(
        len(group) * (group[value].mean() - grand_mean) ** 2
        for _, group in subset.groupby(category)
    )
    total = ((subset[value] - grand_mean) ** 2).sum()
    return (between / total if total else np.nan), len(subset), groups


def scatter_by_category(frame: pd.DataFrame, x: str, y: str, category: str, path: Path) -> None:
    subset = frame[[x, y, category]].dropna()
    fig, axis = plt.subplots(figsize=(8, 6))
    for label, group in subset.groupby(category):
        axis.scatter(group[x], group[y], s=9, alpha=0.55, label=str(label))
    axis.set_xlabel(x)
    axis.set_ylabel(y)
    axis.set_title(f"{x} versus {y}, colored by {category}")
    axis.legend(fontsize=7, frameon=False, markerscale=1.5, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    pcs = pd.read_csv(args.pcs, sep="\t")
    qc = pd.read_csv(args.sample_qc, sep="\t")
    manifest = pd.read_csv(args.sample_manifest)
    pc_columns = [column for column in pcs if column.startswith("PC") and column.endswith("_AVG")]
    if not pc_columns or "IID" not in pcs or "IID" not in qc:
        raise ValueError("PC and sample-QC inputs must contain IID; PC input must contain PC*_AVG")

    lookup = build_manifest_lookup(manifest, args.manifest_id_columns.split(","))
    joined = pcs.merge(lookup, left_on="IID", right_on="CURRENT_IID", how="left")
    joined = joined.merge(qc[["IID", "source", "missingness"]], on="IID", how="left")
    categories = [column for column in args.categorical.split(",") if column in joined]

    summary = []
    means = []
    for category in categories:
        for pc in pc_columns:
            eta2, sample_count, category_count = eta_squared(joined, pc, category)
            summary.append({
                "covariate": category,
                "pc": pc.removesuffix("_AVG"),
                "n": sample_count,
                "categories": category_count,
                "eta_squared": eta2,
            })
        for label, group in joined[[category, *pc_columns]].dropna(subset=[category]).groupby(category):
            for pc in pc_columns:
                values = group[pc].dropna()
                means.append({
                    "covariate": category,
                    "level": label,
                    "pc": pc.removesuffix("_AVG"),
                    "n": len(values),
                    "mean": values.mean(),
                    "sd": values.std(),
                })

    correlations = []
    wgs = joined.loc[joined.source == "WGS"]
    for pc in pc_columns:
        subset = wgs[[pc, "missingness"]].dropna()
        correlations.append({
            "subset": "WGS",
            "metric": "sample_missingness",
            "pc": pc.removesuffix("_AVG"),
            "n": len(subset),
            "pearson_r": subset[pc].corr(subset.missingness),
        })

    pd.DataFrame(summary).to_csv(args.out_dir / "categorical_associations.tsv", sep="\t", index=False)
    pd.DataFrame(means).to_csv(args.out_dir / "category_pc_means.tsv", sep="\t", index=False)
    pd.DataFrame(correlations).to_csv(args.out_dir / "continuous_associations.tsv", sep="\t", index=False)
    if args.write_participant_table:
        joined[["IID", "source", *[c for c in categories if c != "source"], "missingness", *pc_columns]].to_csv(
            args.out_dir / "audit_participants.tsv", sep="\t", index=False
        )

    heatmap = pd.DataFrame(summary).pivot(index="covariate", columns="pc", values="eta_squared")
    heatmap = heatmap.reindex(columns=[pc.removesuffix("_AVG") for pc in pc_columns])
    fig, axis = plt.subplots(figsize=(11, max(3, len(heatmap) * 0.75)))
    image = axis.imshow(heatmap, aspect="auto", cmap="viridis", vmin=0)
    axis.set_xticks(range(len(heatmap.columns)), heatmap.columns)
    axis.set_yticks(range(len(heatmap.index)), heatmap.index)
    axis.set_title("Variance in each PC explained by categorical covariates (unadjusted η²)")
    fig.colorbar(image, ax=axis, label="η²")
    fig.tight_layout()
    fig.savefig(args.out_dir / "categorical_associations.png", dpi=180)
    plt.close(fig)

    for category in categories:
        scatter_by_category(joined, pc_columns[0], pc_columns[1], category, args.out_dir / f"pc1_pc2_by_{category}.png")
        scatter_by_category(joined, pc_columns[1], pc_columns[2], category, args.out_dir / f"pc2_pc3_by_{category}.png")

    print(f"participants={len(joined)} manifest_matched={joined.CURRENT_IID.notna().sum()}")
    print(f"output={args.out_dir}")


if __name__ == "__main__":
    main()
