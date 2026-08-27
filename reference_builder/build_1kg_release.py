#!/usr/bin/env python3
"""Package an existing PLINK2 1000 Genomes PCA and train an ancestry model."""

import argparse
import hashlib
import json
import platform
import shutil
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import sklearn
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.metrics import accuracy_score, confusion_matrix, log_loss


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_pvar(path):
    with path.open() as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                columns = line.lstrip("#").rstrip().split("\t")
                break
        else:
            raise ValueError(f"No #CHROM header in {path}")
        rows = [line.rstrip().split("\t") for line in handle if line.strip()]
    frame = pd.DataFrame(rows, columns=columns)
    required = ["CHROM", "POS", "ID", "REF", "ALT"]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"PVAR missing columns: {', '.join(missing)}")
    return frame[required]


def new_classifier():
    return ExtraTreesClassifier(
        n_estimators=500,
        min_samples_leaf=2,
        class_weight="balanced",
        random_state=20260827,
        n_jobs=-1,
    )


def classifier_json(model, pc_columns, threshold):
    return {
        "model_type": "sklearn_extra_trees",
        "model_file": "classifier.joblib",
        "pc_columns": pc_columns,
        "classes": model.classes_.tolist(),
        "n_estimators": model.n_estimators,
        "min_samples_leaf": model.min_samples_leaf,
        "class_weight": model.class_weight,
        "random_state": model.random_state,
        "scikit_learn_version": sklearn.__version__,
        "assignment_probability_threshold": threshold,
        "uncertain_label": "ADMIXED_OR_UNCERTAIN",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--sample-metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reference-id", default="1kg_grch38_v1")
    parser.add_argument("--probability-threshold", type=float, default=0.8)
    args = parser.parse_args()

    source = args.source.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)

    source_files = {
        "allele_frequencies": source / "06d_pca/ref.acount",
        "eigenvalues": source / "06d_pca/ref.eigenval",
        "loadings": source / "06d_pca/ref.eigenvec.allele",
        "reference_eigenvectors": source / "06d_pca/ref.eigenvec",
        "reference_scores": source / "06d_pca/cohort.sscore",
        "panel_pvar": source / "06b_prune/cohort_pruned.pvar",
        "unrelated_ids": source / "06c_relatedness/unrelateds.king.cutoff.in.id",
    }
    for name, path in source_files.items():
        if not path.is_file():
            raise FileNotFoundError(f"Missing {name}: {path}")

    copied = {}
    for name in (
        "allele_frequencies",
        "eigenvalues",
        "loadings",
        "reference_eigenvectors",
        "reference_scores",
        "unrelated_ids",
    ):
        destination = output / source_files[name].name
        shutil.copy2(source_files[name], destination)
        copied[name] = destination

    panel = read_pvar(source_files["panel_pvar"])
    panel["AMBIGUOUS"] = (
        panel["REF"] + panel["ALT"]
    ).isin(["AT", "TA", "CG", "GC"])
    panel_path = output / "panel_variants.tsv"
    panel.to_csv(panel_path, sep="\t", index=False)
    copied["panel_variants"] = panel_path

    scores = pd.read_csv(source_files["reference_scores"], sep=r"\s+")
    metadata = pd.read_csv(args.sample_metadata, sep=r"\s+")
    labels = metadata[["SampleID", "Population", "Superpopulation"]].rename(
        columns={"SampleID": "IID"}
    )
    joined = scores.merge(labels, left_on="#IID", right_on="IID", validate="one_to_one")
    if len(joined) != len(scores):
        raise ValueError("Not every projected reference sample has metadata")

    pc_columns = [column for column in scores.columns if column.startswith("PC") and column.endswith("_AVG")]
    if not pc_columns:
        raise ValueError("No projected PC columns found")
    X = joined[pc_columns].to_numpy()
    y = joined["Superpopulation"].to_numpy()
    populations = joined["Population"].to_numpy()
    classes = sorted(np.unique(y).tolist())

    predictions = np.empty(len(joined), dtype=object)
    probabilities = np.zeros((len(joined), len(classes)))
    for population in sorted(np.unique(populations)):
        test = populations == population
        train = ~test
        model = new_classifier()
        model.fit(X[train], y[train])
        predictions[test] = model.predict(X[test])
        fold_probabilities = model.predict_proba(X[test])
        for index, label in enumerate(model.classes_):
            probabilities[test, classes.index(label)] = fold_probabilities[:, index]

    validation = joined[["#IID", "Population", "Superpopulation"]].copy()
    validation["PREDICTED"] = predictions
    for index, label in enumerate(classes):
        validation[f"PROB_{label}"] = probabilities[:, index]
    validation["MAX_PROBABILITY"] = probabilities.max(axis=1)
    validation["ASSIGNED"] = np.where(
        validation["MAX_PROBABILITY"] >= args.probability_threshold,
        validation["PREDICTED"],
        "ADMIXED_OR_UNCERTAIN",
    )
    validation_path = output / "classifier_validation_predictions.tsv"
    validation.to_csv(validation_path, sep="\t", index=False)
    copied["classifier_validation_predictions"] = validation_path

    matrix = confusion_matrix(y, predictions, labels=classes)
    per_class = {}
    for index, label in enumerate(classes):
        denominator = matrix[index].sum()
        per_class[label] = {
            "n": int(denominator),
            "recall": float(matrix[index, index] / denominator),
        }
    metrics = {
        "validation": "leave_one_constituent_population_out",
        "n_samples": int(len(joined)),
        "n_populations": int(len(np.unique(populations))),
        "classes": classes,
        "accuracy": float(accuracy_score(y, predictions)),
        "log_loss": float(log_loss(y, probabilities, labels=classes)),
        "coverage_at_threshold": float(
            (validation["MAX_PROBABILITY"] >= args.probability_threshold).mean()
        ),
        "accuracy_at_threshold": float(
            (validation.loc[
                validation["MAX_PROBABILITY"] >= args.probability_threshold,
                "PREDICTED",
            ] == validation.loc[
                validation["MAX_PROBABILITY"] >= args.probability_threshold,
                "Superpopulation",
            ]).mean()
        ),
        "confusion_matrix": matrix.tolist(),
        "per_class": per_class,
    }
    metrics_path = output / "classifier_validation_metrics.json"
    metrics_path.write_text(json.dumps(metrics, indent=2, sort_keys=True) + "\n")
    copied["classifier_validation_metrics"] = metrics_path

    final_model = new_classifier()
    final_model.fit(X, y)
    model_path = output / "classifier.joblib"
    joblib.dump(final_model, model_path)
    copied["classifier_joblib"] = model_path
    portable_model_path = output / "classifier.json"
    portable_model_path.write_text(json.dumps(
        classifier_json(final_model, pc_columns, args.probability_threshold),
        indent=2,
        sort_keys=True,
    ) + "\n")
    copied["classifier"] = portable_model_path

    labels_path = output / "ancestry_labels.tsv"
    labels.to_csv(labels_path, sep="\t", index=False)
    copied["ancestry_labels"] = labels_path

    metadata_path = output / "release_metadata.json"
    release_metadata = {
        "reference_id": args.reference_id,
        "genome_build": "GRCh38",
        "source_dataset": "1000 Genomes 30x high-coverage, 3202 samples",
        "source_directory": str(source),
        "pca_training_samples": int(sum(1 for _ in source_files["unrelated_ids"].open()) - 1),
        "reference_samples": int(len(joined)),
        "panel_variants": int(len(panel)),
        "projected_pcs": pc_columns,
        "plink_source_version": "2.00a3LM (20 Sep 2021)",
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
        "scikit_learn_version": sklearn.__version__,
        "probability_threshold": args.probability_threshold,
    }
    metadata_path.write_text(json.dumps(release_metadata, indent=2, sort_keys=True) + "\n")
    copied["release_metadata"] = metadata_path

    checksums_path = output / "SHA256SUMS"
    manifest = {
        "reference_id": args.reference_id,
        "genome_build": "GRCh38",
        "panel_variants": panel_path.name,
        "allele_frequencies": copied["allele_frequencies"].name,
        "loadings": copied["loadings"].name,
        "reference_scores": copied["reference_scores"].name,
        "ancestry_labels": labels_path.name,
        "classifier": model_path.name,
        "classifier_metadata": portable_model_path.name,
        "checksums": checksums_path.name,
    }
    manifest_json_path = output / "reference_manifest.json"
    manifest_json_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    copied["reference_manifest_json"] = manifest_json_path
    manifest_tsv_path = output / "reference_manifest.tsv"
    pd.DataFrame([manifest]).to_csv(manifest_tsv_path, sep="\t", index=False)
    copied["reference_manifest_tsv"] = manifest_tsv_path

    with checksums_path.open("w") as handle:
        for path in sorted(copied.values(), key=lambda item: item.name):
            handle.write(f"{sha256(path)}  {path.name}\n")


if __name__ == "__main__":
    main()
