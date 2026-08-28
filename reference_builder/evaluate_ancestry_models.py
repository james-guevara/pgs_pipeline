#!/usr/bin/env python3
"""Compare ancestry classifiers with leave-one-population-out validation."""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis, QuadraticDiscriminantAnalysis
from sklearn.ensemble import ExtraTreesClassifier, RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.neighbors import KNeighborsClassifier, NearestCentroid
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC


def models():
    return {
        "nearest_centroid": Pipeline([
            ("scale", StandardScaler()),
            ("model", NearestCentroid()),
        ]),
        "lda_shrinkage": Pipeline([
            ("scale", StandardScaler()),
            ("model", LinearDiscriminantAnalysis(solver="lsqr", shrinkage="auto")),
        ]),
        "qda": Pipeline([
            ("scale", StandardScaler()),
            ("model", QuadraticDiscriminantAnalysis(reg_param=0.1)),
        ]),
        "knn_25_distance": Pipeline([
            ("scale", StandardScaler()),
            ("model", KNeighborsClassifier(n_neighbors=25, weights="distance")),
        ]),
        "rbf_svc": Pipeline([
            ("scale", StandardScaler()),
            ("model", SVC(C=10.0, gamma="scale")),
        ]),
        "random_forest": RandomForestClassifier(
            n_estimators=500, min_samples_leaf=2, class_weight="balanced",
            random_state=20260827, n_jobs=-1,
        ),
        "extra_trees": ExtraTreesClassifier(
            n_estimators=500, min_samples_leaf=2, class_weight="balanced",
            random_state=20260827, n_jobs=-1,
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    scores = pd.read_csv(args.scores, sep=r"\s+")
    metadata = pd.read_csv(args.metadata, sep=r"\s+")
    data = scores.merge(
        metadata[["SampleID", "Population", "Superpopulation"]],
        left_on="#IID", right_on="SampleID", validate="one_to_one",
    )
    pc_columns = [column for column in scores if column.startswith("PC") and column.endswith("_AVG")]
    X = data[pc_columns].to_numpy()
    y = data["Superpopulation"].to_numpy()
    population = data["Population"].to_numpy()

    results = {}
    predictions = data[["#IID", "Population", "Superpopulation"]].copy()
    for name, estimator in models().items():
        predicted = np.empty(len(data), dtype=object)
        for held_out in sorted(np.unique(population)):
            test = population == held_out
            estimator.fit(X[~test], y[~test])
            predicted[test] = estimator.predict(X[test])
        per_population = {}
        for held_out in sorted(np.unique(population)):
            selected = population == held_out
            per_population[held_out] = float(accuracy_score(y[selected], predicted[selected]))
        results[name] = {
            "accuracy": float(accuracy_score(y, predicted)),
            "worst_population_accuracy": float(min(per_population.values())),
            "per_population_accuracy": per_population,
        }
        predictions[name] = predicted

    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "model_comparison.json").write_text(
        json.dumps(results, indent=2, sort_keys=True) + "\n"
    )
    predictions.to_csv(args.output / "model_predictions.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
