#!/usr/bin/env python3
"""Apply the runtime-neutral Extra Trees format emitted by the reference builder."""

import argparse
import csv
import json
from pathlib import Path


def tree_probability(tree, features):
    node = 0
    while tree["children_left"][node] != -1:
        feature = tree["feature"][node]
        if features[feature] <= tree["threshold"][node]:
            node = tree["children_left"][node]
        else:
            node = tree["children_right"][node]
    values = tree["value"][node]
    total = sum(values)
    return [value / total for value in values]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--classifier-metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    model = json.loads(args.model.read_text())
    metadata = json.loads(args.classifier_metadata.read_text())
    if model.get("format") != "pgs_pipeline_extra_trees_v1":
        raise ValueError("Unsupported classifier format")
    if model["features"] != metadata["pc_columns"]:
        raise ValueError("Classifier feature metadata does not match model")
    classes = model["classes"]
    threshold = float(metadata["assignment_probability_threshold"])

    with args.scores.open(newline="") as source, args.output.open("w", newline="") as target:
        reader = csv.DictReader(source, delimiter="\t")
        id_column = "#IID" if "#IID" in reader.fieldnames else "IID"
        output_fields = [id_column, "ANCESTRY", "MAX_PROBABILITY"] + [
            f"PROB_{label}" for label in classes
        ]
        writer = csv.DictWriter(target, fieldnames=output_fields, delimiter="\t")
        writer.writeheader()
        for row in reader:
            features = [float(row[column]) for column in model["features"]]
            probabilities = [0.0] * len(classes)
            for tree in model["trees"]:
                for index, value in enumerate(tree_probability(tree, features)):
                    probabilities[index] += value
            probabilities = [value / len(model["trees"]) for value in probabilities]
            best = max(range(len(classes)), key=probabilities.__getitem__)
            assignment = classes[best] if probabilities[best] >= threshold else metadata["uncertain_label"]
            result = {
                id_column: row[id_column],
                "ANCESTRY": assignment,
                "MAX_PROBABILITY": f"{probabilities[best]:.10g}",
            }
            result.update({f"PROB_{label}": f"{probabilities[index]:.10g}" for index, label in enumerate(classes)})
            writer.writerow(result)


if __name__ == "__main__":
    main()
