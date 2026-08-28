#!/usr/bin/env bash
set -euo pipefail

scores_out=$1
qc_out=$2
shift 2

scores=()
qcs=()
mode=scores
for value in "$@"; do
    if [[ "$value" == "--" ]]; then
        mode=qcs
    elif [[ "$mode" == "scores" ]]; then
        scores+=("$value")
    else
        qcs+=("$value")
    fi
done

if (( ${#scores[@]} == 0 || ${#qcs[@]} == 0 )); then
    echo "At least one score and QC file are required" >&2
    exit 2
fi

reference_ids=$(mktemp)
next_ids=$(mktemp)
next_column=$(mktemp)
trap 'rm -f "$reference_ids" "$next_ids" "$next_column"' EXIT

first=true
for score in "${scores[@]}"; do
    trait=${score%.sscore}
    awk -v trait="$trait" '
        BEGIN { FS=OFS="\t" }
        NR == 1 {
            for (i=1; i<=NF; i++) {
                if ($i == "#IID") iid=i
                if ($i == "SCORE1_AVG") score=i
            }
            if (!iid || !score) exit 2
            print "IID", trait
            next
        }
        { print $iid, $score }
    ' "$score" > "$next_column"

    cut -f1 "$next_column" > "$next_ids"
    if [[ "$first" == true ]]; then
        cp "$next_column" "$scores_out"
        cp "$next_ids" "$reference_ids"
        first=false
    else
        if ! cmp -s "$reference_ids" "$next_ids"; then
            echo "Sample IDs or order in $score differ from the other score files" >&2
            exit 2
        fi
        paste "$scores_out" <(cut -f2 "$next_column") > "${scores_out}.tmp"
        mv "${scores_out}.tmp" "$scores_out"
    fi
done

head -n 1 "${qcs[0]}" > "$qc_out"
for qc in "${qcs[@]}"; do
    tail -n +2 "$qc" >> "$qc_out"
done
