#!/usr/bin/env python3
"""Compare full PGS/PCA runs of the same input in plain and compressed PVAR form.

Usage: python compare_pvar_runs.py EXPERIMENT_ROOT
Expected subdirectories: plain-staged, compressed-staged, compressed-direct.
Each contains the Nextflow results/ directory. Run all three without -resume.
No participant-level data is printed or included in validation.json.
"""
import hashlib
import json
import sys
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def compare(root):
    baseline = root / 'plain-staged/results'
    # Exclude runtime logs, which contain timestamps and work paths. Compare
    # all published scientific tables, filesets, and within-ancestry outputs.
    files = sorted(p.relative_to(baseline) for folder in
                   ('04_summary', '04_score_input', '05_scores', '06_pca', '07_analysis')
                   for p in (baseline / folder).rglob('*')
                   if p.is_file() and p.suffix != '.log')
    assert files, 'No baseline results'
    for required in ('05_scores/combined_scores.tsv', '06_pca/global/global_pcs.tsv',
                     '06_pca/ancestry/ancestry_probabilities.tsv',
                     '06_pca/within_ancestry/status.tsv', '07_analysis/analysis_dataset.tsv'):
        assert Path(required) in files, required
    assert all((baseline / '06_pca/within_ancestry' / group / 'pcs.tsv').is_file()
               for group in ('AFR', 'AMR', 'EAS', 'EUR', 'SAS'))
    checks = []
    for case in ('compressed-staged', 'compressed-direct'):
        results = root / case / 'results'
        for rel in files:
            assert (results / rel).is_file(), (case, str(rel), 'missing')
            assert digest(baseline / rel) == digest(results / rel), (case, str(rel), 'differs')
        checks.append({'case': case, 'identical_scientific_artifacts': len(files), 'status': 'PASS'})
    report = {'status': 'PASS', 'checks': checks}
    (root / 'validation.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    compare(Path(sys.argv[1]))
