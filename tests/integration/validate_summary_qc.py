"""Check target selection, reports, and provenance against native PLINK QC."""
import csv
import json
from pathlib import Path
import subprocess
import sys

root = Path(sys.argv[1])
cases = ['qc_prepare', 'qc_prepare_direct', 'qc_both', 'qc_both_direct',
         'qc_reuse', 'qc_reuse_direct', 'qc_base', 'qc_base_direct']
for target, prefix in [('base', root/'fixture/base'),
                       ('pgs', root/'prepare/results/04_score_input/score_input')]:
    subprocess.run(['plink2', '--pfile', str(prefix), 'vzs', '--missing', '--hardy', '--freq',
                    '--threads', '1', '--memory', '1500', '--out', str(root/('reference_'+target))],
                   check=True, stdout=subprocess.DEVNULL)
for case in cases:
    target = 'base' if case.startswith('qc_base') else 'pgs'
    folder = root/case/'results/04_summary'/target
    assert not (folder.parent/('pgs' if target == 'base' else 'base')).exists(), case
    counts = dict(csv.reader((folder/'summary_counts.txt').open(), delimiter='\t'))
    assert counts['Samples'] == '120' and counts['Variants'] == ('20' if target == 'base' else '18'), case
    provenance = json.loads((folder/'qc_provenance.json').read_text())
    assert provenance['target'] == target and provenance['genome_build'] == 'GRCh38', case
    assert provenance['input_prefix'].endswith('/base' if target == 'base' else '/score_input'), case
    assert provenance['pvar'].endswith('.pvar.zst'), case
    for ext in ['smiss', 'vmiss', 'hardy', 'afreq']:
        assert (folder/('fixture.'+ext)).read_bytes() == (root/('reference_'+target+'.'+ext)).read_bytes(), (case, ext)
    ids = [row.split()[1] for row in (folder/'fixture.afreq').read_text().splitlines()[1:]]
    assert ids == (['22:%d:A:G' % i for i in range(1,21)] if target == 'base' else ['rs%d' % i for i in range(1,19)]), case
    trace = list(csv.DictReader((root/case/'info/trace.txt').open(), delimiter='\t'))
    assert all(row['status'] == 'COMPLETED' for row in trace), case
    tasks = [row['name'] for row in trace]
    qc = [name for name in tasks if 'SUMMARY_QC' in name]
    direct = case in ['qc_reuse_direct', 'qc_base_direct']
    assert len(qc) == 1 and ('SUMMARY_QC_DIRECT' in qc[0]) == direct, case
    assert not any('PCA' in name for name in tasks), case
    if case.startswith(('qc_reuse','qc_base')):
        assert not any('PREPARE_SCORE_PFILE' in name or 'SCORE_TRAIT' in name for name in tasks), case
    if case.startswith('qc_both'):
        assert any('SCORE_TRAIT' in name for name in tasks), case
        assert (root/case/'results/05_scores/combined_scores.tsv').read_bytes() == (root/'both/results/05_scores/combined_scores.tsv').read_bytes(), case
# Disabled QC in the existing preparation suite remains disabled.
for case in ['prepare','prepare_direct','both','reuse','reuse_direct','neither','no_map']:
    assert not (root/case/'results/04_summary').exists(), case
(root/'summary-qc-validation.json').write_text(json.dumps(
    dict(status='PASS', qc_cases=8, existing_preparation_cases=7,
         samples=120, base_variants=20, pgs_variants=18,
         native_report_parity=True, stale_output_ignored=True), indent=2)+'\n')
print('PASS summary QC targets, counts, IDs, native reports, provenance, direct inputs, scoring, stale-output isolation')
