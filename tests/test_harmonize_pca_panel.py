import csv
import importlib.util
import subprocess
import sys
import tracemalloc
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / 'bin/harmonize_pca_panel.py'
spec = importlib.util.spec_from_file_location('panel', SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def fixture(tmp_path):
    panel = tmp_path / 'panel.tsv'
    panel.write_text('CHROM\tPOS\tID\tREF\tALT\n' + ''.join(
        f'1\t{pos}\t{name}\tA\tC\n' for pos, name in [
            (10, 'exact'), (20, 'swap'), (30, 'dup'), (40, 'strand'),
            (50, 'incompatible'), (60, 'absent'), (70, 'exact')]))
    pvar = tmp_path / 'cohort.pvar'
    pvar.write_text('##meta=test\n#CHROM\tPOS\tID\tREF\tALT\n'
        '1\t999\tignore\tA\tC\n'
        '1\t10\ta\tA\tC\n1\t20\tb\tC\tA\n'
        '1\t30\tc\tA\tC\n1\t30\td\tC\tA\n'
        '1\t40\te\tT\tG\n1\t50\tf\tA\tG\n'
        '1\t70\tg\tA\tC\n')
    return panel, pvar


def run(tmp_path, threshold):
    panel, pvar = fixture(tmp_path)
    out = tmp_path / 'out'
    result = subprocess.run([sys.executable, str(SCRIPT), '--panel', str(panel),
        '--cohort-pvar', str(pvar), '--min-overlap', str(threshold),
        '--output-dir', str(out)], capture_output=True, text=True)
    return result, out


def test_matching_and_duplicate_semantics(tmp_path):
    result, out = run(tmp_path, 0.25)
    assert result.returncode == 0, result.stderr
    with (out / 'panel_overlap.tsv').open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    assert [r['STATUS'] for r in rows] == ['exact', 'ref_alt_swapped',
        'duplicate_compatible', 'strand_complement_rejected',
        'incompatible_alleles', 'absent', 'duplicate_panel_id']
    assert (out / 'usable_cohort_ids.txt').read_text() == 'a\nb\n'
    assert (out / 'rename_to_panel_ids.tsv').read_text() == 'exact\ta\nswap\tb\n'
    assert (out / 'reference_alleles.tsv').read_text() == 'exact\tA\nswap\tA\n'


def test_low_overlap_warns_and_preserves_usable_markers(tmp_path):
    result, out = run(tmp_path, 0.90)
    assert result.returncode == 0
    assert 'WARNING: PCA panel overlap' in result.stderr
    assert (out / 'usable_cohort_ids.txt').read_text() == 'a\nb\n'
    assert (out / 'panel_overlap_summary.tsv').exists()


def test_streamed_off_panel_rows_do_not_accumulate(tmp_path):
    pvar = tmp_path / 'large.pvar'
    with pvar.open('w') as f:
        f.write('#CHROM\tPOS\tID\tREF\tALT\n')
        for i in range(200000):
            f.write(f'1\t{i}\tv{i}\tA\tC\n')
    tracemalloc.start()
    rows = list(module.read_pvar(pvar, {('1', '7'), ('1', '199999')}))
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    assert [r['ID'] for r in rows] == ['v7', 'v199999']
    assert peak < 2 * 1024**2


def test_missing_header(tmp_path):
    pvar = tmp_path / 'bad.pvar'
    pvar.write_text('1\t1\ta\tA\tC\n')
    try:
        list(module.read_pvar(pvar, {('1', '1')}))
    except ValueError as e:
        assert 'No #CHROM header' in str(e)
    else:
        raise AssertionError('Missing header accepted')


def test_zero_overlap_still_fails(tmp_path):
    panel, pvar = fixture(tmp_path)
    pvar.write_text('#CHROM\tPOS\tID\tREF\tALT\n2\t999\tunmatched\tA\tG\n')
    result = subprocess.run([sys.executable, str(SCRIPT), '--panel', str(panel),
        '--cohort-pvar', str(pvar), '--min-overlap', '0.9',
        '--output-dir', str(tmp_path/'out')], capture_output=True, text=True)
    assert result.returncode != 0
    assert 'zero usable markers' in result.stderr
