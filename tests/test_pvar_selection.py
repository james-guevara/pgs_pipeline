"""Exercise input format selection without scheduling genotype computations."""
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.skipif(not shutil.which('nextflow'), reason='Nextflow unavailable')
@pytest.mark.parametrize('format,success', [('auto', False), ('pvar', True), ('pvar.zst', True), ('invalid', False)])
def test_ambiguous_input_requires_explicit_format(tmp_path, format, success):
    prefix = tmp_path / 'input'
    for suffix in ('.pvar', '.pvar.zst', '.pgen', '.psam'):
        Path(str(prefix) + suffix).touch()
    result = subprocess.run(
        ['nextflow', 'run', str(ROOT / 'main.nf'), '-ansi-log', 'false',
         '--input_pfile', str(prefix), '--input_pvar_format', format,
         '--direct_inputs', 'true', '--run_summary_qc', 'false'],
        cwd=tmp_path, text=True, capture_output=True, timeout=120)
    output = result.stdout + result.stderr
    assert (result.returncode == 0) == success, output
    if format == 'auto':
        assert 'Both .pvar and .pvar.zst exist' in output
    elif format == 'invalid':
        assert '--input_pvar_format must be' in output
