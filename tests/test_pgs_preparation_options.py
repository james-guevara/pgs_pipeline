"""Real entrypoint validation: invalid combinations fail before tasks launch."""
import shutil
import subprocess
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1]
@pytest.mark.skipif(not shutil.which('nextflow'),reason='Nextflow unavailable')
@pytest.mark.parametrize('extra,message',[
 (['--input_pfile','/tmp/base'],'mutually exclusive'),
 (['--prepare_pgs','true'],'already prepared'),
 (['--score_rsid_map','/tmp/map'],'already prepared'),
 (['--run_pca','true'],'PCA requires the base'),
 (['--run_scores','true'],'--score_sheet is required'),
])
def test_rejects_ambiguous_prepared_input(tmp_path,extra,message):
 r=subprocess.run(['nextflow','run',str(ROOT/'main.nf'),'-ansi-log','false',
  '--input_pgs_pfile','/tmp/prepared','--direct_inputs','true',
  '--run_summary_qc','false']+extra,cwd=tmp_path,text=True,capture_output=True,timeout=120)
 assert r.returncode!=0
 assert message in r.stdout+r.stderr
 assert not list((tmp_path/'work').glob('*/*/.command.sh'))

@pytest.mark.skipif(not shutil.which('nextflow'),reason='Nextflow unavailable')
def test_disabled_preparation_ignores_unused_map(tmp_path):
 r=subprocess.run(['nextflow','run',str(ROOT/'main.nf'),'-ansi-log','false',
  '--input_pfile','/tmp/base','--direct_inputs','true',
  '--score_rsid_map','/tmp/nonexistent-unused-map','--run_summary_qc','false'],
  cwd=tmp_path,text=True,capture_output=True,timeout=120)
 assert r.returncode==0, r.stdout+r.stderr
 assert not list((tmp_path/'work').glob('*/*/.command.sh'))
