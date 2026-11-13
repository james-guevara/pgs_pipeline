import subprocess
from pathlib import Path
import pytest

SCRIPT=Path(__file__).resolve().parents[1]/'bin/run_within_ancestry_pca.sh'

@pytest.mark.parametrize('header,row,expected',[
 ('#IID\tSEX','child1\t0','0\tchild1'),
 ('#FID\tIID\tSEX','family1\tchild1\t0','family1\tchild1'),
 ('#SEX\tIID\tFID','0\tchild1\tfamily1','family1\tchild1')])
def test_keep_list_uses_optional_fid(tmp_path,header,row,expected):
    text=SCRIPT.read_text()
    awk=text.split("awk -F '\\t' -v ancestry=\"$ancestry\" '",1)[1].split("' \"$ancestry_tsv\"",1)[0]
    assignments=tmp_path/'assign.tsv';assignments.write_text('#IID\tANCESTRY\tMOST_LIKELY_ANCESTRY\nchild1\tuncertain\tEUR\nchild2\tAFR\tAFR\n')
    psam=tmp_path/'test.psam';psam.write_text(header+'\n'+row+'\n')
    result=subprocess.check_output(['awk','-F','\t','-v','ancestry=EUR',awk,str(assignments),str(psam)],text=True)
    assert result=='#FID\tIID\n'+expected+'\n'


def test_rejects_unknown_genotype_mode(tmp_path):
    args = ['in', 'a.tsv', '50', '10', '0.0884', '200', '50', '0.2', '0.01', '1', '1000', str(tmp_path / 'out'), '', 'probabilistic']
    result = subprocess.run(['bash', str(SCRIPT), *args], capture_output=True, text=True)
    assert result.returncode == 2 and 'Invalid genotype mode' in result.stderr
    assert not (tmp_path / 'out').exists()


def test_hardcall_mode_erases_dosages_once_and_cleans_up():
    text = SCRIPT.read_text()
    block = text.split('if [[ $genotype_mode == hardcall ]]; then', 1)[1].split('\nfi\n', 1)[0]
    assert '--make-pgen vzs erase-dosage' in block
    assert 'trap' in block and 'rm -f' in block
    assert 'pfile_args=("$pfile" vzs)' in block
    # conversion happens before the per-ancestry loop, so it runs once
    assert text.index('erase-dosage') < text.index('for ancestry in AFR AMR EAS EUR SAS; do')
    # the hard-call copy is written outside the published output directory
    assert 'hardcall_prefix=$PWD/' in block
