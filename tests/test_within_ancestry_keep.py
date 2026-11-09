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
