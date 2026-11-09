import csv
import importlib.util
import subprocess
import sys
import tracemalloc
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / 'bin/explain_score_matching.py'
spec = importlib.util.spec_from_file_location('explain', SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def fixture(root, separator='\t'):
    (root/'weights').write_text('SNP\tA1\tBETA\n' + ''.join(
        x+'\t1\n' for x in ['rsMatch\tA','rsMAF\tA','rsBad\tT','rsOther\tA',
        'rsAbsent\tA','rsUnknown\tA','rsUnmapped\tA','rsMAF\tG',
        'rsOverwritten\tA','rsRenamed\tA']))
    mappings=[('a','rsMatch'),('b','rsMAF'),('c','rsBad'),('d','rsOther'),
        ('missing','rsAbsent'),('z','rsOverwritten'),('z','irrelevant'),('rsRenamed','irrelevant')]
    (root/'map').write_text(''.join(separator.join(pair)+'\n' for pair in mappings))
    header='#CHROM\tPOS\tID\tREF\tALT\n'
    (root/'source').write_text(header+''.join(f'1\t1\t{x}\tA\tG\n' for x in
        ['a','b','c','d','rsUnmapped','z','rsRenamed','ignored']))
    (root/'score').write_text(header+''.join(f'1\t1\t{x}\tA\tG\n' for x in
        ['rsMatch','rsBad','rsOther','ignored']))
    (root/'matched').write_text('rsMatch\nignored\n')
    return ['--trait','test','--weights',str(root/'weights'),'--id-col','1','--allele-col','2',
        '--source-pvar',str(root/'source'),'--score-pvar',str(root/'score'),
        '--rsid-map',str(root/'map'),'--matched-vars',str(root/'matched')]


def test_categories_duplicates_and_whitespace(tmp_path):
    outputs=[]
    for name,separator in [('tabs','\t'),('spaces',' ')]:
        root=tmp_path/name;root.mkdir()
        args=fixture(root,separator)
        subprocess.run([sys.executable,str(SCRIPT)]+args+['--summary',str(root/'summary'),
            '--details',str(root/'details')],check=True)
        outputs.append(((root/'summary').read_bytes(),(root/'details').read_bytes()))
        with (root/'details').open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
        assert [r['category'] for r in rows]==['removed_by_score_view_maf',
            'effect_allele_incompatible','present_but_not_scored_other',
            'absent_from_harmonized_source','no_selected_rsid_mapping_or_absent',
            'source_present_without_selected_mapping','removed_by_score_view_maf',
            'absent_from_harmonized_source','no_selected_rsid_mapping_or_absent']
        assert [r['weight_rows_for_id'] for r in rows if r['variant_id']=='rsMAF']==['2','2']
    assert outputs[0]==outputs[1]


def test_unrequested_variants_and_mappings_do_not_accumulate(tmp_path):
    mapping=tmp_path/'map';pvar=tmp_path/'pvar'
    with mapping.open('w') as m,pvar.open('w') as p:
        p.write('#CHROM\tPOS\tID\tREF\tALT\n')
        for i in range(200000):
            m.write(f'v{i}\trs{i}\n');p.write(f'1\t{i}\trs{i}\tA\tG\n')
    tracemalloc.start()
    result=module.relevant_mapping(mapping,{'rs7'})
    variants=module.pvar(pvar,{'rs7'})
    _,peak=tracemalloc.get_traced_memory();tracemalloc.stop()
    assert result==({'v7':'rs7'},{'rs7'},set(),True)
    assert variants=={'rs7':('A','G')}
    assert peak<2*1024**2


def test_no_map_and_map_outside_requested_are_distinct(tmp_path):
    p=tmp_path/'map';p.write_text('# empty\n')
    assert module.relevant_mapping(p,{'rs7'})==({},set(),set(),False)
    p.write_text('unrelated\tother\n')
    assert module.relevant_mapping(p,{'rs7'})==({},set(),set(),True)
