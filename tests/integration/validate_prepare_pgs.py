"""Validate real Nextflow/PLINK preparation and reuse, without printing genotypes."""
import csv, json, subprocess, sys
from pathlib import Path
root=Path(sys.argv[1])
def trace(case):
 return list(csv.DictReader((root/case/'info/trace.txt').open(),delimiter='\t'))
def names(case):return [r['name'] for r in trace(case)]
for case in ['prepare','prepare_direct','both','reuse','reuse_direct','no_map']:
 assert all(r['status']=='COMPLETED' for r in trace(case)),case
for case in ['prepare','prepare_direct','no_map']:
 assert any('PREPARE_SCORE_PFILE' in n for n in names(case))
 assert not any('SCORE_TRAIT' in n or 'PCA' in n for n in names(case))
for case in ['reuse','reuse_direct']:
 assert not any('PREPARE_SCORE_PFILE' in n for n in names(case))
 assert any('SCORE_TRAIT' in n for n in names(case))
 for f in ['combined_scores.tsv','score_qc_summary.tsv']:
  assert (root/case/'results/05_scores'/f).read_bytes()==(root/'both/results/05_scores'/f).read_bytes(),(case,f)
assert not names('neither')
for case in ['prepare','prepare_direct','both','no_map']:
 folder=root/case/'results/04_score_input'
 for suffix in ['pgen','pvar.zst','psam']:assert (folder/('score_input.'+suffix)).is_file()
 assert not (folder/'score_input.pvar').exists()
 summary=dict(csv.reader((folder/'score_input_filter_summary.tsv').open(),delimiter='\t'))
 assert summary['variants_before']=='20' and summary['variants_after']=='18'
 assert summary['rsid_mapping']==('not_requested' if case=='no_map' else 'applied')
 assert json.loads((folder/'score_input_provenance.json').read_text())['genome_build']=='GRCh38'
 # Export base and prepared genotypes for comparison, ignoring intentional ID changes.
 subprocess.run(['plink2','--pfile',str(folder/'score_input'),'vzs','--threads','1','--memory','1500','--export','A','--out',str(root/case/'genotypes')],check=True,stdout=subprocess.DEVNULL)
 table=list(csv.reader((root/case/'genotypes.raw').open(),delimiter='\t'))
 assert len(table)==121
 if case=='prepare':reference=table
 else:assert table[1:]==reference[1:],case
# Compare retained genotype calls and sample order against the original base.
subprocess.run(['plink2','--pfile',str(root/'fixture/base'),'vzs','--maf','0.01',
                '--threads','1','--memory','1500','--export','A','--out',str(root/'base_filtered')],
               check=True,stdout=subprocess.DEVNULL)
base=list(csv.reader((root/'base_filtered.raw').open(),delimiter='\t'))
assert base[1:]==reference[1:], 'Preparation changed retained genotypes or sample order'
pub=root/'prepare_direct/results/04_score_input/score_input.pgen'
assert any(p.stat().st_ino==pub.stat().st_ino for p in (root/'prepare_direct/work').glob('*/*/score_input.pgen'))
(root/'validation.json').write_text(json.dumps({'status':'PASS','cases':7,'samples':120,'retained_variants':18},indent=2))
print('PASS: preparation, mapping, filtering, compressed PVAR, scoring reuse, direct inputs, hardlink publication')
