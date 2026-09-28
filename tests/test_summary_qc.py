"""Native Nextflow/PLINK checks for QC provenance and invocation routing."""
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]
pytestmark = pytest.mark.skipif(
    not shutil.which('nextflow') or not shutil.which('plink2'),
    reason='Nextflow and native PLINK2 required',
)


@pytest.fixture(scope='module')
def fixture(tmp_path_factory):
    root = tmp_path_factory.mktemp('summary_qc')
    samples = ['s%d' % i for i in range(120)]
    with (root/'input.vcf').open('w') as handle:
        handle.write('##fileformat=VCFv4.2\n##contig=<ID=22,length=100000>\n'
                     '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        handle.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'
                     + '\t'.join(samples) + '\n')
        for i in range(20):
            genotypes = ['0/0', '0/1', '1/1']*40 if i < 18 else ['0/0']*120
            handle.write('22\t%d\t22:%d:A:G\tA\tG\t.\tPASS\t.\tGT\t%s\n'
                         % (i+1, i+1, '\t'.join(genotypes)))
    (root/'rsid.map').write_text(''.join('22:%d:A:G\trs%d\n' % (i, i) for i in range(1, 21)))
    subprocess.run(['plink2', '--vcf', str(root/'input.vcf'), '--make-pgen', 'vzs',
                    '--out', str(root/'base'), '--threads', '1', '--memory', '1500'], check=True)
    subprocess.run(['plink2', '--pfile', str(root/'base'), 'vzs', '--maf', '0.01',
                    '--update-name', str(root/'rsid.map'), '2', '1', '--make-pgen', 'vzs',
                    '--out', str(root/'supplied'), '--threads', '1', '--memory', '1500'], check=True)
    (root/'test.config').write_text('''
process.executor = 'local'
docker.enabled = false
process {
    withLabel: large { cpus = 1; memory = 3.GB }
    withLabel: small { cpus = 1; memory = 3.GB }
}
''')
    return root


@pytest.mark.parametrize('direct', [False, True])
@pytest.mark.parametrize('mode', ['base', 'prepare', 'prepare_no_map', 'supplied'])
def test_native_qc_provenance(fixture, tmp_path, direct, mode):
    target = 'base' if mode == 'base' else 'pgs'
    supplied = mode == 'supplied'
    preparing = mode.startswith('prepare')
    args = ['--input_pgs_pfile' if supplied else '--input_pfile',
            str(fixture/('supplied' if supplied else 'base')),
            '--prepare_pgs', str(preparing).lower()]
    if mode == 'prepare':
        args += ['--score_rsid_map', str(fixture/'rsid.map')]
    if supplied:
        # This is not the historical threshold used for the supplied fileset.
        args += ['--maf', '0.2']
    result = subprocess.run([
        'nextflow', 'run', os.environ.get('PGS_ENTRYPOINT', str(ROOT/'main.nf')), '-c', str(fixture/'test.config'),
        '-ansi-log', 'false', '--run_scores', 'false', '--run_pca', 'false',
        '--run_summary_qc', 'true', '--cohort', 'fixture', '--direct_inputs', str(direct).lower(),
        '--outdir', str(tmp_path/'results'), '--report_dir', str(tmp_path/'info'),
    ] + args, cwd=tmp_path, capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr
    folder = tmp_path/'results/04_summary'/target
    provenance = json.loads((folder/'qc_provenance.json').read_text())
    counts = dict(csv.reader((folder/'summary_counts.txt').open(), delimiter='\t'))
    assert provenance['qc_target'] == provenance['target'] == target
    assert provenance['cohort'] == 'fixture'
    assert provenance['genome_build'] == 'GRCh38'
    assert provenance['sample_count'] == int(counts['Samples']) == 120
    assert provenance['variant_count'] == int(counts['Variants']) == (20 if mode == 'base' else 18)
    if preparing:
        assert provenance['preparation'] == dict(
            source_input=str(fixture/'base'), maf_threshold=0.01,
            rsid_map=str(fixture/'rsid.map') if mode == 'prepare' else None)
    else:
        assert provenance['preparation'] is None
    source = provenance['input_prefix']
    assert provenance['pgen'] == source + '.pgen'
    assert provenance['pvar'] == source + '.pvar.zst'
    assert provenance['psam'] == source + '.psam'
    subprocess.run(['plink2', '--pfile', source, 'vzs', '--missing', '--hardy', '--freq',
                    '--out', str(tmp_path/'native'), '--threads', '1', '--memory', '1500'], check=True)
    for ext in ['smiss', 'vmiss', 'hardy', 'afreq']:
        assert (folder/('fixture.'+ext)).read_bytes() == (tmp_path/('native.'+ext)).read_bytes()
    ids = [line.split()[1] for line in (folder/'fixture.afreq').read_text().splitlines()[1:]]
    assert ids == (['rs%d' % i for i in range(1, 19)] if mode in ['prepare', 'supplied']
                   else ['22:%d:A:G' % i for i in range(1, 21 if mode == 'base' else 19)])
