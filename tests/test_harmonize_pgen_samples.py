import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pgenlib as pg
import pyarrow.parquet as pq


ROOT = Path(__file__).resolve().parents[1]


class HarmonizationTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def fixture(self, name, variants, dosages, samples):
        prefix = self.root / name
        with pg.PgenWriter(str(prefix.with_suffix('.pgen')).encode(), len(samples),
                           variant_ct=len(variants), dosage_present=True,
                           nonref_flags=False) as writer:
            for values in dosages:
                writer.append_dosages(np.asarray(values, dtype=np.float64))
        prefix.with_suffix('.pvar').write_text(
            '#CHROM\tPOS\tID\tREF\tALT\n' + ''.join(
                f'22\t{pos}\t{name}_{pos}\t{ref}\t{alt}\n'
                for pos, ref, alt in variants))
        prefix.with_suffix('.psam').write_text(
            '#IID\tSEX\n' + ''.join(f'{sample}\tNA\n' for sample in samples))
        return prefix

    def test_compact_qc_contract_and_dosage_merge(self):
        first = self.fixture('wgs', [(10, 'A', 'C'), (20, 'G', 'T')],
                             [[0, 1], [2, -9]], ['w1', 'w2'])
        second = self.fixture('imputed', [(10, 'A', 'C'), (30, 'C', 'T')],
                              [[0.3, 1.7], [0, 1]], ['i1', 'i2'])
        manifest = self.root / 'inputs.tsv'
        manifest.write_text(
            f'source\tpgen_prefix\nWGS\t{first}\nImputed\t{second}\n')
        output = self.root / 'output'
        result = subprocess.run(
            [sys.executable, str(ROOT / 'bin/harmonize_pgen_samples.py'),
             '--manifest', str(manifest), '--out-dir', str(output)],
            text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        summary = json.loads((output / 'harmonization_summary.json').read_text())
        self.assertEqual(summary['status'], 'validated')
        self.assertEqual(summary['output_samples'], 4)
        self.assertEqual(summary['union_markers'], 3)
        self.assertEqual(summary['output_markers'], 1)
        self.assertEqual((output / 'sample_qc.tsv').read_text().count('\n'), 5)
        marker = pq.read_table(output / 'marker_qc.parquet').to_pylist()
        self.assertEqual(len(marker), 3)
        retained = [row for row in marker if row['retained']]
        self.assertEqual(len(retained), 1)
        self.assertTrue(retained[0]['wgs_present'])
        self.assertTrue(retained[0]['imputed_present'])
        observed = np.empty(4, dtype=np.float64)
        with pg.PgenReader(str(output / 'harmonized.pgen').encode()) as reader:
            reader.read_dosages(0, observed)
        np.testing.assert_allclose(observed, [0, 1, 0.3, 1.7], atol=1/16384)


if __name__ == '__main__':
    unittest.main()
