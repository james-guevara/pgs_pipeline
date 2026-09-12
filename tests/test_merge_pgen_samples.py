import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import subprocess
import sys

import numpy as np
import pgenlib as pg

SPEC = importlib.util.spec_from_file_location(
    "merge_pgen_samples", Path(__file__).resolve().parents[1] / "bin/merge_pgen_samples.py")
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


class MergeTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def fixture(self, name, variants, calls, ids=None, dosage=False, phased=False):
        prefix = self.root / name
        n = len(calls[0]) // (2 if phased else 1)
        with pg.PgenWriter(str(prefix.with_suffix('.pgen')).encode(), n,
                           variant_ct=len(variants), nonref_flags=False,
                           dosage_present=dosage, hardcall_phase_present=phased) as writer:
            for row in calls:
                if dosage:
                    writer.append_dosages(np.array(row, dtype=np.float64))
                elif phased:
                    writer.append_alleles(np.array(row, dtype=np.int32), all_phased=True)
                else:
                    writer.append_biallelic(np.array(row, dtype=np.int8))
        prefix.with_suffix('.pvar').write_text('#CHROM\tPOS\tID\tREF\tALT\n' + ''.join(
            f'1\t{pos}\t{name}_{i}\t{ref}\t{alt}\n'
            for i, (pos, ref, alt) in enumerate(variants)))
        ids = ids or [f'{name}{i}' for i in range(n)]
        prefix.with_suffix('.psam').write_text('#IID\tSEX\n' + ''.join(f'{iid}\t0\n' for iid in ids))
        return prefix

    def test_reordered_intersection_missingness_metadata_and_counts(self):
        a = self.fixture('a', [(10,'A','C'), (20,'G','T'), (30,'C','T')],
                         [[0,1,-9], [2,-9,0], [1,1,1]])
        b = self.fixture('b', [(20,'G','T'), (40,'A','G'), (10,'A','C'), (30,'T','C')],
                         [[1,2], [0,0], [-9,2], [0,0]])
        out = self.root / 'output'
        summary = MOD.merge(a, b, out)
        self.assertEqual((summary['output_samples'], summary['output_variants']), (5,2))
        self.assertEqual((summary['first_unmatched'], summary['second_unmatched']), (1,2))
        observed = np.empty((2,5), dtype=np.int8)
        with pg.PgenReader(str(out / 'merged.pgen').encode()) as reader:
            reader.read_range(0,2,observed)
        np.testing.assert_array_equal(observed, [[0,1,-9,-9,2],[2,-9,0,1,2]])
        self.assertEqual((out / 'merged.psam').read_text().splitlines(),
                         ['#FID\tIID\tSEX', '0\ta0\t0', '0\ta1\t0',
                          '0\ta2\t0', '0\tb0\t0', '0\tb1\t0'])
        self.assertIn('1:10:A:C', (out / 'merged.pvar').read_text())
        self.assertEqual(json.loads((out / 'merge_summary.json').read_text()), summary)
        with self.assertRaisesRegex(ValueError, 'already exists'):
            MOD.merge(a,b,out)

    def test_dosage_refused_without_partial_output(self):
        a = self.fixture('a', [(10,'A','C')], [[0,1]])
        b = self.fixture('b', [(10,'A','C')], [[0.3,1.7]], dosage=True)
        with self.assertRaisesRegex(ValueError, 'Dosage differs'):
            MOD.merge(a,b,self.root / 'output')
        self.assertFalse((self.root / 'output').exists())

    def test_dosage_preserved_across_sample_axis(self):
        a = self.fixture('a', [(10, 'A', 'C')], [[0, 1]])
        b = self.fixture('b', [(10, 'A', 'C')], [[0.3, 1.7]], dosage=True)
        out = self.root / 'output'
        summary = MOD.merge(a, b, out, dosage_mode='preserve')
        observed = np.empty(4, dtype=np.float64)
        with pg.PgenReader(str(out / 'merged.pgen').encode()) as reader:
            reader.read_dosages(0, observed)
        np.testing.assert_allclose(observed, [0, 1, 0.3, 1.7], atol=1 / 16384)
        self.assertEqual(summary['dosage_policy'], 'preserve')
        self.assertFalse(summary['phase_preserved'])

    def test_phase_preserved(self):
        a = self.fixture('a', [(10,'A','C')], [[1,0,0,1]], phased=True)
        b = self.fixture('b', [(10,'A','C')], [[1,-9]])
        MOD.merge(a,b,self.root / 'output')
        alleles, phase = np.empty(8,dtype=np.int32), np.empty(4,dtype=np.uint8)
        with pg.PgenReader(str(self.root / 'output/merged.pgen').encode()) as reader:
            reader.read_alleles_and_phasepresent(0,alleles,phase)
        np.testing.assert_array_equal(alleles, [1,0,0,1,0,1,-9,-9])
        np.testing.assert_array_equal(phase, [1,1,0,0])

    def test_three_inputs_cli_intersection_and_phase(self):
        a = self.fixture('a', [(10,'A','C'), (20,'G','T'), (30,'C','T')],
                         [[0,1], [2,-9], [1,1]])
        b = self.fixture('b', [(20,'G','T'), (10,'A','C')], [[1],[2]])
        c = self.fixture('c', [(30,'C','T'), (20,'G','T')],
                         [[0,0], [1,0]], phased=True)
        out = self.root / 'output'
        result = subprocess.run([sys.executable, str(SPEC.origin), '--inputs',
                                 str(a), str(b), str(c), '--out-dir', str(out)],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        summary = json.loads(result.stdout)
        self.assertEqual(summary['input_count'], 3)
        self.assertEqual(summary['output_variants'], 1)
        self.assertEqual([row['unmatched'] for row in summary['inputs']], [2,1,1])
        self.assertEqual([row['prefix'] for row in summary['inputs']],
                         [str(p.resolve()) for p in (a,b,c)])
        alleles, phase = np.empty(8,dtype=np.int32), np.empty(4,dtype=np.uint8)
        with pg.PgenReader(str(out / 'merged.pgen').encode()) as reader:
            reader.read_alleles_and_phasepresent(0,alleles,phase)
        np.testing.assert_array_equal(alleles, [1,1,-9,-9,0,1,1,0])
        np.testing.assert_array_equal(phase, [1,0,0,1])
        self.assertEqual((out / 'merged.psam').read_text().splitlines()[1:],
                         ['0\ta0\t0','0\ta1\t0','0\tb0\t0','0\tc0\t0'])
        self.assertIn('1:20:G:T', (out / 'merged.pvar').read_text())

    def test_later_input_validation(self):
        a = self.fixture('a', [(10,'A','C')], [[0]])
        b = self.fixture('b', [(10,'A','C')], [[1]])
        c = self.fixture('c', [(10,'A','C')], [[2]], ids=['b0'])
        out = self.root / 'out'
        with self.assertRaisesRegex(ValueError, 'Overlapping'):
            MOD.merge_many([a,b,c], out)
        c = self.fixture('c', [(20,'A','C')], [[2]])
        with self.assertRaisesRegex(ValueError, 'No matching'):
            MOD.merge_many([a,b,c], out)
        c = self.fixture('c', [(10,'A','C')], [[0.3]], dosage=True)
        with self.assertRaisesRegex(ValueError, 'Dosage differs'):
            MOD.merge_many([a,b,c], out)
        self.assertFalse(out.exists())
        # Missing metadata in input one must not mask a conflict in later inputs.
        for prefix, reference in ((b, 'build1'), (c, 'build2')):
            path = prefix.with_suffix('.pvar')
            path.write_text(f'##reference={reference}\n' + path.read_text())
        with self.assertRaisesRegex(ValueError, 'reference metadata differs'):
            MOD.merge_many([a,b,c], out)

    def test_cli_input_validation_and_legacy_command(self):
        a = self.fixture('a', [(10,'A','C')], [[0]])
        b = self.fixture('b', [(10,'A','C')], [[1]])
        base = [sys.executable, str(SPEC.origin), '--out-dir', str(self.root/'out')]
        for args in ([], ['--inputs', str(a)], ['--first', str(a)],
                     ['--inputs', str(a), str(b), '--first', str(a)]):
            result = subprocess.run(base + args, capture_output=True, text=True)
            self.assertEqual(result.returncode, 2, result.stderr)
        result = subprocess.run(base + ['--first', str(a), '--second', str(b)],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        with self.assertRaisesRegex(ValueError, 'At least two'):
            MOD.merge_many([a], self.root/'other')

    def test_overlap_duplicates_no_intersection_and_count_mismatch(self):
        a = self.fixture('a', [(10,'A','C')], [[0]])
        b = self.fixture('b', [(10,'A','C')], [[1]], ids=['a0'])
        with self.assertRaisesRegex(ValueError, 'Overlapping'):
            MOD.merge(a,b,self.root/'out')
        b = self.fixture('b', [(10,'A','C'),(10,'A','C')], [[0],[1]])
        with self.assertRaisesRegex(ValueError, 'duplicate variant'):
            MOD.merge(a,b,self.root/'out')
        b = self.fixture('b', [(11,'A','C')], [[0]])
        with self.assertRaisesRegex(ValueError, 'No matching'):
            MOD.merge(a,b,self.root/'out')
        b = self.fixture('b', [(10,'A','C')], [[0]])
        b.with_suffix('.psam').write_text('#IID\tSEX\nb0\t0\nb1\t0\n')
        with self.assertRaises(RuntimeError):
            MOD.merge(a,b,self.root/'out')
        self.assertFalse((self.root/'out').exists())

    def test_psam_schema_is_normalized_and_overlap_uses_iid(self):
        a = self.fixture('a', [(10, 'A', 'C')], [[0]], ids=['a0'])
        b = self.fixture('b', [(10, 'A', 'C')], [[1]], ids=['b0'])
        a.with_suffix('.psam').write_text('#IID\tSEX\na0\t1\n')
        b.with_suffix('.psam').write_text('#SID\tIID\tFID\nb_sid\tb0\tfam_b\n')
        out = self.root / 'out'
        MOD.merge(a, b, out)
        self.assertEqual((out / 'merged.psam').read_text().splitlines(), [
            '#FID\tIID\tSEX\tSID',
            '0\ta0\t1\tNA',
            'fam_b\tb0\tNA\tb_sid',
        ])

        c = self.fixture('c', [(10, 'A', 'C')], [[2]], ids=['a0'])
        c.with_suffix('.psam').write_text('#FID\tIID\tSEX\nother_family\ta0\t2\n')
        with self.assertRaisesRegex(ValueError, 'Overlapping'):
            MOD.merge(a, c, self.root / 'overlap')


if __name__ == '__main__':
    unittest.main()
