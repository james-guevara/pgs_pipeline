import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "harmonize_cohorts", ROOT / "bin" / "harmonize_cohorts.py"
)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class HarmonizeCohortsTest(unittest.TestCase):
    def test_eligible_markers_use_maf_and_drop_duplicate_keys(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pvar = root / "test.pvar"
            afreq = root / "test.afreq"
            pvar.write_text(
                "#CHROM\tPOS\tID\tREF\tALT\n"
                "1\t10\t1:10:A:G\tA\tG\n"
                "1\t20\t1:20:C:T\tC\tT\n"
                "1\t30\tdup-a\tG\tA\n"
                "1\t30\tdup-b\tG\tA\n"
            )
            afreq.write_text(
                "#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n"
                "1\t1:10:A:G\tA\tG\t0.2\t100\n"
                "1\t1:20:C:T\tC\tT\t0.005\t100\n"
                "1\tdup-a\tG\tA\t0.2\t100\n"
                "1\tdup-b\tG\tA\t0.2\t100\n"
            )
            result = MODULE.read_eligible(pvar, afreq, 0.01)
            self.assertEqual(result, {("1", 10, "A", "G"): ("1:10:A:G", 0.2)})

    def test_dataset_manifest_requires_unique_ids(self):
        with tempfile.TemporaryDirectory() as directory:
            manifest = Path(directory) / "datasets.tsv"
            manifest.write_text(
                "dataset_id\tpgen_prefix\tpriority\n"
                "same\t/a\t1\n"
                "same\t/b\t2\n"
            )
            with self.assertRaisesRegex(ValueError, "unique"):
                MODULE.read_datasets(manifest)


if __name__ == "__main__":
    unittest.main()
