import csv
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class BuildAnalysisDatasetTest(unittest.TestCase):
    def test_mixed_iid_headers_and_uncertain_argmax_group(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            within = root / "within"
            (within / "EUR").mkdir(parents=True)
            (root / "scores.tsv").write_text("IID\ttrait\ns1\t0.25\n")
            (root / "global.tsv").write_text("#FID\tIID\tPC1_AVG\nf1\ts1\t0.10\n")
            (root / "ancestry.tsv").write_text(
                "#IID\tANCESTRY\tMOST_LIKELY_ANCESTRY\tMAX_PROBABILITY\tPROB_EUR\n"
                "s1\tADMIXED_OR_UNCERTAIN\tEUR\t0.75\t0.75\n"
            )
            (within / "status.tsv").write_text(
                "ancestry\tassigned_samples\tunrelated_training_samples\tpruned_variants\tpcs\tstatus\treliability\treason\n"
                "EUR\t1\t1\t10\t1\tcompleted\tunreliable_small_sample\t.\n"
            )
            (within / "EUR" / "pcs.tsv").write_text(
                "#FID\tIID\tPC1_AVG\nf1\ts1\t0.20\n"
            )
            output = root / "analysis.tsv"
            dictionary = root / "dictionary.tsv"
            result = subprocess.run(
                [
                    "python3", str(ROOT / "bin" / "build_analysis_dataset.py"),
                    "--scores", str(root / "scores.tsv"),
                    "--global-pcs", str(root / "global.tsv"),
                    "--ancestry", str(root / "ancestry.tsv"),
                    "--within-dir", str(within),
                    "--num-global-pcs", "1",
                    "--num-within-pcs", "1",
                    "--output", str(output),
                    "--dictionary", str(dictionary),
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with output.open(newline="") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["IID"], "s1")
            self.assertEqual(row["ANCESTRY"], "ADMIXED_OR_UNCERTAIN")
            self.assertEqual(row["ANCESTRY_MOST_LIKELY"], "EUR")
            self.assertEqual(row["GLOBAL_PC1"], "0.10")
            self.assertEqual(row["WITHIN_PC1"], "0.20")


if __name__ == "__main__":
    unittest.main()
