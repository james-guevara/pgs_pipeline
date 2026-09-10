import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASH_HAS_MAPFILE = subprocess.run(
    ["bash", "-c", "type mapfile >/dev/null 2>&1"], check=False
).returncode == 0


class CollateScoresTest(unittest.TestCase):
    @unittest.skipUnless(BASH_HAS_MAPFILE, "collator requires Bash 4 mapfile")
    def test_accepts_plink_score_header_with_fid_before_iid(self):
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            (work / "alpha.sscore").write_text(
                "#FID\tIID\tALLELE_CT\tSCORE1_AVG\n"
                "fam1\tsample1\t2\t0.25\n"
                "0\tsample2\t2\t-0.50\n"
            )
            (work / "beta.sscore").write_text(
                "#FID\tIID\tALLELE_CT\tSCORE1_AVG\n"
                "fam1\tsample1\t2\t1.25\n"
                "0\tsample2\t2\t1.50\n"
            )
            for trait in ("alpha", "beta"):
                (work / f"{trait}.score_qc.tsv").write_text(
                    "trait\trequested_variants\tmatched_variants\n"
                    f"{trait}\t10\t9\n"
                )

            result = subprocess.run(
                [
                    "bash",
                    str(ROOT / "bin" / "collate_scores.sh"),
                    "combined.tsv",
                    "qc.tsv",
                    "alpha.sscore",
                    "beta.sscore",
                    "--",
                    "alpha.score_qc.tsv",
                    "beta.score_qc.tsv",
                ],
                cwd=work,
                text=True,
                capture_output=True,
                check=False,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(
                (work / "combined.tsv").read_text(),
                "IID\talpha\tbeta\n"
                "sample1\t0.25\t1.25\n"
                "sample2\t-0.50\t1.50\n",
            )


if __name__ == "__main__":
    unittest.main()
