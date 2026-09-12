import csv
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class ExplainScoreMatchingTest(unittest.TestCase):
    def test_categories_account_for_every_weight_row(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "weights.tsv").write_text(
                "SNP\tA1\tBETA\nrs1\tA\t1\nrs2\tG\t1\nrs3\tT\t1\nrs4\tC\t1\nrs5\tA\t1\n"
            )
            (root / "map.tsv").write_text("1:1:A:G\trs1\n1:2:C:G\trs2\n1:3:A:T\trs3\n1:4:A:G\trs4\n")
            (root / "source.pvar").write_text(
                "#CHROM\tPOS\tID\tREF\tALT\n1\t1\t1:1:A:G\tA\tG\n1\t2\t1:2:C:G\tC\tG\n1\t4\t1:4:A:G\tA\tG\n"
            )
            (root / "score.pvar").write_text(
                "#CHROM\tPOS\tID\tREF\tALT\n1\t1\trs1\tA\tG\n1\t4\trs4\tA\tG\n"
            )
            (root / "matched.vars").write_text("rs1\n")
            subprocess.run(
                [
                    "python3", str(ROOT / "bin" / "explain_score_matching.py"),
                    "--trait", "test", "--weights", str(root / "weights.tsv"),
                    "--id-col", "1", "--allele-col", "2",
                    "--source-pvar", str(root / "source.pvar"),
                    "--score-pvar", str(root / "score.pvar"),
                    "--rsid-map", str(root / "map.tsv"),
                    "--matched-vars", str(root / "matched.vars"),
                    "--summary", str(root / "summary.tsv"),
                    "--details", str(root / "details.tsv"),
                ],
                check=True,
            )
            with (root / "summary.tsv").open() as handle:
                observed = {row["category"]: int(row["variants"]) for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(observed["matched"], 1)
            self.assertEqual(observed["removed_by_score_view_maf"], 1)
            self.assertEqual(observed["absent_from_harmonized_source"], 1)
            self.assertEqual(observed["effect_allele_incompatible"], 1)
            self.assertEqual(observed["no_selected_rsid_mapping_or_absent"], 1)
            self.assertEqual(observed["requested_total"], 5)
