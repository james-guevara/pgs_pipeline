import shutil
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(shutil.which("nextflow"), "Nextflow is not installed")
class NextflowStructureTest(unittest.TestCase):
    def test_harmonization_entrypoint_compiles_and_exposes_compact_contract(self):
        source = (ROOT / "workflows" / "harmonization.nf").read_text()
        entrypoint = (ROOT / "harmonize.nf").read_text()
        for output in ("harmonization_summary.json", "marker_qc.parquet", "sample_qc.tsv"):
            self.assertIn(output, source)
        self.assertIn("workflow PGEN_HARMONIZATION_WORKFLOW", source)
        self.assertIn("--input_manifest is required", entrypoint)

    def test_standalone_wrapper_compiles_with_reusable_workflow(self):
        result = subprocess.run(
            [
                "nextflow",
                "run",
                "main.nf",
                "-preview",
                "-ansi-log",
                "false",
                "--input_pfile",
                "/tmp/nonexistent-pgs-structure-fixture",
                "--direct_inputs",
                "true",
                "--run_summary_qc",
                "false",
            ],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertIn("* PREVIEW *", result.stdout + result.stderr)

    def test_named_workflow_exposes_composition_outputs(self):
        source = (ROOT / "workflows" / "pgs.nf").read_text()
        self.assertIn("workflow PGS_WORKFLOW", source)
        for name in (
            "qc_pfile",
            "score_pfile",
            "score_input_summary",
            "combined_scores",
            "global_pcs",
            "ancestry_assignments",
            "within_ancestry",
            "analysis_dataset",
            "analysis_dictionary",
        ):
            self.assertIn(f"{name} =", source)

    def test_scoring_always_uses_a_dedicated_maf_filtered_view(self):
        source = (ROOT / "workflows" / "pgs.nf").read_text()
        self.assertIn("process PREPARE_SCORE_PFILE", source)
        self.assertIn("process PREPARE_SCORE_PFILE_DIRECT", source)
        self.assertIn("--maf ${params.maf} --make-pgen", source)
        self.assertIn("weights.combine(scorePfile)", source)
        self.assertNotIn("weights.combine(qcPfile)", source)
        self.assertIn("params.score_rsid_map", source)
        self.assertIn("--update-name '${rsid_map}' 2 1", source)
        self.assertNotIn("--update-name '${rsid_map}' 1 2", source)
        self.assertIn("{print ${'$'}2}' '${rsid_map}' > mapped_rsid_ids.txt", source)
        self.assertNotIn("mapped_coordinate_ids.txt", source)


if __name__ == "__main__":
    unittest.main()
