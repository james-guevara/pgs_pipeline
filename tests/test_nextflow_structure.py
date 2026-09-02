import shutil
import subprocess
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


@unittest.skipUnless(shutil.which("nextflow"), "Nextflow is not installed")
class NextflowStructureTest(unittest.TestCase):
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
        self.assertIn("workflow ANCESTRY_WORKFLOW", source)
        self.assertIn("ANCESTRY_WORKFLOW(qcPfile, directPfileEnabled)", source)
        self.assertNotIn("${projectDir}/bin/", source)
        self.assertIn("${moduleDir}/../bin/", source)
        for name in (
            "qc_pfile",
            "combined_scores",
            "global_pcs",
            "projection_variants",
            "ancestry_assignments",
            "within_ancestry",
            "analysis_dataset",
            "analysis_dictionary",
        ):
            self.assertIn(f"{name} =", source)

    def test_detachable_ancestry_wrapper_uses_named_workflow(self):
        source = (ROOT / "ancestry.nf").read_text()
        self.assertIn("include { ANCESTRY_WORKFLOW }", source)
        self.assertIn("--input_pfile is required", source)
        self.assertIn("ANCESTRY_WORKFLOW(qcPfile, directPfileEnabled)", source)

        workflow = (ROOT / "workflows" / "pgs.nf").read_text()
        for name in (
            "global_pcs",
            "projection_variants",
            "ancestry_assignments",
            "within_ancestry",
        ):
            self.assertIn(f"{name} =", workflow)


if __name__ == "__main__":
    unittest.main()
