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
        self.assertIn("process EXPLAIN_SCORE_MATCHING", source)
        self.assertIn("variant_match_breakdown.tsv", source)

    def test_vcf_dosage_and_base_maf_are_explicit_and_independent(self):
        entrypoint = (ROOT / "main.nf").read_text()
        source = (ROOT / "workflows" / "pgs.nf").read_text()
        self.assertIn("params.vcf_dosage_field = null", entrypoint)
        self.assertIn("params.base_maf = null", entrypoint)
        self.assertIn('dosage=${params.vcf_dosage_field}', source)
        self.assertIn('--maf ${params.base_maf}', source)
        self.assertNotIn('--maf ${params.maf} \\\n+      --make-pgen \\\n+      --out chr${chr}', source)

    def test_pca_extracts_coordinate_ids_before_renaming_to_panel_ids(self):
        source = (ROOT / "workflows" / "pgs.nf").read_text()
        self.assertEqual(source.count("--out extracted"), 2)
        self.assertEqual(source.count("--pfile extracted"), 2)
        for process_name in ("PREPARE_PCA_PFILE", "PREPARE_PCA_PFILE_DIRECT"):
            block = source.split(f"process {process_name} {{", 1)[1].split("\nprocess ", 1)[0]
            self.assertLess(block.index("--extract"), block.index("--update-name"))
            self.assertLess(block.index("--update-name"), block.index("--ref-allele"))

    def test_within_ancestry_uses_argmax_group_and_fid_iid_keep_files(self):
        classifier = (ROOT / "bin" / "apply_extra_trees.py").read_text()
        within = (ROOT / "bin" / "run_within_ancestry_pca.sh").read_text()
        self.assertIn('"MOST_LIKELY_ANCESTRY": classes[best]', classifier)
        self.assertIn('name == "MOST_LIKELY_ANCESTRY"', within)
        self.assertIn('print "#FID\\tIID"', within)
        self.assertIn('"${pfile}.psam" > "$keep_file"', within)
        self.assertIn('name == "ALT_CTS"', within)
        self.assertIn('name == "OBS_CT"', within)
        self.assertIn('name == "A1"', within)
        self.assertIn('name == "PC1"', within)
        self.assertNotIn('--score "$group_dir/training.eigenvec.allele" 2 5', within)

    def test_analysis_join_resolves_iid_headers_and_argmax_groups(self):
        source = (ROOT / "bin" / "build_analysis_dataset.py").read_text()
        self.assertIn('def sample_id(row):', source)
        self.assertIn('id_candidates = ("IID", "#IID")', source)
        self.assertIn('"MOST_LIKELY_ANCESTRY" in ancestry_fields', source)
        self.assertIn('"ANCESTRY_MOST_LIKELY"', source)


if __name__ == "__main__":
    unittest.main()
