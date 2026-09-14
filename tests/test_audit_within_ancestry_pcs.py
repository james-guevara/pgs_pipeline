import importlib.util
import unittest
from pathlib import Path

import pandas as pd


SCRIPT = Path(__file__).parents[1] / "bin" / "audit_within_ancestry_pcs.py"
SPEC = importlib.util.spec_from_file_location("audit_within_ancestry_pcs", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class AuditWithinAncestryPcsTests(unittest.TestCase):
    def test_manifest_aliases_are_combined(self):
        manifest = pd.DataFrame({
            "Genomic_ID": ["wgs1", None],
            "guid": [None, "imp1"],
            "site": ["A", "B"],
        })
        lookup = MODULE.build_manifest_lookup(manifest, ["Genomic_ID", "guid"])
        self.assertEqual(set(lookup.CURRENT_IID), {"wgs1", "imp1"})

    def test_eta_squared_detects_complete_separation(self):
        frame = pd.DataFrame({"pc": [-1.0, -1.0, 1.0, 1.0], "group": ["A", "A", "B", "B"]})
        eta2, count, groups = MODULE.eta_squared(frame, "pc", "group")
        self.assertEqual((count, groups), (4, 2))
        self.assertAlmostEqual(eta2, 1.0)


if __name__ == "__main__":
    unittest.main()
