"""Integration smoke tests: build the QC layer from real SWC files on disk.

These complement the synthetic unit tests. Synthetic tests assert exact,
hand-computed values; these assert that PreparedMorphology and the metrics
survive contact with real reconstructions (thousands of nodes, real branching,
CCF-registered coordinates, Horta OFFSET headers).
"""
import unittest
from pathlib import Path

from standard_morph.swc_io import read_swc
from standard_morph.preparation import PreparedMorphology
from standard_morph.metrics.local_tortuosity import LocalTortuosityMetric
from standard_morph.models.qc_context import QCContext, Space, MorphologyKind
from standard_morph.policies import get_policy

SWC_DIR = Path(__file__).parent / "swcs"


def _prepared_from(filename):
    df = read_swc(str(SWC_DIR / filename))
    return PreparedMorphology.from_dataframe(df), len(df)


class TestQCIntegration(unittest.TestCase):
    def test_prepared_morphology_structure(self):
        for filename in (
            "N024-648434-CONSENSUS.swc",
            "17109_6601-X5417-Y25287_reg.swc",
            "test_horta_swc.swc",
        ):
            with self.subTest(filename=filename):
                pm, n_rows = _prepared_from(filename)
                self.assertEqual(pm.n, n_rows)
                self.assertGreaterEqual(pm.roots.size, 1)
                self.assertEqual(pm.orphans.size, 0)
                # Every node except roots/orphans should sit in some segment;
                # at minimum a valid tree yields at least one segment.
                self.assertGreater(len(pm.segments), 0)
                # child_counts must reconcile with total non-root edges.
                self.assertEqual(int(pm.child_counts.sum()), n_rows - pm.roots.size)

    def test_tortuosity_runs_on_image_space_file(self):
        pm, _ = _prepared_from("N024-648434-CONSENSUS.swc")
        ctx = QCContext(space=Space.IMAGE_SPACE, morphology_kind=MorphologyKind.MERGED)
        result = LocalTortuosityMetric().evaluate(pm, ctx, get_policy("policy_v1"))

        self.assertIn(result.status, ("pass", "fail"))
        self.assertGreater(result.measurements["n_evaluated"], 0)
        # Full flagged lists stay aligned (no sampling/capping per spec).
        self.assertEqual(len(result.flagged_node_ids), len(result.flagged_node_coordinates))
        self.assertEqual(len(result.flagged_node_ids), result.counts["n_flagged"])

    def test_tortuosity_runs_on_registered_file(self):
        # Tortuosity is applicable in both spaces; it should run on the
        # CCF-registered file and produce a sane result over real branching.
        pm, _ = _prepared_from("17109_6601-X5417-Y25287_reg.swc")
        ctx = QCContext(space=Space.CCF_REGISTERED, morphology_kind=MorphologyKind.MERGED)
        self.assertTrue(LocalTortuosityMetric().is_applicable(ctx))

        result = LocalTortuosityMetric().evaluate(pm, ctx, get_policy("policy_v1"))
        self.assertIn(result.status, ("pass", "fail"))
        self.assertGreater(result.measurements["n_evaluated"], 0)
        self.assertEqual(len(result.flagged_node_ids), result.counts["n_flagged"])


if __name__ == "__main__":
    unittest.main()
