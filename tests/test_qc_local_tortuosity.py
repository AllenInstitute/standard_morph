import unittest

import pandas as pd

from standard_morph.preparation import PreparedMorphology
from standard_morph.metrics.local_tortuosity import LocalTortuosityMetric
from standard_morph.models.qc_context import QCContext, Space, MorphologyKind
from standard_morph.models.qc_policy import Policy


def _df(rows):
    return pd.DataFrame(rows, columns=["node_id", "compartment", "x", "y", "z", "r", "parent"])


def _policy(threshold):
    return Policy(version="test", thresholds={"local_tortuosity": {"tortuosity_threshold": threshold}})


class TestLocalTortuosity(unittest.TestCase):
    def _kinked_path(self):
        # Straight along x except node 3, which sidesteps to y=5.
        # Hand-computed local tortuosity (path/chord over a 3-node window):
        #   node 2 ~ 1.027, node 3 ~ 1.118, node 4 ~ 1.027
        return _df([
            (1, 3, 0.0, 0.0, 0.0, 1.0, -1),
            (2, 3, 10.0, 0.0, 0.0, 1.0, 1),
            (3, 3, 20.0, 5.0, 0.0, 1.0, 2),
            (4, 3, 30.0, 0.0, 0.0, 1.0, 3),
            (5, 3, 40.0, 0.0, 0.0, 1.0, 4),
        ])

    def _ctx(self):
        return QCContext(space=Space.IMAGE_SPACE, morphology_kind=MorphologyKind.MERGED)

    def test_flags_kink_at_tight_threshold(self):
        pm = PreparedMorphology.from_dataframe(self._kinked_path())
        result = LocalTortuosityMetric().evaluate(pm, self._ctx(), _policy(1.1))
        self.assertEqual(result.status, "fail")
        self.assertEqual(result.flagged_node_ids, [3])
        self.assertEqual(result.counts["n_evaluated"], 3)  # interior nodes 2, 3, 4
        self.assertAlmostEqual(result.measurements["max_tortuosity"], 1.11803, places=4)

    def test_passes_at_loose_threshold(self):
        pm = PreparedMorphology.from_dataframe(self._kinked_path())
        result = LocalTortuosityMetric().evaluate(pm, self._ctx(), _policy(10.0))
        self.assertEqual(result.status, "pass")
        self.assertEqual(result.flagged_node_ids, [])

    def test_hairpin_is_infinite_and_flagged(self):
        # Node 2 returns to node 1's coordinate -> zero chord -> inf tortuosity.
        df = _df([
            (1, 3, 0.0, 0.0, 0.0, 1.0, -1),
            (2, 3, 10.0, 0.0, 0.0, 1.0, 1),
            (3, 3, 0.0, 0.0, 0.0, 1.0, 2),
        ])
        pm = PreparedMorphology.from_dataframe(df)
        result = LocalTortuosityMetric().evaluate(pm, self._ctx(), _policy(10.0))
        self.assertEqual(result.flagged_node_ids, [2])

    def test_branch_node_averages_over_children(self):
        # Node 3 is a branch (children 4 and 5); its parent is node 2.
        #   child 4: window (10,0,0)-(20,0,0)-(30,10,0) -> tort ~ 1.0797
        #   child 5: window (10,0,0)-(20,0,0)-(25,0,0)  -> tort = 1.0
        #   branch tortuosity = mean(1.0797, 1.0) ~ 1.0398
        df = _df([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),    # root/soma -> skipped (no parent)
            (2, 3, 10.0, 0.0, 0.0, 1.0, 1),    # reducible -> tort 1.0
            (3, 3, 20.0, 0.0, 0.0, 1.0, 2),    # branch -> averaged tortuosity
            (4, 3, 30.0, 10.0, 0.0, 1.0, 3),   # tip -> skipped
            (5, 3, 25.0, 0.0, 0.0, 1.0, 3),    # tip -> skipped
        ])
        pm = PreparedMorphology.from_dataframe(df)

        # Node 2 (reducible) and node 3 (branch) are measured; tips/root are not.
        flagged = LocalTortuosityMetric().evaluate(pm, self._ctx(), _policy(1.02))
        self.assertEqual(flagged.counts["n_evaluated"], 2)
        self.assertEqual(flagged.status, "fail")
        self.assertEqual(flagged.flagged_node_ids, [3])
        self.assertAlmostEqual(flagged.measurements["max_tortuosity"], 1.03985, places=4)

        # Averaging (not max) over children keeps node 3 under a 1.05 threshold;
        # a max-over-children rule (~1.0797) would wrongly flag it here.
        passed = LocalTortuosityMetric().evaluate(pm, self._ctx(), _policy(1.05))
        self.assertEqual(passed.status, "pass")

    def test_applicable_in_both_spaces(self):
        # Applicability was widened to both coordinate spaces.
        metric = LocalTortuosityMetric()
        self.assertTrue(metric.is_applicable(QCContext(space=Space.IMAGE_SPACE)))
        self.assertTrue(metric.is_applicable(QCContext(space=Space.CCF_REGISTERED)))


if __name__ == "__main__":
    unittest.main()
