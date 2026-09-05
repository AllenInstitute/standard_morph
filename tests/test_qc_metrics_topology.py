import unittest

import pandas as pd

from standard_morph.preparation import PreparedMorphology
from standard_morph.metrics.connected_component import SingleConnectedComponentMetric
from standard_morph.metrics.branch_degree import BranchMaxDegreeMetric
from standard_morph.models.qc_context import QCContext, Space
from standard_morph.models.qc_policy import Policy


def _df(rows):
    return pd.DataFrame(rows, columns=["node_id", "compartment", "x", "y", "z", "r", "parent"])


class TestSingleConnectedComponent(unittest.TestCase):
    def _ctx(self):
        return QCContext(space=Space.IMAGE_SPACE)

    def test_single_tree_passes(self):
        df = _df([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
            (2, 3, 1.0, 0.0, 0.0, 1.0, 1),
            (3, 3, 2.0, 0.0, 0.0, 1.0, 2),
        ])
        pm = PreparedMorphology.from_dataframe(df)
        result = SingleConnectedComponentMetric().evaluate(pm, self._ctx(), Policy("t", {}))
        self.assertEqual(result.status, "pass")
        self.assertEqual(result.measurements["n_components"], 1)

    def test_two_roots_fail(self):
        df = _df([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
            (2, 3, 1.0, 0.0, 0.0, 1.0, 1),
            (3, 1, 5.0, 0.0, 0.0, 1.0, -1),  # second root -> second component
            (4, 3, 6.0, 0.0, 0.0, 1.0, 3),
        ])
        pm = PreparedMorphology.from_dataframe(df)
        result = SingleConnectedComponentMetric().evaluate(pm, self._ctx(), Policy("t", {}))
        self.assertEqual(result.status, "fail")
        self.assertEqual(result.measurements["n_components"], 2)
        self.assertEqual(result.flagged_node_ids, [1, 3])  # the two roots

    def test_orphan_subtree_fails(self):
        df = _df([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
            (2, 3, 1.0, 0.0, 0.0, 1.0, 1),
            (3, 3, 2.0, 0.0, 0.0, 1.0, 999),  # parent missing -> disconnected
        ])
        pm = PreparedMorphology.from_dataframe(df)
        result = SingleConnectedComponentMetric().evaluate(pm, self._ctx(), Policy("t", {}))
        self.assertEqual(result.status, "fail")
        self.assertIn(3, result.flagged_node_ids)


class TestBranchMaxDegree(unittest.TestCase):
    def _policy(self, max_children=2):
        return Policy("t", {"branch_max_degree": {"max_children": max_children}})

    def test_trifurcation_flagged(self):
        # Node 2 has three children (3, 4, 5) -> exceeds max_children=2.
        df = _df([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
            (2, 3, 1.0, 0.0, 0.0, 1.0, 1),
            (3, 3, 2.0, 1.0, 0.0, 1.0, 2),
            (4, 3, 2.0, 0.0, 0.0, 1.0, 2),
            (5, 3, 2.0, -1.0, 0.0, 1.0, 2),
        ])
        pm = PreparedMorphology.from_dataframe(df)
        result = BranchMaxDegreeMetric().evaluate(pm, QCContext(space=Space.IMAGE_SPACE), self._policy())
        # >2 children is unusual but can be real -> flagged for human review.
        self.assertEqual(result.status, "review")
        self.assertEqual(result.flagged_node_ids, [2])
        self.assertEqual(result.measurements["max_children_observed"], 3)

    def test_soma_multifurcation_excluded(self):
        # Soma (node 1) has three children but must NOT be flagged.
        df = _df([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
            (2, 3, 1.0, 1.0, 0.0, 1.0, 1),
            (3, 3, 1.0, 0.0, 0.0, 1.0, 1),
            (4, 2, 1.0, -1.0, 0.0, 1.0, 1),
        ])
        pm = PreparedMorphology.from_dataframe(df)
        result = BranchMaxDegreeMetric().evaluate(pm, QCContext(space=Space.IMAGE_SPACE), self._policy())
        self.assertEqual(result.status, "pass")
        self.assertEqual(result.flagged_node_ids, [])

    def test_stray_type1_node_is_flaggable(self):
        # Only the real soma (type 1 AND root) is exempt. A stray non-root type-1
        # node with too many children IS a branch-degree concern.
        df = _df([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),   # real soma root
            (2, 1, 1.0, 0.0, 0.0, 1.0, 1),    # stray type-1 node with 3 children
            (3, 3, 2.0, 1.0, 0.0, 1.0, 2),
            (4, 3, 2.0, 0.0, 0.0, 1.0, 2),
            (5, 3, 2.0, -1.0, 0.0, 1.0, 2),
        ])
        pm = PreparedMorphology.from_dataframe(df)
        result = BranchMaxDegreeMetric().evaluate(pm, QCContext(space=Space.IMAGE_SPACE), self._policy())
        self.assertEqual(result.status, "review")
        self.assertEqual(result.flagged_node_ids, [2])


if __name__ == "__main__":
    unittest.main()
