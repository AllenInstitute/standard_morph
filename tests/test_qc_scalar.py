"""Every scalar-able metric exposes a canonical `value` + `value_label`."""
import unittest

import numpy as np
import pandas as pd

from standard_morph.preparation import PreparedMorphology
from standard_morph.models.qc_context import QCContext, Space
from standard_morph.models.qc_policy import Policy
from standard_morph.metrics.local_tortuosity import LocalTortuosityMetric
from standard_morph.metrics.connected_component import SingleConnectedComponentMetric
from standard_morph.metrics.branch_degree import BranchMaxDegreeMetric
from standard_morph.metrics.ccf_mesh import NodesOutsideCcfMeshMetric, SomaInsideCcfMeshMetric


def _df(rows):
    return pd.DataFrame(rows, columns=["node_id", "compartment", "x", "y", "z", "r", "parent"])


def _atlas():
    a = np.zeros((10, 10, 10), dtype=np.uint32)
    a[2:8, 2:8, 2:8] = 500
    return a


def _ccf_ctx():
    return QCContext(space=Space.CCF_REGISTERED, ccf_resolution=1,
                     resources={"ccf_annotation": _atlas()})


class TestScalarValues(unittest.TestCase):
    def test_tortuosity_value_is_max(self):
        df = _df([
            (1, 3, 0.0, 0.0, 0.0, 1.0, -1), (2, 3, 10.0, 0.0, 0.0, 1.0, 1),
            (3, 3, 20.0, 5.0, 0.0, 1.0, 2), (4, 3, 30.0, 0.0, 0.0, 1.0, 3),
            (5, 3, 40.0, 0.0, 0.0, 1.0, 4),
        ])
        pm = PreparedMorphology.from_dataframe(df)
        r = LocalTortuosityMetric().evaluate(
            pm, QCContext(space=Space.IMAGE_SPACE),
            Policy("t", {"local_tortuosity": {"tortuosity_threshold": 10}}))
        self.assertEqual(r.value_label, "max_tortuosity")
        self.assertEqual(r.value, r.measurements["max_tortuosity"])
        self.assertAlmostEqual(r.value, 1.11803, places=4)

    def test_connected_component_value(self):
        df = _df([(1, 1, 0.0, 0.0, 0.0, 1.0, -1), (2, 3, 1.0, 0.0, 0.0, 1.0, 1),
                  (3, 1, 5.0, 0.0, 0.0, 1.0, -1)])
        pm = PreparedMorphology.from_dataframe(df)
        r = SingleConnectedComponentMetric().evaluate(pm, QCContext(space=Space.IMAGE_SPACE), Policy("t", {}))
        self.assertEqual(r.value_label, "n_components")
        self.assertEqual(r.value, 2)

    def test_branch_degree_value(self):
        df = _df([(1, 1, 0.0, 0.0, 0.0, 1.0, -1), (2, 3, 1.0, 0.0, 0.0, 1.0, 1),
                  (3, 3, 2.0, 1.0, 0.0, 1.0, 2), (4, 3, 2.0, 0.0, 0.0, 1.0, 2),
                  (5, 3, 2.0, -1.0, 0.0, 1.0, 2)])
        pm = PreparedMorphology.from_dataframe(df)
        r = BranchMaxDegreeMetric().evaluate(
            pm, QCContext(space=Space.IMAGE_SPACE),
            Policy("t", {"branch_max_degree": {"max_children": 2}}))
        self.assertEqual(r.value_label, "max_children_observed")
        self.assertEqual(r.value, 3)

    def test_nodes_outside_value_is_fraction(self):
        df = _df([(1, 1, 5.0, 5.0, 5.0, 1.0, -1), (2, 3, 5.0, 5.0, 6.0, 1.0, 1),
                  (3, 3, 0.0, 0.0, 0.0, 1.0, 2), (4, 3, 9.0, 9.0, 9.0, 1.0, 3)])
        pm = PreparedMorphology.from_dataframe(df)
        r = NodesOutsideCcfMeshMetric().evaluate(pm, _ccf_ctx(), Policy("t", {"nodes_outside_ccf_mesh": {"max_fraction_outside": 0.05}}))
        self.assertEqual(r.value_label, "fraction_outside")
        self.assertAlmostEqual(r.value, 0.5)

    def test_binary_metric_leaves_value_none(self):
        df = _df([(1, 1, 5.0, 5.0, 5.0, 1.0, -1), (2, 3, 5.0, 5.0, 6.0, 1.0, 1)])
        pm = PreparedMorphology.from_dataframe(df)
        r = SomaInsideCcfMeshMetric().evaluate(pm, _ccf_ctx(), Policy("t", {}))
        self.assertIsNone(r.value)  # genuinely binary -> relies on status
        self.assertIn(r.status, ("pass", "fail"))


if __name__ == "__main__":
    unittest.main()
