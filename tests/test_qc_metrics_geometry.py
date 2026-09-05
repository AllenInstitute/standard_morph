import unittest

import pandas as pd

from standard_morph.preparation import PreparedMorphology
from standard_morph.metrics.edge_length import EdgeLengthMetric
from standard_morph.metrics.soma_child_distance import SomaChildDistanceMetric
from standard_morph.models.qc_context import QCContext, Space, MorphologyKind
from standard_morph.models.qc_policy import Policy, PolicyRange


def _df(rows):
    return pd.DataFrame(rows, columns=["node_id", "compartment", "x", "y", "z", "r", "parent"])


def _pm(rows):
    return PreparedMorphology.from_dataframe(_df(rows))


IMG = QCContext(space=Space.IMAGE_SPACE, morphology_kind=MorphologyKind.MERGED)
CCF = QCContext(space=Space.CCF_REGISTERED, morphology_kind=MorphologyKind.MERGED)


def _edge_policy(image=PolicyRange(lo=0.0, hi=30.0), ccf=PolicyRange(lo=0.0, hi=10.0)):
    return Policy("t", {"edge_length": {"max_length_um": {"image_space": image, "ccf_registered": ccf}}})


def _soma_policy(threshold=50.0):
    return Policy("t", {"soma_child_distance": {"max_soma_child_to_soma_um": threshold}})


# A chain along z: soma -> node2 (soma child, 5 um) -> node3 (edge 20) -> node4 (edge 35)
def _chain():
    return _pm([
        (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
        (2, 3, 0.0, 0.0, 5.0, 1.0, 1),    # soma's child (excluded from edge_length)
        (3, 3, 0.0, 0.0, 25.0, 1.0, 2),   # edge 2->3 = 20 um
        (4, 3, 0.0, 0.0, 60.0, 1.0, 3),   # edge 3->4 = 35 um
    ])


class TestEdgeLength(unittest.TestCase):
    def test_excludes_soma_children_and_flags_long_edge(self):
        r = EdgeLengthMetric().evaluate(_chain(), IMG, _edge_policy())
        # the soma->node2 edge is excluded; only 2->3 (20) and 3->4 (35) counted
        self.assertEqual(r.measurements["n_edges"], 2)
        self.assertAlmostEqual(r.measurements["max_edge_length"], 35.0)
        self.assertEqual(r.status, "fail")            # 35 > 30 um
        self.assertEqual(r.flagged_node_ids, [4])

    def test_threshold_keyed_by_space(self):
        pm, pol = _chain(), _edge_policy()
        img = EdgeLengthMetric().evaluate(pm, IMG, pol)   # hi=30 -> only 35 over
        ccf = EdgeLengthMetric().evaluate(pm, CCF, pol)   # hi=10 -> both 20 and 35 over
        self.assertEqual(img.counts["n_over"], 1)
        self.assertEqual(ccf.counts["n_over"], 2)
        self.assertEqual(img.thresholds_used["max_length_um"], {"lo": 0.0, "hi": 30.0})
        self.assertEqual(ccf.thresholds_used["max_length_um"], {"lo": 0.0, "hi": 10.0})

    def test_passes_when_all_within(self):
        r = EdgeLengthMetric().evaluate(_chain(), IMG, _edge_policy(image=PolicyRange(lo=0.0, hi=100.0)))
        self.assertEqual(r.status, "pass")
        self.assertEqual(r.flagged_node_ids, [])

    def test_flags_edge_below_lower_bound(self):
        # Near-zero edge (0.5 um) violates the lower bound of a range.
        pm = _pm([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
            (2, 3, 0.0, 0.0, 5.0, 1.0, 1),   # soma child (excluded)
            (3, 3, 0.0, 0.0, 5.5, 1.0, 2),   # edge 2->3 = 0.5 um  (below lo=1.0)
            (4, 3, 0.0, 0.0, 15.5, 1.0, 3),  # edge 3->4 = 10 um   (within range)
        ])
        pol = _edge_policy(image=PolicyRange(lo=1.0, hi=30.0))
        r = EdgeLengthMetric().evaluate(pm, IMG, pol)
        self.assertEqual(r.status, "fail")
        self.assertEqual(r.flagged_node_ids, [3])
        self.assertEqual(r.counts["n_over"], 1)

    def test_scalar_threshold_backward_compat(self):
        pol = Policy("t", {"edge_length": {"max_length_um": {"image_space": 30.0, "ccf_registered": 10.0}}})
        r = EdgeLengthMetric().evaluate(_chain(), IMG, pol)
        self.assertEqual(r.status, "fail")
        self.assertEqual(r.flagged_node_ids, [4])
        self.assertEqual(r.thresholds_used["max_length_um"], 30.0)
        self.assertEqual(r.value_label, "fraction_edges_over_threshold")

    def test_stray_type1_node_is_not_treated_as_soma(self):
        # Node 3 is type 1 but NOT a root -- a stray soma-typed node, not the soma.
        # Its child (node 4) must be measured as a normal edge, not excluded.
        pm = _pm([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),    # the real soma root
            (2, 3, 0.0, 0.0, 5.0, 1.0, 1),     # soma's child (excluded)
            (3, 1, 0.0, 0.0, 25.0, 1.0, 2),    # stray type-1 node (has a parent)
            (4, 3, 0.0, 0.0, 60.0, 1.0, 3),    # child of the stray -> a real 35 um edge
        ])
        r = EdgeLengthMetric().evaluate(pm, IMG, _edge_policy())
        self.assertEqual(r.measurements["n_edges"], 2)   # only soma->node2 is excluded
        self.assertEqual(r.flagged_node_ids, [4])         # the stray's child is measured


class TestSomaChildDistance(unittest.TestCase):
    def test_measures_only_soma_children(self):
        # node2 (5 um) and node3 (60 um) are soma children; node4 is not.
        pm = _pm([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
            (2, 3, 0.0, 0.0, 5.0, 1.0, 1),    # soma child, 5 um
            (3, 3, 0.0, 0.0, 60.0, 1.0, 1),   # soma child, 60 um (far)
            (4, 3, 0.0, 0.0, 62.0, 1.0, 3),   # child of node3, not the soma
        ])
        r = SomaChildDistanceMetric().evaluate(pm, IMG, _soma_policy(50.0))
        self.assertEqual(r.measurements["n_soma_children"], 2)
        self.assertEqual(r.status, "fail")
        self.assertEqual(r.flagged_node_ids, [3])
        self.assertAlmostEqual(r.value, 60.0)

    def test_passes_within_threshold(self):
        pm = _pm([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
            (2, 3, 0.0, 0.0, 5.0, 1.0, 1),
            (3, 2, 0.0, 0.0, 8.0, 1.0, 1),
        ])
        r = SomaChildDistanceMetric().evaluate(pm, IMG, _soma_policy(50.0))
        self.assertEqual(r.status, "pass")
        self.assertAlmostEqual(r.value, 8.0)

    def test_no_soma_children(self):
        pm = _pm([(1, 2, 0.0, 0.0, 0.0, 1.0, -1), (2, 2, 0.0, 0.0, 5.0, 1.0, 1)])  # no soma
        r = SomaChildDistanceMetric().evaluate(pm, IMG, _soma_policy())
        self.assertEqual(r.status, "pass")
        self.assertIsNone(r.value)


if __name__ == "__main__":
    unittest.main()
