import unittest

import pandas as pd

from standard_morph.preparation import PreparedMorphology
from standard_morph.metrics.axon_origination import AxonOriginationMetric
from standard_morph.metrics.apical_origination import ApicalOriginationMetric
from standard_morph.metrics.compartment_transitions import CompartmentTransitionsMetric
from standard_morph.models.qc_context import QCContext, Space, MorphologyKind
from standard_morph.models.qc_policy import Policy


def _df(rows):
    return pd.DataFrame(rows, columns=["node_id", "compartment", "x", "y", "z", "r", "parent"])


def _pm(rows):
    return PreparedMorphology.from_dataframe(_df(rows))


IMG = QCContext(space=Space.IMAGE_SPACE, morphology_kind=MorphologyKind.MERGED)
AXON_POLICY = Policy("t", {"axon_origination": {"max_axon_origin_to_soma_um": 75.0}})
APICAL_POLICY = Policy("t", {"apical_origination": {"max_origins": 1}})
EMPTY = Policy("t", {})


class TestAxonOrigination(unittest.TestCase):
    def test_valid_axon_from_soma_near(self):
        pm = _pm([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),   # soma
            (2, 3, 10.0, 0.0, 0.0, 1.0, 1),   # basal from soma
            (3, 2, 0.0, 0.0, 15.0, 1.0, 1),   # axon from soma, 15 um away
            (4, 2, 0.0, 0.0, 30.0, 1.0, 3),   # internal axon
        ])
        r = AxonOriginationMetric().evaluate(pm, IMG, AXON_POLICY)
        self.assertEqual(r.status, "pass")
        self.assertEqual(r.measurements["n_axon_origins"], 1)
        self.assertAlmostEqual(r.value, 15.0)

    def test_axon_from_distal_basal_fails_on_soma_distance(self):
        # The key case: axon stems from a VALID parent type (basal) but the basal
        # node is 110 um out, so the origin is far from the SOMA -> fail.
        pm = _pm([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),    # soma
            (2, 3, 0.0, 0.0, 100.0, 1.0, 1),   # distal basal dendrite
            (3, 2, 0.0, 0.0, 110.0, 1.0, 2),   # axon from that basal (110 um from soma)
        ])
        r = AxonOriginationMetric().evaluate(pm, IMG, AXON_POLICY)
        self.assertEqual(r.status, "fail")
        self.assertEqual(r.measurements["n_axon_origins"], 1)   # parent type is fine...
        self.assertAlmostEqual(r.value, 110.0)                  # ...but too far from the soma
        self.assertIn(3, r.flagged_node_ids)

    def test_multiple_origins_fail(self):
        pm = _pm([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
            (2, 2, 0.0, 0.0, 5.0, 1.0, 1),    # axon from soma (origin 1)
            (3, 3, 5.0, 0.0, 0.0, 1.0, 1),    # basal
            (4, 2, 6.0, 0.0, 0.0, 1.0, 3),    # axon from basal (origin 2)
        ])
        r = AxonOriginationMetric().evaluate(pm, IMG, AXON_POLICY)
        self.assertEqual(r.status, "fail")
        self.assertEqual(r.measurements["n_axon_origins"], 2)
        self.assertEqual(sorted(r.flagged_node_ids), [2, 4])

    def test_axon_from_apical_is_invalid_parent(self):
        pm = _pm([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
            (2, 4, 0.0, 0.0, 10.0, 1.0, 1),   # apical from soma
            (3, 2, 0.0, 0.0, 15.0, 1.0, 2),   # axon from apical (invalid parent type)
        ])
        r = AxonOriginationMetric().evaluate(pm, IMG, AXON_POLICY)
        self.assertEqual(r.status, "fail")
        self.assertIn("expected soma", r.message)
        self.assertIn(3, r.flagged_node_ids)

    def test_no_axon_passes(self):
        pm = _pm([(1, 1, 0.0, 0.0, 0.0, 1.0, -1), (2, 3, 1.0, 0.0, 0.0, 1.0, 1)])
        r = AxonOriginationMetric().evaluate(pm, IMG, AXON_POLICY)
        self.assertEqual(r.status, "pass")
        self.assertIsNone(r.value)


class TestApicalOrigination(unittest.TestCase):
    def test_single_apical_trunk_passes(self):
        pm = _pm([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),   # soma
            (2, 4, 0.0, 0.0, 10.0, 1.0, 1),   # apical trunk from soma (one origin)
            (3, 4, 0.0, 0.0, 20.0, 1.0, 2),   # apical from apical (internal)
        ])
        r = ApicalOriginationMetric().evaluate(pm, IMG, APICAL_POLICY)
        self.assertEqual(r.status, "pass")
        self.assertEqual(r.value, 1)

    def test_two_apical_trunks_flagged_for_review(self):
        # Multiple apical trunks are biologically plausible, so this is flagged
        # for human review ("review"), not an objective failure.
        pm = _pm([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),   # soma
            (2, 4, 0.0, 0.0, 10.0, 1.0, 1),   # apical trunk 1
            (3, 4, 10.0, 0.0, 0.0, 1.0, 1),   # apical trunk 2
        ])
        r = ApicalOriginationMetric().evaluate(pm, IMG, APICAL_POLICY)
        self.assertEqual(r.status, "review")
        self.assertEqual(r.value, 2)
        self.assertEqual(sorted(r.flagged_node_ids), [2, 3])

    def test_no_apical_passes(self):
        pm = _pm([(1, 1, 0.0, 0.0, 0.0, 1.0, -1), (2, 3, 0.0, 0.0, 10.0, 1.0, 1)])
        r = ApicalOriginationMetric().evaluate(pm, IMG, APICAL_POLICY)
        self.assertEqual(r.status, "pass")
        self.assertEqual(r.value, 0)


class TestCompartmentTransitions(unittest.TestCase):
    def test_valid_tree(self):
        pm = _pm([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),   # soma
            (2, 3, 1.0, 0.0, 0.0, 1.0, 1),    # basal from soma
            (3, 3, 2.0, 0.0, 0.0, 1.0, 2),    # basal from basal
            (4, 4, 0.0, 1.0, 0.0, 1.0, 1),    # apical from soma
            (5, 4, 0.0, 2.0, 0.0, 1.0, 4),    # apical from apical
            (6, 2, 0.0, 0.0, 1.0, 1.0, 1),    # axon (metric 13's job, not flagged here)
        ])
        r = CompartmentTransitionsMetric().evaluate(pm, IMG, EMPTY)
        self.assertEqual(r.status, "pass")
        self.assertEqual(r.value, 0)

    def test_basal_off_axon_flagged(self):
        pm = _pm([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
            (2, 2, 0.0, 0.0, 10.0, 1.0, 1),   # axon
            (3, 3, 0.0, 0.0, 20.0, 1.0, 2),   # basal off the axon (invalid)
        ])
        r = CompartmentTransitionsMetric().evaluate(pm, IMG, EMPTY)
        self.assertEqual(r.status, "fail")
        self.assertEqual(r.flagged_node_ids, [3])

    def test_apical_off_basal_flagged(self):
        pm = _pm([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
            (2, 3, 1.0, 0.0, 0.0, 1.0, 1),    # basal from soma
            (3, 4, 2.0, 0.0, 0.0, 1.0, 2),    # apical from basal (invalid)
        ])
        r = CompartmentTransitionsMetric().evaluate(pm, IMG, EMPTY)
        self.assertEqual(r.status, "fail")
        self.assertEqual(r.flagged_node_ids, [3])

    def test_multiple_basal_trunks_from_soma_ok(self):
        # Several basal dendrites originating from the soma is normal.
        pm = _pm([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
            (2, 3, 1.0, 0.0, 0.0, 1.0, 1),    # basal trunk 1
            (3, 3, -1.0, 0.0, 0.0, 1.0, 1),   # basal trunk 2
            (4, 3, 0.0, 1.0, 0.0, 1.0, 1),    # basal trunk 3
        ])
        r = CompartmentTransitionsMetric().evaluate(pm, IMG, EMPTY)
        self.assertEqual(r.status, "pass")

    def test_dendrite_off_stray_soma_flagged(self):
        # A stray mid-tree soma (type 1, not the root) is not "the soma": a
        # dendrite hanging off it must be flagged. Both the apical (1-4-4-4-1-...)
        # and basal (1-3-3-3-1-...) variants flag node 6.
        for dendrite in (4, 3):
            pm = _pm([
                (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
                (2, dendrite, 0.0, 0.0, 10.0, 1.0, 1),
                (3, dendrite, 0.0, 0.0, 20.0, 1.0, 2),
                (4, dendrite, 0.0, 0.0, 30.0, 1.0, 3),
                (5, 1, 0.0, 0.0, 40.0, 1.0, 4),        # stray soma mid-tree
                (6, dendrite, 0.0, 0.0, 50.0, 1.0, 5), # dendrite off the stray soma
                (7, dendrite, 0.0, 0.0, 60.0, 1.0, 6),
                (8, dendrite, 0.0, 0.0, 70.0, 1.0, 7),
            ])
            r = CompartmentTransitionsMetric().evaluate(pm, IMG, EMPTY)
            self.assertEqual(r.status, "fail", f"dendrite type {dendrite}")
            self.assertEqual(r.flagged_node_ids, [6], f"dendrite type {dendrite}")


if __name__ == "__main__":
    unittest.main()
