import unittest

import pandas as pd

from standard_morph.preparation import PreparedMorphology
from standard_morph.metrics.single_root import SingleRootNodeMetric
from standard_morph.metrics.soma_first_node import SomaFirstNodeMetric
from standard_morph.metrics.node_identity import NodeIdentityTypesMetric
from standard_morph.metrics.duplicate_coordinates import DuplicateNodeCoordinatesMetric
from standard_morph.models.qc_context import QCContext, Space, MorphologyKind
from standard_morph.models.qc_policy import Policy


def _df(rows):
    return pd.DataFrame(rows, columns=["node_id", "compartment", "x", "y", "z", "r", "parent"])


def _pm(rows):
    return PreparedMorphology.from_dataframe(_df(rows))


IMG = QCContext(space=Space.IMAGE_SPACE, morphology_kind=MorphologyKind.MERGED)
EMPTY_POLICY = Policy("t", {})
NODE_IDENTITY_POLICY = Policy("t", {"node_identity_types": {"allowed_types": [1, 2, 3, 4]}})


class TestSingleRootNode(unittest.TestCase):
    def test_valid_single_root(self):
        pm = _pm([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
            (2, 3, 10.0, 0.0, 0.0, 1.0, 1),
            (3, 2, 20.0, 0.0, 0.0, 1.0, 2),
        ])
        r = SingleRootNodeMetric().evaluate(pm, IMG, EMPTY_POLICY)
        self.assertEqual(r.status, "pass")
        self.assertEqual(r.value, 1)
        self.assertEqual(r.flagged_node_ids, [])

    def test_two_roots_fail(self):
        pm = _pm([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
            (2, 3, 10.0, 0.0, 0.0, 1.0, 1),
            (3, 1, 50.0, 0.0, 0.0, 1.0, -1),  # second root, also a soma
        ])
        r = SingleRootNodeMetric().evaluate(pm, IMG, EMPTY_POLICY)
        self.assertEqual(r.status, "fail")
        self.assertEqual(r.value, 2)
        self.assertEqual(sorted(r.flagged_node_ids), [1, 3])

    def test_extra_type1_node_fail(self):
        pm = _pm([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
            (2, 1, 10.0, 0.0, 0.0, 1.0, 1),  # a second type-1 node (not a root)
        ])
        r = SingleRootNodeMetric().evaluate(pm, IMG, EMPTY_POLICY)
        self.assertEqual(r.status, "fail")
        self.assertEqual(r.measurements["n_type1_nodes"], 2)
        self.assertIn(2, r.flagged_node_ids)

    def test_split_type1_and_root_fail(self):
        # One type-1 node that is NOT the root, and one root that is NOT type 1.
        # Individually n_roots == 1 and n_type1 == 1, but no single node is both,
        # so the soma-root check must fail.
        pm = _pm([
            (1, 3, 0.0, 0.0, 0.0, 1.0, -1),   # root (parent -1) but type 3
            (2, 1, 10.0, 0.0, 0.0, 1.0, 1),   # type 1 (soma) but not the root
        ])
        r = SingleRootNodeMetric().evaluate(pm, IMG, EMPTY_POLICY)
        self.assertEqual(r.status, "fail")
        self.assertEqual(r.measurements["n_soma_roots"], 0)
        self.assertEqual(r.measurements["n_roots"], 1)
        self.assertEqual(r.measurements["n_type1_nodes"], 1)
        self.assertIn("soma root", r.message)
        self.assertEqual(sorted(r.flagged_node_ids), [1, 2])

    def test_root_wrong_id(self):
        # Single soma root, but its node_id is 5 (expected 1). First-line placement
        # is soma_first_node's job, not this metric's.
        pm = _pm([
            (2, 3, 10.0, 0.0, 0.0, 1.0, 5),
            (5, 1, 0.0, 0.0, 0.0, 1.0, -1),
        ])
        r = SingleRootNodeMetric().evaluate(pm, IMG, EMPTY_POLICY)
        self.assertEqual(r.status, "fail")
        self.assertEqual(r.value, 1)  # still exactly one root
        self.assertEqual(r.flagged_node_ids, [5])
        self.assertIn("node_id", r.message)

    def test_wrong_id_no_longer_mentions_first_line(self):
        # single_root_node no longer owns the first-line check.
        pm = _pm([(2, 3, 10.0, 0.0, 0.0, 1.0, 5), (5, 1, 0.0, 0.0, 0.0, 1.0, -1)])
        r = SingleRootNodeMetric().evaluate(pm, IMG, EMPTY_POLICY)
        self.assertNotIn("first line", r.message)

    def test_valid_single_node_soma(self):
        # Just the soma, no children -> still a valid single root.
        pm = _pm([(1, 1, 0.0, 0.0, 0.0, 1.0, -1)])
        r = SingleRootNodeMetric().evaluate(pm, IMG, EMPTY_POLICY)
        self.assertEqual(r.status, "pass")
        self.assertEqual(r.value, 1)
        self.assertEqual(r.measurements["n_soma_roots"], 1)

    def test_no_root_fails(self):
        # A 2-node cycle: no node has parent == -1, so there is no root at all.
        pm = _pm([
            (1, 1, 0.0, 0.0, 0.0, 1.0, 2),
            (2, 3, 10.0, 0.0, 0.0, 1.0, 1),
        ])
        r = SingleRootNodeMetric().evaluate(pm, IMG, EMPTY_POLICY)
        self.assertEqual(r.status, "fail")
        self.assertEqual(r.value, 0)  # n_roots
        self.assertEqual(r.measurements["n_soma_roots"], 0)

    def test_no_soma_fails(self):
        # One topological root, but no type-1 node anywhere.
        pm = _pm([
            (1, 2, 0.0, 0.0, 0.0, 1.0, -1),   # root is axon (type 2), not soma
            (2, 2, 10.0, 0.0, 0.0, 1.0, 1),
        ])
        r = SingleRootNodeMetric().evaluate(pm, IMG, EMPTY_POLICY)
        self.assertEqual(r.status, "fail")
        self.assertEqual(r.measurements["n_roots"], 1)
        self.assertEqual(r.measurements["n_type1_nodes"], 0)
        self.assertEqual(r.measurements["n_soma_roots"], 0)
        self.assertEqual(r.flagged_node_ids, [1])

    def test_extra_non_soma_root_fails(self):
        # A valid soma root PLUS a stray extra root of another type. The
        # soma-root condition itself holds; the extra root is what fails, and
        # each constraint is reported independently.
        pm = _pm([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),   # soma root (valid)
            (2, 2, 10.0, 0.0, 0.0, 1.0, 1),
            (3, 3, 50.0, 0.0, 0.0, 1.0, -1),  # stray extra root, type 3
        ])
        r = SingleRootNodeMetric().evaluate(pm, IMG, EMPTY_POLICY)
        self.assertEqual(r.status, "fail")
        self.assertEqual(r.measurements["n_soma_roots"], 1)  # soma root is fine
        self.assertEqual(r.measurements["n_roots"], 2)       # ...but there are two roots
        self.assertEqual(r.measurements["n_type1_nodes"], 1)
        self.assertEqual(r.value, 2)
        self.assertEqual(sorted(r.flagged_node_ids), [1, 3])


class TestSomaFirstNode(unittest.TestCase):
    def test_soma_first_passes(self):
        pm = _pm([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),   # soma on the first line
            (2, 3, 10.0, 0.0, 0.0, 1.0, 1),
        ])
        r = SomaFirstNodeMetric().evaluate(pm, IMG, EMPTY_POLICY)
        self.assertEqual(r.status, "pass")
        self.assertTrue(r.value)

    def test_soma_not_first_fails(self):
        pm = _pm([
            (2, 3, 10.0, 0.0, 0.0, 1.0, 5),   # a dendrite is the first line
            (5, 1, 0.0, 0.0, 0.0, 1.0, -1),   # soma is the second line
        ])
        r = SomaFirstNodeMetric().evaluate(pm, IMG, EMPTY_POLICY)
        self.assertEqual(r.status, "fail")
        self.assertFalse(r.value)
        self.assertEqual(r.flagged_node_ids, [2])           # the offending first node
        self.assertEqual(r.measurements["first_node_id"], 2)

    def test_first_node_type1_but_not_root_fails(self):
        # First row is type 1 but has a parent -> not a soma root, so not "the soma".
        pm = _pm([
            (1, 1, 0.0, 0.0, 0.0, 1.0, 2),    # type 1 but parented (not a root)
            (2, 1, 0.0, 0.0, 10.0, 1.0, -1),  # the actual soma root
        ])
        r = SomaFirstNodeMetric().evaluate(pm, IMG, EMPTY_POLICY)
        self.assertEqual(r.status, "fail")

    def test_survives_topology_failure(self):
        # requires_topology=False, so a duplicate-id file still runs this check.
        from standard_morph.engine import run_qc
        df = _df([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
            (2, 3, 1.0, 0.0, 0.0, 1.0, 1),
            (2, 3, 9.0, 9.0, 9.0, 1.0, 1),   # duplicate id -> TOPOLOGY failure
        ])
        report = run_qc(df, IMG, metrics=["soma_first_node", "single_root_node"])
        by_name = {r.name: r.status for r in report.results}
        self.assertNotEqual(by_name["soma_first_node"], "skipped")  # coordinate-only
        self.assertEqual(by_name["single_root_node"], "skipped")    # needs topology


class TestNodeIdentityTypes(unittest.TestCase):
    def test_valid_merged(self):
        pm = _pm([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
            (2, 2, 10.0, 0.0, 0.0, 1.0, 1),
            (3, 3, 20.0, 0.0, 0.0, 1.0, 1),
        ])
        r = NodeIdentityTypesMetric().evaluate(pm, IMG, NODE_IDENTITY_POLICY)
        self.assertEqual(r.status, "pass")
        self.assertEqual(r.value, 0)
        self.assertEqual(r.measurements["unique_types"], [1, 2, 3])

    def test_unexpected_type_flagged(self):
        pm = _pm([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
            (2, 2, 10.0, 0.0, 0.0, 1.0, 1),
            (3, 3, 20.0, 0.0, 0.0, 1.0, 1),
            (4, 7, 30.0, 0.0, 0.0, 1.0, 3),  # type 7 is not allowed
        ])
        r = NodeIdentityTypesMetric().evaluate(pm, IMG, NODE_IDENTITY_POLICY)
        self.assertEqual(r.status, "fail")
        self.assertEqual(r.value, 1)
        self.assertEqual(r.flagged_node_ids, [4])

    def test_merged_missing_axon_fail(self):
        pm = _pm([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
            (2, 3, 10.0, 0.0, 0.0, 1.0, 1),  # only soma + dendrite, no axon
        ])
        r = NodeIdentityTypesMetric().evaluate(pm, IMG, NODE_IDENTITY_POLICY)
        self.assertEqual(r.status, "fail")
        self.assertEqual(r.value, 0)  # no out-of-set types, but an identity is missing
        self.assertIn("axon", r.message)


class TestDuplicateNodeCoordinates(unittest.TestCase):
    def test_no_duplicates(self):
        pm = _pm([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
            (2, 3, 10.0, 0.0, 0.0, 1.0, 1),
            (3, 3, 20.0, 0.0, 0.0, 1.0, 2),
        ])
        r = DuplicateNodeCoordinatesMetric().evaluate(pm, IMG, EMPTY_POLICY)
        self.assertEqual(r.status, "pass")
        self.assertEqual(r.value, 0)

    def test_duplicates_flagged(self):
        pm = _pm([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
            (2, 3, 10.0, 10.0, 10.0, 1.0, 1),
            (3, 3, 10.0, 10.0, 10.0, 1.0, 2),  # duplicate of node 2
        ])
        r = DuplicateNodeCoordinatesMetric().evaluate(pm, IMG, EMPTY_POLICY)
        self.assertEqual(r.status, "fail")
        self.assertEqual(r.value, 2)
        self.assertEqual(r.measurements["n_duplicate_groups"], 1)
        self.assertEqual(sorted(r.flagged_node_ids), [2, 3])


if __name__ == "__main__":
    unittest.main()
