import unittest

import numpy as np
import pandas as pd

from standard_morph.engine import run_qc
from standard_morph.metrics.integrity import (
    ValidParentReferencesMetric,
    AcyclicMetric,
    ParentBeforeChildMetric,
    CastableColumnsMetric,
    _find_cyclic_nodes,
)
from standard_morph.models.qc_context import QCContext, Space
from standard_morph.models.qc_policy import Policy

EMPTY = Policy("t", {})
IMG = QCContext(space=Space.IMAGE_SPACE)


def _df(rows):
    return pd.DataFrame(rows, columns=["node_id", "compartment", "x", "y", "z", "r", "parent"])


class TestValidParentReferences(unittest.TestCase):
    def test_clean_passes(self):
        df = _df([(1, 1, 0., 0., 0., 1., -1), (2, 3, 0., 0., 10., 1., 1)])
        r = ValidParentReferencesMetric().evaluate(df, IMG, EMPTY)
        self.assertEqual(r.status, "pass")
        self.assertEqual(r.value, 0)

    def test_dangling_parent_flagged(self):
        # node 3's parent (99) does not exist in the file.
        df = _df([(1, 1, 0., 0., 0., 1., -1), (2, 3, 0., 0., 10., 1., 1), (3, 3, 0., 0., 20., 1., 99)])
        r = ValidParentReferencesMetric().evaluate(df, IMG, EMPTY)
        self.assertEqual(r.status, "fail")
        self.assertEqual(r.value, 1)
        self.assertEqual(r.flagged_node_ids, [3])
        self.assertIn("99", r.message)


class TestAcyclic(unittest.TestCase):
    def test_clean_passes(self):
        df = _df([(1, 1, 0., 0., 0., 1., -1), (2, 3, 0., 0., 10., 1., 1), (3, 3, 0., 0., 20., 1., 2)])
        r = AcyclicMetric().evaluate(df, IMG, EMPTY)
        self.assertEqual(r.status, "pass")

    def test_self_parent_is_a_cycle(self):
        df = _df([(1, 1, 0., 0., 0., 1., -1), (2, 3, 0., 0., 10., 1., 2)])  # node 2 parents itself
        r = AcyclicMetric().evaluate(df, IMG, EMPTY)
        self.assertEqual(r.status, "fail")
        self.assertIn(2, r.flagged_node_ids)

    def test_two_node_cycle(self):
        # 2 -> 3 -> 2, with no root at all.
        df = _df([(2, 3, 0., 0., 0., 1., 3), (3, 3, 0., 0., 10., 1., 2)])
        r = AcyclicMetric().evaluate(df, IMG, EMPTY)
        self.assertEqual(r.status, "fail")
        self.assertEqual(sorted(r.flagged_node_ids), [2, 3])

    def test_tail_into_cycle_flags_only_the_cycle(self):
        # 1(root) -> ... but 3 -> 4 -> 5 -> 3 is a cycle; 6 hangs off the cycle.
        df = _df([
            (1, 1, 0., 0., 0., 1., -1),
            (3, 3, 0., 0., 0., 1., 5),
            (4, 3, 0., 0., 0., 1., 3),
            (5, 3, 0., 0., 0., 1., 4),
            (6, 3, 0., 0., 0., 1., 5),  # leads into the cycle but is not on it
        ])
        r = AcyclicMetric().evaluate(df, IMG, EMPTY)
        self.assertEqual(r.status, "fail")
        self.assertEqual(sorted(r.flagged_node_ids), [3, 4, 5])  # not 6, not 1

    def test_helper_directly(self):
        # a clean tree plus a disjoint 3-cycle
        self.assertEqual(_find_cyclic_nodes({1: -1, 2: 1, 3: 2}), set())
        self.assertEqual(_find_cyclic_nodes({7: 8, 8: 9, 9: 7}), {7, 8, 9})


class TestParentBeforeChild(unittest.TestCase):
    def test_in_order_passes(self):
        df = _df([(1, 1, 0., 0., 0., 1., -1), (2, 3, 0., 0., 10., 1., 1), (3, 3, 0., 0., 20., 1., 2)])
        r = ParentBeforeChildMetric().evaluate(df, IMG, EMPTY)
        self.assertEqual(r.status, "pass")

    def test_child_before_parent_flagged(self):
        # node 2 (row 0) references parent 3, which is written later (row 1).
        df = _df([(2, 3, 0., 0., 10., 1., 3), (3, 3, 0., 0., 20., 1., 1), (1, 1, 0., 0., 0., 1., -1)])
        r = ParentBeforeChildMetric().evaluate(df, IMG, EMPTY)
        self.assertEqual(r.status, "fail")
        self.assertIn(2, r.flagged_node_ids)

    def test_topological_but_not_bfs_or_dfs_passes(self):
        # The README counterexample: 1->2, 1->3, 2->4, 4->6, 3->5 in order 1,2,4,3,6,5.
        # Every parent precedes its child, though the order is neither BFS nor DFS.
        df = _df([
            (1, 1, 0., 0., 0., 1., -1),
            (2, 3, 0., 0., 0., 1., 1),
            (4, 3, 0., 0., 0., 1., 2),
            (3, 3, 0., 0., 0., 1., 1),
            (6, 3, 0., 0., 0., 1., 4),
            (5, 3, 0., 0., 0., 1., 3),
        ])
        r = ParentBeforeChildMetric().evaluate(df, IMG, EMPTY)
        self.assertEqual(r.status, "pass")


class TestCastableInf(unittest.TestCase):
    def test_inf_coordinate_flagged(self):
        df = _df([(1, 1, 0., 0., 0., 1., -1), (2, 3, np.inf, 0., 10., 1., 1)])
        r = CastableColumnsMetric().evaluate(df, IMG, EMPTY)
        self.assertEqual(r.status, "fail")
        self.assertEqual(r.measurements["problems_by_column"].get("x"), 1)


class TestAlwaysRunAndBlocking(unittest.TestCase):
    def test_new_checks_run_on_every_call(self):
        df = _df([(1, 1, 0., 0., 0., 1., -1), (2, 3, 0., 0., 10., 1., 1)])
        report = run_qc(df, IMG, metrics=["node_identity_types"])
        names = [r.name for r in report.integrity_results]
        self.assertIn("valid_parent_references", names)
        self.assertIn("acyclic", names)

    def test_cycle_blocks_topology_but_not_coordinate_metrics(self):
        # A self-parent cycle at node 2; single_root_node needs topology (skipped),
        # node_identity_types is coordinate/attribute-only (still runs).
        df = _df([(1, 1, 0., 0., 0., 1., -1), (2, 3, 0., 0., 10., 1., 2)])
        report = run_qc(df, IMG, metrics=["single_root_node", "node_identity_types"])
        self.assertFalse(report.integrity_ok)
        statuses = {r.name: r.status for r in report.results}
        self.assertEqual(statuses["single_root_node"], "skipped")
        self.assertNotEqual(statuses["node_identity_types"], "skipped")


if __name__ == "__main__":
    unittest.main()
