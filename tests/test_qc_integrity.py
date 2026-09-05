import unittest

import pandas as pd

from standard_morph import run_qc, QCContext, Space, MorphologyKind
from standard_morph.metrics.integrity import (
    RequiredColumnsMetric,
    NonEmptyMetric,
    CastableColumnsMetric,
    UniqueNodeIdsMetric,
)
from standard_morph.metrics.base import EvaluationPhase, BlockScope
from standard_morph.models.qc_policy import Policy
import numpy as np

COLS = ["node_id", "compartment", "x", "y", "z", "r", "parent"]
IMG = QCContext(space=Space.IMAGE_SPACE, morphology_kind=MorphologyKind.MERGED)
EMPTY_POLICY = Policy("t", {})


def _df(rows, cols=COLS):
    return pd.DataFrame(rows, columns=cols)


def _clean_tree():
    return _df([
        (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
        (2, 2, 10.0, 0.0, 0.0, 1.0, 1),
        (3, 3, 5.0, 5.0, 0.0, 1.0, 1),
    ])


# --------------------------------------------------------------- unit tests
class TestIntegrityMetrics(unittest.TestCase):
    def test_required_columns_pass_and_phase(self):
        m = RequiredColumnsMetric()
        self.assertEqual(m.evaluation_phase, EvaluationPhase.INPUT_INTEGRITY)
        self.assertEqual(m.blocks_on_failure, BlockScope.BUILD)
        r = m.evaluate(_clean_tree(), IMG, EMPTY_POLICY)
        self.assertEqual(r.status, "pass")
        self.assertEqual(r.value, 0)

    def test_required_columns_missing_radius_fails(self):
        df = _df([(1, 1, 0.0, 0.0, 0.0, -1)],
                 cols=["node_id", "compartment", "x", "y", "z", "parent"])  # no 'r'
        r = RequiredColumnsMetric().evaluate(df, IMG, EMPTY_POLICY)
        self.assertEqual(r.status, "fail")
        self.assertEqual(r.measurements["missing_columns"], ["r"])

    def test_non_empty_fails_on_empty(self):
        r = NonEmptyMetric().evaluate(_df([]), IMG, EMPTY_POLICY)
        self.assertEqual(r.status, "fail")
        self.assertEqual(r.value, 0)

    def test_unique_node_ids_flags_duplicates(self):
        df = _df([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
            (2, 3, 1.0, 0.0, 0.0, 1.0, 1),
            (2, 3, 9.0, 9.0, 9.0, 1.0, 1),  # duplicate id 2
        ])
        r = UniqueNodeIdsMetric().evaluate(df, IMG, EMPTY_POLICY)
        self.assertEqual(r.status, "fail")
        self.assertEqual(r.value, 1)
        self.assertEqual(r.flagged_node_ids, [2])

    def test_unique_node_ids_pass(self):
        r = UniqueNodeIdsMetric().evaluate(_clean_tree(), IMG, EMPTY_POLICY)
        self.assertEqual(r.status, "pass")

    def test_castable_columns_pass(self):
        m = CastableColumnsMetric()
        self.assertEqual(m.blocks_on_failure, BlockScope.BUILD)
        r = m.evaluate(_clean_tree(), IMG, EMPTY_POLICY)
        self.assertEqual(r.status, "pass")
        self.assertEqual(r.value, 0)

    def test_castable_columns_flags_non_numeric_coordinate(self):
        df = _df([
            (1, 1, "abc", 0.0, 0.0, 1.0, -1),  # non-numeric x
            (2, 3, 1.0, 0.0, 0.0, 1.0, 1),
        ])
        r = CastableColumnsMetric().evaluate(df, IMG, EMPTY_POLICY)
        self.assertEqual(r.status, "fail")
        self.assertIn("x", r.measurements["problems_by_column"])

    def test_castable_columns_flags_null_int(self):
        df = _df([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
            (np.nan, 3, 1.0, 0.0, 0.0, 1.0, 1),  # null node_id
        ])
        r = CastableColumnsMetric().evaluate(df, IMG, EMPTY_POLICY)
        self.assertEqual(r.status, "fail")
        self.assertIn("node_id", r.measurements["problems_by_column"])

    def test_castable_columns_flags_non_integer_id(self):
        df = _df([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
            (2.5, 3, 1.0, 0.0, 0.0, 1.0, 1),  # non-integral node_id
        ])
        r = CastableColumnsMetric().evaluate(df, IMG, EMPTY_POLICY)
        self.assertEqual(r.status, "fail")
        self.assertIn("node_id", r.measurements["problems_by_column"])


# ------------------------------------------------------ engine two-phase flow
class TestTwoPhaseEngine(unittest.TestCase):
    def test_clean_input_runs_both_phases(self):
        report = run_qc(_clean_tree(), IMG, suite_name="default_pre_registration_tests")
        self.assertTrue(report.morphology_evaluated)
        self.assertTrue(report.integrity_ok)
        # buildability always runs
        self.assertEqual(
            [r.name for r in report.integrity_results],
            ["required_columns", "non_empty", "castable_columns", "unique_node_ids",
             "valid_parent_references", "acyclic"],
        )
        self.assertTrue(all(r.status == "pass" for r in report.integrity_results))
        self.assertGreater(len(report.results), 0)
        self.assertNotIn("skipped", [r.status for r in report.results])

    def test_topology_failure_skips_only_topology_metrics(self):
        # Duplicate id -> TOPOLOGY scope. Coordinate/attribute-only metrics still
        # run; only the topology-dependent one is skipped.
        df = _df([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
            (2, 3, 1.0, 0.0, 0.0, 1.0, 1),
            (2, 3, 9.0, 9.0, 9.0, 1.0, 1),  # duplicate id 2
        ])
        report = run_qc(df, IMG, metrics=[
            "node_identity_types",          # requires_topology=False -> runs
            "duplicate_node_coordinates",   # requires_topology=False -> runs
            "single_connected_component",   # requires_topology=True  -> skipped
        ])
        self.assertTrue(report.morphology_evaluated)   # built fine; only topology suspect
        self.assertFalse(report.integrity_ok)          # unique_node_ids failed
        self.assertFalse(report.passed)

        by_name = {r.name: r for r in report.results}
        self.assertNotEqual(by_name["node_identity_types"].status, "skipped")
        self.assertNotEqual(by_name["duplicate_node_coordinates"].status, "skipped")
        self.assertEqual(by_name["single_connected_component"].status, "skipped")
        self.assertIn("topology", by_name["single_connected_component"].message)

    def test_ccf_fraction_runs_despite_duplicate_ids(self):
        # The motivating case: fraction of nodes outside the CCF mesh is a pure
        # coordinate lookup, so it runs even though the ids are duplicated.
        ann = np.zeros((10, 10, 10), dtype=np.uint32)
        ann[2:8, 2:8, 2:8] = 500
        ctx = QCContext(space=Space.CCF_REGISTERED, ccf_resolution=1,
                        resources={"ccf_annotation": ann})
        df = _df([
            (1, 1, 5.0, 5.0, 5.0, 1.0, -1),
            (2, 3, 5.0, 5.0, 6.0, 1.0, 1),
            (2, 3, 0.0, 0.0, 0.0, 1.0, 1),  # duplicate id 2, and outside the mesh
        ])
        report = run_qc(df, ctx, metrics=["nodes_outside_ccf_mesh", "single_root_node"])
        by_name = {r.name: r for r in report.results}
        # coordinate metric ran and produced a real scalar
        self.assertNotEqual(by_name["nodes_outside_ccf_mesh"].status, "skipped")
        self.assertIsNotNone(by_name["nodes_outside_ccf_mesh"].value)
        # topology metric was skipped
        self.assertEqual(by_name["single_root_node"].status, "skipped")

    def test_build_failure_skips_even_coordinate_metrics(self):
        # A BUILD-scope failure (missing column) skips everything -- there are no
        # arrays for even a coordinate-only metric to run on.
        df = _df([(1, 1, 0.0, 0.0, 0.0, -1)],
                 cols=["node_id", "compartment", "x", "y", "z", "parent"])  # no 'r'
        report = run_qc(df, IMG, metrics=["node_identity_types"])
        self.assertFalse(report.morphology_evaluated)
        # required_columns blocks (BUILD), so later integrity checks are not run
        self.assertEqual([r.name for r in report.integrity_results], ["required_columns"])
        self.assertEqual(report.integrity_results[0].status, "fail")
        self.assertEqual(report.results[0].status, "skipped")  # even the coordinate metric
        self.assertIn("built", report.results[0].message)

    def test_buildability_runs_even_for_custom_metric_list(self):
        report = run_qc(_clean_tree(), IMG, metrics=["single_connected_component"])
        self.assertEqual(
            [r.name for r in report.integrity_results],
            ["required_columns", "non_empty", "castable_columns", "unique_node_ids",
             "valid_parent_references", "acyclic"],
        )
        self.assertEqual([r.name for r in report.results], ["single_connected_component"])

    def test_non_numeric_value_reports_instead_of_crashing(self):
        # Previously this raised ValueError inside the build; now it is a report.
        df = _df([
            (1, 1, "abc", 0.0, 0.0, 1.0, -1),
            (2, 3, 1.0, 0.0, 0.0, 1.0, 1),
        ])
        report = run_qc(df, IMG, metrics=["single_connected_component"])
        self.assertFalse(report.morphology_evaluated)
        cast = next(r for r in report.integrity_results if r.name == "castable_columns")
        self.assertEqual(cast.status, "fail")
        self.assertEqual(report.results[0].status, "skipped")

    def test_null_node_id_reports_instead_of_silent_garbage(self):
        # A NaN node_id used to cast to a garbage integer with only a warning.
        df = _df([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
            (np.nan, 3, 1.0, 0.0, 0.0, 1.0, 1),
        ])
        report = run_qc(df, IMG, metrics=["single_connected_component"])
        self.assertFalse(report.morphology_evaluated)
        cast = next(r for r in report.integrity_results if r.name == "castable_columns")
        self.assertEqual(cast.status, "fail")

    def test_read_swc_does_not_crash_on_bad_values(self):
        import os
        import tempfile
        from standard_morph.swc_io import read_swc

        content = "1 1 abc 0 0 1 -1\n2 3 1 0 0 1 1\n"
        with tempfile.NamedTemporaryFile("w", suffix=".swc", delete=False) as f:
            f.write(content)
            path = f.name
        try:
            df = read_swc(path)  # must not raise
            self.assertTrue(df["x"].isna().any())  # 'abc' -> NaN
        finally:
            os.remove(path)

    def test_report_with_skips_is_json_serialisable(self):
        import json

        df = _df([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
            (2, 3, 1.0, 0.0, 0.0, 1.0, 1),
            (2, 3, 9.0, 9.0, 9.0, 1.0, 1),
        ])
        report = run_qc(df, IMG, suite_name="default_pre_registration_tests")
        d = json.dumps(report.to_dict())
        self.assertIn("morphology_evaluated", d)
        self.assertIn("integrity_ok", d)


if __name__ == "__main__":
    unittest.main()
