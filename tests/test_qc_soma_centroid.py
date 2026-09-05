import os
import tempfile
import unittest

import pandas as pd

from standard_morph.engine import run_qc
from standard_morph.exceptions import IncompatibleMetricContextError, MissingPolicyValueError
from standard_morph.metrics.soma_at_centroid import SomaAtCentroidMetric
from standard_morph.preparation import PreparedMorphology
from standard_morph.models.qc_context import QCContext, Space
from standard_morph.models.qc_policy import Policy

POLICY = Policy("t", {"soma_at_centroid": {"max_offset_fraction": 0.5}})
EMPTY = Policy("t", {})  # missing max_offset_fraction -> MissingPolicyValueError


def _df(rows):
    return pd.DataFrame(rows, columns=["node_id", "compartment", "x", "y", "z", "r", "parent"])


def _pm(rows):
    return PreparedMorphology.from_dataframe(_df(rows))


def _ctx(resources):
    return QCContext(space=Space.IMAGE_SPACE, resources=resources)


def _soma_at(x, y, z):
    # A one-soma tree whose soma sits at (x, y, z).
    return _pm([
        (1, 1, x, y, z, 1.0, -1),
        (2, 3, x, y, z + 10.0, 1.0, 1),
    ])


REF = {"image_soma_xyz": (0.0, 0.0, 0.0), "image_soma_radius_xyz": (100.0, 50.0, 200.0)}


class TestSomaAtCentroid(unittest.TestCase):
    def test_the_worked_example_fails(self):
        # ref (0,0,0), radius (100,50,200), swc (99,40,100)
        # -> per-axis (0.99, 0.80, 0.50); worst axis 0.99 > 0.5 -> fail.
        pm = _soma_at(99.0, 40.0, 100.0)
        r = SomaAtCentroidMetric().evaluate(pm, _ctx(dict(REF)), POLICY)
        self.assertEqual(r.status, "fail")
        self.assertAlmostEqual(r.value, 0.99)
        self.assertEqual(r.flagged_node_ids, [1])
        frac = r.measurements["per_axis_offset_fraction"]
        self.assertAlmostEqual(frac["x"], 0.99)
        self.assertAlmostEqual(frac["y"], 0.80)
        self.assertAlmostEqual(frac["z"], 0.50)
        # x and y exceed 0.5; z is exactly 0.5 (not over).
        self.assertIn("x", r.message)
        self.assertIn("y", r.message)

    def test_near_centroid_passes(self):
        # swc (40, 20, 50) -> (0.40, 0.40, 0.25); worst 0.40 <= 0.5 -> pass.
        pm = _soma_at(40.0, 20.0, 50.0)
        r = SomaAtCentroidMetric().evaluate(pm, _ctx(dict(REF)), POLICY)
        self.assertEqual(r.status, "pass")
        self.assertAlmostEqual(r.value, 0.40)

    def test_missing_threshold_raises(self):
        pm = _soma_at(60.0, 0.0, 0.0)
        with self.assertRaises(MissingPolicyValueError):
            SomaAtCentroidMetric().evaluate(pm, _ctx(dict(REF)), EMPTY)

    def test_multiple_somas_error(self):
        pm = _pm([
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
            (2, 3, 0.0, 0.0, 10.0, 1.0, 1),
            (3, 1, 500.0, 500.0, 500.0, 1.0, -1),  # a second soma root
        ])
        r = SomaAtCentroidMetric().evaluate(pm, _ctx(dict(REF)), POLICY)
        self.assertEqual(r.status, "error")
        self.assertEqual(r.measurements["n_soma"], 2)

    def test_nonpositive_radius_error(self):
        res = {"image_soma_xyz": (0.0, 0.0, 0.0), "image_soma_radius_xyz": (100.0, 0.0, 200.0)}
        r = SomaAtCentroidMetric().evaluate(_soma_at(1.0, 1.0, 1.0), _ctx(res), POLICY)
        self.assertEqual(r.status, "error")
        self.assertIn("positive", r.message)

    def test_bad_shape_error(self):
        res = {"image_soma_xyz": (0.0, 0.0), "image_soma_radius_xyz": (100.0, 50.0, 200.0)}
        r = SomaAtCentroidMetric().evaluate(_soma_at(1.0, 1.0, 1.0), _ctx(res), POLICY)
        self.assertEqual(r.status, "error")


class TestSomaCentroidPlumbing(unittest.TestCase):
    def test_missing_resources_is_incompatible(self):
        with self.assertRaises(IncompatibleMetricContextError):
            run_qc(_df([(1, 1, 0.0, 0.0, 0.0, 1.0, -1)]),
                   QCContext(space=Space.IMAGE_SPACE), metrics=["soma_at_centroid"])

    def test_wrong_space_is_incompatible(self):
        # Image-based reference -> not applicable in CCF space.
        with self.assertRaises(IncompatibleMetricContextError):
            run_qc(
                _df([(1, 1, 0.0, 0.0, 0.0, 1.0, -1), (2, 3, 0.0, 0.0, 10.0, 1.0, 1)]),
                QCContext(space=Space.CCF_REGISTERED, resources=dict(REF)),
                metrics=["soma_at_centroid"],
            )

    def test_runs_through_engine(self):
        df = _df([(1, 1, 40.0, 20.0, 50.0, 1.0, -1), (2, 3, 40.0, 20.0, 60.0, 1.0, 1)])
        report = run_qc(df, _ctx(dict(REF)), metrics=["soma_at_centroid"])
        r = next(x for x in report.results if x.name == "soma_at_centroid")
        self.assertEqual(r.status, "pass")


class TestSomaMipDegradation(unittest.TestCase):
    def test_bad_zarr_path_records_error_but_still_evaluates(self):
        # A bogus zarr path (or missing optional deps) must not fail the metric --
        # the offset is still computed and the failure is recorded.
        with tempfile.TemporaryDirectory() as d:
            res = dict(REF)
            res["image_zarr_path"] = os.path.join(d, "does_not_exist.zarr")
            res["soma_mip_path"] = os.path.join(d, "out.png")
            r = SomaAtCentroidMetric().evaluate(_soma_at(40.0, 20.0, 50.0), _ctx(res), POLICY)
        self.assertEqual(r.status, "pass")
        self.assertIn("soma_mip_error", r.measurements)
        self.assertEqual(r.artifacts, [])

    def test_zarr_path_without_output_path_records_error(self):
        res = dict(REF)
        res["image_zarr_path"] = "s3://bucket/img.zarr"
        r = SomaAtCentroidMetric().evaluate(_soma_at(40.0, 20.0, 50.0), _ctx(res), POLICY)
        self.assertEqual(r.status, "pass")
        self.assertIn("soma_mip_path missing", r.measurements["soma_mip_error"])


if __name__ == "__main__":
    unittest.main()
