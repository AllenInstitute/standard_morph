import os
import tempfile
import unittest

import pandas as pd

from standard_morph.engine import run_qc
from standard_morph.exceptions import IncompatibleMetricContextError, MissingPolicyValueError
from standard_morph.metrics.filename_format import FilenameFormatMetric, _is_valid_filename
from standard_morph.models.qc_context import QCContext, Space
from standard_morph.models.qc_policy import Policy
from standard_morph.preparation import PreparedMorphology

AIND = Policy("t", {"filename_format": {"name_format": "AIND"}})
AIBS = Policy("t", {"filename_format": {"name_format": "AIBS"}})
EMPTY = Policy("t", {})  # missing name_format -> MissingPolicyValueError


def _valid_df():
    return pd.DataFrame(
        [
            (1, 1, 0.0, 0.0, 0.0, 1.0, -1),
            (2, 3, 0.0, 0.0, 10.0, 1.0, 1),
        ],
        columns=["node_id", "compartment", "x", "y", "z", "r", "parent"],
    )


def _write_swc(directory, name):
    path = os.path.join(directory, name)
    with open(path, "w") as f:
        f.write("1 1 0.0 0.0 0.0 1.0 -1\n2 3 0.0 0.0 10.0 1.0 1\n")
    return path


class TestFilenamePattern(unittest.TestCase):
    def test_valid_aind_names(self):
        for name in (
            "N123-000000.swc",
            "N1_123456_ABC.swc",
            "N42-654321-consensus.swc",
            "N7-123456-axon-XY.swc",
            "N7_123456_dendrite_consensus.swc",
        ):
            self.assertTrue(_is_valid_filename(name, "AIND"), name)

    def test_invalid_aind_names(self):
        for name in (
            "neuron.swc",             # no N-prefixed id
            "123-000000.swc",         # missing leading N
            "N123-00000.swc",         # 5-digit block, needs 6
            "N123-000000-TOOLONG.swc",  # tag is neither 2-3 letters nor consensus
            "N123-000000.txt",        # wrong extension
        ):
            self.assertFalse(_is_valid_filename(name, "AIND"), name)

    def test_aibs_is_todo_and_always_passes(self):
        self.assertTrue(_is_valid_filename("anything at all", "AIBS"))
        self.assertTrue(_is_valid_filename("N123-000000.swc", "AIBS"))

    def test_unknown_format_raises(self):
        with self.assertRaises(ValueError):
            _is_valid_filename("N123-000000.swc", "MADEUP")


class TestFilenameMetric(unittest.TestCase):
    def _run(self, filename, policy):
        ctx = QCContext(space=Space.IMAGE_SPACE, resources={"filename": filename})
        return FilenameFormatMetric().evaluate(None, ctx, policy)  # swc_df is ignored

    def test_valid_passes(self):
        r = self._run("N123-000000.swc", AIND)
        self.assertEqual(r.status, "pass")
        self.assertTrue(r.value)

    def test_invalid_fails(self):
        r = self._run("garbage.swc", AIND)
        self.assertEqual(r.status, "fail")
        self.assertFalse(r.value)

    def test_full_path_is_reduced_to_basename(self):
        r = self._run("/some/nested/dir/N123-000000.swc", AIND)
        self.assertEqual(r.status, "pass")
        self.assertEqual(r.measurements["filename"], "N123-000000.swc")

    def test_aibs_passes_regardless(self):
        r = self._run("not a valid name", AIBS)
        self.assertEqual(r.status, "pass")

    def test_missing_name_format_raises(self):
        with self.assertRaises(MissingPolicyValueError):
            self._run("garbage.swc", EMPTY)


class TestFilenamePlumbing(unittest.TestCase):
    def test_missing_filename_is_incompatible(self):
        # DataFrame input with no resources["filename"] -> fail-fast, like a
        # CCF metric with no ccf_resolution.
        with self.assertRaises(IncompatibleMetricContextError):
            run_qc(_valid_df(), QCContext(space=Space.IMAGE_SPACE), metrics=["filename_format"])

    def test_missing_filename_is_incompatible_for_prepared_morph(self):
        pm = PreparedMorphology.from_dataframe(_valid_df())
        with self.assertRaises(IncompatibleMetricContextError):
            run_qc(pm, QCContext(space=Space.IMAGE_SPACE), metrics=["filename_format"])

    def test_prepared_morph_skips_all_integrity_metrics(self):
        # PreparedMorphology has no raw table; all integrity metrics are skipped.
        pm = PreparedMorphology.from_dataframe(_valid_df())
        ctx = QCContext(space=Space.IMAGE_SPACE, resources={"filename": "N123-000000.swc"})
        report = run_qc(pm, ctx, metrics=["filename_format"])
        fr = next(r for r in report.integrity_results if r.name == "filename_format")
        self.assertEqual(fr.status, "skipped")

    def test_prepared_morph_integrity_skip_message_is_clear(self):
        pm = PreparedMorphology.from_dataframe(_valid_df())
        ctx = QCContext(space=Space.IMAGE_SPACE, resources={"filename": "garbage.swc"})
        report = run_qc(pm, ctx, metrics=["filename_format"])
        fr = next(r for r in report.integrity_results if r.name == "filename_format")
        self.assertEqual(fr.status, "skipped")
        self.assertIn("PreparedMorphology", fr.message)

    def test_filename_from_resources(self):
        ctx = QCContext(space=Space.IMAGE_SPACE, resources={"filename": "N123-000000.swc"})
        report = run_qc(_valid_df(), ctx, metrics=["filename_format"])
        fr = next(r for r in report.integrity_results if r.name == "filename_format")
        self.assertEqual(fr.status, "pass")
        self.assertTrue(report.integrity_ok)

    def test_filename_auto_derived_from_path(self):
        with tempfile.TemporaryDirectory() as d:
            path = _write_swc(d, "N123-000000.swc")
            report = run_qc(path, QCContext(space=Space.IMAGE_SPACE), metrics=["filename_format"])
        fr = next(r for r in report.integrity_results if r.name == "filename_format")
        self.assertEqual(fr.status, "pass")
        self.assertEqual(fr.measurements["filename"], "N123-000000.swc")

    def test_bad_name_auto_derived_fails_the_report(self):
        with tempfile.TemporaryDirectory() as d:
            path = _write_swc(d, "not_a_valid_name.swc")
            report = run_qc(path, QCContext(space=Space.IMAGE_SPACE), metrics=["filename_format"])
        fr = next(r for r in report.integrity_results if r.name == "filename_format")
        self.assertEqual(fr.status, "fail")
        # A bad name does not block the build, but it does fail the overall run.
        self.assertFalse(report.integrity_ok)
        self.assertFalse(report.passed)

    def test_explicit_filename_overrides_path(self):
        # An explicit resources["filename"] wins over the path's basename.
        with tempfile.TemporaryDirectory() as d:
            path = _write_swc(d, "not_a_valid_name.swc")
            ctx = QCContext(space=Space.IMAGE_SPACE, resources={"filename": "N123-000000.swc"})
            report = run_qc(path, ctx, metrics=["filename_format"])
        fr = next(r for r in report.integrity_results if r.name == "filename_format")
        self.assertEqual(fr.status, "pass")


if __name__ == "__main__":
    unittest.main()
