"""Tests for PolicyRange — the inclusive [lo, hi] policy threshold type."""
import pytest

from standard_morph.models.qc_policy import Policy, PolicyRange
from standard_morph.exceptions import MissingPolicyValueError
from standard_morph.metrics.base import threshold_for_space


# ---------------------------------------------------------------------------
# PolicyRange construction and validation
# ---------------------------------------------------------------------------

class TestPolicyRangeConstruction:
    def test_valid_range(self):
        r = PolicyRange(lo=5.0, hi=50.0)
        assert r.lo == 5.0
        assert r.hi == 50.0

    def test_equal_bounds_allowed(self):
        r = PolicyRange(lo=10.0, hi=10.0)
        assert r.lo == r.hi

    def test_inverted_bounds_raise(self):
        with pytest.raises(ValueError, match="must not exceed"):
            PolicyRange(lo=50.0, hi=5.0)

    def test_frozen(self):
        r = PolicyRange(lo=1.0, hi=2.0)
        with pytest.raises(Exception):  # FrozenInstanceError
            r.lo = 99.0


# ---------------------------------------------------------------------------
# PolicyRange.contains
# ---------------------------------------------------------------------------

class TestPolicyRangeContains:
    def setup_method(self):
        self.r = PolicyRange(lo=10.0, hi=50.0)

    def test_value_below_lo_fails(self):
        assert not self.r.contains(9.99)

    def test_value_at_lo_passes(self):
        assert self.r.contains(10.0)

    def test_value_in_middle_passes(self):
        assert self.r.contains(30.0)

    def test_value_at_hi_passes(self):
        assert self.r.contains(50.0)

    def test_value_above_hi_fails(self):
        assert not self.r.contains(50.01)

    def test_negative_range(self):
        r = PolicyRange(lo=-20.0, hi=-5.0)
        assert r.contains(-10.0)
        assert not r.contains(0.0)
        assert not r.contains(-21.0)

    def test_zero_width_range(self):
        r = PolicyRange(lo=7.0, hi=7.0)
        assert r.contains(7.0)
        assert not r.contains(7.001)
        assert not r.contains(6.999)


# ---------------------------------------------------------------------------
# PolicyRange as a policy value — Policy.__getitem__
# ---------------------------------------------------------------------------

class TestPolicyRangeInPolicy:
    def setup_method(self):
        self.policy = Policy(
            version="test",
            thresholds={
                "my_metric": {
                    "range_key": PolicyRange(lo=5.0, hi=50.0),
                    "scalar_key": 30.0,
                }
            },
        )

    def test_getitem_returns_range(self):
        val = self.policy["my_metric", "range_key"]
        assert isinstance(val, PolicyRange)
        assert val.lo == 5.0
        assert val.hi == 50.0

    def test_getitem_returns_scalar(self):
        val = self.policy["my_metric", "scalar_key"]
        assert val == 30.0

    def test_missing_key_raises(self):
        with pytest.raises(MissingPolicyValueError):
            _ = self.policy["my_metric", "nonexistent"]


# ---------------------------------------------------------------------------
# PolicyRange as a per-space value inside threshold_for_space
# ---------------------------------------------------------------------------

class TestThresholdForSpaceWithRange:
    def setup_method(self):
        self.policy = Policy(
            version="test",
            thresholds={
                "my_metric": {
                    # per-space dict whose values are themselves PolicyRanges
                    "range_per_space": {
                        "image_space": PolicyRange(lo=5.0, hi=30.0),
                        "ccf_registered": PolicyRange(lo=1.0, hi=10.0),
                    },
                    # plain scalar per-space (existing behavior, must still work)
                    "scalar_per_space": {
                        "image_space": 30.0,
                        "ccf_registered": 10.0,
                    },
                    # top-level PolicyRange (same for all spaces)
                    "range_flat": PolicyRange(lo=0.0, hi=100.0),
                }
            },
        )

    def test_per_space_range_image(self):
        val = threshold_for_space(self.policy, "my_metric", "range_per_space", "image_space")
        assert isinstance(val, PolicyRange)
        assert val.lo == 5.0
        assert val.hi == 30.0

    def test_per_space_range_ccf(self):
        val = threshold_for_space(self.policy, "my_metric", "range_per_space", "ccf_registered")
        assert isinstance(val, PolicyRange)
        assert val.lo == 1.0
        assert val.hi == 10.0

    def test_per_space_scalar_still_works(self):
        val = threshold_for_space(self.policy, "my_metric", "scalar_per_space", "image_space")
        assert val == 30.0

    def test_flat_range_returned_regardless_of_space(self):
        val = threshold_for_space(self.policy, "my_metric", "range_flat", "image_space")
        assert isinstance(val, PolicyRange)
        assert val.contains(50.0)

    def test_missing_space_raises(self):
        with pytest.raises(MissingPolicyValueError):
            threshold_for_space(self.policy, "my_metric", "range_per_space", "unknown_space")


# ---------------------------------------------------------------------------
# Typical metric-side usage pattern
# ---------------------------------------------------------------------------

class TestTypicalMetricPattern:
    """Demonstrate the idiomatic isinstance branch a metric would use."""

    @staticmethod
    def _check(threshold, computed) -> bool:
        if isinstance(threshold, PolicyRange):
            return threshold.contains(computed)
        return computed <= threshold

    def test_scalar_below_passes(self):
        assert self._check(50.0, 30.0)

    def test_scalar_above_fails(self):
        assert not self._check(50.0, 51.0)

    def test_range_inside_passes(self):
        assert self._check(PolicyRange(lo=10.0, hi=50.0), 25.0)

    def test_range_outside_fails(self):
        assert not self._check(PolicyRange(lo=10.0, hi=50.0), 5.0)

    def test_range_at_bounds_passes(self):
        assert self._check(PolicyRange(lo=10.0, hi=50.0), 10.0)
        assert self._check(PolicyRange(lo=10.0, hi=50.0), 50.0)


# ---------------------------------------------------------------------------
# Public API export
# ---------------------------------------------------------------------------

def test_policy_range_importable_from_top_level():
    from standard_morph import PolicyRange as PR  # noqa: F401
    assert PR is PolicyRange


def test_policy_range_importable_from_policies():
    from standard_morph.policies import PolicyRange as PR  # noqa: F401
    assert PR is PolicyRange
