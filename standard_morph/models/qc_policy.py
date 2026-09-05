"""Versioned threshold policy object.

A policy is a versioned, read-only mapping of ``metric_name -> {key: value}``.
Metrics read thresholds from the policy rather than receiving them as call-time
arguments, and the policy version is recorded in every run report so results
are traceable.
"""
from dataclasses import dataclass, field

from standard_morph.exceptions import MissingPolicyValueError


@dataclass(frozen=True)
class PolicyRange:
    """An inclusive numeric range ``[lo, hi]`` for a policy threshold.

    Use instead of a plain scalar when a metric should accept any value within
    a band rather than everything below (or above) a single cutoff.

    A metric checks membership with ``threshold.contains(computed_value)``,
    which returns ``True`` when ``lo <= value <= hi``.

    Parameters
    ----------
    lo : float
        Lower bound (inclusive).
    hi : float
        Upper bound (inclusive).

    Examples
    --------
    In a policy definition::

        "some_metric": {"some_key": PolicyRange(lo=5.0, hi=50.0)}

    In a metric's ``evaluate()``::

        threshold = policy[self.name, "some_key"]
        if isinstance(threshold, PolicyRange):
            passed = threshold.contains(computed_value)
        else:
            passed = computed_value <= threshold   # plain upper-bound scalar
    """

    lo: float
    hi: float

    def __post_init__(self):
        if self.lo > self.hi:
            raise ValueError(
                f"PolicyRange lower bound {self.lo} must not exceed upper bound {self.hi}"
            )

    def contains(self, value: float) -> bool:
        """Return True when ``lo <= value <= hi``."""
        return self.lo <= value <= self.hi


@dataclass(frozen=True)
class Policy:
    """Versioned, read-only threshold map for a QC run.

    A ``Policy`` is the single source of truth for all numeric thresholds and
    configuration values used during a run. Metrics read from it rather than
    accepting thresholds as call-time arguments, so every result is traceable
    to an exact, named policy version.

    Parameters
    ----------
    version : str
        The policy version name (e.g. ``"policy_v1"``). Recorded in every
        :class:`~standard_morph.models.qc_run.RunReport` so results can be
        reproduced later.
    thresholds : dict
        Nested mapping of ``{metric_name: {key: value}}``. Values may be plain
        scalars, :class:`PolicyRange` instances, strings, or lists — whatever
        the metric declares via ``required_policy_keys``.
    """

    version: str
    thresholds: dict = field(default_factory=dict)

    def __getitem__(self, item):
        """Return ``policy[metric_name, key]``, raising if absent."""
        metric_name, key = item
        try:
            return self.thresholds[metric_name][key]
        except KeyError:
            raise MissingPolicyValueError(
                f"Policy '{self.version}' has no value for {metric_name}.{key}"
            ) from None

    def for_metric(self, metric_name):
        """Return the full threshold dict for a metric (empty if none)."""
        return dict(self.thresholds.get(metric_name, {}))
