"""Metric result object.

Every metric returns a `MetricResult`. Quantitative metrics always populate
`measurements`; failing metrics return the *full* flagged node lists (no
sampling or capping) per the architecture spec.

``status`` is one of:

* ``"pass"``    -- the check ran and the morphology is fine on this axis.
* ``"fail"``    -- the check ran and found an objective defect.
* ``"review"``  -- the check ran and flagged something unusual but possibly
  valid, which a human must judge (e.g. multiple apical origins). Not a failure:
  it makes the run's ``overall_status`` ``"incomplete"``, not ``"fail"``. A
  metric declares this via ``Metric.violation_severity`` (see
  ``standard_morph.metrics.base.Severity``).
* ``"error"``   -- the check could not complete (e.g. a precondition it needs
  was absent).
* ``"skipped"`` -- the check did not run. Morphology-quality metrics are skipped
  when a blocking input-integrity metric failed, so the morphology could not be
  built (see ``standard_morph.metrics.base``).
"""
from dataclasses import dataclass, field
from typing import List, Optional, Tuple


@dataclass
class MetricResult:
    """The output of a single metric's ``evaluate()`` call.

    Every metric returns exactly one ``MetricResult``. The ``status`` field is
    the primary verdict; the remaining fields carry evidence for dashboards,
    debugging, and downstream review.

    ``flagged_node_ids`` and ``flagged_node_coordinates`` are always
    parallel — element *i* of the coordinates list is the position of node
    *i* in the id list. Metrics must maintain this invariant.

    ``value`` is a single canonical scalar suitable for trending across many
    cells (e.g. a maximum, count, or fraction). Leave it ``None`` for metrics
    that are genuinely binary and have no meaningful scalar summary.
    """

    name: str
    status: str = "pass"  # "pass" | "fail" | "review" | "error" | "skipped"
    message: str = ""
    # Canonical headline scalar for the metric (e.g. a max, count, or fraction),
    # for dashboards/trends. None for metrics that are genuinely binary.
    value: Optional[float] = None
    value_label: Optional[str] = None  # what `value` measures, e.g. "max_tortuosity"
    thresholds_used: dict = field(default_factory=dict)
    measurements: dict = field(default_factory=dict)
    flagged_node_ids: List[int] = field(default_factory=list)
    flagged_node_coordinates: List[Tuple[float, float, float]] = field(default_factory=list)
    counts: dict = field(default_factory=dict)
    #: Generated QC artifacts (e.g. images), each a dict {type, path, description}.
    artifacts: List[dict] = field(default_factory=list)
    runtime_ms: float = 0.0
