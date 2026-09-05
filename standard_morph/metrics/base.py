"""Metric base class, evaluation phases, and applicability declaration.

QC runs in **two phases**, reflecting two genuinely different kinds of "wrong":

1. ``INPUT_INTEGRITY`` -- is the raw SWC *table* even well-formed? (unique node
   ids, required columns, non-empty, parseable). These metrics read the raw
   ``pandas.DataFrame`` before any graph is built, because a malformed table
   cannot be turned into a faithful morphology in the first place.
2. ``MORPHOLOGY_QUALITY`` -- given a well-formed table, does the *neuron* have
   problems? (roots, connectivity, tortuosity, position in the CCF, ...). These
   metrics read a fully-built ``PreparedMorphology``.

The engine runs the integrity phase first. If an integrity metric that is
marked ``blocks_morphology_phase`` fails, the morphology phase is not run (you
cannot do graph checks on a graph that could not be built); those metrics are
reported with status ``"skipped"`` instead of crashing. This keeps a batch of
thousands of files robust: a malformed file yields a normal report saying
exactly what is wrong, never a stack trace.

A metric therefore declares which phase it belongs to (input-integrity vs morphology). The default is
``MORPHOLOGY_QUALITY`` so that ordinary metrics need no extra boilerplate.
"""
from dataclasses import dataclass, field
from enum import Enum

from standard_morph.exceptions import IncompatibleMetricContextError, MissingPolicyValueError


class EvaluationPhase(str, Enum):
    """Which phase of a QC run a metric belongs to, and hence what data it reads.

    * ``INPUT_INTEGRITY`` -- ``evaluate`` receives the raw ``pandas.DataFrame``
      (an unvetted SWC table, not yet a graph).
    * ``MORPHOLOGY_QUALITY`` -- ``evaluate`` receives a ``PreparedMorphology``
      (a built, validated graph).
    """

    INPUT_INTEGRITY = "input_integrity"
    MORPHOLOGY_QUALITY = "morphology_quality"


class BlockScope(str, Enum):
    """What an input-integrity metric's *failure* invalidates. 
    For example, if a swc file is missing a node ID, some morphology metrics, specifically 
    topology based ones cannot run because the graph cannot be built. The engine will skip 
    those metrics and report them as "skipped" instead of crashing.

    Set on ``Metric.blocks_on_failure`` (only meaningful for ``INPUT_INTEGRITY``
    metrics). ``None`` means a failure blocks nothing -- it is informational.

    * ``BUILD``    -- the morphology cannot be constructed at all (e.g. a
      required column is missing, or the table is empty). *Every* morphology
      metric is skipped, because there are not even coordinates to run on.
    * ``TOPOLOGY`` -- the arrays build fine, but the tree topology
      (``parent`` / ``children`` / ``roots`` / ``segments``) is untrustworthy
      (e.g. duplicate node ids make the id->index map ambiguous). Only metrics
      that *need* topology are skipped; coordinate/attribute-only metrics still
      run, since their per-node data (``xyz``, ``compartment``) is intact.
    """

    BUILD = "build"
    TOPOLOGY = "topology"


class Severity(str, Enum):
    """What a metric's *violation* means for the report.

    Set on ``Metric.violation_severity``; it decides the ``status`` a metric
    reports when its check is violated. Both severities still *run* and report
    normally -- the difference is only how the result is graded and rolled up.

    * ``FAIL``   -- an objective defect: the file is malformed or the tree is
      impossible. Reports status ``"fail"`` and makes the run's
      ``overall_status`` ``"fail"``.
    * ``REVIEW`` -- unusual but possibly valid: a threshold or convention that a
      human must judge (e.g. multiple apical origins, high branch degree).
      Reports status ``"review"``; the run's ``overall_status`` becomes
      ``"incomplete"`` (needs human oversight) rather than ``"fail"``.

    A metric's bucket is declared once, at the top of its class, so moving a
    check between ``FAIL`` and ``REVIEW`` is a one-line change with no engine or
    framework edits.
    """

    FAIL = "fail"
    REVIEW = "review"


@dataclass(frozen=True)
class Applicability:
    """Declares the contexts in which a metric may run.

    Parameters
    ----------
    spaces : frozenset[Space]
        Coordinate spaces the metric supports.
    morphology_kinds : frozenset[MorphologyKind]
        Morphology kinds the metric supports.
    required_resources : frozenset[str]
        Keys that must be present in ``context.resources`` for the metric to run
        (e.g. ``"image_path"``, ``"ccf_atlas_path"``).
    """

    spaces: frozenset
    morphology_kinds: frozenset
    required_resources: frozenset = field(default_factory=frozenset)


class Metric:
    """Base class for QC metrics.

    Subclasses set the ``name``, ``display_name``, and ``applicability`` class
    attributes and implement :meth:`evaluate`. They may also override the two
    phase-related attributes below; the defaults make a metric an ordinary
    morphology-quality check.

    Attributes
    ----------
    evaluation_phase : EvaluationPhase
        Which phase this metric runs in. Defaults to ``MORPHOLOGY_QUALITY``.
        Determines what :meth:`evaluate` is handed (a DataFrame for
        ``INPUT_INTEGRITY``, a ``PreparedMorphology`` for ``MORPHOLOGY_QUALITY``).
    blocks_on_failure : BlockScope or None
        Only meaningful for ``INPUT_INTEGRITY`` metrics: what this metric's
        *failure* invalidates -- ``BUILD`` (can't build anything), ``TOPOLOGY``
        (topology untrustworthy), or ``None`` (failure blocks nothing).
    requires_topology : bool
        Only meaningful for ``MORPHOLOGY_QUALITY`` metrics: whether the check
        needs valid tree topology. Coordinate/attribute-only metrics (that read
        just ``xyz`` / ``compartment``) set this ``False`` so they still run when
        only the topology is unreliable. Defaults to ``True``.
    violation_severity : Severity
        What this metric's *violation* means -- ``FAIL`` (an objective defect) or
        ``REVIEW`` (unusual, needs human oversight). Defaults to ``FAIL``.
        Metrics report ``self.violation_severity.value`` as their status when
        their check is violated, so switching a metric between the two buckets is
        a single-line change here. See :class:`Severity`.
    """

    name = None
    display_name = None
    applicability = None  # type: Applicability
    required_policy_keys = frozenset()

    evaluation_phase = EvaluationPhase.MORPHOLOGY_QUALITY
    blocks_on_failure = None  # BlockScope | None; see class docstring
    requires_topology = True
    violation_severity = Severity.FAIL  # Severity.FAIL | Severity.REVIEW

    def validate_context(self, context):
        """Raise if this metric is incompatible with ``context``.

        Incompatibility is terminal per the architecture spec: the caller is
        expected to let this propagate and halt the run.
        """
        reasons = []
        if context.space not in self.applicability.spaces:
            allowed = sorted(s.value for s in self.applicability.spaces)
            reasons.append(f"space '{context.space.value}' not in allowed {allowed}")
        if context.morphology_kind not in self.applicability.morphology_kinds:
            allowed = sorted(m.value for m in self.applicability.morphology_kinds)
            reasons.append(
                f"morphology_kind '{context.morphology_kind.value}' not in allowed {allowed}"
            )
        missing = [r for r in self.applicability.required_resources if r not in context.resources]
        if missing:
            reasons.append(f"missing required resources: {sorted(missing)}")

        if reasons:
            raise IncompatibleMetricContextError(
                f"Metric '{self.name}' is incompatible with context: " + "; ".join(reasons)
            )

    def is_applicable(self, context):
        """Return True if the metric can run in ``context`` (no raise)."""
        try:
            self.validate_context(context)
            return True
        except IncompatibleMetricContextError:
            return False

    def evaluate(self, data_under_test, context, policy):
        """Compute the metric and return a `MetricResult`.

        Parameters
        ----------
        data_under_test :
            The data being interrogated. Its type depends on
            ``evaluation_phase``: a ``pandas.DataFrame`` (raw SWC table) for
            ``INPUT_INTEGRITY`` metrics, a ``PreparedMorphology`` for
            ``MORPHOLOGY_QUALITY`` metrics. Concrete metrics name this parameter
            more specifically -- ``swc_df`` or ``prepared_morph`` respectively.
        context : QCContext
            The run context (space, morphology kind, resources, policy version).
        policy : Policy
            The active threshold policy.
        """
        raise NotImplementedError


def threshold_for_space(policy, name, key, space_value):
    """Return a threshold that may vary by coordinate space.

    Use this instead of ``policy[name, key]`` for thresholds whose optimal
    value differs between coordinate spaces (e.g. edge lengths before vs.
    after resampling).

    The policy value may be either:

    * A **scalar** -- the same threshold applies in every space.
    * A **per-space dict** ``{space_value: threshold, ...}`` -- the threshold
      is looked up by ``space_value`` (a ``Space.value`` string such as
      ``"image_space"`` or ``"ccf_registered"``).

    The resolved value (after any space lookup) may itself be a plain scalar
    **or** a :class:`~standard_morph.models.qc_policy.PolicyRange`. Callers
    should branch on ``isinstance(resolved, PolicyRange)`` to choose between
    ``resolved.contains(value)`` and a direct comparison.

    Raises
    ------
    MissingPolicyValueError
        If the top-level key is absent from the policy, or if the value is a
        dict that does not contain an entry for ``space_value``.
    """

    val = policy[name, key]
    if isinstance(val, dict):
        if space_value not in val:
            raise MissingPolicyValueError(
                f"Policy has no value for {name}.{key}[{space_value}]"
            )
        return val[space_value]
    return val


