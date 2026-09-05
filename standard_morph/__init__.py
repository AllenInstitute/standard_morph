"""standard_morph: modular quality control for neuron SWC reconstructions.

Typical use::

    from standard_morph import run_qc, QCContext, Space

    context = QCContext(space=Space.IMAGE_SPACE)
    report = run_qc("cell.swc", context, suite_name="default_pre_registration_tests")
    print(report.summary)

See docs/qc_overhaul_architecture_spec.md for the design.
"""
from standard_morph.preparation import PreparedMorphology
from standard_morph.swc_io import read_swc
from standard_morph.atlas import load_ccf_annotation, clear_atlas_cache
from standard_morph.models.qc_context import (
    QCContext,
    Space,
    MorphologyKind,
    ALL_MORPHOLOGY_KINDS,
)
from standard_morph.metrics.base import (
    EvaluationPhase,
    BlockScope,
    Severity,
    Metric,
    Applicability,
    threshold_for_space,
)
from standard_morph.models.qc_result import MetricResult
from standard_morph.models.qc_policy import Policy, PolicyRange
from standard_morph.models.qc_run import RunReport, SCHEMA_VERSION
from standard_morph.policies import get_policy, available_policies
from standard_morph.suites import resolve_suite, available_suites, BUILTIN_SUITES
from standard_morph.registry import REGISTRY, register
from standard_morph.engine import run_qc
from standard_morph.exceptions import (
    QCError,
    IncompatibleMetricContextError,
    MissingPolicyValueError,
    MissingPolicyValuesError,
)

__all__ = [
    "run_qc",
    "PreparedMorphology",
    "read_swc",
    "load_ccf_annotation",
    "clear_atlas_cache",
    "QCContext",
    "Space",
    "MorphologyKind",
    "ALL_MORPHOLOGY_KINDS",
    "EvaluationPhase",
    "BlockScope",
    "Severity",
    "Metric",
    "Applicability",
    "threshold_for_space",
    "register",
    "MissingPolicyValueError",
    "MissingPolicyValuesError",
    "MetricResult",
    "Policy",
    "PolicyRange",
    "RunReport",
    "SCHEMA_VERSION",
    "get_policy",
    "available_policies",
    "resolve_suite",
    "available_suites",
    "BUILTIN_SUITES",
    "REGISTRY",
    "QCError",
    "IncompatibleMetricContextError",
]
