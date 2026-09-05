"""Metric 14 -- Multiple apical dendrite origination points.

An apical dendrite normally leaves the soma as a single trunk. Two or more apical
origination points -- apical nodes (type 4) whose parent is *not* itself apical --
are flagged for interrogation. This can be biologically plausible for some cell
types, but it is unusual and often a tracing error.

Closely related to ``axon_origination`` (metric 13): same "origination = a node
whose parent is a different compartment" idea. The differences: apical dendrite
is optional (zero is fine), and the concern is *more than one* trunk rather than
exactly one. Whether each apical origin stems from a valid parent (soma/apical)
is ``compartment_transitions``' (metric 15) job; this metric only counts trunks.

Note: AIND tracings do not currently distinguish basal (3) from apical (4)
dendrite, so this only fires on files that carry apical (type 4) annotations.
"""
import time

import numpy as np

from standard_morph.metrics.base import Metric, Applicability, Severity
from standard_morph.models.qc_context import ALL_COORDINATE_SPACES, ALL_MORPHOLOGY_KINDS
from standard_morph.models.qc_result import MetricResult
from standard_morph.registry import register

APICAL = 4


class ApicalOriginationMetric(Metric):
    name = "apical_origination"
    display_name = "Multiple apical dendrite origination points"
    metric_number = 14

    # Multiple apical trunks can be biologically plausible, so a violation is
    # flagged for human review rather than treated as an objective failure.
    violation_severity = Severity.REVIEW

    applicability = Applicability(
        spaces=ALL_COORDINATE_SPACES,
        morphology_kinds=ALL_MORPHOLOGY_KINDS,
        required_resources=frozenset(),
    )
    required_policy_keys = frozenset({"max_origins"})

    def evaluate(self, prepared_morph, context, policy):
        t0 = time.perf_counter()
        max_origins = policy[self.name, "max_origins"]

        parent = prepared_morph.parent
        compartment = prepared_morph.compartment

        result = MetricResult(name=self.name, status="pass")
        result.thresholds_used = {"max_origins": max_origins}
        result.value_label = "n_apical_origins"

        apical = np.flatnonzero(compartment == APICAL)
        if apical.size == 0:
            result.value = 0
            result.measurements = {"n_apical_origins": 0}
            result.counts = {"n_apical_origins": 0}
            result.message = "No apical dendrite nodes present."
            result.runtime_ms = (time.perf_counter() - t0) * 1000
            return result

        # Origins: apical nodes whose parent is not apical (incl. root/orphan).
        p = parent[apical]
        parent_is_apical = (p >= 0) & (compartment[np.where(p >= 0, p, 0)] == APICAL)
        origins = apical[~parent_is_apical]
        n_origins = int(origins.size)
        
        result.value = n_origins
        result.measurements = {"n_apical_origins": n_origins}
        result.counts = {"n_apical_origins": n_origins}

        if n_origins > max_origins:
            result.status = self.violation_severity.value
            result.flagged_node_ids = [int(prepared_morph.node_id[o]) for o in origins]
            result.flagged_node_coordinates = [
                tuple(float(v) for v in prepared_morph.xyz[o]) for o in origins
            ]
            result.message = (
                f"{n_origins} apical dendrite origination points (expected at most "
                f"{max_origins}); may be biologically plausible -- flag for review."
            )
        else:
            result.message = f"{n_origins} apical dendrite origination point(s)."

        result.runtime_ms = (time.perf_counter() - t0) * 1000
        return result


register(ApicalOriginationMetric())
