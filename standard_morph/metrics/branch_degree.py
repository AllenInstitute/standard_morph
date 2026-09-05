"""Metric 17 -- Branch points with more than two children.

Bifurcations are expected; branch points with more than ``max_children``
children are flagged for interrogation. The soma is excluded, since it
legitimately has many primary neurites radiating from it.
"""
import time

import numpy as np

from standard_morph.metrics.base import Metric, Applicability, Severity
from standard_morph.models.qc_context import ALL_COORDINATE_SPACES, ALL_MORPHOLOGY_KINDS
from standard_morph.models.qc_result import MetricResult
from standard_morph.registry import register


class BranchMaxDegreeMetric(Metric):
    name = "branch_max_degree"
    display_name = "Branch points with more than two children"
    metric_number = 17

    # A branch point with >2 children is unusual but can be biologically real, so
    # a violation is flagged for human review rather than an objective failure.
    violation_severity = Severity.REVIEW

    # Topological check; valid in either coordinate space, any morphology kind.
    applicability = Applicability(
        spaces=ALL_COORDINATE_SPACES,
        morphology_kinds=ALL_MORPHOLOGY_KINDS,
        required_resources=frozenset(),
    )
    required_policy_keys = frozenset({"max_children"})

    def evaluate(self, prepared_morph, context, policy):
        t0 = time.perf_counter()
        max_children = policy[self.name, "max_children"]

        result = MetricResult(name=self.name, status="pass")
        result.thresholds_used = {"max_children": max_children}

        # The soma (a soma-typed root) legitimately branches many times, so it is
        # excluded -- but only the actual soma, not any stray type-1 node.
        soma_roots = prepared_morph.soma_roots
        over = prepared_morph.child_counts > max_children
        over[soma_roots] = False
        flagged_idx = np.flatnonzero(over)

        # Max degree excluding the soma, for context in the report.
        non_soma = np.ones(prepared_morph.n, dtype=bool)
        non_soma[soma_roots] = False
        non_soma_counts = prepared_morph.child_counts[non_soma]
        result.measurements = {
            "max_children_observed": int(non_soma_counts.max()) if non_soma_counts.size else 0,
        }
        result.value = result.measurements["max_children_observed"]
        result.value_label = "max_children_observed"
        result.flagged_node_ids = [int(prepared_morph.node_id[i]) for i in flagged_idx]
        result.flagged_node_coordinates = [
            tuple(float(v) for v in prepared_morph.xyz[i]) for i in flagged_idx
        ]
        result.counts = {"n_flagged": int(flagged_idx.size)}

        if flagged_idx.size:
            result.status = self.violation_severity.value
            result.message = (
                f"{flagged_idx.size} branch point(s) have more than "
                f"{max_children} children."
            )
        else:
            result.message = f"No branch points exceed {max_children} children."

        result.runtime_ms = (time.perf_counter() - t0) * 1000
        return result


register(BranchMaxDegreeMetric())
