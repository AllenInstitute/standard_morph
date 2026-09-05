"""Metric 4/5 -- Edge length threshold.

The length of each edge (a node to its parent) should stay within a threshold:
roughly 30 um before resampling (image space) and 10 um after resampling (CCF
space). One metric covers both stages; the threshold comes from the policy,
keyed by coordinate space.

Edges from the soma to its immediate children are **excluded** here -- the soma
is modelled as a single point with a radius, but real somata are not spheres, so
a subtree's origin legitimately sits farther out than a normal edge. That
soma-to-child distance is checked separately, with its own larger threshold, by
``soma_child_distance``.
"""
import time

import numpy as np

from standard_morph.metrics.base import Metric, Applicability, threshold_for_space
from standard_morph.models.qc_context import ALL_COORDINATE_SPACES, ALL_MORPHOLOGY_KINDS
from standard_morph.models.qc_policy import PolicyRange
from standard_morph.models.qc_result import MetricResult
from standard_morph.registry import register


class EdgeLengthMetric(Metric):
    name = "edge_length"
    display_name = "Edge length within threshold"
    metric_number = 4  # 4 (pre-resampling) / 5 (post-resampling)

    applicability = Applicability(
        spaces=ALL_COORDINATE_SPACES,
        morphology_kinds=ALL_MORPHOLOGY_KINDS,
        required_resources=frozenset(),
    )
    required_policy_keys = frozenset({"max_length_um"})

    def evaluate(self, prepared_morph, context, policy):
        t0 = time.perf_counter()
        threshold = threshold_for_space(
            policy, self.name, "max_length_um", context.space.value
        )

        parent = prepared_morph.parent
        xyz = prepared_morph.xyz

        # Every real edge is a node with a valid parent index ...
        child = np.flatnonzero(parent >= 0)
        parent_idx = parent[child]
        # ... except edges whose parent is the soma (its immediate children).
        # "The soma" is a soma-typed root, not any type-1 node (see soma_roots).
        keep = ~np.isin(parent_idx, prepared_morph.soma_roots)
        child, parent_idx = child[keep], parent_idx[keep]

        is_range = isinstance(threshold, PolicyRange)
        threshold_repr = {"lo": threshold.lo, "hi": threshold.hi} if is_range else threshold

        result = MetricResult(name=self.name, status="pass")
        result.thresholds_used = {"max_length_um": threshold_repr, "space": context.space.value}
        result.value_label = "fraction_edges_out_of_range" if is_range else "fraction_edges_over_threshold"

        if child.size == 0:
            result.value = 0.0
            result.measurements = {"n_edges": 0, "n_over": 0, "fraction_over": 0.0, "max_edge_length": None}
            result.counts = {"n_over": 0, "n_edges": 0}
            result.message = "No non-soma edges to measure."
            result.runtime_ms = (time.perf_counter() - t0) * 1000
            return result

        lengths = np.linalg.norm(xyz[child] - xyz[parent_idx], axis=1)

        if is_range:
            violations = (lengths < threshold.lo) | (lengths > threshold.hi)
        else:
            violations = lengths > threshold

        over_idx = child[violations]

        result.value = float(violations.mean())
        result.measurements = {
            "n_edges": int(lengths.size),
            "n_over": int(violations.sum()),
            "fraction_over": float(violations.mean()),
            "max_edge_length": float(lengths.max()),
        }
        result.flagged_node_ids = [int(prepared_morph.node_id[i]) for i in over_idx]
        result.flagged_node_coordinates = [tuple(float(v) for v in xyz[i]) for i in over_idx]
        result.counts = {"n_over": int(violations.sum()), "n_edges": int(lengths.size)}

        if violations.any():
            result.status = self.violation_severity.value
            if is_range:
                result.message = (
                    f"{int(violations.sum())}/{lengths.size} edge(s) outside "
                    f"[{threshold.lo}, {threshold.hi}] um "
                    f"(min {lengths.min():.1f}, max {lengths.max():.1f} um)."
                )
            else:
                result.message = (
                    f"{int(violations.sum())}/{lengths.size} edge(s) exceed {threshold} um "
                    f"(max {lengths.max():.1f} um)."
                )
        else:
            if is_range:
                result.message = (
                    f"All {lengths.size} edge(s) within [{threshold.lo}, {threshold.hi}] um."
                )
            else:
                result.message = f"All {lengths.size} edge(s) within {threshold} um."

        result.runtime_ms = (time.perf_counter() - t0) * 1000
        return result


register(EdgeLengthMetric())
