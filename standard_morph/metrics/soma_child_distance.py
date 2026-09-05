"""Metric 4/5 companion -- Soma-to-child distance.

The soma's immediate children -- the origin of each subtree -- sit farther from
the soma point than a normal edge would, because the soma is modelled as a
single point plus a radius while real somata are not spheres. So they are
excluded from ``edge_length``. They should still not be *too* far, though, which
is what this metric checks, with its own (larger) threshold.

"Soma children" are nodes whose parent's compartment is soma (type 1); this does
not depend on there being exactly one soma (that is ``single_root_node``'s job).
"""
import time

import numpy as np

from standard_morph.metrics.base import Metric, Applicability
from standard_morph.models.qc_context import ALL_COORDINATE_SPACES, ALL_MORPHOLOGY_KINDS
from standard_morph.models.qc_result import MetricResult
from standard_morph.registry import register


class SomaChildDistanceMetric(Metric):
    name = "soma_child_distance"
    display_name = "Soma-to-child distance within threshold"
    metric_number = None  # a companion to metrics 4/5, not a numbered inventory item

    applicability = Applicability(
        spaces=ALL_COORDINATE_SPACES,
        morphology_kinds=ALL_MORPHOLOGY_KINDS,
        required_resources=frozenset(),
    )
    required_policy_keys = frozenset({"max_soma_child_to_soma_um"})

    def evaluate(self, prepared_morph, context, policy):
        t0 = time.perf_counter()
        threshold = policy[self.name, "max_soma_child_to_soma_um"]

        parent = prepared_morph.parent
        xyz = prepared_morph.xyz

        child = np.flatnonzero(parent >= 0)
        parent_idx = parent[child]
        # Children of a soma root (type 1 AND parent -1), not of any type-1 node.
        is_soma_child = np.isin(parent_idx, prepared_morph.soma_roots)
        soma_child, soma_idx = child[is_soma_child], parent_idx[is_soma_child]

        result = MetricResult(name=self.name, status="pass")
        result.thresholds_used = {"max_soma_child_to_soma_um": threshold}
        result.value_label = "max_soma_child_distance_um"

        if soma_child.size == 0:
            result.value = None
            result.measurements = {"n_soma_children": 0, "max_soma_child_distance_um": None}
            result.counts = {"n_soma_children": 0, "n_over": 0}
            result.message = "No soma children to measure."
            result.runtime_ms = (time.perf_counter() - t0) * 1000
            return result

        distances = np.linalg.norm(xyz[soma_child] - xyz[soma_idx], axis=1)
        over = distances > threshold
        over_idx = soma_child[over]

        result.value = float(distances.max())
        result.measurements = {
            "n_soma_children": int(soma_child.size),
            "n_over": int(over.sum()),
            "max_soma_child_distance_um": float(distances.max()),
            "soma_child_distances": [float(d) for d in distances],
        }
        result.flagged_node_ids = [int(prepared_morph.node_id[i]) for i in over_idx]
        result.flagged_node_coordinates = [tuple(float(v) for v in xyz[i]) for i in over_idx]
        result.counts = {"n_soma_children": int(soma_child.size), "n_over": int(over.sum())}

        if over.any():
            result.status = self.violation_severity.value
            result.message = (
                f"{int(over.sum())}/{soma_child.size} soma child(ren) farther than "
                f"{threshold} um (max {distances.max():.1f} um)."
            )
        else:
            result.message = f"All {soma_child.size} soma child(ren) within {threshold} um."

        result.runtime_ms = (time.perf_counter() - t0) * 1000
        return result


register(SomaChildDistanceMetric())
