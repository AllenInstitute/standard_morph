"""Metric 11 -- Duplicate node coordinates.

Flags nodes that share an identical ``(x, y, z)`` coordinate with another node.
Exact equality is used because duplicates arise as literally repeated rows.
"""
import time

import numpy as np

from standard_morph.metrics.base import Metric, Applicability
from standard_morph.models.qc_context import ALL_COORDINATE_SPACES, ALL_MORPHOLOGY_KINDS
from standard_morph.models.qc_result import MetricResult
from standard_morph.registry import register


class DuplicateNodeCoordinatesMetric(Metric):
    name = "duplicate_node_coordinates"
    display_name = "Duplicate node coordinates"
    metric_number = 11
    requires_topology = False  # groups on xyz only; reads no tree structure

    applicability = Applicability(
        spaces=ALL_COORDINATE_SPACES,
        morphology_kinds=ALL_MORPHOLOGY_KINDS,
        required_resources=frozenset(),
    )

    def evaluate(self, prepared_morph, context, policy):
        t0 = time.perf_counter()
        result = MetricResult(name=self.name, status="pass")

        # Group nodes by exact coordinate. Any coordinate held by 2+ nodes is a
        # duplicate group; every node in such a group is flagged.
        _, inverse, counts = np.unique(
            prepared_morph.xyz, axis=0, return_inverse=True, return_counts=True
        )
        inverse = inverse.ravel()  # numpy version-robust (some return an (N,1) column)
        duplicated_group = counts > 1
        dup_idx = np.flatnonzero(duplicated_group[inverse])

        n_dup_nodes = int(dup_idx.size)
        n_dup_groups = int(duplicated_group.sum())

        result.value = n_dup_nodes
        result.value_label = "n_duplicate_nodes"
        result.measurements = {
            "n_duplicate_nodes": n_dup_nodes,
            "n_duplicate_groups": n_dup_groups,
        }
        result.flagged_node_ids = [int(prepared_morph.node_id[i]) for i in dup_idx]
        result.flagged_node_coordinates = [
            tuple(float(v) for v in prepared_morph.xyz[i]) for i in dup_idx
        ]
        result.counts = {"n_duplicate_nodes": n_dup_nodes, "n_duplicate_groups": n_dup_groups}

        if n_dup_nodes:
            result.status = self.violation_severity.value
            result.message = (
                f"{n_dup_nodes} node(s) in {n_dup_groups} group(s) share coordinates "
                "with another node."
            )
        else:
            result.message = "No duplicate node coordinates."

        result.runtime_ms = (time.perf_counter() - t0) * 1000
        return result


register(DuplicateNodeCoordinatesMetric())
