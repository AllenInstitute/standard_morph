"""Metric 6 -- Expected node identity types.

Every node's compartment code should be one of the expected SWC types
(soma=1, axon=2, basal dendrite=3, apical dendrite=4). A *merged* reconstruction
should additionally carry soma, axon, and dendrite identities together; a missing
identity can indicate an axon-only / dendrite-only / soma-only file, or axon and
dendrite traced with the same type.
"""
import time

import numpy as np

from standard_morph.metrics.base import Metric, Applicability
from standard_morph.models.qc_context import (
    ALL_COORDINATE_SPACES,
    ALL_MORPHOLOGY_KINDS,
    MorphologyKind,
)
from standard_morph.models.qc_result import MetricResult
from standard_morph.registry import register

SOMA, AXON, BASAL_DENDRITE, APICAL_DENDRITE = 1, 2, 3, 4


class NodeIdentityTypesMetric(Metric):
    name = "node_identity_types"
    display_name = "Expected node identity types"
    metric_number = 6
    requires_topology = False  # reads only the compartment column

    applicability = Applicability(
        spaces=ALL_COORDINATE_SPACES,
        morphology_kinds=ALL_MORPHOLOGY_KINDS,
        required_resources=frozenset(),
    )
    required_policy_keys = frozenset({"allowed_types"})

    def evaluate(self, prepared_morph, context, policy):
        t0 = time.perf_counter()
        allowed = set(policy[self.name, "allowed_types"])

        result = MetricResult(name=self.name, status="pass")

        types, counts = np.unique(prepared_morph.compartment, return_counts=True)
        present = {int(t) for t in types}
        type_counts = {int(t): int(c) for t, c in zip(types, counts)}

        # Nodes whose compartment code is outside the allowed set.
        unexpected_idx = np.flatnonzero(~np.isin(prepared_morph.compartment, list(allowed)))

        result.value = int(unexpected_idx.size)
        result.value_label = "n_unexpected_type_nodes"
        result.measurements = {
            "unique_types": sorted(present),
            "type_counts": type_counts,
            "n_unexpected_type_nodes": int(unexpected_idx.size),
        }

        problems = []
        if unexpected_idx.size:
            problems.append(f"unexpected node type(s) present: {sorted(present - allowed)}")

        # A merged file must carry soma, axon, and a dendrite identity.
        if context.morphology_kind == MorphologyKind.MERGED:
            missing = []
            if SOMA not in present:
                missing.append("soma(1)")
            if AXON not in present:
                missing.append("axon(2)")
            if BASAL_DENDRITE not in present and APICAL_DENDRITE not in present:
                missing.append("dendrite(3/4)")
            if missing:
                problems.append(
                    "merged file missing expected identities: "
                    + ", ".join(missing)
                    + " (axon and dendrite may have been traced with the same type)"
                )

        result.flagged_node_ids = [int(prepared_morph.node_id[i]) for i in unexpected_idx]
        result.flagged_node_coordinates = [
            tuple(float(v) for v in prepared_morph.xyz[i]) for i in unexpected_idx
        ]
        result.counts = {"n_unexpected_type_nodes": int(unexpected_idx.size)}

        if problems:
            result.status = self.violation_severity.value
            result.message = "; ".join(problems)
        else:
            result.message = f"Node types {sorted(present)} are all expected."

        result.runtime_ms = (time.perf_counter() - t0) * 1000
        return result


register(NodeIdentityTypesMetric())
