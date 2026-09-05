"""Metric 15 -- Compartment type consistency (dendrites).

Every dendrite node's parent must be a compatible compartment -- a transition
may only happen where it makes biological sense:

* basal dendrite (type 3)  -> parent is the **soma root** or another basal (3)
* apical dendrite (type 4)  -> parent is the **soma root** or another apical (4)

"Soma root" is the canonical soma (type 1 **and** a root), not just any type-1
node. That distinction matters: a dendrite hanging off a *stray* mid-tree
type-1 node (e.g. a second soma erroneously placed inside a dendrite) is an
error, and it must be caught here -- especially for basal dendrites, where
multiple trunks from the real soma are normal, so there is no "too many origins"
count metric to catch it (unlike ``apical_origination``, metric 14).

Checking every dendrite edge also enforces two things for free: dendrites
originate from valid places, and the axon tree stays pure -- a dendrite hanging
off an axon fails its own rule. The axon is validated separately by
``axon_origination`` (metric 13) and the soma by ``single_root_node`` (metric 1),
so together the three cover every edge.
"""
import time

import numpy as np

from standard_morph.metrics.base import Metric, Applicability
from standard_morph.models.qc_context import ALL_COORDINATE_SPACES, ALL_MORPHOLOGY_KINDS
from standard_morph.models.qc_result import MetricResult
from standard_morph.registry import register

BASAL, APICAL = 3, 4


class CompartmentTransitionsMetric(Metric):
    name = "compartment_transitions"
    display_name = "Compartment type consistency"
    metric_number = 15

    applicability = Applicability(
        spaces=ALL_COORDINATE_SPACES,
        morphology_kinds=ALL_MORPHOLOGY_KINDS,
        required_resources=frozenset(),
    )

    def evaluate(self, prepared_morph, context, policy):
        t0 = time.perf_counter()
        result = MetricResult(name=self.name, status="pass")

        parent = prepared_morph.parent
        compartment = prepared_morph.compartment
        has_parent = parent >= 0
        soma_roots = prepared_morph.soma_roots

        flagged = []
        for dendrite_type in (BASAL, APICAL):
            idx = np.flatnonzero((compartment == dendrite_type) & has_parent)
            if idx.size:
                parent_idx = parent[idx]
                # Valid parent: the canonical soma root, or a node of the same
                # dendrite type. A type-1 node that is NOT the soma root does not
                # count -- a dendrite off a stray mid-tree soma is an error.
                valid = np.isin(parent_idx, soma_roots) | (compartment[parent_idx] == dendrite_type)
                flagged.extend(int(i) for i in idx[~valid])

        flagged_idx = sorted(flagged)
        result.value = len(flagged_idx)
        result.value_label = "n_invalid_transitions"
        result.measurements = {"n_invalid_transitions": len(flagged_idx)}
        result.flagged_node_ids = [int(prepared_morph.node_id[i]) for i in flagged_idx]
        result.flagged_node_coordinates = [
            tuple(float(v) for v in prepared_morph.xyz[i]) for i in flagged_idx
        ]
        result.counts = {"n_invalid_transitions": len(flagged_idx)}

        if flagged_idx:
            result.status = self.violation_severity.value
            result.message = (
                f"{len(flagged_idx)} dendrite node(s) originate from an invalid "
                "compartment (not the soma root or same-type dendrite)."
            )
        else:
            result.message = "All dendrite compartment transitions are valid."

        result.runtime_ms = (time.perf_counter() - t0) * 1000
        return result


register(CompartmentTransitionsMetric())
