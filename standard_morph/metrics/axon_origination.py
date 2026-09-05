"""Metric 13 -- Axon origination.

The axon should enter the reconstruction at a single, well-placed point. This
checks three things about where the axon originates -- a node whose parent is
not itself an axon:

1. there is **exactly one** such origination point;
2. its parent is the **soma or a basal dendrite** (type 1 or 3); and
3. the origination sits **near the soma**.

The distance in (3) is the straight-line distance from the axon-origin node to
the *soma point* -- not to its immediate parent. The axon may legitimately stem
from a basal dendrite node, but that node (and hence the origin) should still be
close to the cell body, not out on a distal dendrite.

Note this ``max_axon_origin_to_soma_um`` threshold is its own policy value, larger
than ``soma_child_distance``'s: an axon origin may sit a hop or few out from the
soma, whereas a soma *child* is immediate.
"""
import time

import numpy as np

from standard_morph.metrics.base import Metric, Applicability
from standard_morph.models.qc_context import ALL_COORDINATE_SPACES, ALL_MORPHOLOGY_KINDS
from standard_morph.models.qc_result import MetricResult
from standard_morph.registry import register

AXON, SOMA, BASAL = 2, 1, 3


class AxonOriginationMetric(Metric):
    name = "axon_origination"
    display_name = "Axon origination"
    metric_number = 13

    applicability = Applicability(
        spaces=ALL_COORDINATE_SPACES,
        morphology_kinds=ALL_MORPHOLOGY_KINDS,
        required_resources=frozenset(),
    )
    required_policy_keys = frozenset({"max_axon_origin_to_soma_um"})

    def evaluate(self, prepared_morph, context, policy):
        t0 = time.perf_counter()
        max_distance = policy[self.name, "max_axon_origin_to_soma_um"]

        parent = prepared_morph.parent
        compartment = prepared_morph.compartment
        xyz = prepared_morph.xyz

        result = MetricResult(name=self.name, status="pass")
        result.thresholds_used = {"max_axon_origin_to_soma_um": max_distance}
        result.value_label = "axon_origin_distance_to_soma_um"

        axon = np.flatnonzero(compartment == AXON)
        if axon.size == 0:
            result.value = None
            result.measurements = {"n_axon_origins": 0}
            result.counts = {"n_axon_origins": 0}
            result.message = "No axon nodes present."
            result.runtime_ms = (time.perf_counter() - t0) * 1000
            return result

        # Axon origins: axon nodes whose parent is NOT an axon (incl. root/orphan).
        p = parent[axon]
        parent_is_axon = (p >= 0) & (compartment[np.where(p >= 0, p, 0)] == AXON)
        origins = axon[~parent_is_axon]
        n_origins = int(origins.size)

        problems = []
        flagged = set()

        # (1) exactly one origination point
        if n_origins != 1:
            problems.append(f"expected exactly one axon origination point, found {n_origins}")
            flagged.update(int(i) for i in origins)

        # (2) valid parent type: soma or basal
        for o in origins:
            po = int(parent[o])
            if po < 0:
                problems.append("axon origination has no parent (should stem from soma or basal)")
                flagged.add(int(o))
            elif int(compartment[po]) not in (SOMA, BASAL):
                problems.append(
                    f"axon originates from compartment {int(compartment[po])}, expected soma(1) or basal(3)"
                )
                flagged.add(int(o))

        # (3) distance from each origin to the soma point (needs a unique soma)
        soma_roots = prepared_morph.soma_roots
        origin_distances = None
        if soma_roots.size == 1 and n_origins:
            soma_xyz = xyz[int(soma_roots[0])]
            distances = np.linalg.norm(xyz[origins] - soma_xyz, axis=1)
            origin_distances = [float(d) for d in distances]
            result.value = float(distances.max())
            for o, d in zip(origins, distances):
                if d > max_distance:
                    problems.append(f"axon origination is {d:.1f} um from the soma, exceeds {max_distance} um")
                    flagged.add(int(o))

        result.measurements = {
            "n_axon_origins": n_origins,
            "origin_distances_to_soma_um": origin_distances,
        }
        flagged_idx = sorted(flagged)
        result.flagged_node_ids = [int(prepared_morph.node_id[i]) for i in flagged_idx]
        result.flagged_node_coordinates = [
            tuple(float(v) for v in xyz[i]) for i in flagged_idx
        ]
        result.counts = {"n_axon_origins": n_origins, "n_flagged": len(flagged_idx)}

        if problems:
            result.status = self.violation_severity.value
            result.message = "; ".join(problems)
        else:
            result.message = "Single axon origination from the soma/basal, near the soma."

        result.runtime_ms = (time.perf_counter() - t0) * 1000
        return result


register(AxonOriginationMetric())
