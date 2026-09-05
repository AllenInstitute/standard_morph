"""Metric 1 -- Single root node validation.

A valid SWC has exactly one root -- a node with ``parent == -1`` -- and that root
is the soma (node type 1), with no other node of type 1, and by convention its
``node_id`` is 1. (The soma's *first-line* placement is a separate concern,
checked by ``soma_first_node``.)
"""
import time

import numpy as np

from standard_morph.metrics.base import Metric, Applicability
from standard_morph.models.qc_context import ALL_COORDINATE_SPACES, ALL_MORPHOLOGY_KINDS
from standard_morph.models.qc_result import MetricResult
from standard_morph.registry import register


class SingleRootNodeMetric(Metric):
    name = "single_root_node"
    display_name = "Single root node validation"
    metric_number = 1

    applicability = Applicability(
        spaces=ALL_COORDINATE_SPACES,
        morphology_kinds=ALL_MORPHOLOGY_KINDS,
        required_resources=frozenset(),
    )

    def evaluate(self, prepared_morph, context, policy):
        t0 = time.perf_counter()
        result = MetricResult(name=self.name, status="pass")

        # The soma root (type 1 AND parent -1) is the canonical property; this
        # metric is its validator, so it also checks the two components -- that
        # there is exactly one root of any type, and exactly one type-1 node.
        soma_root_idx = prepared_morph.soma_roots
        root_idx = prepared_morph.roots
        type1_idx = np.flatnonzero(prepared_morph.compartment == prepared_morph.SOMA_COMPARTMENT)

        n_soma_roots = int(soma_root_idx.size)
        n_roots = int(root_idx.size)
        n_type1 = int(type1_idx.size)

        result.value = n_roots
        result.value_label = "n_roots"
        result.measurements = {
            "n_soma_roots": n_soma_roots,
            "n_roots": n_roots,
            "n_type1_nodes": n_type1,
        }

        problems = []
        flagged = set()

        # Core requirement, checked directly: exactly one node that is BOTH the
        # soma (type 1) and the root (parent == -1). A separate type-1 node and a
        # separate parent==-1 node therefore cannot slip through as two nodes.
        if n_soma_roots != 1:
            problems.append(
                f"expected exactly one soma root (type 1 with parent == -1), found {n_soma_roots}"
            )
            flagged.update(int(i) for i in root_idx)
            flagged.update(int(i) for i in type1_idx)

        # No extra roots of any type ...
        if n_roots != 1:
            problems.append(f"expected exactly one root (parent == -1), found {n_roots}")
            flagged.update(int(i) for i in root_idx)

        # ... and no extra type-1 (soma) nodes anywhere.
        if n_type1 != 1:
            problems.append(f"expected exactly one type-1 (soma) node, found {n_type1}")
            flagged.update(int(i) for i in type1_idx)

        # Canonical id of the (unique) soma root. Its first-line placement is
        # checked separately by `soma_first_node`.
        if n_soma_roots == 1:
            r = int(soma_root_idx[0])
            if int(prepared_morph.node_id[r]) != 1:
                problems.append(f"soma root node_id is {int(prepared_morph.node_id[r])}, expected 1")
                flagged.add(r)

        flagged_idx = sorted(flagged)
        result.flagged_node_ids = [int(prepared_morph.node_id[i]) for i in flagged_idx]
        result.flagged_node_coordinates = [
            tuple(float(v) for v in prepared_morph.xyz[i]) for i in flagged_idx
        ]
        result.counts = {
            "n_roots": n_roots,
            "n_type1_nodes": n_type1,
            "n_flagged": len(flagged_idx),
        }

        if problems:
            result.status = self.violation_severity.value
            result.message = "; ".join(problems)
        else:
            result.message = "Single valid root: one soma-type root, node_id 1."

        result.runtime_ms = (time.perf_counter() - t0) * 1000
        return result


register(SingleRootNodeMetric())
