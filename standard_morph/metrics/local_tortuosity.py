"""Metric 10 -- Large sidestep / local tortuosity deviations.

Measures local tortuosity as the ratio of local path length to the straight
chord over a 3-node window: one node above (the parent) and one node below.

* Reducible nodes (exactly one child): a single ``(parent, node, child)`` window.
* Branch nodes (two or more children): one window per child --
  ``(parent, node, child_k)`` -- averaged into a single tortuosity for the node.
* Tip nodes (no children) and root/orphan nodes (no valid parent -- which
  includes the soma) are not measured.

A straight run scores 1.0; a sharp sidestep or hairpin scores higher. Nodes
whose tortuosity exceeds the policy threshold are flagged.

v1 measures geometry only. The ``use_image_reference`` flag is reserved for a
later version that will consult image data to avoid flagging genuine sharp
turns in real neurons; when False the metric is purely geometric.
"""
import time

import numpy as np

from standard_morph.metrics.base import Metric, Applicability
from standard_morph.models.qc_context import ALL_COORDINATE_SPACES, ALL_MORPHOLOGY_KINDS
from standard_morph.models.qc_result import MetricResult
from standard_morph.registry import register

EPS = 1e-9


def _window_tortuosity(prev_xyz, node_xyz, next_xyz):
    """Local tortuosity for each 3-point window, vectorized over many windows.

    A "window" is three consecutive points along the neuron path: the node being
    scored (``node_xyz``) together with its neighbour one step above
    (``prev_xyz``) and one step below (``next_xyz``). Each argument is an
    ``(N, 3)`` array of x/y/z coordinates, so ``N`` windows are evaluated in a
    single call -- row ``i`` of every array belongs to window ``i``.

    Tortuosity is the ratio of the distance actually travelled along the path to
    the straight-line shortcut between the two outer points::

        tortuosity = (|prev -> node| + |node -> next|) / |prev -> next|

    A perfectly straight run scores 1.0 (path length == chord). The further the
    middle node juts off the line between its neighbours -- a "sidestep" -- the
    larger the ratio grows.

    Returns
    -------
    numpy.ndarray
        ``(N,)`` array of tortuosity values, one per window.
    """
    prev_xyz = np.asarray(prev_xyz, dtype=float)
    node_xyz = np.asarray(node_xyz, dtype=float)
    next_xyz = np.asarray(next_xyz, dtype=float)

    # The two legs of the path actually walked, as per-row (per-window) lengths.
    leg_in = np.linalg.norm(node_xyz - prev_xyz, axis=1)   # |prev -> node|
    leg_out = np.linalg.norm(next_xyz - node_xyz, axis=1)  # |node -> next|
    path_length = leg_in + leg_out

    # The straight-line shortcut that skips the middle node.
    euclid_dist = np.linalg.norm(next_xyz - prev_xyz, axis=1)  # |prev -> next|

    # When prev and next coincide the chord is ~0 and the ratio is undefined /
    # blows up: this is a hairpin (the path folds back on itself). Treat it as
    # infinite tortuosity so those nodes are always flagged, and silence the
    # divide-by-zero warning that np.where would otherwise still emit.
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(euclid_dist > EPS, path_length / euclid_dist, np.inf)


class LocalTortuosityMetric(Metric):
    name = "local_tortuosity"
    display_name = "Large sidestep / local tortuosity deviations"
    metric_number = 10

    applicability = Applicability(
        spaces=ALL_COORDINATE_SPACES,
        morphology_kinds=ALL_MORPHOLOGY_KINDS,
        required_resources=frozenset(),
    )
    required_policy_keys = frozenset({"tortuosity_threshold"})

    def __init__(self, use_image_reference: bool = False):
        self.use_image_reference = use_image_reference

    def evaluate(self, prepared_morph, context, policy):
        t0 = time.perf_counter()
        threshold = policy[self.name, "tortuosity_threshold"]

        result = MetricResult(name=self.name, status="pass")
        result.thresholds_used = {"tortuosity_threshold": threshold}
        result.value_label = "max_tortuosity"
        # `measurements` starts empty (dataclass default) and is filled with
        # measurements below. `use_image_reference` is metric *config*, not a
        # measurement, so it does not belong here.

        xyz = prepared_morph.xyz
        parent = prepared_morph.parent
        child_counts = prepared_morph.child_counts
        # A node needs a node above (valid parent) and at least one below (child).
        # `parent >= 0` excludes roots (-1) and orphans (-2) -- including the soma.
        has_parent = parent >= 0

        measured_idx = []   # middle-node indices that got a tortuosity value
        measured_tort = []  # matching tortuosity values (same order)

        # -- Reducible nodes: exactly one (parent, node, child) window each. --
        reducible = np.flatnonzero((child_counts == 1) & has_parent)
        if reducible.size:
            # A reducible node's single child is the unique node whose parent is it.
            child_of = np.full(prepared_morph.n, -1, dtype=np.int64)
            valid = np.flatnonzero(has_parent)
            child_of[parent[valid]] = valid  # reducible parents have exactly one child
            nxt = child_of[reducible]
            measured_idx.append(reducible)
            measured_tort.append(
                _window_tortuosity(xyz[parent[reducible]], xyz[reducible], xyz[nxt])
            )

        # -- Branch nodes: one window per child, averaged into one value. --
        branches = np.flatnonzero((child_counts >= 2) & has_parent)
        if branches.size:
            prev_i, mid_i, child_i, owner = [], [], [], []
            for j, b in enumerate(branches):
                b = int(b)
                kids = prepared_morph.children(b)
                prev_i.extend([int(parent[b])] * len(kids))
                mid_i.extend([b] * len(kids))
                child_i.extend(kids)
                owner.extend([j] * len(kids))
            window_tort = _window_tortuosity(xyz[prev_i], xyz[mid_i], xyz[child_i])
            owner = np.asarray(owner)
            # Mean tortuosity across each branch's per-child windows.
            sums = np.bincount(owner, weights=window_tort, minlength=branches.size)
            counts = np.bincount(owner, minlength=branches.size)
            measured_idx.append(branches)
            measured_tort.append(sums / counts)

        if not measured_idx:
            result.message = "No measurable nodes available to evaluate tortuosity."
            result.measurements.update(
                {"n_evaluated": 0, "max_tortuosity": None, "mean_tortuosity": None}
            )
            result.counts = {"n_flagged": 0, "n_evaluated": 0}
            result.runtime_ms = (time.perf_counter() - t0) * 1000
            return result

        mid = np.concatenate(measured_idx)
        tort = np.concatenate(measured_tort)

        flagged_mask = tort > threshold
        flagged_idx = mid[flagged_mask]
        flagged_tort = tort[flagged_mask]

        finite = tort[np.isfinite(tort)]
        result.measurements.update(
            {
                "n_evaluated": int(tort.size),
                "max_tortuosity": (float(np.max(finite)) if finite.size else float("inf")),
                "mean_tortuosity": (float(np.mean(finite)) if finite.size else None),
                "flagged_tortuosities": [float(v) for v in flagged_tort],
            }
        )
        result.value = result.measurements["max_tortuosity"]
        result.flagged_node_ids = [int(prepared_morph.node_id[i]) for i in flagged_idx]
        result.flagged_node_coordinates = [
            tuple(float(v) for v in prepared_morph.xyz[i]) for i in flagged_idx
        ]
        n_flagged = int(flagged_mask.sum())
        result.counts = {"n_flagged": n_flagged, "n_evaluated": int(tort.size)}

        if n_flagged:
            result.status = self.violation_severity.value
            result.message = (
                f"{n_flagged} node(s) exceed local tortuosity threshold {threshold}."
            )
        else:
            result.message = (
                f"All {tort.size} evaluated node(s) within tortuosity threshold {threshold}."
            )

        result.runtime_ms = (time.perf_counter() - t0) * 1000
        return result


register(LocalTortuosityMetric())
