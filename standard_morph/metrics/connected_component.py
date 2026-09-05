"""Metric 16 -- Single connected component.

Each SWC must be a single connected tree. This fails when the reconstruction
splits into multiple components -- typically extra root nodes (more than one
``parent == -1``) or orphan subtrees whose parent id is missing from the file.
"""
import time

from standard_morph.metrics.base import Metric, Applicability
from standard_morph.models.qc_context import Space, ALL_MORPHOLOGY_KINDS
from standard_morph.models.qc_result import MetricResult
from standard_morph.registry import register


class SingleConnectedComponentMetric(Metric):
    name = "single_connected_component"
    display_name = "Single connected component"
    metric_number = 16

    # Topological check; valid in either space, any morphology kind.
    applicability = Applicability(
        spaces=frozenset({Space.IMAGE_SPACE, Space.CCF_REGISTERED}),
        morphology_kinds=ALL_MORPHOLOGY_KINDS,
        required_resources=frozenset(),
    )

    def evaluate(self, prepared_morph, context, policy):
        t0 = time.perf_counter()
        result = MetricResult(name=self.name, status="pass")

        n_components = prepared_morph.n_components
        roots = prepared_morph.roots
        orphans = prepared_morph.orphans

        result.measurements = {
            "n_components": n_components,
            "n_roots": int(roots.size),
            "n_orphans": int(orphans.size),
        }
        result.value = n_components
        result.value_label = "n_components"
        result.counts = {
            "n_components": n_components,
            "n_roots": int(roots.size),
            "n_orphans": int(orphans.size),
        }

        if n_components > 1:
            # Report the component "seeds": every root plus every orphan node.
            flagged = sorted(set(roots.tolist()) | set(orphans.tolist()))
            result.status = self.violation_severity.value
            result.flagged_node_ids = [int(prepared_morph.node_id[i]) for i in flagged]
            result.flagged_node_coordinates = [
                tuple(float(v) for v in prepared_morph.xyz[i]) for i in flagged
            ]
            result.message = (
                f"Reconstruction has {n_components} connected components "
                f"({roots.size} root(s), {orphans.size} orphan(s)); expected 1."
            )
        else:
            result.message = "Reconstruction is a single connected component."

        result.runtime_ms = (time.perf_counter() - t0) * 1000
        return result


register(SingleConnectedComponentMetric())
