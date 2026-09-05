"""Metric -- Soma is the first node.

By SWC convention the soma is written as the first data line of the file. This is
its own check rather than part of ``single_root_node`` because it is a distinct
concern (row *order*, not soma *uniqueness*) and it is coordinate/attribute-only:
it reads just the first row's ``compartment`` and ``parent``, so it sets
``requires_topology=False`` and still runs when the tree wiring is untrustworthy
(duplicate ids, a cycle, ...), unlike ``single_root_node``.

Index 0 is the first data row of the SWC: ``PreparedMorphology`` preserves the
input row order, and ``read_swc`` preserves file order (comments stripped).
"""
import time

from standard_morph.metrics.base import Metric, Applicability
from standard_morph.models.qc_context import ALL_COORDINATE_SPACES, ALL_MORPHOLOGY_KINDS
from standard_morph.models.qc_result import MetricResult
from standard_morph.registry import register


class SomaFirstNodeMetric(Metric):
    name = "soma_first_node"
    display_name = "Soma is the first node"
    metric_number = None
    requires_topology = False  # reads only the first row's compartment + parent

    applicability = Applicability(
        spaces=ALL_COORDINATE_SPACES,
        morphology_kinds=ALL_MORPHOLOGY_KINDS,
        required_resources=frozenset(),
    )

    def evaluate(self, prepared_morph, context, policy):
        t0 = time.perf_counter()
        result = MetricResult(name=self.name, status="pass")

        # The soma is the first line iff row 0 is a soma root (type 1, parent -1).
        first_is_soma = bool(
            (prepared_morph.compartment[0] == prepared_morph.SOMA_COMPARTMENT)
            and (prepared_morph.parent[0] == prepared_morph.ROOT)
        )
        first_id = int(prepared_morph.node_id[0])
        first_type = int(prepared_morph.compartment[0])

        result.value = first_is_soma
        result.value_label = "soma_is_first_node"
        result.measurements = {
            "first_node_id": first_id,
            "first_node_compartment": first_type,
            "soma_is_first_node": first_is_soma,
        }

        if first_is_soma:
            result.message = f"First node (id {first_id}) is the soma."
        else:
            result.status = self.violation_severity.value
            result.flagged_node_ids = [first_id]
            result.flagged_node_coordinates = [tuple(float(v) for v in prepared_morph.xyz[0])]
            result.message = (
                f"First node (id {first_id}, type {first_type}) is not the soma; "
                "the soma (type 1, root) should be the first line of the SWC."
            )

        result.runtime_ms = (time.perf_counter() - t0) * 1000
        return result


register(SomaFirstNodeMetric())
