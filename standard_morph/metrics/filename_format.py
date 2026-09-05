"""Metric 7 -- Filename format.

Checks that an SWC file is named according to the expected convention for its
source. The convention is selected by policy (``name_format``); currently:

* ``AIND`` -- Allen Institute for Neural Dynamics. Validated against the naming
  pattern ported from the legacy ``tools.has_valid_name``.
* ``AIBS`` -- Allen Institute for Brain Science. TODO: the pattern is not yet
  defined, so this always passes.

Unlike every other metric, this one inspects the *filename*, not the SWC data --
so it reads the name from ``context.resources["filename"]`` and ignores the
table it is handed. When ``run_qc`` is given a file path, the engine populates
``resources["filename"]`` automatically from the path's basename. For DataFrame
or PreparedMorphology inputs the caller must supply it in ``resources`` (or omit
this metric); if it is absent the metric is reported as incompatible with the
context (fail-fast) -- the same contract the CCF metrics use for
``ccf_resolution``.
"""
import os
import re
import time

from standard_morph.metrics.base import Metric, Applicability, EvaluationPhase
from standard_morph.models.qc_context import ALL_COORDINATE_SPACES, ALL_MORPHOLOGY_KINDS
from standard_morph.models.qc_result import MetricResult
from standard_morph.registry import register

#: AIND SWC naming pattern (ported from ``tools.has_valid_name``). Matches, e.g.:
#:   N123-000000.swc
#:   N1_123456_ABC.swc                (2-3 letter initials, or "consensus")
#:   N42-654321-consensus.swc
#:   N7-123456-axon-XY.swc            (axon/dendrite tag + initials/consensus)
#:   N7_123456_dendrite_consensus.swc
AIND_PATTERN = (
    r"^N\d+[-_]\d{6}"
    r"(?:[-_](?:[A-Za-z]{2,3}|consensus)"
    r"|[-_](?:axon|dendrite)[-_](?:[A-Za-z]{2,3}|consensus))?"
    r"\.swc$"
)


def _is_valid_filename(filename, name_format):
    """Return True if ``filename`` matches the ``name_format`` convention."""
    if name_format == "AIND":
        return re.match(AIND_PATTERN, filename, re.IGNORECASE) is not None
    if name_format == "AIBS":
        # TODO: define and validate the AIBS naming pattern. Passes for now.
        return True
    raise ValueError(
        f"Unknown name_format {name_format!r}. Supported formats: 'AIND', 'AIBS'."
    )


class FilenameFormatMetric(Metric):
    name = "filename_format"
    display_name = "Filename format"
    metric_number = 7

    evaluation_phase = EvaluationPhase.INPUT_INTEGRITY
    # A bad filename is worth flagging, but it does not prevent building or
    # trusting the morphology, so its failure blocks nothing.
    blocks_on_failure = None

    applicability = Applicability(
        spaces=ALL_COORDINATE_SPACES,
        morphology_kinds=ALL_MORPHOLOGY_KINDS,
        required_resources=frozenset({"filename"}),
    )
    required_policy_keys = frozenset({"name_format"})

    def evaluate(self, swc_df, context, policy):
        t0 = time.perf_counter()
        name_format = policy[self.name, "name_format"]
        filename = os.path.basename(context.resources["filename"])

        result = MetricResult(name=self.name, status="pass")
        is_valid = _is_valid_filename(filename, name_format)

        result.value = bool(is_valid)
        result.value_label = "filename_valid"
        result.measurements = {"filename": filename, "name_format": name_format}
        if is_valid:
            result.message = f"Filename {filename!r} matches the {name_format} convention."
        else:
            result.status = self.violation_severity.value
            result.message = (
                f"Filename {filename!r} does not match the {name_format} naming convention."
            )
        result.runtime_ms = (time.perf_counter() - t0) * 1000
        return result


register(FilenameFormatMetric())
