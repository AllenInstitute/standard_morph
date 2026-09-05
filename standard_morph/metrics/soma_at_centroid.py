"""Metric -- Soma at image centroid.

Compares the reconstruction's soma position against a user-provided, image-based
soma centroid and per-axis soma radii. The SWC soma should sit NEAR the centroid,
not merely inside the soma extent -- so the offset is measured as a fraction of
the soma radius on each axis, and the *worst axis* must stay within
``max_offset_fraction`` (default 0.5 -> within half a radius on every axis).

    per_axis_fraction = |swc_soma - image_soma| / soma_radius   (elementwise)
    value = max(per_axis_fraction)                              # worst axis

Inputs are per-file, supplied via ``context.resources`` (required -> the metric
is incompatible with the context, i.e. fail-fast, if either is absent -- the
same contract the CCF metrics use for ``ccf_resolution``):

* ``image_soma_xyz``         -- (x, y, z) image-based soma centroid
* ``image_soma_radius_xyz``  -- (rx, ry, rz) soma radius per axis (all > 0)

Optionally, if ``image_zarr_path`` (and ``soma_mip_path``) are provided, a soma
MIP QC image is rendered via ``standard_morph.imaging.render_soma_mip`` and
recorded in ``result.artifacts``. A rendering failure is recorded in
``measurements["soma_mip_error"]`` but never fails the metric -- the image is a
bonus, the offset is the check.
"""
import time

import numpy as np

from standard_morph.metrics.base import Metric, Applicability
from standard_morph.models.qc_context import Space, ALL_MORPHOLOGY_KINDS
from standard_morph.models.qc_result import MetricResult
from standard_morph.registry import register

_AXES = ("x", "y", "z")


class SomaAtCentroidMetric(Metric):
    name = "soma_at_centroid"
    display_name = "Soma at image centroid"

    # Image-based reference, so only meaningful in image space. Reads only the
    # soma's coordinate + compartment, so it survives an untrustworthy topology.
    requires_topology = False

    applicability = Applicability(
        spaces=frozenset({Space.IMAGE_SPACE}),
        morphology_kinds=ALL_MORPHOLOGY_KINDS,
        required_resources=frozenset({"image_soma_xyz", "image_soma_radius_xyz"}),
    )
    required_policy_keys = frozenset({"max_offset_fraction"})

    def evaluate(self, prepared_morph, context, policy):
        t0 = time.perf_counter()
        result = MetricResult(name=self.name, status="pass")
        max_offset = policy[self.name, "max_offset_fraction"]
        result.thresholds_used = {"max_offset_fraction": max_offset}

        # -- Validate the SWC soma: exactly one canonical soma root. --
        soma_idx = prepared_morph.soma_roots
        if soma_idx.size != 1:
            return self._error(
                result, t0,
                f"Expected exactly one soma (type 1, root); found {soma_idx.size}.",
                {"n_soma": int(soma_idx.size)},
            )

        # -- Validate the user-provided reference centroid + radii. --
        ref = np.asarray(context.resources["image_soma_xyz"], dtype=float)
        radius = np.asarray(context.resources["image_soma_radius_xyz"], dtype=float)
        if ref.shape != (3,) or radius.shape != (3,):
            return self._error(
                result, t0,
                "image_soma_xyz and image_soma_radius_xyz must each have 3 values (x, y, z).",
                {},
            )
        if np.any(radius <= 0):
            return self._error(
                result, t0,
                f"Soma radius must be positive on every axis; got {radius.tolist()}.",
                {},
            )

        # -- Offset as a fraction of the soma radius, per axis. --
        i = int(soma_idx[0])
        swc_soma = prepared_morph.xyz[i].astype(float)
        per_axis = np.abs(swc_soma - ref) / radius
        max_frac = float(np.max(per_axis))

        result.value = max_frac
        result.value_label = "max_soma_offset_fraction"
        result.measurements = {
            "swc_soma_xyz": [float(v) for v in swc_soma],
            "image_soma_xyz": [float(v) for v in ref],
            "image_soma_radius_xyz": [float(v) for v in radius],
            "per_axis_offset_fraction": {a: float(f) for a, f in zip(_AXES, per_axis)},
            "max_offset_fraction": max_frac,
            "mean_offset_fraction": float(np.mean(per_axis)),
            "ellipsoid_offset_fraction": float(np.sqrt(np.sum(per_axis ** 2))),
        }

        if max_frac > max_offset:
            over = ", ".join(a for a, f in zip(_AXES, per_axis) if f > max_offset)
            result.status = self.violation_severity.value
            result.flagged_node_ids = [int(prepared_morph.node_id[i])]
            result.flagged_node_coordinates = [tuple(float(v) for v in swc_soma)]
            result.message = (
                f"SWC soma is {max_frac:.2f} soma-radii from the image centroid on its "
                f"worst axis ({over}), exceeding {max_offset:.2f}."
            )
        else:
            result.message = (
                f"SWC soma is within {max_offset:.2f} soma-radii of the image centroid on "
                f"every axis (worst {max_frac:.2f})."
            )

        self._maybe_render_mip(context, ref, swc_soma, result)

        result.runtime_ms = (time.perf_counter() - t0) * 1000
        return result

    @staticmethod
    def _error(result, t0, message, measurements):
        result.status = "error"
        result.message = message
        result.measurements = measurements
        result.counts = measurements
        result.runtime_ms = (time.perf_counter() - t0) * 1000
        return result

    @staticmethod
    def _maybe_render_mip(context, image_soma_xyz, swc_soma_xyz, result):
        """Render the soma MIP if requested; never raise (bonus artifact only)."""
        zarr_path = context.resources.get("image_zarr_path")
        if not zarr_path:
            return
        mip_path = context.resources.get("soma_mip_path")
        if not mip_path:
            result.measurements["soma_mip_error"] = (
                "image_zarr_path given but soma_mip_path missing; no MIP rendered."
            )
            return
        try:
            from standard_morph.imaging import render_soma_mip

            out = render_soma_mip(
                zarr_path=zarr_path,
                output_path=mip_path,
                image_soma_xyz=image_soma_xyz,
                swc_soma_xyz=swc_soma_xyz,
                crop_size=context.resources.get("soma_mip_crop_size", 128),
                mip_depth=context.resources.get("soma_mip_depth", 10),
            )
            result.artifacts.append({
                "type": "soma_mip",
                "path": out,
                "description": "Soma MIP with image centroid (green) and SWC soma (red) marked.",
            })
        except Exception as exc:  # noqa: BLE001 -- imaging must never fail the metric
            result.measurements["soma_mip_error"] = f"{type(exc).__name__}: {exc}"


register(SomaAtCentroidMetric())
