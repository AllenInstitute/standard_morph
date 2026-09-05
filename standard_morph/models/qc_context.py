"""QC context object and shared context enums.

The context carries the runtime metadata a metric needs to decide whether it
is applicable and how to evaluate. The user provides `space` explicitly; the
pipeline stage (pre-/post-registration) is derived from it, never entered
directly.
"""
from dataclasses import dataclass, field
from enum import Enum


class Space(str, Enum):
    """Coordinate space the morphology lives in."""

    IMAGE_SPACE = "image_space"
    CCF_REGISTERED = "ccf_registered"


class MorphologyKind(str, Enum):
    """What the SWC represents."""

    AXON = "axon"
    DENDRITE = "dendrite"
    MERGED = "merged"


#: Convenience set for metrics that apply to any morphology kind.
ALL_MORPHOLOGY_KINDS = frozenset(MorphologyKind)
ALL_COORDINATE_SPACES = frozenset(Space)


@dataclass
class QCContext:
    """Runtime metadata required to evaluate metric applicability.

    Parameters
    ----------
    space : Space
        Coordinate space, ``image_space`` or ``ccf_registered`` (required).
    morphology_kind : MorphologyKind
        Whether this is an axon, dendrite, or merged reconstruction.
    resources : dict
        Optional external inputs keyed by name, e.g. ``{"image_path": ...}``
        or ``{"ccf_atlas_path": ...}``. Metrics declare which keys they need.
    ccf_resolution : int
        Microns per voxel used to convert micron coordinates to atlas voxel
        indices. Defaults to 10 (the bundled Allen CCF atlas). Only needs to
        be set explicitly when supplying a custom atlas via
        ``resources["ccf_atlas_path"]`` or ``resources["ccf_annotation"]``
        with a different resolution. Only meaningful when
        ``space == ccf_registered``.
    policy_version : str
        Selected threshold policy version, e.g. ``"policy_v1"``.
    """

    space: Space
    morphology_kind: MorphologyKind = MorphologyKind.MERGED
    resources: dict = field(default_factory=dict)
    ccf_resolution: int = 10
    policy_version: str = "policy_v1"

    def __post_init__(self):
        # Coerce plain strings into enums so callers can pass either.
        self.space = Space(self.space)
        self.morphology_kind = MorphologyKind(self.morphology_kind)
        if self.resources is None:
            self.resources = {}
