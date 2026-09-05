"""CCF annotation atlas access for post-registration QC metrics.

The atlas is the Allen CCF annotation volume: a 3-D ``uint`` array where each
voxel holds a structure id and ``0`` means "outside the brain". A CCF-registered
SWC stores coordinates in microns, so a node's voxel index is
``floor(coord / resolution)`` and the array is indexed directly as
``annotation[x, y, z]`` (no axis swap). A node is inside the brain when its voxel
is in bounds and non-zero.

Conventions and coordinate handling mirror ``morph_utils.ccf``.

The real 10um volume decompresses to ~4.8 GB, so it is loaded once and cached.
Metrics may also be handed a preloaded array to avoid any disk load (used in
tests and to reuse an atlas across many cells).
"""
import numpy as np

# Cache keyed by (resolution, path-or-"bundled") -> annotation array.
_ATLAS_CACHE = {}


def _bundled_atlas_path(resolution):
    from importlib.resources import files

    return files("standard_morph") / "data" / f"annotation_{resolution}.nrrd"


def load_ccf_annotation(resolution=10, atlas_path=None):
    """Load (and cache) the CCF annotation volume as a 3-D array.

    Parameters
    ----------
    resolution : int
        Microns per voxel; selects the bundled ``annotation_{resolution}.nrrd``
        when ``atlas_path`` is not given.
    atlas_path : str, optional
        Explicit path to a ``.nrrd`` atlas, overriding the bundled file.
    """
    key = (resolution, str(atlas_path) if atlas_path else "bundled")
    if key in _ATLAS_CACHE:
        return _ATLAS_CACHE[key]

    try:
        import nrrd
    except ImportError as e:  # pragma: no cover - exercised only without the extra
        raise ImportError(
            "Reading a CCF atlas requires the 'nrrd' package. "
            "Install the ccf extra:  pip install '.[ccf]'"
        ) from e

    path = atlas_path if atlas_path is not None else _bundled_atlas_path(resolution)
    array, _ = nrrd.read(str(path))
    _ATLAS_CACHE[key] = array
    return array


def clear_atlas_cache():
    """Drop all cached atlas arrays (frees the ~4.8 GB 10um volume from RAM)."""
    _ATLAS_CACHE.clear()


def coordinates_to_voxels(coords, resolution):
    """Convert an ``(n, 3)`` array of micron coordinates to integer voxels."""
    return np.floor(np.asarray(coords, dtype=float) / resolution).astype(int)


def in_brain_mask(voxels, annotation):
    """Boolean ``(n,)`` mask: True where a voxel is in bounds and non-zero.

    Out-of-bounds voxels (including negative indices) are treated as outside the
    brain rather than wrapping around the array.
    """
    voxels = np.asarray(voxels)
    n = len(voxels)
    in_bounds = np.ones(n, dtype=bool)
    shape = annotation.shape
    for d in range(3):
        in_bounds &= (voxels[:, d] >= 0) & (voxels[:, d] < shape[d])

    result = np.zeros(n, dtype=bool)
    if in_bounds.any():
        v = voxels[in_bounds]
        result[in_bounds] = annotation[v[:, 0], v[:, 1], v[:, 2]] != 0
    return result
