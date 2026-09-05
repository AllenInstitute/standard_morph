"""Optional image-rendering helpers for QC artifacts.

These pull in heavy, optional dependencies (``zarr``, ``s3fs``, ``imageio``,
``scikit-image``) *lazily*, so the core QC framework never depends on them.
Callers must wrap invocations in ``try/except``: a rendering failure means a
missing QC *bonus*, never a reason to fail a metric.

Adapted from the legacy ``tools.get_soma_mip``. Differences: it supports a local
OME-Zarr path as well as an anonymous ``s3://`` URL, it centres the crop on the
image-based soma centroid, and it marks both the image centroid and the SWC soma
so a reviewer can see the offset the metric measured.

NOTE: this has not been exercised against a live OME-Zarr in the test suite
(the optional deps are not installed there); validate against a real image
before relying on the rendered output.
"""
import math


def _open_ome_zarr(zarr_path):
    """Open an OME-Zarr root group from a local path or an ``s3://`` URL."""
    import zarr

    path = str(zarr_path)
    if path.startswith("s3://"):
        import s3fs

        s3 = s3fs.S3FileSystem(anon=True, client_kwargs={"region_name": "us-west-2"})
        store = s3fs.S3Map(root=path[len("s3://"):], s3=s3, check=False)
        return zarr.open(store, mode="r")
    return zarr.open(path, mode="r")


def _draw_cross(rgb, py, px, color, arm=4):
    """Draw a small crosshair at (py, px) in an (H, W, 3) uint8 array."""
    h, w = rgb.shape[:2]
    if not (0 <= py < h and 0 <= px < w):
        return  # marker falls outside the crop; skip silently
    for d in range(-arm, arm + 1):
        if 0 <= px + d < w:
            rgb[py, px + d] = color
        if 0 <= py + d < h:
            rgb[py + d, px] = color


def render_soma_mip(zarr_path, output_path, image_soma_xyz, swc_soma_xyz,
                    crop_size=128, mip_depth=10):
    """Render a soma-centred MIP PNG with the centroid and SWC soma marked.

    The crop is centred on ``image_soma_xyz`` (the reference centroid). The image
    centroid is marked green and the SWC soma red. Returns ``output_path``.

    Raises if the optional dependencies are missing or the image cannot be read
    or cropped -- the caller is expected to catch and record that, not propagate.
    """
    import imageio
    import numpy as np
    from skimage.exposure import rescale_intensity

    group = _open_ome_zarr(zarr_path)
    arr = group["0"]
    # OME-Zarr scale metadata: microns/voxel per axis, ordered [t, c, z, y, x].
    scale = (
        group.attrs.asdict()["multiscales"][0]["datasets"][0]
        ["coordinateTransformations"][0]["scale"]
    )

    def to_voxel(xyz):
        x, y, z = xyz
        return (
            int(round(z / scale[2])),
            int(round(y / scale[3])),
            int(round(x / scale[4])),
        )

    cz, cy, cx = to_voxel(image_soma_xyz)
    s = math.ceil(crop_size / 2)
    y0, y1 = cy - s, cy + s
    x0, x1 = cx - s, cx + s

    crop = arr[0, 0, cz - mip_depth:cz + mip_depth, y0:y1, x0:x1].max(axis=0)
    img = rescale_intensity(np.asarray(crop), out_range=(0, 255)).astype("uint8")

    # Grayscale -> RGB so the markers can be coloured.
    rgb = np.stack([img, img, img], axis=-1)
    for xyz, color in ((image_soma_xyz, (0, 255, 0)), (swc_soma_xyz, (255, 0, 0))):
        _, vy, vx = to_voxel(xyz)
        _draw_cross(rgb, vy - y0, vx - x0, color)

    imageio.imwrite(output_path, rgb)
    return output_path
