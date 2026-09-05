"""SWC file reading for the QC framework.

Standalone loader so the QC layer does not depend on the legacy ``Standardizer``.
Handles the Janelia/Horta ``# OFFSET X Y Z`` header by folding the offset into
the coordinates, and returns a DataFrame with the canonical SWC columns.

Loading is deliberately *lenient*: malformed values are coerced to ``NaN`` rather
than crashing, and an empty file is returned as an empty DataFrame rather than
raising. Reporting those problems is the job of the input-integrity metrics
(``non_empty``, ``castable_columns``), so a bad file produces a clean QC report
instead of an exception before the run even starts.
"""
import pandas as pd

SWC_COLUMN_NAMES = ["node_id", "compartment", "x", "y", "z", "r", "parent"]


def _read_offset(path):
    """Return the (x, y, z) offset from a ``# OFFSET`` header, or zeros."""
    with open(path, "r") as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith("# OFFSET"):
                parts = stripped.split()
                # Expected: "# OFFSET X Y Z" -> 5 tokens.
                if len(parts) != 5:
                    raise ValueError(
                        f"Invalid OFFSET header in {path!r}; expected exactly 3 "
                        f"numbers after '# OFFSET', got: {stripped!r}"
                    )
                try:
                    return tuple(float(v) for v in parts[2:5])
                except ValueError:
                    raise ValueError(
                        f"Invalid OFFSET header in {path!r}; expected numeric "
                        f"values after '# OFFSET', got: {stripped!r}"
                    )
            if not stripped.startswith("#") and stripped:
                break  # reached data before any OFFSET line
    return (0.0, 0.0, 0.0)


def read_swc(path):
    """Load an SWC file into a canonical DataFrame.

    Parameters
    ----------
    path : str
        Path to the SWC file. Comment lines (``#``) are ignored; a Horta
        ``# OFFSET`` header, if present, is applied to the coordinates.

    Returns
    -------
    pandas.DataFrame
        Columns ``node_id, compartment, x, y, z, r, parent``. Every column is
        coerced to numeric -- values that are not numeric become ``NaN`` (the
        ``castable_columns`` integrity metric reports them). The result may be
        empty (the ``non_empty`` metric reports that); this function does not
        raise on malformed content, only on an unreadable file or a malformed
        ``# OFFSET`` header.
    """
    offset_x, offset_y, offset_z = _read_offset(path)

    swc_df = pd.read_csv(
        path,
        sep=r"\s+",
        comment="#",
        header=None,
        names=SWC_COLUMN_NAMES,
    )

    # Coerce to numeric so a malformed value becomes NaN instead of crashing the
    # read (or, for coordinates, poisoning the OFFSET arithmetic below). We do
    # NOT hard-cast to int here -- that would turn a NaN into a garbage integer;
    # the integrity phase validates types and the morphology build casts once
    # the values are known to be clean.
    for col in SWC_COLUMN_NAMES:
        swc_df[col] = pd.to_numeric(swc_df[col], errors="coerce")

    swc_df["x"] = swc_df["x"] + offset_x
    swc_df["y"] = swc_df["y"] + offset_y
    swc_df["z"] = swc_df["z"] + offset_z
    return swc_df
