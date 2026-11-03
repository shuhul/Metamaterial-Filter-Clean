"""Place reusable cells by a named bounding-box anchor or a local point."""

import gdstk


def place(cell, *, origin=(0, 0), anchor="bottom_left"):
    """Return a reference with anchor located at origin in parent coordinates.

    anchor can also be an explicit (x, y) point in the component's coordinates.
    Geometry remains unchanged, so one cell can be reused in many placements.
    """
    (x0, y0), (x1, y1) = cell.bounding_box()
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    anchors = {
        "bottom_left": (x0, y0), "bottom_right": (x1, y0),
        "top_left": (x0, y1), "top_right": (x1, y1),
        "center": (cx, cy), "top": (cx, y1), "bottom": (cx, y0),
        "left": (x0, cy), "right": (x1, cy),
    }
    if isinstance(anchor, str):
        if anchor not in anchors:
            raise ValueError(f"Unknown anchor {anchor!r}; choose from {tuple(anchors)}")
        anchor = anchors[anchor]
    return gdstk.Reference(cell, origin=(origin[0] - anchor[0], origin[1] - anchor[1]))
