"""The top bar and two side legs of the filter_v6 capacitor conductor."""

from dataclasses import dataclass
from math import isfinite

import gdstk


@dataclass(frozen=True)
class CapParams:
    """Grouped capacitor dimensions in um; w and l are outer dimensions."""

    w: float = 302
    l: float = 600
    t_top: float = 14
    t_side: float = 6

    def __post_init__(self):
        for name in ("w", "l", "t_top", "t_side"):
            value = getattr(self, name)
            if not isfinite(value) or value <= 0:
                raise ValueError(f"cap.{name} must be finite and positive")
        if 2 * self.t_side >= self.w:
            raise ValueError("cap.w must exceed 2 * cap.t_side")
        if self.t_top >= self.l:
            raise ValueError("cap.l must exceed cap.t_top")


def make_cap(cap: CapParams, *, name="CAP", layer=1, datatype=0):
    """Return a reusable cell with one continuous Pi-shaped metal polygon.

    Origin: lower-left outer corner. +x is right; +y is up.
    Placement belongs to the parent cell, using a gdstk.Reference.
    """
    w, l = cap.w, cap.l
    side, top = cap.t_side, cap.t_top
    points = [
        (0, 0), (side, 0), (side, l - top),
        (w - side, l - top), (w - side, 0), (w, 0),
        (w, l), (0, l),
    ]
    cell = gdstk.Cell(name)
    cell.add(gdstk.Polygon(points, layer=layer, datatype=datatype))
    return cell
