"""Assemble a capacitor, centered meander, and two identical pads."""

from dataclasses import dataclass, field
from math import isfinite

import gdstk

if __package__:
    from .capacitor import CapParams, make_cap
    from .meander import MeanderParams, make_meander
    from .pad import PadParams, make_pad
    from .placement import place
else:
    from capacitor import CapParams, make_cap
    from meander import MeanderParams, make_meander
    from pad import PadParams, make_pad
    from placement import place


@dataclass(frozen=True)
class UnitCellParams:
    cap: CapParams = field(default_factory=CapParams)
    meander: MeanderParams = field(default_factory=MeanderParams)
    pad: PadParams = field(default_factory=PadParams)
    d_top: float = 8
    d_pad: float = 12

    def __post_init__(self):
        if not isfinite(self.d_top) or self.d_top < 0:
            raise ValueError("d_top must be finite and nonnegative")
        if not isfinite(self.d_pad) or self.d_pad <= 0:
            raise ValueError("d_pad must be finite and positive")
        inner_half = self.cap.w / 2 - self.cap.t_side
        if inner_half - self.d_pad - self.pad.w <= self.meander.w / 2:
            raise ValueError("Pads must fit between cap legs and meander with clearance")
        if self.pad.l >= self.meander.l:
            raise ValueError("Pads must end below the underside of the cap top bar")


def make_unit_cell(params: UnitCellParams, *, name="UNIT_CELL", layer=1, datatype=0):
    """Return a hierarchical cell, origin centered d_top above the cap top.

    Cap and meander meet at the underside of the top bar. Pads are separated
    from the inner cap-leg edges by d_pad and aligned to the meander bottom.
    """
    p = params
    cap = make_cap(p.cap, name=f"{name}_CAP", layer=layer, datatype=datatype)
    meander = make_meander(p.meander, name=f"{name}_MEANDER", layer=layer, datatype=datatype)
    pad = make_pad(p.pad, name=f"{name}_PAD", layer=layer, datatype=datatype)
    cap_top = -p.d_top
    meander_top = cap_top - p.cap.t_top
    bottom = meander_top - p.meander.l
    pad_outer = p.cap.w / 2 - p.cap.t_side - p.d_pad
    cell = gdstk.Cell(name)
    cell.add(
        place(cap, origin=(0, cap_top), anchor="top"),
        place(meander, origin=(0, meander_top), anchor="top"),
        place(pad, origin=(-pad_outer, bottom), anchor="bottom_left"),
        place(pad, origin=(pad_outer, bottom), anchor="bottom_right"),
    )
    return cell
