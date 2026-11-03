"""Reusable rectangular pad; dimensions in microns."""

from dataclasses import dataclass
from math import isfinite

import gdstk


@dataclass(frozen=True)
class PadParams:
    w: float = 24
    l: float = 586

    def __post_init__(self):
        for name in ("w", "l"):
            if not isfinite(getattr(self, name)) or getattr(self, name) <= 0:
                raise ValueError(f"pad.{name} must be finite and positive")


def make_pad(pad: PadParams, *, name="PAD", layer=1, datatype=0):
    cell = gdstk.Cell(name)
    cell.add(gdstk.rectangle((0, 0), (pad.w, pad.l), layer=layer, datatype=datatype))
    return cell
