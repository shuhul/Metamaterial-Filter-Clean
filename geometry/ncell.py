"""Parameterized array, two end cells, feeds, and planar ground geometry."""

from dataclasses import dataclass, field
from math import isfinite

import gdstk

if __package__:
    from .unit_cell import UnitCellParams, make_unit_cell
else:
    from unit_cell import UnitCellParams, make_unit_cell


@dataclass(frozen=True)
class NCellParams:
    unit: UnitCellParams = field(default_factory=UnitCellParams)
    n_cells: int = 4  # Interior resonators; two end cells are added.
    d_cell: float = 60  # Outer cap edge to outer cap edge.
    t_ground: float = 4
    d_end_sep: float = 14
    w_side: float = 412
    l_top_ground: float = 448
    l_sim: float = 1500
    t_feed: float = 10
    d_feed: float = 6  # Feed underside to side-ground top.

    def __post_init__(self):
        if isinstance(self.n_cells, bool) or not isinstance(self.n_cells, int) or self.n_cells < 1:
            raise ValueError("n_cells must be a positive integer")
        for name in ("d_cell", "t_ground", "d_end_sep", "w_side", "l_top_ground",
                     "l_sim", "t_feed", "d_feed"):
            if not isfinite(getattr(self, name)) or getattr(self, name) <= 0:
                raise ValueError(f"{name} must be finite and positive")
        if self.t_ground >= self.d_cell:
            raise ValueError("t_ground must be smaller than d_cell")
        if self.l_top_ground + self.l_total >= self.l_sim:
            raise ValueError("l_sim must leave space for the bottom ground plane")
        if self.t_feed > self.unit.cap.t_top:
            raise ValueError("t_feed must fit within the cap top bar")
        if self.unit.d_top <= 0:
            raise ValueError("unit.d_top must separate caps and feeds from top ground")
        if self.unit.cap.l >= self.unit.cap.t_top + self.unit.meander.l:
            raise ValueError("Cap legs must end above bottom ground")
        if self.unit.d_top + self.t_feed + self.d_feed >= self.l_total:
            raise ValueError("Side ground must extend above the bottom ground")

    @property
    def total_cells(self):
        return self.n_cells + 2

    @property
    def pitch(self):
        return self.unit.cap.w + self.d_cell

    @property
    def l_total(self):
        """Ground-strip length from unit origin to meander/pad bottom."""
        return self.unit.d_top + self.unit.cap.t_top + self.unit.meander.l

    @property
    def w_sim(self):
        return (2 * (self.w_side + self.d_end_sep) + self.total_cells * self.unit.cap.w
                + (self.total_cells - 1) * self.d_cell)

    @property
    def cell_centers(self):
        first = self.w_side + self.d_end_sep + self.unit.cap.w / 2
        return tuple(first + i * self.pitch for i in range(self.total_cells))

    @property
    def box(self):
        return ((0, self.l_top_ground - self.l_sim), (self.w_sim, self.l_top_ground))


def make_ncell(params: NCellParams, *, name="NCELL", layer=1, datatype=0):
    """Return a hierarchical layout. y=0 is the bottom edge of top ground.

    Every unit origin is at (cell_center, 0). GDS contains planar metal only;
    the intended simulation box is available as params.box.
    """
    p = params
    unit = make_unit_cell(p.unit, name=f"{name}_UNIT", layer=layer, datatype=datatype)
    cell = gdstk.Cell(name)
    for x in p.cell_centers:
        cell.add(gdstk.Reference(unit, origin=(x, 0)))

    ground = gdstk.Cell(f"{name}_GROUND")
    rects = []

    def rect(x0, y0, x1, y1):
        rects.append(gdstk.rectangle((x0, y0), (x1, y1)))

    bottom = p.box[0][1]
    side_top = -p.unit.d_top - p.t_feed - p.d_feed
    rect(0, 0, p.w_sim, p.l_top_ground)
    rect(0, bottom, p.w_sim, -p.l_total)
    rect(0, bottom, p.w_side, side_top)
    rect(p.w_sim - p.w_side, bottom, p.w_sim, side_top)
    for left, right in zip(p.cell_centers, p.cell_centers[1:]):
        x = (left + right) / 2
        rect(x - p.t_ground / 2, -p.l_total, x + p.t_ground / 2, 0)
    ground.add(*gdstk.boolean(rects, [], "or", precision=.001, layer=layer, datatype=datatype))
    cell.add(gdstk.Reference(ground))

    feeds = gdstk.Cell(f"{name}_FEEDS")
    top = -p.unit.d_top
    for x0, x1 in ((0, p.cell_centers[0]), (p.cell_centers[-1], p.w_sim)):
        feeds.add(gdstk.rectangle((x0, top-p.t_feed), (x1, top), layer=layer, datatype=datatype))
    cell.add(gdstk.Reference(feeds))
    return cell
