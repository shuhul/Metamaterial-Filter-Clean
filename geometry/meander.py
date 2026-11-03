"""Square-turn meander with centered top and bottom leads; dimensions in um."""

from dataclasses import dataclass
from math import isfinite

import gdstk


@dataclass(frozen=True)
class MeanderParams:
    w: float = 190
    l_top: float = 10
    l: float = 596
    n: int = 72  # Total horizontal runs, including the two partial end runs.
    t: float = 4
    d: float = 4  # Clear edge-to-edge gap between horizontal runs.

    def __post_init__(self):
        for name in ("w", "l_top", "l", "t", "d"):
            if not isfinite(getattr(self, name)) or getattr(self, name) <= 0:
                raise ValueError(f"meander.{name} must be finite and positive")
        if isinstance(self.n, bool) or not isinstance(self.n, int) or self.n < 3:
            raise ValueError("meander.n must be an integer of at least 3")
        if self.w <= 3 * self.t:
            raise ValueError("meander.w must exceed 3 * meander.t")
        if self.l_bot <= 0:
            raise ValueError("Meander body and top lead must leave a positive l_bot")

    @property
    def l_body(self):
        """Outer height of the horizontal runs and U-turns."""
        return self.n * self.t + (self.n - 1) * self.d

    @property
    def l_bot(self):
        """Clear bottom lead length below the final horizontal run."""
        return self.l - self.l_top - self.l_body


def make_meander(meander: MeanderParams, *, name="MEANDER", layer=1, datatype=0):
    """Build one continuous square-turn conductor, anchored at its lower left.

    The first run goes from the centered top lead to the right edge. The last
    run returns from its U-turn to the centered bottom lead. Intermediate runs
    span the full width. l_top and l_bot exclude the horizontal run thickness.
    """
    m = meander
    parts = []

    def rect(x0, y0, x1, y1):
        parts.append(gdstk.rectangle((x0, y0), (x1, y1)))

    center = m.w / 2
    body_top = m.l - m.l_top
    for row in range(m.n):
        y_top = body_top - row * (m.t + m.d)
        x0, x1 = 0, m.w
        if row == 0:
            x0 = center - m.t / 2
        elif row == m.n - 1:
            # Odd final rows enter from the right; even rows from the left.
            if row % 2:
                x0 = center - m.t / 2
            else:
                x1 = center + m.t / 2
        rect(x0, y_top - m.t, x1, y_top)
        if row < m.n - 1:
            x0 = m.w - m.t if row % 2 == 0 else 0
            rect(x0, y_top - 2 * m.t - m.d, x0 + m.t, y_top)

    rect(center - m.t / 2, body_top - m.t, center + m.t / 2, m.l)
    rect(center - m.t / 2, 0, center + m.t / 2, m.l_bot + m.t)
    # Merge overlaps so all bars, turns, and leads form one conductor.
    polygons = gdstk.boolean(parts, [], "or", precision=0.001,
                             layer=layer, datatype=datatype)
    cell = gdstk.Cell(name)
    cell.add(*polygons)
    return cell
