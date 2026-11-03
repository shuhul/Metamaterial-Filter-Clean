"""Export a single meander and a labeled preview. Run directly or with -m."""

from pathlib import Path

import gdstk
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

if __package__:
    from .meander import MeanderParams, make_meander
else:
    from meander import MeanderParams, make_meander


def save_preview(cell, m, path):
    fig, ax = plt.subplots(figsize=(7, 9))
    for polygon in cell.polygons:
        ax.add_patch(Polygon(polygon.points, facecolor="steelblue"))
    ax.annotate("", xy=(0, m.l + 20), xytext=(m.w, m.l + 20),
                arrowprops={"arrowstyle": "<->"})
    ax.text(m.w / 2, m.l + 30, f"w = {m.w:g} um", ha="center")
    for label, xy, text_y in (
        (f"l_top = {m.l_top:g} um", (m.w / 2, m.l - m.l_top / 2), m.l - 30),
        (f"l_bot = {m.l_bot:g} um (derived)", (m.w / 2, m.l_bot / 2), 40),
    ):
        ax.annotate(label, xy=xy, xytext=(m.w + 20, text_y),
                    arrowprops={"arrowstyle": "->"}, fontsize=9)
    ax.text(m.w + 20, m.l / 2,
            f"n = {m.n} horizontal runs\n{m.n - 1} U-turns\n"
            f"t = {m.t:g} um\nd = {m.d:g} um\n"
            f"l = {m.l:g} um", va="center")
    ax.set(xlim=(-25, m.w + 210), ylim=(-20, m.l + 60),
           xlabel="x (um)", ylabel="y (um)", title="Single square-turn meander")
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main():
    # Independent parameters. l_bot is calculated to preserve l.
    meander = MeanderParams(w=190, l_top=10, l=596, n=72, t=4, d=4)
    cell = make_meander(meander)
    library = gdstk.Library(unit=1e-6, precision=1e-9)
    library.add(cell)
    top = library.new_cell("TOP")
    top.add(gdstk.Reference(cell))
    output_dir = Path(__file__).resolve().parent / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    gds_path = output_dir / "single_meander.gds"
    library.write_gds(str(gds_path))
    preview_path = output_dir / "single_meander.png"
    save_preview(cell, meander, preview_path)
    print(f"Generated: {gds_path}")
    print(f"Preview: {preview_path}")
    print(f"Parameters (um except n): {meander}")
    print(f"Derived: l_body = {meander.l_body:g} um; l_bot = {meander.l_bot:g} um")
    print(f"Bounding box (um): {top.bounding_box()}")
    print("TOP -> MEANDER; layer 1, datatype 0; units 1 um, precision 1 nm")


if __name__ == "__main__":
    main()
