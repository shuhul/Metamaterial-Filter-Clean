"""Generate the assembled unit cell and its annotated preview."""

from pathlib import Path

import gdstk
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

if __package__:
    from .unit_cell import UnitCellParams, make_unit_cell, CapParams, MeanderParams, PadParams
else:
    from unit_cell import UnitCellParams, make_unit_cell, CapParams, MeanderParams, PadParams


def save_preview(cell, params, path):
    fig, ax = plt.subplots(figsize=(8, 10))
    for ref, color, label in zip(cell.references,
                                 ("steelblue", "darkorange", "seagreen", "seagreen"),
                                 ("Cap", "Meander", "Pads", None)):
        for i, polygon in enumerate(ref.get_polygons()):
            ax.add_patch(Polygon(polygon.points, facecolor=color,
                                 label=label if i == 0 else None))
    ax.plot(0, 0, "+", color="black", markersize=12)
    ax.annotate(f"Origin (0, 0)\nd_top = {params.d_top:g} um", xy=(0, 0),
                xytext=(35, 28), arrowprops={"arrowstyle": "->"}, fontsize=9)
    x_inner = -params.cap.w / 2 + params.cap.t_side
    gap_y = -params.d_top - params.cap.l / 2
    ax.annotate("", xy=(x_inner, gap_y), xytext=(x_inner + params.d_pad, gap_y),
                arrowprops={"arrowstyle": "<->"})
    ax.annotate(f"d_pad = {params.d_pad:g} um\n(both sides)",
                xy=(x_inner + params.d_pad/2, gap_y), xytext=(-270, gap_y + 65),
                arrowprops={"arrowstyle": "->"}, fontsize=9)
    bottom = -params.d_top - params.cap.t_top - params.meander.l
    ax.axhline(bottom, color="gray", linestyle="--", linewidth=.7)
    ax.text(0, bottom - 20, "Pads and meander share bottom edge", ha="center", fontsize=9)
    ax.set(xlim=(-280, 245), ylim=(bottom - 50, 70), xlabel="x (um)", ylabel="y (um)",
           title="Unit cell: cap + meander + two pads")
    ax.set_aspect("equal")
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main():
    params = UnitCellParams(
        cap=CapParams(w=302, l=600, t_top=14, t_side=6),
        meander=MeanderParams(w=190, l=596, l_top=10, n=72, t=4, d=4),
        pad=PadParams(w=24, l=586),
        d_top=8,
        d_pad=12,
    )
    cell = make_unit_cell(params)
    library = gdstk.Library(unit=1e-6, precision=1e-9)
    library.add(cell, *cell.dependencies(True))
    top = library.new_cell("TOP")
    top.add(gdstk.Reference(cell))
    output_dir = Path(__file__).resolve().parent / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    gds_path = output_dir / "unit_cell.gds"
    library.write_gds(str(gds_path))
    preview_path = output_dir / "unit_cell.png"
    save_preview(cell, params, preview_path)
    print(f"Generated: {gds_path}")
    print(f"Preview: {preview_path}")
    print(f"Parameters (um except n): {params}")
    print(f"Meander derived l_bot: {params.meander.l_bot:g} um")
    print(f"Origin: cap horizontal center, {params.d_top:g} um above cap top")
    for label, ref in zip(("Cap", "Meander", "Left pad", "Right pad"), cell.references):
        print(f"{label} bounds (um): {ref.bounding_box()}")
    print("Layer 1, datatype 0; units 1 um; precision 1 nm")


if __name__ == "__main__":
    main()
