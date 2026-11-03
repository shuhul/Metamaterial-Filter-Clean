"""Generate one capacitor conductor and an annotated preview."""

from pathlib import Path

import gdstk
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

# Support both direct execution (including an editor's Run button) and -m.
if __package__:
    from .capacitor import CapParams, make_cap
else:
    from capacitor import CapParams, make_cap


def save_preview(cell, cap, path):
    """Draw the actual generated polygon; annotations are not GDS geometry."""
    fig, ax = plt.subplots(figsize=(6, 8))
    ax.add_patch(Polygon(cell.polygons[0].points, facecolor="steelblue"))

    def dimension(start, end, label, text_position):
        ax.annotate("", xy=start, xytext=end,
                    arrowprops={"arrowstyle": "<->", "color": "black"})
        ax.text(*text_position, label, ha="center", va="center", fontsize=10)

    dimension((0, cap.l + 30), (cap.w, cap.l + 30),
              f"cap.w = {cap.w:g} um", (cap.w / 2, cap.l + 48))
    dimension((-30, 0), (-30, cap.l),
              f"cap.l = {cap.l:g} um", (-100, cap.l / 2))
    ax.annotate(f"cap.t_top = {cap.t_top:g} um",
                xy=(cap.w / 2, cap.l - cap.t_top / 2),
                xytext=(cap.w / 2, cap.l - 100), ha="center",
                arrowprops={"arrowstyle": "->"})
    ax.annotate(f"cap.t_side = {cap.t_side:g} um (both legs)",
                xy=(cap.t_side / 2, cap.l / 3),
                xytext=(cap.w / 2, cap.l / 3 - 70), ha="center",
                arrowprops={"arrowstyle": "->"})
    ax.set(xlim=(-180, cap.w + 50), ylim=(-35, cap.l + 75),
           xlabel="x (um)", ylabel="y (um)", title="Single capacitor conductor")
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main():
    # Edit these grouped parameters; every dimension is in microns.
    cap = CapParams(w=302, l=600, t_top=14, t_side=6)
    cap_cell = make_cap(cap)

    library = gdstk.Library(unit=1e-6, precision=1e-9)
    library.add(cap_cell)
    top = library.new_cell("TOP")
    top.add(gdstk.Reference(cap_cell))

    output_dir = Path(__file__).resolve().parent / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    gds_path = output_dir / "single_cap.gds"
    library.write_gds(str(gds_path))
    preview_path = output_dir / "single_cap.png"
    save_preview(cap_cell, cap, preview_path)

    print(f"Generated: {gds_path}")
    print(f"Preview: {preview_path}")
    print(f"Parameters (um): {cap}")
    print("TOP -> CAP: one conductor, layer 1, datatype 0")
    print(f"Bounding box (um): {top.bounding_box()}")
    print("Units: 1 um; precision: 1 nm")


if __name__ == "__main__":
    main()
