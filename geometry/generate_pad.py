"""Export one rectangular pad and a labeled preview. Run directly or with -m."""

from pathlib import Path

import gdstk
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

if __package__:
    from .pad import PadParams, make_pad
else:
    from pad import PadParams, make_pad


def save_preview(cell, pad, path):
    fig, ax = plt.subplots(figsize=(5, 9))
    ax.add_patch(Polygon(cell.polygons[0].points, facecolor="seagreen"))
    ax.annotate("", xy=(0, pad.l + 20), xytext=(pad.w, pad.l + 20),
                arrowprops={"arrowstyle": "<->"})
    ax.text(pad.w / 2, pad.l + 35, f"pad.w = {pad.w:g} um", ha="center")
    ax.annotate("", xy=(-25, 0), xytext=(-25, pad.l),
                arrowprops={"arrowstyle": "<->"})
    ax.text(-40, pad.l / 2, f"pad.l = {pad.l:g} um", rotation=90,
            ha="center", va="center")
    ax.set(xlim=(-100, pad.w + 100), ylim=(-25, pad.l + 65),
           xlabel="x (um)", ylabel="y (um)", title="Single pad")
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main():
    pad = PadParams(w=24, l=586)
    cell = make_pad(pad)
    library = gdstk.Library(unit=1e-6, precision=1e-9)
    library.add(cell)
    top = library.new_cell("TOP")
    top.add(gdstk.Reference(cell))
    output_dir = Path(__file__).resolve().parent / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    gds_path = output_dir / "single_pad.gds"
    library.write_gds(str(gds_path))
    preview_path = output_dir / "single_pad.png"
    save_preview(cell, pad, preview_path)
    print(f"Generated: {gds_path}")
    print(f"Preview: {preview_path}")
    print(f"Parameters (um): {pad}")
    print(f"Bounding box (um): {top.bounding_box()}")
    print("TOP -> PAD; layer 1, datatype 0; units 1 um, precision 1 nm")


if __name__ == "__main__":
    main()
