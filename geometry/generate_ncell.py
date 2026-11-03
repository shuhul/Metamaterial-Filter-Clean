"""Generate the n-cell layout, review image, and intended box metadata."""

from dataclasses import asdict
from pathlib import Path
import argparse
import json

import gdstk
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

if __package__:
    from .ncell import NCellParams, make_ncell
else:
    from ncell import NCellParams, make_ncell


def save_preview(cell, p, path):
    fig, ax = plt.subplots(figsize=(max(12, p.w_sim / 250), 7))
    for i, ref in enumerate(cell.references):
        if i < p.total_cells:
            for component, color in zip(ref.cell.references,
                                        ("steelblue", "darkorange", "seagreen", "seagreen")):
                for polygon in component.get_polygons():
                    points = polygon.points + ref.origin
                    ax.add_patch(Polygon(points, facecolor=color))
        else:
            for polygon in ref.get_polygons():
                ax.add_patch(Polygon(polygon.points, facecolor="lightgray" if i == p.total_cells else "steelblue"))
    for i, x in enumerate(p.cell_centers):
        label = "End L" if i == 0 else "End R" if i == p.total_cells-1 else f"Cell {i}"
        ax.text(x, 40, label, ha="center", fontsize=9)
        ax.plot(x, 0, "+", color="black", markersize=4)
    for left, right in zip(p.cell_centers, p.cell_centers[1:]):
        ax.annotate("", xy=(left+p.unit.cap.w/2, -p.l_total-45),
                    xytext=(right-p.unit.cap.w/2, -p.l_total-45),
                    arrowprops={"arrowstyle": "<->", "linewidth": .7})
    ax.text(p.w_sim/2, -p.l_total-100,
            f"d_cell = {p.d_cell:g} um   |   t_ground = {p.t_ground:g} um   |   "
            f"ground length = {p.l_total:g} um", ha="center", fontsize=9)
    ax.text(p.w_side/2, -p.l_total/2, f"Side ground\nw = {p.w_side:g}", ha="center", fontsize=9)
    ax.text(p.w_sim-p.w_side/2, -p.l_total/2,
            f"d_end_sep\n= {p.d_end_sep:g} um", ha="center", fontsize=9)
    ax.text(p.w_sim/2, p.l_top_ground/2, f"Top ground: l = {p.l_top_ground:g} um", ha="center")
    ax.set(xlim=(-25, p.w_sim+25), ylim=(p.box[0][1]-25,p.box[1][1]+25),
           xlabel="x (um)", ylabel="y (um)",
           title=f"{p.n_cells} interior cells + 2 end cells | intended box {p.w_sim:g} x {p.l_sim:g} um")
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main(n_cells=4):
    # Unit-cell dimensions are grouped in params.unit; defaults match unit_cell.py.
    params = NCellParams(n_cells=n_cells, d_cell=60, t_ground=4,
                         d_end_sep=14, w_side=412, l_top_ground=448, l_sim=1500)
    cell = make_ncell(params)
    library = gdstk.Library(unit=1e-6, precision=1e-9)
    library.add(cell, *cell.dependencies(True))
    top = library.new_cell("TOP")
    top.add(gdstk.Reference(cell))
    output = Path(__file__).resolve().parent / "output"
    output.mkdir(parents=True, exist_ok=True)
    stem = output / f"ncell_{params.n_cells}"
    library.write_gds(str(stem.with_suffix(".gds")))
    save_preview(cell, params, stem.with_suffix(".png"))
    metadata = {"parameters": asdict(params), "units": "um", "total_cells": params.total_cells,
                "w_sim": params.w_sim, "l_total": params.l_total,
                "intended_box": params.box, "cell_centers": params.cell_centers,
                "note": "Planar GDS geometry only; Sonnet box, ports, stackup and simulation are not configured."}
    stem.with_suffix(".json").write_text(json.dumps(metadata, indent=2)+"\n", encoding="utf-8")
    print(f"Generated: {stem.with_suffix('.gds')}")
    print(f"Preview: {stem.with_suffix('.png')}")
    print(f"Metadata: {stem.with_suffix('.json')}")
    print(f"Cells: {params.n_cells} interior + 2 end = {params.total_cells}")
    print(f"Intended box: {params.w_sim:g} x {params.l_sim:g} um; bounds {params.box}")
    print(f"Ground-strip length: {params.l_total:g} um; pitch: {params.pitch:g} um")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-cells", type=int, default=4)
    main(parser.parse_args().n_cells)
