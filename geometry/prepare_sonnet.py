"""Prepare an isolated generated Sonnet project; never run the solver."""

from hashlib import sha256
from pathlib import Path
import argparse
import json
import re
import subprocess
import xml.etree.ElementTree as ET

import gdstk
import numpy as np

if __package__:
    from .validation import require, validate_gds_metadata
else:
    from validation import require, validate_gds_metadata

BASE = Path(__file__).resolve().parents[1]
DEST = BASE / "sonnet_files" / "generated_ncell_4"
SOURCE = BASE / "sonnet_files" / "filter_v6_cells_n=4.sonx"
GDS = BASE / "geometry" / "output" / "ncell_4.gds"
BIN = Path(r"C:\Program Files\Sonnet Software\19.52.2025.261\bin")
NS = "https://www.sonnetsoftware.com/schema/project"
ET.register_namespace("", NS)
ET.register_namespace("xsi", "http://www.w3.org/2001/XMLSchema-instance")


def tag(name):
    return f"{{{NS}}}{name}"


def write_xml(tree, path):
    ET.indent(tree, space="    ")
    tree.write(path, encoding="utf-8", xml_declaration=True)


def points(poly):
    return np.array([(float(x), float(y)) for x, y in re.findall(
        r"\(([^,]+),([^\)]+)\)", " ".join(p.text or "" for p in poly.findall(tag("Points"))))])


def run(command, log):
    result = subprocess.run([str(x) for x in command], cwd=log.parent,
                            capture_output=True, text=True)
    log.write_text(result.stdout + result.stderr, encoding="utf-8")
    print(result.stdout + result.stderr)
    result.check_returncode()


def save_review(actual, width, height, feed_y, dest=DEST, n_cells=4):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon

    fig, ax = plt.subplots(figsize=(14, 8))
    for poly in actual:
        ax.add_patch(Polygon(poly.points, facecolor="seagreen"))
    for number, x, direction in ((1, 0, 1), (2, width, -1)):
        ax.plot(x, feed_y, "ro", markersize=4)
        ax.annotate("", xy=(x+direction*400,feed_y), xytext=(x,feed_y),
                    arrowprops={"arrowstyle": "->", "color": "black"})
        ax.text(x+direction*80,feed_y-35,f"Port {number}", ha="center")
    ax.set(xlim=(-35,width+35), ylim=(height+35,-35),
           xlabel="Sonnet x (um)", ylabel="Sonnet y (um, downward)",
           title=f"Imported {n_cells}-cell project: 2 x 2 um mesh; 400 um reference planes")
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(dest / "imported_review.png", dpi=180)
    plt.close(fig)


def main(n_cells=4, *, gds_path=None, output_dir=None, geometry_validator=None):
    dest = Path(output_dir).resolve() if output_dir is not None else BASE / "sonnet_files" / f"generated_ncell_{n_cells}"
    gds = Path(gds_path).resolve() if gds_path is not None else BASE / "geometry" / "output" / f"ncell_{n_cells}.gds"
    dest.mkdir(parents=True, exist_ok=True)
    before = sha256(SOURCE.read_bytes()).hexdigest()
    metadata = json.loads(gds.with_suffix(".json").read_text())
    width = metadata["w_sim"]
    height = metadata["parameters"]["l_sim"]
    if metadata["parameters"]["n_cells"] != n_cells:
        raise ValueError("GDS metadata cell count does not match requested cell count")
    grid = 2.0
    if width % grid or height % grid:
        raise ValueError("Box dimensions must fit the 2 um mesh")

    # Keep the design GDS intact. Snap a separate import artifact in box coordinates.
    # Ties go toward positive x/y; the 4 um center leads shift +1 um in x.
    _, original = (geometry_validator or validate_gds_metadata)(gds, metadata)
    top_ground = metadata["parameters"]["l_top_ground"]
    expected = []
    changed = 0
    lib = gdstk.Library(unit=1e-6, precision=1e-9)
    top = lib.new_cell("TOP")
    for poly in original:
        coords = poly.points.copy()
        coords[:, 1] += height - top_ground
        snapped = np.floor(coords / grid + .5) * grid
        changed += int(np.count_nonzero(np.any(snapped != coords, axis=1)))
        require((np.max(np.abs(snapped-coords)) <= 1.000001), 'Validation failed: np.max(np.abs(snapped-coords)) <= 1.000001')
        require((np.all((snapped >= 0) & (snapped <= (width, height)))), 'Validation failed: np.all((snapped >= 0) & (snapped <= (width, height)))')
        out = gdstk.Polygon(snapped, layer=1, datatype=0)
        top.add(out)
        expected.append(out)
    snapped_gds = dest / f"ncell_{n_cells}_grid2.gds"
    lib.write_gds(str(snapped_gds))
    expected_components = metadata.get('expected_metal_components', 1) if geometry_validator else 1
    require(len(gdstk.boolean(expected, [], 'or', precision=.001)) == expected_components,
            'Unexpected disconnected metal in GDS')

    template = ET.parse(SOURCE)
    metal = template.find(f".//{tag('PlanarTechLayer')}")
    require((metal.get("MaterialName") == "SUPERCONDUCTOR"), 'Validation failed: metal.get("MaterialName") == "SUPERCONDUCTOR"')
    mapping = ET.Element(tag("ImportExportMapping"), GDSStreamDataTypes="1:0")
    metal.insert(0, mapping)
    mapped_template = dest / "mapped_template.sonx"
    write_xml(template, mapped_template)
    raw = dest / "mapped_import.sonx"
    layer_file = dest / "superconductor.lay"
    layer_file.write_text('Metal-1 metal_boundary 1 0 0 "SUPERCONDUCTOR" on 0\n')
    run([BIN / "gds.exe", "-v", "-i", mapped_template, "-g", "-l", layer_file, "-p", "TOP",
         "-u", "um", "-s", f"{width},{height}", "-t", "min", "-C", "2,2",
         "-o", raw, snapped_gds], dest / "import.log")

    tree = ET.parse(raw)
    root = tree.getroot()
    geometry = root.find(tag("Geometry"))
    require((geometry.get("YDirection") == "Negative"), 'Validation failed: geometry.get("YDirection") == "Negative"')
    layers = geometry.findall(f".//{tag('PlanarTechLayer')}")
    imported = []
    for layer in layers:
        polys = layer.findall(f".//{tag('PlanarPolygon')}")
        if polys:
            require((layer.get("MaterialName") == "SUPERCONDUCTOR"), 'Validation failed: layer.get("MaterialName") == "SUPERCONDUCTOR"')
            imported.extend(polys)
    actual = [gdstk.Polygon(points(poly)) for poly in imported]
    # Sonnet's top-left coordinate convention flips GDS y about box height.
    target = [gdstk.Polygon([(x, height-y) for x,y in poly.points]) for poly in expected]
    require((not gdstk.boolean(actual, target, "xor", precision=.001)), 'Validation failed: not gdstk.boolean(actual, target, "xor", precision=.001)')
    require((all(np.allclose(points(poly)/grid, np.round(points(poly)/grid)) for poly in imported)), 'Validation failed: all(np.allclose(points(poly)/grid, np.round(points(poly)/grid)) for poly in imported)')

    size = geometry.find(f"{tag('Box')}/{tag('Size')}")
    mesh = geometry.find(f"{tag('Box')}/{tag('NumCells')}")
    require((float(size.get('X')) == width and float(size.get('Y')) == height), "Validation failed: float(size.get('X')) == width and float(size.get('Y')) == height")
    require((float(mesh.get('X')) == width/grid and float(mesh.get('Y')) == height/grid), "Validation failed: float(mesh.get('X')) == width/grid and float(mesh.get('Y')) == height/grid")

    # Add ports on the exact feed-edge segments, using the baseline XML structure.
    cap_top = top_ground + metadata["parameters"].get("d_feed_top", metadata["parameters"]["unit"]["d_top"])
    feed_bottom = cap_top + metadata["parameters"]["t_feed"]
    require((not geometry.findall(f".//{tag('Port')}")), 'Validation failed: not geometry.findall(f".//{tag(\'Port\')}")')
    for number, edge_x in ((1, 0), (2, width)):
        matches = []
        for poly in imported:
            xy = points(poly)
            for i, (a,b) in enumerate(zip(xy, np.roll(xy, -1, axis=0))):
                if (a[0] == b[0] == edge_x and
                    sorted((a[1],b[1])) == [cap_top, feed_bottom]):
                    matches.append((poly,i))
        require((len(matches) == 1), 'Validation failed: len(matches) == 1')
        poly, vertex = matches[0]
        objects = poly.find(tag("Objects"))
        if objects is None:
            objects = ET.SubElement(poly, tag("Objects"))
        port = ET.SubElement(objects, tag("Port"), MacroID=f"portGenerated{number}",
                             Number=str(number), Vertex=str(vertex))
        ET.SubElement(port, tag("Type"), Type="BOX")
        ET.SubElement(port, tag("Impedance"), Resistance="50", Reactance="0", Inductance="0", Capacitance="0")
        ET.SubElement(port, tag("Center")).text = f"({edge_x},{(cap_top+feed_bottom)/2})"

    planes = geometry.findall(f"{tag('Box')}/{tag('RefPlane')}")
    require((len(planes) == 2), 'Validation failed: len(planes) == 2')
    for plane in planes:
        plane.set("RefLength", "400")
    require((metadata["parameters"]["w_side"] > 400), 'Validation failed: metadata["parameters"]["w_side"] > 400')
    root.find(tag("OutputFiles")).set("Folder", ".")
    # Leave legacy optimization goals disabled; they point to an unrelated external file.
    for goal in root.findall(f"{tag('Optimize')}/{tag('Goal')}"):
        goal.set("On", "FALSE")
    final = dest / f"generated_ncell_{n_cells}.sonx"
    write_xml(tree, final)

    # Check simulation settings against baseline without running em.
    baseline = ET.parse(SOURCE)
    def signature(element):
        return (element.tag, element.attrib, (element.text or '').strip(),
                [signature(child) for child in element])
    for name in ("Units", "Control", "Sweeps"):
        require((signature(root.find(tag(name))) == signature(baseline.find(tag(name)))), 'Validation failed: signature(root.find(tag(name))) == signature(baseline.find(tag(name)))')
    for name in ("Dielectric", "Conductor", "DielectricMaterialModel", "BoxCover"):
        require(([signature(x) for x in geometry.findall(f'.//{tag(name)}')] == [
            signature(x) for x in baseline.findall(f'.//{tag(name)}')]), "Validation failed: [signature(x) for x in geometry.findall(f'.//{tag(name)}')] == [ signature(x) for x in baseline.findall(f'.//{tag(name)}')]")
    require((before == sha256(SOURCE.read_bytes()).hexdigest()), 'Validation failed: before == sha256(SOURCE.read_bytes()).hexdigest()')
    report = {"project": str(final), "baseline_sha256": before,
              "source_gds_sha256": sha256(gds.read_bytes()).hexdigest(),
              "source_metadata_sha256": sha256(gds.with_suffix('.json').read_bytes()).hexdigest(),
              "project_sha256": sha256(final.read_bytes()).hexdigest(),
              "box_um": [width,height], "mesh_um": [2,2], "polygons": len(imported),
              "snapped_vertex_count": changed, "max_snap_um": 1,
              "ports": 2, "reference_planes_um": [400,400],
              "material": "SUPERCONDUCTOR", "geometry_xor_empty": True,
              "settings_preserved": ["stackup", "materials", "sweeps", "control"],
              "solver_run": False}
    (dest / "validation.json").write_text(json.dumps(report, indent=2)+"\n")
    run([BIN / "soncmd.exe", "-Read", final], dest / "sonnet_read.log")
    save_review(actual, width, height, (cap_top+feed_bottom)/2, dest, n_cells)
    print(f"Ready for GUI review: {final}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-cells", type=int, default=4)
    main(parser.parse_args().n_cells)
