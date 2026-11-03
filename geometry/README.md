For the editable, no-arguments runner, open [run_filter.py](run_filter.py) and
see [RUN_FILTER.md](RUN_FILTER.md). All geometry and stage flags are in its settings block.

Read the [complete workflow and handoff](WORKFLOW_HANDOFF.md) and
[current state JSON](workflow_state.json) for the selected design and iteration procedure.

# Parameterized geometry

The [workflow audit](../WORKFLOW_AUDIT.md) records validation and limitations.
Sonnet preparation now rejects stale GDS/metadata pairs and dimensions that
cannot safely use the fixed 2 um mesh. Dimension parameters must be multiples
of 2 um; centered lead placement can still shift 1 um during grid snapping.
Keep at least 2 um clearance on either side of ground strips and between pads
and meanders. The existing 400 um reference planes require `w_side > 400`.

## First Sonnet import test

For the generated single-cell ideal-capacitor workflow and inductance sweep, see
[generated single-cell instructions](../sonnet_files/generated_ncell_1/README.md).
Both `generate_ncell.py` and `prepare_sonnet.py` accept `--n-cells` (default 4).

`python geometry/prepare_sonnet.py` (using the project `.venv`) prepares the
four-cell GUI-review project from the existing `ncell_4.gds` and JSON metadata.
It creates a separate GDS snapped to a 2 um grid, translates with Sonnet's
`gds.exe`, maps the superconducting metal, attaches two ports, and validates the
result without running the EM solver. See
[the review instructions](../sonnet_files/generated_ncell_4/README.md).

## Full n-cell layout

Run `uv run --locked python -m geometry.generate_ncell` or
`.\.venv\Scripts\python.exe geometry\generate_ncell.py` from the project root.
Edit `NCellParams` in `geometry/generate_ncell.py`. Outputs are
`geometry/output/ncell_4.gds`, `.png`, and `.json` (the number follows `n_cells`).

| Parameter | Default | Meaning (dimensions in um) |
| --- | --- | --- |
| `n_cells` | 4 | Interior resonators; two identical end cells are added. |
| `unit` | `UnitCellParams()` | Grouped cap, meander, pad and unit-placement dimensions. |
| `d_cell` | 60 | Clear outer-cap-edge spacing between every neighboring pair. |
| `t_ground` | 4 | Width of a ground strip centered in each inter-cell gap. |
| `d_end_sep` | 14 | Horizontal gap from outer end-cap edge to side ground. |
| `w_side` | 412 | Width of each side-ground rectangle. |
| `l_top_ground` | 448 | Top-ground rectangle height. |
| `l_sim` | 1500 | Intended simulation height. |
| `t_feed` | 10 | Feedline thickness; aligned to the cap top edge. |
| `d_feed` | 6 | Vertical gap from feed underside to side-ground top. |

The starting `d_cell=60`, `t_ground=4`, `t_feed=10`, and `d_feed=6` follow
the source layout. User dimensions take priority over its small asymmetries.
For custom components, pass `unit=UnitCellParams(...)` to `NCellParams`.

Derived dimensions:

```text
total_cells = n_cells + 2
pitch = unit.cap.w + d_cell
l_total = unit.d_top + unit.cap.t_top + unit.meander.l
w_sim = 2*(w_side + d_end_sep) + total_cells*unit.cap.w + (total_cells-1)*d_cell
```

`l_total` is the ground-strip length, not a renamed meander parameter. It is
618 um = 8 + 14 + 596 by default: the meander begins below the cap top bar.
There are `total_cells - 1` ground strips. They connect the top and bottom
ground planes. The bottom ground begins at the common meander/pad bottom,
connecting all pads and bottom meander leads. Cap legs remain 10 um above it.
The centered ground strips leave 28 um clearance to each neighboring cap.

Coordinate convention: x=0 is the left box edge; y=0 is the top-ground bottom
edge, where the unit-cell origins lie. The intended box is
`(0, l_top_ground-l_sim)` to `(w_sim, l_top_ground)`. For N=4 it is
`(0, -1052)` to `(2964, 448)`, giving 2964 x 1500 um. The first cap begins
at x=426 um; centers are 577, 939, 1301, 1663, 2025, 2387 um.
Negative y coordinates are intentional.

Feedlines extend from the box edges to the two end-cell centers and overlap
their cap bars. Side grounds end 6 um below the feeds. The top and bottom
planes span the full width; all ground rectangles are merged before export.
The cells remain hierarchical references, including the two end cells.

The GDS is planar metal on layer 1 / datatype 0 with 1 nm precision. It does
**not** configure a Sonnet simulation box, materials, ports, reference planes,
or analysis. The JSON records intended box coordinates and dimensions for a
later import step; the PNG shows the complete geometry for review.

Reusable components, assembly, and export scripts live here. `generate_unit_cell.py`
exports a cap, a meander, and two pads based on the central cell of
`sonnet_files/filter_v6.sonx`. Standalone cap and meander generators remain available.

## Naming and units

All dimensions are in **microns**. Local coordinates use +x to the right and
+y upward. Dimension names describe the component before rotation or placement.

| Name | Meaning |
| --- | --- |
| `w` | Width: horizontal extent (x). |
| `l` | Length: vertical extent (y). |
| `d` | Edge-to-edge distance between structures, unless explicitly documented otherwise. |
| `t` | In-plane thickness of a metal feature, measured across that feature. |

Use suffixes to identify a feature, such as `t_top` and `t_side`. Here `t` does
not mean deposited film thickness in z; that belongs to the simulation stackup.
For the meander, `d` is the clear gap between neighboring horizontal runs.

Group parameters by component instead of encoding the component in every name:

```python
from geometry import CapParams, make_cap

cap = CapParams(w=302, l=600, t_top=14, t_side=6)
cell = make_cap(cap)
print(cap.t_top)
```

Thus the requested `w_cap` and `l_cap` are `cap.w` and `cap.l`; the top
thickness is `cap.t_top`. `CapParams` is an immutable dataclass. For a variant,
use `dataclasses.replace(cap, t_top=16)` and build a new cell with a unique name.

## Capacitor definition

The conductor has one top bar joining two equal side legs, with an open bottom.
The lower-left outer corner is (0, 0).

| Parameter | Default (um) | Measurement |
| --- | --- | --- |
| `cap.w` | 302 | Full outer width, including both side legs. |
| `cap.l` | 600 | Full outer length, including the top bar. |
| `cap.t_top` | 14 | Vertical thickness of the top bar. |
| `cap.t_side` | 6 | Horizontal thickness of each side leg. |

The inner opening is `cap.w - 2 * cap.t_side` = 290 um wide and
`cap.l - cap.t_top` = 586 um tall. Changing thickness preserves the outer
dimensions. The top bar occupies y = 586 to 600 um.

The source Sonnet file uses downward y coordinates and separate polygons for
the bar and legs (central polygon IDs 22, 23, 24). Their outer bounds are about
302 x 600 um, with a top thickness about 15 um. This generator uses the
requested **14 um** thickness and joins the three features into one polygon.

This is the capacitor conductor only, not a complete electrical capacitor or
resonator: surrounding conductors, gaps, ground, inductor, and ports come later.

## Generate and view

From the project root:

```powershell
uv run --locked python -m geometry.generate_cap
```

Or use the existing local environment directly:

```powershell
.\.venv\Scripts\python.exe -m geometry.generate_cap
```

Direct script execution also works, including an editor's Run button when
the project `.venv` interpreter is selected:

```powershell
.\.venv\Scripts\python.exe geometry\generate_cap.py
```

Outputs are saved in `geometry/output/`, beside the geometry code:

- `geometry/output/single_cap.gds`: TOP references one CAP cell, containing
  one polygon on layer 1, datatype 0. Units are 1 um; precision is 1 nm.
- `geometry/output/single_cap.png`: annotated view of the generated polygon.
  Labels are only in the preview.

`make_cap` builds local geometry; the parent assembles components using
`gdstk.Reference` for translation and rotation. Components do not write files.
The export script owns library settings, placement, and output paths. This keeps
the same component usable in a single cell or a repeated filter assembly.

## Meander

`MeanderParams` groups the dimensions of a square-turn meander:

```python
from geometry import MeanderParams, make_meander

meander = MeanderParams(w=190, l_top=10, l=596, n=72, t=4, d=4)
cell = make_meander(meander)
print(meander.l_bot)  # 14 um, derived
```

| Parameter | Default | Definition |
| --- | --- | --- |
| `w` | 190 um | Full outer width of the meander body. |
| `l_top` | 10 um | Top lead length above the outer edge of the first run. |
| `l` | 596 um | Full outer vertical extent, including both leads. |
| `n` | 72 | Total horizontal runs, including the first and last partial runs. |
| `t` | 4 um | Trace thickness, also used for turns and leads. |
| `d` | 4 um | Clear edge-to-edge vertical gap between adjacent runs. |

There are `n - 1` square U-turns. There is no bend radius. Derived dimensions:

```text
pitch = t + d
l_body = n * t + (n - 1) * d
l_bot = l - l_top - l_body
```

The defaults give a 572 um body and 14 um bottom lead. Increasing `n` by one
reduces `l_bot` by 8 um while preserving the full height. Invalid combinations
that leave `l_bot <= 0` are rejected. At least three runs are required to span
both sides of the specified width; odd and even counts are supported.

Local +y is upward. Both leads are centered at x = w/2. From the top lead,
the first run heads right, alternating across the body, with the final run
returning to the bottom lead. The body topology, 72 runs, and 4 um gaps follow
central polygon 305 in `filter_v6.sonx`. The generated leads are exactly centered;
the source's separate connector polygons are not copied.

Run from the project root (direct execution and module execution both work):

```powershell
.\.venv\Scripts\python.exe geometry\generate_meander.py
# Or:
uv run --locked python -m geometry.generate_meander
```

Edit independent parameters in `geometry/generate_meander.py`. Outputs are
`geometry/output/single_meander.gds` and `geometry/output/single_meander.png`.
The GDS has TOP referencing MEANDER, on layer 1 / datatype 0, with 1 um units
and 1 nm precision. It contains only the meander, not an assembled capacitor
and inductor. GDS export may split the continuous conductor into touching
polygons to comply with the writer's vertex limit.

## Unit cell and pads

Edit the grouped parameters in `geometry/generate_unit_cell.py`:

```python
params = UnitCellParams(
    cap=CapParams(w=302, l=600, t_top=14, t_side=6),
    meander=MeanderParams(w=190, l=596, l_top=10, n=72, t=4, d=4),
    pad=PadParams(w=24, l=586),
    d_top=8,
    d_pad=12,
)
```

All components use `w` for outer width and `l` for outer length. The former
meander `l_total` parameter is now `meander.l`. Both pads reuse one rectangular
PAD cell, with dimensions `params.pad.w` and `params.pad.l`.

The unit-cell origin is at the cap's horizontal center, `d_top` above its
**outer top edge**, not above the center of its bounding box. +y is upward,
so the geometry lies below the origin. The centered meander's top lead touches
the underside of the cap top bar. Pads sit inside the cap legs, with a clear
horizontal `d_pad` gap to each adjacent leg's inner edge. Their bottom edges
align with the meander bottom. With the defaults:

| Component | x bounds (um) | y bounds (um) |
| --- | --- | --- |
| Cap | -151 to 151 | -608 to -8 |
| Meander | -95 to 95 | -618 to -22 |
| Left pad | -133 to -109 | -618 to -32 |
| Right pad | 109 to 133 | -618 to -32 |

The pads are 12 um from the cap legs and 14 um from the meander's outer sides.
The cap and meander form one connected conductor; the two pads remain separate.
Ground, ports, stackup, and simulation configuration are not added here.

```powershell
.\.venv\Scripts\python.exe geometry\generate_unit_cell.py
# Or:
uv run --locked python -m geometry.generate_unit_cell
```

Outputs: `geometry/output/unit_cell.gds` and `geometry/output/unit_cell.png`.
TOP references UNIT_CELL, which references its cap, meander, and shared pad cell.
All metal is on layer 1 / datatype 0, with 1 um units and 1 nm precision.

## Component origins and placement

For a standalone pad GDS and dimensioned PNG, run:

```powershell
.\.venv\Scripts\python.exe geometry\generate_pad.py
# Or:
uv run --locked python -m geometry.generate_pad
```

Edit `PadParams(w=24, l=586)` in `geometry/generate_pad.py`. Outputs are
`geometry/output/single_pad.gds` and `geometry/output/single_pad.png`.
The standalone generator and unit-cell generator each specify their own
parameters; editing a standalone example does not change the unit cell.

Primitive cells keep their lower-left bounding-box corner at (0, 0). Use `place`
to position any chosen anchor at a parent coordinate, without modifying shared
geometry:

```python
from geometry import place

cap_ref = place(cap_cell, origin=(0, -8), anchor="top")
pad_ref = place(pad_cell, origin=(-133, -618), anchor="bottom_left")
# An arbitrary component-local point can be used as the anchor as well:
custom_ref = place(cap_cell, origin=(0, 0), anchor=(151, 608))
```

Named anchors: `center`, `top`, `bottom`, `left`, `right`, `top_left`,
`top_right`, `bottom_left`, `bottom_right`. Center means bounding-box center.
`origin` is where the selected anchor lands, in microns. For a complete unit
cell, `gdstk.Reference(unit_cell, origin=(x, y))` places its designed origin
(the point above the cap) at (x, y) without reanchoring its bounding box.

## Four-cell check with tuned meanders

Run `.\.venv\Scripts\python.exe geometry\test_ncell_186.py` from the project root.
This creates a separate array with all six physical meanders at 186 um, uses
the existing 2 um mesh and 4.5-5.5 GHz/500-interval sweep, and fits C0/Cg with
the measured L0 from `output_data/meander_tuning_ref400/best_refined/results.json`.
Projects go to `sonnet_files/generated_ncell_4_w186/`, GDS/preview/metadata to
`geometry/output/ncell_4_w186.*`, and fit results to
`output_data/generated_ncell_4_w186/`. Existing experiment folders are not
overwritten. The earlier 190 um results and shared parameter defaults stay intact.

The same runner accepts `--d-cell` for isolated spacing trials, for example
`.\.venv\Scripts\python.exe geometry\test_ncell_186.py --d-cell 58`.
Nondefault spacing adds `_d58` (or the requested value) to output folders.
The meander span and saved inductance remain fixed; fit JSON includes the Cg
target error and a `Cg_within_2_percent` result.

To repeat the selected design in fresh files, run:

```powershell
.\.venv\Scripts\python.exe -u geometry\test_ncell_186.py --d-cell 56 --run-name manual_01
```

Use a new run name for each repeat. This runs the full array workflow using the
saved selected inductance; it does not rerun the single-cell inductance sweep.

## First taper experiments

See [taper/README.md](taper/README.md) for the modular mixed taper/regular layout,
measured trial records, and separate strong/weak coupling targets. Targets are
in [taper/first_taper_targets.json](taper/first_taper_targets.json). Original
regular and test projects are preserved.
