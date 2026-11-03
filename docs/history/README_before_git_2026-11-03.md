# Purcell filter parameter extraction

For an editable script with all geometry and simulation flags, use
[geometry/run_filter.py](geometry/run_filter.py). See the [run guide](geometry/RUN_FILTER.md).

Start with the [complete workflow and handoff](geometry/WORKFLOW_HANDOFF.md) and
[machine-readable current state](geometry/workflow_state.json). The selected design
uses 186 um meanders and d_cell=56 um; generic defaults still describe older geometry.

Fit circuit models to Sonnet two-port S-parameters, compare them with the
**4.8?5.2 GHz** target band, and save S21 plots. The component targets are
**L0 = 4.9194 nH, C0 = 190.4 fF, Cg = 8.2865 fF**.

## Setup with uv

Install [uv](https://docs.astral.sh/uv/getting-started/installation/) if needed,
then run from this directory:

```powershell
uv sync --locked
uv run --locked extract_circuit_params.py --no-simulate
uv run --locked extract_circuit_params_ncell.py --no-simulate
```

`pyproject.toml` declares dependencies; `uv.lock` records the resolved versions.
uv manages `.venv/`, so manual activation is unnecessary. `.python-version` selects
Python 3.14; the project declares Python 3.11 or later. See the
[uv project guide](https://docs.astral.sh/uv/guides/projects/) for dependency management.

The single-cell script currently has `RUN_SIMULATION = True` and
`RUN_INDUCTANCE = True`. Use `--no-simulate` for saved-data fitting. To run Sonnet,
check `RUN_SIMULATION` in the appropriate script and `EM_EXE` in `constants.py`,
and omit `--no-simulate`.
A working Sonnet installation/license is required only for simulations.

## Parameterized GDS geometry

Run `uv run --locked python -m geometry.generate_ncell` to generate the modular
n-cell layout, two end cells, feeds, ground planes, and preview in `geometry/output/`. See
[geometry/README.md](geometry/README.md) for grouped parameters, dimension
conventions, and the component API.

## Taper exploration

The active ground-placement correction and **0.025 fixed-target S21 RMS**
criterion are in the [V3 campaign](sonnet_files/taper_exploration/three_tapers_v3_2026-10-26/README.md),
run by `geometry/taper/three_tapers_v3.py`. The regular junction strip is centered;
there is no strip in the taper input gap. Earlier V2 acceptance used 0.05.
The selected T1 now passes the revised limit with target RMS **0.01888848**;
its complete case is in
`sonnet_files/taper_exploration/final_designs/new_T1_plain_cw344_w220_p532_g56_j66_c590_tg6_final/`.

The active replacement specification is T1 = 5.97 nH / 148 fF / 15.3 fF,
T2 = 4.33 nH / 149 fF / 153 fF, and T3 = 2.77 nH / 251 fF / 261 fF.
See the [corrected T1 campaign](sonnet_files/taper_exploration/three_tapers_v3_2026-10-26/README.md)
and [active taper state](geometry/taper/three_tapers_status.json).
It uses one input electrode, with no input-end T1 resonator. Earlier selected
taper designs below are historical and are not the active T1 geometry.

[geometry/specifications.py](geometry/specifications.py) centralizes the unit-cell,
taper-1, taper-2 and taper-3 targets, units, grid and tolerance. See the
[taper-2 exploration notes](sonnet_files/taper_exploration/TAPER_2_README.md)
for geometry, calibration, trial history and model limitations. Editable scripts
are `geometry/taper/run_taper.py` and `geometry/taper/run_second_taper.py`.
Accepted native projects and their shareable plots belong under
`sonnet_files/taper_exploration/final_designs/`; exploratory trial folders are
not automatically accepted designs.

## Where to change settings

| Location | Settings |
| --- | --- |
| `geometry/specifications.py` | Unit-cell and both taper electrical targets, units, tolerance and geometry grid. |
| `constants.py` | Unit-target aliases imported from the specification, band targets, derived target Cg/C0 ratio, paths, `EM_EXE`, and 50-ohm reference impedance. |
| `extract_circuit_params.py` | Single-cell filename (`filter_v6`), derived paths, simulation/inductance flags, C_add sweep, small C_add for capacitance extraction, fit window, and evaluation limit. |
| `extract_circuit_params_ncell.py` | `N_CELLS` (default 4), derived filenames/paths, simulation flag, fit margin, and evaluation limit. Included cell counts: 2, 3, 4. |
| `circuit_common.py` | Shared fitting bounds/seeds, peak thresholds, dB clipping, and plot settings. |

The original single-cell workflow uses `TARGET_L0_NH` when not extracting inductance.
The multi-cell CLI also defaults to that target, but accepts explicit `--l0-nh`. The
multi-cell fit uses `TARGET_C0_FF` and `TARGET_CG_FF` directly as its first seed
and plot references; there are no redundant fixed-L0 or reference aliases.
A fitted result can differ from its target. An extracted single-cell L0 is **not**
automatically written back to `constants.py`; keep targets intact and pass the measured
value with `--l0-nh` when fitting subsequent multi-cell data.

Paths are anchored to the project directory. Units are explicit in setting names:
Hz, nH, fF, or pF. Model calculations convert to SI units.

## Single-cell workflow and C_add

The included single-cell inputs are `sonnet_files/filter_v6.sonx` and
`sonnet_files/filter_v6.snp`. The export still identifies its original project
name as `filter_new_v6_test2` in its header; it has been provided under the new
filename. The project export rule is `$BASENAME.snp`, so future runs of
`filter_v6.sonx` should produce `filter_v6.snp`.

For an inductance sweep, set both `RUN_SIMULATION = True` and
`RUN_INDUCTANCE = True` in `extract_circuit_params.py`, then run:

```powershell
uv run --locked extract_circuit_params.py
```

C_add is an ideal capacitor intended to sit **in parallel with the resonator's
C0**, adding a known amount to its capacitance. For an ideal isolated parallel LC:

```text
C_total = C0 + C_add
omega_res = 1 / sqrt(L0 * (C0 + C_add))
1 / omega_res^2 = L0 * C_add + L0 * C0
```

A line fitted to `1/(2*pi*f_res)^2` versus C_add (in farads) therefore has slope
L0 (in henries). Constant effective capacitance/loading contributes to the
intercept. A sufficiently large C_add can dominate the intrinsic capacitance,
but **domination is not required for the ideal linear relation**. The configured
0.02?0.05 pF sweep is 20?50 fF, compared with a C0 target of 190.4 fF, so this is
a capacitive perturbation rather than a dominating capacitance.

For the actual symmetric lossless `Cg - shunt LC - Cg` model with equal
50-ohm terminations, setting the reflection coefficient to zero gives the
transmission-peak condition (derived from `circuit_models.py`):

```text
1 / omega_peak^2 = L0 * [C0 + C_add + 2*Cg / (1 + (omega_peak*Cg*Z0)^2)]
```

Here `omega*Cg*Z0` is about 0.013 near 5 GHz, so the coupling contribution is
approximately the constant `2*Cg`. Thus C0 and coupling mainly affect the
intercept; C_add need not dominate either one for the slope method. What matters
is that L0 and the effective loading stay stable, and the resonance shift is
resolved. This calculation assumes the dummy/feed structures are represented by
the model's ideal coupling capacitors.

Using the target values, the model predicts 4.763, 4.661, 4.566, and 4.476 GHz
for 20, 30, 40, and 50 fF. Its continuous-frequency slope differs from the target
L0 by only about 0.001% due to this coupling correction. These are analytical
predictions, not a new Sonnet sweep or an error bound on the physical device.
Adding 500 fF predicts about 2.699 GHz, outside the project's currently enabled
4-6 GHz adaptive sweep. Larger steps therefore require changing the frequency
coverage and may probe a different effective inductance. The current code takes
the largest sampled S21 value, rejects endpoint or ambiguous peaks, and checks
that the resonance decreases as C_add increases before fitting the slope.

The implementation edits the project's `C_add` variable at each point, runs
Sonnet, takes the maximum sampled S21 magnitude as the resonance, and performs
linear regression. The actual transmission peak includes coupling/loading;
frequency-dependent loading, changing modes, and insufficient frequency
resolution can bias this approximation. Check the linearity and R?, confirm the
same resonance is tracked, and make sure the Sonnet frequency sweep covers its
shift. Increasing C_add can move it outside the original frequency range.

After the sweep, the script sets C_add to **0.0001 pF = 0.1 fF**, reruns Sonnet,
and fits complex S21 (real and imaginary parts) to determine C0 and Cg. The model
neglects that small added capacitance. With `RUN_INDUCTANCE = False`, it skips the
sweep and uses the target L0. It writes `inductance_fit.png` only when sweeping.

`--no-simulate` overrides simulation and skips the sweep, project edits, and
Sonnet execution. It uses the target L0 and whatever data are already in the
export. Files must correspond to the intended geometry and small C_add.

Simulation **edits the single-cell project in place**. The original C_add value
is not restored; a successful full run leaves it at 0.0001 pF. Interrupted runs
can leave it at a sweep value. All C_add-specific handling lives in the single-cell script; generic project
variable editing lives in `circuit_common.py`.
The multi-cell workflow neither changes C_add nor extracts inductance.

## Multi-cell topology: one capacitor between cells

The implemented ABCD cascade is:

```text
M = M_Cg * (M_LC * M_Cg)^N

N = 2:
input -- Cg -- node 1 -- Cg -- node 2 -- Cg -- output
                 |              |
              L0 || C0       L0 || C0
                 |              |
               ground         ground
```

There are **N+1 coupling capacitors and N resonators**: exactly one Cg between
neighboring resonators. The code starts with one input Cg, then appends one shunt
LC and one Cg for each cell. It does not concatenate complete
`Cg | LC | Cg` blocks, which would double the internal series capacitors.
The topology has been inspected and left unchanged.

In the supplied N=2 screenshot, the two middle cells are the resonators; the two
outer dummy structures provide the terminal coupling structures, as described
by the designer. They are not counted as two extra resonators in this model.
The screenshot supports that intended interpretation, but cannot by itself
establish that the dummy structures are electrically equivalent to ideal Cg
capacitors over the whole band.

Multi-cell fitting uses S21 magnitude in dB, clipped to ?80?0 dB. It tries target,
peak-based, and generic Cg seeds and retains the best residual. The fit window
extends 0.5 GHz beyond the outermost peaks, falling back to the dominant peak
when fewer than two are detected. Single-cell fitting retains its complex-S21
objective and ?1.5 GHz resonance window.

Both workflows use the **same plot function**: colors, ?60?5 dB y-axis, 0.5 GHz
margin beyond the fit window, target-band shading, and C0/Cg target comparisons.
Their displayed frequency ranges can differ because their fit windows differ.
Plots are saved to `output_data/<project>/s21_fit.png`; reruns replace the image.

## Interpreting the model

The reported band estimates are:

```text
f_high = 1 / (2*pi*sqrt(L0*C0))
f_low  = 1 / (2*pi*sqrt(L0*(C0 + 4*Cg)))
```

These describe the ideal uniform-chain band limits, not measured ?3 dB edges of
a finite filter. The fit assumes identical lossless resonators and equal nearest-
neighbor coupling capacitors, including the end couplings. It does not include
longer-range mutual capacitance, mutual inductance, or separate dummy/feed models.

For the current regular unit-cell layouts through N=4, negligible additional
parasitic coupling is the designer's working assumption. This applies to these
geometries, not as a general consequence of the cell count. The taper is a separate,
unfinished design and is not covered by that assumption. Its unequal cells and
couplings will need a nonuniform model. The current multi-cell objective fits
magnitude only; agreement does not independently establish phase agreement.

## Future geometry automation ? investigation only

Modular cap, meander, pad, unit-cell, and n-cell generation is implemented in
`geometry/`. Command-line Sonnet import, grid checks, port attachment, and a
generated single-cell ideal-capacitor extraction workflow are working. See
[geometry instructions](geometry/README.md) and [the workflow audit](WORKFLOW_AUDIT.md).
Automated target optimization has not started.

Sonnet 19.52 supports importing into the current project and replacing its
polygons by clearing **Keep polygons in present project**. A project template
can carry stackup, import settings, frequencies, and EM options. GDS port-marker
rectangles on a mapped port layer can recreate numbered ports using cells such
as `SonPort_Box_1`. See the
[19.52 Import Wizard documentation](https://www.sonnetsoftware.com/support/doc_19.52/Topics/Interfaces_Translators/Import_Wizard.htm).

For this layout, I recommend generating geometry and explicit port markers from
one parameter set, then translating into a **copy of a configured Sonnet template**.
Keep coordinate conventions, mesh, stackup, and reference planes controlled.
The installed `gds.exe -h` confirms project-input and geometry append options;
Sonnet also documents its [command-line translator and macro runner](https://www.sonnetsoftware.com/support/doc_19.52/Topics/Appendix/Command_Line_Programs.htm).
An actual import needs validation before relying on it in a sweep.

An alternative is to keep the port-bearing feed polygons fixed and replace only
the interior geometry in a copied `.sonx`. Inspection of these local XML files
shows ports nested inside polygons with vertex references. Consequently, blindly
replacing the entire metal layer can delete the ports or invalidate their edge
attachments. Selective XML editing is a possible custom approach, not a verified
replacement API. Polygon IDs, parameter references, grounding, and connections
would need to remain consistent.

Increasing N also changes the footprint. The current right-hand box port changes
position between projects; preserving every coordinate and box dimension may
therefore be inappropriate. Either use a sufficiently large fixed box with
controlled feed extensions, or update the box/port positions and de-embedding
consistently. For a fixed N, the existing named geometry variables are a simpler
first automation step than regenerating every polygon. Exporting GDS for fabrication
can remain a separate output of the same geometry definition.

## File guide

| File/path | Purpose |
| --- | --- |
| `README.md` | Setup, settings, model assumptions, workflows, and automation investigation. |
| `pyproject.toml` | uv project metadata and Python dependency declarations. |
| `uv.lock` | uv-managed dependency lockfile; keep with the project. |
| `.python-version` | Default Python version for uv. |
| `.gitignore` | Excludes local environments, caches, and bytecode from version control. |
| `constants.py` | Only genuinely shared configuration. |
| `extract_circuit_params.py` | Single-cell configuration, C_add/L0 sweep, complex capacitance fit, and entry point. |
| `extract_circuit_params_ncell.py` | Multi-cell configuration, fixed-L0 dB capacitance fit, and entry point. |
| `circuit_models.py` | Single-cell S21, resonance seed, and N-cell ABCD equations. |
| `circuit_common.py` | Generic Sonnet project-variable editing, execution, export reading, fitting helpers, reports, and plots. |
| `Screenshot 2026-10-01 002323.png` | N=2 geometry reference with two interior cells and two outer dummy structures. |
| `sonnet_files/filter_v6.sonx` / `.snp` | Single-cell Sonnet project and saved two-port export. |
| `sonnet_files/filter_v6_cells_n=2.sonx` / `.snp` | Two-cell project and export. |
| `sonnet_files/filter_v6_cells_n=3.sonx` / `.snp` | Three-cell project and export. |
| `sonnet_files/filter_v6_cells_n=4.sonx` / `.snp` | Four-cell project and export. |
| `output_data/filter_v6/s21_fit.png` | Single-cell comparison plot. |
| `output_data/filter_v6_cells_n=2/s21_fit.png` | Two-cell comparison plot. |
| `output_data/filter_v6_cells_n=3/s21_fit.png` | Three-cell comparison plot. |
| `output_data/filter_v6_cells_n=4/s21_fit.png` | Four-cell comparison plot. |
| `output_data/filter_v6/inductance_fit.png` | Created by the optional inductance sweep. |

`sonnet_files/sondata/` contains solver-managed folders for `filter_v6`,
`filter_new_v6_test2` (the original single-cell name), and N=2, 3, 4. Not every
folder contains every file below. Extraction reads `.snp` exports, not these caches.

| Solver file | Purpose |
| --- | --- |
| `info.log` | Run summary metadata. |
| `log_composite.log` | Combined solver log. |
| `log_errors.log` | Diagnostics/errors and run context. |
| `log_response.log` | Response/progress log. |
| `log_timing.log` | Timing information. |
| `geo_response_d.sid`, `response_d.sid` | Sonnet de-embedded response cache files. |
| `geo_response_nd.sid`, `response_nd.sid` | Non-de-embedded response cache files. |
| `geo_interpdata_d.sid`, `interpdata_d.sid` | Interpolated response cache files. |

`.venv/` is uv's generated environment; `__pycache__/` contains Python bytecode.
The scripts print fitted parameters to the console rather than saving CSV/JSON.
If fitting fails, check input paths, frequency coverage, L0, peak selection, and
bounds. If simulation fails, check `EM_EXE`, the Sonnet license, and export settings.
