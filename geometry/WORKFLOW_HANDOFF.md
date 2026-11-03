# Filter workflow: implementation, iteration record, and handoff

Updated: 2026-10-13. Read this before continuing the project.
Machine-readable companion: [workflow_state.json](workflow_state.json).
Paths are relative to the project root unless stated otherwise.

This document describes the implemented workflow and saved evidence. Completed
actions, assumptions, and future work are distinguished. No new simulations
were run to write it.

## 1. Current decision and files to use

The latest task was to bring Cg within 2% of its target, retaining acceptable C0
and the tuned meander. It is complete. Stop iterating until the user requests
another objective.

| Quantity | Target | Current result | Meaning |
| --- | --- | --- | --- |
| L0 | 4.9194 nH | 4.922723945 nH, +0.067568% | Native filter_v6 single-cell capacitor-sweep extraction |
| C0 | 190.4 fF | 188.878116386 fF, about -0.799% | Selected N=4 fit; user accepted C0 as good enough |
| Cg | 8.2865 fF | 8.241620220 fF, -0.541601% | Selected N=4 fit; meets +/-2% |
| meander.w | Tuned geometry | 186 um | Horizontal span of the runs |
| d_cell | Tuned geometry | 56 um | Gap between neighboring outer cap edges |

Selected artifacts:

- [Design GDS](output/ncell_4_w186_d56.gds), [parameters](output/ncell_4_w186_d56.json), [design preview](output/ncell_4_w186_d56.png).
- [Final Sonnet project](../sonnet_files/generated_ncell_4_w186_d56/generated_ncell_4.sonx).
- [Simulation GDS](../sonnet_files/generated_ncell_4_w186_d56/ncell_4_grid2.gds).
- [Actual simulated S-parameters](../sonnet_files/generated_ncell_4_w186_d56/generated_ncell_4.snp).
- [S21 and fitted circuit](../output_data/generated_ncell_4_w186_d56/s21_fit.png).
- [Fit and hashes](../output_data/generated_ncell_4_w186_d56/fit_parameters.json).
- [Spacing summary](../output_data/cg_tuning/summary.json), [Cg plot](../output_data/cg_tuning/cg_vs_spacing.png), [S21 comparison](../output_data/cg_tuning/s21_comparison.png).
- [Inductance evidence used](../output_data/meander_tuning_ref400/best_refined/results.json).

**Generic defaults have not been promoted to the selected design.**
MeanderParams still defaults to w=190 and NCellParams to d_cell=60.
The array-test runner explicitly uses w=186; its default spacing is still 60.
Use the selected metadata or pass --d-cell 56. A plain default generator run
otherwise recreates earlier geometry. TARGET constants remain targets, not
measured values, and were not overwritten.

## 2. What is automated and what the assistant chooses

```mermaid
flowchart LR
    P[Grouped geometry parameters] --> G[Design GDS, JSON, preview]
    G --> I[Snap to 2 um grid and import]
    I --> S[Sonnet: materials, ports, calibration]
    S --> EM[One adaptive EM analysis]
    EM --> SP[Touchstone S-parameters]
    L[Saved single-cell L0] --> F[Python circuit fit]
    SP --> F
    F --> E[Target error and acceptance]
    E -->|Outside tolerance: choose next parameter| P
    E -->|Inside tolerance| D[Save selection and stop]
```

Python generates geometry, calls native Sonnet command-line programs, reads
exports, fits the circuit, and saves results. The assistant chooses the objective,
parameter to vary, next candidate, and stopping point.

The meander tuner implements a bounded sensitivity-based search in Python.
For spacing, the assistant selected 58 then 56 um and ran the same script for
each. There is no general-purpose autonomous global optimizer yet. Sonnet's
legacy optimizer goals are disabled. SciPy optimizes the fitted C0/Cg; it does
not change the GDS.

```text
Cg_error_percent = 100 * (Cg_fitted / 8.2865 - 1)
accept if abs(Cg_error_percent) <= 2
```

Acceptance is 8.12077 to 8.45223 fF. This applies to an effective fitted circuit
parameter, not an independent capacitance measurement or the entire S21 curve.

## 3. Environment and entry points

Windows workspace: C:\Users\shuhul\Downloads\filter_new_automation_clean.
Sonnet bin: C:\Program Files\Sonnet Software\19.52.2025.261\bin.
Installed release: Sonnet 19.52.2025.261. uv manages the local .venv;
.python-version selects Python 3.14.

```powershell
# From the project root, restore declared dependencies if needed:
uv sync --locked
# Use the local interpreter:
.\.venv\Scripts\python.exe geometry\test_ncell_186.py --help
```

uv-created environments may have no pip; that is normal. Do not install
project dependencies with global pip. gdstk is declared in pyproject.toml and
uv.lock. Other dependencies are NumPy, SciPy, matplotlib, and scikit-rf.
gdsfactory was not introduced. No CAD drawing conversion or GUI-driving
simulation automation is involved.

The Windows execution sandbox previously failed to initialize with error 1223;
the approved alternate execution path worked. Sonnet may need network access
for its license. These execution issues are distinct from native solver errors.
Keep edits within the authorized workspace.

| File | Responsibility |
| --- | --- |
| geometry/capacitor.py, meander.py, pad.py | Primitive geometry and grouped dimension dataclasses |
| geometry/placement.py, unit_cell.py, ncell.py | Anchors, composition, array, grounds, feeds |
| geometry/generate_*.py | Standalone examples and GDS/preview exporters |
| geometry/validation.py | GDS/metadata identity and supported-grid checks |
| geometry/prepare_sonnet.py | Import, material mapping, ports, reference planes, native read check |
| geometry/prepare_inductance.py | Ideal capacitor on the generated N=1 layout |
| geometry/tune_meander.py | Native filter_v6 copies and bounded run-span tuning |
| extract_generated_inductance.py | Archive C_add sweeps, regress L0, optionally fit single-cell capacitances |
| geometry/test_ncell_186.py | Isolated w=186 array, optional d_cell, one EM run, fixed-L0 fit |
| extract_circuit_params_ncell.py | Fit explicit N-cell export and save S21 overlay |
| circuit_models.py, circuit_common.py, constants.py | Equations, shared I/O/fitting, targets/paths |
| audit_workflow.py, tests/test_workflow.py | Historical artifact audit and offline regression checks |

Supporting detail: [geometry conventions](README.md), [historical audit](../WORKFLOW_AUDIT.md),
[meander notes](MEANDER_TUNING.md), [speed diagnosis](../output_data/performance_audit/README.md).

## 4. Geometry definitions and selected dimensions

Dimensions are in microns. GDS uses unit=1e-6 meters and precision=1e-9 meters:
coordinates are um with 1 nm database resolution. Database precision is separate
from the 2 um simulation mesh. Metal is layer 1, datatype 0; TOP is the single
top-level cell.

w is horizontal span; l is vertical span; d is edge-to-edge gap; t is in-plane
bar thickness/trace width, not film thickness. Parameters are grouped, such as
unit.cap.t_top, unit.meander.w, unit.pad.l. Native in_width means meander.w.

| Group | Selected values (um except counts) | Construction |
| --- | --- | --- |
| cap | w=302, l=600, t_top=14, t_side=6 | Pi-shaped top bar and two legs, open bottom |
| meander | w=186, l=596, l_top=10, n=72, t=4, d=4 | Square turns, centered top/bottom leads |
| pad | w=24, l=586 | One primitive referenced twice per unit |
| unit | d_top=8, d_pad=12 | Origin above cap, pads separated from inner cap legs |
| array | n_cells=4, d_cell=56, t_ground=4 | Four interior units plus two feed-connected end units |
| surroundings | d_end_sep=14, w_side=412, l_top_ground=448, l_sim=1500 | Side/top/bottom grounds |
| feeds | t_feed=10, d_feed=6 | Existing feed width and feed-to-side-ground gap retained |

n=72 counts all horizontal runs, including partial first/last runs: 71 turns.

```text
l_body = n*t + (n-1)*d = 572
l_bot = l - l_top - l_body = 596 - 10 - 572 = 14
```

Changing l alone with n/t/d fixed changes the bottom lead, not the number of
runs. The L experiment changed horizontal w, not vertical l.

Primitive origins are the lower-left bounding-box corner. place(cell, origin,
anchor) positions a named bounding-box anchor or explicit local point at a
parent coordinate without changing the shared cell. The unit's designed origin
is at the cap's horizontal center, 8 um above its outer top. Unit-local bounds:

- Cap x=-151..151, y=-608..-8.
- Meander x=-93..93, y=-618..-22; its top lead touches the cap underside.
- Pads x=-133..-109 and 109..133, y=-618..-32, aligned to the meander bottom.

```text
total_cells = n_cells + 2 = 6
pitch = cap.w + d_cell = 358
l_total = d_top + cap.t_top + meander.l = 618
w_sim = 2*(w_side+d_end_sep) + total_cells*cap.w + (total_cells-1)*d_cell
      = 2*(412+14) + 6*302 + 5*56 = 2944
```

Unit origins have y=0 at the top ground's underside. x centers are 577, 935,
1293, 1651, 2009, 2367. The design box runs from (0,-1052) to (2944,448).
Inter-cell ground strips extend from y=0 to -618; bottom ground joins at -618.
Spacing automatically changes pitch, positions, horizontal ground dimensions,
box width, and right port position. These are derived consequences, not extra
tuned inputs.

The feed's 10 um width differs from the end cap's 14 um top-bar thickness.
The user chose to retain this difference. Pads join bottom ground in the full
array. Meanders ground cap nodes inductively: DC connectivity does not make
all RF nodes equivalent.

## 5. What N means physically and in the circuit model

N=4 means four interior resonators in the fitted circuit. Sonnet contains six
physical structures: four interior units plus two end units connected to feeds.
It solves fields/current on all six, including mutual coupling and end loading.
The reduced model has N identical shunt LC resonators and N+1 series Cg:

```text
port -- Cg -- node1 -- Cg -- node2 -- ... -- nodeN -- Cg -- port
                 |               |                 |
              L0 || C0        L0 || C0          L0 || C0
                 |               |                 |
               ground          ground            ground
```

The end units are not separate LC branches in this model. Their influence is
approximated by terminal coupling/loading. Four visible peaks do not prove the
end units have no modes. Black reference-plane arrows do not exclude physical
structures from EM analysis. Earlier conversational mode explanations were
heuristic, not a validated modal analysis.

## 6. GDS import and deliberate Sonnet preparation changes

prepare_sonnet.py uses sonnet_files/filter_v6_cells_n=4.sonx as the common
settings/material template, including for generated N=1. Originals are preserved.

1. Read GDS/JSON, rebuild geometry from parameters, require empty polygon XOR,
   correct units/layer/TOP and consistent derived dimensions. Reject unsupported
   grid dimensions and collapsed clearances.
2. Add l_sim-l_top_ground=1052 to y to put the GDS box in positive coordinates.
3. Snap with 2*floor(coord/2+0.5), ties positive, at most 1 um movement. Preserve
   the design GDS; save a separate ncell_4_grid2.gds for actual simulation.
4. Explicitly map layer 1/datatype 0 to Metal-1 and SUPERCONDUCTOR. This fixes the
   unknown-metal assignment seen in an early diagnostic import.
5. Run native gds.exe using TOP, mapping, box dimensions and -C 2,2. Sonnet uses
   downward y. Verify imported polygons against the transformed/snapped geometry
   by XOR, and verify every imported vertex lies on the 2 um grid.
6. Find the unique feed edge on each box wall and add a BOX port with 50 ohm
   real impedance. Selected port centers: (0,461), (2944,461); feed y=456..466.
   Attachment indices come from the actual imported polygon edges.
7. Set both fixed reference planes to 400 um; OutputFiles folder to '.'; disable
   legacy Optimize goals pointing to an unrelated external target file.
8. Compare units, controls, sweeps, dielectric/conductor definitions, layer
   thicknesses and covers with the template. Confirm source hash unchanged.
   Save validation.json and imported_review.png.
9. Run soncmd.exe -Read to validate native readability. This is not an EM run.

A typical native import has this structure:

```text
gds.exe -v -i mapped_template.sonx -g -l superconductor.lay -p TOP
        -u um -s 2944,1500 -t min -C 2,2 -o mapped_import.sonx ncell_4_grid2.gds
```

Use -o for output. An early attempt with two positional paths failed.

Even dimensions do not guarantee grid-aligned positions. Centered 4 um leads
at odd centers need a +1 um x shift. At d_cell=58, centered ground strips also
round and acquire unequal side clearances. At selected d_cell=56, strip centers
fit the grid; lead shifts remain. The simulation belongs to the snapped geometry.
A fabrication export must explicitly choose design versus simulated geometry.

The array runner sets the Touchstone filename explicitly to generated_ncell_4.snp
and refreshes the final-project hash after that change. Native project copies
previously retained an expanded template basename, exporting filter_v6.snp rather
than the expected trial filename. That issue was corrected and completed data
were reused instead of repeating their EM solve.

### The 600-to-400 um reference-plane change

Original projects use 600 um shifts. The native filter_v6 control emitted warnings
that the planes could extend past uniform feeds. The generated side-ground region
is 412 um wide; 400 um ends before its first discontinuity. 600 um can extend
into the end-cell region. Generated projects and successful native tuning copies
use 400 um. Original source projects were not edited.

This changes the calibration/reference of reported S-parameters. It does not
remove metal, move the geometrical box-wall ports, or suppress end-cell modes.
Geometric plausibility and zero warnings are established; no reference-plane
sensitivity study has been run. See Sonnet's
[reference-plane documentation](https://www.sonnetsoftware.com/support/doc_19.52/Topics/Analysis_Setup/De-embedding/Shifting_Reference_Planes.htm).

### Preserved physical/material settings

Stored XML levels, in order:

| Level | Material | Thickness | Relative permittivity | Loss tangent |
| --- | --- | --- | --- | --- |
| 0 | Vacuum | 1500 um | 1 | 0 |
| 1 | Silicon (intrinsic) | 500 um | 11.9 | 5e-6 |
| 2 | Air | 5000 um | 1 | 0 |

Metal-1 uses SUPERCONDUCTOR on level 0. Box Z=7000 um; top cover FREE_SPACE,
bottom cover PERFECT_CONDUCTOR; symmetry off. SUPERCONDUCTOR uses
SURFACE_IMPEDANCE with XML Ls=0.0388 and Rdc=Rrf=Xdc=0.
**Ls is in pH/square**, not the global circuit inductance unit nH. See
[Sonnet General Loss Model](https://www.sonnetsoftware.com/support/help-18/users_guide/GeneralLossModel.html).
These values were inherited, not recalibrated to actual fabrication. Extracted
L0 includes the configured material and EM geometry effects.

Physical feed impedance in logs is about 56.89 ohm, while terminations/export
normalization are 50 ohm. Normalization does not make the physical feed matched.

## 7. Inductance extraction: a separate multi-run loop

Do not confuse these geometry families:

- Generated N=1 plus copied ideal capacitor, prepared by prepare_inductance.py:
  initial w=190 extraction gave 5.015715851 nH.
- Native copies of original filter_v6.sonx, prepared/tuned by tune_meander.py:
  selected w=186 extraction gave 4.922723945 nH.

Both contain one interior unit plus two end units. C_add is an ideal Sonnet
circuit component from the interior cap node to top ground, not fabricated GDS
metal. The generated version copies the original CAP component/internal pins
and bridges y=448 to y=456 at snapped central x=940. Native copies retain the
original component location. Component terminals 3/4 are not additional external
feed ports. Generated geometry and pin contact are checked before simulation.

```text
C_add values in Sonnet project (pF): 0.02, 0.03, 0.04, 0.05
Actual added capacitance (fF):       20,   30,   40,   50
```

For each value: change C_add, snapshot the project, run em, check completion and
zero errors, verify exported C_add header, read S21, take the largest sampled
magnitude's frequency. Reject boundary maxima, ambiguous prominent peaks, or a
resonance that fails to decrease with increasing C_add. Archive each project,
SNP and log separately. Fit:

```text
x = C_add in farads
y = 1/(2*pi*f_peak)^2 in seconds squared
y = L0*x + intercept
L0 = slope in henries
```

For an isolated parallel LC, y=L0*(C0+C_add). With ideal symmetric weak coupling,
loading adds approximately constant 2*Cg to effective capacitance, mostly affecting
the intercept. C_add need not dominate C0. Distributed/frequency-dependent loading
can bias the slope. High R squared measures linearity, not absolute accuracy.

Selected refined frequencies: 4.8165, 4.7115, 4.6130, 4.5205 GHz.
L0=4.922723945 nH, R squared=0.999999973462, regression standard
error=0.000567052 nH. This standard error is not the physical uncertainty.

Changing only the ideal capacitor can reuse the unchanged EM geometry solution;
Sonnet evaluates its component ports in circuit analysis. Four capacitor values
do not generally require four complete EM geometry solves, although every
invocation still incurs processing/export overhead.
[Sonnet components](https://www.sonnetsoftware.com/support/doc_19.52/Topics/Components/Components.htm)

Early sweeps used 4-6 GHz, Target=20000, actual 0.1 MHz output. For speed, native
screening changed to 4.35-4.95 GHz, Target=500, actual 1 MHz output. Selected
refinement used Target=1500, actual 0.5 MHz output. Sonnet chooses convenient
spacing: Target does not guarantee exactly (stop-start)/Target. Read net.f.
This narrower range covers the test-cap resonances for these candidates, not
necessarily the final 0.1 fF resonance. The older default prepare_inductance.py
still creates the broad/dense sweep; do not blindly rerun it for quick tuning.

The selected native w186 project was restored to C_add=0.0001 pF=0.1 fF and
original 4-6 GHz/500 settings. **No final 0.1 fF EM response was generated for
that selected native copy.** Its last raw response is explicitly named
last_sweep_cadd_50ff.snp; refined sweeps live under best_refined. Do not use the
50 fF response in a fit assuming negligible C_add.

## 8. Meander iteration history and selection

Native dimensional editing materializes moved vertices:

```text
soncmd.exe -ParamFile parameters.txt -ReadWrite filter_v6.sonx trial.sonx
parameters.txt: in_width=186
```

The tuner verifies that only the three linked meander polygons move horizontally,
their y coordinates stay fixed, requested spans are achieved, and other polygon
coordinates, ports, component, box and material settings are preserved. All three
native physical meanders share the changed span.

| Trial | Geometry family | Span (um) | L0 (nH) | Role |
| --- | --- | --- | --- | --- |
| Initial single cell | Generated | 190 | 5.015715851 | Earlier extraction, separate family |
| Native control | filter_v6 copy, ref=400 | 190 | 5.014942955 | High-resolution control |
| Native screen | Same | 186 | 4.926571148 | Coarse output |
| Native bracket | Same | 184 | 4.875503242 | Coarse output, below target |
| Selected refinement | Same | 186 | 4.922723945 | Final L0 used subsequently |

The tuner measures a control, perturbs span by 4 um, estimates local slope from
two measurements, predicts target-crossing span, and chooses an untried 2 um-grid
candidate. After 190/186 it tried 184 to bracket the target. It selected 186 as
closest and refined it. Sampled frequency resolution explains the small change
between screening and refinement.

The 600 um reference-plane attempt was stopped after warnings and preserved in
meander_tuning/. Successful experiments are in meander_tuning_ref400/.
Other runs were stopped while switching dense output to fast screening. A stale
lock was released only after our solver was stopped. A completed point with the
wrong export basename was recovered. Interrupted attempts are not selected results.

## 9. N=4 simulation and capacitance extraction: once per candidate

For each array trial, test_ncell_186.py generates GDS/JSON/preview, imports into
a new folder, sets the export name, and calls **em once**. It then reads selected
single-cell L0 JSON. It does not run an inductance sweep or add C_add to N=4.

Controls: adaptive ABS 4.5-5.5 GHz, Target=500, DOUBLE precision, MEDIUM_SPEED,
de-embedding on, STOP_RESTART ABS cache, 2 um mesh. Sonnet used 36 threads.
Exports contain 501 samples at 2 MHz spacing; these are interpolated, not 501
independent EM solves. The selected trial required six EM frequencies; others
required seven.

read_snp temporarily copies .snp to .s2p for scikit-rf. It checks two ports,
finite values, positive increasing frequencies and 50 ohm normalization.

The fitter holds L0 fixed, varies C0/Cg, detects prominent S21 peaks and sets a
window from first/last peaks plus 0.5 GHz on each side, intersecting available
data. The selected trial uses all 501 points. It fits S21 magnitude in dB,
clipped to [-80,0], reducing the influence of very deep notches. It does not
fit N=4 phase. The ABCD model is:

```text
Zg = 1/(j*omega*Cg)
Ylc = j*omega*C0 + 1/(j*omega*L0)
M_Cg = [[1,Zg],[0,1]]
M_LC = [[1,0],[Ylc,1]]
M_total = M_Cg * (M_LC * M_Cg)^N
S21 = 2 / (A + B/50 + C*50 + D)
```

s21_n_cells accepts C0/Cg in fF and L0 in henries. The single-cell closed form
uses SI capacitances; N=1 agrees numerically with the cascade.

SciPy tries target and resonance-derived seeds, then Cg seeds 0.5/1/2/5/10/20 fF.
Both parameters are bounded 0.1..2000 fF; max evaluations 8000 per seed. Lowest
sum of squared clipped-dB residuals wins. For the selected case the target seed
converged near Cg=8.532 fF, but another seed gave 8.241620 fF with lower residual.
Multiple seeds therefore matter.

Selected magnitude-fit RMS residual: 4.687275 dB. Effective C0/Cg are useful for
tuning, but the circuit does not reproduce all distributed effects. The earlier
single-cell capacitance fit used real/imaginary S21 instead; its values are not
interchangeable with N=4 magnitude-fit values.

Printed ideal-chain band edges:

```text
f_high = 1/(2*pi*sqrt(L0*C0))
f_low = 1/(2*pi*sqrt(L0*(C0+4*Cg)))
```

Selected implied edges: 4.8161-5.2195 GHz. They are distinct from the finite
array's four transmission peaks. Meeting the Cg criterion does not certify the
4.8-5.2 GHz target band or exact peak loss.

## 10. Spacing iteration history and choice of 56 um

L0 remained 4.922723945 nH. All six meanders remained 186 um. No L sweep was
repeated for these spacing tests.

| d_cell | C0 (fF) | Cg (fF) | Cg error | EM analysis |
| --- | --- | --- | --- | --- |
| 60 um baseline | 188.465916 | 7.940202 | -4.179% | 7 frequencies, 12 min 18 sec during slowdown |
| 58 um trial | 188.719106 | 8.083948 | -2.444% | 7 frequencies, 5 min 28 sec |
| 56 um selected | 188.878116 | 8.241620 | -0.542% | 6 frequencies, 4 min 43 sec |

Closer cells were expected to increase Cg, but sensitivity was measured:

```text
slope = (8.083948-7.940202)/(58-60) = about -0.071873 fF/um
predicted exact-target spacing = 58 + (8.2865-8.083948)/slope = about 55.18 um
```

56 um is a supported grid choice predicted to meet 2%; its actual fit did.
C0 changed slightly. A physical spacing change cannot promise that only one
circuit parameter changes. The user accepted C0 and requested Cg tuning, so cap
dimensions were held fixed. The loop stopped after two new spacing trials.

Selected external solver wall time: about 288 seconds; Sonnet's reported analysis:
283 seconds. Preparation: about 3 seconds; final read/fit/save stage: about 12
seconds. Timer boundaries differ; use stage_timings.json with simulation.log.

## 11. Speed diagnosis: established facts versus hypotheses

Earlier generated w190 N=4 took 5 min 56 sec; w186/d60 took 12 min 18 sec with
identical settings, seven frequencies and 12,489 subsections. Matrix-fill totals
increased 200 to 346 sec; matrix-solve totals 81 to 232 sec; ABS totals 42 to 75
sec. Most extra time was per-frequency EM computation, not Python, hidden L
sweeps, additional frequencies, or denser output.

An estimate-only em -N check found original GUI geometry had 12,491 subsections
versus generated 12,489. Imported polygon fragmentation did not materially
increase the solver problem size. Short diagnostics observed background Photos
CPU activity, a separate GUI solver during the slow run, 18 physical cores/36
logical processors, Balanced power plan, about 76 GB free RAM and negligible
paging/disk queue. These do not prove which historical process or scheduling
condition caused the slowdown. No throttling or 18-vs-36-thread benchmark was
run; global power/thread settings were not changed.

After the user closed background processes, spacing runs returned to normal
speed. That is consistent with contention, not a controlled causal experiment.

The unchanged-project cache-only test took 0.805 sec and zero new EM frequencies.
Its original export bytes/timestamp were restored. This shows how cached GUI
analysis can seem much faster than fresh geometry analysis; earlier manual GUI
cache state is unknown.

Efficient future policy:

- Refit saved data without Sonnet when only fitting assumptions change.
- Preserve the correct sondata cache. Never transplant caches between different
  geometries or attribute an old result to a new candidate.
- Use narrow/coarse screening with adequate mode coverage; refine only the
  selected case when the acceptance margin requires it.
- Use measured sensitivity for the next grid-valid candidate and stop at tolerance.
- Avoid repeating L sweeps for each small spacing trial. Re-extract when meander,
  trace/gap, material, grounding or substantial loading changes invalidate transfer.
- Avoid overlapping heavy jobs. Treat 18 physical-core threads as an untested
  benchmark candidate, not an established faster setting.
- Save stage timers and before/after CPU counters. Their failure should not block
  a valid simulation; the runner now makes these snapshots best-effort.

Sonnet's -AbsNoDiscrete forbids new discrete solves with existing ABS cache.
Insufficient cache may produce no usable export: verify completion and data.
-N estimates subsections without EM solves. See
[CLI documentation](https://www.sonnetsoftware.com/support/help-18/users_guide/EmCommandLineforBatch.html)
and [thread controls](https://www.sonnetsoftware.com/support/help-18/Sonnet_suites/ThreadControlWin.html).

## 12. Implemented validation, provenance, and corrections

- Explicit ValueError validation replaced preparation asserts, surviving Python -O.
  Supported-grid dimensions and minimum clearances are enforced.
- Metadata/GDS identity and imported geometry are checked by rebuilding/XOR.
- Material mapping, ports, mesh, box, covers and inherited settings are checked.
- Return code/completion/errors and C_add header values protect simulation exports.
- Boundary/ambiguity/monotonicity guards protect resonance slope extraction.
- SNP input checks cover ports, normalization, finiteness and frequency order.
- Offline tests cover primitives, placement, N=1/4/10 geometry, grid rejection,
  stale metadata, N=1 model equivalence and simulation working-directory handling.
  Ten tests passed in normal and optimized Python during the audit.
- Selected saved networks passed reciprocity/passivity checks within tolerance.
- Newer runs include SHA-256 input/output hashes and running/failed/complete state.
- Shared simulation runner uses project.parent as cwd for relative outputs.
- Multi-cell fitting supports explicit SNP, N, L0 and output folder. Default L0
  remains target: omitting --l0-nh can choose an unintended fixed assumption.
- Generated extraction gained fit_capacitance=False to skip unnecessary final
  capacitance runs in L screening. Its default retains historical full extraction.
- Preparation accepts explicit GDS/output folder for isolated trials. The array
  runner gained --d-cell, --run-name, Cg error/tolerance fields and timing/resource snapshots.
- Export filename correction updates the final-project hash after the edit.

| Artifact | Meaning |
| --- | --- |
| Geometry JSON | Parameter record, not proof EM ran |
| validation.json | Import/preparation checks; solver_run:false describes preparation even after later EM |
| run_status.json | Script completion; confirm with logs and hashes |
| fit_parameters.json/results.json | Values, fixed assumptions, available provenance |
| simulation.log/cadd_*.log | Actual frequencies, duration, warnings/errors |
| mapped_template.sonx/mapped_import.sonx | Intermediate import files, not final configured project |
| sondata/ | Native cache tied to project/geometry/settings |
| PNG | Review only, not source geometry or numerical data |

Older PowerShell logs may be UTF-16; Python logs are plain text. Decode correctly
before parsing. A failed regex on wrong encoding does not mean zero analysis time.
Historical missing hashes/status are not retroactively fabricated.

An important earlier fit correction: generated w190/d60 fitted C0=193.540455,
Cg=8.241376 fF when fixed to TARGET L0=4.9194 nH. Refit of the same data with its
extracted L0=5.015715851 gave C0=189.820129, Cg=8.084707 fF. Both calculations
have explicit assumptions; changing fixed L0 changes fitted C, not saved S21.
Never describe the target-L fit as a measured-L fit.

## 13. Remaining assumptions and accuracy limits

1. The 2 um mesh was retained at the user's request. A 4 um trace/gap spans two
   cells. No convergence study establishes absolute L/C accuracy.
2. Selected L0 comes from native filter_v6; N=4 is generated geometry. Earlier
   w190 extractions were close, but w186 transfer was not independently checked
   using a generated single cell or a separate N=4 inductance measurement.
3. Spacing changes can affect effective loading/mutual inductance. Fixed L0 is
   the fit assumption here, not proof every array has identical inductance.
4. The identical-cell circuit approximates end loading and distributed fields.
   Its approximately 4.69 dB residual limits parameter interpretation.
5. ABS exports are interpolated. N=4's 2 MHz spacing can miss narrow peak maxima;
   do not claim precise insertion loss from sampled peak heights.
6. Reference-plane plausibility and zero warnings are established; calibration
   sensitivity was not measured.
7. Material values were inherited, not independently verified for film thickness,
   temperature or fabrication.
8. Component tolerances do not certify exact finite-array passband/phase/loss,
   ten-cell behavior, yield, or measured hardware performance.

These identify what a stronger signoff would need. They are not instructions to
repeat broad expensive validation for every small authorized tuning step.

## 14. Commands: saved-data fitting versus fresh EM simulation

Run PowerShell from the project root.

### Refit selected saved data without Sonnet

```powershell
.\.venv\Scripts\python.exe extract_circuit_params_ncell.py --no-simulate --snp sonnet_files\generated_ncell_4_w186_d56\generated_ncell_4.snp --n-cells 4 --l0-nh 4.922723945068658 --output-dir output_data\selected_d56_refit
```

This writes a separate fit folder. It does not recompute L0 or EM fields.

### New spacing trial

```powershell
# Example only: run if this new candidate is actually justified/requested:
.\.venv\Scripts\python.exe -u geometry\test_ncell_186.py --d-cell 54
```

This launches fresh EM analysis. The runner refuses existing experiment folders. For repeated studies, pass a fresh
--run-name, for example:

```powershell
.\.venv\Scripts\python.exe -u geometry\test_ncell_186.py --d-cell 56 --run-name manual_01
```

This creates separate files/folders ending in _d56_manual_01. A repeated identical
run name raises FileExistsError. Do not delete saved results to bypass this guard. The runner fixes N=4 and w=186; a different
N/w requires an explicit extension, not a renamed folder.

For custom geometry, instantiate grouped NCellParams, export GDS and metadata
together, and call prepare_sonnet.main(n_cells, gds_path=..., output_dir=...).
Then run em from that project folder and fit with explicit SNP/N/L0.

### Inductance tuning is a distinct multi-run experiment

```powershell
.\.venv\Scripts\python.exe -u geometry\tune_meander.py --max-trials 3 --refine-best
```

This describes historical logic, not an instruction to run now. The tuner has
fixed experiment roots and limited resume logic. Completed trial reuse checks
hashes, but repeated refinement can collide with best_refined or
last_sweep_cadd_50ff.snp. Use fresh roots for a new study. It is not yet a fully
resumable production optimizer.

### Cached processing or estimate only

```powershell
& 'C:\Program Files\Sonnet Software\19.52.2025.261\bin\em.exe' -v -AbsNoDiscrete sonnet_files\generated_ncell_4_w186_d56\generated_ncell_4.sonx
& 'C:\Program Files\Sonnet Software\19.52.2025.261\bin\em.exe' -N -v sonnet_files\generated_ncell_4_w186_d56\generated_ncell_4.sonx
```

Cached processing can rewrite export headers/bytes, making recorded hashes stale.
Preserve the original or identify the rerun separately; recompute provenance
when intentionally replacing data. Verify current exports and completion.
-AbsNoDiscrete cannot replace a new-geometry solution.

### Offline checks when code changes justify them

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s tests -p test_workflow.py -v
.\.venv\Scripts\python.exe audit_workflow.py
```

The artifact audit covers earlier default generated N=1/N=4 and initial L data.
It does not exhaustively audit every later spacing trial. Validate a new trial's
own metadata, logs, status, hashes and network.

## 15. Protocol for a future assistant

1. Read this document, workflow_state.json, selected metadata, fit/L0 JSON and
   relevant runner code. Confirm referenced files/hashes still agree.
2. Establish the requested target/tolerance and authorized variable. Current Cg
   task is complete; do not optimize further without a new user objective.
3. Distinguish targets from fixed measured inputs, native from generated geometry,
   fresh EM from fitting, and dense exports from independent EM frequencies.
4. Choose an informative candidate from measured sensitivity. Keep units/grid
   explicit; save all parameters and allocate a new experiment directory.
5. Reuse extraction only when geometry/material/loading assumptions remain
   appropriate. Explain any transfer assumption or calibration change.
6. Generate, validate and review; run only necessary simulation. Capture time,
   native frequency count, warnings and provenance.
7. Fit consistently with correct N, SNP and L0. Inspect seeds, plausibility,
   peaks, residual and reciprocity/passivity as appropriate.
8. Compare error to requested tolerance. If met, stop and save selection. If not,
   use measured sensitivity for the next supported-grid trial.
9. Report values, changed inputs, actual time, file links and material limitations.
   Do not assert an unmeasured slowdown cause or precise loss from coarse peaks.
10. Update this document/state when selection, goals, material, calibration or
    extraction method changes. Preserve earlier evidence.

Future production improvements: central parameter/config input, configurable
experiment roots, full resume/cache identity checks, and model validation across
geometries/array sizes. These are future work. The prototype already executes
the stages demonstrated above.


## 16. Editable user entry point added after this handoff

The preferred user-facing runner is now [run_filter.py](run_filter.py), with all
geometry and stage flags edited directly in the script. See [RUN_FILTER.md](RUN_FILTER.md)
for flag combinations, known/fresh L0 selection, output paths and reuse rules.
It defaults to selected w186/d56 geometry but does not change dataclass defaults
or earlier selected artifacts. prepare_inductance.main now accepts optional
GDS/project/output paths and sweep settings; extract_generated_inductance.main
accepts optional project/output paths and C_add values. Existing default entry
points and native tuning overrides remain compatible. No new EM simulations
were run to verify this addition; offline orchestration tests and a real fit of
existing selected S-parameters were used.
