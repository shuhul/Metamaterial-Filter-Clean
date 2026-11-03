# First taper: selected design and iteration handoff

> This document describes the historical 6.06 nH / 100 fF / 56.7 fF design.
> Active targets are in `three_taper_targets.json`; the corrected unit-cell-based
> T1 workflow is `three_tapers_v3.py`, with status in `three_tapers_status.json`.
> Its input gap has no ground strip, the regular-junction strip is centered,
> and the active fixed-target complex S21 RMS limit is 0.025.
> Do not use this old selected T1 for the replacement three-taper specification.

Current T2 work runs `python -B -m geometry.taper.taper2_idc` from the project
root (using `.venv\Scripts\python.exe` on Windows). Edit the grouped `SETTINGS`
and unique `CASE` in `taper2_idc.py`; it delegates geometry and simulation to
`taper2_plain.py` and the common response analysis in `taper2_v3.py`.
The frozen accepted T1 is in
`taper2_T1_reference.json` and measured progress/blocker in `taper2_progress.json`.
The detailed current campaign is
`../../sonnet_files/taper_exploration/three_tapers_v3_T2_2026-10-26/README.md`.
These files distinguish fitted-curve error from the fixed-target S21 score;
T2 is not yet accepted. License connectivity was restored on October 27;
parameter-only Sonnet iterations are active with one four-thread job at a time.
The campaign README records completed trials and distinguishes prepared layouts
from simulated and accepted designs.

The default driver remains serial. The user also authorized a two-job trial:
`parallel_trials.run` launches two isolated worker processes with separate native
projects and output logs. Each solver is capped at four threads; the shared
two-slot lease excludes a third solver and excludes serial callers while either
slot is occupied. A failed concurrent attempt is archived and retried serially
after both original workers stop. Batch evidence is under
`investigations/parallel_t2_2026-10-28`; this is optional, not an unlimited solver
pool. Shared identical-body inductance caches are completed before concurrent
readers use them. Only the parent can select/finalize a passing design.

`coupling.w_input_stub` optionally controls the horizontal feed span before the
input rail independently of T2 cap width. Its default `None` reproduces historical
geometry and metadata. Changing it translates the IDC and every downstream
component together; their relative geometry remains fixed. The external ports
remain50ohm normalized and retain400um reference planes. This control tests
residual feed loading without changing the electrical acceptance model. The
current trial ledger is
`output_data/taper_exploration/three_tapers_v3_T2_plain_2026-10-28/thin_rail_trials.json`
relative to the project root; per-trial `response/proposal.json` records the reason.

### Explicit-capacitor inductance extraction restored October 28

The active T2 driver now uses `explicit_inductance.py`, which calls the original
`extract_generated_inductance.py` with real ideal CAP components in the Sonnet
project. One fresh native CAP calculation is followed by three cache-only
capacitor changes (`-AbsNoDiscrete`); each subsequent log must confirm zero new
EM frequencies. A separate CAP cache is created instead of inheriting the
four-port cache, which crashed during the attempted conversion. Original
fixtures and the failed-run logs were preserved.

Current checks, using the same added capacitances as the prior sparse method:

| Body | Explicit CAP L (nH) | Previous sparse L (nH) | Change | Target error |
|---|---:|---:|---:|---:|
| Frozen T1 | 6.00141179 | 6.00627210 | -0.081% | +0.526% |
| T2, w148/pad200 | 4.31103489 | 4.31481655 | -0.088% | -0.438% |

T1 evidence is in `output_data/taper_exploration/explicit_L_checks_2026-10-28/T1`
relative to the repository root; its native CAP project is under the matching
`sonnet_files/taper_exploration/explicit_L_checks_2026-10-28/T1` directory.
T2 central evidence is `output_data/taper_exploration/explicit_L_cache/0eb5a2f4f41f1395`;
its native project is `sonnet_files/taper_exploration/explicit_L_cache/0eb5a2f4f41f1395/Lcal/generated_ncell_1.sonx`.
Each result includes capacitor values, sampled resonances, regression, hashes,
timing, and new-EM-frequency counts; all four native projects/exports are saved.

The T1 check used 0.1 MHz exports and took365.3 seconds including extraction.
Routine T2 uses 1 MHz exports:133.3 seconds total, including111.3 seconds for the
first native call and approximately7.0 seconds for each cached call. These are
completed run times, not a controlled benchmark against the sparse method.
Both started with fresh CAP caches. No claim of comparative speedup is made.

Identical T2 calibration metal reuses the checked central result across assembly
coupling changes. The body fingerprint, source hashes and GDS XOR are checked;
a body change needs a new extraction. Case-specific copies are in `Lcal_explicit`.
Historical `Lcal/results.json` sparse estimates remain for comparison and are
not the active driver's inductance source. Sparse four-port interpolation is
retained only as a diagnostic. The original audit evidence is in
`investigations/inductance_evidence_review_2026-10-28`.

Both independently measured taper inductances and all regular-cell values stay
fixed in the five-capacitance heuristic fit. Acceptance still uses the fixed
nominal circuit over4.65-5.40GHz with7501 points and one matched phase delay;
complex S21 RMS must be at most0.025. Replacing extracted L does not change that
target score for an existing EM export. No T2 design is accepted yet.

The first taper benchmark meets all four targets within the inherited 2% screening
criterion. This is a result from the specified Sonnet setup and equivalent-circuit
extraction, not a claim of 2% absolute accuracy in a fabricated device.
Targets are retained in `first_taper_targets.json`; the exact selection, paths,
measurements, errors, and audit are in `selected_design.json`.

| Parameter | Target | Selected | Error |
|---|---:|---:|---:|
| L (nH) | 6.06 | 6.014106 | -0.757% |
| C0 (fF) | 100 | 101.492048 | +1.492% |
| Next-taper/input coupling (fF) | 56.7 | 56.151540 | -0.967% |
| Taper-to-regular coupling (fF) | 8.2865 | 8.229325 | -0.690% |

## Geometry and parameter meanings

All geometry dimensions are microns. `TaperParams` groups the cap dimensions;
`taper.meander` is the reusable existing `MeanderParams` object. Widths are
horizontal (`w`), lengths vertical (`l`), separations clear edge-to-edge (`d`),
and metal widths/thicknesses are `t`. Meander `n` counts horizontal runs, including
the two partial end runs, rather than U-turns.

- Cap width 328, top thickness 14, left/right thickness 24/6.
- Left leg length 584, right leg length 310. Leg lengths include the top bar.
- Meander width 240, total length 596, top lead 10, 58 runs, trace 4, clear gap 6.
  The body is 574 high and the derived bottom lead is 12. Bends remain square.
- No grounded pads inside the taper. The regular cell retains its two pads,
  cap 302 x 600, meander 186 x 596 with 72 runs and 4 um trace/gaps.
- Input feed-coupling finger: thickness 20, outer vertical length 536,
  input separation 6. Its top bar is 328 x 14. It has no meander.
- Taper-to-regular separation 56, with a 4 um grounded strip. Regular-to-output-end
  separation is also 56. Side grounds are 412 wide with 14 um end separation.
- Cap tops share y=-8; the entire upper ground boundary is y=0. There is no
  step at the interface. All meanders end at y=-618. The common extent is
  8 + 14 + 596 = 618, rather than 8 + 596.
- The box is 2230 x 1500. Top ground height 448. The actual Sonnet mesh is 2 x 2,
  with 400 um reference planes and 50 ohm normalized ports.

The unit origin is centered 8 um above the cap top, with GDS +y upward. Original
fixture coordinates span y=-1052..448. The import GDS uses box coordinates
0..2230 by 0..1500; Sonnet uses downward y. Regular center-lead edges can start
on odd-micron coordinates and are snapped by up to 1 um, using the same positive
half-grid rule as the existing workflow. Imported metal is checked by XOR.

## Saved files and how to run

- `output/first_taper_unit.gds` and `.png`: the standalone reusable taper.
- `output/first_taper_fixture.gds` and `.json`: the full benchmark source.
- `output/first_taper_fixture_grid2.gds`: the exact snapped import artifact.
  Use the source with its metadata in Python; the grid2 artifact is for direct
  import/inspection and does not use the unsnapped source metadata.
- `output/first_taper_fixture.png`: preview of actual Sonnet metal and planes.
- `../../sonnet_files/taper_exploration/final_designs/finger536_t20_right310_input6_j56_final/generated_ncell_2.sonx`:
  the actual simulated project. Open this file to review settings and geometry.
- `../../output_data/taper_exploration/first_taper_s21_fit.png`: main fit plot.
- The selected case folder contains `sparameters_comparison.png`, `trial.json`,
  `fit.json`, `audit.json`, `response.snp`, `project.sonx`, native logs, and the
  saved geometry/response identity. The `_isolation` result folder contains L
  sweep snapshots and regression plots. Its last sweep project/export is 50 fF;
  it is a calibration fixture, not the final mixed-response project.

From the project root:

```powershell
.\.venv\Scripts\python.exe geometry\taper\run_taper.py
```

The default flags only refit the saved selected result, with no new simulation.
All settings live inside `run_taper.py`; there are no CLI arguments. For a new
geometry, choose a new `CASE_NAME`, set `SIMULATE_MIXED=True`, and set
`SIMULATE_INDUCTANCE=True` when changing the taper unit itself. Leave
`FIT_RESPONSE=True` to get the plots and fitted values. If only the input finger
or cell spacing changes, keep the matching `INDUCTANCE_RESULT_FILE` and
`SIMULATE_INDUCTANCE=False`. Parameter checks reject an unrelated L reference.
`FIT_RESPONSE=False` permits simulation/export without fitting. An inductance-only
run is also possible with `SIMULATE_INDUCTANCE=True`, `SIMULATE_MIXED=False`, and
`FIT_RESPONSE=False`. Completed trial names are protected; use a new name.

## What the code does

1. `layout.py` constructs the asymmetric cap and places `make_meander` at the
   cap underside. `build` assembles the input coupling fixture, taper, regular
   resonator and regular output end, feeds, common grounds and weak-gap strips.
   It writes one TOP GDS library with 1 um units, 1 nm precision, layer 1:datatype 0.
2. `prepare_sonnet.main(..., geometry_validator=layout.validate)` rebuilds the
   expected source geometry and validates the GDS. It snaps an import copy to
   the 2 um mesh and invokes Sonnet `gds.exe`. Materials, substrate, covers and
   controls come from the existing regular-cell Sonnet template. Native polygons
   must XOR-match the snapped source. Ports are attached to the unique box-edge
   feed segments; vertex numbers are discovered rather than hardcoded.
3. `experiment.measure_unit` builds a three-unit isolation fixture and adds the
   original ideal shunt component using `sonnet.add_capacitor`. Component pins
   bridge the top-ground/cap gap. This ideal component appears only in the
   isolation project, never in the final physical GDS or mixed-response project.
4. `extract_generated_inductance.main` sweeps 20/30/40/50 fF C_add and tracks each
   sampled S21 maximum. It regresses 1/(2*pi*f)^2 against C_add in farads; the
   slope in henries is L. The first point solves EM; subsequent points reuse
   the unchanged geometry cache. Every point preserves project, Touchstone,
   log, peak and hashes. L for this selection is 6.014105932 nH,
   with R_squared=0.9999955055 and regression standard error 0.009016 nH.
   This statistical error does not include mesh/material/model errors.
5. `sonnet.configure` sets the short ABS sweep; `sonnet.simulate` runs em and
   records completion and input/output identity. Mixed screening uses 4..5.8 GHz,
   target 900, producing 901 samples at 2 MHz. The selected fresh run used 7 actual
   EM frequencies in about 1 minute 47 seconds. Each isolation geometry needed
   about 1 minute 10 seconds, followed by cheap cached C_add points.
6. `fitting.py` builds a nonuniform lossless ABCD chain:
   input Cg -> taper shuntLC -> junction Cg -> regular shuntLC -> output Cg.
   L_taper comes from its own sweep, and L_regular remains the previously measured
   4.922723945 nH. Both self capacitances and all three coupling capacitances are
   fitted independently. C0 is never pinned to 100, and junction/output Cg are
   never tied together. The fitted regular values are C0=188.450589 fF
   and output Cg=8.312406 fF, close to the retained regular design.
7. `audit.py` rechecks source/actual metal, units, stack/materials, mesh, planes,
   ports, archived sweep hashes, Touchstone peak frequencies and L regression.
   It checks actual completed logs for zero errors AND zero warnings. Final
   project/response hashes must match the recorded run. `index.py` refreshes the
   trial index in `output_data/taper_exploration/summary.json`.

The measured circuit L is not entered as Sonnet metal inductance. The original
SUPERCONDUCTOR material's Ls=0.0388 pH/square is inherited unchanged, as are the
substrate and other material settings. The 6.06 nH target is used for comparison;
the fitted circuit uses measured L, not the target value.

## Extraction model change: residual fixture phase

After switching to the original taper-unit's finger-feed topology, fitting raw
complex S21 biased C0 relative to magnitude-only fitting. The circuit nodes are
not identical to the 400 um reference-plane locations, and the feed/coupling
fixture adds smooth phase that the minimal lumped chain omits. We added bounded,
nonnegative matched port delay sections as an effective fixture correction.
This is an approximation to residual fixture phase, not a measured geometric
feed length or a hidden capacitor adjustment.

Total delay is fitted with complex S21; its division between ports is then
chosen using S11/S22, without changing the extracted L/C values. Selected total
is 4.985216 ps, with left/right
4.395758/0.589458 ps. Broad positive
capacitance bounds and multiple starting points are used. All original and
alternative fits remain saved, so this change is explicit and reviewable.

| Fit | C0 fF | Input Cg fF | Junction Cg fF | Complex S21 RMS |
|---|---:|---:|---:|---:|
| amplitude | 101.0984 | 56.5114 | 8.3372 | 0.024440 |
| complex | 101.3001 | 56.7450 | 8.1043 | 0.013892 |
| complex_with_delay | 101.4920 | 56.1515 | 8.2293 | 0.006472 |

The primary delay-corrected fit and the magnitude-only cross-check both meet
all four targets within 2%. Reflection discrepancies remain: complex RMS is
0.06852 for S11 and 0.07627 for S22. This is
an effective circuit model for screening, not a complete model of the feed ends.
The strong coupling is measured against a feed finger representing the facing
next-taper interface. A complete next taper resonator and the whole multi-stage
filter still require a joint EM check when those stages are designed.

## What was tried and learned

- Original taper isolation measured L 6.0484 nH, C0 99.759 fF and symmetric
  coupling 57.184 fF, already near the requested values. Its uniform left-side
  ground is only 384 um wide, so the baseline used 360 um reference planes rather
  than the invalid original 600 um setting. Generated fixtures use 412/400.
- The native taper-unit project has a 26 um ground step: left ground ends at
  y=422, while an extra rectangle extends the right ground to y=448. The cap
  tops themselves were aligned. Flattening the ground raised fitted C0.
- A fully modular flat control gave about 115..117 fF. Trying fewer meander runs
  and wider gaps changed L but barely changed C0 after restoring L with width.
  This empirical result redirected C0 tuning to the capacitor legs.
- 58 runs with 6 um gaps needed width 240 to restore L near 6.06 nH. Width 228 gave
  only 5.747 nH. This valid wider-gap meander was retained in the selection.
- Shortening both legs to 460 gave C0 104.44 fF; 420 gave C0 99.98 fF, but strong
  coupling became too small. The 2 um input gap trial was not selected.
- Initial modular trials copied a full input end cell. Inspection of the original
  taper-unit fixture showed that its input is instead a feed bar and finger,
  without a meander. Later trials use this topology and tune the input finger
  independently. Early/full-end and later/finger results are not interchangeable.
- The selected asymmetry keeps the left leg long for strong coupling, shortens
  the right leg for low C0, and adjusts the ground-screened junction gap for weak
  coupling. Input finger overlap and thickness provide an independent strong
  coupling knob. This avoids a displaced top-ground boundary.
- The last unsuccessful trial, right 330 and junction 60, passed L/C0/input but
  missed junction by 2.29%. Right 310 and junction 56 brought every target inside.

Old trials and the original `tests/` projects were retained. Early logs from the
shared L helper print its old regular-cell target; their taper result JSON was
corrected to 6.06 without changing the regression. The helper now accepts an
explicit `target_l0_nh`, so later taper logs use the correct comparison target.

## Cache/output-grid check and stability

Cache-only requests for a denser grid did NOT change Sonnet's existing 901-point
export in this installation. We verified the actual count/range, restored the
original 900-target project and exact response snapshot, and explicitly recorded
`grid_refinement_effective=false`. No finer-grid convergence result is claimed.
Those cache calls solved zero new EM frequencies. See the
[Sonnet cache-only option documentation](https://www.sonnetsoftware.com/support/help-18/users_guide/EmCommandLineforBatch.html).
`sonnet.simulate` now rejects an export that does not match the requested
frequency range/grid, instead of treating an old export as a new result.

As an independent fitting-sampling check, even and odd 901-point subsets were
fitted separately at 4 MHz spacing. Their C0, input Cg and junction Cg differ
from the full fit by less than 0.001 fF. These results are in the selected case's
`sampling_check/`. This checks fit sampling stability, not physical mesh convergence.

Offline verification: 20 existing workflow tests plus 7 taper tests passed,
including energy conservation/reciprocity, the uniform-chain limit, synthetic
recovery with fixture phase, stale GDS rejection, geometry grounding/alignment,
input-finger topology, and refusal of a stale cache-only export. Running the
editable script with default flags was verified to refit saved data without EM.

Shared changes outside this folder are limited to an optional geometry validator
in `prepare_sonnet.py`, resolved explicit paths and an optional report target in
`extract_generated_inductance.py`, the geometry README link, and offline tests.
Existing regular geometry/defaults, selected regular results, and original test
projects are preserved.
