# Editable filter workflow

Open [run_filter.py](run_filter.py), edit the block marked EDIT YOUR SETTINGS
HERE, then run it. No command-line arguments are required.

```powershell
cd C:\Users\shuhul\Downloads\filter_new_automation_clean
.\.venv\Scripts\python.exe geometry\run_filter.py
```

The defaults reproduce the selected geometry: 4 interior cells plus 2 end cells,
186 um meanders, d_cell=56 um. The default run name is my_filter_01. All geometry
parameters are explicitly grouped in PARAMS. Targets remain in constants.py.
The validated 2 um grid and 400 um reference planes are preserved; materials and
stackup are inherited from the original template. This is an editable runner,
not an automatic geometry optimizer.

## Main switches

| Flag | What True does |
| --- | --- |
| GENERATE_GEOMETRY | Write array GDS, JSON and preview from PARAMS |
| PREPARE_SONNET | Import GDS, configure ports and simulation, save project |
| SIMULATE_INDUCTANCE | Generate/use matching N=1 companion, run C_add sweep, extract L0 |
| SIMULATE_ARRAY | Run the N-cell adaptive EM analysis once |
| FIT_S21 | Fit C0/Cg with selected L0 and save s21_fit.png |

Use a fresh RUN_NAME for changed geometry. Existing geometry/projects are not
silently overwritten. Flags are genuine stage switches: when a simulation flag
is False, that simulation is not launched.

## Common modes

| Mode | Generate | Prepare | Inductance | Array EM | Fit |
| --- | --- | --- | --- | --- | --- |
| New design with known L0 (default) | True | True | False | True | True |
| New design with fresh L0 | True | True | True | True | True |
| GDS and preview only | True | False | False | False | False |
| Prepare a project without simulation | True | True | False | False | False |
| Refit a saved run | False | False | False | False | True |
| Rerun a saved array | False | False | False | True | True |
| Extract L0 only for a new design | True | True | True | False | False |

For saved-run modes, keep the same RUN_NAME and PARAMS. They apply to projects
created by this runner; historical trial folders have a different layout and
lack its simulation_record.json. The script checks GDS/parameters, native
geometry and recorded project/export hashes before fitting stored data.
Changing a geometry parameter while fitting old data is rejected. Changing the
chosen fixed L0 is allowed as an explicit circuit-fit assumption.

Array frequency settings take effect when preparing or simulating a project.
Fit-only mode uses the actual saved frequency vector; it does not create new
frequency data. Repeated analyses save a new analysis_<UTC timestamp> subfolder,
so earlier fits are preserved.

## Which inductance is used?

- With SIMULATE_INDUCTANCE=True, the script creates a generated N=1 companion
  using the same unit, spacing, feeds, grounds and material template. It adds
  the ideal parallel capacitor and runs C_add=20/30/40/50 fF. The slope gives
  L0, which is passed to the array fit in that invocation. No final tiny-cap
  capacitance fit is performed on the companion; its last test value is 50 fF.
- With the flag False, USE_SAVED_INDUCTANCE=True and INDUCTANCE_RESULT_FILE=None, it reuses this RUN_NAME's
  previously extracted L0 if one exists. Otherwise it uses FIXED_L0_NH, which
  defaults to the selected known 4.922723945 nH.
- Set INDUCTANCE_RESULT_FILE to an explicit extraction results.json to select
  another saved result. New generated results check complete parameter identity
  and template hash. Historical native results can only check their stored
  meander span; remaining transfer assumptions are explicitly printed.

Re-extract when meander/trace/gap/material or substantial loading changes make
fixed L0 inappropriate. Small spacing sweeps may deliberately hold L0 fixed,
as the earlier experiment did. A fixed value is an assumption, not a new
measurement. To ignore previously extracted L0 and force your manual value,
set USE_SAVED_INDUCTANCE=False and edit FIXED_L0_NH. A newly enabled sweep
always supplies its measured result regardless of this reuse flag.

Default capacitor-sweep range: 4.35-4.95 GHz, Target=1500. This suits the current
geometry's test-cap resonances. Expand/change it for significantly different
geometry; boundary/ambiguity checks reject out-of-range/unresolved peaks.
Default array range: 4.5-5.5 GHz, Target=500. Target controls output density;
actual Sonnet spacing may differ from range divided by Target. C_add is in pF;
L0 is in nH; geometric dimensions are um. The n=72 meander count includes the
partial first/last horizontal runs; bottom lead length is derived.

## Where files appear

For RUN_NAME="my_filter_01":

```text
geometry/output/my_filter_01.gds
geometry/output/my_filter_01.json
geometry/output/my_filter_01.png
sonnet_files/my_filter_01/array/generated_ncell_4.sonx
sonnet_files/my_filter_01/array/generated_ncell_4.snp
output_data/my_filter_01/s21_fit.png
output_data/my_filter_01/fit_parameters.json
output_data/my_filter_01/configuration.json
output_data/my_filter_01/run_status.json
output_data/my_filter_01/stage_timings.json
```

The console prints the absolute plot/results paths. For repeated saved-data
analysis, plots go into analysis_<timestamp> under the result folder. The root
latest_run.json points to the most recent successful analysis. It records which
stages ran. Exact response/project snapshots preserve data across reruns.

When extracting L0, additional files are:

```text
geometry/output/my_filter_01_single.gds
sonnet_files/my_filter_01/inductance/generated_inductance.sonx
output_data/my_filter_01/inductance/results.json
output_data/my_filter_01/inductance/cadd_*.sonx / .snp / .log
```

On repeated analyses the extraction outputs belong to that analysis subfolder;
inductance_record.json points to the current extraction. A stopped run is marked
failed. There is no automatic full-stage resume for a partly generated run;
use a new run name or deliberately reuse valid stages with their flags disabled.

## Checks and limits

The new runner was checked offline: real saved-data fitting reproduces the
selected Cg, geometry-only mode avoids native execution, optional extraction
uses the measured L0, and stale inputs/overwrites are rejected. Native execution
was mocked for those orchestration tests; no new EM run was launched for code
verification. Existing native helpers retain their previous default interfaces.

The same physical limits remain: 2 um convergence is untested, single-cell L0
is transferred to N-cell fitting, and the reduced circuit is approximate. See
[WORKFLOW_HANDOFF.md](WORKFLOW_HANDOFF.md) for detailed assumptions and history.
