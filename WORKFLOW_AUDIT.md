# Workflow audit — 2026-10-09

The saved geometry, imported Sonnet projects, simulation exports, and inductance extraction are internally consistent. The automated checks pass. This establishes a working prototype workflow, not a convergence study or proof that the equivalent circuit captures every electromagnetic effect.

This audit reviewed the Python workflow and existing artifacts, added validation and regression tests, and refit saved data. It did not run new EM simulations, change the 2 x 2 um mesh, change design targets, or overwrite historical fit results.

## Results reproduced

| Quantity | Result | Interpretation |
| --- | --- | --- |
| Single-cell extracted L0 | 5.015715851 nH | Slope from the 20, 30, 40, 50 fF added-capacitor sweep |
| Sweep regression R squared | 0.999999291 | Strong linearity over this sweep; not total physical accuracy |
| Single-cell fitted C0 | 188.401905 fF | Complex S21 fit using extracted L0 |
| Single-cell fitted Cg | 8.111116 fF | Same fit |
| Four-cell fitted C0 | 189.820129 fF | New fit of existing data using extracted L0 |
| Four-cell fitted Cg | 8.084707 fF | Same fit |

The earlier four-cell C0 = 193.540455 fF and Cg = 8.241376 fF used the **target** L0 = 4.9194 nH. Those results are not a calculation error, but must not be described as fits using the extracted inductance. The new fit is saved separately in `output_data/workflow_audit/generated4_extracted_l0/`.

The single-cell and four-cell fits use different objectives and loading environments. Their fitted capacitances need not agree exactly. The four-cell magnitude fit still has approximately 5.1 dB RMS residual over the fitting samples, so it is not an exact circuit representation of the EM response.

## What was checked

| Stage | Evidence and outcome |
| --- | --- |
| Parameterized geometry | Tests cover capacitor dimensions/area, meander runs and gaps, placement anchors, and arrays with 1, 4, and 10 interior cells. These are geometry tests, not EM simulations of all three sizes. |
| GDS and metadata | Regenerated polygons agree with saved metadata and GDS. Units are microns, precision is 1 nm, and metal uses layer 1/datatype 0. |
| Grid conversion | Imported polygons agree with the snapped GDS by polygon XOR; vertices lie on the 2 um grid. The original centered leads required a 1 um horizontal shift in the simulation geometry. The design GDS and snapped simulation GDS are intentionally separate artifacts. |
| Sonnet import | The saved N=1 and N=4 projects have the expected box sizes, geometry, two feed ports, and 2 um mesh. Stackup, material definitions, and relevant analysis settings were checked against the original template. |
| Feed calibration | Port attachments are on the intended feed edges. The 400 um reference-plane shifts end within the 412 um uniform side-feed region. This is a geometrical check, not a sensitivity study. |
| Added ideal capacitor | All five saved snapshots use the same geometry and the requested C_add values. The component terminals connect the top ground and the interior cap. Saved logs report completed runs without EM errors or warnings. |
| Touchstone data | Seven datasets passed finite-value, frequency-order, 50 ohm normalization, reciprocity, and passivity checks within numerical tolerance. Maximum singular values are below one; maximum reciprocity discrepancy is about 1.5e-10. |
| Inductance extraction | Independent regression reproduces 5.015715851 nH. Resonances decrease monotonically: 4.7441, 4.6420, 4.5460, 4.4556 GHz for 20–50 fF. |
| Circuit equations | The single-cell closed-form response and N=1 ABCD response agree numerically to approximately 2.5e-13. |

Detailed machine-readable evidence and export hashes: [checks.json](output_data/workflow_audit/checks.json).

## Corrections made during the audit

- Added `geometry/validation.py`: reject stale parameter/GDS combinations, stale derived metadata, unsupported grid dimensions, and insufficient clearances before import. The inductance preparation also checks that the source Sonnet geometry matches the current metadata and that component pins lie on the intended metal edges.
- Replaced preparation-script assertions with explicit validation errors, so running Python with `-O` cannot disable these checks.
- Made the shared simulation runner use the project directory as its working directory, keeping relative Sonnet export paths consistent.
- Added Touchstone input validation and a single-resonance guard. Inductance extraction now rejects ambiguous or boundary peaks, nondecreasing sweep frequencies, and nonphysical slopes.
- Added completion/error checks, artifact hashes, and a running/failed/complete status file for future generated inductance runs. These additions do not retroactively establish provenance that was not saved by older runs.
- Added explicit `--snp`, `--n-cells`, `--l0-nh`, and `--output-dir` controls to the multi-cell fitter, with a saved fit-parameter JSON. Defaults still use the target inductance; specify measured L0 explicitly when that is intended.
- Updated workflow documentation and added ten offline regression tests plus a repeatable saved-artifact audit.

## Limits that remain

1. **Mesh convergence is untested.** The requested 2 um grid is retained. A 4 um trace/gap spans only two cells; agreement with the previous layout does not establish discretization accuracy.
2. **Dense exports are not dense independent EM solves.** Sonnet ABS interpolation produced the 0.1 MHz single-cell export spacing. The regression standard error of about 0.003 nH excludes mesh, interpolation, calibration, and model errors. Independent frequency checks near resonances would be needed to quantify interpolation error.
3. **The four-cell export is too coarse for reliable peak insertion loss.** Its 2 MHz spacing can miss narrow peak maxima. The observed baseline/generated peak shifts of approximately 8–12 MHz are useful comparisons, but the sampled peak heights are not a precise loss comparison.
4. **All six physical structures participate in the four-interior-cell EM simulation.** End cells are connected to the feeds and load the interior array. Four visible transmission peaks do not prove the end cells have no modes or are electromagnetically absent. Reference-plane arrows do not remove their modes. The reduced four-resonator model is an approximation that needs validation across array sizes.
5. **Transferring single-cell L0 to longer arrays is an assumption.** Mutual coupling and end loading can affect fitted effective parameters. The four-cell fit is a consistency check, not an independent measurement of its inductance.
6. **Calibration and matching need separate assessment.** Reference planes are geometrically plausible, but their sensitivity has not been measured. The saved physical feed impedance is approximately 56.89 ohm, while S-parameters are normalized to 50 ohm. This is allowed, but normalization does not make the physical feed matched.
7. **Automation is modular but not yet an optimization controller.** Standalone component generators contain example parameters; changing one does not automatically change the array. Array metadata is the parameter record for import. Before optimization, each candidate should have its own parameter/artifact manifest and results directory, with explicit propagation of extracted L0 to subsequent fits. Existing folders are reusable working locations and are not an immutable experiment database.

Sonnet documents that shifted reference planes require a suitable uniform transmission-line region: [Shifting Reference Planes](https://www.sonnetsoftware.com/support/doc_19.52/Topics/Analysis_Setup/De-embedding/Shifting_Reference_Planes.htm). Its ideal components are incorporated through circuit analysis and internal component ports, which supports reusing an unchanged EM solution when only the ideal capacitor value changes: [Components](https://www.sonnetsoftware.com/support/doc_19.52/Topics/Components/Components.htm).

## Reproduce the checks

Run from the repository directory in PowerShell:

```powershell
.\.venv\Scripts\python.exe -m unittest discover -s tests -p test_workflow.py -v
.\.venv\Scripts\python.exe -O -m unittest discover -s tests -p test_workflow.py
.\.venv\Scripts\python.exe audit_workflow.py
```

Both normal and optimized Python runs passed all ten tests. Compilation checks passed for the geometry package and audited workflow scripts. Tests mock native execution where applicable; the newly added preflight code was not followed by a new native EM run in this audit.

To reproduce the four-cell fit with the extracted inductance, without simulating:

```powershell
.\.venv\Scripts\python.exe extract_circuit_params_ncell.py --no-simulate --snp sonnet_files\generated_ncell_4\generated_ncell_4.snp --n-cells 4 --l0-nh 5.015715850762367 --output-dir output_data\workflow_audit\generated4_extracted_l0
```
