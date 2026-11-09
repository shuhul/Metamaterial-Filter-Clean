# Current investigation: selective ground meshing and shorted-circuit extraction

Acceptance: all physical C0/Cg/L within 1% of targets, circuit complex-S21 RMS below 0.01. The user relaxed accelerated-EM screening from 0.005 to **0.01** complex-S21 RMS of the original numerical reference and authorized rough screening followed by fine verification. These are separate tests.

## Files to inspect

- Corrected geometry, original mesh: `geometry/input_IDC_trimmed/project.sonx`.
- Corrected geometry, selective coarse-ground preview: `geometry/input_IDC_trimmed_coarse_ground_preview/project.sonx`.
- Close-up visually checked: `geometry/input_IDC_trimmed_coarse_ground_preview/IDC_closeup.png`.
- Mesh regions, not an actual subsection plot: `mesh_benchmarks/coarse_ground_interior/mesh_region_plan.png`.
- Live sequential worker: `campaign_status.json`.
- Plan and primary research links: `PLAN.md`.

**The mesh benchmarks deliberately retain the original untrimmed metal**, so the mesh effect is compared independently of the IDC geometry correction. They are not the corrected design to use for subsequent design work. The actual extraction runs use the trimmed geometry. Neither original selected files nor old simulation caches are modified.

## Completed numerical checks

Fresh independent direct calculations at 4.85, 5.0 and 5.15 GHz reproduce the archived original-mesh S21 interpolant with RMS **0.00000545**. These are checks at three frequencies, not a full-band mesh-convergence claim.

| Setting | Subsections | Three-point elapsed time | Complex S21 RMS versus fresh baseline | Decision |
|---|---:|---:|---:|---|
| Original medium-speed mesh, edges on | 14,538 | 811.08 s | reference | retained reference |
| Edge mesh disabled globally | 3,328 | 158.36 s | 0.247035 | rejected |
| Coarse distant ground interiors, fine active edges | 14,251 | 258 s for first point | 0.11525 at first point | rejected; remaining pilot solves stopped |
| Remote box/port extension by 10 um, same physical reference plane | 14,543 | 269 s for first point | about 0.000018 from printed first-point rows | accurate at pilot point but no timing benefit; deprioritized |
| Native polygon merge, 261 editor polygons to four | 14,538 | mesh estimate only | not solved | same subsection count; deprioritized |

Global edge removal was 5.12 times faster, but fails accuracy badly; that speedup cannot be used. Initial per-polygon minimum-size tests with the fine global setting increased subsection count (30,616 / 26,314 / 20,493 for minima 2 / 4 / 8 cells) rather than improving on the existing middle global setting. Increasing maximum subsection size to 1,000 cells only reduced the count to 14,420. These last numbers are mesh estimates, not solved timing improvements.

The distant-ground partition reduces the subsection count by only 287 (1.97%). This suggests that the original automatic mesh already avoids spending most of its unknowns in those broad interiors. Its first-point error alone makes the intended three-point RMS at least 0.06654, even if the other points had zero error, so the failed pilot was stopped early. This specific partition is not an acceptable acceleration; it does not prove that every possible selective ground mesh would fail.

Another controlled experiment extended only the remote right port/box boundary by 10 um, changing X cell count from 2,155 (5 x 431) to 2,160 (2^4 x 3^3 x 5). Right reference length increased from 400 to 410 um, retaining the same physical reference location. This is a boundary experiment, not identical-box meshing. First-point S21 agreed closely, but elapsed time showed no useful improvement; remaining points were stopped and the candidate deprioritized. This is not a completed three-point/full-band validation. Source rationale is Sonnet's Japanese distributor's FFT/cell-count timing guidance linked in `PLAN.md`.

Native `soncmd -MergePolygons` reduced the number of editor polygons from 261 to four while preserving canonical metal XOR, port locations and impedances. Native mesh estimation still reports exactly 14,538 subsections and the same memory as baseline, so a new full timing sweep was not justified. The solver appears already to account for those internal boundaries.

No acceleration has passed both accuracy and benefit checks. The original mesh was selected and the original archived full-band response reused as its reference; the fresh three-point independent checks above verify it without a redundant 20-minute same-settings resimulation. Trimmed unshorted/shorted geometry experiments now proceed on that mesh.

## Circuit baseline and shorting rationale

An exploratory fit of the existing original response, with every capacitance constrained to nominal +/-1%, all calibrated inductances fixed, and port delays free, reached S21 RMS **0.05070**. All ten capacitances hit their target bounds. This is not a successful extraction or a proof of a global minimum; it establishes that merely tightening the old fit is insufficient. Saved evidence: `joint_fits/configs_0/S21_target_bounds_1pct/fit.json` and `joint_S21.png`.

Prepared shorted variants retain every original resonator and pad and add metal bridges across the first one, two or three series couplings. Shorting a coupling does not delete earlier shunt tanks. The joint model retains those shunts and shares physical capacitances across the configurations.

Ideal-circuit synthetic Jacobians suggest that adding the first two shorted variants improves conditioning substantially. This is only experiment-design evidence: physical bridges can alter local fields and add inductance, so actual shared-parameter predictions, residuals and parameter profiles must verify transferability. Synthetic conditioning assumes the ideal circuit is correct and is not a confidence interval.

## Automation scope

The serial four-thread worker completes conservative mesh screens, validates a usable setting against the saved full-band short-assembly response, then simulates trimmed unshorted and shorted configurations. It performs shared-capacitance fits with calibrated L fixed, broad diagnostic bounds and nominal +/-1% bounds, saving both S21-only and all-S fits. It stops for evidence review before geometry tuning; this is not a claim that the final geometry already meets the new targets. No new full 10-cell sweep is launched.

The first sandboxed solver attempt could not reach the license server. Re-running the authorized solver with license-server access succeeded; the completed timings above are from the successful run. The originally stopped 114-point full-filter run remains stopped and preserved.

## Prepared later geometry controls

`geometry_controls/manifest.json` lists +/-2 um input-IDC span and T3 pad-width changes. These are native-read checked and local-XOR checked, with all meanders, old downstream metal and corrected IDC termination preserved. They have **not** been simulated. They are available for measured local sensitivities after shared extraction is tested, not a substitute for that test.

## Further solver research

The installed release also documents a Beowulf solver, which requires a separate license. It can benefit large models, while setup overhead can make small models slower. It has not been benchmarked or claimed available here. [Official release-19 solver documentation](https://www.sonnetsoftware.com/support/doc_19.52/Topics/Analysis_Setup/EM_Options/EM_Solvers.htm).

## Expanded speed tests after the user's 0.01 screening allowance

The eight-thread, same-metal, same-mesh direct 4.85 GHz pilot completed successfully: native frequency time **156 s versus 270 s** for the fresh four-thread reference, about **1.73x faster**. Complex S21 differs by **9.30e-11**; all-S RMS is 3.05e-10. Native total including setup was 164 s. The wrapper's single-point postprocessor initially rejected one sample, so native timing is authoritative; approximate wall recovered from the last log write is 170.91 s. The completed EM calculation was retained and only postprocessing repaired, without repeating it. Evidence: `mesh_benchmarks/eight_threads/status.json`, `em.log`, `project.snp`. A cold full-band corrected-geometry comparison is queued before treating this as sweep-wide verification.

**Full-band eight-thread verification is now complete.** Independent cold folders and byte-identical `.sonx` files, same four directly solved frequencies: **1094.40 s (18 min 14 s) at four threads versus 628.99 s (10 min 29 s) at eight threads**, a measured **1.740x wall speedup**. Native totals are 18 min 7 s versus 10 min 20 s. Dense S21 and all-S exports agree exactly at their saved precision (RMS 0.0); this is not a claim of exact arithmetic or continuum accuracy. Native errors/warnings are zero. This thread-count improvement preserves the original numerical method. See `speed_followup/eight_threads_full_band/comparison.json`, `S21_comparison.png`, and both run logs. The direct pilot above resolves smaller numerical differences than this exported full-band comparison.

These trials use an investigation-local capped runner (`capped_solver.py`), not an unrestricted global thread increase. It reserves both solver slots, caps its own child at eight declared threads/ten core equivalents, and stops its own child on sustained host load above 85%. Other processes are not targeted.

Intermediate mesh estimates on the corrected metal: conformal full coils fail native meshing (box-edge touching); conformal straight coil parts increase the count to 15,824; coils alone without edge refinement reduce it to 7,187. The latter completed a pilot in 65 s per frequency (76.16 s wrapper wall) but complex S21 error is **0.08086 at 4.861 GHz**, rejecting it as a direct substitute. Larger interior minimum sizes with edges retained give 11,971 / 11,343 / 11,029 subsections for minima 100 / 200 / 500 cells. Actual S-parameter pilots follow; counts alone are not measured speedups. See `intermediate_meshes/estimates.json`, `speed_followup/additional_estimates.json`, and `speed_followup/pilot_results.json`.

The middle slider setting shown by the user is already the baseline. Sonnet documents only three slider positions, so intermediate tradeoffs need polygon/Tech Layer controls. Fewer exported ABS samples alone do not materially speed up the expensive EM calculations; selective direct screening is useful, but must not be presented as a converged full curve. [Release-19 subsection controls](https://www.sonnetsoftware.com/support/doc_19.52/Topics/Subsectioning/Subsectioning.htm), [ABS sweep documentation](https://www.sonnetsoftware.com/support/doc_19.52/Topics/AnalysisPlan/Adaptive_Band_Synthesis.htm).

## First real shorted-data transfer check

Corrected unshorted and first-shorted native sweeps are both complete, with no EM errors, using four direct frequencies each. Shared calibrated-L, broad-capacitance fitting gives S21 RMS **0.02566 / 0.03317** for unshorted/shorted, above 0.01. The excellent unshorted-only broad fit (0.00296) does not predict the shorted configuration: frozen prediction RMS 0.16894, or 0.14820 with only phase delay aligned. This disproves treating that particular unshorted fit as a uniquely verified extraction. It does not prove that the shorting idea is unusable; omitted fixture/distributed effects or fields changed by the short must be checked.

An additional shared input-electrode ground-shunt diagnostic only improves the joint S21 errors to 0.02277 / 0.03090. It is insufficient and has not been adopted as a verified circuit. Evidence: `joint_fits/configs_0_1/`, `heldout_predictions/short_1_from_previous_configs_bounds30/`, `input_fixture_diagnostic/`. All actual calibrated inductances remain fixed, no target-compliant geometry is claimed.

Widening the shared capacitance diagnostic bounds to +/-80% gives **0.00612 / 0.00590** S21 RMS for the same two real configurations. This succeeds as a shared-response fit, but requires markedly different capacitances (for example T3-to-T2 Cg 57.81 fF, not 153 fF), with output-fixture C0 hitting the lower bound. It does **not** establish the values as physically unique. Local profiles still permit +/-5% changes in T3 C0, input Cg and T3-to-T2 Cg while keeping both fits below 0.01. More shorted configurations and held-out predictions are needed. Files: `joint_fits/wide_shared_diagnostic/configs_0_1/S21_target_bounds_80pct/`.

Seeding the all-S fit from the good S21 fit (with several left/right delay splits) avoids poor local minima from random nominal starts. This gives S21 RMS 0.00664 / 0.00955, S11 RMS 0.04153 / 0.03353 and S22 RMS 0.00154 / 0.00266 under the same broad bounds. This is better optimization evidence, not a new physical EM result or a proof of uniqueness. `joint_fit.py` now records the additional seed source. A separate known-feedline-impedance diagnostic retains external ports at 50 ohms; it is not a port-renormalization change and is not the adopted model.

## Completed speed-stage decision

Use eight threads with the original mesh for the remaining shorted identification runs. The fully verified wall speedup is 1.740x, with identical saved full-band exports. The next two-shorted configuration began at 15:29:55 PST; the controller then runs the third and saves joint fits. `resume_extraction.py` owns these jobs, and `reflection_followup.py` adds all-S fits and held-out predictions without launching EM. The original full filter stays stopped.

**A useful rough mesh was found:** retain the original feedlines, ground, capacitors, 2 um grid and all edge refinement, increasing minimum interior subsection size to 100 cells only on coil polygons. Subsection count falls from 14,528 to 12,586. Full-band complex S21 RMS is **0.00535839**, maximum pointwise error 0.0186366. This passes the new RMS 0.01 screening allowance. The maximum is not subject to that RMS threshold. Evidence: `speed_followup/coils_min100_edges/status.json`, `abs_validation/em.log`, `S21_comparison.png`, `intermediate_meshes/coils_min100_edges/geometry_checks.json`.

The rough-mesh pilot used 130 s per frequency versus 150 s for the original-mesh eight-thread resonance-region point, a modest per-point saving. The complete pilot-plus-validation campaign consumed **692.13 s of native-runner time** across **five** directly calculated frequencies, including one extra pilot. Thus this campaign itself was not faster than the four-frequency original-mesh eight-thread sweep. Endpoint and adaptive validation stages total about 548 s, but this is a split, cache-reusing measurement, not a separately measured cold single-run ABS benchmark. Do not advertise a larger full-sweep speedup from subsection counts alone.

During validation, an initial controller changed numeric sweep limits without converting the discrete sweep type. The two endpoint calculations completed normally and were retained. `complete_local_validation.py` preserved the logs in `direct_endpoint_screen/`, converted the sweep to ABS with interpolated export, and reused the completed field cache. The final ABS calculation added 5.0975 and 4.86 GHz; no completed endpoint was recalculated. The corrected `local_edge_pilot.py` performs this conversion explicitly for future use. Full-band validation, rather than the initial pilot, determines acceptance.

The sixteen-thread same-geometry direct pilot completed in **100 s per frequency**, 109 s native total, **117.24 s actual runner elapsed**, with S21 difference 9.30e-11 versus the fresh four-thread reference. This is 2.70x faster per frequency than the 270 s four-thread pilot. Its 525.59 s wrapper timer includes waiting behind another solver and must not be reported as solve time. It is **one-point evidence only**, not a full-band sixteen-thread benchmark. Files: `mesh_benchmarks/sixteen_threads/status.json`, `em.log`, `em.log.resources.json`.

Resource samples show average host load 33.2% / 44.4% / 52.6% for the four-thread full sweep, eight-thread full sweep, and sixteen-thread pilot respectively. Other host activity was not controlled; these are measured workstation timings rather than a dedicated-machine benchmark. Sixteen threads briefly reached 86.7% total host load, below the sustained-load stop duration; eight threads provides more headroom for unattended identification. `speed_followup/resource_comparison.json` records these observations. No CPU guard interrupted these successful runs.

These mesh-comparison RMS values describe numerical agreement between Sonnet settings. They are **not circuit-fit RMS**, and do not certify extracted capacitances to 1%. Rough screening must be followed by the original-mesh parameter and response checks.
