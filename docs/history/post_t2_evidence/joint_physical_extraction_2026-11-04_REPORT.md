# Shared physical extraction and geometry-prediction investigation

**No parameter set has yet passed the requested verification.** The completed saved-data tests demonstrate that an excellent taper-section S21 fit can still give a poor full-filter prediction. The independent native interventions needed to distinguish physical parameters and test a geometry change are prepared, read-checked, and queued behind the running full sweep. They have not been simulated yet.

The objective remains an actual filter response around **4.8–5.2 GHz**. Component values are useful only to the extent they describe transferable physical behavior and tell us how to change the approved geometry.

## Completed evidence

### Geometry and circuit definitions were audited first

`results/geometry_and_settings_audit.json` verifies native metal, rather than relying on file names:

- The full filter's entire left taper/fixture region exactly matches the accepted short test, by empty Boolean XOR.
- The right taper region is an exact x reflection of the left.
- All ten regular-cell bodies match the short test's regular-cell body.
- The regular-only reference's right fixture and adjacent regular body exactly match the corresponding region of the short test after translation.
- Materials, conductor, covers, layer thicknesses, control settings, mesh and port definitions are recorded. Mesh is 2 um; external reference impedances are 50 ohms.

The short test is **T3–T2–T1–one regular resonator–a directly feed-connected output-end resonator**. Its output-end C0 is a separate fixture parameter; there is no further series coupling after that last shunt LC. The full filter instead has ten regular resonators and mirrored tapers. The custom network evaluator reproduces the existing short and full ABCD implementations to better than 1e-10 in the original topology. These implementation checks are distinct from whether that topology fits EM.

### The short-section fit alone is demonstrably insufficient

With independently calibrated L fixed and all five C0 plus five couplings free within half to twice their extraction anchors, the short-only complex S21 fit reaches **0.00286435 RMS**. It is visually excellent in transmission. Its full-filter **magnitude** prediction error is **0.273951 RMS**. Its short reflection errors are also large: complex S11 RMS 0.22443 and S22 RMS 0.13546.

The full-filter phase nuisance parameter is untrained in a short-only fit, so its large unaligned complex prediction error must not be interpreted as component error alone. The magnitude error is independent of that phase ambiguity and already rejects this parameter set as transferable.

The fitted output-end C0 reaches its upper bound. That is another reason not to interpret this low transmission residual as a unique physical measurement.

Reoptimizing all other parameters while fixing selected values at +/-20% gives these explicit alternatives:

| Fixed change from the best short-only fit | Reoptimized short complex S21 RMS |
|---|---:|
| T3 C0 -20% | 0.003179 |
| T3 C0 +20% | 0.003155 |
| T2 C0 +20% | 0.002878 |
| Input–T3 Cg +20% | 0.002884 |
| T3–T2 Cg +20% | 0.002981 |

All these listed alternative optimizations converged. Thus **fixed L plus RMS below 0.005 still does not constrain those capacitances to a few percent in this short measurement**. In contrast, the T1–regular coupling is more sensitive in this model: +/-20% produces approximately 0.026–0.031 short-fit RMS. These are deterministic feasibility/profile tests, not noise-based confidence intervals, and a sensitive parameter is not automatically an accurate one if the model is wrong.

`results/short_parameter_profiles.json` contains the completed profiles of all ten short capacitances, including convergence flags and unchanged-parameter full-filter predictions. Out-of-bound profile points are omitted. T2–T1 coupling and output-fixture C0 also admit substantial changes with low short-section residual.

### Shared fitting was tried, not just proposed

The first series fits the first 72 completed full EM calculations and holds out the next 27, with the completed short response as a second configuration. A second series uses all **99** direct full points captured at the investigation snapshot. The latter has no independent full-frequency holdout. The short dense export is sampled uniformly for efficiency, but it comes from four native direct EM fields and must not be counted as hundreds of independent measurements.

The completed search covers **20 model/objective configurations and 142 multistart fits**, including:

- short-only, full-only, and shared short/full original LC models;
- complex S21, magnitude-only, and all-S fitting objectives;
- shared positive shunt loss;
- additional positive series leads as a limited distributed-trace diagnostic;
- last-resort fitted L variation of +/-1%, +/-3%, and +/-10%, after fixed-L tests;
- distinct boundary and interior regular-cell capacitances/couplings;
- a third regular-only configuration with separately modeled left fixture.

None established a shared response fit near 0.005. For example, shared loss plus series leads with calibrated L fixed gives short/full complex S21 RMS **0.07415/0.09586**. Distinguishing boundary and bulk values improves those to **0.04770/0.08150** in the corresponding lead model, still insufficient. Allowing L to move +/-10% gives **0.07747/0.08987** in the tested shared-loss model. Adding the regular-only reference still fails the required accuracy in the tried model families. These finite searches do **not** prove a global lower bound or rule out every possible physical circuit. They do mean these resulting values must not be promoted as verified extraction.

See [COMPARISONS.md](COMPARISONS.md), `comparison_summary.json`, and per-start JSON files for every result and objective choice. The current three-configuration diagnostic uses the regular reference's saved 4.5–5.5 GHz range; the queued intervention analysis explicitly restricts that reference to the common 4.65–5.40 GHz band. Magnitude-only fits are not evidence of correct phase or reflection. Extra leads, boundary parameters, or fitted inductances are explicitly diagnostic additions, rather than silently relabeled original C0/Cg/L values. No calibration or selected design was overwritten.

## A physically relevant distinction: boundary versus bulk

Identical resonator bodies do not imply identical surrounding metal. At native y=700 um:

- T1–first-regular ground trace occupies x=3186–3188: **2 um** wide, **36 um** from the first regular electrode at x=3224.
- A bulk regular–regular ground trace occupies x=3552–3556: **4 um** wide, **26 um** from the next regular electrode at x=3582.

Consequently, transferring the short test's single regular C0 to all ten full-filter cells is an assumption about the environment. The boundary/bulk extension tests that assumption and improves some fits, but **does not by itself resolve the mismatch**. Finite ground impedance, nonlocal coupling, distributed trace/meander response, and fixture modeling remain possible contributors. Their individual causes have not been isolated. The regular-only input and output fixtures are not exact geometric mirrors, so the three-configuration test gives the left fixture independent capacitance/coupling/delay parameters instead of forcing an unsupported symmetry.

## The useful existing intervention evidence

An older T3 geometry already has explicit native capacitor responses at 0.1, 100, 250 and 400 fF, plus two fresh direct-EM checks of the 100 fF case. I reanalyzed these saved results.

At each frequency, a single ideal capacitor in a fixed linear EM environment has a common-denominator response of the form `S(C)=(A+B*C)/(1+D*C)`. Fitting the first three capacitor settings predicts the withheld 400 fF response to numerical precision (S21 RMS 2.39e-13). That extraordinarily small number is an **algebraic/cache consistency check**, not independent EM accuracy: the settings share the native EM/cache calculation. The two separately calculated native frequencies at 5.0069 and 5.1427 GHz differ from the cached 100 fF response by about **5.2e-6 in complex S21**, independently supporting native capacitor termination at those two frequencies.

The nearly-open 0.1 fF component case differs from the unprobed historical circuit by 0.00430 complex S21 RMS before alignment. This demonstrates why unchanged metal alone is not enough: adding a component's internal terminal/calibration structure also needs a check. The older geometry is not the accepted current geometry, so its coefficients are **not** transferred to the new design.

This method supplies a measured response to a known local capacitance change without forcing the EM environment into an inadequate original LC ladder. It does not establish the absolute original C0 or tell us that widening a pad is exactly equivalent to adding an ideal capacitor. The new experiments explicitly test those distinctions.

## Native verification now queued

`intervention_results/status.json` records the live worker. It is waiting for **all** of these conditions before launching any new EM:

1. The healthy current full-filter sweep finishes every frequency.
2. Its exported response passes integrity checks and the controller completes full analysis.
3. The already queued saved-data refinement worker finishes or records a failure.
4. No other `em.exe` is running.

It will then execute serially with the existing four-thread guarded runner. It never stops or restarts the current full solver.

### A. Known capacitors on the current accepted metal

Two read-checked `.sonx` files attach a native CAP across the existing top-ground gap, at the actual T2 or T3 top node. **Every metal vertex, material, external port, mesh and reference plane is preserved.** Added capacitor values are 0.1, 25, 50 and 100 fF. The first three train the intervention model; 100 fF is withheld. The nearly-open case must agree with the unprobed circuit before its extraction can be used. Two fresh direct frequencies additionally check the 100 fF native calculation independently of cache interpolation.

Normal native ABS is allowed to add fields when capacitor changes move resonances; an unverified cache-only interpolation is not forced. Cache reuse is a possible efficiency, not an assumed measured speed improvement. The numeric capacitor responses come from Sonnet, not from the earlier sparse four-port approximation.

### B. A real geometry change, not another circuit fit

The T3 pads are prepared at widths **176, 184 (baseline), and 192 um**. Their tops and outer edges remain fixed; the inner edges move by +/-8 um. The accepted topology, flat top ground, clearances to other components, regular metal, T1 and T2 remain intact.

Only **four existing vertices** move in each short case, or **eight** across both mirrored full-filter T3 cells. Original polygon partitions, IDs, remote vertices and port indices are retained, avoiding an unnecessary whole-layout reimport/mesh repartition. Boolean XOR checks confine every change to the intended pad strips. All four projects pass Sonnet's native read check, and the before/after geometry was visually inspected.

Short geometry cases receive a completed sweep. Full geometry cases initially use six **direct** baseline-matched frequencies: **4.8292, 4.9305, 4.9776, 5.0248, 5.0758 and 5.1818 GHz**. These test the response-change direction without another hundred-frequency adaptive run for each probe. They do not certify the whole passband.

### C. What will count as evidence

`analyze_interventions.py` will:

- fit one shared underlying physical set across the short section, regular-only reference, full filter and training capacitor loads;
- predict both withheld 100 fF responses without changing those physical values;
- check complex S11/S21/S22, not merely transmission magnitude;
- perform reoptimized +/-5% profiles if response gates pass;
- test whether each actual short pad change is well represented by the measured local capacitor sensitivity;
- compare its predicted full-filter response change with the six independently simulated full geometry points;
- report actual improvement toward the original target at those points separately from extraction error.

The working response gates are S21 complex RMS below 0.005 and all-S RMS below 0.01 on each configuration. Those are investigation criteria, not a guarantee of absolute physical accuracy. Passing them still requires reviewing profile spread, calibration consistency and geometry-change prediction. Nothing is automatically promoted to the production design or declared uniquely measured.

If a gate fails, the worker preserves the evidence and reports the failure rather than quietly loosening definitions or declaring the least-bad fit correct. A geometry step will only have a verified tuning direction after its measured response agrees with the prediction. Absolute certainty is not warranted by finite simulations; the practical requirement is a quantified prediction that survives independent interventions.

## Files

- `campaign.py`, `three_config.py`: completed multiconfiguration model searches.
- `audit.py`: exact native metal/environment/fixture checks.
- `profile_short.py`: reoptimized short-parameter ambiguity tests.
- `check_historical_loads.py`: historical native capacitor audit and common-denominator operator check.
- `prepare_known_loads.py`, `prepared/manifest.json`: current-geometry known-load files, **prepared, not simulated**.
- `prepare_pad_interventions.py`, `prepared/geometry_interventions.json`: small physical pad interventions, **prepared, not simulated**.
- `run_after_full.py`: active completion-gated serial study worker.
- `analyze_interventions.py`: subsequent parameter and geometry-prediction validation, not yet executed on new EM data.
- `results/`, `results_all_direct/`, `results_three_config/`: frozen data, fits, profiles, plots and provenance.
- `intervention_results/status.json`: live next-stage state; `REPORT.md` will be generated there after the new EM study completes.

No original Markdown, project, calibration, published selected design, or running native full sweep was altered by this investigation.
