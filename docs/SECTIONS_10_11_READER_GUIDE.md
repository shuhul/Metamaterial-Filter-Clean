# Reader guide: the assembled taper and the full filter

This guide explains Sections 10 and 11 of the Filter Automation report. It starts with the purpose, then follows the actual calculation from Sonnet geometry to the accepted taper and the full-filter test. It separates measurements made from simulations, numerical fits, design comparisons, and what remains unverified. No new simulation was run to write this guide.

## 1. What we are trying to build

The verified project goal from the design discussions is a manufacturable superconducting microwave filter with a passband around **4.8–5.2 GHz**, connected to 50-ohm ports. Ten repeated regular resonators form the middle. Three taper resonators on each side transition between the ports and the regular chain:

`port — T3 — T2 — T1 — 10 regular cells — T1 — T2 — T3 — port`

There are sixteen resonators in the full filter. The taper changes the resonator and coupling values near the ends to reduce mismatch. Good matching permits energy to enter the regular array instead of reflecting back toward the source. The regular-cell electrical scale puts the desired collective band near 5 GHz, but the coupled array's band edges cannot be calculated from an isolated LC frequency alone.

The user supplied the target circuit values and Frank's paper as a reference, and repeatedly identified the 4.8–5.2 GHz filter response as the practical objective. The material available here does **not** establish why that band was selected for a particular qubit, transition, readout channel, experiment, or fabrication package. Those possible applications must not be invented. Likewise, the 50-ohm port impedance is the wave normalization/termination, not a promise that every piece of physical transmission line has exactly 50-ohm characteristic impedance.

The supplied AWR inputs are preserved unchanged in `reference_inputs/awr/`. The screenshot shows an enabled goal for the magnitude of Full_Filter S(1,1), measurement less than 0.001, from **4.82 to 5.18 GHz**, with Difference cost, default L=2 and weight 1. That is a reflection objective in the interior of the requested band, leaving 20 MHz at each band edge. It is not a screenshot of an achieved EM response. The user reports that minimizing S11 worked better than directly maximizing S21 in their AWR optimization; that experience is recorded without claiming universal optimizer superiority or that this goal was already used by our Sonnet geometry search.

For a lossless two-port with the second port matched, energy conservation gives

`|S11|² + |S21|² = 1`.

Therefore, in that ideal case an S11 magnitude below 0.001 corresponds to reflected power below one part in a million and S21 magnitude above about 0.9999995. A real EM model can also dissipate or radiate power, so low reflection alone is not sufficient: a matched absorber can have low S11 and low S21. We must still inspect transmission and energy loss. The threshold is a linear magnitude; 0.001 corresponds to −60 dB, not −0.001 dB.

## 2. What Section 10 actually simulates

The accepted taper fixture is smaller than the full filter. It contains T3, T2, T1, a regular resonator, and an output-end resonator/fixture connected to the second port. Including nearby metal captures the loading and coupling that would be missed by isolated component tests. This fixture is useful for faster design iterations, but its two-port S21 is not the full sixteen-resonator filter's S21.

The canonical input is `final_designs/final_T3/project.sonx`; its completed response is `project.snp`. The continuous-ground drawing is `continuous_ground/project.sonx`. The original and merged drawings have identical metal coverage, but the merged taper drawing is not a separate new EM validation. The full-filter geometry uses merged ground and has its own EM calculation.

Open `s21_with_target_parameters.png` to see the native response and the bounded circuit comparison. Open `s21_verified.png` to see the independent network prediction as well. These are three distinct curves: the geometry simulated by Sonnet, a model of that geometry, and a nearby target circuit.

## 3. What Sonnet and S-parameters tell us

Sonnet solves the electromagnetic problem for the metal, dielectric stack, enclosure/boundaries and ports. An S-parameter describes outgoing waves produced by an incoming wave, with other ports terminated in their reference impedances. S21 is transmission from port 1 to port 2; S11 is reflection at port 1. These are complex numbers: magnitude measures amplitude, phase records timing/phase shift.

Two external ports provide only a limited view of many coupled internal structures. Several different allocations of capacitance and inductance can reproduce the same external response. An excellent fit is useful, but it does not by itself tell us which individual physical electrode has a uniquely measured C0 or Cg.

Before this method, the regular/T1/T2 workflow calibrated L with known capacitors and fitted a small LC ladder. For the standalone T2 selection, three new capacitances and a fixture delay varied while older values and L were fixed. This is a conditional extraction: the returned T2 values depend on the fixed downstream assumptions. Those historical selections remain under `final_regular`, `final_T1` and `final_T2`; the later embedded T1/T2 geometry in `final_T3` is not identical to every standalone baseline.

## 4. Why add internal ports, and what does “shorting” mean here?

The current method adds five internal GAP ports at existing coil-return connections, retaining the two external ports. This produces a seven-port network. Polygon splitting exposes the connection while preserving the total metal union, checked by a Boolean XOR test. The ports give independent electrical access to the internal connections rather than attempting to infer everything from external transmission alone.

These ports are **not automatically nonloading voltage probes**. Looking at S21 with all extra ports terminated in 50 ohms would describe a different circuit. Instead, the saved multiport response is reduced by numerically restoring the internal connections as shorts.

This differs from the earlier physical bridges across series capacitors. Those bridges changed the connection topology and left the earlier shunt inductors present. They did not remove unwanted resonators. Here the internal-port network is an instrumented representation, and its restoration is tested against the original uninstrumented two-port calculation.

The mathematical reduction is straightforward. Separate the impedance matrix into external (`e`) and internal (`i`) blocks:

`Ve = Zee Ie + Zei Ii`

`Vi = Zie Ie + Zii Ii`.

Restoring ideal shorts sets `Vi=0`, so `Ii = −inv(Zii) Zie Ie`. Substitution gives

`Zrestored = Zee − Zei inv(Zii) Zie`.

This is the Schur complement: it accounts for internal currents that flow when the connection is closed. With equal real 50-ohm references,

`Z = 50 (I+S) inv(I−S)` and `S = (Z−50I) inv(Z+50I)`.

The implementation also checks the equivalent direct S-matrix reduction and applies the recorded internal-port polarity convention. It does not infer physical correctness merely from unchanged metal. For the accepted geometry, restored-versus-original maximum error across all four external S entries is about **0.00011846**, below its 0.001 restoration gate. See `final_T3/local_model/analysis.json` and `manifest.json`.

## 5. What is extracted from the seven-port calculation?

There are two related representations.

**A low-frequency terminal capacitance matrix.** At 0.5, 1 and 2 GHz, the selected six-terminal imaginary impedance is represented as

`Im Z(r) = −A/r + B r + D r³`, where `r=f/(1 GHz)`.

The inverse-frequency term is capacitive. In the script's normalization, the corresponding matrix in fF is `C = inv(A) × 10^6/(2π)`. A separate 1.5 GHz point is held out, and a two-frequency extraction checks sensitivity to the extrapolation. The accepted held-out relative impedance error is 8.25e−6.

Each diagonal matrix entry includes connections from a terminal to ground and to other terminals. Negative off-diagonal entries represent mutual capacitances. The row sum is the capacitance to the represented ground. These are definitions in the terminal model, not names that can be substituted directly into the original ideal ladder.

**A frequency-dependent operating-band network.** The entire 7×7 complex Z matrix is interpolated using powers `r^-1, r, r³, r⁵, r⁷`, now with `r=f/(5 GHz)`. Training frequencies are 4.65, 4.861, 5.0, 5.25 and 5.4 GHz. The separate 5.096 GHz native calculation checks a point that did not determine the coefficients. Its restored S21 error is **1.45418e−6**. Across the dense two-port export, S21 RMS is **8.11746e−6** and all-S RMS is **3.44610e−5**.

This operating-band model retains distributed effects that a single static LC per cell can miss. It is validated for the stated geometry and band; it is not asserted to be a globally passive/causal model for all frequencies. The 1,501 samples in the comparison are dense exported values, not 1,501 independent native EM calculations.

## 6. A real example: terminal capacitance versus ladder C0

The accepted six-terminal matrix lists terminals as input, T3, T2, T1, regular and output end. In the T3 row, its diagonal is about 644.0767 fF. Couplings to the other represented terminals have magnitudes about 202.0960, 102.3464, 0.3601, 0.0606 and 0.0155 fF. Thus the ground row sum is

`644.0767 − 202.0960 − 102.3464 − 0.3601 − 0.0606 − 0.0155 ≈ 339.1981 fF`.

Yet the normal ladder comparator has T3 C0=**251.6275 fF**, exactly 0.25% above its 251 fF target. These are not inconsistent measurements of one physical quantity. They belong to different decompositions of the response: a coupled terminal matrix with distributed background, versus a restricted ladder with a single shunt C0 at each node.

Calling 251.6275 fF “the independently measured capacitance of the T3 pad” would be incorrect. It is an optimized comparator value. Calling 339.1981 fF the original circuit C0 would also be incorrect. This distinction is the reason for keeping both representations and their evidence.

## 7. How the target ladder is fitted

For a series coupling capacitor, `Zg=1/(jωCg)`. For a shunt tank, `Y=jωC0+1/(jωL0)`. Their ABCD matrices are `[1,Zg;0,1]` and `[1,0;Y,1]`. Multiplying them along the chain gives entries A,B,C,D. Between equal 50-ohm ports,

`S21 = 2/(A + B/50 + 50C + D)`.

The final comparison permits ten capacitances and five inductances within **±0.25% of their nominal targets**. It also permits a total reference delay within ±100 ps. The fitted delay multiplies S21 by `exp(−j2πfτ)`; it changes phase, not magnitude or resonance locations. The accepted delay is 3.58105 ps.

The residual is complex S21 error with equal frequency weighting:

`R = sqrt(mean(|S21_EM − S21_ladder exp(−j2πfτ)|²))`.

The accepted R is **0.00249186334**, just below 0.0025. This states that one tightly bounded normal circuit matches the taper fixture's response. It does not establish uniquely extracted physical L/C values. Final L values are comparison parameters, not fresh explicit-capacitor measurements of every final meander. Parameter bounds are design constraints, not confidence intervals.

The terminal-model RMS of 8.12e−6 and this ladder-comparison RMS of 0.002492 measure different things. A tiny modeling error alone does not make a geometry meet the target. A good target comparison alone does not prove that its individual fitted components are uniquely physical.

## 8. How we decide which geometry to change

The useful step is measuring how the complete network changes with a physical dimension. We alter a pad width, rail spacing, finger span or coil parameter, run the corresponding native model, and estimate a local derivative `dZ/ddk`. For a small update:

`Znew ≈ Zbase + Σ Δdk (dZ/ddk)`.

This includes changes in distributed loading as well as static capacitance. It avoids assuming that widening one pad changes only one ideal C0. The derivative model is a local approximation; it is not a guarantee for arbitrarily large changes.

The practical iteration is:

1. Verify a native baseline and its restored terminal network.
2. Measure controlled geometry changes and construct local derivatives.
3. Search the permitted geometry step range and bounded ladder values together.
4. Round/enumerate changes on the 2 µm manufacturing/simulation grid.
5. Save the predicted response before the candidate's EM run.
6. Compare the completed candidate response with the prediction and target.
7. Accept actual improvement; otherwise reduce the step range or remeasure the derivatives.

Sparse difference models have their own held-out/dense response gates. The recorded 0.0003 all-S threshold is an error in the **49 complex scattering entries**, not an error in 49 fitted physical parameters. Full baseline-network validation and restoration use separate gates.

A concrete final step is recorded in the cover-converged history. At pad-width increment +32 µm the actual target RMS was **0.00250163429**, narrowly failing. Increasing that existing pad-width control by another **2 µm**, to +34 µm, gave a completed native result of **0.00249186334** and passed. That is a real measured geometry result, not just a circuit prediction. It does not prove that every further width increment would improve it.

The opposite example is informative too: a regular-pad proposal predicted roughly 0.003164 RMS but produced roughly 0.008012 natively. It was not the retained solution. These cases explain why simulation verifies the geometry proposal rather than the optimizer's prediction being treated as a result.

## 9. Why the upper air region matters

The physical layout is planar, but electric and magnetic fields extend into the substrate and above the metal. The simulation therefore needs a vertical environment and a boundary condition at its top. The native model's top cover is labeled `FREE_SPACE`; it is not an added metal ground plate. Its finite distance can nevertheless influence the calculated response.

The original upper air height was 1.5 mm. With unchanged metal, moving it to 12 mm changed complex S21 by RMS **0.00966824**, and reduced maximum port-1 power deficit `1−|S11|²−|S21|²` from **0.177283** to **0.00104456**. That was large relative to the desired design error. Refinement therefore used the farther boundary.

Increasing from 12 to 30 mm on the accepted taper changed S21 by RMS **3.14856e−9**. This supports top-height convergence for these settings. It does not identify exactly how much of the earlier deficit was physical radiation, dielectric loss, or interaction with the finite numerical boundary. It does not establish convergence with respect to every mesh or enclosure parameter.

The retained stack has 12 mm upper air, 500 µm silicon, 5 mm lower air, and a bottom perfect-conductor cover. A fabricated package or nearby metal lid could change the real device. A future package model should use its actual surroundings instead of assuming this boundary test proves every enclosure is equivalent. The native surface-model `Ls=0.0388` is preserved; its field-specific unit convention remains unconfirmed here.

## 10. What Section 11 verifies

The full filter replaces the small fixture with ten regular cells and mirrored three-stage tapers. Its conditional circuit prediction repeats the regular values and mirrors the taper components, omitting the small fixture's output-end resonator and phase delay. That prediction is not a full-filter EM calculation.

Before committing to a long sweep, we solved one full-filter frequency, **4.92 GHz**, where the circuit prediction could be checked directly:

| Quantity | Circuit prediction | Full-filter EM |
|---|---:|---:|
| Linear magnitude | 0.957403974 | 0.972139673 |
| 20 log10(magnitude) | −0.378095480 dB | −0.245426656 dB |

EM is higher by **0.014735699** in magnitude, or **0.132668825 dB**. The calculation took 867.274 seconds. Phase conventions were not aligned, so this comparison is deliberately magnitude-only. It justified proceeding, but says nothing by itself about all ripple extrema or the band edges.

The full sweep uses 4.5–5.5 GHz ABS, a 2 µm grid, medium-speed edge meshing, double precision and de-embedding. Both ports are 50 ohms with 400 µm reference-plane shifts. The run uses 36 threads, selected by its launcher; the single downloadable project does not itself guarantee that runtime thread policy on another machine. ABS chooses direct frequencies adaptively, so its output interval target is not a known count of expensive solves. The supplied native project requests Touchstone S/MA output as `project.snp`, with `UseInterpolatedData=FALSE`; do not mistake a partial export for a converged full-band result.

**Final full-filter results placeholder:** add the completed native response, convergence/completion evidence, and passband assessment only after the sweep finishes and validation succeeds. Live progress belongs in a separate timestamped status record, not in this stable results section. No completion or passband-wide success is claimed here.

## 11. Files and reproduction

Current source tree:

```text
geometry/                 geometry, import, stage helpers, frozen-layout replay
tools/                    publication checks and saved-response replay
docs/                     this guide, technical method, history and audits
final_designs/
    final_regular/        historical standalone regular selection
    final_T1/             historical standalone T1 selection
    final_T2/             historical standalone T2 selection
    final_T3/             current complete taper fixture and response
    final_full_filter/    current full-band input
        verification_4p92GHz/  completed independent single-point check
final_old/                superseded publication copies and old indexes
reference_inputs/awr/     supplied target circuit and optimizer screenshot
```

The project root also has `pyproject.toml`, `uv.lock`, `reference_taper_geometry.png` and `reference_passband_plot.png`. Local experiment directories `investigations/`, `sonnet_files/` and `output_data/` are intentionally excluded from Git. Moving the curated publication folders did not move the active solver or its files.

From the repository root:

```powershell
uv sync --locked
uv run python tools/check_published_files.py
uv run python -m geometry.rebuild_final --verify-only
uv run python tools/verify_current_response.py
```

The first verifier checks artifact hashes. The second checks native/GDS metal equality, the grid and the ports for each selected input. The third replays the saved terminal coefficients and normal comparator against the accepted native response, without EM. It should reproduce 0.0024918633384848348 and 8.117461191373584e−6 for the two RMS definitions above.

For native frozen-layout regeneration:

```powershell
uv run python -m geometry.rebuild_final --output-dir build/final_designs
```

This needs Sonnet and its license for GDS import/read. It does not launch EM. It is exact frozen-layout replay, not a promise that every old exploratory optimizer is portable. The complete `.sonx` is already sufficient to inspect and simulate the frozen design; GDS, Python scripts and old solver caches are not needed merely to open it. A new simulation needs a writable output location and a working Sonnet installation/license.

Technical source snapshots are in `docs/history/post_t2_evidence/`: `reachable_network.py.txt`, `reachable_iterate.py.txt`, `terminal_profiles.py.txt`, `quarter_target_optimizer.py.txt` and related files. Their original local sources are under `investigations/filter_refinement_2026-11-05/`. The portable checker is `tools/verify_current_response.py`, and model inputs are `final_T3/local_model/network.json`, `project.snp` and `result.json`. `PUBLICATION_PATHS.json` records the directory changes; `SOURCE_PATHS.json` traces copied results to working originals.

The dependable conclusion is that we have a response-qualified taper and an encouraging single-point full-filter test. Unique original-component identification, finer mesh convergence, the role of an actual fabricated package, and completed full-band performance remain separate questions.
