# Connecting the full-filter response model to physical design

The unchanged 72-point passive model predicts all 26 subsequent completed Sonnet calculations with **complex S21 RMS 0.00383870** and magnitude RMS **0.00174428**. This validates useful forward prediction for the existing metal. It does not supply a derivative with respect to metal dimensions. The actual filter still has substantial passband ripple and does not yet reproduce the original target.

The native adaptive full-filter job was left running throughout. No native geometry, calibration, selected design, solver setting, or original rational coefficient was changed. All new outputs are in this investigation directory. The current evidence contains 98 completed direct frequencies; unfinished dense native exports were not used.

## 1. New-point validation

`frozen_model_validation.json` contains every held-out point and the unchanged coefficient hash. The script checked that the original 72 training rows still match the frozen training snapshot, including complex S parameters.

| Quantity | Result |
|---|---:|
| Original training complex S21 RMS | 0.00134715 |
| New points, never refitted | 26 |
| New-point complex S21 RMS | 0.00383870 |
| New-point magnitude S21 RMS | 0.00174428 |
| Worst new complex S21 error | 0.0123649 at 4.8252 GHz |
| New-point RMS across all four S entries | 0.00753368 |

`unseen_validation_and_target.png` was visually inspected. The new circles closely track the frozen response prediction. The largest errors cluster around the sharp lower skirt. The model predicts the existing passband ripple accurately; the orange target is substantially flatter. The dense lines are model predictions, not completed dense EM. Agreement with this model does not establish mesh convergence or correct nominal filter behavior. Deep-stopband relative/dB accuracy is poorer than absolute complex accuracy.

## 2. A concrete explanation of nonuniqueness

This is more than an optimizer being unreliable. Even the original ideal LC ladder admits an exact internal voltage-scaling ambiguity when its inductances are free.

In nodal coordinates, write `Y(s) = s C + K/s`, where `K = inverse(L)` and the off-diagonal entries of `C` are minus the coupling capacitances. For a diagonal matrix `D` whose terminal entries equal one, replace both matrices by `C' = D C D` and `K' = D K D`. Eliminating the internal nodes gives exactly the same terminal admittance, because the internal scaling cancels in its Schur complement. This transformation can retain positive physical capacitors and inductors.

For the current reference ladder, scaling both hidden T2 nodes by 1.1 gives this example:

| Parameter | Original | Alternative |
|---|---:|---:|
| T2 L (nH) | 4.34479 | 3.59073 |
| T2 C0 (fF) | 141.550 | 188.863 |
| T3 C0 (fF) | 238.450 | 223.915 |
| T1 C0 (fF) | 149.044 | 147.591 |
| T3–T2 Cg (fF) | 145.350 | 159.885 |
| T2–T1 Cg (fF) | 14.5350 | 15.9885 |

All other values remain unchanged. The maximum difference across **all four S parameters at 10,001 frequencies is 4.12e-14**. Neither more terminal frequencies nor adding S11 to the same unmodified circuit distinguishes these two parameter sets. `exact_nonuniqueness.json` preserves the complete values and numerical check. This is an algebraic example, not two new EM geometries or an assertion that the alternative matches independent L calibration.

**Fixing independently calibrated L removes this particular ambiguity.** With L fixed, the local capacitance Jacobian has all nine independent columns at the tested reference. Its condition number is still approximately **24,725**, so sensitivity is very uneven. Full local rank is not a guarantee of accurate extraction in the presence of a wrong circuit model, finite numerical error, or other distant minima. The arbitrary pole/residue macromodel has even less of a unique mapping to named physical elements.

## 3. Ways to obtain physical information that is missing

### Known internal loads: retain the explicit-capacitor approach

A controlled internal load changes the information content, unlike another terminal-frequency point on unchanged metal. I tested this on the original circuit equations, using T2 shunt additions of 0, 50, and 100 fF and keeping the same underlying L/C values across the three responses. The joint complex-S21 Jacobian has rank 13 for the 13 shared L/C parameters at the tested anchor, versus 11 at the relative 1e-6 sensitivity threshold for the unloaded data. The two exactly equivalent example circuits above produce **0.16766 complex S21 RMS difference** when both receive the same known 100 fF load. The load breaks their equivalence.

This is a **synthetic identifiability test**, not a newly completed loaded Sonnet experiment, and not a global uniqueness proof. The practical version is the previously used explicit Sonnet capacitor calibration: identical metal, common component locations, identical mesh/material/port definitions, several known added capacitors, and a joint fit sharing the original physical parameters. The added capacitor must not be allowed to refit away. Its connection and mode assignment must be verified. These experiments cannot repair an incorrect distributed-to-lumped model by themselves.

### Internal admittance: useful only with explicit reference definitions

If calibrated internal observations are needed, physical node access can constrain a local admittance model. After accounting for ports, ground references and feed sections, the lumped approximation gives `C_ij = -Cg_ij` off-diagonal and `C0_i = sum_j C_ij` for a complete set of represented signal nodes. The inductive coefficient is `K`; it is diagonal only when mutual inductance is negligible. Fit several directly simulated frequencies and check held-out frequencies and port-attachment effects.

The external two-port admittance is a Schur complement over hidden nodes; its entries are not the original local C0/Cg values. Simply converting external S to Y does not solve this problem. Sonnet documents calibrated internal ports and warns that port structures introduce parasitics requiring calibration: [Introduction to Ports](https://www.sonnetsoftware.com/support/doc_19.52/Topics/Ports/Introduction_to_Ports.htm), [Co-calibrated Internal Ports](https://www.sonnetsoftware.com/support/doc_19.52/Topics/Ports/Co-calibrated_Internal_Ports.htm), [De-embedding Overview](https://www.sonnetsoftware.com/support/doc_19.52/Topics/Analysis_Setup/De-embedding/De-embedding_Overview.htm). This is an optional follow-up, not a new four-port/sparse-interpolation calibration adopted without validation.

## 4. What is wrong with the present filter response?

At the worst completed sample between 4.85 and 5.15 GHz, **5.0248 GHz**, the directly computed power balance is:

- transmitted: **28.07%**;
- reflected: **69.68%**;
- missing from the two ports: **2.25%**.

Thus reflection dominates this particular central dip. It is productive to investigate taper matching and loading before assuming dissipation explains the ripple. Other frequencies have different power balances; this does not rule out radiation, distributed modes, or loss elsewhere. The missing-power mechanism has not been isolated.

The frozen response estimates central dips around -4.6 to -5.9 dB, while the target is almost flat centrally. It also predicts a lower useful edge near 4.824 GHz, versus the target's -3 dB edge near 4.803 GHz. Its uniform-band magnitude discrepancy from the target is **0.25446** over 4.65–5.40 GHz. This measures the difference in shape; it is not the much smaller error of the model against EM.

The nominal infinite regular-chain dispersion is:

`f_low = 1 / (2*pi*sqrt(L*(C0+4*Cg)))`

`f_high = 1 / (2*pi*sqrt(L*C0))`

The requested regular-cell targets give **4.79933–5.20033 GHz**, confirming that they are consistent with the intended band. An archived regular-only extraction gives an ideal-chain estimate of **4.81607–5.21947 GHz**. That is a plausible lead for the displaced band, not a measured Bloch dispersion or proof that regular cells cause the full-filter shift. The historical extraction itself has model error. It also illustrates why sub-percent component errors can shift a sharp edge by tens of MHz. Regular metal remains unchanged, as requested. If its actual bulk band is displaced, taper matching alone cannot correct the infinite-chain dispersion; that question should be measured rather than hidden by a taper refit.

## 5. What the geometry history actually establishes

`geometry_evidence.py` found **32 completed usable T3 assemblies** and **9 pairs differing in one recorded numeric T3 setting**, with the same recorded downstream reference. The completed raw responses, rather than their fitted capacitor labels, were compared on a uniform 4.7–5.3 GHz grid. These are small-assembly comparisons; native material/port fingerprints and transfer to the current cleanup/full assembly have not been independently re-audited. No family provided three usable single-control values for a held-out interpolation check.

Examples of measured small-assembly effects:

| Recorded change | Magnitude RMS to small target: before → after |
|---|---:|
| Input fingers 28 → 40, original T2 junction family | 0.11101 → 0.07199 |
| Input fingers 40 → 70, another junction family | 0.06595 → 0.08733 |
| T3–T2 fingers 22 → 40 | 0.06523 → 0.05328 |
| Pad width 212 → 220 um, same other settings | 0.05309 → 0.05328 |

This supports using input/junction coupling as influential controls. It does **not** establish that increasing either coupling improves the current full filter. The input response can be undercoupled or overcoupled depending on the surrounding structure. Fitting each trial into a different degenerate capacitor set is not a reliable substitute for comparing its actual S response.

The ideal-target circuit sensitivity also shows appreciable response to T1 C0 and the T1–regular junction. It would be unjustified to assume only the outer T3 IDC can repair every ripple. Those are circuit-only directions, not measured geometry gradients. Preserve the regular metal and change approved taper dimensions locally; no new ground notches, new circuit topology, or speculative extra metal is implied.

## 6. Recommended route to an actual working geometry

Use **forward geometry sensitivity** to guide the changes, with the original LC model retained as an interpretable guide and independent calibration check. Unique extraction is helpful, but it is not a prerequisite for measuring how a dimension affects S21.

1. Allow the current full sweep to finish and pass response integrity checks. Keep the frozen model and validate additional points against it before refitting. The accurate rational model is a baseline response/interpolation tool only.
2. Preserve the approved geometry family, flat top ground, regular metal, centered components, and 2 um grid. Keep measured L as a physical constraint. Start a local experiment around this exact accepted cleanup geometry, rather than transferring a large derivative from a different historic layout.
3. Probe a small number of independent taper controls: input IDC overlap, T3–T2 IDC overlap, T3 pad extent, and (if the matching response calls for them) T1 C0 pad extent and T1–regular gap/overlap. Use one changed control per comparison, with modest positive/negative, grid-resolved steps and unchanged remaining settings. Increasing overlap usually increases intended coupling but also changes parasitics; measure the entire S response.
4. Use the shorter assembly for fast screening, but validate the proposed direction on the full filter at informative passband troughs and both skirts before relying on the shorter model. It has a different termination/loading, so small-circuit success alone is insufficient. Choose actual direct frequencies from the completed full-filter valleys/edges; confirm shift-sensitive points with nearby frequencies rather than trusting one point. Do not start another large adaptive sweep for every trial.
5. Form the measured finite-difference response Jacobian with respect to those physical dimensions. Compute a bounded, regularized step that improves magnitude against the original target on a uniform frequency grid. Separately report band edges, worst central transmission, ripple, reflection, and stopband rejection. Averages alone can hide a narrow failed region. Do not introduce an arbitrary new ripple/rejection acceptance threshold without a specification.
6. Simulate the candidate, compare the actual response change against the predicted change, and shrink the step/rebuild the derivative if prediction fails. A predicted low objective is not an achieved geometry. Finish with the complete mirrored 10-cell full-filter verification.

`../geometry_response_step.py` implements the bounded local step from **measured** plus/minus EM responses. It does not claim the required baseline derivatives already exist. It neither invokes EM nor edits geometry; only the command-line/input path was exercised here. New native controlled perturbations and their validation remain necessary. `geometry_step_manifest.example.json` shows its input contract and deliberately fails its setup-verification gate until real evidence is supplied.

This route avoids requiring an arbitrary equivalent-circuit inversion to be unique before choosing a metal change. It also prevents the opposite mistake of treating a good interpolation fit as a geometry optimizer.

## 7. Additional conditional fit check and limits

`fixed_L_cap_fits.json` records three-start original-ladder fits of the 98 points with independently calibrated L fixed and all nine shared capacitances bounded within +/-20% of the extraction anchor. The best complex-S21-only trial reached 0.14964; fitting all S entries reached S21 RMS 0.21698 and all-S RMS 0.21208. Some starts hit evaluation limits or different local minima. These limited searches do **not** establish a global error floor, nor prove that no equivalent circuit could fit. They supply no credible new physical extraction and were not promoted. The earlier broader circuit-family results remain separate.

The exact nonuniqueness example, frozen-model validation, and direct power audit do not depend on these optimizers. The new practical recommendation rests on those tests and measured response differences, not on interpreting failed/local fits as uniquely measured capacitances.

## Files and reproducibility

- `../validate_and_design.py`: frozen-model validation, exact equivalence, loaded synthetic Jacobian, target sensitivity, dispersion and direct power audit.
- `../geometry_evidence.py`: completed historical single-setting response comparisons.
- `../fixed_L_check.py`: the limited conditional fitting diagnostic.
- `../geometry_response_step.py`: future measured-geometry response step; no solver or geometry actions.
- `direct_snapshot.npz`, `em_snapshot.log`: the 98-point read-only evidence snapshot.
- `frozen_model_validation.json`, `unseen_validation_and_target.png`: independent prediction check.
- `exact_nonuniqueness.json`, `local_identifiability.json`, `known_load_identifiability.json`: physical identifiability evidence, with scope labels.
- `direct_power_audit.json`, `filter_performance.json`, `regular_dispersion.json`, `regular_band_hypothesis.json`: mismatch diagnostics and their limits.
- `historical_geometry_evidence.json`: source paths and before/after raw-response metrics.
- `ideal_circuit_sensitivity.json`: circuit-only parameter perturbation guide.

Replay with `.venv/Scripts/python.exe -B investigations/full_filter_fit_72_2026-11-03/validate_and_design.py` and the analogous `geometry_evidence.py` command. The validation rerun takes a fresh completed-direct-data snapshot but never refits the frozen coefficients. Existing snapshot/output files in this new investigation subfolder will be refreshed; production inputs and original calibration evidence remain immutable.
