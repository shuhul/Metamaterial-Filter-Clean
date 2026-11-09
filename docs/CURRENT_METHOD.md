# Current taper method and evidence (November 8, 2026)

The accepted taper matches a normal ladder comparator with all ten capacitances and five inductances within 0.25% of nominal, at complex S21 RMS 0.00249186334 over 4.65–5.40 GHz. These are **optimized comparator values, not uniquely extracted physical L/C values**. An independently simulated terminal network predicts the taper response at RMS 0.00000811746. These two errors answer different questions.

Current artifacts are in `final_designs/current_2026_11_08/`. The older `taper_3/` and `full_filter/` packages remain historical. `SOURCE_PATHS.json` maps every copied original artifact to its new repository path and hash. No working simulation was moved.

## What changed after the selected T2 stage

The earlier T2 procedure fixed calibrated inductances and adopted downstream capacitances, fitting T2 C0, input Cg, and T2–T1 Cg against complex S21 with a reference delay. That was conditional extraction: fixing older values removed fitting freedoms, but did not demonstrate that those values remained physically correct in the assembled layout. A small two-port residual did not establish uniqueness or transfer to the full filter.

The first accepted T3 selection (`geometry/taper/selected_taper3.json`) used width 130 µm, 68 horizontal runs, wire 6 µm, gap 2 µm, pads 184 × 594 µm. It recorded RMS 0.0200526, C0 256.1465 fF, input Cg 263.0294 fF, L 2.7923246 nH, and junction Cg 156.6004 fF. Older capacitances could move 3%; inductances were fixed. This is a historical result, not the current 0.25% checkpoint. Staged 1/2/3% refits, retained-junction tests, wider-band and loaded-probe diagnostics are implemented in `geometry/taper/taper3_*.py`; selection files and their source paths preserve the actual decisions.

Subsequent work corrected IDC rail ends and centered the T2 coil/pads. Full-filter fits and joint short-fixture fits then showed why matching one response was insufficient. In the saved shorting investigation, an unshorted-only fit of 0.00296 predicted the first shorted fixture with RMS 0.16894 (0.14820 with delay alignment). Joint broad fits gave 0.02566/0.03317. Widening capacitance bounds to 80% reduced both errors to 0.00612/0.00590, but required substantially different parameters and still admitted compensation. Shorting series couplings retains all shunt inductors and changes local fields; it does not isolate a lone unit cell.

Those numbers are from `investigations/shorting_speed_2026-11-05/REPORT.md`, preserved in the mailbox evidence package. They concern earlier geometries/settings and must not be pooled with current acceptance statistics. The eight-family full-filter fit and joint physical-extraction investigations likewise establish model flexibility, not unique original circuit parameters. Their reports retain failed transfer checks.

## Terminal network construction and restoration

The retained method in `reachable_network.py` adds five internal GAP ports at existing coil return connections, retaining the two external ports. It splits polygons at a shared edge, checks empty XOR of the total metal union, and verifies native reading. It does not treat extra 50 Ω terminations as nonloading probes. The internal connections are numerically short-restored before comparison with the original two-port EM result.

With external/internal partitions e/i, the short-restored impedance is

`Zeff = Zee - Zei inv(Zii) Zie`.

For equal real 50 Ω references, `Z = 50 (I+S) inv(I-S)` and `S = (Z-50I) inv(Z+50I)`. The equivalent direct scattering closure is `See - Sei inv(I+Sii) Sie`. The code checks these two closures against each other and uses a polarity transform `diag(1,1,-1,-1,-1,-1,-1)` before terminal analysis. This is a validated series-port restoration, not an assertion that internal ports are intrinsically nonperturbing or a general known-capacitor calibration result.

For the accepted geometry, restored-port versus original native data has maximum all-S error 0.000118456, below the 0.001 reconstruction gate. See `taper/local_model/analysis.json`. The source project, seven-port response, manifest, completion record and network coefficients are included there. No additional known-capacitor connection validation is claimed by this report.

## Static capacitance and distributed response

At low frequencies, a selected six-terminal imaginary-Z block is represented as `-A/r + B*r + D*r^3`, with `r=f/(1 GHz)`. Solving at 0.5, 1, and 2 GHz gives `C_fF = inv(A) * 1e6/(2*pi)`. A two-frequency variant checks extrapolation sensitivity; 1.5 GHz is held out. Ground capacitances are matrix row sums and mutual capacitances are negative off-diagonal entries. The accepted low-frequency held-out relative Z error is 8.24754e-6.

This matrix is a terminal representation; its entries are not relabeled as original ladder C0/Cg. For example, the accepted T3 ground row sum is about 339.20 fF, whereas the normal comparator T3 C0 is 251.6275 fF. They have different definitions. Distributing frequency-dependent loading between a static capacitance and an effective inductance is model dependent.

For the operating band, the full 7×7 Z matrix is interpolated in powers `r^-1, r, r^3, r^5, r^7`, now with `r=f/(5 GHz)`. Training frequencies are 4.65, 4.861, 5.0, 5.25, and 5.4 GHz. The independent 5.096 GHz point gives restored S21 error 1.45418e-6. Against the dense native two-port export, RMS is 8.11746e-6 for S21 and 3.44610e-5 for all S. The 1,501 export samples are not 1,501 independent EM solves.

These powers accommodate capacitive low-frequency behavior and smooth distributed corrections. They are an empirical local approximation, not a globally causal/passive broadband model or independent measurements of five physical inductances. The saved validation supports the stated band and geometry only. Final normal L values were optimized in tight bounds; they are not fresh explicit-capacitor calibrations of the modified meanders.

## Identifiability and actual geometry updates

`terminal_profiles.py` tests six ground terms and fifteen mutual terms. For seven primary terms, it fixes a change of ±0.5%, ±1%, or ±5%, then refits the other twenty positive factors within 0.2–5.0 using all 49 complex S entries. The impedance change is `B.T @ (inv(Cnew)-inv(Cref)) @ B / (j*omega)`. The mapping has rank six. The distributed background is held fixed. These are deterministic compensation profiles, not statistical confidence intervals or proof of unique original ladder values. Full profiles are packaged in `terminal_profiles.json`.

Crucially, geometry prediction uses changes in the **whole terminal Z**, not only static C. Native ±dimension controls provide finite-difference derivatives. Local proposals use `Znew ~= Zbase + sum(delta_d[k] * dZ/dd[k])`. Sparse difference models are checked at held-out resonance-region frequencies and against a separate native two-port simulation; inadequate models trigger remeasurement. The difference-model gate is 0.0003 for held-out S21, local all-S RMS, and dense S21 RMS, plus restoration. Early failed directions are not accepted merely because their surrogate predicts improvement.

`reachable_iterate.py` and its fine/control extensions store proposals before EM, compare predicted and actual changes, shrink trust regions after failures, and update derivatives after accepted measurements. The final search uses existing geometry dimensions on the 2 µm grid, enumerating nearby manufacturable rail choices. Later analytic Jacobians accelerate numerical optimization; they do not change the EM model or objective. Historical trajectories are in the copied history evidence rather than presented as universally reliable gradients.

`quarter_target_optimizer.py` minimizes the equal-frequency, unweighted real/imaginary residual to a normal ladder. There are 15 multiplicative component factors bounded to [0.9975,1.0025], and total reference delay bounded to ±100 ps. SciPy `least_squares` uses several delay starts. Reported RMS is `sqrt(sum(abs(S21native-S21model)^2)/N)`, equivalently `sqrt(2*cost/N)`. Delay multiplies S21 by `exp(-j*2*pi*f*tau)`; it does not alter magnitude. The accepted total delay is 3.58104935 ps. Low residual and tight bounds establish existence of a nearby comparator, not physical parameter uniqueness or guaranteed full-array transfer.

## Final geometry and boundary correction

Compared with the recorded refinement baseline, the accepted changes include T3 pad width +34 µm inward on each side; T3 meander width 130→122 µm and horizontal runs 68→72; T3–T2 IDC span +16 µm; T2–T1 facing rail width +16 µm; T2/T1 pad-height changes −20/−32 µm; and tuned T1–regular facing rails. Coil vertical-position and input lower-ground shifts are recorded in the geometry JSONs. The regular body and pads were retained in the final accepted variant. The upper ground edge remains straight. Do not infer the final geometry from one old selection file: use the complete packaged native project and metadata.

`n` counts horizontal runs, including partial end runs (`geometry/meander.py`); calling 68→72 an increase in *turn count* is imprecise. Both external native reference-plane shifts are 400 µm. The two-port reference impedance is 50 Ω, which is distinct from the computed CPW characteristic impedance.

The original 1.5 mm free-space top boundary was too close for boundary-independent results. A same-metal comparison showed S21 RMS change 0.00966824 on increasing it to 12 mm, with maximum power deficit dropping from 0.177283 to 0.00104456. These are observed boundary sensitivities, not proof that all missing power was physical radiation or metal loss. The final 12/30 mm comparison gives S21 RMS difference 3.14856e-9 and maximum all-S difference 9.65177e-8. It validates top-height convergence for these settings, not continuum mesh convergence or a fabricated enclosure.

The accepted stack uses 12 mm upper air, 500 µm silicon, 5 mm lower air, and the saved covers/materials. Grid 2×2 µm, MEDIUM_SPEED, edge meshing, double precision, de-embedding, 36 threads. Earlier coarse-mesh screens were not substituted for final evidence. A 36-thread same-metal pilot took 70 s per frequency versus 95 s at 16 threads, with S21 difference 1.08e-12; that is one-point evidence, not a universal runtime guarantee. Earlier four/eight-thread full-band timing was 1094.40/628.99 s with identical saved exports.

The isolated short-input check reduced clearance 740→440 µm and worsened comparator RMS to 0.00885226; original-versus-shortened magnitude RMS was 0.00774039. The user retained the original clearance. This rejected test is not the input used by the running full filter.

## Full-filter evidence is separate

Ten regular cells plus three mirrored tapers at each end make sixteen resonators. The conditional ladder prediction repeats the regular values, omits the output fixture, mirrors taper components, and omits fixture phase delay. At 4.92 GHz it predicts |S21|=0.957403974 (−0.378095480 dB). Direct EM gives 0.972139673 (−0.245426656 dB), a +0.132668825 dB difference, in 867.274 s. Phase was not aligned, so this is a magnitude sanity check at one frequency only.

The separate ABS sweep covers 4.5–5.5 GHz with 2,000 output intervals and 36 threads. Its input is copied under `full_filter_sweep_input/`; no partial response is published as final. Direct solve count and ABS SNR are progress observations, not percent complete or stopband attenuation. See the timestamped mailbox snapshot for live status. Completing the taper comparator criterion does not automatically validate passband ripple, loss, or band edges of the full device.

## Reproduction and remaining limits

Run `uv sync --locked`, `uv run python tools/check_published_files.py`, and `uv run python -m geometry.rebuild_final --verify-only`. The design registry includes the three current projects; GDS/native equality, grid and ports are checked without Sonnet. Native regeneration uses `--output-dir build/final_designs` and imports/reads projects without launching EM. Frozen layout replay is self-contained; it does not claim every historical campaign is portable.

`uv run python tools/verify_current_response.py` independently replays the saved network coefficients and bounded ladder comparison without EM. Original source paths and provenance are retained. The old root design folders are not silently overwritten. Solver caches, investigations, logs and environments stay excluded by directory-based Git policy.

Unresolved: unique physical C0/Cg/L identification; extrapolation beyond the validated band; physical material uncertainty; finer-than-2 µm convergence; effect of alternative valid external reference planes; and complete full-filter verification. Native `Ls=0.0388` is recorded, but the specific surface-model unit convention still needs authoritative confirmation; a general project inductance unit of NH alone does not settle the Ls field units. Do not repeat the old pH/square assertion as verified.
