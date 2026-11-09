# T2-to-T3 history addendum

Publication paths in this dated record were subsequently simplified; see `final_designs/PUBLICATION_PATHS.json` for current locations. Original campaign paths are unchanged.

Paths below are relative to the project root unless noted. This supplement answers mailbox request 0007 using saved records, not new simulations. `T3DATA` abbreviates `output_data/taper_exploration/three_tapers_v3_T3_2026-10-30/`.

## Initial T3 extension and failed constrained fits

The starting point was `geometry/taper/selected_taper2_cap5.json`. T3 replaced the former T2 input electrode with a resonator while retaining downstream metal, verified by translated metal XOR. The initial design had a 480×590 µm asymmetric U-capacitor, an 84 µm wide meander with 72 horizontal runs and 4 µm wire/gaps, and 70-finger input/junction IDCs. The intended sequence was explicit-capacitor L calibration, direct feed without the new input IDC, then the complete fixture. Targets were 2.77 nH, 251 fF, 261 fF input coupling, and 153 fF T3–T2 coupling. At that time acceptance meant new capacitances within 5%, calibrated L within 2%, and complex fit RMS ≤0.025.

The first measured L was 2.719061755 nH (−1.8389%), using explicit 500/700/900/1100 fF loads. Passing that L check did not establish a working taper. Source: `T3DATA/README.txt`; implementation: `geometry/taper/taper3.py` and `run_taper3.py`. Clearance, pad height/width, coil run count and IDC dimensions were subsequently explored; ledgers retain individual cases rather than a single geometry parameter varying independently.

The October 31 `CURRENT_T3_STATUS.txt` records a particularly informative failure, `regular_w80_n72_dl476_dr500_padL460_R24_top2_i84ni33_stub200_FIXED`: calibrated L=2.762630600 nH, fitted C0=254.5111 fF, input Cg=253.2206 fF, fixed junction=156.6004 fF, but RMS=0.0548982. Component estimates were near target while the fit failed. It was explicitly not accepted.

That status also records why simple explanations were insufficient: the older T2 S21 RMS was 0.00624769 but S11/S22 errors were 0.14371872/0.11671336; redistributing delay improved S22 but left S11=0.1123423. An older-EM-block diagnostic still gave RMS 0.04319382. Magnitude-only fitting reached 0.03108, so phase alone was not a cure. A compact-coil operating-band L check changed 2.74255310 to 2.65660795 nH, but the corresponding fit remained about 0.052. These were diagnostics, not adopted calibrations or proof of a single cause.

## Retaining the old junction

`geometry/taper/taper3_fixed_junction.py` restored the exact accepted T2 IDC: 70 fingers, 22 µm rail separation, 4 µm fingers/gaps, 2 µm tips and 8 µm rails. Its metal ROI XOR was checked. The previous extracted junction value, 156.6003999 fF, was fixed; only new T3 C0, input Cg and delay varied.

Recorded examples in `T3DATA/fixed_junction_ledger.json`: `w80n73_dr512_pad104x594_r404_i40_j22` achieved RMS 0.0217268 but C0=107.2927 fF, failing the target. The `...i84ni28_originalT2junction_FIXED` case gave RMS 0.0212773, C0=89.5002 fF and input Cg=213.2866 fF. Input-finger reduction 40→28 changed input Cg from about 357.832 to 213.287 fF in the expected direction, but that did not solve C0 and fit quality together. Retaining an IDC's metal is not proof that its earlier equivalent-circuit value transfers unchanged.

## Staged and historical bounded refits

`taper3_staged_bounds.py` tried 1%, 2%, 3% old-capacitance bounds anchored to the immutable selected T2 values, not accumulated from previous refits. It freed T2 first, then T2+T1, then regular parameters; the output-end terms remained fixed in this particular staged test. All L remained fixed and T3 C0/input Cg/delay were free. On the failed regular-w80 case, the broadest scope gave RMS 0.044314, 0.038975, 0.036381 for 1%, 2%, 3% respectively: none passed 0.025. At 3%, C0/input Cg were 310.202/183.732 fF. Source: that case's `response/staged_bounds/SUMMARY.txt` and `summary.json`.

The subsequent symmetric-pad historical refit expanded scope to the output end as well. `T3DATA/historical_staged_bounds/results.json` records 25 geometries, 300 fits, and 811.21 s of saved-data processing; no new EM. No case passed both fit and component criteria in that comparison. The lowest RMS was 0.004817 for `w130_t6g2n69_dl342_dr340_pad220x450_i104_j28n18`, but C0/input Cg were 361.724/105.044 fF. The closest target case under 1% old bounds and RMS<0.025 was `w130_t6g2n71_dl342_dr340_pad220x594_i56_j28n40`: RMS 0.0228033, C0=287.2401 (+14.44%), input Cg=250.8271 (−3.90%). This motivated geometry refinement rather than selecting the lowest residual. Sources: `taper3_historical_bounds.py`, `taper3_historical_report.py`, and the saved `REPORT.txt`.

## Wider band and known-capacitor probes

`taper3_wide_band.py` tested window dependence on `w130_t6g2n66_dl342_dr340_pad80x590_mc24_i134ni40_j28n16`. The saved 4.0–5.4 GHz fit had RMS 0.0189655 and C0=321.4291 fF; the narrow-window refit had RMS 0.0257070 and C0=321.4358 fF. A lower broad-band average did not make C0 meet target. Files: that case's `response_wide_4to5p4/fit/fit.json` and `fit_narrow/fit.json`.

An older `w84n73_dr256_pad104x594_i36_j22` geometry has completed explicit native capacitor settings 0.1, 100, 250 and 400 fF under `known_Cadd/`. `taper3_loaded_probe.py` produced fresh 100 fF direct checks at 5.1427 and 5.0069 GHz; complex S21 differences from the cached result were 5.1180e-6 and 5.2629e-6. The later joint investigation fit a common-denominator S(C) operator and predicted withheld 400 fF at RMS 2.39e-13, an algebra/cache consistency check, not independent broadband EM accuracy. The almost-open component differed from the unprobed geometry by RMS 0.00430, warning that port/component insertion needs validation.

The derivative diagnostic tried fitting d(1/S21)/dCadd so T3 L/C0 cancel algebraically. It gave input/junction couplings 89.649/125.457 fF, normalized derivative RMS 0.0184594, and a strongly frequency-dependent, complex inferred C0 (real percentiles −18.09/292.32/602.08 fF). It was not an accepted static extraction. Source: `known_Cadd/derivative_diagnostic/result.json`; independently summarized in `investigations/joint_physical_extraction_2026-11-04/REPORT.md`.

## First accepted T3 and cleanup

After the historical comparison, narrower tall pads and adjusted run counts were tested. `pad_refinement_saved_results.json` records 212 µm pads with 71 runs (3% bounds: C0 279.682, input Cg 235.273 fF) and 192 µm pads with 69 runs (3%: 270.443/248.185 fF); those entries were rejected. The eventually selected `w130_t6g2n68_dl342_dr340_pad184x594_i54_j28n40_pad_tune` is explicitly accepted in `geometry/taper/selected_taper3.json`: RMS 0.0200526, C0 256.1465 fF, input Cg 263.0294 fF, calibrated L 2.7923246 nH, fixed junction 156.6004 fF, with 3% older bounds. This supersedes the earlier status file saying no T3 had been accepted.

Next the user requested removing the IDC rail tails and centering T2's coil/pads. `t2_cleanup.py` cuts the 40-finger junction rails from y=−598 to −340 µm, retaining fingers. It translates the existing T2 coil and the two pads by −56, −52 and −54 µm respectively, with localized checks. `run_t2_cleanup_and_full.py` evaluates separate and combined variants, including explicit L checks where needed.

`output_data/taper_exploration/T2_cleanup_2026-11-01/cleanup_summary.json` records RMS: baseline 0.0204577; trimmed IDC 0.0183814; centered T2 0.0158044; combined 0.0203123. All passed the then-current 5%-bounded capacitance/0.025 fit rule. The combined version was selected to satisfy both requested geometry cleanups, deliberately not because it had the lowest RMS. The exact recorded selection reason and all source fit paths are in that JSON. These are still conditional equivalent-circuit fits, not unique measurements.

## Prepared interventions versus completed tests

The later joint investigation proposed new known capacitors at T2/T3 and pad-width interventions at 176/184/192 µm on the then-current accepted geometry. The old report described them as queued behind the full sweep. Their definitive later status is `investigations/joint_physical_extraction_2026-11-04/intervention_results/status.json`: `stopped_by_user`, `EM_started=false`, UTC 2026-11-05T21:23:22.8793533Z. Only that status file is present in the results directory. Thus those particular prepared projects did not complete their proposed EM study; do not confuse them with the older completed known-C data above.

The later terminal-network and measured geometry-derivative campaign superseded this as the active design workflow. That is an inference from the subsequent executed scripts/results, not evidence that the cancelled interventions were retrospectively validated. The current terminal method has its own restoration, held-out and native-geometry checks described in CURRENT_METHOD.md.

## Older full sweeps and the current sweep

The pre-cleanup `output_data/taper_exploration/full_filter/full_filter_10_regular_3_tapers_each_side/` run began October 31 10:54 UTC and ended November 1 00:54 UTC with failure return code 4294967295, after 50,428.62 s. Its log ends during 5.1744 GHz. The inspected status/resource files do not establish why that process ended or whether someone intentionally stopped it; do not invent a cause.

The cleanup sweep began November 1 11:20:14 UTC and was explicitly stopped by user request November 5 21:23:22 UTC after **106 h 3 min 8 s**, **114 completed frequencies**, with 4.994 GHz interrupted. It had not converged. `STOP_TIMING.md`, `stop_timing_record.json`, and `run_status.json` in its `...T2_cleanup_2026-11-01/` output folder document this. Both waiting workers were also stopped to prevent automatic follow-ups.

The current November 8 36-thread sweep is a different project in `investigations/filter_refinement_2026-11-05/quarter_percent_full_sweep_2026-11-08/`. It remains untouched by this publication work. Current counts belong in a timestamped status snapshot, not merged into the older timing record.

## Repository and temporary-file audit

The published tree contains `geometry/` (generators, import, circuit and stage helpers), `tools/` (verification/publication utilities), `docs/` (current method plus dated history), `final_designs/` (curated projects/results), `tests/`, `pyproject.toml`, `uv.lock`, and the two root reference images. Local campaigns live in ignored `sonnet_files/`, `output_data/`, `investigations/`; caches/environments/builds are excluded.

Main current entry points: `tools/check_published_files.py`, `geometry.rebuild_final --verify-only`, and `tools/verify_current_response.py`. Native replay with `geometry.rebuild_final --output-dir build/final_designs` imports/reads inputs without EM. Historical exploratory scripts may reference excluded campaigns; they are not advertised as portable workflow drivers.

At the audit inspection, `git ls-files '*tmp*' '*temp*'` returned no tracked temporary paths, and a recursive mailbox search found no remaining `*.tmp` files. In particular no `.terminal_profiles.json.tmp` is committed. A same-directory temporary name can appear briefly during an atomic write; readers should ignore dotfiles and `.tmp` until the final filename appears. No unrelated temporary files were deleted.
