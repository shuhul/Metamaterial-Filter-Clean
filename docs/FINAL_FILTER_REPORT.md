# Completed full-filter simulation

The completed design contains ten regular resonators between mirrored three-stage tapers, connected to two 50-ohm external ports. Its native Sonnet response has a broad, high-transmission band near the intended 4.8-5.2 GHz range. The final full-filter simulation contains no added internal delta-gap ports. It independently tests the complete metal geometry rather than relying on a cascade of the small-taper model.

The original input is `final_designs/final_full_filter/project.sonx`. Sonnet completed 123 adaptive direct frequency solves over 4.5-5.5 GHz on November 10, 2026 at 04:08:02 Pacific, reporting no errors or warnings and final ABS SNR 40.9246 dB against a 40 dB goal. The elapsed time to that native completion was about 41 hours 19 minutes 5 seconds. Resource-monitor elapsed time was 41 hours 20 minutes 13 seconds, including the subsequent process shutdown interval. These are different timestamps, not contradictory solve durations.

The workflow initially flagged the output as failed because its export validator required at least 2000 samples, while the native project was configured to export only the directly solved frequencies. The 123-point file is valid and agrees with the printed solver results to their precision. An isolated copy of the completed cache was then used to enable native interpolated output, with `-AbsNoDiscrete` prohibiting new direct solves. Sonnet produced the 2001-point curve with zero new frequencies and zero errors or warnings. The original working input and cache were preserved. The original wrapper error is retained in the evidence folder rather than erased.

## Main results and figure captions

Use `final_designs/final_full_filter/figures/annotated_main_passband_linear.png` as the main performance figure. Caption: **Completed full-filter Sonnet transmission magnitude.** The vertical blue lines mark the requested 4.8 and 5.2 GHz band edges. The green region marks the main continuous interval where |S21| is at least 0.8. The dotted orange lines are its estimated boundaries, 4.8260 and 5.1681 GHz, giving approximately 342.1 MHz bandwidth. The curve is Sonnet's native ABS interpolation, exported on a 0.5 MHz grid; threshold crossings are interpolated between those exported samples. The cutoff is amplitude, equivalent to -1.94 dB and 64% transmitted power for the equal-reference ports. The entire requested 4.8-5.2 GHz interval does not satisfy this threshold.

Use `figures/reference_vs_completed_sonnet.png` as the main qualitative reference comparison. Caption: **Reference-paper panel (left) and the completed full-filter simulation (right).** Both use 0 to -80 dB transmission, but their frequency ranges differ: the right panel is 4.7-5.3 GHz to give approximately comparable visual proportions. The reference image retains both its no-taper and with-taper curves; our panel contains only the tapered simulation. Both tapered curves show a flat central band, upper-edge structure, a notch and a finite stopband tail. This is a visual comparison, not a quantitative equality of specifications, bandwidth or rejection. No untapered version of our geometry was simulated for this figure.

Use `figures/sonnet_interpolated_s21.png` for a main or supporting comparison with the ideal circuit. Caption: **Completed Sonnet full-filter response versus the conditional ladder prediction.** Panels show linear magnitude, transmission in dB, and passband detail. Red is the native 2001-sample ABS curve; red markers show 123 direct solves; blue is the previously generated conditional circuit prediction using the accepted comparator parameters, repeated regular cells and mirrored tapers. The full-filter curve is not fitted to these data in this figure. Agreement is close through much of the central band; edge locations and depths differ. The circuit prediction is not a measurement and the dense samples are not 2001 independent EM solves.

The dB representation uses `20*log10(|S21|)`. On the wide 0 to -80 dB view, the simulated tapered response resembles the reference's tapered shape. That resemblance supports the design outcome without establishing equal stopband rejection or identical edge behavior. The reference itself also contains upper-edge dips, a recovery and a notch; those features should not be described as uniquely defective behavior in our filter.

## All requested plots and placement

| File under figures/ | Axes / definition | Recommended placement |
| --- | --- | --- |
| annotated_main_passband_linear.png / .pdf | GHz vs linear transmission magnitude; main contiguous >=0.8 interval only | Main; user's final annotation preference |
| reference_vs_completed_sonnet.png / .pdf | Side-by-side paper panel and ours; dB, different GHz ranges | Main, qualitative comparison |
| sonnet_interpolated_s21.png | Linear, dB and passband zoom; native ABS plus prediction and direct points | Main/supporting |
| sonnet_reference_style.png | Our native ABS curve alone, 4.7-5.3 GHz, 0 to -80 dB | Alternative main figure |
| completed_convergence.png | ABS SNR in dB versus SNR-check number | Appendix; this is interpolation convergence, not stopband rejection |
| full_filter_s21.png | Earlier completed direct-point plot with connecting guide lines and circuit prediction | Appendix/history; superseded for continuous curve display |
| annotated_achieved_passband.png / .pdf | Earlier dB annotation with all above-cutoff intervals | Appendix/history only; user replaced it with linear main-band-only annotation |

All original requested completed plots are preserved. The last two rows are historical presentations, not the recommended final report figures. Their companion validation JSON files describe the source samples. The narrow peak outside the main achieved band is not assigned a reliable bandwidth from one above-cutoff sample on the 0.5 MHz export grid.

## Numerical settings and remaining limits

The native full-filter calculation used a 2-by-2 micrometre grid, medium-speed edge meshing, double precision, de-embedding, 36 threads, and 400 micrometre reference-plane shifts at both 50-ohm external ports. It required 41,299 subsections and reported 14 GB peak memory. The stack has 12 mm upper air under a free-space top boundary, 500 micrometre silicon, 5 mm lower air and a perfect-conductor bottom. Materials and exact settings are embedded in the native file; no external component models are required to open it. Licensed Sonnet and a writable directory are needed for a new run. The .sonx alone does not enforce the workstation launcher's 36-thread policy.

ABS convergence supports the interpolation of this simulation. It does not prove finer-mesh convergence, exact fabricated material parameters, behavior in an actual package, or unique physical extraction of each ladder component. The full-filter performance is now supported by the native simulation itself; earlier taper-fit bounds were a design criterion, not a guarantee of full-array performance.

## Evidence and replay

`metrics.json` records thresholds, runtime definitions, counts and response/input hashes. `response_direct.snp` contains the original 123 samples; `response_interpolated.snp` contains the native 2001-sample export. `interpolated_output/project.sonx` differs from the original only in enabling interpolated Touchstone output and is the convenient project for a future dense-output run. It is not a newly solved geometry. Completion, cache-only export and wrapper-error evidence are in `evidence/`.

From the repository root, `uv run python tools/analyze_final_full_filter.py` validates the saved response against the logs, recalculates metrics and regenerates the annotated main-band and convergence plots. It does not invoke Sonnet. `uv run python tools/check_published_files.py` checks packaged hashes; `uv run python -m geometry.rebuild_final --verify-only` checks frozen native/GDS geometry and ports. Selected comparison figures are preserved rendered artifacts; the script does not claim to regenerate every historical presentation.
