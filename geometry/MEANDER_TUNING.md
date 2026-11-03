# Native filter_v6 meander tuning

Run from the project root:

```powershell
.\.venv\Scripts\python.exe -u geometry\tune_meander.py --max-trials 3 --refine-best
```

Copies of the original `sonnet_files/filter_v6.sonx` vary `in_width`, the
horizontal span of each meander run (`meander.w`). The native Sonnet parameter
moves the three linked meander polygons, including the end structures.
Vertical geometry, runs, trace/gaps, caps, pads, box, materials and 2 um mesh
are preserved. Trial reference planes are shortened from 600 to 400 um because
the original setting produces warnings about extending beyond uniform feeds.
The original project and shared design targets are preserved.

Screening uses adaptive 4.35-4.95 GHz with Target=500, matching the original
interval target over a narrower range. Sonnet selects a convenient actual
export spacing; inspect each result's export_frequency_spacing_MHz.
The final inductance check uses the same range with Target=1500. Control
resampling experiments are saved in spacing_check.json; these estimate
sampling error, not mesh convergence or total physical uncertainty.

Each trial archives 20/30/40/50 fF capacitor sweep snapshots, Touchstone files,
solver logs, fitted L0 and regression plots. The initial 190 um control was
already completed at higher density and also has a final 0.1 fF capacitance fit.
Subsequent trials focus only on inductance to save time. The selected project
is restored to C_add=0.0001 pF and the original 4-6 GHz/500 interval sweep
after refinement. The last raw export is renamed last_sweep_cadd_50ff.snp; no final small-capacitor S21
simulation is claimed for the selected geometry. Run a new small-capacitor
simulation before using that project's raw export for capacitance fitting.

Projects: sonnet_files/meander_tuning_ref400/w*/filter_v6_w*.sonx.
Results: output_data/meander_tuning_ref400/w*/.
Final refined extraction: output_data/meander_tuning_ref400/best_refined/.
Overview: summary.json and length_vs_inductance.png in the results root.

Completed trial reuse requires matching source/result hashes. Interrupted
trials may resume from validated archived sweep points. The interrupted
600 um calibration experiment remains separate under meander_tuning/ and is
not a completed inductance result. These copies use the original filter_v6,
not the generated layout. No N=4 geometry is changed by this experiment.
