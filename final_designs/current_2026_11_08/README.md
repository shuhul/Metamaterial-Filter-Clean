# Current accepted checkpoint

| Directory | Status | Open |
|---|---|---|
| taper | Completed taper fixture, all 15 comparator L/C within 0.25%, complex RMS 0.00249186334 | [project](taper/project.sonx), [continuous ground](taper/continuous_ground/project.sonx), [fit and values](taper/s21_with_target_parameters.png) |
| full_filter_single_point | Completed 4.92 GHz direct EM only | [project](full_filter_single_point/project.sonx), [comparison](full_filter_single_point/comparison.png), [result](full_filter_single_point/complete.json) |
| full_filter_sweep_input | Snapshot of input to ongoing ABS run; no completed response claimed | [project](full_filter_sweep_input/project.sonx) |

The terminal model is independently validated, but normal comparator values are not unique physical extractions. Full-filter passband performance is not established by the taper fit or the one-frequency check. Original 740 um input clearance retained; shortened-input variant rejected.

All working files remain in their original locations. [SOURCE_PATHS.json](SOURCE_PATHS.json) records old-to-new copy mappings and hashes. Relative files in this directory are the published inputs; absolute paths inside historical JSON remain provenance, not portable dependencies. The response in taper/project.snp corresponds to taper/project.sonx. The merged continuous-ground copy preserves metal coverage; it is not a separate taper EM verification.

See [current method](../../docs/CURRENT_METHOD.md) for equations, settings, history, limits, and verification commands. The earlier final_designs/taper_3 and full_filter folders are preserved historical checkpoints.

Figure captions: taper/geometry.png shows the accepted fixture; s21_with_target_parameters.png compares native EM to the bounded normal ladder; s21_verified.png also shows the independent terminal model. conditional_ten_cell_s21.png is a circuit-only prediction in linear magnitude and dB, not a native full-filter curve. full_filter_single_point/comparison.png compares the one completed EM point to that prediction. The root reference_passband_plot.png is the reference-paper crop, not our simulation.
