# Completed full filter

Open `project.sonx` for the exact final input used in the completed native sweep. `interpolated_output/project.sonx` preserves the same geometry/settings but enables dense native Touchstone output for a future run. Neither needs a source GDS or Python dependency merely to open; licensed Sonnet and a writable directory are needed to simulate.

- `response_direct.snp`: 123 directly solved samples.
- `response_interpolated.snp`: 2001 native Sonnet ABS samples across 4.5-5.5 GHz.
- `metrics.json`: verified completion, runtime, hashes and main-band threshold calculation.
- `figures/`: all user-requested completed plots, including superseded presentations retained for history.
- `evidence/`: native completion, zero-new-solve interpolation export, original wrapper error and resource/timing summaries.
- `verification_4p92GHz/`: earlier independent one-frequency check.
- `COMPLETED_SOURCE_PATHS.json`: original response/figure/evidence sources.

Recommended report figures: annotated_main_passband_linear.png, reference_vs_completed_sonnet.png and sonnet_interpolated_s21.png. See ../../docs/FINAL_FILTER_REPORT.md for full captions, placement and limitations. Speed claims are audited in ../../docs/SIMULATION_SPEED_EVIDENCE.md.

Replay from repository root: `uv run python tools/analyze_final_full_filter.py`. No EM is launched. The original wrapper status says failed because it expected a dense export count from a direct-only export; native EM completed successfully. Do not use that stale wrapper flag as the final simulation verdict.
