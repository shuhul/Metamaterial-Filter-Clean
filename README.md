# Metamaterial Filter Automation

Filter geometry, Sonnet simulation automation, explicit-capacitor inductance calibration, and equivalent-circuit fitting.

## Selected deliverables

All publishable design artifacts live under **[final_designs/](final_designs/README.md)**. Each design contains its Sonnet project, approved GDS geometry, parameters, and saved results where available. Calibration inputs and supporting responses are included inside that design's folder.

```text
final_designs/
    regular_cell/
    taper_1/
    taper_2/
    taper_3/                 # current cleanup assembly
    full_filter/             # ten regular cells, tapers at both ends
    references/              # the pre-cleanup T3 baseline
reference_taper_geometry.png
```

**Full-filter EM verification is still running as of November 3, 2026.** The full-filter folder contains the input project and layout, not a completed response. Current refinement goals are RMS below 0.005 and component errors below 1%; these are not achieved results.

| Cell | L0 (nH) | C0 (fF) | Input Cg (fF) |
|---|---:|---:|---:|
| Regular | 4.9194 | 190.4 | 8.2865 |
| T1 | 5.97 | 148 | 15.3 |
| T2 | 4.33 | 149 | 153 |
| T3 | 2.77 | 251 | 261 |

These are nominal targets. Measured inductances and fitted capacitances are recorded separately. The current staged cleanup fit has complex S21 RMS 0.020312 with capacitances bounded within 5% of nominal targets. Several values reach bounds; the fit does not establish unique physical capacitances.

## Rebuild from a fresh clone

Install Python and [uv](https://docs.astral.sh/uv/), then:

```powershell
uv sync --locked
uv run python tools/check_published_files.py
uv run python -m geometry.rebuild_final --verify-only
```

To regenerate the native projects from the frozen approved GDS inputs, with a licensed Sonnet installation:

```powershell
uv run python -m geometry.rebuild_final --output-dir build/final_designs
```

This imports geometry and reads the resulting projects; it **does not launch an EM calculation**. An alternate Sonnet installation can be supplied with `--sonnet-bin`. See [the rebuild guide](docs/REBUILD.md) for what is regenerated and how to run new simulations deliberately.

## Documentation and source

- [Selected designs](docs/SELECTED_DESIGNS.md) and [publication policy](docs/PUBLISHING.md).
- [Geometry guide](geometry/README.md) and [editable regular-cell runner](geometry/RUN_FILTER.md).
- [Electrical specifications](geometry/specifications.py) and [inductance measurement](docs/MEASURING_INDUCTANCE.md).
- [Sonnet process control](tools/README.md).
- [Reference taper geometry](reference_taper_geometry.png).

The reusable source includes older helper modules used during development. Historical campaign drivers can depend on excluded trial data and original Windows paths. The documented `geometry.rebuild_final` path is the self-contained replay route for approved designs; it does not depend on those campaigns. Original provenance JSON is retained unchanged, including historical paths and hashes.

Normal folder-based `.gitignore` rules keep campaigns, generated builds, environments, logs, and solver caches local. There is no explicit allowlist of tracked files. Sonnet regenerates `sondata` when a project is simulated; those caches are not required to reconstruct its metal, materials, ports, mesh, or sweep settings.
