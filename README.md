# Metamaterial Filter Automation

Filter geometry, Sonnet simulation automation, explicit-capacitor inductance calibration, and equivalent-circuit fitting.

## Selected deliverables

All publishable design artifacts live under **[final_designs/](final_designs/README.md)**. Each design contains its Sonnet project, approved GDS geometry, parameters, and saved results where available. Calibration inputs and supporting responses are included inside that design's folder.

```text
final_designs/
    regular_cell/
    taper_1/
    taper_2/
    taper_3/                 # historical cleanup assembly
    full_filter/             # historical full-filter input
    current_2026_11_08/      # accepted taper and separate full-filter evidence
    references/              # the pre-cleanup T3 baseline
reference_taper_geometry.png
```

**Current checkpoint (November 8):** [current_2026_11_08](final_designs/current_2026_11_08/README.md) contains the accepted taper, completed 4.92 GHz full-filter check, and separate ongoing-sweep input. Taper complex S21 RMS is 0.00249186334 to a comparator whose 15 L/C values are within 0.25% of nominal. These are not uniquely extracted physical values. The full-band filter sweep remains preliminary. See [the current method and evidence](docs/CURRENT_METHOD.md).

| Cell | L0 (nH) | C0 (fF) | Input Cg (fF) |
|---|---:|---:|---:|
| Regular | 4.9194 | 190.4 | 8.2865 |
| T1 | 5.97 | 148 | 15.3 |
| T2 | 4.33 | 149 | 153 |
| T3 | 2.77 | 251 | 261 |

These are nominal targets. The older staged cleanup fit (RMS 0.020312, 5% bounds) remains in the historical folders. The current checkpoint supersedes it; measured, fitted, and optimized comparator values must be distinguished.

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
