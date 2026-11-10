# Rebuilding the selected designs

The canonical registry is `final_designs/designs.json`. It contains final_regular, final_T1, final_T2, final_T3, final_full_filter, and its verification_4p92GHz subfolder. The standalone regular/T1/T2 selections are staged baselines; final_T3 contains the later embedded geometry. Superseded publication copies are in final_old and are not part of the default replay.

Each registered input has project.sonx (complete native settings and metal), layout.gds (the same frozen metal), design.json, and parameters.json. Finished response evidence is included where available. The full-filter folder contains the completed direct and interpolated full-band responses, evidence, metrics and the earlier one-frequency verification subfolder. `uv run python tools/analyze_final_full_filter.py` validates completed evidence and regenerates the principal annotated/convergence plots without EM. Original absolute paths in provenance JSON remain historical; they are not required by the replay.

## Check without Sonnet

```powershell
uv sync --locked
uv run python tools/check_published_files.py
uv run python -m geometry.rebuild_final --verify-only
uv run python tools/verify_current_response.py
```

Checks cover current and archived artifact hashes, current native/GDS metal equality, grid and external ports, and saved taper network/comparator replay. Native regeneration of current designs uses their exact frozen layouts rather than relying on old campaign directories. The legacy block assembler remains in the source for the archived cleanup filter only; it is not used to rebuild the new full-filter geometry.

## Import and read native projects

```powershell
uv run python -m geometry.rebuild_final --output-dir build/final_designs
```

With another Sonnet installation, add `--sonnet-bin 'C:\path\to\Sonnet\bin'`. The tool imports the frozen GDS, reattaches native port edges, preserves materials/stack/mesh/reference planes/controls/sweeps, compares the metal, and asks Sonnet to read the result. It does not start em.exe. It needs a compatible Sonnet installation and license. Current inputs were produced with Sonnet 19.52.2025.261.

## Simulate deliberately

The existing project.sonx is sufficient to open the frozen design in Sonnet; no source GDS or old cache is required. Use a writable independent run directory and a licensed installation for any new EM analysis. Sonnet regenerates sondata; solver caches remain excluded from Git. Published native files do not automatically enforce the workstation launcher's 36-thread policy.

`final_designs/PUBLICATION_PATHS.json` maps the former publication paths to current names. `SOURCE_PATHS.json` maps copies to original working sources. Running solver directories were not moved. See [the reader guide](SECTIONS_10_11_READER_GUIDE.md) for the two different fit errors, the single-point/full-band distinction and limitations.
