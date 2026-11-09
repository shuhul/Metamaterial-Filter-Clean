# Metamaterial Filter Automation

Python geometry generation, Sonnet electromagnetic simulation, and verified circuit-response comparison for a superconducting metamaterial filter.

## Selected designs

The selected files are in [final_designs](final_designs/README.md):

- [final_regular](final_designs/final_regular/): standalone regular-cell baseline.

- [final_T1](final_designs/final_T1/): standalone T1 baseline.

- [final_T2](final_designs/final_T2/): standalone T2 baseline.

- [final_T3](final_designs/final_T3/): accepted assembled T3/T2/T1/regular taper fixture.

- [final_full_filter](final_designs/final_full_filter/): current full-filter input, with completed 4.92 GHz verification in its subfolder.

The accepted taper has complex S21 RMS 0.00249186334 to a normal circuit whose ten capacitances and five inductances are within 0.25% of nominal. Those are comparator values, not unique physical extractions. **Final full-band filter results are pending.** A one-frequency check is not full-band validation.

Superseded publication copies are preserved in [final_old](final_old/README.md). [PUBLICATION_PATHS.json](final_designs/PUBLICATION_PATHS.json) records the renames; [SOURCE_PATHS.json](final_designs/SOURCE_PATHS.json) links current copies to original working artifacts. Running simulations and their source directories were not moved.

## Read and reproduce

Start with the [Sections 10?11 reader guide](docs/SECTIONS_10_11_READER_GUIDE.md), [current technical method](docs/CURRENT_METHOD.md), and [rebuild guide](docs/REBUILD.md). The supplied AWR target project and optimizer screenshot are in [reference_inputs/awr](reference_inputs/awr/README.md).

```powershell

uv sync --locked

uv run python tools/check_published_files.py

uv run python -m geometry.rebuild_final --verify-only

uv run python tools/verify_current_response.py

```

To regenerate the frozen native projects from the packaged GDS layouts:

```powershell

uv run python -m geometry.rebuild_final --output-dir build/final_designs

```

The replay imports and reads projects; it does not start EM. Native project inputs include geometry, materials, mesh, ports, and sweep settings. Sonnet and a working license are required for native import or simulation. Caches can be regenerated and are excluded from Git.

Source is under geometry/ and tools/, documentation under docs/, and nominal electrical targets in geometry/specifications.py. Local trial campaigns and solver outputs stay in ignored investigations/, sonnet_files/, and output_data/. Historical scripts may depend on those campaigns; the frozen replay and saved-response verifier above are self-contained.

Reference images: [taper geometry](reference_taper_geometry.png), [reference passband](reference_passband_plot.png).
