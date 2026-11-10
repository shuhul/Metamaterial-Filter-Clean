# Publication policy

Publish source, documentation, reference artwork, supplied `reference_inputs/`, the curated `final_designs/`, and preserved publication history in `final_old/`. Use normal folder-based ignore rules: new files placed inside a final design are eligible for Git automatically unless they are ordinary logs, caches, or build products.

Do not add local `sonnet_files/`, `output_data/`, `investigations/`, virtual environments, duplicate handoff bundles, generated `build/` outputs, or Sonnet `sondata` caches. Those directories remain untouched on disk. No original campaign files are moved or deleted to prepare this repository.

All selected native `.sonx` projects, GDS layouts, saved `.snp` responses, calibration inputs, settings, previews, and the reference taper image are published. Completed full-filter responses and concise completion evidence are published. Large diagnostic fit ensembles, solver caches and raw resource time-series remain local. The selected-design index identifies the current assembly and its baselines.

`final_designs/manifest.json` records integrity hashes for current designs, archived publication copies, and supplied reference inputs. It is **not a Git allowlist** and does not control which files Git can track. Refresh it after intentionally updating a final design:

```powershell

uv run python tools/prepare_git_publication.py

uv run python tools/check_published_files.py

```

Inspect `git diff --stat` and `git status` before committing. `.gitattributes` preserves native/evidence file bytes so their hashes are stable across clones. Git LFS is unnecessary for the current selected package; its largest response files are only a few megabytes.

Historical documentation describes some older methods and paths. Follow the root README, `docs/REBUILD.md`, and `final_designs/README.md` for the current replay procedure. The old root README is retained under `docs/history/`.
