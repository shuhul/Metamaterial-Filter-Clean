# Selected filter designs

All paths below are relative to this directory. Open `project.sonx` to inspect a design in Sonnet. `layout.gds` is its approved frozen metal; `parameters.json` preserves its original geometry settings.

| Folder | Role | Saved response |
|---|---|---|
| [regular_cell](regular_cell/) | Selected regular-cell array | Yes |
| [taper_1](taper_1/) | T1 baseline assembly | Yes |
| [taper_2](taper_2/) | T2 baseline assembly | Yes |
| [taper_3](taper_3/) | Current T3/T2/T1/regular cleanup assembly | Yes |
| [full_filter](full_filter/) | Ten regular cells and mirrored three-taper sections | **Still running; input only** |
| [references/taper_3_before_cleanup](references/taper_3_before_cleanup/) | Pre-cleanup baseline used to establish the current assembly | Yes |

Each completed assembly has `response.snp`, its recorded fit or selection, and an S21 plot. `calibration/` includes available native calibration inputs and supporting results. The current T3 cleanup shortens the T3-T2 IDC rails to the finger ends and centers T2's coil and pads. Its staged bounded fit RMS is 0.020312; it is not proof of unique extracted capacitances or a completed full-filter target match.

The full-filter native input uses a 2 um mesh and 50 ohm reference impedances. The right taper section is mirrored. Its EM verification is unfinished, so no incomplete full-filter response is presented as final.

See [the rebuild guide](../docs/REBUILD.md). Frozen native/GDS replay is verified for all six projects, and the full-filter geometry regenerates from the current taper assembly and repeated regular body. Original provenance retains old paths and hashes; the replay tool uses only these repository-contained files.
