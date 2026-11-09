# Selected-design index

**Current checkpoint:** see [current_2026_11_08](../final_designs/current_2026_11_08/README.md). The older entries below are retained as historical baselines; their earlier acceptance criteria are superseded.
The publishable files are organized under [final_designs/](../final_designs/README.md).

| Design | Native input | Plot |
|---|---|---|
| Regular array | [project.sonx](../final_designs/regular_cell/project.sonx) | [S21](../final_designs/regular_cell/s21_fit.png) |
| T1 baseline | [project.sonx](../final_designs/taper_1/project.sonx) | [S21](../final_designs/taper_1/s21_fit.png) |
| T2 baseline | [project.sonx](../final_designs/taper_2/project.sonx) | [S21](../final_designs/taper_2/s21_fit.png) |
| Historical cleanup assembly | [project.sonx](../final_designs/taper_3/project.sonx) | [S21](../final_designs/taper_3/s21_fit.png) |
| Full-filter input, EM pending | [project.sonx](../final_designs/full_filter/project.sonx) | [Geometry](../final_designs/full_filter/geometry.png) |

The T3 pre-cleanup baseline is retained in `final_designs/references/`. The original reference drawing is [reference_taper_geometry.png](../reference_taper_geometry.png). Each design's parameters, metadata, calibration inputs, and saved response (where complete) are next to its project. Follow [REBUILD.md](REBUILD.md) for verification and native project regeneration from a fresh clone.

Nominal specifications are centralized in `geometry/specifications.py`. Current staged results and future full-filter verification are different evidence: the unfinished full filter must not be labeled as having achieved its target response.
