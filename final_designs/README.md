# Selected final designs

| Folder | Role | Verification |

|---|---|---|

| final_regular | Selected standalone regular array | Historical calibration/response |

| final_T1 | Selected standalone T1 fixture | Historical calibration/response |

| final_T2 | Selected standalone T2 fixture | Historical calibration/response |

| final_T3 | Current assembled taper fixture | Completed native response and independent terminal model |

| final_full_filter | Current 10-cell filter with mirrored tapers | Input; final full-band results pending |

| final_full_filter/verification_4p92GHz | Same full-filter geometry, single-frequency test | Completed 4.92 GHz EM |

Open project.sonx in a selected folder. layout.gds is its frozen metal, and design.json/parameters.json record its role. Older baseline responses are named response.snp; the current taper and single-point exports are project.snp. A standalone T1/T2 baseline is not identical to its later embedded version in final_T3.

For a continuous-ground taper-only drawing, use final_T3/continuous_ground/project.sonx. Its metal union matches the accepted original, but it has not been separately simulated as a taper validation. The full filter uses merged ground and has its own verification.

Final full-filter result placeholder: **await completed native sweep and validation**. No live point count or provisional full response is represented as final here.

PUBLICATION_PATHS.json records the publication-only rename mapping. SOURCE_PATHS.json records original source paths and hashes. Original solver/campaign paths are unchanged. Earlier alternative packages are archived in ../final_old/. See ../docs/SECTIONS_10_11_READER_GUIDE.md and ../docs/REBUILD.md.
