# Recorded geometry-refinement trajectory

These are saved native iteration records, not newly run tests. A low predicted error is not an achieved result. Separate campaigns use different bounds/settings; compare only within the stated campaign. Reasons not recorded here must not be inferred from directory names.

| Campaign | Records | First actual RMS | Last actual RMS | Last predicted RMS | Last prediction error |
|---|---:|---:|---:|---:|---:|
| half_percent | 2 | 0.017308389589291494 | 0.017029941912746274 | 0.01678102271618207 | 0.0014903079257567048 |
| half_percent_all_LC | 2 | 0.00959965576012396 | 0.009421518967794924 | 0.009192247119719476 | 0.0012634315867117085 |
| half_percent_rail_spacing | 1 | 0.004633961997885624 | 0.004633961997885624 | not recorded | not recorded |
| half_percent_T1_pad | 3 | 0.007756049794268594 | 0.006534072119510343 | 0.0064595511138115175 | 0.0009752857738644839 |
| iteration_fine | 4 | 0.028263076493604145 | 0.007901471555329307 | 0.007861929874156633 | 0.0009469941631200972 |
| quarter_percent | 4 | 0.0061085971297451256 | 0.0037412945115759513 | 0.003767522631831491 | 0.0009413718435652373 |
| quarter_percent_body_position | 4 | 0.0037469454454327456 | 0.0037115067077685697 | 0.003675183128327214 | not recorded |
| quarter_percent_converged_cover | 8 | 0.002902406797511282 | 0.0024918633384848348 | not recorded | not recorded |
| quarter_percent_direct | 1 | 0.0038689902881285415 | 0.0038689902881285415 | 0.003613651226326912 | not recorded |
| quarter_percent_ground_width | 2 | 0.004387751410411329 | 0.004049451235034293 | 0.003661137226544254 | not recorded |
| quarter_percent_input_ground | 3 | 0.004122908959244891 | 0.003615890682690924 | 0.003603810516615021 | not recorded |
| quarter_percent_meanders | 1 | 0.00357635377438462 | 0.00357635377438462 | 0.003569575164908222 | not recorded |
| quarter_percent_regular_pads | 1 | 0.008011650596368151 | 0.008011650596368151 | 0.003164082036201824 | not recorded |

The final cover-converged sequence ends at 0.0024918633384848348, with native terminal-network and 30 mm boundary checks. The earlier regular-pad control overshot its surrogate prediction; its last recorded actual RMS is about 0.00801, and it was not the retained final geometry. The final retained regular body/pads are unchanged.

For the historical T3 selection, staged bounded fitting, shorted-response transfer failures, and subsequent change in objective, see CURRENT_METHOD.md and the unmodified reports in this directory. The history JSON includes exact per-step source projects and recorded reasons where available. This is an auditable sequence of recorded results, not a claim that every historical experimental intent has been reconstructed.
