# Shared-parameter model comparisons

| Dataset choice / model | Short complex RMS | Full complex RMS | Full magnitude RMS | Held-out full complex RMS |
|---|---:|---:|---:|---:|
| results/full_only_fixedL | 0.092998 | 0.202914 | 0.182434 | 0.366418 |
| results/joint_fixedL | 0.076069 | 0.200689 | 0.180509 | 0.360321 |
| results/joint_fixedL_allS | 0.069996 | 0.283701 | 0.081235 | 0.487580 |
| results/short_only_fixedL | 0.002864 | 1.035175 | 0.273951 | 1.131173 |
| results_all_direct/boundary_fixedL | 0.054422 | 0.085609 | 0.048758 | — |
| results_all_direct/boundary_fixedL_allS | 0.061573 | 0.145244 | 0.058041 | — |
| results_all_direct/boundary_fixedL_magnitude | 0.049747 | 1.388606 | 0.026067 | — |
| results_all_direct/boundary_series_leads_fixedL | 0.047700 | 0.081498 | 0.045878 | — |
| results_all_direct/joint_fixedL | 0.075823 | 0.101286 | 0.057767 | — |
| results_all_direct/joint_fixedL_magnitude | 0.231047 | 1.312396 | 0.028989 | — |
| results_all_direct/joint_L10pct | 0.077466 | 0.089875 | 0.051212 | — |
| results_all_direct/joint_L1pct | 0.074644 | 0.099784 | 0.058037 | — |
| results_all_direct/joint_L3pct | 0.075170 | 0.097314 | 0.056313 | — |
| results_all_direct/joint_series_leads_fixedL | 0.074149 | 0.095858 | 0.055441 | — |
| results_all_direct/joint_series_leads_fixedL_allS | 0.070986 | 0.178399 | 0.056627 | — |
| results_all_direct/joint_shared_loss_fixedL | 0.074397 | 0.101064 | 0.058944 | — |
| results_three_config/three_config_allS | 0.057524 | 0.156554 | 0.067504 | — |
| results_three_config/three_config_complex | 0.057379 | 0.092245 | 0.054377 | — |
| results_three_config/three_config_leads_complex | 0.047470 | 0.084120 | 0.048863 | — |
| results_three_config/three_config_magnitude | 0.061391 | 1.252227 | 0.045703 | — |

Short-only/full-only cross-configuration values do not refit electrical parameters; separate full phase nuisance parameters are untrained in a short-only fit. Magnitude error is independent of that phase ambiguity. `results` trains the first 72 full points and holds out 27. `results_all_direct` fits all 99 captured direct full points and has no independent frequency holdout. `results_three_config` adds the regular-only reference with separately modeled input fixture.

Low magnitude RMS alone does not imply correct complex/reflection response. No table entry is accepted solely for minimizing an objective.
