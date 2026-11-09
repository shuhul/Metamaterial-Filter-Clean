# Supplied AWR target-design inputs

`Filter4_Final_Version.emp` and `AWR_optimizer.png` were supplied with mailbox request 0010. Their original bytes are preserved. The native project is an AWR file; it was not executed or optimized during this publication task. Opening/reproducing its optimization requires a compatible licensed AWR environment, and any project dependencies should be checked there. It is not an input dependency of the standalone Sonnet project.

The inspected screenshot shows the enabled `Full_Filter:|S(1,1)|` goal: measurement less than 0.001, 4.82–5.18 GHz, Difference, default L=2, weight 1. These are optimizer settings, not proof the circuit or geometry achieves them. The user reports S11 minimization worked better than directly maximizing S21. No controlled comparison of those optimizer choices was performed here.

The requested filter band remains 4.8–5.2 GHz; the pictured goal covers its interior. For an ideal lossless matched-reference two-port, reflection and transmission powers sum to one. With dissipation or radiation, low S11 alone does not guarantee high S21. See `docs/SECTIONS_10_11_READER_GUIDE.md` for the distinction between this supplied circuit objective and the existing taper geometry objective.
