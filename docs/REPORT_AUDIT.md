# Audit of the November 9 reference-only Notion drafts

Audited against source code and saved results; this is not an edit of either Notion page. The supplied references were `revision-main-reference.md` and `revision-appendix-reference.md` in the lab-share `results/local-reference/` folder. Preserve their older results as explicitly dated history.

| Draft section | Required correction or qualification |
|---|---|
| Main 1 | Target table remains applicable. Distinguish nominal targets, independently modeled terminal values, and optimized normal-comparator values. |
| Main 2 | The listed calibration-then-three-C fit describes earlier stages. Replace the current-method summary with the validated terminal network and measured geometry-derivative workflow in CURRENT_METHOD.md. |
| Main 3 | Confirmed: meander `n` counts horizontal runs, including partial first/last runs. It does not count complete turns. |
| Main 4 | Current native external planes are both 400 µm; both external impedances are 50 Ω. Final top air is 12 mm, checked against 30 mm. The current taper has 1,501 output samples, not the older 7,501. Full sweep has 2,000 intervals, not that number of direct solves. Keep historical mesh/reference-plane choices dated. |
| Main 5 | Explicit-capacitor calibration was used for earlier L estimates. Final tightly bounded normal inductances are optimized comparator values, not newly measured inductances after every meander modification. High regression R² does not prove absolute physical accuracy. |
| Main 6 | ABCD and complex-RMS formulas remain valid for the normal ladder, but no longer describe the whole extraction/geometry procedure. The static terminal matrix and distributed seven-port model must be described separately. The final objective is existence of a bounded comparator response, not unique fitted C0/Cg. |
| Main 7–9 | Retain regular/T1/T2 selection numbers as historical staged results. Their 5%/0.025 criteria, frozen older values and source paths are not the current assembled-taper acceptance rule. |
| Main 10 | Most result numbers already agree with evidence. Change “turn count 68 to 72” to “horizontal runs 68 to 72.” Explain that the independent model's RMS and bounded comparator's RMS differ by purpose. The merged taper copy is coverage-equivalent, not separately EM-verified. |
| Main 11 | Completed one-point values agree. Update the snapshot timestamp/count from the mailbox status JSON when editing. Keep full-sweep completion unclaimed until native final output confirms it. |
| Main 12 | The nominal sensitivity numbers concern the historical T2 circuit and are not universal sensitivities of the final T3 fixture. Do not call its 5%/0.025 rule “current.” Explain measured geometry derivatives, held-out checks and conditional terminal-capacitance profiles. Tight bounds are not confidence intervals. |
| Main 13 | Historical runners are useful development provenance but are not the reproduction entry point for the latest frozen geometry. Use the current design registry, `geometry.rebuild_final`, `tools/check_published_files.py`, and `tools/verify_current_response.py`. |
| Appendix E | The sparse four-port isolated-inductance method is distinct from the later seven-port assembled terminal network. Do not claim that reverting to explicit capacitors ruled out all multiport methods. Controlled timing superiority and universal interpolation error remain unproven. |
| Appendix F–I | Extend past the T2 cutoff: T3 pad/IDC revisions, bounded refits, cleanup, full-filter/joint-fit transfer failures, shorting tests, terminal network, geometry derivatives, boundary correction and stricter checkpoints. |
| Appendix J | 400 µm planes verified in current input; an alternative-reference-plane extraction study remains unverified. 12/30 mm top-air convergence does not resolve finer-mesh convergence. |
| Appendix K–L | “T3 still in development” is superseded by the accepted taper checkpoint. Full-band filter validation remains outstanding. Earlier stopped sweeps and current running sweep must not be conflated. |

## NEEDS JUSTIFICATION disposition

- Historical regular-array end-structure mismatch: not isolated by this audit. Keep the old hypothesis qualified; a later accurate different-fixture network does not resolve the original four-versus-six-resonator discrepancy.
- Original four-port calibration timing/accuracy: not a controlled same-geometry comparison established here. Preserve as unresolved, not as proof of superiority of either method.
- Downstream-value uncertainty in the old three-C T2 fit: later bounded refits, shared-fixture failures and terminal compensation tests demonstrate why fixed older values need qualification. They do not provide a universal confidence interval for those three historical values.
- External reference planes: source settings verified at 400 µm. Independence of extracted parameters from all valid plane choices remains unverified. Fitted delay handles smooth phase offsets, not every fixture error.
- Ls units: native surface-impedance conductor has `Ls="0.0388"`; general project inductance units read `NH`. The earlier documentation says pH/square, but that assertion is not confirmed by those XML fields. Keep the exact recorded value and mark the field's unit convention unresolved pending authoritative native UI/help evidence.
- Finer-than-2 µm convergence: unresolved. The available 2/4 µm and subsection screens do not establish continuum convergence.
- Whether the old 5% criterion guarantees a working full filter: unresolved as a general claim, and superseded as the current design criterion. The accepted taper now meets a 0.25%-bounded comparator with RMS below 0.0025. The full-filter single point is encouraging but not a full-band test.

## Evidence and presentation rules

`result.json` explicitly sets `physical_parameter_uniqueness_claim=false`; `local_model/analysis.json` similarly sets uniqueness false. Do not turn the normal comparator table into “measured physical parameters.” Its older `full_filter_EM_verified_for_this_geometry=false` field was written before the subsequent one-point test: qualify it with the newer `full_filter_single_point/complete.json`, while preserving both originals.

The reference-only main draft already includes a timestamped 21-solve snapshot and correct one-point values. Those are not errors; they are observations that can become stale. Use the newer timestamped mailbox snapshot separately instead of rewriting an old count without a time.

The repository's earlier `final_designs/taper_3` and `full_filter` folders remain historical and are not the latest designs. Current published paths are under `final_designs/current_2026_11_08/`. The working original campaign paths are unchanged. SOURCE_PATHS.json provides exact mappings.

The report must not assert an exhaustive causal reconstruction of every abandoned trial merely from a directory name. The supplied chronology/index records saved evidence; cases without a contemporaneous reason should be labeled “reason not recorded/verified” rather than inventing intent.
