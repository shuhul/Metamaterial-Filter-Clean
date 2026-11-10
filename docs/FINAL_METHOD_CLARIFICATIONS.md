# Sections 10 and 11: final method clarification

## Suggested main-report explanation

The initial objective was to turn a target LC ladder into metal whose simulated transmission has the intended filter response. Earlier stages calibrated inductance and fitted a few capacitances to two-port transmission. As the assembly grew, low residuals did not establish unique physical component values, and some fitted circuits did not transfer reliably to changed fixtures.

The later implementation therefore built a numerical model of the simulated geometry and measured how small geometry changes changed its response. Its chosen representation was a seven-port network: two external ports plus five delta-gap ports at coil returns. Sonnet supplied the full matrix of responses between these ports. Closing the five internal interfaces mathematically recovered the ordinary two-port response, which was checked against a separate uninstrumented simulation.

Two descriptions were calculated from that data. A low-frequency fit provided a diagnostic terminal-capacitance representation; it was not the final optimization objective. A frequency-dependent model near 5 GHz predicted the geometry's response. The optimizer used changes in that full network under controlled geometry perturbations to predict combinations of pad, rail, IDC and coil adjustments without running EM for every mathematical candidate. Each predicted network was reduced to external S21 and compared with an ideal ladder whose component values were constrained near the requested targets. Promising realizable dimensions were then generated and tested with ordinary two-port Sonnet simulations.

The useful general strategy was measured geometry sensitivity, inexpensive local prediction and native verification. Internal ports were used in this implementation, but their advantage was not established. A simpler approach could obtain sensitivities directly from ordinary two-port S matrices, score the predicted S21 and verify candidates in exactly the same manner. No comparative test showed that the internal representation was faster, more accurate or necessary. It also did not uniquely extract the original C0, Cg and L values.

The accepted taper had native complex-S21 RMS 0.00249186334 against a ladder with all ten capacitances and five inductances within 0.25% of the nominal values. This means the response matches a tightly constrained target-like circuit. Those bounds are not physical-parameter measurement uncertainties. Separately, the terminal model reproduced the existing taper response with RMS 8.11746e-6. The full sixteen-resonator filter was subsequently simulated without the internal ports and has its own completed response evidence.

## Appendix: what the evidence establishes

- Polygon union checks show unchanged metal coverage when preparing the internal interfaces. They do not alone prove identical electrical response.
- Short closure is a Schur-complement elimination of internal currents under zero interface voltage. Its direct S and Z formulations agree, checking implementation. The original two-port comparison is the separate physical restoration check; maximum all-S difference was approximately 0.00011846.
- Added ports, calibration and frequency interpolation can create discrepancies. Exact closure algebra does not remove inaccurate input data. Final candidate acceptance used uninstrumented native responses.
- The low-frequency representation selects ports 1,3,4,5,6,2 from polarity-adjusted Z and omits port 7. In that Z relationship the omitted current is zero. The selected coordinates are not automatically common-ground circuit-node voltages; row sums cannot simply be called uniquely measured pad capacitances.
- The five-frequency operating-band expansion interpolates the training matrices. Held-out 5.096 GHz and native two-port comparisons check its response accuracy over the tested band. This does not establish broadband passivity/causality or optimality of the chosen basis.
- Candidate prediction uses the complete operating-band matrix changes. It is not a step of diagnosing “one C0 is wrong, therefore change this pad.” The low-frequency C matrix did not drive the final objective.
- The raised-boundary refinement initially reused older geometry sensitivities. It checked predictions, updated the model and restricted steps; it did not newly measure every derivative at every candidate.
- The final +32 to +34 micrometre pad-control step changed that geometry control only between those two candidates. The target comparator could be refitted within the permitted bounds. Earlier changes to meanders and other dimensions remained in the accepted baseline.

## NEEDS JUSTIFICATION if claimed

| Potential claim | Current evidence and proper qualification |
| --- | --- |
| Internal ports were required or improved the optimization | Not demonstrated. Requires a controlled comparison with a direct two-port sensitivity optimizer. |
| Original physical C0/Cg/L are uniquely known to 0.25% | Not established. Bounds are imposed on comparator values; no unique physical decomposition follows. |
| Six-terminal row sums equal physical pad-to-ground capacitances | Requires a justified mapping from selected differential port coordinates and termination conditions to physical node voltages. Do not assert automatically. |
| Odd-power frequency basis is uniquely best or globally physical | No comparison/guarantee established; evidence is local held-out response agreement. |
| Final L values are new independent calibrated inductances | They are bounded comparator parameters, not fresh calibrations of every final coil. |
| ABS 40 dB proves mesh or fabricated-device accuracy | It assesses interpolation convergence, not these other uncertainties. |
| A measured EM threading speedup guarantees the same ratio for the full filter | Earlier controlled trials have narrower geometry/frequency scope; final run has no matched lower-thread full-array benchmark. |

These are limits on claims, not assertions that the achieved native full-filter response is invalid. The completed filter evidence stands separately from whether the chosen optimization machinery was the simplest possible route.

## Screenshot captions

Use `final_designs/final_T3/port_review/gap-ports-overview.png` to illustrate temporary ports 3-7 in the small taper fixture. Use `gap-port-detail.png` to show port 6 at the regular coil return. The visible negative mark above and positive mark below indicate local reference orientation, not actual charges or a positively biased ground. Neither image is the full ten-unit-cell filter. Detailed captions and the native location are recorded in the adjacent README.
