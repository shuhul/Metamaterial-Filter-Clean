# Simulation speed: measured changes and limits

The final full-filter run used one workflow-managed Sonnet job with 36 threads. Its resource record reports 36 logical processors and a peak child usage of about 34.87 core equivalents. This establishes the available logical-processor count and observed use, not 36 physical cores. An investigation-local launcher acquired an exclusive serial solver lease and set thread limits for its child. The lease coordinates these workflows, not arbitrary unrelated GUI processes. The user temporarily authorized a separate low-priority short-input test during the long sweep; no claim of a completely idle benchmark machine is made.

## What was actually measured

| Experiment | Timing | Accuracy evidence | Scope |
| --- | --- | --- | --- |
| Four versus eight threads, identical native project bytes, independent cold folders | 1094.40 versus 628.99 seconds; 1.740x wall speedup | Saved dense exports identical at stored precision | Earlier taper full-band ABS case with four direct solves; not the final 16-resonator filter |
| Sixteen versus 36 threads at 4.65 GHz | 95 versus 70 seconds per frequency; about 1.36x | S21 difference 1.08e-12; max all-S difference 1.67e-10 | One-point same-metal pilot; not a full-band speed guarantee |
| Edge refinement removed globally | 811.08 versus 158.36 seconds for three points | S21 RMS 0.247035 | Fast but rejected; cannot claim as usable acceleration |
| Coarse distant ground interiors | Subsections 14,538 to 14,251; first point 258 seconds | First-point S21 error 0.11525 | Rejected; no demonstrated useful timing/accuracy improvement |
| Coil minimum interiors 100 cells, edges retained | Pilot 130 versus 150 seconds per point | Full-band S21 RMS 0.00535839, max pointwise 0.0186366 | Accepted rough screening under the then-authorized 0.01 RMS setting; not final reference settings |

The rough-mesh pilot plus validation cost 692.13 seconds across five native points. A cited approximately 548-second endpoint/adaptive subtotal reused cache and was not an independent cold full-sweep benchmark. The original-mesh eight-thread reference used four points. Those totals must not be used to claim a controlled full-sweep speedup from the rough mesh. Editor polygon merging changed 261 polygons to four but left the estimated 14,538 subsections unchanged; no measured solve-speed benefit was established.

The detailed historical report and controlled 4/8-thread and 36-thread records are packaged in `final_designs/final_full_filter/evidence/`. The historical report records evolving decisions, so its early “no acceleration passed” statement precedes later successful threading and rough-screening tests. It must be read chronologically.

## Final numerical settings

The final run retained 2-by-2 micrometre cells, medium-speed edge meshing, double precision and de-embedding. Its native log reports 41,299 subsections, 13 GB estimated memory and 14 GB peak memory. The complete 4.5-5.5 GHz ABS run solved 123 frequencies in about 41 hours 19 minutes to native completion. The monitor's process-lifetime total is about 41 hours 20 minutes. No equivalent full-filter run at another thread count was performed as a controlled benchmark, so there is no measured final-filter speedup factor.

A large regular ground region need not be meshed as finely as the minimum base cell everywhere: subsectioning already combines cells. Therefore changing a broad interior or merging editor polygons need not reduce matrix size much. Fewer unknowns can help, but counts alone are not timing measurements, and coarse changes must pass response checks.

## ABS, caching and output are different operations

ABS chooses expensive direct frequencies adaptively and constructs a frequency-response interpolant. The 2000 requested output intervals did not mean 2000 field solves. Exporting fewer samples does not by itself eliminate the expensive direct calculations. The original output setting saved 123 direct points, even though the native ABS analysis had converged. A wrapper incorrectly treated that count as failure.

After completion, enabling interpolated output on an isolated project/cache copy produced 2001 samples with `-AbsNoDiscrete`. Its log confirms zero new frequencies, zero errors and zero warnings. That operation reuses a completed calculation; it is not evidence that a fresh full-filter EM calculation can be performed in seconds. Caches depend on compatible geometry/settings; changing metal does not generally make the old EM solution valid for the new geometry.

## Numerical search speedups are not EM speedups

The geometry optimizer uses a local response approximation to evaluate many candidate combinations without a fresh EM solve for every candidate. Analytical derivatives of the circuit and network calculations reduce numerical optimization effort. The early finite-difference geometry runs and the selected candidate's native verification still require EM. No controlled end-to-end runtime comparison established how much the internal-port representation improved search speed, or whether it improved it at all relative to a simpler two-port surrogate.

Sparse frequency checks and small-fixture simulations can screen ideas before a dense export or full-array run, but they do not independently prove full-band or full-filter accuracy. A shorter response interpolation error is not a capacitance uncertainty, and a setting-to-setting numerical RMS is not a circuit-fit RMS.

## Explanations that remain unsubstantiated

- Linear speedup with thread count, or a specific faster final-filter runtime at another count.
- A usable large speedup from the rejected ground-coarsening or global-edge-removal trials.
- An EM speedup from merely reducing plotted/output samples or merging editor polygons.
- A demonstrated computational advantage of internal ports over a direct two-port sensitivity model.
- Finer-mesh or fabricated-device accuracy from ABS convergence alone.
