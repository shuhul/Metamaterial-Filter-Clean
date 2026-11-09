# Full-filter fit investigation: first 72 direct EM points

A stable passive rational two-port response fits the 72 direct Sonnet points at **complex S21 RMS 0.00134715**, with maximum sampled complex error about 0.0035. It passes scikit-rf's global passivity check and has stable poles. This is an accurate response model, not a unique extraction of the original cell capacitances and inductances.

The native full simulation was left running. No project, geometry, calibration, selected design, or solver setting was changed. Only the first 72 completed direct calculations were fitted; no unfinished dense ABS export was used.

## Validation

- Fit all 72 points, then enforce passivity: RMS 0.00134715.
- Refit 58 points, leaving 14 out, then enforce passivity: held-out RMS 0.00279649; training RMS 0.00413794. Frequency-decimated validation is exploratory because the samples are clustered.
- New independent direct point at 4.8428 GHz: actual |S21| 0.986419, predicted 0.985856; complex error 0.00062802. This point was not fitted or used to select the model.

The best response was visually inspected in magnitude, dB, phase, and complex residuals. It tracks the measured passband ripple and skirts. Predictions between samples remain model predictions, not completed EM calculations. Agreement at one subsequent point does not establish accuracy at every uncomputed frequency, nor does it prove mesh convergence.

## Circuit-family comparison

Every circuit test freed its L and C parameters; no original calibration values were overwritten. Resonance-frequency/impedance coordinates reduce the numerical correlation of L and C. Repeated-cell values were initially shared; later tests allowed position dependence and additional parasitic couplings. All modeled capacitors and inductors remained positive; modeled losses remained nonnegative. Nonlocal mutual-L tests enforce a positive-definite inductance matrix. Each circuit also has fitted reference-plane phase terms. These are diagnostic equivalent circuits, not validated component measurements.

| Model | Free parameters | Complex S21 RMS | Magnitude RMS |
|---|---:|---:|---:|
| A_shared_lossless | 15 | 0.0621521 | 0.0466768 |
| B_shared_shunt_loss | 19 | 0.0582671 | 0.0440197 |
| C_shared_loss_feed | 21 | 0.0523781 | 0.0276428 |
| D_shared_loss_ESR_feed | 26 | 0.0474693 | 0.0244561 |
| E_mirror_position_loss | 38 | 0.0148335 | 0.0108324 |
| F_individual_position_loss | 62 | 0.0137718 | 0.0103032 |
| G_shared_nonlocal_passive | 35 | 0.0263533 | 0.0204870 |
| H_position_nonlocal_passive | 47 | 0.0142836 | 0.0108606 |

Finite multistart searches do not prove a global optimum. Some spatially expanded fits reached iteration limits; this is recorded per start. Good response agreement must not be presented as unique physical extraction. Untying repeated cells and adding parasitics changes the model family, rather than establishing that the nominal identical cells physically have those fitted values.

## Why varying only the old LC values can be insufficient

The saved two-port data is not lossless/unitary: the apparent missing port-1 power, computed as `1-|S11|^2-|S21|^2`, has median 3.79% and maximum 29.46%. The metal model has explicit Rdc=Rrf=0, and the project top cover is FREE_SPACE. Sonnet documents the free-space cover as a radiation/absorption boundary. This makes omitted loss/loading a plausible contributor; it does **not** isolate the source of all missing power. No cover-change or mesh comparison simulation was run here.

Official sources: [Sonnet 19.52 analysis box](https://www.sonnetsoftware.com/support/doc_19.52/Topics/Getting_Started/The_Analysis_Box.htm), [Sonnet antennas and radiation](https://www.sonnetsoftware.com/support/doc_19.54/Topics/Antennas_and_Radiation/Antennas_and_Radiation.htm).

The tests show that the dataset admits an accurate stable passive model. They have not established a unique set of the original C0/Cg/L parameters, or that the selected physical geometry now meets the nominal target filter response. A fitted macromodel of the full filter cannot be decomposed into uniquely measured cells, or used unchanged to predict a smaller circuit.

## Saved artifacts and replay

- `results/best_s21_comparison.png`: readable magnitude, dB, and residual comparison.
- `results/rational_passive_fit.png`: response/phase fit inspected during the investigation.
- `results/rational_coefficients.json`: complete 32-state stable passive two-port model parameters.
- `evaluate_model.py`: evaluate those parameters without refitting or simulation.
- `results/passive_model_s21.csv`: predicted response on a 7501-point grid, clearly a model.
- `results/direct_72_points.s2p`: the actual frozen direct samples in an exchange format.
- `results/snapshot.npz`, `em_snapshot.log`, `snapshot_metadata.json`: evidence and provenance.
- `results/passive_rational_holdout.json`, `rational_passivity.json`: validation and rejected model variants.
- `results/summary.json` and per-model JSON/plots: circuit parameter sets and fit diagnostics.
- `campaign.py`, `parasitic_fit.py`, `parasitic_position.py`, `passivity_check.py`, `summarize.py`: analysis code. Initial anchors and prior seeds are frozen alongside the evidence.

Replay the selected response from the project root:

```powershell
uv run python investigations/full_filter_fit_72_2026-11-03/evaluate_model.py
```

The rational equation is a sum of fitted poles/residues plus constant/proportional terms, with conjugate terms for complex poles. Coefficient ordering is documented in the JSON. The coefficient replay was independently checked against scikit-rf's evaluator to numerical precision. The fit band is 4.65-5.40 GHz; out-of-band predictions are not experimentally validated.

## Visual limitations

The flexible component fit predicts a sharp, unverified peak around 5.32 GHz between direct samples. The rational model gives a smooth skirt there. Current data alone does not conclusively resolve every feature between samples. The raw full sweep remains running. The absolute complex-RMS criterion also does not guarantee accurate deep-stopband dB values: the passive model differs substantially in relative terms where the actual transmission is very small. Those differences remain visible in the dB plot.

See `results/best_circuit_parameters.md` for a readable node-by-node table of the best component fit. The original calibrated and nominal parameter sets were left unchanged.
