## 5. Measuring inductance independently

We extract the inductance before fitting the capacitances. A single resonance frequency depends on both L and C, so it cannot determine them separately. Adding several known capacitances gives us a way to measure L from how the resonance moves.

In the calibration fixture, an ideal capacitor, `C_add`, is connected between the resonator node and ground. The metal geometry stays fixed while the added capacitance changes. Under a parallel LC approximation,

```text
omega_res ≈ 1 / sqrt[L × (C_eff + C_add)]

1 / omega_res² ≈ L × C_add + L × C_eff

omega_res = 2 × pi × f_res
```

Plotting `1 / omega_res²` against `C_add` gives an approximately straight line. With capacitance in farads, its slope is inductance in henries. The intercept includes the resonator’s background capacitance and the loading from the calibration fixture.

We do not need to know C₀ or Cg beforehand. Their approximately constant contribution is absorbed into the intercept, allowing L to be found from the slope. This depends on the loading remaining sufficiently constant as the resonance moves. Frequency-dependent coupling or a change in the resonant mode can introduce curvature and bias the extracted inductance.

### Choosing the added capacitances

The original workflow used added capacitances of 20, 30, 40 and 50 fF. These values can still be useful when C₀ is around 150 fF. The added capacitance only needs to produce a clearly measurable frequency shift. For example, ignoring coupling loading, changing `C_add` from 20 to 50 fF with C₀ = 150 fF lowers the resonance frequency by about 8%.

Increasing the loads further gives a larger shift, but also moves the measurement to lower frequencies. The useful range depends on the fixture and the available simulation bandwidth. We want well-separated resonances that remain within the calibration range and follow the same mode.

### The current calibration workflow

The current taper workflow uses a separate isolation fixture containing the resonator being measured and calibration feed electrodes. The resonator body matches the corresponding response geometry. Two external ports provide the input and output connections, and two additional component ports provide the connection where `C_add` will be applied.

Sonnet calculates this four-port network directly at 1, 2, 2.5, 3 and 4 GHz. The exported network contains the response needed to connect different ideal capacitors afterward in Python. Each capacitor value therefore uses the same EM solution.

The code converts the network into an admittance representation and constructs a frequency-dependent model using the results at 1, 2, 3 and 4 GHz. The directly simulated point at 2.5 GHz is held out and used to check how well the model predicts an intermediate frequency. A separate numerical check confirms that two ways of connecting the capacitor give the same response.

The capacitor values are then selected automatically. The code calculates positive loads intended to place the loaded resonances near 3.8, 3.4, 3.0 and 2.6 GHz. This keeps the measurement within the simulated frequency range and gives useful separation between the resonances without requiring prior fitted values of C₀ and Cg.

For each load, Python calculates the loaded transmission response, locates its maximum and refines the peak frequency numerically. It then fits the line relating `1 / omega_res²` to `C_add`.

For the accepted current T1, the saved calibration used the following values.

| Added capacitance | Resonance frequency |
| --- | --- |
| 128.57 fF | 3.8 GHz |
| 201.60 fF | 3.4 GHz |
| 305.42 fF | 3.0 GHz |
| 460.43 fF | 2.6 GHz |

The resulting inductance was 6.00627 nH, with a regression R² of approximately 0.9999988. The inductances calculated from adjacent pairs of loads ranged from about 5.984 to 6.017 nH, providing an additional check on how much the slope changes across the calibration range.

### Checks and interpretation

The workflow rejects loads that are nonpositive or fail to increase, peaks at the search boundaries, resonance frequencies that fail to decrease with increasing capacitance, and excessive error at the held-out frequency. Geometry and file hashes are also checked before reusing saved calibration data.

Once L has been extracted, it is held fixed during the capacitance fit to the assembled fixture. This separates the inductance calibration from the later capacitance extraction. The result is an effective inductance for the isolation fixture and calibration frequency range. Its accuracy still depends on the mesh, material model, fixture loading and validity of the parallel LC approximation. A high R² supports linearity over the chosen load range, but does not quantify all of these errors.
