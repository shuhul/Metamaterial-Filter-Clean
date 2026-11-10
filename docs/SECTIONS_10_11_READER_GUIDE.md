# Reader guide: the assembled taper and the full filter

This guide explains report Section 10 (Building an Independent Network Model) and Section 11 (Matching the Target Circuit and Updating Geometry). Full-filter verification belongs to report Section 13. Chapter numbers below are local to this guide. It starts with the purpose, then follows the actual calculation from Sonnet geometry to the accepted taper and the full-filter test. It separates measurements made from simulations, numerical fits, design comparisons, and what remains unverified. No new simulation was run to write this guide.

**Implementation clarification:** the optimizer used the internal seven-port network to predict geometry changes, but its objective was external S21 agreement with a bounded target ladder. It did not use uniquely extracted internal C0/Cg/L values to determine those changes. The same sensitivity-based search could be built from ordinary two-port responses, without internal ports or short restoration. No controlled comparison demonstrated that the internal-port approach was necessary, faster or more accurate. It is the implementation we used, not a proven advantage over that simpler alternative. Final acceptance used ordinary two-port EM, and the full filter was simulated without internal ports.

## 1. What we are trying to build

The verified project goal from the design discussions is a manufacturable superconducting microwave filter with a passband around **4.8–5.2 GHz**, connected to 50-ohm ports. Ten repeated regular resonators form the middle. Three taper resonators on each side transition between the ports and the regular chain:

`port — T3 — T2 — T1 — 10 regular cells — T1 — T2 — T3 — port`

There are sixteen resonators in the full filter. The taper changes the resonator and coupling values near the ends to reduce mismatch. Good matching permits energy to enter the regular array instead of reflecting back toward the source. The regular-cell electrical scale puts the desired collective band near 5 GHz, but the coupled array's band edges cannot be calculated from an isolated LC frequency alone.

The user supplied the target circuit values and Frank's paper as a reference, and repeatedly identified the 4.8–5.2 GHz filter response as the practical objective. The material available here does **not** establish why that band was selected for a particular qubit, transition, readout channel, experiment, or fabrication package. Those possible applications must not be invented. Likewise, the 50-ohm port impedance is the wave normalization/termination, not a promise that every piece of physical transmission line has exactly 50-ohm characteristic impedance.

The supplied AWR inputs are preserved unchanged in `reference_inputs/awr/`. The screenshot shows an enabled goal for the magnitude of Full_Filter S(1,1), measurement less than 0.001, from **4.82 to 5.18 GHz**, with Difference cost, default L=2 and weight 1. That is a reflection objective in the interior of the requested band, leaving 20 MHz at each band edge. It is not a screenshot of an achieved EM response. The user reports that minimizing S11 worked better than directly maximizing S21 in their AWR optimization; that experience is recorded without claiming universal optimizer superiority or that this goal was already used by our Sonnet geometry search.

For a lossless two-port with the second port matched, energy conservation gives

`|S11|² + |S21|² = 1`.

Therefore, in that ideal case an S11 magnitude below 0.001 corresponds to reflected power below one part in a million and S21 magnitude above about 0.9999995. A real EM model can also dissipate or radiate power, so low reflection alone is not sufficient: a matched absorber can have low S11 and low S21. We must still inspect transmission and energy loss. The threshold is a linear magnitude; 0.001 corresponds to −60 dB, not −0.001 dB.

## 2. What Section 10 actually simulates

The accepted taper fixture is smaller than the full filter. It contains T3, T2, T1, a regular resonator, and an output-end resonator/fixture connected to the second port. Including nearby metal captures the loading and coupling that would be missed by isolated component tests. This fixture is useful for faster design iterations, but its two-port S21 is not the full sixteen-resonator filter's S21.

The canonical input is `final_designs/final_T3/project.sonx`; its completed response is `project.snp`. The continuous-ground drawing is `continuous_ground/project.sonx`. The original and merged drawings have identical metal coverage, but the merged taper drawing is not a separate new EM validation. The full-filter geometry uses merged ground and has its own EM calculation.

Open `s21_with_target_parameters.png` to see the native response and the bounded circuit comparison. Open `s21_verified.png` to see the independent network prediction as well. These are three distinct curves: the geometry simulated by Sonnet, a model of that geometry, and a nearby target circuit.

## 3. What Sonnet and S-parameters tell us

Sonnet solves the electromagnetic problem for the metal, dielectric stack, enclosure/boundaries and ports. An S-parameter describes outgoing waves produced by an incoming wave, with other ports terminated in their reference impedances. S21 is transmission from port 1 to port 2; S11 is reflection at port 1. These are complex numbers: magnitude measures amplitude, phase records timing/phase shift.

Two external ports provide only a limited view of many coupled internal structures. Several different allocations of capacitance and inductance can reproduce the same external response. An excellent fit is useful, but it does not by itself tell us which individual physical electrode has a uniquely measured C0 or Cg.

Before this method, the regular/T1/T2 workflow calibrated L with known capacitors and fitted a small LC ladder. For the standalone T2 selection, three new capacitances and a fixture delay varied while older values and L were fixed. This is a conditional extraction: the returned T2 values depend on the fixed downstream assumptions. Those historical selections remain under `final_regular`, `final_T1` and `final_T2`; the later embedded T1/T2 geometry in `final_T3` is not identical to every standalone baseline.

## 4. Why add internal ports, and what does “shorting” mean here?

### What a GAP port actually is

In our Sonnet files, **GAP means a delta-gap port**: an internal, two-terminal electrical access point placed where two metal polygons meet. One terminal belongs to each polygon. Sonnet can excite that connection and calculate its voltage/current relationship to every other port. Unlike the external input and output ports, this access point sits inside the circuit.

Despite the name, it does **not require a visible physical gap**. The two polygons abut: they meet along an edge. The port introduces an electrical interface at that edge. Its voltage is the difference between its two terminals, not necessarily a voltage relative to a universal ground. Polarity therefore matters. See [Sonnet's delta-gap port documentation](https://www.sonnetsoftware.com/support/doc_19.52/Topics/Ports/Delta_Gap_Ports.htm).

For our coil-return connections, the idea is:

```text
Original connection:        meander ---- continuous metal ---- return metal
Instrumented connection:    meander -------- [ + port - ] ----- return metal
Restored in Python:         meander -------- [  short  ] ----- return metal
```

The drawing is schematic: the port interface is at the shared polygon edge, not a new finite-length device. These are temporary simulation access points, not extra connectors or components to fabricate.

### Where we put them and what information they provide

The current method adds **five internal GAP ports at existing coil-return connections**, retaining external ports 1 and 2. Internal ports 3, 4, 5, 6 and 7 correspond to T3, T2, T1, the regular resonator and the output-end resonator, respectively. This produces a seven-port network. We split the metal polygons at the chosen connections while preserving their total metal coverage. A Boolean XOR check verifies that the instrumented and original layouts cover the same area; it does not, by itself, prove that adding ports has no electrical effect.

Sonnet calculates a 7-by-7 S-parameter matrix at each frequency. Conceptually, it excites one port at a time and records the waves at all seven ports, with the other ports matched to their reference impedances. This gives independent access to internal connections as well as the input and output. We can then convert that matrix into relationships between terminal voltages and currents. Those additional relationships constrain the terminal-network model more than external S21 alone can.

A port does not directly report “this capacitor is 150 fF.” Capacitances and frequency-dependent network coefficients are inferred from these relationships in the following sections. More access improves the available information, but does not automatically make the original ladder's C0, Cg and L uniquely identifiable.

### Why we must close the internal ports again

These ports are **not automatically nonloading voltage probes**. Reading the raw seven-port S21 with the five extra ports matched to 50 ohms describes a circuit with those additional terminations. That is not the original continuous-metal filter. The 50-ohm reference used to express S-parameters does not mean we want five physical 50-ohm resistors in the final device.

Instead, Python uses the complete multiport response to reconnect each internal interface with an ideal short. **The intention is to make that interface behave as it did before the port was inserted: a continuous electrical connection.** A short means **zero voltage across that interface, while allowing current to flow**. It does not mean zero current, zero voltage everywhere in the resonator, or removal of the meander. The inductor is still in the circuit; only its original return connection has been restored. Python calculates this reconnection from the saved network; it does not need to run another EM simulation to impose it.

This is why the workflow has two separate steps: first obtain the richer seven-port description, then calculate the external response with the internal connections closed. The resulting two-port response can be compared directly with the original, uninstrumented Sonnet simulation.

This differs from the earlier physical bridges across series capacitors. Those bridges changed the connection topology and left the earlier shunt inductors present. They did not remove unwanted resonators. Here the internal-port network is an instrumented representation, and its restoration is tested against the original uninstrumented two-port calculation.

### How the mathematical short works

Start with Sonnet's **7-by-7 S matrix at a particular frequency**. The following calculation is repeated at each frequency.

**Step 1: Convert the S matrix to an impedance matrix, Z.** S relates incoming and outgoing waves. Z relates terminal voltages and currents: `V = Z I`. The conversion makes a short easy to express as a zero-voltage condition. For our equal, real 50-ohm port references:

`Z = 50 (I + S) inv(I − S)`.

Here `I` without a subscript is the identity matrix, not a current. The resulting Z matrix is still 7-by-7; no ports have been removed yet.

**Step 2: Separate external and internal variables.** Use `e` for the two input/output ports and `i` for the five added ports. The matrix equation becomes:

`Ve = Zee Ie + Zei Ii`

`Vi = Zie Ie + Zii Ii`.

`Ve` and `Ie` each contain two values; `Vi` and `Ii` each contain five. The blocks have these sizes:

| Block | Size | What it describes |
| --- | --- | --- |
| `Zee` | 2-by-2 | External voltages due to external currents |
| `Zei` | 2-by-5 | External voltages due to internal currents |
| `Zie` | 5-by-2 | Internal voltages due to external currents |
| `Zii` | 5-by-5 | Internal voltages due to internal currents |

**Step 3: Impose the five shorts.** Require `Vi = 0`, meaning zero voltage difference across each internal interface. The second equation then determines the internal currents:

`Ii = −inv(Zii) Zie Ie`.

These currents are generally nonzero. They are exactly the currents required to satisfy the closed-connection condition for a given external excitation.

**Step 4: Include those currents in the external response.** Substitute that expression for `Ii` into the first equation:

`Ve = (Zee − Zei inv(Zii) Zie) Ie`.

Therefore the effective two-port impedance is

`Zrestored = Zee − Zei inv(Zii) Zie`.

**This is the Schur complement: the effective matrix that remains after eliminating one set of variables from a linear system.** Here we eliminate the five internal currents as independent variables, while retaining their effects on the external ports. We do not eliminate the physical resonators. `Zrestored` is now 2-by-2.

**Step 5: Convert the effective Z matrix back to S.**

`Srestored = (Zrestored − 50 I) inv(Zrestored + 50 I)`.

This 2-by-2 matrix gives S11, S21, S12 and S22 with all five internal connections restored. The inverse notation expresses the algebra; a numerical implementation can solve the corresponding linear systems without explicitly forming an inverse.

Simply deleting the five internal rows and columns from the original S matrix would not work. The remaining `See` block describes the external response with the internal ports matched to their reference impedances, rather than shorted. The coupling terms in the reduction account for the change in termination.

**Illustrative arithmetic, not filter data:** with one external and one internal variable, suppose Zee=100 ohms, Zei=Zie=20 ohms and Zii=10 ohms. For Ie=1 ampere, the short condition gives `0 = 20*1 + 10*Ii`, so Ii=-2 amperes. Then Ve=100*1+20*(-2)=60 volts. Equivalently the effective impedance is `100 - 20*(1/10)*20 = 60 ohms`. The negative current simply flows opposite the chosen reference direction. Keeping only Zee would incorrectly give 100 ohms.

### How we checked that this represents the original circuit

There are three distinct checks; they answer different questions.

**An equivalent calculation directly in S checks the implementation.** We can impose the same shorts using incoming waves `a` and outgoing waves `b`, without converting to Z. At each ideal short, the returning wave has the opposite voltage amplitude: `ai = −bi`. Starting with `b = S a` and eliminating the internal waves gives

`Srestored = See − Sei inv(I + Sii) Sie`.

This should agree with the Z-based calculation above. This is what “equivalent direct S-matrix reduction” means: two algebraic routes to the same terminated network. Their agreement checks the conversion, indexing and reduction code; it is not a second independent EM validation.

**The port-polarity convention keeps the signs consistent.** Each delta-gap port has positive and negative terminals. Its voltage is measured from negative to positive, with a corresponding current direction. The implementation uses the recorded orientations consistently when representing and comparing the terminal network. Reversing a port's orientation requires corresponding changes to both its voltage/current definitions and the matrix entries that couple it to other ports. Polarity is bookkeeping, not a freely fitted parameter. Although a short always has zero voltage under either orientation, inconsistent signs in the rest of the network would corrupt the calculations.

**Comparison with the original two-port EM simulation checks physical restoration.** The decisive check is against an independent Sonnet simulation of the original connected metal without the five added ports. For the accepted geometry, restored-versus-original maximum error across all four external S entries is about **0.00011846**, below its 0.001 restoration gate. See `final_designs/final_T3/local_model/analysis.json` and `manifest.json`.

That agreement supports using the instrumented network to describe this geometry under the tested conditions. It does not establish perfect calibration at every frequency or unique physical values for every original lumped component. Those are separate questions from whether closing the added ports recovers the original external response.

At this point, we have added five internal access points, obtained the complete response between all seven ports, and checked that electrically closing the added interfaces reproduces the original two-port circuit. This is more information than five local voltage readings: it records how excitation at each port affects every other port. We retain the seven-port data for the extraction below; the restored two-port response is its validation and the externally observable prediction. No individual C0 or Cg has been extracted merely by performing this reduction.

## 5. What is extracted from the seven-port calculation?

### Why do another fit after obtaining the seven-port response?

Section 4 gives us a matrix at each simulated frequency: a detailed electrical description of the geometry. It does not directly label particular matrix entries as the original circuit's C0 or Cg. We use this data in two different ways:

1. **Interpret the capacitive behavior:** use low-frequency data to estimate an effective terminal capacitance matrix. This is a diagnostic description of loading and coupling.
2. **Predict the response where the filter operates:** build a frequency-dependent terminal-network model around 5 GHz. This is the model used to predict responses and assess small geometry changes.

**The low-frequency fit is not an attempt to reduce the stopband.** Neither fit changes the metal or improves transmission by itself. They describe the existing simulated geometry. Improvement requires choosing a physical geometry change and verifying its response, as described in Section 7. The low-frequency capacitance extraction is not a prerequisite for constructing the operating-band model.

### A. Low-frequency fitting: interpreting capacitance

An ideal capacitor has impedance `ZC = 1/(j 2π f C)`, whose imaginary part is `−1/(2π f C)`. An ideal inductor has impedance `ZL = j 2π f L`. Thus capacitive and inductive contributions have different frequency dependencies: inverse frequency versus frequency. We use that distinction to help separate the network's capacitive behavior from other contributions.

For the selected six-terminal representation, the low-frequency calculation fits

`Im Z(x) = −A/x + B x + D x³`, with `x = f/(1 GHz)`.

**What is x?** It is just frequency divided by a reference frequency, so it has no units. At 0.5 GHz, x=0.5; at 1 GHz, x=1; at 2 GHz, x=2. Earlier versions of this guide called this variable r. It is not resistance, a circuit element, or another unknown to fit. Normalizing frequency keeps the numerical values convenient.

`Im Z` means the imaginary part of the impedance matrix. A, B and D are coefficient matrices, not three individual physical components:

- **A** describes the inverse-frequency capacitive contribution.
- **B** captures the linear-frequency contribution, including inductive behavior.
- **D** supplies a higher-order frequency-dependent correction.

The code applies polarity signs `[1,1,-1,-1,-1,-1,-1]`, converts seven-port S to Z, then selects the principal Z submatrix with zero-based indices `[0,2,3,4,5,1]`: ports 1, 3, 4, 5, 6 and 2, omitting port 7. In Z-parameter terms, this selected relationship assumes zero current at the omitted port. It is not the five-port short closure used for final transmission. The saved six-terminal labels are input, T3, T2, T1, regular and output end. These are labels for selected terminal coordinates; they do not prove that every voltage is a node voltage relative to a common physical ground. The coefficients describe this representation rather than isolated physical capacitors.

We use results at **0.5, 1 and 2 GHz** to determine the coefficients. In this model and normalization, the effective capacitance matrix in femtofarads is

`C = inv(A) × 10^6/(2π)`.

The matrix inverse is needed because capacitive impedance is proportional to the inverse capacitance matrix; taking the reciprocal of each entry separately would not give the same result. The numerical factor converts the 1 GHz normalization and farads into femtofarads. This interpretation depends on the chosen terminal representation and the adequacy of the frequency expansion.

A separate **1.5 GHz** result is held out: it is not used to determine A, B or D, but tests whether the fitted model predicts additional data. The accepted held-out relative impedance error is **8.25e−6**. A two-frequency extraction also checks sensitivity to the extrapolation. These are checks of the representation, not proof of uniquely identified original ladder components.

In a common-reference nodal capacitance matrix, diagonal entries include the total attached capacitance, negative off-diagonal entries describe mutual capacitance, and row sums give capacitance to the reference. The saved calculation uses this interpretation and records row sums. However, the selected GAP-port coordinates require care: a row sum is not by itself a uniquely measured pad-to-ground capacitance. These diagnostic values cannot simply replace the ladder C0 or Cg values, and they were not the objective in the final geometry search.

### B. Operating-band fitting: predicting the response near 5 GHz

A static capacitance matrix alone does not capture the complete behavior near the filter band. The metal has inductive, capacitive and distributed electromagnetic effects. For accurate response prediction, we therefore fit the **entire complex 7-by-7 Z matrix** near the operating band, using the frequency functions

`u^-1, u, u³, u⁵, u⁷`, with `u = f/(5 GHz)`.

Here u is another dimensionless frequency variable: at 5 GHz, u=1; at 4.8 GHz, u=0.96. Earlier text used r for both normalizations; this guide now uses x and u to distinguish them. These powers form a numerical approximation of the frequency dependence, not a list of five identifiable physical element types.

Each frequency function has a complex 7-by-7 coefficient matrix. The code solves five linear equations for the five coefficients of each matrix entry. This interpolates the training data; the held-out checks, rather than the training agreement itself, test prediction. The inverse-frequency term accommodates capacitive behavior and the odd powers supply additional smooth frequency dependence. No controlled comparison established this basis as uniquely optimal. Here ?independent? means constructed without forcing the target ladder values; it still depends on Sonnet, the geometry, port definitions and numerical approximation.

Training frequencies are **4.65, 4.861, 5.0, 5.25 and 5.4 GHz**. Once the coefficients are determined, the model can evaluate Z at other frequencies in the band. We then apply the internal-port shorting calculation from Section 4 and convert back to two-port S to predict transmission and reflection.

The separate **5.096 GHz** native calculation checks a point that did not determine the coefficients. Its restored S21 error is **1.45418e−6**. Across the dense two-port export, S21 RMS is **8.11746e−6** and all-S RMS is **3.44610e−5**. These errors measure prediction of the existing geometry, not its agreement with the desired filter target.

This operating-band model retains distributed effects that a single static LC per cell can miss. It is validated for the stated geometry and band; it is not asserted to be a globally passive/causal model for all frequencies. The 1,501 samples in the comparison are dense exported values, not 1,501 independent native EM calculations.

The practical motivation is to obtain a trustworthy description of the geometry before asking how to change it. The low-frequency matrix helps interpret its capacitances; the operating-band network predicts the response that matters for the filter. Section 7 then measures how physical dimensions change that network and uses those measured sensitivities to propose geometry updates.

## 6. Connecting the terminal network to the target circuit (report Section 11)

After Section 5 we have a model that predicts the existing geometry's electrical response near 5 GHz. We also have your desired circuit, specified by its coupling capacitors Cg, shunt capacitors C0 and inductors L. These are two descriptions that we can compare through their predicted S21:

```text
Geometry -> seven-port EM data -> operating-band Z model
         -> close internal ports -> predicted geometry S21
                                            |
                                      compare the curves
                                            |
Target Cg, C0, L -> ideal ladder calculation -> target-circuit S21
```

### Do the internal ports directly give each coupling and tank impedance?

**No.** The terminal network describes the coupled geometry as a whole. It includes mutual coupling and distributed loading, so a particular Z-matrix entry is not automatically the impedance of one coupling capacitor or one shunt tank. For example, Z at one terminal depends on the rest of the network and the specified terminal conditions. The port positions also do not isolate every capacitor as a separate two-terminal component.

The low-frequency capacitance matrix helps interpret that loading, but its entries cannot simply replace the C0 and Cg values in the original ladder. The operating-band fit similarly gives a reliable network response, not a unique decomposition into the original ladder components. This is the distinction that connects Section 5 to the design procedure: **we can predict how the geometry behaves without first assigning an independently measured C0, Cg and L to every feature.**

### How the target circuit supplies the comparison

For the ideal ladder, the component definitions are known. At angular frequency `omega = 2*pi*f`:

- A coupling capacitor has impedance `Zcoupling = 1/(j*omega*Cg)`.
- A parallel shunt capacitor and inductor have admittance `Ytank = j*omega*C0 + 1/(j*omega*L)`; the tank impedance is `1/Ytank`.

Here j is the imaginary unit and admittance means inverse impedance. These equations calculate the intended component impedances from the chosen ideal values; they do not extract them directly from a Sonnet matrix.

To calculate the entire ideal circuit, we combine the series and shunt sections using their ABCD matrices, `[1, Zcoupling; 0, 1]` and `[1, 0; Ytank, 1]`. Multiplying them in circuit order gives A, B, C and D for the complete ladder. With 50-ohm external ports:

`S21_ladder = 2/(A + B/50 + 50*C + D)`.

We now have two curves at matching frequencies: the restored terminal-model response of the geometry and the ladder response. Comparing them tells us how close that geometry is to the desired behavior.

### What the final bounded ladder fit means

The final comparison allows ten capacitances and five inductances in the ladder to move within **+/-0.25% of their nominal targets**. It also allows a total reference delay within +/-100 ps. The accepted delay is 3.58105 ps. Multiplying by `exp(-j*2*pi*f*delay)` adjusts phase; it cannot change the transmission magnitude or resonance positions.

The saved parameter ordering is explicit:

| Entries | Order | Nominal values |
| --- | --- | --- |
| C vector 0-4 (fF) | T3 C0, T2 C0, T1 C0, regular C0, output-end C0 | 251, 149, 148, 190.4, 190.4 |
| C vector 5-9 (fF) | input-T3, T3-T2, T2-T1, T1-regular, regular-output | 261, 153, 15.3, 8.2865, 8.2865 |
| L vector 0-4 (nH) | T3, T2, T1, regular, output end | 2.77, 4.33, 5.97, 4.9194, 4.9194 |

The optimizer's vector holds these 15 component multipliers followed by total delay in picoseconds. Their accepted values are in `final_designs/final_T3/result.json`. The output-end fixture is not an extra regular cell to retain in the final full-filter prediction.

This defines a small family of acceptable target-like circuit responses. For a given geometry, the fit finds which member of that family best matches its response. It does not freely move the ideal circuit until any geometry appears successful.

The comparison error is

`R = sqrt(mean(abs(S21_geometry - S21_ladder*exp(-j*2*pi*f*delay))**2))`.

The numerical least-squares residual concatenates the real and imaginary errors. If there are N frequency points and its cost is half the sum of squared residuals, the reported RMS is `sqrt(2*cost/N)`, not division by 2N. As an illustrative example, errors 0.01 and 0.03j at two frequencies give RMS `sqrt((0.01**2+0.03**2)/2) = 0.02236`. It uses complex S21, equally weighted over the comparison frequencies. The final native taper result gives **R = 0.00249186334**, below 0.0025, while the comparator's components remain within their allowed bounds. This establishes agreement with a target-like circuit; it does not prove that every physical component was uniquely extracted to 0.25%. In particular, the final comparator inductances are not fresh independent calibrations of every final meander.

There are therefore two different errors to keep separate:

| Question | Comparison | Accepted result |
| --- | --- | --- |
| Does our terminal model predict the existing geometry accurately? | Model S21 versus native Sonnet S21 | RMS about 0.00000812 |
| Does the geometry behave like an allowed target circuit? | Native Sonnet S21 versus bounded ladder S21 | RMS about 0.002492 |

The first makes the model useful for prediction. The second is the design objective. The next section explains how we change the geometry to improve the second while checking that the first remains reliable.

## 7. From a response model to an actual geometry change (report Section 11)

### The missing information: how does changing a dimension change the network?

The baseline frequency fits tell us how **one geometry** behaves as frequency changes. They do not, on their own, tell us what happens when a pad is widened. To learn that, we need additional native simulations of deliberately changed geometries.

For example, keep the port definitions and all other geometry controls fixed, increase one pad width by a small known amount, and calculate the changed seven-port network. Compare it with the baseline at matching frequencies. This measures the effect of that physical change on the complete network, including changes to coupling, shunt loading and distributed effects together.

### A sensitivity is simply a change per unit of geometry

Use a descriptive name rather than the earlier confusing `dZ/ddk` notation:

`pad_width_sensitivity(f) = [Z_wider_pad(f) - Z_baseline(f)] / width_step`.

This is a matrix at each frequency. If width_step is measured in micrometres, each entry has units of ohms per micrometre. It describes how much that entry changes for a small increase in pad width.

For illustration only, suppose the baseline pad is 100 micrometres wide and the test pad is 102 micrometres wide. The measured sensitivity is `(Z_102 - Z_100)/2`. A proposed 1 micrometre increase would then predict approximately half the observed matrix change. This numerical example explains the calculation; it is not a recorded design or a permissible off-grid fabrication proposal.

An actual early test in `reachable_geometry.py` measured input IDC span +2 micrometres and T3 pad width +2 micrometres separately, using forward differences from the baseline divided by 2. It saved a prediction for the combined +2/-2 micrometre change before the candidate EM result. Later search models reused and updated derivatives, rather than independently remeasuring every control at every iteration. Central finite-difference checks in `quarter_analytic_response.py` validated analytical derivatives of the mathematical response calculation; those checks were not new plus/minus geometry EM experiments.

We measure analogous sensitivities for controls such as rail spacing, finger span and coil geometry. A discrete control, such as meander run count, uses the measured effect of a finite step rather than an infinitesimal derivative. The sparse frequency models help represent these changed networks over the operating band, with their own accuracy checks.

### Predict candidate geometries without solving every candidate in Sonnet

For small changes, add the measured contributions:

```text
Z_candidate(f) approximately equals
    Z_baseline(f)
  + pad_width_change * pad_width_sensitivity(f)
  + rail_spacing_change * rail_spacing_sensitivity(f)
  + finger_span_change * finger_span_sensitivity(f)
  + ...
```

All terms refer to compatible terminal definitions and the same frequency. This is a local, first-order approximation. It may miss nonlinear changes or interactions between simultaneous geometry adjustments, so the allowed steps are kept small.

For each proposed set of dimensions, we then perform the calculations already explained:

1. Evaluate the predicted seven-port Z matrix across the operating band.
2. Close the internal ports using the Section 4 reduction.
3. Convert the resulting two-port Z matrix to S21.
4. Compare that predicted S21 with the bounded target-circuit family from Section 6.

This is how a change in a pad width becomes a predicted change in filter response. We do not have to guess that widening the pad changes only C0, nor force the complicated network into a unique ladder decomposition first. The low-frequency capacitance interpretation can suggest which controls to investigate, but the measured operating-band network changes determine the predicted response.

### Choose a candidate, then check it with a native simulation

The optimizer searches small allowed geometry changes and the bounded ladder values together, looking for lower comparison error. Geometry controls remain subject to the existing topology and the 2 micrometre simulation/design grid. The procedure is:

1. Validate a baseline network against the native two-port response.
2. Simulate controlled geometry perturbations and measure their sensitivities.
3. Search small candidate changes using the inexpensive network predictions.
4. Round or enumerate realizable choices on the 2 micrometre grid, and save the chosen candidate's predicted response.
5. Generate that candidate's actual geometry and run Sonnet.
6. Compare the actual response both with the saved prediction and with the target-like ladder response.
7. Accept actual improvement. If prediction fails or the candidate is worse, reject it, shrink the step range or remeasure the sensitivities. Re-establish the model around an updated geometry as needed.

The optimizer therefore produces an explicit change to a geometry control, not just a desired change in an abstract capacitance. Its proposal remains provisional until the native simulation checks it.

Sparse difference models have separate held-out/dense response checks. The recorded 0.0003 all-S threshold concerns the **49 complex scattering entries** of a seven-port response, not 49 physical fitted parameters. Baseline-model accuracy and restoration of the original circuit have separate checks.

### Were internal ports necessary for this search?

No necessity or advantage was demonstrated. They were genuinely used: `refine_quarter_coils.py::propose_network` forms a candidate seven-port Z matrix from a baseline plus measured matrix changes, calls the short-closure response calculation, and scores external S21 against the bounded ladder. `refine_quarter_raised_cover.py` verifies proposed geometries with native two-port Sonnet responses and updates the network model. There is no additional step in this optimizer where a uniquely extracted internal capacitor identifies the required geometry change. The low-frequency capacitance matrix is not its objective.

The essential strategy is **simulate a baseline, perturb dimensions, measure response sensitivities, cheaply search predicted candidates, then verify the selected geometry in Sonnet**. We could implement it using the ordinary 2-by-2 S matrix instead:

```text
S_candidate(f) approximately equals
    S_baseline(f)
  + pad_width_change * two_port_pad_sensitivity(f)
  + finger_span_change * two_port_finger_sensitivity(f)
  + ...
```

Here each sensitivity comes from ordinary two-port simulations of the baseline and perturbed metal. The optimizer can score the predicted S21 directly. This requires neither internal ports, short closure, nor a capacitance extraction. It also supports predicting combinations of changes without a fresh EM solve for every combination.

The two approaches need not give identical finite-step predictions. Approximating internal Z and then reducing it is a different approximation from approximating external S directly. Both require local steps and native validation. We did not perform a controlled comparison establishing which was better. A direct two-port sensitivity model would be a simpler first approach for an external-S21 objective; we must not retrospectively claim that it was the implementation used.

Short closure itself is exact linear algebra for the supplied network. Additional ports, their calibration, and fitting the frequency dependence can introduce numerical discrepancies. The restored-versus-original maximum all-S difference of about 0.00011846 measures one such check; it is not proof of improved accuracy. Direct two-port modeling avoids the port-restoration discrepancy, though its own geometry/frequency approximations still need validation. The accepted designs were checked with uninstrumented two-port EM, so the final filter result does not rely solely on the internal model.

The final raised-cover search initially transferred sensitivities measured under the older upper-boundary setting while recalibrating its baseline at 12 mm air height. `refine_quarter_raised_cover.py` records that transfer, checks native prediction error, updates the matrix sensitivities from observed steps and can request fresh derivatives. Its trust steps shrink after inaccurate or unsuccessful proposals; only certain validated pad-only moves permit expansion. Reuse saves measurements but does not demonstrate that the old sensitivities remain exact under the changed boundary.

### What happened in an actual final iteration

In the recorded cover-converged sequence, the pad-width increment control at +32 micrometres gave a native target-comparison RMS of **0.00250163429**, just above the acceptance threshold. Increasing that existing control by one more 2 micrometre step, to +34 micrometres, produced a native RMS of **0.00249186334** and passed. These numbers are changes relative to the geometry control's baseline, not absolute pad widths. The saved geometry vectors differ only at the pad-width control between these two candidates; the ladder comparator is fitted again within its bounds, so its numerical values need not remain fixed. The +34 record identifies this as one extra 2 micrometre step beyond the former search bound; it should not be portrayed as an independently proven internal-capacitance diagnosis.

There were also failed predictions: a regular-pad proposal predicted about 0.003164 RMS but gave about 0.008012 in the native simulation. It was rejected. This is why the measured geometry sensitivities guide the search, while the actual Sonnet response decides whether a candidate succeeds.

## 8. Why the upper air region matters

The physical layout is planar, but electric and magnetic fields extend into the substrate and above the metal. The simulation therefore needs a vertical environment and a boundary condition at its top. The native model's top cover is labeled `FREE_SPACE`; it is not an added metal ground plate. Its finite distance can nevertheless influence the calculated response.

The original upper air height was 1.5 mm. With unchanged metal, moving it to 12 mm changed complex S21 by RMS **0.00966824**, and reduced maximum port-1 power deficit `1−|S11|²−|S21|²` from **0.177283** to **0.00104456**. That was large relative to the desired design error. Refinement therefore used the farther boundary.

Increasing from 12 to 30 mm on the accepted taper changed S21 by RMS **3.14856e−9**. This supports top-height convergence for these settings. It does not identify exactly how much of the earlier deficit was physical radiation, dielectric loss, or interaction with the finite numerical boundary. It does not establish convergence with respect to every mesh or enclosure parameter.

The retained stack has 12 mm upper air, 500 µm silicon, 5 mm lower air, and a bottom perfect-conductor cover. A fabricated package or nearby metal lid could change the real device. A future package model should use its actual surroundings instead of assuming this boundary test proves every enclosure is equivalent. The native surface-model `Ls=0.0388` is preserved; its field-specific unit convention remains unconfirmed here.

## 9. Full-filter verification (report Section 13)

The full filter replaces the small fixture with ten regular cells and mirrored three-stage tapers. Its conditional circuit prediction repeats the regular values and mirrors the taper components, omitting the small fixture's output-end resonator and phase delay. That prediction is not a full-filter EM calculation.

Before committing to a long sweep, we solved one full-filter frequency, **4.92 GHz**, where the circuit prediction could be checked directly:

| Quantity | Circuit prediction | Full-filter EM |
|---|---:|---:|
| Linear magnitude | 0.957403974 | 0.972139673 |
| 20 log10(magnitude) | −0.378095480 dB | −0.245426656 dB |

EM is higher by **0.014735699** in magnitude, or **0.132668825 dB**. The calculation took 867.274 seconds. Phase conventions were not aligned, so this comparison is deliberately magnitude-only. It justified proceeding, but says nothing by itself about all ripple extrema or the band edges.

The full sweep uses 4.5–5.5 GHz ABS, a 2 µm grid, medium-speed edge meshing, double precision and de-embedding. Both ports are 50 ohms with 400 µm reference-plane shifts. The run uses 36 threads, selected by its launcher; the single downloadable project does not itself guarantee that runtime thread policy on another machine. ABS chooses direct frequencies adaptively, so its output interval target is not a known count of expensive solves. The supplied native project requests Touchstone S/MA output as `project.snp`, with `UseInterpolatedData=FALSE`; do not mistake a partial export for a converged full-band result.

**Completed full-filter result (November 10, 2026):** Sonnet finished 123 direct frequencies with final ABS SNR 40.9246 dB against a 40 dB goal and zero reported errors or warnings. The wrapper then incorrectly rejected the 123-point direct export because it expected at least 2000 points. An isolated copy of the completed cache, with interpolated output enabled and `-AbsNoDiscrete`, exported Sonnet's 2001-point curve with zero new EM frequencies. Thus the wrapper error was not an EM failure.

The verified native interpolated result is local at `investigations/filter_refinement_2026-11-05/quarter_percent_full_sweep_2026-11-08/interpolated_export/project.snp`; figures are in the sibling `completed_plots/` folder. At the chosen amplitude cutoff |S21| >= 0.8, the main continuous achieved band is approximately 4.8260-5.1681 GHz (342.1 MHz), estimated from that native 0.5 MHz export grid. This is a full uninstrumented-filter EM result, not an inference from the internal-port model. The completed responses, native input, metrics and figures are now published in `final_designs/final_full_filter/`; the original working files remain at the local paths above.

## 10. Files and reproduction

Current source tree:

```text
geometry/                 geometry, import, stage helpers, frozen-layout replay
tools/                    publication checks and saved-response replay
docs/                     this guide, technical method, history and audits
final_designs/
    final_regular/        historical standalone regular selection
    final_T1/             historical standalone T1 selection
    final_T2/             historical standalone T2 selection
    final_T3/             current complete taper fixture and response
    final_full_filter/    completed full-band input, response and evidence
        verification_4p92GHz/  completed independent single-point check
final_old/                superseded publication copies and old indexes
reference_inputs/awr/     supplied target circuit and optimizer screenshot
```

The project root also has `pyproject.toml`, `uv.lock`, `reference_taper_geometry.png` and `reference_passband_plot.png`. Local experiment directories `investigations/`, `sonnet_files/` and `output_data/` are intentionally excluded from Git. Moving the curated publication folders did not move the active solver or its files.

From the repository root:

```powershell
uv sync --locked
uv run python tools/check_published_files.py
uv run python -m geometry.rebuild_final --verify-only
uv run python tools/verify_current_response.py
```

The first verifier checks artifact hashes. The second checks native/GDS metal equality, the grid and the ports for each selected input. The third replays the saved terminal coefficients and normal comparator against the accepted native response, without EM. It should reproduce 0.0024918633384848348 and 8.117461191373584e−6 for the two RMS definitions above.

For native frozen-layout regeneration:

```powershell
uv run python -m geometry.rebuild_final --output-dir build/final_designs
```

This needs Sonnet and its license for GDS import/read. It does not launch EM. It is exact frozen-layout replay, not a promise that every old exploratory optimizer is portable. The complete `.sonx` is already sufficient to inspect and simulate the frozen design; GDS, Python scripts and old solver caches are not needed merely to open it. A new simulation needs a writable output location and a working Sonnet installation/license.

Technical source snapshots are in `docs/history/post_t2_evidence/`: `reachable_network.py.txt`, `reachable_iterate.py.txt`, `terminal_profiles.py.txt`, `quarter_target_optimizer.py.txt` and related files. Their original local sources are under `investigations/filter_refinement_2026-11-05/`. The portable checker is `tools/verify_current_response.py`, and model inputs are `final_T3/local_model/network.json`, `project.snp` and `result.json`. `final_designs/PUBLICATION_PATHS.json` records directory changes. The current source mapping is `final_designs/SOURCE_PATHS.json`; the preserved pre-move mapping is `final_old/publication_2026_11_08/SOURCE_PATHS.json`. Both exist and serve different historical roles.

The completed native full-filter response is now available. Unique original-component identification, finer mesh convergence and the behavior of an actual fabricated package remain separate questions. Internal ports were used in the local prediction model, but a benefit over ordinary two-port sensitivity optimization has not been demonstrated.


### Implementation reading map

| Calculation | Source/function | Saved evidence |
| --- | --- | --- |
| Port construction, polarity, closure and low/band fits | `reachable_network.py`: prepare, analyze, closez | final_T3/local_model/manifest.json, analysis.json, network.json |
| Early +2 micrometre geometry differences | `reachable_geometry.py`: main | original reachable_network controls and prospective validation |
| Bounded comparator and joint proposals | `quarter_target_optimizer.py`: distance, propose | final_T3/result.json |
| Derivatives of closure and ladder response | `quarter_analytic_response.py`: network_with_jac, target_with_jac | original analytic_derivative_validation.json |
| Actual final proposal and grid search | `refine_quarter_coils.py`: propose_network | candidate history and geometry vectors |
| Raised-boundary verification and model updates | `refine_quarter_raised_cover.py`: main, finish | final_T3/iteration_history.json, parameters.json, cover_30mm_verification/ |

Original scripts are under `investigations/filter_refinement_2026-11-05/`. Selected historical snapshots are under `docs/history/post_t2_evidence/`; they represent different stages and should not all be assigned the final bounds. `tools/verify_current_response.py` is the portable saved-response replay. Earlier T2 fixed-L capacitance fitting, later multiport modeling, bounded comparator fitting and geometry optimization are distinct operations.

### Short glossary

- **Native EM:** Sonnet solving the fields for a specified physical geometry.
- **Terminal network:** voltage/current or wave relationships between the defined ports.
- **Surrogate:** a cheaper numerical response predictor built from simulation data.
- **Training frequency:** a frequency used to determine model coefficients.
- **Held-out frequency:** a separately checked frequency not used for those coefficients.
- **Comparator:** the ideal ladder whose components are constrained near the requested targets.
- **Sensitivity:** measured response change per unit geometry change.
- **Trust step:** the limited range of geometry changes over which a local prediction is used.
- **Native verification:** testing the proposed metal rather than accepting the predictor's score alone.
