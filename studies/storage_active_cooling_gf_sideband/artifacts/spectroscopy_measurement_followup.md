# Experiment-Facing Spectroscopy Follow-Up

## Main answers
- To locate the Step A storage-cooling transition, the primary observable should be the population in `|f,0_r,n_s-1>` after a weak long storage-sideband probe. In practice this means read out `P_f`, not an `X/Y` projection, during the first spectroscopy pass.
- A direct transmon `g-f` carrier probe is **not** the cooling transition because it keeps the storage photon number fixed. It is still useful as an auxiliary calibration for qutrit readout and for preparing `|f,0_r,n_s-1>` when calibrating Step B.
- If direct `f` discrimination is weak, the best fallback is a dump-assisted witness: scan Step A, then apply the calibrated Step B dump pulse and read out the resulting readout occupation or integrated ringdown.
- `X/Y` projections in the `g-f` manifold are needed only after the resonance is found, for Ramsey/tomography-style checks of phase, Stark shift, and axis control.
- Neighboring `n_s`-resolved Step A and Step B lines are only about `5.680842 MHz` apart, so the first search should use weak, long probes rather than the short high-fidelity control pulses.

## Exact transition search workflow
1. Prepare `|g,0_r,n_s>` for the desired storage manifold.
2. Coarse scan the Step A sideband around the predicted line with a weak square probe of duration `400 ns` and amplitude `0.30/sqrt(n_s) MHz`.
3. Measure either direct `P_f` or the dump-assisted readout witness.
4. Fit the peak and then repeat a fine scan over `+-2 MHz` with smaller steps.
5. At the fitted frequency, perform a duration sweep to calibrate the Rabi period and then switch to the faster high-fidelity pulse family from the main study.
6. For Step B, prepare `|f,0_r,n_s-1>` using either the calibrated Step A pulse or the direct transmon `g-f` carrier, then scan the readout sideband and monitor readout occupation.
7. Only after those `Z`-axis calibrations are stable should you run `g-f` Ramsey or analysis-pulse tomography to project onto `X_{gf}` and `Y_{gf}`.

## Predicted exact lines and recommended observables
| n_s | Direct g-f carrier (GHz) | Step A sideband (GHz) | Step B dump (GHz) | Preferred Step A readout | Preferred Step B readout |
|---|---:|---:|---:|---|---|
| 1 | 12.039366992 | 6.804115034 | 3.448825278 | `P_f` or dump witness | `\langle n_r \rangle` or integrated ringdown |
| 2 | 12.033686150 | 6.798434192 | 3.443144436 | `P_f` or dump witness | `\langle n_r \rangle` or integrated ringdown |
| 3 | 12.028005308 | 6.792753350 | 3.437463594 | `P_f` or dump witness | `\langle n_r \rangle` or integrated ringdown |
| 4 | 12.022324466 | 6.787072508 | 3.431782752 | `P_f` or dump witness | `\langle n_r \rangle` or integrated ringdown |

## Simulated weak-drive spectroscopy peaks
| n_s | Step A peak detuning (MHz) | Step A peak signal | Step B peak detuning (MHz) | Step B peak signal |
|---|---:|---:|---:|---:|
| 1 | +0.00 | 0.4674 | +0.00 | 0.4679 |
| 2 | +0.00 | 0.4673 | +0.00 | 0.4678 |
| 3 | +0.00 | 0.4673 | +0.00 | 0.4678 |
| 4 | +0.00 | 0.4672 | +0.00 | 0.4677 |

## What to measure in each experiment
- **Step A coarse spectroscopy:** measure `P_f`. This is a `Z_{gf}` measurement, not an `X/Y` projection.
- **Step A fallback when `f` readout is weak:** append the calibrated Step B pulse and detect the induced readout photon as a dump witness.
- **Step A Rabi calibration:** stay with `P_f` versus pulse duration.
- **Step A Stark-shift and axis check:** switch to `g-f` Ramsey with a final `pi/2` analysis pulse so that `X_{gf}` or `Y_{gf}` is mapped onto readout.
- **Step B spectroscopy:** measure readout occupation or integrated homodyne ringdown, because the target state already lives in `|g,1_r,n_s-1>`.
- **Full cooling validation:** compare storage number-splitting or Wigner/tomography before and after repeated cycles, not just transmon population.

## Most useful experiment sequence
The cleanest lab sequence is:
1. Calibrate qutrit transmon readout (`g/e/f`) and the direct `g-f` carrier first.
2. Use the storage-sideband probe to find Step A at each `n_s`.
3. Use the direct `g-f` carrier from `|g,0_r,n_s-1>` to prepare `|f,0_r,n_s-1>` and calibrate Step B independently.
4. Reconnect the two calibrated pieces into the full cooling primitive.

The central practical point is that you should not begin with full `g-f` tomography. Begin with population spectroscopy in the `g/f` manifold, then add `X/Y` projections only after the line center and pulse area are known.
