# Waveform Optimization for the `|f,n-1> <-> |g,n>` Sideband Interaction in cQED

## Executive Summary
- **Best selective waveform in the closed-system truth model:** wide Gaussian (`sigma_fraction = 0.24`) for both storage and readout sidebands.
- **Best fast/unselective waveform in the closed-system truth model:** square for both storage and readout sidebands.
- **Fastest closed-system selective duration:** storage `300 ns`, readout `220 ns`.
- **Fastest closed-system unselective duration:** storage `12 ns`, readout `12 ns`.
- **Best gate-oriented family within the present analytic family set:** Gaussian for storage and Blackman for readout, but neither is a phase-clean coherent SWAP across `n = 1,2,3`.
- **Main practical caveat:** once transmon decoherence is included using the matched local reference (`T1 = 9.813 us`, `T2,Ramsey = 6.325 us`, `Tphi = 9.332 us`), no pulse family continues to satisfy the original selective threshold, and only the storage fast-square pulse remains threshold-valid.
- **Simultaneous two-tone extension:** for `|g,0,1> -> |g,1,0>`, a resonant continuous two-tone drive reaches `0.999998` closed-system target transfer in `29.5 ns` and still peaks near `0.953` with noise, while the detuned Raman-like variant suppresses the intermediate `|f>` population to `0.121` but slows to `348 ns` and falls to about `0.556` once readout decay is included.

## Scope and Parameter Provenance
The study used the editable local `cqed_sim` device example `cQED_simulation/examples/sequential_sideband_reset.py` for the Hamiltonian and the mode noise terms. The key device parameters actually used in the simulations were:

| Parameter | Value |
|---|---:|
| Readout frequency | `8.596222556 GHz` |
| Qubit frequency | `6.150358764 GHz` |
| Storage frequency | `5.240932800 GHz` |
| Readout linewidth `kappa_r / 2pi` | `4.156 MHz` |
| Transmon anharmonicity `alpha / 2pi` | `-255.669695 MHz` |
| Storage dispersive shift `chi_s / 2pi` | `-2.840421 MHz` |
| Readout dispersive shift `chi_r / 2pi` | `-3.000000 MHz` |
| Storage `T1` | `250 us` |
| Storage `T2,Ramsey` | `150 us` |

The sideband-reset example does **not** export transmon coherence values. To study transmon decoherence explicitly, the open-system extension used a matched local reference from `cQED_simulation/examples/workflows/simulate_fock_tomo_and_sqr_calibration.py`, which shares the same storage/readout frequencies and transmon anharmonicity scale:

| Transmon reference parameter | Value |
|---|---:|
| `T1` | `9.813 us` |
| `T2,Ramsey` | `6.325 us` |
| Derived `Tphi` | `9.332 us` |

This transmon tuple is treated as a **clearly labeled sensitivity anchor**, not as a claimed device-specific measurement from the sideband-reset example itself.

## Model and Definitions
The truth model is the full multilevel dispersive transmon-storage-readout Hamiltonian replayed through the native `cqed_sim` sideband-control operator. The two target processes were:

- storage sideband: `|f,0_r,n_s-1> <-> |g,0_r,n_s>`
- readout sideband: `|f,n_r-1,0_s> <-> |g,n_r,0_s>`

The waveform families compared were:

- square
- Gaussian
- cosine
- flat-top cosine
- flat-top Gaussian
- smooth compact-support bump
- Blackman

Operational thresholds:

- **Selective pulse**
  - `P_target >= 0.99`
  - `P_leak <= 0.02`
  - `P_neighbor^max <= 0.01`
- **Fast/unselective pulse**
  - `P_target >= 0.985`
  - `P_leak <= 0.03`
- **Gate-oriented diagnostic**
  - same transfer/leakage floor as the fast metric plus `P_neighbor^max <= 0.02`
  - then rank by projected `2 x 2` SWAP fidelity

## Closed-System Ranking
### Storage Sideband
The selective winner is the wide Gaussian family:

| `n_s` | Duration | Amplitude | Target transfer | Max neighboring transfer |
|---|---:|---:|---:|---:|
| 1 | `220 ns` | `1.193 MHz` | `0.99405` | `0.00434` |
| 2 | `300 ns` | `0.619 MHz` | `0.99405` | `0.00619` |
| 3 | `300 ns` | `0.505 MHz` | `0.99407` | `0.00053` |

The fast/unselective winner is square:

| `n_s` | Duration | Amplitude | Target transfer | Max neighboring transfer |
|---|---:|---:|---:|---:|
| 1 | `12 ns` | `21.875 MHz` | `0.99695` | `0.54192` |
| 2 | `12 ns` | `15.468 MHz` | `0.99695` | `0.82029` |
| 3 | `12 ns` | `12.630 MHz` | `0.99696` | `0.92718` |

### Readout Sideband
The closed-system selective winner is also the wide Gaussian family:

| `n_r` | Duration | Amplitude | Target transfer | Max neighboring transfer |
|---|---:|---:|---:|---:|
| 1 | `220 ns` | `1.193 MHz` | `0.99405` | `0.00812` |
| 2 | `220 ns` | `0.844 MHz` | `0.99407` | `0.00894` |
| 3 | `220 ns` | `0.689 MHz` | `0.99408` | `0.00736` |

The fast/unselective winner is again square at `12 ns` for all three manifolds, with target transfer `0.99695-0.99696` and neighboring-manifold response `0.54-0.93`.

## Why These Families Win
The ranking follows the usual sideband bandwidth tradeoff, but the numbers make the tradeoff explicit in this device:

- the same-mode line spacing is only `5.68 MHz` for the storage sideband and `6.00 MHz` for the readout sideband
- square pulses win when speed is the only priority because they maximize pulse area per unit time
- the price is very broad spectral content and large neighboring-manifold response
- wide Gaussian pulses suppress high-frequency sidelobes enough to reach the selective threshold, even though they are roughly `18-25x` slower than the fast square winner

## Open-System Results
### Mode Noise Only
If the open-system model includes only the noise terms exported directly by the sideband-reset example, the storage recommendations survive while the readout selective regime already fails:

| Scenario | Mode | Regime | Mean noisy target |
|---|---|---|---:|
| mode-only | storage | selective Gaussian | `0.99225` |
| mode-only | storage | fast square | `0.99542` |
| mode-only | readout | selective Gaussian | `0.33126` |
| mode-only | readout | fast square | `0.92441` |

The storage selective pulse remains strong. The readout selective pulse is destroyed by the readout linewidth on the `220 ns` timescale.

### Including the Matched Local Transmon Reference
Adding transmon decoherence changes the practical recommendations further:

| Scenario | Mode | Regime | Mean noisy target | Worst noisy target |
|---|---|---|---:|---:|
| transmon reference | storage | selective Gaussian | `0.96075` | `0.95655` |
| transmon reference | storage | fast square | `0.99391` | `0.99385` |
| transmon reference | readout | selective Gaussian | `0.32240` | `0.09952` |
| transmon reference | readout | fast square | `0.92300` | `0.87740` |

Two conclusions matter:

1. The long storage selective Gaussian is still the best **structural** selective family, but under the original strict threshold it no longer qualifies once the matched transmon decoherence is included.
2. The fast storage square pulse remains comfortably above the fast/unselective threshold.

The reranked transmon-inclusive winner table is therefore extremely stark: **only the storage fast-square family remains threshold-valid**. No storage selective family survives, and no readout family survives.

## Simultaneous Two-Tone Transfer Extension
The study was extended to ask whether a continuous storage-sideband drive plus a simultaneous readout-sideband drive can directly move a single storage photon into the readout mode. For the single-photon manifold, the dominant reduced ladder is

- `|g,0,1> <-> |f,0,0> <-> |g,1,0>`

which supports two useful limits.

1. **Resonant bright-state transfer**
  - Equal `12 MHz` leg couplings and zero common detuning produce a near-perfect closed-system transfer peak of `0.999998` at `29.5 ns`.
  - The price is large intermediate-state participation: the peak `|f,0,0>` probability is about `0.500`.
  - With mode noise only, the target still peaks at `0.9554`; with the matched transmon reference included, it peaks at `0.9530`.

2. **Detuned Raman-like transfer**
  - Equal `4 MHz` leg couplings with `20 MHz` common detuning reduce the peak intermediate-state population to `0.121`.
  - The transfer then slows to `348 ns` in the closed system.
  - Under the noisy replay, the target peak drops to `0.5606` with mode noise only and `0.5560` with the matched transmon reference, so readout decay removes most of the practical benefit of the darker path.

The extension also shows that this clean reduced three-state picture is specific to `n = 1`. For `n = 2,3`, the same constant simultaneous drive opens a longer conversion chain, and the best resonant target probability stalls around `0.61-0.69` even before noise is added. In this parameter regime, the continuous two-tone protocol is therefore a plausible fast single-photon transfer primitive, but not a clean general transfer solution for arbitrary storage Fock states.

## Gate-Oriented Diagnostic
The added gate-oriented ranking confirms that the present study family set should be interpreted as state-transfer control, not as coherent SWAP synthesis.

| Mode | Best family under gate-oriented diagnostic | Conservative duration | Mean projected SWAP fidelity | Mean absolute phase asymmetry |
|---|---|---:|---:|---:|
| Storage | Gaussian | `300 ns` | `0.716` | `0.944 rad` |
| Readout | Blackman | `300 ns` | `0.765` | `0.839 rad` |

These are not phase-clean coherent gates. Even the best gate-oriented family within the analytic ansatz set remains far from a uniformly high-fidelity `2 x 2` SWAP across `n = 1,2,3`.

## Validation
Three checks were repeated after the extended study:

1. **Analytic sanity**
   - the exact `cqed_sim` rotating-frame sideband frequencies agree with the dispersive formulas to numerical precision for both modes and all `n = 1,2,3`
2. **Timestep convergence**
   - finalist baselines used `0.25 ns`
   - using `0.5 ns` changes target transfer by only `1.8e-4` to `1.4e-3`
   - using `1.0 ns` changes target transfer by `1.5e-3` to `3.8e-3`
3. **Truncation convergence**
   - increasing from `(n_tr, n_s, n_r) = (4,5,5)` to `(5,6,6)` changes finalist target transfer by at most `5.7e-5`
4. **Two-tone reduced-model validation**
  - for the selected single-photon resonant case, the full simulator and reduced three-state model agree to within `4.1e-7` in peak target probability and exactly in peak time
  - the selected single-photon cases change by at most `2.5e-3` in peak target probability when rerun at `1.0 ns` instead of `0.25 ns`
  - the selected single-photon cases change by less than `2.3e-6` in peak target probability when the truncation is enlarged to `(5,6,6)`

## Practical Recommendation
The recommendation now depends on which layer of realism matters.

- **Closed-system control design:** use a wide Gaussian for selective sideband transfer and a square pulse for the fastest strong transfer.
- **Mode-noise-only open-system replay:** storage selective Gaussian remains viable; readout selective Gaussian does not.
- **Transmon-inclusive practical replay using the matched local reference:** only the storage fast-square pulse remains threshold-valid under the original metrics.

That means the current study supports a careful experimental recommendation:

- use storage-sideband Gaussian pulses if you are designing a spectrally selective control primitive and the transmon coherence is much better than the matched local reference
- use storage-sideband square pulses if you want the most robust fast transfer in the current model
- do not rely on long readout-sideband selective pulses in this hardware regime
- a resonant simultaneous two-tone pulse is viable if the goal is a fast single-photon storage-to-readout transfer and substantial temporary `|f>` occupation is acceptable
- a detuned constant two-tone pulse is not competitive in this device because readout decay dominates before the slower Raman-like transfer can complete
- do not interpret any of the present analytic families as a finished coherent SWAP gate without additional phase compensation or unitary-level optimization

## Limitations and Future Work
- The sideband remains an **effective control operator**, not a microscopic pump Hamiltonian.
- The sideband-reset example still does not export transmon coherence values, so the transmon-noise extension is a matched local sensitivity study rather than an exact device-verification pass.
- The gate-oriented extension added a ranking diagnostic, but not a full unitary-level optimal-control search.
- The simultaneous two-tone extension only tests constant square pulses, so it does not yet address shaped counter-intuitive overlap, STIRAP-like timing, or photon-number-selective transfer design.
- A next study should combine:
  - a microscopic pump-aware model,
  - direct unitary optimization,
  - shaped or optimized simultaneous two-tone control for the single-photon transfer problem,
  - and device-verified transmon `T1/T2` values on the exact sideband hardware instance.
