# Fast and Robust Active Cooling / Vacuum Reset in a Transmon-Storage-Readout cQED System

## Problem Class
DES, ANA, OPT

## Motivation
This study compares active-cooling architectures for a storage-transmon-readout cQED device, with the goal of identifying which protocol is fastest, which is most robust to transmon decoherence, and which offers the best overall speed-versus-robustness tradeoff on the local device parameters exposed through `cqed_sim`.

## Goals
1. Compare pulsed ladder cooling, continuous resonant bright-state cooling, continuous detuned Raman-like cooling, and an auxiliary effective autonomous-cooling benchmark.
2. Quantify cooling time, residual storage occupation, residual transmon excitation, leakage, and robustness to `T1`, `T2`, dephasing, detuning error, amplitude error, reset error, and thermal loading.
3. Determine how much benefit the readout resonator provides as the dominant dissipative dump.
4. End with a decisive recommendation for the best current experimental target and the best longer-term autonomous-cooling direction.

## Methods
- Native `cqed_sim` three-mode multilevel Lindblad replay for all physical schemes.
- Reuse of the validated storage and readout `g-f` sideband pulse recommendations from the earlier local studies.
- Targeted coupling, detuning, ringdown, decoherence, calibration, and readout-linewidth sweeps.
- A small reduced-model benchmark only for the idealized effective `L \propto a_s` limit.

## Analytic Preliminary
The numerical results support the basic first-principles expectation: schemes that rely on large real transmon occupation can be very fast, but they become fragile under decoherence and calibration error; schemes that keep the transmon more virtual are more robust, but slower unless the dissipative channel is made stronger.

## cqed_sim Gap Analysis
| Functionality | Needed? | Available in cqed_sim? | Plan |
|---|---|---|---|
| Full three-mode driven Lindblad replay | Yes | Yes | Use native model directly |
| Pulsed sequential ladder with readout dump | Yes | Yes | Use native sideband channels |
| Continuous simultaneous storage/readout sidebands | Yes | Yes | Use overlapping native channels |
| Effective storage-only jump operator `L \propto a_s` | Helpful | No | Use a clearly labeled reduced benchmark only |

## Assumptions
- The local sideband-reset example is the authoritative device tuple.
- The matched local tomography workflow provides the transmon coherence sensitivity anchor.
- The continuous protocols are centered on the `n=1` sideband line, so higher-photon performance is a real diagnostic of scalability rather than an optimized many-photon control result.

## Compute & Resource Strategy
The study reuses earlier validated pulse choices and spends the compute budget on cross-scheme robustness sweeps instead of re-optimizing waveform families. The latest full artifact refresh completed in `2510.6 s` on CPU, followed by a targeted initial-state artifact refresh using the saved selected settings.

## Expected Outcomes
The final outputs now include ranked scheme recommendations, saved figures and machine-readable artifacts, a technical report, and a reproducibility notebook.

## Known Limitations
- The sideband control layer is still effective rather than pump-microscopic.
- The autonomous benchmark is auxiliary and not a direct device replay.
- The continuous schemes are not yet optimized with shaped counter-intuitive timing or optimal control.
- The continuous-scheme end-of-window residual occupations remain time-step sensitive between `0.5 ns` and `1.0 ns`, so headline final-state numbers are quoted from the finer baseline while robustness sweeps still use the coarser grid.

## Validation
- [x] Sanity checks
- [ ] Convergence
- [ ] Literature comparison (if applicable)

## Status
ACTIVE
