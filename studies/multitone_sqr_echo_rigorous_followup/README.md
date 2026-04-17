# Rigorous Echo-Ansatz Follow-Up for Multitone SQR in Dispersive cQED

## Problem Class
ANA | REP | DES

## Motivation
This follow-up reruns the strict simultaneous shared-line multitone-SQR study with a narrower and more rigorous question about the echoed construction. The previous study established a controlled no-go result for the no-detuning amplitude-plus-azimuth ansatz and found that two simple echoed replays did not rescue the gate. The remaining concern is whether that echo verdict was under-optimized.

This study therefore keeps the same strict physical model and no-detuning control restriction, but upgrades the echo section from inherited half-pulse replay to explicit sequence-level echo optimization. The goal is to decide, on stronger evidence, whether the echoed ansatz can actually rescue the gate or only reshuffle the error budget.

All Hilbert-space objects follow the `cqed_sim` qubit-first convention: qubit tensor cavity, with logical basis ordered as `(|g,0>, |e,0>, |g,1>, |e,1>, ...)`.

## Study Composition
| Component | Source | Prior role | Role here |
|---|---|---|---|
| Strict no-detuning shared-line baseline | `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed` | Established the controlled no-go statement and baseline numerics | Reproduced here as the physical baseline and report comparison point |
| Sequence-level echo optimization pattern | `studies/ideal_sqr_direct_vs_echoed_multitone` | Broader echoed optimizer with extra controls | Reused only as implementation prior art for composite-sequence optimization structure |
| Targeted-subspace diagnostics | `cqed_sim.calibration.targeted_subspace_multitone` | Framework validation metrics | Primary metric engine in this follow-up |

## Goals
1. Reproduce the strict no-detuning simultaneous-shared-line baseline on the same physical model.
2. Replace the old echoed replay with genuinely optimized echoed ansatz families.
3. Make every reported metric explicit, including which ones are phase-sensitive and which ones can hide failure modes.
4. Test whether the echoed sequence fails because of under-optimization, finite refocusing pulses, or a deeper structural limitation.
5. Apply the most important prior future-work items: fair duration-matched comparison and a manifold-aware refocusing alternative.
6. Regenerate the full report, figures, machine-readable artifacts, and reproducibility notebook.

## Methods
- `cqed_sim.core.DispersiveTransmonCavityModel` and `FrameSpec` for the shared-line dispersive model.
- `cqed_sim.calibration.conditioned_multitone` for tone generation and waveform construction.
- `cqed_sim.calibration.targeted_subspace_multitone.evaluate_targeted_subspace_multitone` and `analyze_targeted_subspace_operator` for strict restricted-subspace diagnostics.
- `cqed_sim.sequence.SequenceCompiler` and `cqed_sim.sim.prepare_simulation` for exact full-sequence replay when finite refocusing pulses are included.
- Local helpers only where no public `cqed_sim` API exists:
  - exact reduced blockwise replay,
  - decoupled-block comparator,
  - logical probe-state metrics,
  - echoed sequence wrappers and joint echo optimizers.

## Analytic Preliminary
The first-principles baseline remains the dispersive block-resolved Hamiltonian
```text
H_0 = sum_n (omega_n / 2) sigma_z tensor |n><n|
```
driven through one shared qubit line by simultaneous resonant tones placed exactly at the addressed block frequencies. In the interaction frame of `H_0`, each block sees its resonant transverse term plus all spectator tones from the other addressed manifolds. The old no-go argument already shows that the spectator tones generate effective blockwise `Z` generators that are not generically removable with amplitude and azimuth knobs alone.

For the echoed construction, the ideal toggling-frame algebra says that
```text
U_echo = X_pi U_2 X_pi U_1
```
can cancel a blockwise `Z` contribution only if the inserted `X_pi` really acts as the same `pi` rotation in every addressed block and the two halves retain the required transverse action. This follow-up explicitly tests how much survives when that echoed ansatz is actually optimized rather than inherited from a half-target replay.

Controlled approximations:
- The baseline analytical no-go remains a dispersive block-resolved argument rather than an all-regime theorem.
- Ideal instantaneous refocusing is treated as an upper-bound thought experiment.
- Finite refocusing is tested with exact full-sequence propagation rather than extra rotating-wave simplifications beyond those already built into the dispersive model.

## cqed_sim Gap Analysis
| Functionality | Needed? | Available in cqed_sim? | Plan |
|---|---|---|---|
| Strict shared-line multitone propagation | Yes | Yes | Use `evaluate_targeted_subspace_multitone` and exact replay |
| Sequence-level echoed ansatz optimization with multiple segments | Yes | Partial | Compose local wrappers around `cqed_sim` evaluations |
| Restricted-subspace fidelity, leakage, and block-phase diagnostics | Yes | Yes | Use `analyze_targeted_subspace_operator` |
| Fair direct-versus-echo total-duration comparison | Yes | Yes | Re-run direct optimization at the matched total duration |
| Manifold-aware refocusing pulse benchmark | Yes | Partial | Optimize a shared-line multitone `X_pi` pulse locally using the same strict controls |
| Probe-state fidelity diagnostics across cavity superpositions | Yes | No public helper | Implement locally and cross-check against framework metrics |

## Assumptions
- Dispersive block-resolved qubit-cavity model with shared qubit control line.
- No artificial per-tone detuning in any optimized construction.
- Amplitude and azimuth corrections only for all multitone segments, including the manifold-aware refocusing variant.
- Primary active window sizes: `N_active = 2, 3`.
- Primary duration grid: `|chi| T / 2pi = 3, 5`.
- Model variants: `chi_only` and `chi_plus_chiprime`.
- Target families: aligned `x` as the best-case echo setting and structured `XY` as the generic non-aligned setting.

## Compute & Resource Strategy
- Expected bottleneck: joint echoed optimization with exact full-sequence replay for finite refocusing pulses.
- Planned strategy:
  - use a smaller but more carefully chosen case grid than the previous broad sweep,
  - run a pilot case first to calibrate optimization budgets,
  - keep ideal-instantaneous echo optimization operator-composed where possible,
  - use exact compiled full-sequence replay only where finite refocusing pulses make it necessary.
- Realized strategy:
  - production budget `direct_starts=2`, `direct_maxiter=10`, `segment_maxiter=8`, `echo_ideal_maxiter=8`, `echo_gaussian_maxiter=5`, `refocus_maxiter=8`,
  - exact reduced replay and decoupled-block checks remained inexpensive,
  - the full production run completed in about `1810 s` on CPU without extra packages, GPU backends, or multiprocessing.

## Expected Outcomes
- A stronger answer to whether the echo sequence truly fails under the strict no-detuning shared-line model.
- Clear separation between:
  - ideal instantaneous echoed upper bound,
  - finite vacuum-calibrated echoed replay,
  - optimized finite echoed sequence,
  - manifold-aware refocusing alternative.
- Explicit metric definitions showing which success claims are or are not supported.

## Known Limitations
- This follow-up still studies the strict no-detuning shared-line family, not richer segmented controls with explicit detuning or blockwise `Z` compensation.
- Sequence-level optimization remains numerical; a closed-form no-go for every echoed variant may still be unavailable.
- The manifold-aware multitone-refocusing benchmark was limited to the representative harder subset with `chi_plus_chiprime` and `|chi|T/2pi = 5`, because even the refocusing pulse itself was already clearly poor under the strict ansatz.

## Findings
1. The strict shared-line direct pulse remained the strongest physical baseline, with mean restricted average gate fidelity `0.7133`.
2. The exact reduced blockwise replay matched the full strict result to machine precision, confirming that the failure is already present in the block-resolved shared-line dynamics.
3. The replayed ideal instantaneous echo strongly reduced residual-`Z` while still failing badly as a gate: mean max residual-`Z` `0.0098 rad`, mean restricted average gate fidelity `0.2006`, mean explicit probe fidelity `0.0886`.
4. Jointly optimizing the ideal instantaneous echo improved it substantially and beat the plain direct pulse in `8/16` cases, but only at `|chi|T/2pi = 5`; the mean remained below the direct baseline and no exact gate was found.
5. The finite Gaussian echoed sequence still failed as a practical rescue: mean restricted average gate fidelity `0.3176`, no cases better than the active-duration direct pulse, and only `3/16` cases better than the total-duration-matched direct comparator.
6. The manifold-aware shared-line multitone refocusing benchmark did not change the physical verdict: the refocusing pulse itself had mean fidelity `0.1805`, and the echoed construction built from it stayed below the plain direct baseline on the tested hard subset.

## Validation
- [x] Sanity checks
  - Exact reduced replay versus the full shared-line pulse matched with minimum reduced-versus-full restricted process fidelity `1.0`.
  - The stronger decoupled-block comparator reproduced the ideal target with fidelity `1.0` in every tested case.
- [x] Convergence
  - Representative higher-budget direct and echoed reruns remained consistent with the production conclusion; no finite echoed construction approached unit fidelity under the stricter budgets.
- [x] Literature comparison (if applicable)
  - The underlying dispersive and echo-frame reasoning remained consistent with standard cQED and spin-echo expectations.
  - The local repository audit still shows that earlier nearby positive-looking cases used extra resources such as `d_omega` or broader waveform families, so they are not direct counterexamples to the strict claim.

## Status
COMPLETE
