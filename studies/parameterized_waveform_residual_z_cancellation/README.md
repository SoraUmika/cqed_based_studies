# Parameterized Waveform Design for Fock-Conditional Qubit Rotations with Residual-Z Cancellation

## Problem Class
OPT | ANA | DES

## Motivation
The completed baseline study `multitone_sqr_arbitrary_fock_conditional_rotations` established a clear negative result for the standard Gaussian multitone SQR ansatz: arbitrary block-diagonal qubit rotations on a finite active Fock window are not realized at high fidelity, and the dominant coherent error is not explained by spectator-population crosstalk alone. The next question is narrower and more constructive: can richer but still structured waveform families suppress the residual blockwise qubit-Z error while preserving the intended conditional qubit rotations?

This study keeps the same ideal dispersive transmon+cavity model with `chi` and optional `chi'` only, no cavity self-Kerr in the main comparison, and the same qubit-first `cqed_sim` convention. Hilbert-space operators are interpreted as qubit `otimes` cavity, and the logical basis is ordered as `(|g,0>, |e,0>, |g,1>, |e,1>, ...)`.

## Goals
1. Reuse the strict active-subspace diagnostics from the baseline arbitrary-block study as the common comparison layer for richer waveform families.
2. Implement and compare five parameterized waveform families: baseline Gaussian multitone, symmetric two-segment multitone, explicit echoed multitone, complex-envelope multitone, and basis-expanded multitone.
3. Quantify whether the richer families reduce the coherent error unitary's blockwise Z component relative to the baseline ansatz.
4. Measure the tradeoff between residual-Z suppression, overall active-subspace fidelity, and leakage/cavity-block preservation.
5. Compare `chi`-only and `chi + chi'` models on a larger grid spanning `N_active = 2, 3, 4`, `chi T / 2pi = 3, 5`, one structured target family, and a small random ensemble.
6. Save machine-readable artifacts, figures, and a study summary that support follow-on analysis rather than a pilot-only snapshot.

## Methods
- `cqed_sim.core.DispersiveTransmonCavityModel` and `FrameSpec` for the dispersive qubit-cavity Hamiltonian.
- `cqed_sim.core.frequencies.manifold_transition_frequency` for conditioned qubit transition frequencies.
- `cqed_sim.calibration.conditioned_multitone` for tone construction and sampled multitone pulse building.
- `cqed_sim.calibration.targeted_subspace_multitone` for the baseline targeted-subspace multitone optimizer and for the logical-subspace diagnostics reused across all families.
- `cqed_sim.sequence.SequenceCompiler` and `cqed_sim.sim.prepare_simulation` for direct replay of arbitrary sampled pulse sequences.
- Local study code for family-specific envelope parameterizations, local optimization over sampled envelopes, residual-Z error decomposition, and figure/report generation.

## Expected Outcomes
- If the baseline failure is mostly a blockwise Z-coherence problem, symmetric or complex-envelope families should reduce the mean residual-Z error angle while improving average gate fidelity on the same targets.
- If the failure is instead a deeper controllability or bandwidth limitation, richer families may only weakly improve the fidelity and leave large non-Z coherent error components.
- `chi'` is expected to make residual-Z cancellation harder by breaking simple symmetry relations across active manifolds.
- A decisive outcome would separate two regimes: smooth structured targets where single-segment richer families may be enough, and random block targets where echoed structure or still-richer protocols may be necessary.

## Current Findings
- The expanded comparison grid is complete: `48` cases and `240` waveform-family rows covering `N_active = 2, 3, 4`, `chi T / 2pi = 3, 5`, structured target family `C`, and three random `D` seeds for each model/duration/active-window combination.
- On the smooth structured targets, the single-segment richer families still help only marginally. The best overall point is the basis-expanded family at `0.874158` average gate fidelity for `chi + chi'`, `N_active = 2`, `chi T / 2pi = 5`, compared with `0.873359` for the baseline.
- On structured targets, the explicit echoed family is decisively the wrong protocol in its current form: its mean fidelity drops to `0.164414` versus `0.763677` for the baseline, while both residual-Z and transverse coherent error increase.
- On random targets, the picture reverses. The echoed family wins `23` of the `36` random-target cases, raising mean fidelity from `0.247163` to `0.321155`, reducing mean residual-Z error from `0.980429` rad to `0.627334` rad, and reducing mean transverse error from `1.790404` rad to `1.590729` rad.
- The complex-envelope and basis-expanded families remain the best single-segment refinements: they improve mean fidelity only slightly over the baseline, but preserve the good structured-target performance while modestly lowering transverse coherent error.
- Fidelity still degrades strongly with larger active windows. For the structured family, the best fidelity falls from about `0.874` at `N_active = 2` to about `0.780` at `N_active = 3` and about `0.674` at `N_active = 4`. For the harder random ensemble, median fidelities remain well below `0.5` across the full grid.

## Known Limitations
- The richer families are implemented locally because `cqed_sim` does not currently expose echoed or basis-expanded targeted-subspace multitone optimizers.
- The present echoed family uses a single global mid-sequence Gaussian `X_pi` pulse calibrated on the `n = 0` manifold. It is not a manifold-selective echo and is not jointly optimized with a second corrective qubit operation.
- The pilot uses a two-level qubit model (`n_tr = 2`) to stay aligned with the targeted-subspace logical validation workflow.
- The main comparison remains closed-system and does not include open-system noise in the optimization objective.
- No formal convergence sweep over timestep, cavity truncation, or optimizer budget has been completed yet; the present conclusions rely on internal consistency checks and cross-family comparison at fixed numerics.
- The random ensemble remains modest at three seeds per configuration, so the hard-target statistics are informative but not yet production-level.

## Suggested Upstreaming
- Add a public `cqed_sim` API for targeted-subspace optimization over arbitrary sampled multitone envelopes rather than only the built-in Gaussian family.
- Add built-in logical diagnostics that decompose the blockwise error unitary into Z-like and transverse components.
- Add support for multi-segment conditioned multitone sequences with shared tone specifications and mirrored or echoed envelope constraints.

## Status
COMPLETE