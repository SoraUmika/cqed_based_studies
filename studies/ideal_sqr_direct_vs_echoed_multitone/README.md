# Ideal SQR Feasibility with Direct and Echoed Multitone Waveforms

## Problem Class
OPT | ANA | DES

## Motivation
The two earlier studies in this repository answered adjacent but not identical questions.

1. `multitone_sqr_arbitrary_fock_conditional_rotations` optimized a single Gaussian multitone waveform for arbitrary block-diagonal SU(2) targets on an active Fock window. That is a harder and more general control problem than an ideal SQR gate, but it is not the same target.
2. `parameterized_waveform_residual_z_cancellation` compared richer waveform families, but its `echoed_multitone` family used only one mid-sequence `X_pi` pulse rather than the requested
   \[
   \text{half-SQR} \rightarrow \pi \rightarrow \text{half-SQR} \rightarrow \pi,
   \]
   and it optimized a scalarized fidelity-plus-residual-Z objective instead of a pure ideal-SQR gate target.

This follow-up study therefore re-centers the question on the narrower and physically cleaner target:

\[
U_{\mathrm{SQR}}^{\mathrm{ideal}}(\{\theta_n\}) =
\sum_{n=0}^{N_{\mathrm{active}}-1} |n\rangle\!\langle n| \otimes R_x(\theta_n),
\qquad
R_x(\theta)=e^{-i\theta \sigma_x/2}.
\]

The phase is fixed to the qubit \(x\) axis so that the echoed construction can be tested exactly in the form requested by the user. All simulation objects follow the `cqed_sim` qubit-first convention: Hilbert-space operators are interpreted as qubit \(\otimes\) cavity, and the logical basis is ordered as \((|g,0\rangle, |e,0\rangle, |g,1\rangle, |e,1\rangle, \ldots)\).

## Goals
1. Audit the two earlier studies and record the specific ways in which they do and do not answer the ideal-SQR question.
2. Define a direct multitone Gaussian-envelope waveform ansatz for an ideal SQR target with per-manifold angles \(\theta_n\) and common \(x\)-axis phase.
3. Define and test the echoed multitone construction
   \[
   \text{half-SQR} \rightarrow X_\pi \rightarrow \text{half-SQR} \rightarrow X_\pi
   \]
   in three scenarios:
   - fully symmetric echoed,
   - equal-duration independently corrected,
   - weakly asymmetric echoed.
4. Derive low-\(n\) analytical conditions showing when direct or echoed parameterizations can cancel first-order residual \(Z\)-type errors in principle.
5. Quantify whether any improvement from the echoed construction is symmetry-protected or only obtained through extra asymmetry or fine tuning.
6. Save machine-readable artifacts, figures, validation results, a reproducibility notebook, and a unified report.

## Methods
- `cqed_sim.core.DispersiveTransmonCavityModel` and `FrameSpec` for the dispersive transmon-cavity Hamiltonian with `chi` and optional `chi'`.
- `cqed_sim.core.frequencies.manifold_transition_frequency` and `carrier_for_transition_frequency` for conditioned qubit-tone and echo-pulse carriers.
- `cqed_sim.calibration.conditioned_multitone.build_conditioned_multitone_tones` and `build_conditioned_multitone_waveform` for the direct multitone and half-SQR waveform ansatz.
- `cqed_sim.calibration.targeted_subspace_multitone.optimize_targeted_subspace_multitone`, `build_block_rotation_target_operator`, `build_spanning_state_transfer_set`, and `analyze_targeted_subspace_operator` for direct ideal-SQR optimization and logical-subspace diagnostics.
- `cqed_sim.sequence.SequenceCompiler` and `cqed_sim.sim.prepare_simulation` for replay of the full echoed pulse schedule.
- Local study code only for the composite echoed optimization layer, analytic low-\(n\) reduction, figure generation, report writing, and notebook generation.

Analytical starting point:
- The direct multitone ansatz is treated in each active manifold as an effective driven qubit with desired \(x\)-rotation area \(\Theta_n\) and unwanted residual \(Z\)-phase \(\Phi_n\).
- For the ideal echoed schedule with \(x\)-axis half-rotations, \(X_\pi\) commutes with the desired \(x\)-rotation while toggling \(Z \mapsto -Z\), so identical halves cancel first-order \(Z\)-type phase if the two halves experience the same residual \(Z\) accumulation.
- The analysis explicitly records the obstruction that the echo loses this clean property once the finite \(X_\pi\) pulses become manifold dependent.

Numerical study grid:
- model variants: `chi_only`, `chi_plus_chiprime`
- target families: `smooth_x`, `staggered_x`
- active windows: `N_active = 2, 3`
- durations: \(|\chi|T/2\pi = 3, 5`

## Expected Outcomes
- The direct waveform should be competitive when the ideal SQR target is smooth and the active window is small, but it may leave manifold-dependent residual \(Z\)-type error because the simultaneous spectator tones induce coherent Stark-like phase shifts.
- The fully symmetric echoed construction should cancel first-order \(Z\)-type error in the reduced analytical model, but its practical success depends on whether the inserted finite \(X_\pi\) pulses are sufficiently manifold independent.
- If the independently corrected or weakly asymmetric echoed cases are the only ones that improve materially, that will indicate the apparent gain is not a clean symmetry argument but a more fragile fine-tuning effect.

## Known Limitations
- `cqed_sim` does not currently expose a public optimizer for multi-segment targeted-subspace multitone sequences, so the direct single-segment optimization stays inside the public API while the echoed composite optimization is implemented locally on top of the same simulation stack.
- The echoed sequence uses finite Gaussian \(X_\pi\) pulses, not instantaneous ideal kicks. Their manifold dependence is part of the tested physics, not an external correction.
- The current follow-up targets ideal \(x\)-axis SQR gates rather than arbitrary XY-axis conditional rotations. That matches the user’s \(\theta_n\)-only specification and the cleanest echo algebra, but it does not exhaust every possible ideal-SQR convention.
- The main grid is intentionally focused on \(N_{\mathrm{active}}=2,3\) to keep the direct-vs-echoed comparison analytically interpretable and computationally tractable.

## Suggested Upstreaming
- Add a public `cqed_sim` helper for multi-segment targeted-subspace conditioned multitone optimization with shared or tied correction parameters.
- Add a public echoed-SQR builder that composes `ConditionedMultitoneWaveform` segments with calibrated refocusing pulses and exposes a toggling-frame diagnostic.
- Add logical-subspace diagnostics that separate residual \(Z\)-type error, transverse coherent error, and manifold-dependent echo-pulse error directly in the public API.

## Status
COMPLETE
