# Can an Ideal Multitone SQR Realize Arbitrary Fock-Conditional Qubit Rotations?

## Problem Class
OPT | ANA | DES

## Motivation
The original version of this study tested only a single standard Gaussian multitone SQR pulse against arbitrary block-diagonal conditional qubit-rotation targets. That was enough to show that the single-pulse ansatz performs poorly on strict arbitrary-control metrics, but it was not enough to justify a broader impossibility claim about echoed or composite SQR constructions.

This extension adds the explicit echoed sequence requested by the user:

\[
\text{half SQR} \rightarrow X_\pi \rightarrow \text{half SQR} \rightarrow X_\pi.
\]

The purpose of the extension is to determine whether the earlier negative result was primarily a failure of the single Gaussian multitone ansatz, or whether the explicit echoed composite gate materially changes the conclusion.

Throughout this study we follow the `cqed_sim` qubit-first convention: Hilbert-space objects are interpreted as qubit \(\otimes\) cavity, and the restricted logical basis is ordered as \((|g,0\rangle, |e,0\rangle, |g,1\rangle, |e,1\rangle, \ldots)\).

## Goals
1. State clearly what the previous report did and did not test.
2. Preserve the original single-pulse Gaussian multitone SQR results as the baseline comparison.
3. Implement the explicit echoed schedule `half SQR -> X_pi -> half SQR -> X_pi` and test it directly.
4. Compare single-pulse and echoed SQR under two fairness conventions:
   - fixed total gate duration,
   - fixed active SQR duration.
5. Quantify whether any echoed improvement is mainly due to cancellation of block-dependent residual-Z error or to deeper changes in controllability.
6. Report exact best-fit parameters, waveform plots, machine-readable artifacts, validation outputs, and a revised report that does not overclaim.

## Methods
- `scripts/run_study.py` preserves the original 192-case single-pulse Gaussian multitone baseline across target families A/B/C/D.
- `scripts/run_echo_comparison.py` reuses the validated single-pulse baseline rows for families C and D, then adds two explicit echoed branches:
  - `echoed_fixed_total`: total gate duration held equal to the baseline single pulse, so each half-SQR segment uses \((T - 2\tau_\pi)/2\).
  - `echoed_fixed_active`: total active SQR time held equal to the baseline single pulse, so each half-SQR segment uses \(T/2\) and the total gate becomes \(T + 2\tau_\pi\).
- The echoed construction uses identical first and second half-SQR waveforms. The second half is not phase-shifted or independently reparameterized.
- The inserted `X_pi` pulses are finite Gaussian qubit pulses, not instantaneous ideal kicks. They use:
  - axis: `x`,
  - phase: `0`,
  - duration: `40 ns`,
  - same qubit drive channel as the multitone pulse,
  - carrier chosen from the `n = 0` qubit transition frequency in the study frame.
- For the ideal blockwise design model, the half-pulse target blocks are chosen as
  \[
  W_n = X_\pi^\dagger \sqrt{R_n},
  \]
  so that the time-ordered composite sequence reproduces the requested final block \(R_n\) under an ideal echoed picture.
- The echoed extension focuses on:
  - family C: the strongest structured target family from the baseline study,
  - family D: the decisive random-target stress test.
- Echo extension sweep:
  - `N_active = 2, 3, 4`
  - `|chi| T / 2 pi = 1, 3, 5`
  - model variants: `chi_only`, `chi_plus_chiprime`
  - family D random ensemble: the same 4 seeds per operating point as the baseline subset.
- Combined comparison dataset size:
  - 90 baseline single-pulse rows,
  - 180 echoed rows,
  - 270 total comparison rows.

## Expected Outcomes
- If the earlier negative result was mainly a refocusable blockwise phase problem, the explicit echoed sequence should reduce residual-Z error and improve fidelity on matched cases.
- If the limitation is deeper than residual-Z accumulation, echo may reduce residual phase on some cases without fixing the transverse coherent mismatch that limits full-block fidelity.
- A scientifically honest final conclusion should distinguish:
  - failure of the single Gaussian multitone ansatz,
  - failure of the explicitly tested echoed construction,
  - untested richer composite or independently parameterized waveform families.

## Known Limitations
- The echoed construction tested here is the exact repeated-half schedule requested by the user. The two half-SQR pulses are identical; they are not independently optimized.
- The inserted `X_pi` pulses are finite Gaussian pulses, not ideal instantaneous rotations. Any residual manifold dependence of those pulses is part of the tested physical implementation.
- The echo extension focuses on the representative structured family C and the decisive random family D; it does not resweep families A and B.
- The strict optimization loop still uses a two-level qubit model (`n_tr = 2`) and omits cavity Kerr, open-system noise, and transmon leakage in the main optimization.
- Failure of both the single-pulse and this specific echoed construction still does not rule out richer segmented or independently parameterized composite waveforms.

## Suggested Upstreaming
- Add a public `cqed_sim` helper for time-ordered composite targeted-subspace sequences, not only single conditioned multitone waveforms.
- Add a reusable echoed-target helper that maps a desired block operator `R_n` to a half-sequence target `W_n` under a specified refocusing pulse convention.
- Expose direct blockwise residual-Z and transverse coherent-error decomposition in the public targeted-subspace analysis API.

## Status
COMPLETE