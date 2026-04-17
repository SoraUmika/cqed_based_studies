# Strong Validation of SQR / CPSQR for Arbitrary Fock-Conditional Qubit Rotations

## Problem Class
OPT | ANA | DES

## Motivation
This study asks a stronger question than the earlier multitone arbitrary-rotation report: can the SQR / CPSQR control families available in the patched current `cqed_sim` realize arbitrary blockwise qubit rotations on an active set of cavity Fock manifolds, and if not, what fails first?

The study is designed around explicit operator-action validation rather than aggregate fidelity alone. In particular, it distinguishes strict blockwise operator agreement from a relaxed per-block `Z`-gauge equivalence class, and it validates candidate controls on basis states, qubit superpositions, and at least one cross-manifold superposition.

Throughout, the Hilbert-space convention is qubit tensor cavity, and the active logical basis is ordered as
\[
(|g,0\rangle,\ |e,0\rangle,\ |g,1\rangle,\ |e,1\rangle,\ \ldots).
\]

## Goals
1. Audit and document the patched `cqed_sim` conventions relevant to arbitrary Fock-conditional control.
2. Define strict SQR and relaxed CPSQR target classes for arbitrary blockwise qubit rotations.
3. Compare multiple waveform families:
   1. direct single-pulse Gaussian multitone,
   2. independent-half echoed multitone,
   3. symmetry-constrained echoed variants,
   4. segmented composite conditional control,
   5. a higher-expressivity benchmark family.
4. Separate failure modes into blockwise angle error, axis error, residual conditional `Z`, off-block contamination, and addressed-subspace leakage.
5. Validate every highlight case on a spanning qubit input set within each addressed block and on a cross-block cavity superposition.
6. Produce publication-quality figures, saved artifacts, a reproducibility notebook, and a final report with conservative claims.

## Methods
- Use only the patched current `cqed_sim` checkout at `../cQED_simulation/`.
- Reuse `cqed_sim.calibration.targeted_subspace_multitone` as the optimization and validation backbone.
- Build arbitrary blockwise SU(2) targets explicitly, then optimize SQR-style multitone control families against:
  - strict restricted-operator fidelity,
  - relaxed per-block `Z`-gauge fidelity,
  - state-transfer validation on a spanning state set.
- Organize the study in stages:
  1. convention audit and targeted package tests,
  2. structured-grid screening,
  3. stress and random-target validation,
  4. stronger-model replay on highlight cases,
  5. report and notebook generation.

### Analytic Preliminary
The dispersive multitone drive on manifold `n` produces an effective qubit Hamiltonian of the form
\[
H_n^{\mathrm{eff}}(t)
\approx
\frac{1}{2}\Omega_{x,n}(t) X
+
\frac{1}{2}\Omega_{y,n}(t) Y
+
\frac{1}{2}\Delta_n(t) Z
\]
with block-dependent transverse drive and residual detuning / Stark-like `Z` terms. Single-segment SQR control is naturally aligned with in-plane rotations, while arbitrary blockwise SU(2) control generally requires either cancellation or absorption of the residual `Z` structure, or additional segmentation that turns the accumulated `Z` terms into a usable gauge freedom.

This motivates the three headline criteria used in the study:
1. strict arbitrary blockwise SU(2) agreement,
2. relaxed per-block `Z`-gauge agreement (CPSQR),
3. state-transfer-only success on a finite validation set.

## Expected Outcomes
- Easier in-plane target subclasses may be reachable with direct SQR-style families.
- General structured SU(2) targets may require segmentation or more expressive envelopes.
- If the benchmark family significantly outperforms Gaussian SQR on the same Hamiltonian, the limitation is primarily ansatz rigidity.
- If even the benchmark family fails on a target subclass, the evidence points toward a deeper limitation under the modeled Hamiltonian and duration constraints.

## Known Limitations
- The study uses closed-system optimization in the main sweep. Stronger-model replays are limited to selected highlight cases.
- The higher-expressivity benchmark is a richer constrained waveform ansatz rather than a full GRAPE implementation.
- Random-target statistics are concentrated on one representative regime rather than every operating point.
- Echo fairness is enforced primarily by fixed active multitone duration; not every family is benchmarked under every fairness convention.
- The saved dataset intentionally emphasizes the strongest informative slices rather than a brute-force Cartesian sweep over every family/target/model combination.

## Suggested Upstreaming
- A public helper for arbitrary blockwise SU(2) target specification, not only conditioned `(theta_n, phi_n)` seeds.
- A reusable relaxed per-block `Z`-gauge metric in `cqed_sim.calibration.targeted_subspace_multitone`.
- Built-in composite-sequence helpers for echoed and segmented targeted-subspace optimization.

## Status
COMPLETE
