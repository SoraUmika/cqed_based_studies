# Corrected SQR Optimization with a Fock-Resolved Effective-Qubit Metric

## Problem Class
OPT | ANA | DES

## Motivation
Recent `cqed_sim` changes corrected the SQR tensor-ordering, carrier-sign, and waveform-amplitude conventions. That makes it necessary to revisit pulse optimization with a metric that matches the actual control objective. The present question is not whether the full joint qubit-cavity unitary matches an ideal block-diagonal operator, but whether each relevant cavity Fock label `n` induces the desired effective qubit rotation parameters `(theta_n, phi_n)`.

The current `cqed_sim.calibration.conditioned_multitone` source now uses the same `phi_n` rotation-axis convention as `cqed_sim.core.ideal_gates.qubit_rotation_xy`. That fixed an earlier subtlety in which the reduced layer had effectively behaved like a Bloch-azimuth target. This study nevertheless keeps the reduced effective-unitary metric as the primary SQR objective, because matching the image of `|g,n>` is still weaker than matching the full effective qubit rotation on manifold `n`.

## Goals
1. Audit the corrected SQR conventions now used by `cqed_sim`.
2. Define a reduced optimization metric that depends only on the Fock-resolved effective qubit action.
3. Compare that metric against narrower alternatives based on final-state Bloch angles.
4. Optimize the corrected Gaussian multitone waveform over multiple pulse durations.
5. Report per-manifold target versus achieved `(theta_n, phi_n)` values, residual errors, and optimized waveform parameters.
6. Assess whether this direct multitone pulse family is expressive enough for the corrected SQR objective.

## Methods
- Primary simulation framework: `cqed_sim`.
- Convention audit sources:
  - `cqed_sim.core.ideal_gates`
  - `cqed_sim.pulses.calibration`
  - `cqed_sim.calibration.conditioned_multitone`
  - `cqed_sim.calibration.sqr`
  - `documentations/physics_conventions.md`
  - `inconsistency/20260327_hilbert_space_ordering_audit.md`
- Analytic preliminary:
  - derive the small-angle reduced multitone kernel
    `beta_n ~= -i sum_m K_nm a_m`
    with
    `K_nm = int_0^T g(t/T) exp[i(Delta_n - Delta_m)t] dt`,
    where `g` is the normalized Gaussian envelope and `a_m` are the complex tone amplitudes;
  - use the kernel rank / conditioning to test first-order reachability before numerical optimization;
  - document the fixed mapping between rotation-axis phase and final-state Bloch azimuth for `|g>` input.
- Numerical workflow:
  - optimize the direct Gaussian multitone waveform
    `epsilon(t) = g(t/T) sum_n A_n exp{i[(omega_n + d_omega_n)t + (phi_n + d_alpha_n)]}`
    with
    `A_n = theta_n/(2T) + lambda0 d_lambda_n`,
    `lambda0 = pi/(2T)`;
  - sweep several durations `|chi| T / 2pi`;
  - compare baseline, reduced-state optimization, and reduced-unitary optimization.
- Framework gap that requires local study code:
  - `conditioned_multitone` has a reduced sector-by-sector state evaluator but does not expose a public reduced **multitone unitary** extractor.
  - `sqr.py` exposes effective-unitary extraction for single-manifold calibration and full multitone block extraction from the full Hilbert-space propagator, but not the reduced multitone unitary metric requested here.
  - This study therefore adds a local helper that reuses the `conditioned_multitone` compiled waveform and evolves the corresponding two-level reduced Hamiltonian with `qutip.propagator`.

## Expected Outcomes
- The reduced effective-unitary metric should be stricter and more aligned with the requested SQR objective than the reduced final-state metric.
- The direct multitone Gaussian ansatz should perform best on smooth low-dimensional targets and degrade with larger active Fock windows or shorter durations.
- If the pulse family is fundamentally limited, state-based optimization may look better than unitary-based optimization because it ignores wrong action on qubit superpositions.

## Known Limitations
- The main optimization remains closed-system and uses a two-level qubit model in the reduced metric.
- The study focuses on the corrected direct Gaussian multitone parameterization rather than richer segmented envelopes.
- The reduced-unitary helper is local study code because the exact public API is not yet available in `cqed_sim`.
- The duration sweep is finite, so the reported best point is only within the scanned operating window.

## Suggested Upstreaming
- Add a public `cqed_sim.calibration.conditioned_multitone` helper that returns the reduced per-sector effective qubit unitary, not only the final state from `|g>`.
- Add an explicit conversion helper between SQR axis phase and Bloch-azimuth diagnostics.
- Add a built-in reduced multitone optimizer against effective-unitary process fidelity.

## Status
COMPLETE
