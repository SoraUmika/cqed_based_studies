# GRAPE-Like Optimization of cQED Readout Pulses

## Problem Class

OPT | ANA | DES

## Motivation

Dispersive readout in circuit QED is the canonical method for measuring the state of a
superconducting qubit. The readout pulse drives the coupled resonator, whose response
differs depending on the qubit state due to the dispersive shift χ. Optimizing the
readout pulse shape is therefore directly relevant to measurement fidelity, speed, and
the QND character of the readout process.

This study systematically addresses what constitutes the most physically meaningful and
practically useful formulation for readout pulse optimization — comparing analytical
pulse families, piecewise-constant GRAPE-like control, multi-objective formulations,
open-system effects, and robustness to parameter uncertainty.

## Goals

1. Derive and validate the linear dispersive resonator model for qubit-state-conditioned
   field evolution α_g(t) and α_e(t) under arbitrary drive ε(t).
2. Establish the optimal drive frequency and its dependence on χ/κ.
3. Compare five pulse families (square, Gaussian, Hann, cosine-rise, spline) across a
   sweep of χ/κ ∈ {0.25, 0.5, 1.0, 2.0, 4.0} for integrated SNR² and endpoint separation.
4. Implement GRAPE-like adjoint-gradient optimization for piecewise-constant pulses and
   compare against the best pulse-family baseline.
5. Study the joint readout + active depletion optimization and quantify residual photon
   suppression.
6. Assess the effect of qubit T1 relaxation and measurement-induced dephasing on net
   assignment fidelity as a function of readout duration.
7. Evaluate robustness of the nominally optimal pulse under ±20% uncertainty in χ and κ,
   and determine whether robust optimization is warranted.
8. Produce a final ranked recommendation of control formulation, waveform family, and
   optimization strategy for practical use.

## Methods

### cqed_sim modules used

- `cqed_sim.measurement.ReadoutResonator` — steady-state amplitude α_ss, time-domain
  `response_trace()` for constant-ε validation, and `gamma_meas()` for dephasing rate.
- `cqed_sim.measurement.ReadoutChain` — integrated I/Q trace and `iq_centers()` for
  steady-state I/Q separation validation.
- `cqed_sim.pulses.envelopes` — `square_envelope`, `gaussian_envelope`,
  `cosine_rise_envelope` imported for cross-checking against standalone implementations.

### Standalone extensions (documented gap in cqed_sim)

**Gap 1**: `cqed_sim.measurement.ReadoutResonator.response_trace()` only supports
constant-amplitude drives (solves analytically α_ss + (α_0 − α_ss) e^{−λt}).
Arbitrary time-varying envelopes ε(t) require a numerical ODE integrator.
→ `common.py::integrate_readout_ode()` implements a 4th-order Runge-Kutta solver.

**Gap 2**: `cqed_sim.optimal_control` (GrapeSolver) is designed for unitary/state-
transfer objectives, not readout-specific objectives (SNR, pointer separation, residual
photons). → `common.py::grape_adjoint_gradient()` implements the analytical adjoint
gradient for the linear dispersive model.

Both gaps are candidates for upstreaming into cqed_sim (see Suggested Upstreaming).

## Assumptions

- Physical model: single-mode dispersive readout resonator (linear dispersive limit).
  Rotating frame at the drive frequency ω_d; EOM: dα/dt = −(κ/2 + iΔ)α − iε(t).
- Dispersive shift: χ/(2π) ∈ {−1.25, −2.5, −5.0, −10.0, −20.0} MHz for χ/κ ∈ {0.25,
  0.5, 1.0, 2.0, 4.0} at fixed κ/(2π) = 5.0 MHz.
- Maximum drive amplitude: |ε_max|/(2π) = 5.0 MHz (≤ 5 steady-state photons on resonance).
- Qubit decoherence: T1 = 30 μs, T2 = 20 μs (pure dephasing T_φ = 60 μs).
- ODE integration time step: dt = 2 ns; GRAPE segments: N = 60.
- Hilbert space truncation: not applicable (linear model has no truncation).
- Phase 4 open-system: qubit T1 and measurement-induced dephasing treated analytically.
- Convergence criterion: SNR² stable to 0.1% when dt halved.

## Expected Outcomes

- Square pulse at midpoint drive (ω_d = ω_r + χ/2) is already near-optimal for SNR²
  among simple families; Hann/Gaussian shapes add <5% improvement.
- GRAPE improves SNR² by ~10–20% over square pulse at intermediate κT (≈ 2–5).
- Joint readout + depletion optimization reduces residual photons by >10× compared with
  readout-only.
- Qubit T1 limits useful readout duration to T ≲ 0.3 T1; optimal T_read ≈ 3–5/κ.
- Nominal GRAPE pulse is fragile to >10% χ uncertainty; robust GRAPE recovers ≈70% of
  nominal gain with 5-point ensemble.
- The most practically useful formulation: time-integrated SNR² objective, midpoint
  drive, Hann or square envelope for simple use, GRAPE for precision readout.

## Status

COMPLETE — Validated on 2026-03-20. In the present linear, peak-limited model, the recommended strategy is a midpoint-driven square pulse with $\kappa T \approx 5$--$6$; GRAPE matches but does not exceed square, and the current robust-GRAPE and active-depletion formulations are not recommended.

---

## Suggested Upstreaming

1. **`ReadoutResonator.simulate_drive(epsilon_t, tlist)`**: general ODE solver for
   time-varying drive, returning conditioned field traces α_g(t), α_e(t).
2. **Readout GRAPE module**: `solve_readout_grape(resonator, T, N_seg, objective)`
   implementing the adjoint gradient for SNR²/multi-objective readout costs.
3. **`ReadoutChain.snr2(alpha_g, alpha_e, tlist)`**: convenience function for
   computing the heterodyne SNR² from conditioned field traces.
