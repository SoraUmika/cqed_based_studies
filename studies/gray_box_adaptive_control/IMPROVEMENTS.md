# Improvement Log: Gray-Box Adaptive Control for cQED Systems

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- **Chi Ramsey probe not in cqed_sim**: The multi-Fock dispersive Ramsey probe (measuring qubit oscillation frequency shift vs Fock number) is implemented analytically in `scripts/probe_library.py`, not in `cqed_sim.calibration_targets`. Upstreaming needed for reuse.
- **Single-parameter control update**: The probe fit returns both χ and χ', but the production gray-box loop only feeds back χ into the learner model. A fully consistent update should jointly incorporate χ', Kerr, and eventually noise parameters; otherwise the omission study still leaves a residual ~0.1-0.2% fidelity gap.

## Recommended Improvements (P2)
- **Bayesian parameter estimation**: Replace the point-estimate L-BFGS-B inference with online Bayesian inference (e.g., sequential Monte Carlo). This provides uncertainty quantification and adaptive measurement scheduling.
- **Multi-mode extension**: Extend to systems with readout cavity and/or second storage mode. The probe protocol needs generalization for cross-talk between modes.
- **Noisy GRAPE objective**: Current GRAPE optimizes coherent fidelity. Optimizing under a Lindbladian objective would make gray-box correction more impactful since it would start from a stronger baseline.

## Nice-to-Haves (P3)
- Hardware validation on a real cQED device.
- Comparison with randomized benchmarking-based adaptive control (Kelly et al. 2016 approach).
- Integration with `cqed_sim.optimal_control.GrayBoxAdaptiveConfig` dataclass for clean API.

## Open Questions
- At what χ mismatch level does the gray-box strategy become insufficient and model-free optimization becomes necessary?
- How does probe efficiency scale with cavity truncation (n_cav > 4)?
- Can the multi-Fock Ramsey data be used to simultaneously track slow drift and trigger recalibration adaptively?

## What Was Tried and Did Not Work
- **Black-box (model-free) differential evolution**: significantly less data-efficient than gray-box approach. Requires many more truth-model evaluations to converge.
- **Omitting chi_higher from probe model**: at chi_higher/chi ~ 1% level, omission causes systematic bias in chi inference that propagates to ~0.5% fidelity degradation.
- **Feeding back inferred chi_higher without Kerr**: a 30% mismatch spot check during repository validation inferred χ'/2π ≈ -21.4 kHz correctly, but re-optimizing on a learner with `(χ_hat, χ'_hat, Kerr = 0)` reduced truth-model fidelity relative to the production `(χ_hat, χ' = 0, Kerr = 0)` update. Partial higher-order correction is therefore not adopted as a drop-in fix; χ', Kerr, and the control objective likely need to be upgraded together.

## Compute & Resource Notes
- GRAPE: 16 time slices × 10 ns = 160 ns total gate. 12-dim Hilbert space (n_cav=4, n_tr=3). Fast per evaluation.
- Phase 5 parameter sweeps: systematic sweeps over readout confusion, probe budget, drift rate, and Hamiltonian omissions. Each sweep is modest compute.
- Data files: `phase4_results.npz`, `phase5_*.npz` in `data/`.

## Resolved
- **Study status is ACTIVE**: resolved during repository-wide consolidation validation. README updated to `COMPLETE`, a dedicated `scripts/validate_results.py` harness now checks artifact completeness, archived-data consistency, truncation sensitivity, and multistart sensitivity, and writes `data/validation_summary.json`.
