# Execution Plan

## Task Summary
- Study: `studies/waveform_level_gate_realization_dispersive_cqed`
- Run: `task_runs/waveform_level_gate_realization_dispersive_cqed`
- Prompt date: 2026-04-01
- Problem class: `ANA`, `DES`, `OPT`
- Objective: determine how well different cavity-control strategies approximate `I_q \otimes D(alpha)` in a realistic dispersive cQED model, and identify the best practical protocol.

## Protocol Order
1. Establish the baseline failure with a naive single-tone pulse.
2. Test the "make it fast" strategy using short square, Gaussian, and cosine pulses.
3. Test direct compensation with a two-tone branch-matched pulse.
4. Test an echoed displacement using inserted qubit `pi` pulses.
5. Run a bounded hardware-aware optimal-control benchmark.
6. Compare all families on a common state-test set and produce a final recommendation.

## Model Hierarchy
1. Minimal dispersive model with `chi` only
2. Higher-order model with `chi` and `chi'`
3. Full model with `chi`, `chi'`, and cavity self-Kerr `K`

## Core Metrics
- `delta_alpha = |alpha_g - alpha_e|`
- Superposition-state entanglement entropy and branch overlap
- State fidelities against the ideal displaced target on the explicit cavity-state test set
- Wigner-function agreement for representative success and failure cases
- Duration, bandwidth, and implementation-complexity comparisons

## Validation Gates
- Sanity: recover the nearly unconditional short-pulse limit and verify that `chi` dominates the long-pulse failure.
- Convergence: rely on the already validated propagator settings for this study family and keep the same default truncation and step size.
- Literature alignment: ensure the observed breakdown tracks the expected `|chi| T` criterion.

## Deliverable Checklist
- `README.md` and `IMPROVEMENTS.md` updated for unconditional displacement
- `SCIENCE_DIRECTIVE.md` written
- unconditional-displacement artifacts and figures generated
- `report.tex` and `report.pdf` refreshed around the new protocol comparison
- reproducibility notebook regenerated against the new artifacts
- `EXECUTION_SUMMARY.md` and `REVIEW_REQUEST.md` updated for review handoff
