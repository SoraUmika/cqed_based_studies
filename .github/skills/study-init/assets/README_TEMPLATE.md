# <Study Title>

## Problem Class

<!-- OPT | REP | DES | ANA — pick one or more -->

## Motivation

<!-- Why this study matters. Link to paper if REP class. -->

## Goals

<!-- Numbered, concrete, falsifiable goals. Example:
1. Achieve gate fidelity > 99.9% for a SNAP gate on the storage cavity.
2. Reproduce Fig. 3 of [Author et al., PRL 2024] within 2% agreement.
-->

## Methods

<!-- Which cqed_sim modules/functions will be used. Example:
- `DispersiveReadoutTransmonStorageModel` for system setup
- `simulate_sequence` for time-domain evolution
- `UnitarySynthesizer` for gate optimization
-->

## Analytic Preliminary

<!-- Start from the first-principles model (Hamiltonian, equations of motion, symmetry argument, or conservation law) whenever feasible.
State the closed-form or limiting-case reasoning attempted before numerics.
List every controlled approximation introduced to simplify the problem and why it is valid.
If no useful analytic answer exists, explain why the problem must be numerical. -->

## cqed_sim Gap Analysis

<!-- Table template:
| Functionality | Needed? | Available in cqed_sim? | Plan |
|---|---|---|---|
| Example feature | Yes | Partial | Use cqed_sim + local helper |
-->

## Assumptions

<!-- Physical assumptions, parameter ranges, convergence criteria. Example:
- Transmon: ω_q = 5.0 GHz, α = -200 MHz
- Storage: ω_s = 7.5 GHz, κ = 10 kHz
- Hilbert space: N_transmon = 4, N_storage = 20
- Convergence: fidelity stable to 1e-4 when doubling N_storage
-->

## Compute & Resource Strategy

<!-- Upfront cost estimate and acceleration plan. Example:
- Single-point simulation: ~10 s
- Parameter sweep: ~20 min over 16 parallel workers
- Long-run bottleneck: GRAPE at 4 pulse durations
- Planned acceleration: joblib.Parallel, JAX backend unavailable
-->

## Expected Outcomes

<!-- What success looks like — quantitative where possible. Example:
- Fidelity ≥ 99.5% for target unitary
- Chi shift matches analytic prediction within 1%
-->

## Known Limitations

<!-- Updated throughout the study. What approximations are being made?
What is constrained by compute time or framework capability? -->

## Validation

<!-- Keep this current as the study progresses:
- [ ] Sanity checks
- [ ] Convergence
- [ ] Literature comparison (if applicable)
-->

## Status

ACTIVE
<!-- ACTIVE | COMPLETE | BLOCKED — update as work progresses -->
