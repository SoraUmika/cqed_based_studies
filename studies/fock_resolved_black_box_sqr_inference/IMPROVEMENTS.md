# Improvement Log: Fock-Resolved Qubit-State Inference for a Black-Box SQR Gate

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- **[P1 | HIGH]** The stated measurement set `{qubit tomography, cavity displacement, dispersive wait}` has an analytic nullspace: it cannot identify `p_n` and `Z_n` separately. Any future extension that claims full Fock-resolved qubit-state recovery must add a new measurement primitive or an explicit prior that removes this degeneracy.
- **[P1 | HIGH]** The desired inverse problem is not a first-class `cqed_sim` workflow. The study will need a local inference layer for the wait/displacement tomography model unless a dedicated upstream helper is added.
- **[P1 | HIGH]** Combined-protocol fits must not be interpreted literally when the residual is large. In the completed study, coherent and leakage-prone black-box outputs produced large combined residuals and large weighted-sector errors precisely because the diagonal inverse model had failed.

## Recommended Improvements (P2)
- **[P2 | MEDIUM]** Upstream a reusable `black_box_fock_inference` helper into `cqed_sim.tomo` or `cqed_sim.calibration_targets`, including kernels for `D(alpha) -> wait` and a documented statement of the recoverable subspace.
- **[P2 | MEDIUM]** Add Bayesian model comparison between the diagonal model and a low-rank cavity-coherence surrogate, so combined-protocol residuals become a quantitative coherence witness rather than only a heuristic flag.
- **[P2 | MEDIUM]** Add an optional auxiliary measurement primitive, such as cavity parity or a Fock-selective qubit tag pulse, and repeat the identifiability analysis to show exactly which null directions disappear.

## Nice-to-Haves (P3)
- **[P3 | LOW]** Expand the pulse-level black-box case library with more waveform families and duration points once the core identifiability conclusion is stable.
- **[P3 | LOW]** Add a benchmark against `cqed_sim.tomo.run_fock_resolved_tomo` as an oracle-with-extra-control comparison, clearly labeled as a stronger measurement setting than the black-box protocol under study.

## Open Questions
- If independent cavity-population information is supplied from a separate calibration, does the remaining wait/displacement protocol become sufficient to estimate `Z_n`, or is an additional qubit-population-coupling operation still required?
- What is the minimum combined-protocol alpha-grid needed to make cavity coherences reliably detectable in the presence of realistic calibration drift?
- How much of the pulse-level failure budget in the imperfect black-box cases comes from leakage versus coherent in-subspace error?

## What Was Tried and Did Not Work
- **Displacement-only tagging intuition**: a cavity displacement applied immediately before tracing out the cavity leaves the reduced qubit state invariant. Scratch calculations confirmed that `t = 0` makes every displacement-only setting exactly identical to ordinary qubit tomography.
- **Direct use of `run_conditioned_multitone_validation(..., simulation_mode=\"full\")` with `n_tr = 3`**: scratch testing failed because the helper compares 2x2 target states against 3x3 reduced transmon states. This is not a blocker for the current study because waveform construction still works and pulse-level validation can be handled with `n_tr = 2` or with explicit qubit-subspace truncation in local code.

## Compute & Resource Notes
- Study initialization and API review completed on 2026-03-31.
- Scratch timing: a small full-mode conditioned-multitone optimization for a 4-sector target completed in under one minute and reached weighted fidelity essentially equal to unity on an `n_tr = 2`, `n_cav = 4` model.
- Planned runtime split: one optimized pulse-level case (slow), one nominal pulse-level replay (moderate), all identifiability / MLE sweeps (fast, vectorized).
- Full study runtime (`STUDY_PROFILE=full`): approximately 105 s on the active workstation, including figure generation and one optimized multitone pulse-level case.
- Validation summary: all automated checks passed. The exact wait-only diagonal benchmark and the exact combined diagonal benchmark fit to machine precision; the coherent combined benchmark produced a large residual (`0.1387`), confirming the intended coherence-witness behavior.

## Resolved
- **`cqed_sim` access path**: resolved locally. The active installation is `C:\Users\jl82323\Box\Shyam Shankar Quantum Circuits Group\Users\Users_JianJun\cQED_simulation`; the stale `C:\Users\dazzl\...` path in some documentation is not present in this environment.
- **Wait-phase sign convention in the recoverable kernel**: resolved during implementation. The exact joint-state simulator agreed with the analytic inverse model only after using the transverse phase factor `e^{-i\phi_n(t)}` rather than `e^{+i\phi_n(t)}`.
- **Full-state non-identifiability demonstration**: resolved with two complementary checks. Random-restart full-state MLE produced widely separated population assignments, and an explicit gauge-family construction produced multiple exact sector decompositions with identical recoverable data.
