---
description: "Pre-study checklist — run before starting any new study to verify all prerequisites from AGENTS.md are satisfied. Ensures: API Reference consulted, first-principles analytic preliminary attempted, controlled approximations stated, compute cost estimated, gap analysis written, Hilbert space dims justified."
---

# Pre-Study Checklist

Before proceeding with implementation, answer each question below. All items must be addressed — mark N/A only with justification.

## 1. API Reference Consulted?

Have you read the [cqed_sim API Reference](https://github.com/SoraUmika/qubox_cQEDsim/blob/main/API_REFERENCE.md) to confirm:

- [ ] The required simulation functionality exists in `cqed_sim`
- [ ] You know which modules and functions to use
- [ ] Any gaps are documented in the README under `## cqed_sim Gap Analysis`

**Which modules/functions will you use?**
> (list them here)

## 2. Analytic Preliminary Attempted?

Before running numerics, have you:

- [ ] Started from the first-principles model, Hamiltonian, or governing equations when feasible
- [ ] Attempted a closed-form or limiting-case solution
- [ ] Recorded the analytic argument (or why it fails) in README `## Analytic Preliminary`
- [ ] Listed the controlled approximations used and the regime in which they are valid
- [ ] Identified whether the central question is achievable in principle

**Summary of analytic result:**
> (brief description of the first-principles picture, the approximations used, or "no useful analytic result because...")

## 3. Compute Cost Estimated?

- [ ] Estimated wall-clock time for the main computation
- [ ] Decided on acceleration strategy (GPU, multiprocessing, batching) if > few minutes
- [ ] Recorded the plan in README `## Compute & Resource Strategy`

**Estimated cost:** (e.g., "GRAPE sweep over 6 durations, ~10 min per duration, parallelized with joblib")

## 4. Hilbert Space Dimensions Justified?

- [ ] Chosen truncation dimensions for each mode (qubit, cavity, readout)
- [ ] Justified the choice (max photon number, convergence from prior studies, etc.)
- [ ] Planned a convergence check to verify truncation is sufficient

**Dimensions:** (e.g., "N_qubit=4, N_cavity=15, N_readout=N/A")

## 5. Physical Assumptions Stated?

- [ ] Listed all physical assumptions in README `## Assumptions`
- [ ] Identified which approximations are used (RWA, dispersive, etc.)
- [ ] Stated parameter ranges and convergence criteria

## 6. Study Folder Initialized?

- [ ] `studies/<name>/` created with full folder structure
- [ ] README.md has all mandatory sections
- [ ] IMPROVEMENTS.md created with placeholder headings
- [ ] study_state.json created with initial state

If any item is not satisfied, address it before proceeding to implementation.
