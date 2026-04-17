# Improvement Log: Native / Rich Multitone Feasibility for Ideal SQR and CPSQR

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- **[P1 | MEDIUM] Best-case higher-level validation not yet included in the main grid**: The principal study uses `n_tr = 2`. Follow the strongest strict-SQR and CPSQR cases with `n_tr = 3` and, if warranted, cavity Kerr to check whether the conclusion is robust.
- **[P1 | MEDIUM] No dedicated duration refinement for `N_active = 4`**: The comparison grid is strong enough to show that strict full-joint SQR degrades substantially at `N_active = 4`, especially with `chi_plus_chiprime`, but the current study does not include a fine duration threshold sweep there.

## Recommended Improvements (P2)
- **[P2 | MEDIUM] Add open-system follow-up for best cases only**: Once the closed-system answer is stable, test whether the best direct and echoed families survive realistic `T1/Tphi`.
- **[P2 | MEDIUM] Add multi-restart statistics for the richest sampled-envelope families**: The current workflow uses bounded restart counts to keep the study tractable. More restarts would strengthen negative claims about local minima.
- **[P2 | LOW] Upstream the corrected helper layer**: The study will carry a local corrected run-config helper because several earlier study-local wrappers passed frame-shifted `fock_fqs_hz` values into the patched absolute-frequency API.

## Nice-to-Haves (P3)
- **[P3 | LOW] Add a stronger operator-norm proxy**: The current plan uses spectral/Frobenius proxies on the addressed subspace. A tighter diamond-like bound would sharpen the full-joint conclusion.
- **[P3 | LOW] Add richer random-x ensembles**: The stage-2 random target set is intentionally small to keep the workflow interpretable.

## Open Questions
- Does the best CPSQR case remain predominantly a post-rotation manifold-dependent `Z` correction, or does it hide a more general near-SU(2) mismatch?
- Are any apparent echoed gains symmetry-protected, or are they mainly optimization-budget dependent once finite `X_pi` pulses are modeled?
- In the structured `N_active = 2` screen cases at `|chi| T / 2pi = 3, 5`, several echoed variants reached near-unit CPSQR joint fidelity while their strict joint fidelity remained poor. Determine whether this is a clean CPSQR realization or an artifact of the current echo-optimizer parameterization.
- Can a composite objective that explicitly optimizes the strict inter-manifold phase relations recover the echoed CPSQR benefit while also lifting strict joint fidelity above the current `~0.87` ceiling?

## What Was Tried and Did Not Work
- The fully symmetric echoed replay of the strict half-SQR solution performs poorly as a strict-SQR construction in the early structured screen cases. It consistently underperforms the direct strict families on strict addressed-subspace fidelity, so symmetry alone is not enough once finite inserted `X_pi` pulses are included.
- The richer sampled-envelope direct families are not automatically better than the native strict family. In the early `N_active = 2` smooth-target screen, the symmetric two-segment family helped at `|chi| T / 2pi = 1`, but the complex-envelope and basis-expanded families were noticeably weaker than the strict direct optimizer at `|chi| T / 2pi = 3, 5`.

## Compute & Resource Notes
- Upstream convention-audit regression tests run successfully on March 27, 2026:
  - `tests/test_25_tensor_product_convention.py`
  - `tests/test_23_sqr_additive_amplitude_correction.py`
  - `tests/test_20_gaussian_iq_convention.py`
- The audit also confirmed that passing frame-relative `fock_fqs_hz` values into the patched `build_sqr_tone_specs(...)` helper shifts the internal carriers back to the bare qubit frequency and must not be reused in the new study workflow.
- The staged screen run started on March 27, 2026 under `python -u studies/native_rich_multitone_sqr_cpsqr_feasibility/scripts/run_study.py` so progress could be monitored live from the log.
- Runtime observation from the screen stage: direct families are inexpensive enough for the staged grid, while independently corrected echoed variants are the main wall-clock bottleneck. The full screen remains tractable, but brute-force multi-restart echo studies would likely need a narrower case selection.
- The final artifact set was completed with `python -u studies/native_rich_multitone_sqr_cpsqr_feasibility/scripts/resume_study.py` plus a separate echoed-duration completion run under `python -u studies/native_rich_multitone_sqr_cpsqr_feasibility/scripts/run_echoed_duration.py`.
- Total saved workload of the final study: 120 screen evaluations, 240 comparison evaluations, 72 direct duration evaluations, and 36 echoed duration evaluations.

## Resolved
- **Runtime import mismatch for background runs**: Added `scripts/runtime_compat.py` and imported it early in the study entrypoints so background or resumed runs can still discover the sibling patched `cqed_sim` checkout.
- **Stale `fock_fqs_hz` helper mismatch**: The study-local run-config helper now leaves `fock_fqs_hz=None` and documents that patched `cqed_sim` interprets the override as an absolute transition frequency. The issue is also recorded in `BLOCKERS.md`.
