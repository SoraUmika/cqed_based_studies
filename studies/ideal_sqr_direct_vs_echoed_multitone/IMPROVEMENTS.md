# Improvement Log: Ideal SQR Feasibility with Direct and Echoed Multitone Waveforms

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- **[P1 | MEDIUM] No public multi-segment targeted-subspace optimizer in `cqed_sim`**: the direct single-segment ideal-SQR optimization is fully inside the framework, but the equal-duration and asymmetric echoed cases require a local composite optimization wrapper. If a public composite optimizer lands upstream, replace the local Powell layer immediately.
- **[P1 | MEDIUM] Finite echo-pulse manifold dependence may dominate the echoed verdict**: the echoed construction is only a clean symmetry argument if the inserted \(X_\pi\) pulses act uniformly across the active manifolds. The current study must track whether echo failure is due to the multitone halves or due to the refocusing pulses themselves.

## Recommended Improvements (P2)
- **[P2 | MEDIUM] Add manifold-selective or jointly optimized echo pulses**: if the symmetric echo fails mainly because the global \(X_\pi\) pulses are not uniform across manifolds, the next pass should upgrade the echo pulses before concluding that the composite construction itself is fundamentally weak.
- **[P2 | MEDIUM] Extend the ideal-SQR target family beyond fixed \(x\)-axis rotations**: the present follow-up chooses the \(x\)-axis convention because it matches the user’s \(\theta_n\)-only specification and gives the cleanest analytic echo algebra. A broader follow-up should test arbitrary common XY-axis phases.
- **[P2 | LOW] Push the active-window grid to \(N_{\mathrm{active}}=4\)**: the current main grid stops at \(N_{\mathrm{active}}=3\) for runtime and interpretability. A larger-window extension would test whether the same conclusions persist once spectral crowding is stronger.

## Nice-to-Haves (P3)
- **[P3 | LOW] Export a compact per-case analytic summary table**: the report will state the small-\(n\) formulas, but a machine-readable table linking each numerical case to its direct and echoed analytic expectations would make future comparisons easier.
- **[P3 | LOW] Bundle representative waveform overlays by case family**: storing side-by-side direct/symmetric/independent/asymmetric overlays for the best and worst cases would help future agents inspect qualitative differences quickly.

## Open Questions
- If a symmetric echoed sequence still underperforms after first-order residual-\(Z\) cancellation is analytically available, which term dominates numerically: finite-\(X_\pi\) nonuniformity, second-order commutators, or simple loss of available active-drive time?
- When independent corrections help, do they mainly repair the multitone halves or do they compensate errors introduced by the echo pulses?
- If weak timing asymmetry helps only marginally, is that evidence for robustness, or evidence that the chosen parameterization is already too constrained to exploit the extra degree of freedom?

## What Was Tried and Did Not Work
- **Fully symmetric echoed construction on the full focused grid**: it underperformed the direct multitone waveform in every tested case. Mean average fidelity was `0.248211` versus `0.595142` for the direct construction, and its mean residual-Z error was also worse (`0.318233` rad versus `0.058152` rad). The clean first-order echo algebra did not survive the finite-pulse implementation strongly enough to pay off numerically.
- **Equal-duration independently corrected echoed construction**: allowing separate correction vectors for the two halves changed the answer only in the fifth decimal place. It improved over the symmetric echo in `9/16` cases, but the mean fidelity difference was negligible (`0.2482169` versus `0.2482112`), and it never overtook the direct waveform.
- **Weakly asymmetric echoed construction**: the Powell layer never chose a nonzero asymmetry within numerical tolerance on the focused grid. Every saved best-fit case returned `duration_asymmetry_eta = 0.0`, so the extra timing freedom provided no practical gain in this parameterization.
- **Residual-Z-only echo narrative**: no echoed construction reduced residual-Z error relative to the direct waveform on any saved best case. The smallest residual-Z penalty among echoed cases was still positive (`+0.107005` rad relative to direct on `chi_plus_chiprime_smooth_x_na3_chiT5p0`).

## Compute & Resource Notes
- Full focused grid (`16` cases, `64` saved construction rows) completed in about `407 s` wall-clock on the current workstation with system Python `3.12.10`.
- Validation rerun (sanity + dt convergence + cavity-padding convergence) took about `8 s`.
- The heaviest cases were the `N_active = 3` runs with `chi + chi'`, especially the direct multitone optimization and the echoed local Powell refinements.
- The composite echoed Powell stages routinely exhausted the configured evaluation budget (`70` or `80` function evaluations) without discovering a materially better valley; this is evidence that the current echoed parameterization is flat or genuinely unhelpful, not merely under-optimized by one or two steps.
