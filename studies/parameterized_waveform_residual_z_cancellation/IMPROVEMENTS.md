# Improvement Log: Parameterized Waveform Design for Fock-Conditional Qubit Rotations with Residual-Z Cancellation

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- **[P1 | MEDIUM] No upstream arbitrary-envelope targeted-subspace optimizer**: `cqed_sim` can optimize only the built-in Gaussian conditioned multitone family directly. The richer families in this study must therefore be optimized locally while still replaying through the same simulation stack.

## Recommended Improvements (P2)
- **[P2 | MEDIUM] Replace the current global echo with a manifold-aware echoed protocol**: the implemented echoed family uses one global mid-sequence `X_pi` pulse and helps many random targets, but it catastrophically degrades the smooth structured family. The next version should use manifold-selective echo structure or jointly optimized auxiliary qubit rotations.
- **[P2 | MEDIUM] Promote the objective from scalarized residual-Z suppression to a fidelity/error Pareto study**: the current local objective still emphasizes weighted targeted-subspace loss plus residual-Z regularization. The expanded sweep shows that reducing residual-Z alone does not explain the best families; transverse coherent error must be treated as a first-class competing objective.
- **[P2 | LOW] Extend the random-target ensemble and active-window scaling**: the current grid reaches `N_active = 4` with three random seeds per configuration. The next production pass should push to `N_active = 5` and at least `10` random seeds per duration/model point.
- **[P2 | MEDIUM] Add formal convergence and noise validation**: the current conclusions rely on fixed numerical settings, state-level checks, and cross-family comparison. Dedicated timestep/truncation sweeps and open-system follow-up are still missing.

## Nice-to-Haves (P3)
- **[P3 | LOW] Store compressed aggregate waveform archives**: per-case `NPZ` exports now exist, but a bundled archive keyed by family and representative case would make downstream plotting cheaper.

## Open Questions
- Why does the explicit echoed protocol help most of the hard random targets while failing so badly on the smooth structured family? The current data suggest that the answer is not just over-rotation of the echo pulse.
- Are the echo-driven gains on random targets robust, or are they specific to the present small random ensemble and manifold-0 echo calibration?
- Does `chi'` mainly degrade performance through spectral crowding, or through a symmetry-breaking mechanism that specifically frustrates Z cancellation?

## What Was Tried and Did Not Work
- **Global echoed family on structured targets**: the explicit mid-sequence global `X_pi` echo is not a universal fix. On the structured family `C`, its mean fidelity falls to `0.164414` versus `0.763677` for the baseline, while residual-Z and transverse error both worsen. This version of the echo should not be used as the default structured-target protocol.
- **Residual-Z-only interpretation of the hard random cases**: the expanded grid confirms that residual-Z reduction is only part of the story. The echoed family lowers the random-target mean residual-Z error by about `0.353` rad relative to the baseline and also lowers mean transverse error by about `0.200` rad, while the best single-segment families change both metrics only weakly.
- **Symmetric two-segment family as a generic upgrade**: this family still gives slightly lower residual-Z on average, but its fidelity is consistently a bit worse than the best single-segment refinements and much worse than the echoed family on the hard random cases.

## Compute & Resource Notes
- Single expanded structured case (`5` waveform families): about `18` to `60` seconds wall clock depending on `N_active` and model variant.
- Expanded random cases are slower and dominate the runtime: about `18` to `67` seconds per case in the present `N_active = 2, 3, 4` sweep.
- Full expanded grid (`48` cases, `240` waveform-family rows): about `35` minutes wall clock on the current workstation.
- Promoting the residual/transverse diagnostics into the baseline study did not require re-running the old optimization campaign; regenerating summaries and figures from the saved artifacts is fast enough for routine comparison updates.

## Resolved
- **[P2 | MEDIUM] True echoed multitone protocol added**: implemented as `echoed_multitone` with an explicit mid-sequence Gaussian `X_pi` pulse. The current implementation is informative but not yet the final protocol.
- **[P2 | LOW] Larger active-window sweep completed**: the study now covers `N_active = 2, 3, 4` with three random seeds per duration/model point.
- **[P2 | MEDIUM] Baseline residual/transverse diagnostics promoted**: `multitone_sqr_arbitrary_fock_conditional_rotations` now records and plots the same coherent-error decomposition for direct comparison.
- **[P2 | MEDIUM] Residual-Z and transverse error separated in aggregates**: the updated figures and machine summaries now track both components explicitly.
- **[P3 | LOW] Waveform spectrum plots added**: the study now saves representative baseband spectrum comparisons across the waveform families.
- **[P3 | LOW] Pulse samples exported in `NPZ`**: each case/family artifact now has a companion waveform file in `artifacts/waveforms/`.
- **Initial optimizer-bounds warning**: resolved by moving the symmetric two-segment initial guess inside its amplitude-scale bounds.