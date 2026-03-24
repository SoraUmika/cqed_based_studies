# Improvement Log: Open-System Deep Dive for Selective Qubit Rotation in Dispersive cQED

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- **None currently at the archive-validity level**: the replay blocker was resolved by
  switching GRAPE evaluation to the simulator-backed control-schedule path. Remaining work is
  quantitative scope expansion rather than a correctness blocker.

## Recommended Improvements (P2)
- **Extend the three-mode study to nonzero `chi_sr`**: the current conservative baseline
  uses `chi_sr = 0`, so the spectator-readout result is likely a lower bound on direct
  storage-readout coupling effects.

## Nice-to-Haves (P3)
- **Persist separate reduced-fidelity and coherence convergence checkpoints**: the three-mode
  coherence-definition fix did not affect the reduced-fidelity path, but the current
  partial-file schema forced both metrics to be invalidated together.
- **Store explicit GRAPE replay probe tables**: save the representative replay-versus-substeps
  scan to a dedicated data file instead of relying on log output.

## Open Questions
- How much of the earlier GRAPE replay instability was due to pulse-export semantics versus
  the coarser replay grid itself?
- Would a hardware-smoothed parameterization or held-sample export materially reduce replay
  cost while preserving the optimized closed-system objective?

## What Was Tried and Did Not Work
- **Validator bug fix alone was not sufficient**: the original validator incorrectly compared
  a fresh re-optimization against the archived control. After switching to replay of the
  archived control itself, the representative noisy-fidelity delta remained `6.50e-2`
  between `16` and `32` substeps per slice, confirming a real replay-stability problem.
- **Higher replay resolution without changing the replay method did not converge cleanly**:
  the representative archived control produced noisy target fidelities `0.8320`, `0.8970`,
  `0.9137`, `0.9077`, `0.9179`, and `0.9189` at `16`, `32`, `64`, `96`, `128`, and `256`
  substeps per slice.
- **Removing the extra end padding did not fix the instability**: replaying the same control
  with `t_end = duration` instead of `duration + 4 dt` still yielded `0.8206`, `0.8940`,
  `0.9130`, and `0.9177` at `16`, `32`, `64`, and `128` substeps per slice.
- **Forcing the solver `max_step` to match the replay grid did not materially change the
  result**: the no-padding replay with `SimulationConfig(max_step = dt)` matched the
  unstable sequence above to numerical precision.

## Compute & Resource Notes
- Full `run_grape_noisy_replay.py` regeneration with the stabilized
  `evaluate_control_with_simulator(...)` replay path completed successfully for all four
  representative durations in about `95 s`.
- The corrected three-mode validation baseline and `dt` refinement were successfully
  recomputed and checkpointed; the larger truncation rerun was much more expensive in
  wall-clock time than the other convergence legs.
- The earlier three-mode truncation reduced-fidelity scalar `0.7260030664` was preserved
  from the pre-fix partial state because the coherence-definition change did not affect that
  reduced-fidelity code path. The corrected truncation-point coherence ratio was not
  recomputed before the study hit the GRAPE blocker and is not part of the pass/fail gate.
- The stabilized representative GRAPE replay at `|chi_s| T / 2pi = 2` gives noisy target
  fidelities `0.958318`, `0.958278`, and `0.958358` at `512`, `768`, and `1024` replay
  substeps per slice, respectively, so the saved archive now comfortably meets the
  `5e-4` convergence target.

## Resolved
- **GRAPE noisy replay was not numerically converged**: resolved by replacing the old
  pulse-export replay path with `evaluate_control_with_simulator(...)` and regenerating the
  GRAPE archive plus convergence artifacts. The representative `512` to `1024` substep
  replay delta is now `3.97e-5`.
- **The GRAPE archive and figure were stale with respect to the stable replay method**:
  resolved by rerunning `run_grape_noisy_replay.py`, `validate_results.py`, and
  `plot_results.py` after the replay-path fix.
