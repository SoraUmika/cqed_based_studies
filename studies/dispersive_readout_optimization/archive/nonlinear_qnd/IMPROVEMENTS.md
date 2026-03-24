# Improvement Log: Nonlinear-QND Limits and Hardware-Realistic Optimization of Dispersive Readout Pulses in cQED

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- **Microscopic strong-drive backaction is still not natively exposed in `cqed_sim`**: the
  package now has reusable phenomenological strong-readout disturbance helpers for arbitrary
  envelopes, but the physics remains operational rather than microscopic.
- **No theorem-tight QND-aware discrimination bound**: the planned realistic benchmark will
  be model-based and optimization-backed rather than an information-theoretic proof.
- **The executed frontier is still sparse in duration**: the stable saved frontier covers
  `96`, `240`, and `496` ns. That is enough to answer the control question qualitatively,
  but not enough to claim a fully resolved duration-phase diagram.

## Recommended Improvements (P2)
- **Validate the effective mixing layer against a more microscopic beyond-dispersive model**:
  this remains the highest-priority follow-on if the control conclusion is going to inform
  hardware calibration policy. The present rich-model ranking is plausible, but the mixing
  law is still phenomenological.
- **Search a constrained hybrid family around `ring_hold`**: the rich frontier is led by
  `ring_hold` at all three executed durations, which suggests that a small number of extra
  smooth degrees of freedom may recover more detector-limited fidelity without the QND
  penalty of fully free procedural segments.
- **Promote the hardware-profile presets into reusable shared utilities**: the present
  study-local settings are likely reusable across future readout studies.

## Nice-to-Haves (P3)
- **Add a dedicated three-mode optimization sweep**: useful once the two-mode-plus-mixing
  hierarchy is validated.
- **Add trajectory-based repeated-readout diagnostics**: would separate detector noise from
  state disturbance more cleanly.
- **Run a denser duration sweep around `200-350` ns**: this is the region where the control
  ranking changes from legacy `procedural_segments` dominance to rich-model `ring_hold`
  dominance.

## Open Questions
- Which part of the richer model changes the control conclusion most: occupancy-driven
  mixing, slew-driven mixing, hardware low-pass filtering, or IQ distortion?
- Does the nulling tail lose the frontier mainly because finite bandwidth destroys the exact
  cancellation condition, or because the richer model now rewards smoother occupancy
  trajectories more than exact emptying?
- Can a lightly augmented `ring_hold` family close the remaining detector gap at `240 ns`
  without triggering the high induced-transition cost seen in `procedural_segments`?

## What Was Tried and Did Not Work
- **Initial effective-mixing activation law**: the first version used a conservative onset
  and scaling, and the rich-model smoke test produced effectively zero induced transition
  and zero synthesized mixing amplitude. The onset had to be lowered and the activation
  rescaled before the richer model actually generated family-dependent QND penalties.
- **Powell local refinement in the physical regimes**: local search was stable enough in the
  linear problem, but it was unreliable for the richer physical sweeps, especially for
  `procedural_segments`. Direct random evaluations were fine; the unstable path was the
  local optimizer. The stable workaround was warm-start plus random search only
  (`n_local = 0`).
- **Long silent stage execution**: the larger stagewise runners appeared to stall or be cut
  off when a stage produced no output for a while. Adding heartbeat prints per duration and
  per family was enough to stabilize long runs in this environment.

## Compute & Resource Notes
- The stable saved frontier is in `data/frontier_summary.json` and covers `96`, `240`, and
  `496` ns.
- The representative deep-dive is in `data/study_summary.json` and
  `data/representative_traces.npz`; it carries the `240 ns` regime breakdown, reference
  slice, tradeoff slice, QND stress suite, robustness sweep, and convergence checks.
- Figure generation is cheap once the JSON and NPZ artifacts exist. The expensive part is
  the rich-model multistart optimization, not the plotting.
- The representative rich-model validation passed with `dt`-refinement fidelity change
  `1.43e-3` and truncation change `4.73e-5`.

## Resolved
- **Hardware replay activation was ambiguous during early testing**: this was resolved by
  switching to native `SequenceCompiler(hardware=...)` replay and by adding an explicit
  transport-analysis diagnostic.
- **The richer model initially failed to separate families by induced transition**: this was
  resolved by retuning the effective mixing onset and scales until the nominal and stress
  sweeps produced a nonzero, family-dependent measurement-induced transition channel.
- **Purely study-local effective mixing support**: resolved at the package level by adding
  reusable arbitrary-waveform strong-readout helpers, including higher-ladder continuation,
  to `cqed_sim.measurement.strong_readout`.
