# Blockers: Native / Rich Multitone Feasibility for Ideal SQR and CPSQR

## Active Blockers
- None.

## Resolved
- **Resolved on 2026-03-27: stale `fock_fqs_hz` helper mismatch in earlier study-local wrappers**
  - Observed: multiple older study `make_run_config(...)` helpers populated `fock_fqs_hz` from `manifold_transition_frequency(..., frame=frame) / 2pi`, i.e. already frame-shifted frequencies.
  - Expected: patched `cqed_sim.pulses.calibration.build_sqr_tone_specs(...)` interprets `fock_fqs_hz` as absolute qubit transition frequencies in Hz and subtracts the frame internally.
  - Specific call path: `ConditionedMultitoneRunConfig(..., fock_fqs_hz=<frame-relative list>)` -> `build_conditioned_multitone_tones(...)` -> `build_sqr_tone_specs(...)`.
  - Suspected cause: earlier local helper code predated the patched absolute-frequency interpretation and preserved the old frame-relative override pattern.
  - Resolution used in this study: do not reuse those wrappers; the new study constructs corrected run configs with `fock_fqs_hz=None` unless an explicit absolute-frequency override is intentionally needed and audited.
