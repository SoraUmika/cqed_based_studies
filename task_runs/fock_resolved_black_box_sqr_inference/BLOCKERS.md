# Blockers

## Active Blockers
- None.

## Resolved Blockers
- **2026-03-31 — `conditioned_multitone` full-mode helper on `n_tr = 3`**:
  scratch execution hit a dimension mismatch because the helper compares 2x2
  target density matrices against 3x3 reduced transmon states. This is not an
  active blocker because the study can still use the waveform-construction path
  and can evaluate pulse-level cases either on `n_tr = 2` models or with an
  explicit local truncation to the qubit subspace.
