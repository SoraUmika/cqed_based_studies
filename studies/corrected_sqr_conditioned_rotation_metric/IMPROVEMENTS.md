# Improvement Log: Corrected SQR Optimization with a Fock-Resolved Effective-Qubit Metric

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- **[P1 | MEDIUM] Public reduced-unitary API missing**: the convention fix aligned `conditioned_multitone` phase handling with `R_phi(theta)`, but there is still no public reduced multitone effective-unitary extractor. This study therefore kept a local helper.

## Recommended Improvements (P2)
- **[P2 | MEDIUM] Expand target family coverage**: only one smooth four-level corrected-SQR profile and its prefix subwindows were scanned here. Add sparse and random corrected-SQR targets.
- **[P2 | MEDIUM] Add open-system replay**: the current conclusions are closed-system only.

## Nice-to-Haves (P3)
- **[P3 | LOW] Add public reporting of axis-z contamination**: the reduced unitary rerun found residual extracted axis-z tilt to be a useful diagnostic.

## Open Questions
- How much further can segmented or basis-modulated multitone waveforms reduce the residual axis-z contamination without needing an echoed protocol?
- Does a richer duration sweep reveal a sharper optimum than the present `|chi|T/2pi = 1, 3, 5` grid?

## What Was Tried and Did Not Work
- Optimizing only the reduced final-state fidelity remained insufficient even after the phase-convention fix. The resulting pulses routinely overestimated success relative to the stricter reduced effective-unitary metric.

## Compute & Resource Notes
- Main rerun: `scripts/run_study.py`.
- Optimization grid: 12 cases (`N_active = 1..4`, `|chi|T/2pi = 1,3,5`), each with baseline, analytic warm start, state-based optimization, and unitary-based optimization.
- Reduced objective stayed cheap because every evaluation used only sector-by-sector two-level propagators.
