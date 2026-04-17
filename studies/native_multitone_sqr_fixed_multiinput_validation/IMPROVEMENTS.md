# Improvement Log: Native Multitone SQR Fixed Multi-Input Validation

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- **[P1 | MEDIUM] Full historical re-optimization not repeated**: this study re-simulates representative saved native multitone pulses under the fixed package, but it does not re-run every original optimization grid. A future extension should re-optimize the older studies end-to-end with the corrected `fock_fqs_hz` handling and compare the resulting parameter shifts.

## Recommended Improvements (P2)
- **[P2 | MEDIUM] Expand beyond representative cases**: add a larger sample of hard/random cases from the arbitrary-target studies to determine whether the multi-input conclusions are representative of the full saved ensembles.
- **[P2 | LOW] Add reduced effective-unitary extraction everywhere**: the quartet probe set is strong, but a uniform reduced process-tomography API across all studies would simplify direct comparisons.

## Nice-to-Haves (P3)
- **[P3 | LOW] Add cavity-superposition probes**: the present study focuses on single-manifold qubit inputs because that is the user's requested reduced objective. A supplementary pass on cavity-superposition states could show where reduced and full logical-subspace validation begin to diverge.

## Open Questions
- Do the older native multitone pulses fail the fixed multi-input test mainly because of the now-incorrect `fock_fqs_hz` override, or because the target families themselves were too ambitious for the chosen pulse durations?
- How much of the gap between single-state and quartet validation is genuine control limitation versus hidden phase freedom left unconstrained by the earlier objectives?

## What Was Tried and Did Not Work
- Initial direct shell-to-Python JSON inspection using piped here-strings on Windows did not pass stdin to `python -`; switched to `python -c` for reliable probing of saved artifacts.

## Compute & Resource Notes
- The four representative fixed-package reevaluations plus figure/report/notebook generation completed in a few seconds on the current workstation; the actual pulse propagation was not the runtime bottleneck.
- The report compiled successfully with `pdflatex`; `bibtex` reports that there are no citation commands, which is expected for this internal validation note.
