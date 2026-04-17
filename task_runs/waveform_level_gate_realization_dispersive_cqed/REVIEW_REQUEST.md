# Review Request - Unconditional Cavity Displacement in Dispersive cQED

**Study:** `studies/waveform_level_gate_realization_dispersive_cqed`  
**Run:** `task_runs/waveform_level_gate_realization_dispersive_cqed`  
**Date:** 2026-04-02

## Status
The execution package is ready for critical review.

## What Changed In This Iteration
- Reframed the study around unconditional cavity displacement rather than the earlier mixed displacement-plus-qubit-rotation scope.
- Added a unified protocol benchmark covering:
  - naive single-tone
  - fast broadband single-tone
  - two-tone branch-compensated displacement
  - echoed displacement
  - hardware-aware optimal control
- Added a structured multiplex follow-up covering full-duration multicarrier fits, segmented branch-resonant fits, and a jointly optimized shaped two-tone family.
- Generated new `unconditional_*` artifacts and figures.
- Updated the report, summary files, and reproducibility notebook to use the new artifact chain and the corrected broad-state ranking.

## Key Scientific Claims To Review
1. The naive single-tone pulse fails as an unconditional displacement once `|chi| T` becomes appreciable, and this failure is severe on qubit superposition inputs.
2. `chi` is the dominant mechanism at the representative operating point; `chi'` and self-Kerr only weakly perturb the main conclusion.
3. Short two-tone branch compensation is the best simple, physically interpretable protocol in this study.
4. After the structured multiplex follow-up, the short `20 ns` two-tone pulse is also the best tested protocol overall on the explicit broad state-test set.
5. The bounded hardware-aware optimal-control pulse is best only in the narrower class of bounded sampled-waveform references.
6. The practical echoed protocol tested here is not competitive because the inserted qubit inversion is itself manifold dependent.
7. Full-duration multicarrier fits, segmented branch-resonant fits, and a low-parameter jointly optimized shaped-two-tone family all remain negative results on the broad state set.

## Suggested Review Focus
- Evidence-to-claim mapping in the unconditional-displacement report
- Whether the distinction between "best overall tested protocol" and "best sampled-waveform benchmark" is clearly and fairly justified
- Whether the report clearly separates vacuum branch matching from broader state-set performance
- Whether the limitations around decoherence, broader optimal-control search, coherent-state robustness, and the now-negative structured multiplex follow-up are stated strongly enough

## Main Files
- Report: `studies/waveform_level_gate_realization_dispersive_cqed/report/report.pdf`
- Summary: `task_runs/waveform_level_gate_realization_dispersive_cqed/EXECUTION_SUMMARY.md`
- Science directive: `task_runs/waveform_level_gate_realization_dispersive_cqed/SCIENCE_DIRECTIVE.md`
- Notebook: `studies/waveform_level_gate_realization_dispersive_cqed/scripts/reproducibility_notebook.ipynb`

## Review Invocation
`study=studies/waveform_level_gate_realization_dispersive_cqed run=task_runs/waveform_level_gate_realization_dispersive_cqed phase=review`
