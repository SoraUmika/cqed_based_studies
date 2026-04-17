# Task Checklist

## Status Summary
- Study: `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed`
- Run: `task_runs/multitone_sqr_no_go_decoupled_echo_dispersive_cqed`
- Problem class: `ANA`, `REP`, `DES`
- Current state: science-director-approved report polished; updated LaTeX source verified via `report_build.pdf`; `report.pdf` remains locked by another process and could not be refreshed during this pass
- Next phase: `complete`
- Accepted baseline: archived report, figures, data, notebook, and validation package remain in force unless targeted extension regressions contradict them

## Extension Pass
- [x] Re-read the study context and archived iteration handoff files
- [x] Write the extension-scoped science directive
- [x] Map the strict two-block tuned cancellation set and select exact-simulation checkpoints
- [x] Run focused shared-line exact checks on tuned and nearby off-tuned points
- [x] Test symmetry-aligned and, if supported, manifold-aware echoed refocusing against the archived echo baseline
- [x] Validate the extension-specific results and prepare the preserved report inputs for an appended extension section

## Initialize And Plan
- [x] Read AGENTS instructions and repository conventions
- [x] Audit nearby SQR studies and identify scope mismatches
- [x] Initialize the study folder, README, and improvement log
- [x] Write the study-scoped science directive

## Implement
- [x] Build strict no-detuning shared-line multitone helpers on top of `cqed_sim`
- [x] Implement exact logical-subspace diagnostics and blockwise `X/Y/Z` decomposition
- [x] Implement the stronger decoupled-block reduced model
- [x] Implement the echoed `half-SQR -> pi -> half-SQR -> pi` replay and analysis
- [x] Generate machine-readable artifacts and report-quality figures

## Validate
- [x] Confirm the analytical two-block and many-block no-go arguments against numerical trends
- [x] Check convergence with respect to truncation, timestep, and optimization budget
- [x] Compare representative conclusions against prior repository claims and relevant literature where needed

## Report And Reproducibility
- [x] Update the README and improvement log with final findings
- [x] Write `report/report.tex`
- [x] Verify an updated report build (`report_build.pdf`; in-place overwrite of `report.pdf` was blocked by a Windows file lock)
- [x] Create `scripts/reproducibility_notebook.ipynb`
- [x] Write `EXECUTION_SUMMARY.md`
- [x] Write `REVIEW_REQUEST.md`

## Review Handoff
- [x] Prepare the study package for science-director review

## Polish
- [x] Re-read the APPROVE directive and current run/state files before editing
- [x] Polish the report flow, method scanability, and reproducibility appendix without reopening the scientific scope
- [x] Verify the polished PDF build via `report_build.pdf`
- [x] Re-test the canonical `report.pdf` refresh and record the persistent Windows file lock
- [x] Write `POLISH_COMPLETE.md` and advance the study/run metadata to `COMPLETE`
