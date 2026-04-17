# Task Checklist

## Status Summary
- Study: `studies/fast_robust_storage_vacuum_reset_comparison`
- Run: `task_runs/fast_robust_storage_vacuum_reset_comparison`
- Problem class: `DES`, `ANA`, `OPT`
- Current state: iteration 3 repair pass implemented; metric contradictions are fixed, but notebook execution has not been re-verified and the study remains `ACTIVE` pending a resubmission decision

## Initialize And Plan
- [x] Read AGENTS instructions and repo conventions
- [x] Create the study folder, state file, README, and improvement log
- [x] Write the study-scoped science directive
- [x] Record the first-principles model and `cqed_sim` gap analysis

## Implement
- [x] Implement the comparative architecture runner
- [x] Reuse the validated pulsed ladder and two-tone sideband control settings
- [x] Add mixed-state, coherent-state, and higher-Fock initial-state tests
- [x] Run robustness sweeps over decoherence, detuning, amplitude, reset error, thermal loading, and readout linewidth
- [x] Export machine-readable artifacts and report-quality figures

## Validate
- [x] Demonstrate timestep and truncation convergence on representative cases
- [x] Quantify speed-versus-robustness tradeoffs across schemes
- [ ] Demonstrate analytic-versus-numeric agreement for the reduced effective models
- [ ] Check whether the conclusions remain stable under the matched transmon-reference noise model

## Report And Reproducibility
- [x] Update the README validation section and improvement log
- [x] Write `report/report.tex` (iteration 2: full manuscript rebuild)
- [x] Compile `report/report.pdf`
- [x] Create `scripts/reproducibility_notebook.ipynb`
- [x] Update `EXECUTION_SUMMARY.md` and `REVIEW_REQUEST.md`

## Iteration 2 Follow-up Actions (from REVIEW_DIRECTIVE)
- [x] [REPORT] Rebuild report into complete scientific manuscript
- [x] [VALIDATE] Surface convergence data quantitatively
- [x] [UNCERTAINTY] Add sensitivity/uncertainty treatment
- [x] [FIGURE] Fix broken summary figure and integrate omitted evidence
- [x] [CLAIMS] Repair or qualify unsupported claims
- [x] [REPRODUCIBILITY] Expand reproducibility appendix
- [x] [STATE] Reconcile task-run bookkeeping

## Iteration 3 Follow-up Actions (from REVIEW_DIRECTIVE)
- [x] [METRICS] Define canonical metric dictionary and crosswalk across report numbers and saved artifacts
- [x] [ARTIFACTS] Reconcile or regenerate scheme summary, initial-state, convergence, and study-results artifacts
- [x] [REPORT] Rewrite contradictory abstract, results, validation, discussion, and conclusion claims after reconciliation
- [x] [VALIDATE] Correct the dt-stability narrative so it matches the verified convergence values
- [x] [STATE] Align README, TASK_CHECKLIST, study_state, and notebook verification status with actual completion state
