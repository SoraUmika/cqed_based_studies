# Task Checklist

## Status Summary
- Study: studies/fock_resolved_black_box_sqr_inference
- Run:   task_runs/fock_resolved_black_box_sqr_inference
- Loop iteration: 0

## Bootstrap
- [x] B0.1 Initialize study directory and state files
- [x] B0.2 Write the first study-specific SCIENCE_DIRECTIVE.md
- [x] B0.3 Replace boilerplate README and IMPROVEMENTS scaffolding with study-specific content

## Phase 1: Planning
- [x] P1.1 Review `AGENTS.md`, the local `cqed_sim` API reference, and relevant prior studies
- [x] P1.2 Perform an analytic identifiability check before writing simulation code
- [x] P1.3 Classify the protocol gaps between existing `cqed_sim` functionality and the requested inverse problem
- [x] P1.4 Finalize the experiment matrix and success criteria in `study_state.json`

## Phase 2: Implementation
- [x] I2.1 Implement runtime/path compatibility helpers for the study scripts
- [x] I2.2 Implement the single-qubit forward simulator and Cholesky MLE baseline
- [x] I2.3 Implement the Fock-diagonal forward kernels for wait-only, displacement-only, and combined protocols
- [x] I2.4 Implement constrained least squares and binomial MLE on the recoverable transverse-sector parameters
- [x] I2.5 Implement the diagnostic full-state Cholesky-plus-softmax MLE to expose non-identifiability
- [x] I2.6 Implement black-box case generation: ideal diagonal, coherent, pulse-level near-ideal, imperfect/leaky, and noisy cases
- [x] I2.7 Run the main study and serialize data / artifact summaries
- [x] I2.8 Generate all required figures in both PNG and PDF formats

## Phase 3: Validation
- [x] V3.1 Sanity-check the analytic invariance claims against direct numerical simulation
- [x] V3.2 Verify convergence with respect to time-grid density, alpha-grid density, and cavity truncation
- [x] V3.3 Verify that the pulse-level case ranking is stable under a modest transmon-truncation check
- [x] V3.4 Compare least squares versus binomial MLE on the recoverable transverse parameters
- [x] V3.5 Quantify full-state MLE non-uniqueness across random restarts on a representative diagonal case
- [x] V3.6 Review and finalize `IMPROVEMENTS.md`

## Phase 4: Reporting
- [x] R4.1 Write `EXECUTION_SUMMARY.md`
- [x] R4.2 Write `report/report.tex` with main text, validation section, limitations, appendix, and reproducibility section
- [x] R4.3 Compile `report/report.pdf` and confirm no layout issues
- [x] R4.4 Create `scripts/reproducibility_notebook.ipynb`
- [x] R4.5 Write `REVIEW_REQUEST.md`
