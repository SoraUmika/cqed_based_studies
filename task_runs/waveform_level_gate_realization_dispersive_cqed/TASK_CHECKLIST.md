# Task Checklist

## Status Summary
- Study: `studies/waveform_level_gate_realization_dispersive_cqed`
- Run: `task_runs/waveform_level_gate_realization_dispersive_cqed`
- Problem class: `ANA`, `DES`, `OPT`
- Current state: study complete; archival package updated through the structured multiplex follow-up

## Initialize And Plan
- [x] Read AGENTS instructions and repo conventions
- [x] Reuse the existing `cqed_sim` stack rather than creating ad hoc simulation code
- [x] Write the study-scoped science directive and execution plan
- [x] Confirm required deliverables: report, artifacts, figures, notebook, summary files

## Implement
- [x] Extend shared helpers in `scripts/common.py` for shaped displacement, branch-frequency extraction, and distance metrics
- [x] Implement the unconditional-displacement driver covering naive, fast, two-tone, echoed, and optimal-control families
- [x] Generate machine-readable artifacts for each protocol family
- [x] Generate report-quality figures for branch mismatch, entanglement, filtering, `chi` scaling, protocol comparison, and Wigner diagnostics

## Validate
- [x] Confirm the minimal analytic picture matches the simulated naive-pulse failure
- [x] Compare minimal, higher-order, and full Hamiltonian variants to isolate the dominant error mechanism
- [x] Reuse the previously validated propagator settings for the same model family and maintain the same default truncation and time step
- [x] Compare the observed regime boundary against the expected `1 / |chi|` timescale

## Report And Reproducibility
- [x] Update the study README and improvement log
- [x] Rewrite `report.tex` around the unconditional-displacement study
- [x] Compile `report.pdf`
- [x] Regenerate `scripts/reproducibility_notebook.ipynb` against the `unconditional_*` artifacts
- [x] Update the execution summary and review request files
- [x] Add the multiplex-drive interpretation and follow-up guidance to the report and study docs
- [x] Run and archive an explicit multiplex-drive follow-up benchmark
- [x] Test segmented and jointly optimized structured multiplex refinements and fold the final negative result into the archival package

## Review Handoff
- [x] Prepare the study package for science-director review
