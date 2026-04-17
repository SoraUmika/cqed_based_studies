# Task Checklist

## Status Summary
- Study: studies/microwave_component_thermalization_cqed
- Run: task_runs/microwave_component_thermalization_cqed
- Loop iteration: 0
- Execution status: READY_FOR_REVIEW

## Bootstrap
- [x] B0.1 Initialize study directory and state files
- [x] B0.2 Write initial SCIENCE_DIRECTIVE.md

## Phase 1: Planning
- [x] P1.1 Classify the problem and write the analytic preliminary
- [x] P1.2 Audit the `cqed_sim` API and document gaps
- [x] P1.3 Define observables, validation gates, and compute strategy

## Phase 2: Implementation
- [x] I2.1 Implement shared `cqed_sim` thermal study utilities
- [x] I2.2 Run thrust A thermometer simulations
- [x] I2.3 Run thrust B coherence and dephasing simulations
- [x] I2.4 Run thrust C multimode auxiliary-mode maps
- [x] I2.5 Run thrust D transient temperature-step simulations

## Phase 3: Validation
- [x] V3.1 Zero-temperature and weak-coupling sanity checks
- [x] V3.2 Hilbert-space truncation checks
- [x] V3.3 Analytic thermal-scaling comparisons

## Phase 4: Reporting
- [x] R4.1 Generate figures and machine-readable artifacts
- [x] R4.2 Write markdown research memo
- [x] R4.3 Write LaTeX report and compile PDF
- [x] R4.4 Create reproducibility notebook
- [x] R4.5 Write execution summary and review request
