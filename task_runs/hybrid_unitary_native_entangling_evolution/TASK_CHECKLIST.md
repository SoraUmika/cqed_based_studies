# Task Checklist

## Status Summary
- Study: studies/hybrid_unitary_native_entangling_evolution
- Run: task_runs/hybrid_unitary_native_entangling_evolution
- Loop iteration: 2

## Bootstrap
- [x] B0.1 Initialize study directory and core files
- [x] B0.2 Audit the earlier hybrid-unitary study and candidate archive
- [x] B0.3 Document the new entangling-weighted study scope

## Phase 1: Planning
- [x] P1.1 Confirm cqed_sim support for synthesis, replay, Bloch, and Wigner diagnostics
- [x] P1.2 Define the entangling-aware cost model and candidate ranking strategy
- [x] P1.3 Bootstrap the archive-derived candidate frontier

## Phase 2: Implementation
- [x] I2.1 Add a native-block composition search for full U_target candidates
- [x] I2.2 Generate Phase 2 figures and machine-readable outputs
- [x] I2.3 Record the best native-heavy candidates in study_state.json and EXECUTION_SUMMARY.md

## Phase 3: Validation
- [x] V3.1 Re-verify replay support for ConditionalPhaseSQR and FreeEvolveCondPhase in the installed cqed_sim build
- [x] V3.2 Run depth-resolved X/Y/Z diagnostics for shortlisted candidates
- [x] V3.3 Run Wigner-versus-depth diagnostics for shortlisted candidates

## Phase 5: Runtime Validation
- [x] V5.1 Replace the exact qubit Hadamard with a replayable local decomposition
- [x] V5.2 Build replayable or surrogate-backed runtime candidates for the required shortlist
- [x] V5.3 Generate runtime artifacts, symbolic-versus-runtime comparison tables, and truncation sweeps at n_cav = 10, 12, 14
- [x] V5.4 Run nominal-noise replay for the runtime finalists at n_cav = 12, n_tr = 3
- [x] V5.5 Save Bloch and Wigner diagnostics for the replay-backed finalists

## Phase 4: Reporting
- [x] R4.1 Draft report/report.tex with the candidate table and Phase 2 conclusions
- [x] R4.2 Compile report.pdf after validation
- [x] R4.3 Revise the report to reflect the runtime-validated conclusion and saved artifacts