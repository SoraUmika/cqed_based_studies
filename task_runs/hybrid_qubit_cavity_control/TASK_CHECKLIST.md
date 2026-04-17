# Task Checklist

## Status Summary
- Study: studies/hybrid_qubit_cavity_control
- Run: task_runs/hybrid_qubit_cavity_control
- Loop iteration: 14

## Bootstrap
- [x] B0.1 Resume existing study and create task-run state files
- [x] B0.2 Audit existing hybrid and SQR study context

## Phase 1: Planning
- [x] P1.1 Review cqed_sim selective-control support and waveform-bridge limits
- [x] P1.2 Define the speed-limit / feasibility extension scope
- [x] P1.3 Record the extended plan in study outputs
- [x] P1.4 Identify the missing second-pass gaps: `B_local` fairness, sparse-variant selection, optimizer refinement, and qutrit replay

## Phase 2: Implementation
- [x] I2.1 Add speed-limit and physical-feasibility analysis scripts
- [x] I2.2 Generate ideal decomposition frontier outputs with t_R=16 ns and t_D=48 ns
- [x] I2.3 Generate SQR/CPSQR/SNAP duration and sparse-tone frontier outputs
- [x] I2.4 Generate strategy-ranking figures and summary tables
- [x] I2.5 Add fair `B_local` physical compilation using the optimized gate-set archive
- [x] I2.6 Add adaptive sparse-variant selection across `full`, `top4`, and `top2`
- [x] I2.7 Run fixed-duration optimizer-backed refinement for the selected sparse variants
- [x] I2.8 Replay refined waveforms on an `n_tr=3` model and record qutrit leakage diagnostics

## Phase 3: Validation
- [x] V3.1 Run sanity checks on decomposition formulas and gate realizability outputs
- [x] V3.2 Run convergence / stability spot checks for the new physical study
- [x] V3.3 Update validation summary artifacts
- [x] V3.4 Validate refined-entangler dt stability and qutrit leakage bounds

## Phase 4: Reporting
- [x] R4.1 Update README and IMPROVEMENTS with the new extension results
- [x] R4.2 Revise report.tex with separated speed-limit sections and practical recommendations
- [x] R4.3 Compile report.pdf
- [x] R4.4 Refresh study_state and task-run handoff files for the second-pass conclusion

## Phase 5: Extension Pass (Improvement Studies)
- [x] E5.1 Create `run_extension_pass.py` implementing P2.1, P2.2, P2.3, P1.3, P1.1
- [x] E5.2 P2.1: SNAP-equivalent pi-SQR unification (F=0.5321, matches B_ent_S)
- [x] E5.3 P2.2: Stage-2 refinement diagnostics (all 9 gates ABNORMAL, flat landscape)
- [x] E5.4 P2.3: B_local_S1 bottleneck isolation (complexity ratio 8.8×, saturated)
- [x] E5.5 P1.3: Joint duration-plus-correction optimization (L1d_S1 +0.34)
- [x] E5.6 P1.1: Full compiled sequence replay (B_ent=0.691, B_local=0.338, L1d=0.004)
- [x] E5.7 Fix Subspace API bug (from_standard_block → qubit_cavity_block)
- [x] E5.8 Fix fidelity computation bug (pass target_subspace to simulate_sequence)
- [x] E5.9 Write unified extended report (23 pages, clean PDF)
- [x] E5.10 Update IMPROVEMENTS.md with resolved items and new findings
- [x] E5.11 Update README.md, study_state.json, and task-run documents

## Phase 6: Extension Pass 2 (GRAPE and Comprehensive Analysis)
- [x] E6.1 Create `run_extension_pass_2.py` implementing GRAPE, fine grid, hardware, SU(4)
- [x] E6.2 GRAPE optimization for 3 bottleneck gates (B_local_S1→0.9910, B_ent_S→0.9919, L1d_S1→0.9888)
- [x] E6.3 Fine 20-point duration grid for L1d_S1 and B_ent_S
- [x] E6.4 Hardware constraint analysis (9/9 pass)
- [x] E6.5 Random SU(4) test (partially resolved; model-based blocked by framework error)
- [x] E6.6 Fix GRAPE channel specs and target matrix ordering bugs
- [x] E6.7 Fix TargetUnitary, UnitarySynthesizer, and gate constructor API bugs
- [x] E6.8 Update report.tex with Extension Pass 2 section (29 pages)
- [x] E6.9 Fix LaTeX math-mode error and recompile PDF (clean, no errors)
- [x] E6.10 Update IMPROVEMENTS.md with GRAPE resolved items and post-GRAPE priorities
- [x] E6.11 Update study_state.json, PROGRESS_LOG, EXECUTION_SUMMARY, TASK_CHECKLIST

## Phase 7: Upstream Fixes and Study Extension (Iterations 5-14)
- [x] U7.1 Upstream fix: `bloch_xyz_from_qubit_state` accepts n_tr > 2 via `truncate_to_qubit_subspace()` (extractors.py)
- [x] U7.2 Upstream fix: ConditionalPhaseSQR waveform bridge (io/gates.py + waveform_bridge.py)
- [x] U7.3 Investigate `_expand_target_matrix` — confirmed NOT a bug; model-based UnitarySynthesizer works correctly (objective=0.530195)
- [x] U7.4 Re-verify all prior study results under updated cqed_sim (24/24 checks passed)
- [x] U7.5 Upstream fix: `_conditioned_metric` truncates rho_q for n_tr>2 (targeted_subspace_multitone.py). n_tr=3 inner-loop study: 9/9 validation + 3/3 re-optimization complete. Max |f> leak < 8.1e-7.
- [x] U7.6 L2d/L2c CPSQR full-sequence waveform replay (enabled by U7.2). L2d: ideal=0.9441, pulse=0.2001, leak_avg=0.280. L2c: ideal=0.9558, pulse=0.1319, leak_avg=0.275. Full five-strategy comparison complete.
- [x] U7.7 Full-sequence GRAPE for B_local (F_nom=0.993, F_strict=0.775, leak=0.239) and L2d (F_nom=0.774, F_strict=0.361, leak=0.245). Proves B_local unitary achievable at 2440 ns. Strict-replay gap indicates simulation-GRAPE disagreement, not physics limit.
- [x] U7.8 Iter 12: Finer-dt (4 ns) sequence GRAPE — dt hypothesis DISPROVEN. B_local: F_nom=0.999999, F_strict=0.703 (WORSE). L2d: F_nom=0.999969, F_strict=0.754. Propagator model mismatch confirmed.
- [x] U7.9 Iter 13: 5 cross-validation tests. V3 PASS (n_cav=8 converged, delta=0). V2 confirms gap accumulates with sequence length (single-gate gap=0.058 vs multi-gate gap=0.218-0.297). V1 shows identity is non-trivial in rotating frame. All data files intact.
- [x] U7.10 Iter 14: Final report extension with Iterations 5-13 content (upstream fixes, CPSQR replay, sequence GRAPE, propagator analysis, validation). 25-page PDF compiled with 5 new main-text sections + 1 appendix section. All figures referenced.