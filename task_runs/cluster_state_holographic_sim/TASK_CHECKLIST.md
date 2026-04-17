# Task Checklist: cluster_state_holographic_sim

## Phase 1: Target Verification
- [x] T1.1: Construct canonical cluster-state transfer matrix
- [x] T1.2: Compare with make_target('cluster', 1) — compute fidelity
- [x] T1.3: Extract MPS tensors from both conventions and verify equivalence
- [x] T1.4: Generate N-site cluster state via holographic protocol and verify

## Phase 2: Ideal Observables
- [x] T2.1: Compute single-site Pauli expectations ⟨Xi⟩, ⟨Yi⟩, ⟨Zi⟩
- [x] T2.2: Compute stabilizer expectations ⟨Ki⟩
- [x] T2.3: Compute ZXZ correlators
- [x] T2.4: Compute string-order correlators
- [x] T2.5: Tabulate and verify against analytic predictions

## Phase 3: Decomposition
- [x] T3.1: Decompose target into QubitRotation + Displacement + SNAP — **BLOCKED: SNAP is cavity-only, cannot entangle qubit-cavity**
- [x] T3.2: Report decomposition fidelity — **F=0.5 (SNAP limitation)**
- [ ] T3.3: Tabulate SNAP gate parameters — **N/A: decomposition not valid**

## Phase 4: SNAP Pulse Optimization → GRAPE Full-Target Optimization
- [x] T4.1: GRAPE-optimize full target unitary at multiple durations (50–400 ns)
- [x] T4.2: Report optimized fidelity vs duration: best F=0.9990 at 400 ns
- [x] T4.3: Check for leakage at short durations — all zero (ideal mode)

## Phase 5: Timing Analysis
- [x] T5.1: GRAPE-dominated timing summary (single-unitary approach)
- [ ] T5.2: Full sequence timing with 16 ns rotations, 48 ns displacements — **N/A: D-R-SNAP decomposition not valid**

## Phase 6: Waveform-Level Simulation
- [ ] T6.1: Build pulse-level sequence — **Deferred (P3): requires GRAPE pulse extraction**
- [ ] T6.2: Compare decomposition vs pulse fidelity — **Deferred**
- [ ] T6.3: Compare observables at pulse level — **Deferred**

## Phase 7: Validation & Report
- [x] T7.1: Sanity checks (limiting cases: unitarity, stabilizer values, transfer matrix eigenvalues)
- [x] T7.2: Convergence analysis (documented GRAPE seed/iteration sensitivity)
- [x] T7.3: Write report.tex and compile PDF
- [x] T7.4: Finalize IMPROVEMENTS.md

## Phase 8: Multi-Decomposition Comparison (Iteration 2)
- [x] T8.1: Gather context from hybrid_qubit_cavity_control and hybrid_unitary_native_entangling_evolution
- [x] T8.2: Discover and validate cqed_sim UnitarySynthesizer API signatures
- [x] T8.3: Write decomposition_comparison.py (5 strategies, A-E)
- [x] T8.4: Run Strategy A (D+R+SNAP): F=0.500 confirmed
- [x] T8.5: Run Strategy B (D+SQR+CP, 1-4 blocks): F=1.000 at 2 blocks
- [x] T8.6: Run Strategy C (SQR+CP, no D): F=0.500 confirmed
- [x] T8.7: Run Strategy D (D+R+FE, 2-4 blocks): F=0.9999 at 2 blocks
- [x] T8.8: Generate 5 decomposition comparison figures (PNG+PDF)
- [x] T8.9: Save artifacts (target_unitary.npz, decomposition_best.json, decomposition_comparison.json)
- [x] T8.10: Extend report.tex with decomposition comparison + reproducibility appendix
- [x] T8.11: Compile extended PDF (10 pages, 529285 bytes)
- [x] T8.12: Update IMPROVEMENTS.md, README.md, study_state.json

## Phase 9: AGENTS.md & Skills Update (Iteration 2)
- [x] T9.1: Add artifacts/ directory to AGENTS.md study structure
- [x] T9.2: Add Reproducibility appendix template to AGENTS.md LaTeX template
- [x] T9.3: Add Reproducibility Requirements section to AGENTS.md
- [x] T9.4: Update "Is the study complete?" decision tree in AGENTS.md
- [x] T9.5: Update latex-report SKILL.md with mandatory reproducibility step
- [x] T9.6: Update validate-results SKILL.md with artifacts check

## Phase 10: Model-Based Verification & Decoherence (Iteration 3)
- [x] T10.1: Write SCIENCE_DIRECTIVE.md for Iteration 3
- [x] T10.2: Create iteration3_evaluate.py (model-based evaluation at N_cav=8,12,15)
- [x] T10.3: Run evaluation — CRITICAL FINDING: N_cav=2 results are truncation artifacts (F=1.000→0.094)
- [x] T10.4: FE wait-time analysis: FE0=FE1=177.0 ns = 1.005 × τ_CZ — validated
- [x] T10.5: Create and run iteration3_coherence_budget.py (T1=30μs, T2=20μs)
- [x] T10.6: Coherence budget result: GRAPE 200ns achieves best combined F=0.987
- [x] T10.7: Verify GRAPE extended sweep data (500,600,800 ns already in results_combined.json)
- [x] T10.8: Create and run iteration3_figures.py (3 new figures)
- [x] T10.9: Run iteration3_validation.py (6/6 checks passed)
- [x] T10.10: Extend report.tex with Iteration 3 section
- [x] T10.11: Compile PDF (12 pages, 649 KB)
- [x] T10.12: Update IMPROVEMENTS.md, study_state.json, TASK_CHECKLIST.md
