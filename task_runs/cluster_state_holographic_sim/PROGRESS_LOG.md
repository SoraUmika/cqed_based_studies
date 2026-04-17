# Progress Log: cluster_state_holographic_sim

## 2026-03-24T00:00:00Z - Study initialized
- Created study folder structure: README.md, IMPROVEMENTS.md, study_state.json
- Problem class: DES + ANA
- Identified cqed_sim API: make_target('cluster', 1) builds SWAP·CZ·(H⊗I) per-site MPS isometry
- Key API modules: unitary_synthesis.targets, quantum_algorithms.holographic_sim, optimal_control

## Phase 1: Target Verification — COMPLETE
- Ran `make_target('cluster', 1)`: 4×4 unitary confirmed
- Unitarity error: ||U†U - I|| = 4.46e-16
- Matched canonical construction SWAP·CZ·(H⊗I)
- MPS tensors: A^0 = I/√2, A^1 = Z/√2
- Files: `scripts/run_study_v3.py` (Phase 1)

## Phase 2: Ideal Observables — COMPLETE
- N=6 cluster state via transfer-matrix formalism
- All single-site Pauli expectations: |⟨σ⟩| < 3e-32 (numerically zero)
- All stabiliser expectations: ⟨Ki⟩ = 1.0 (within machine precision)
- String-order correlator: +1.0 at nearest-neighbor separation (|j-i|=2), 0 otherwise
- Files: `scripts/run_study_v3.py` (Phase 2), `data/results.json`

## Phase 3: D-R-SNAP Decomposition — BLOCKED
- **CRITICAL FINDING**: SNAP gate in cqed_sim is cavity-only (I_q ⊗ diag(e^{iθ_n}))
- Does NOT apply different phases to |g,n⟩ vs |e,n⟩ — cannot create entanglement
- 2-SNAP synthesis: fidelity stuck at 0.500 (consistent with separable-only constraint)
- 3-SNAP synthesis: did not converge after >600 s
- Pivoted to GRAPE full-target optimization
- Files: `scripts/run_study_v3.py` (Phase 3)

## Phase 4: GRAPE Full-Target Optimization — COMPLETE
- Ran GRAPE (3 seeds × 300 iterations, N_cav=8, N_tr=2) at 6 durations
- Results:
  - 50 ns: F=0.6337 (75 s)
  - 100 ns: F=0.9494 (76 s)
  - 150 ns: F=0.9561 (100 s)
  - 200 ns: F=0.9966 (130 s)
  - 300 ns: F=0.9957 (214 s) — non-monotonic, GRAPE stochasticity
  - 400 ns: F=0.9990 (275 s)
- Files: `scripts/complete_study.py`, `scripts/run_study_v3.py`

## Holographic Channel Analysis — COMPLETE
- HolographicChannel.from_unitary: bond_dim=2
- Transfer matrix eigenvalues: {1, 0, 0, 0}
- Correlation length: infinite (characteristic of SPT-ordered cluster state)
- Files: `scripts/run_study_v3.py` (Phase 6), `data/results.json`

## Figures — COMPLETE
- Generated 7 figure types in PNG + PDF:
  stabilisers, string_order, grape_fidelity, pauli, string_order_decay, transfer_eigenvalues, grape_infidelity
- Files: `scripts/make_figures.py`, `figures/`

## Report — COMPLETE
- Wrote `report/report.tex` (6 pages, revtex4-2 two-column format)
- 15 references in `report/references.bib`
- Compiled via pdflatex → bibtex → pdflatex → pdflatex
- Includes mandatory appendix with full GRAPE data, observables, and transfer matrix analysis
- File: `report/report.pdf`

## IMPROVEMENTS.md — FINALIZED
- P1: SNAP cavity-only limitation, GRAPE stochasticity
- P2: Extended sweep, decoherence, SQR decomposition, Hilbert space convergence
- P3: Waveform simulation, sequential measurement protocol, GRAPE pulse visualization
- Documented failed approaches: D-R-SNAP, SQR+Powell, naive holographic state construction

## Study Status: COMPLETE

---

## Iteration 2: Multi-Decomposition Strategy Comparison

### Context Gathering — COMPLETE
- Read hybrid_qubit_cavity_control study: README.md, IMPROVEMENTS.md, phase2_idealized.py, common.py
- Read hybrid_unitary_native_entangling_evolution study: README.md, common.py
- Learned UnitarySynthesizer API: GateSequence pattern, gate constructor signatures
- Key discovery: SNAP has NO drift_model parameter; SQR, CP, FE all require drift_model

### API Discovery — COMPLETE
- Verified all gate constructor signatures via Python inspect.signature()
- Confirmed make_target returns ndarray directly (not dict)
- Identified correct GateSequence construction pattern from hybrid study scripts

### Multi-Decomposition Comparison — COMPLETE
- Wrote `scripts/decomposition_comparison.py`: Full comparison of 5 strategies (A-E)
- Wrote `scripts/decomposition_comparison_fast.py`: Streamlined version with precomputed results
- Wrote `scripts/decomposition_analysis.py`: Final analysis with hardcoded results, 5 figures, 3 artifacts
- Results (ideal mode, N_cav=2, chi=K=0):
  - A (D+R+SNAP): F=0.500 — SNAP cavity-only, confirmed limitation
  - B (D+SQR+CP, 1 block): F=0.707 — insufficient depth
  - B (D+SQR+CP, 2 blocks): F=1.000 — **PERFECT ideal decomposition** (9 gates, 117s)
  - B (D+SQR+CP, 3-4 blocks): F=1.000 — confirms convergence
  - C (SQR+CP, no D): F=0.500 — no Fock mixing without Displacement
  - D (D+R+FE, 2 blocks): F=0.9999 — native chi-wait, 8 gates
  - D (D+R+FE, 3-4 blocks): F=1.000 — converges with more depth

### Figure Generation — COMPLETE
- decomposition_comparison.{png,pdf} — bar chart of all strategies
- grape_fidelity_vs_duration.{png,pdf} — GRAPE sweep with infidelity
- strategy_ranking.{png,pdf} — ranked comparison
- depth_scaling_D_SQR_CP.{png,pdf} — scaling with circuit depth
- infidelity_comparison.{png,pdf} — log-scale infidelity

### Artifacts — COMPLETE
- artifacts/target_unitary.npz — 4×4 target unitary matrix
- artifacts/decomposition_best.json — best fidelity per gate-set family
- data/decomposition_comparison.json — full comparison data

### Report Extension — COMPLETE
- Backed up report.tex → report.tex.bak
- Inserted Extension section: Multi-Decomposition Strategy Comparison
- Added Reproducibility appendix (6 subsections: optimized parameters, waveforms, gate sequences, modeling assumptions, reproduction procedure, saved artifacts)
- Compiled PDF: 10 pages, 529285 bytes (up from 6 pages, 394471 bytes)

### State Files Updated — COMPLETE
- study_state.json: Added decomposition_comparison results, updated loop_iteration=2
- IMPROVEMENTS.md: Updated P1/P2 items, added compute notes, moved resolved items
- README.md: Updated problem class, goals, key results table, known limitations
- TASK_CHECKLIST.md: Added Phase 8 tasks, all checked complete

### AGENTS.md & Skills Updated — COMPLETE
- AGENTS.md: Added artifacts/ directory, reproducibility appendix template, reproducibility requirements section, updated decision tree
- .github/skills/latex-report/SKILL.md: Added Step 4 "Write Reproducibility Appendix (MANDATORY)"
- .github/skills/validate-results/SKILL.md: Added "Check 4 — Reproducibility Artifacts"

## Study Status: COMPLETE (Iteration 2)

## Iteration 3: Model-Based Verification & Decoherence Budget

### Science Director Phase
- Reviewed IMPROVEMENTS.md: 6 active P1/P2 items identified
- Classified: [P1|MEDIUM] ideal-mode verification, [P2|LOW] extended GRAPE, [P2|MEDIUM] decoherence, [P2|MEDIUM] FE wait-time, [P2|MEDIUM] hybrid strategy, [P2|LOW] Hilbert space convergence
- Write SCIENCE_DIRECTIVE.md with 5 numbered experiments

### Execution Phase
- Created iteration3_evaluate.py: evaluates ideal-mode parameters at N_cav=8,12,15
- **CRITICAL FINDING**: N_cav=2 ideal-mode results are truncation artifacts:
  - Strategy B (D+SQR+CP, 2 blocks): F=1.000 at N_cav=2 → F=0.094 at N_cav=8 (93% leakage)
  - Strategy D (D+R+FE, 2 blocks): F=0.9999 at N_cav=2 → F=0.188 at N_cav=8 (74% leakage)
  - Mechanism: Displacement gates populate n≥2 Fock states; SQR/CP have only n=0,1 parameters
- FE wait-time validation: FE0=FE1=177.0 ns = 1.005 × τ_CZ (π/|χ|=176.1 ns) — excellent match
- Hilbert space convergence: |F(N_cav=12) - F(N_cav=15)| = 0.0002 — converged
- GRAPE extended sweep (500,600,800 ns) already in results_combined.json — no new computation needed
- Coherence budget: T1=30μs, T2=20μs analytical estimate → GRAPE 200ns: F_comb=0.987 (best)
- Generated 3 new figures: hilbert_space_validity, hilbert_convergence, combined_fidelity_budget

### Validation Phase
- 6/6 validation checks passed (FE wait-time, leakage consistency, GRAPE truncation immunity, coherence monotonicity, N_cav convergence, GRAPE duration convergence)

### Report Phase
- Extended report.tex with new section: "Extension: Model-Based Verification and Decoherence Analysis"
- 5 subsections: truncation artifacts, FE validation, coherence budget, validation, revised conclusions
- 3 new figures, 2 new tables
- Compiled PDF: 12 pages, 649 KB

### State Update
- Updated: study_state.json (loop_iteration=3), IMPROVEMENTS.md (6 resolved items, 2 new P1/P2), TASK_CHECKLIST.md (Phase 10)
- Remaining P1: Re-optimise parametric strategies at N_cav=8 (HIGH difficulty)
- Remaining P2: Lindblad GRAPE, fine GRAPE resolution, warm-started GRAPE

## Study Status: COMPLETE (Iteration 3)
