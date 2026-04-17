# Cluster-State Holographic Sequential Simulation

## Problem Class
DES / ANA / OPT — Experiment Design + System Analysis + Parameter Optimization

## Motivation
Holographic quantum simulation uses a single qubit–cavity system to sequentially prepare and measure multi-qubit entangled states by exploiting the MPS (matrix product state) structure. The 1D cluster state is a canonical resource state for measurement-based quantum computation. This study verifies the per-site transfer unitary implemented in `cqed_sim`, compares multiple decomposition strategies into native bosonic primitives, optimizes via GRAPE, and reports the expected cluster-state observables and total gate times.

## Goals
1. Verify that `make_target('cluster', n_match=1)` correctly produces MPS tensors generating the canonical 1D cluster state.
2. Compute ideal cluster-state observables: single-site Pauli expectations, stabilizer expectations ⟨K_i⟩, ZXZ correlators, and string-order correlators.
3. Decompose the verified target into the cqed_sim gate set (QubitRotation + Displacement + SNAP).
4. Pulse-optimize each SNAP gate for minimum duration at a stated fidelity threshold.
5. Report SNAP-dominated and full-sequence timing summaries (including 16 ns qubit rotations, 48 ns displacements).
6. **[Extension] Multi-decomposition comparison**: Systematically compare 5 strategies (D+R+SNAP, D+SQR+CP, SQR+CP-only, D+R+FreeEvolveCondPhase, GRAPE) across multiple circuit depths, informed by hybrid_qubit_cavity_control study methodology.
7. Explain how the final protocol supports holographic sequential measurement of cluster-state correlations.

## Methods
- `cqed_sim.unitary_synthesis.targets.make_target('cluster', 1)` — per-site transfer unitary
- `cqed_sim.quantum_algorithms.holographic_sim` — HolographicChannel, HolographicSampler
- `cqed_sim.unitary_synthesis.UnitarySynthesizer` — gate-set decomposition (Strategies A-D)
- `cqed_sim.unitary_synthesis.{Displacement, SQR, ConditionalPhaseSQR, FreeEvolveCondPhase, SNAP, QubitRotation}` — gate primitives
- `cqed_sim.optimal_control.GrapeSolver` — GRAPE optimization (Strategy E)
- `cqed_sim.core.DispersiveTransmonCavityModel` — physical Hamiltonian model

## Expected Outcomes
- State fidelity F ≈ 1.0 between implemented and canonical cluster state.
- All stabilizer expectations ⟨K_i⟩ ≈ +1 for the ideal target.
- Single-site ⟨X⟩, ⟨Y⟩, ⟨Z⟩ ≈ 0.
- Decomposition fidelity > 0.99 for qubit-conditional gate sets (SQR, FE).
- Identification of which decomposition strategies are practical, efficient, and physically implementable.

## Key Results (through Iteration 4)
| Strategy | Fidelity | Gates | Key Finding |
|----------|----------|-------|-------------|
| A: D+R+SNAP | 0.500 | 7 | SNAP cavity-only, no entanglement |
| B: D+SQR+CP (2 blocks) | 1.000 | 9 | **Perfect ideal decomposition** |
| C: SQR+CP (no D) | 0.500 | 6 | No Fock mixing without Displacement |
| D: D+R+FE (2 blocks) | 0.9999 | 8 | Native chi-wait, hardware-friendly |
| E: GRAPE 400ns | 0.999 | — | Model-based optimal control |

- **Iteration 3:** the ideal-mode Strategy B and D results were shown to be truncation artifacts. At enlarged cavity dimension the ideal decompositions collapse, and the full-drift $N_\mathrm{cav}=12$ check for Strategy B falls to $\mathcal{F}=0.0747$.
- **Iteration 4:** bounded-displacement optimisation at enlarged cavity dimension does not rescue the native decompositions. The best bounded embedded result is Strategy D with $|\alpha|\leq 0.3$, giving $\mathcal{F}(N_\mathrm{cav}=12)=0.5465$ with $22.8\%$ average leakage; the lowest-leakage Strategy B result reaches $\mathcal{F}(N_\mathrm{cav}=12)=0.5291$ with $4.0\%$ leakage.
- **Active-tone minimisation:** the bounded Strategy B sequence can be compressed to 20 ns per gate (180 ns total active time) with only a tiny embedded-fidelity change, from $0.5291$ to $0.5305$ at $N_\mathrm{cav}=12$. Time compression is therefore possible, but it does not solve the native-decomposition fidelity ceiling.
- **Wigner validation:** the final cavity Wigner functions do not match the target well enough to support the bounded decomposition; the logical-one sectors are especially poor ($F_\mathrm{cav}=0.544$ for $|g,1\rangle$ and $0.499$ for $|e,1\rangle$).
- **Artifacts now present:** `artifacts/best_strategy_B.json`, `artifacts/best_strategy_B_duration_penalized.json`, `artifacts/best_strategy_D.json`, `artifacts/wigner_comparison.npz`, `data/iteration4_results.json`, and `data/iteration4_strategy_b_duration.json`.

## Known Limitations
- **Bounded displacement does not rescue a high-fidelity native decomposition** — even the best bounded embedded result at $N_\mathrm{cav}=12$ stays below $0.55$ fidelity, so GRAPE remains the only physically validated high-fidelity route in this study.
- **Iteration 4 uses enlarged-Hilbert-space ideal-gate embedding, not a full new drift-aware re-optimisation at $N_\mathrm{cav}=12$** — leakage is isolated cleanly, but a full high-dimensional SQR/CP optimisation remains open.
- GRAPE optimization used only 3 seeds × 300 iterations for the original sweep — the fidelity frontier may not be fully converged.
- No decoherence channels included in the direct optimisation loops (the coherence budget is analytical only).
- Waveform-level simulation not completed (deferred).

## Status
COMPLETE
