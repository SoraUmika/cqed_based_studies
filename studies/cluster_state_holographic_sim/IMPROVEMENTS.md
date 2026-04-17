# Improvement Log: Cluster-State Holographic Sequential Simulation

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)

- **[P1 | HIGH] SNAP gate is cavity-only in cqed_sim ideal mode**: `SNAP(phases).ideal_unitary(n_cav, n_tr)` applies identical phases to |g,n⟩ and |e,n⟩, meaning `I_q ⊗ diag(e^{iθ_n})`. The Displacement–Rotation–SNAP gate set therefore **cannot generate any entangling unitary** in ideal mode. Either cqed_sim must implement qubit-conditional SNAP, or a different entangling primitive must be used. **UPDATE (Iteration 2)**: This is now well-understood and documented. The D+SQR+CP and D+R+FE gate sets bypass this limitation entirely with F=1.000 and F=0.9999 respectively.

- **[P1 | HIGH] Re-optimise parametric strategies at N_cav=8+ with full Fock-level support**: Iteration 3 revealed that ideal-mode results (N_cav=2) are truncation artifacts. Strategy B drops from F=1.000 to F=0.094 at N_cav=8 (93% leakage). Re-optimisation with extended SQR/CP parameters (all 8+ Fock levels, ~90 parameters) is needed. Preliminary attempts (multistart=4, maxiter=300) did not converge within compute budget. Gradient-based methods may be required. **UPDATE (Iteration 4)**: The bounded-displacement hypothesis was tested directly at enlarged cavity dimension. Tight bounds help leakage but do not rescue the target unitary: at N_cav=12, B2 with |α|≤0.3 gives F=0.5291 with 3.99% leakage, while D3 with |α|≤0.3 gives F=0.5465 with 22.8% leakage. Relaxing to |α|≤0.5 or unbounded recovers higher ideal-mode fidelity but drives leakage back to 40–55%. Conclusion: low-population operation alone is insufficient; the gate family needs more expressivity, not just smaller displacements.

## Recommended Improvements (P2)

- **[P2 | LOW] Extend GRAPE sweep to longer durations**: Only 6 duration points (50–400 ns) were computed. Points at 500, 600, 800 ns would reveal the asymptotic fidelity limit and confirm convergence. **UPDATE (Iteration 3)**: RESOLVED — 9 points (50–800 ns) already exist in results_combined.json. F>0.99999 for ≥500 ns.

- **[P2 | MEDIUM] Add decoherence channels**: All simulations are unitary (no T1, T2, Tphi). Including Lindblad dissipators would give realistic fidelity estimates and identify the optimal gate time balancing coherence loss vs control fidelity. **UPDATE (Iteration 3)**: First-order coherence budget computed analytically (T1=30μs, T2=20μs). GRAPE 200ns has best combined fidelity F_comb=0.987. Full Lindblad integration deferred to P2.

- **[P2 | MEDIUM] SQR-based decomposition study**: SQR (Selective Qubit Rotation) IS qubit-conditional and could serve as the entangling primitive. **UPDATE (Iteration 2)**: RESOLVED — D+SQR+CP with 2 blocks achieves F=1.000 (exact ideal decomposition, 9 gates). See `scripts/decomposition_comparison.py`.

- **[P2 | MEDIUM] FreeEvolveCondPhase wait-time analysis**: Strategy D (D+R+FE) achieves F=0.9999 at 2 blocks. The optimised wait times should be compared to τ=π/χ≈176 ns (the CZ-equivalent interaction time) to verify physical consistency. **UPDATE (Iteration 3)**: RESOLVED — FE0=FE1=177.0 ns = 1.005 × τ_CZ. Perfect match.

- **[P2 | MEDIUM] Hybrid analytical + GRAPE strategy**: Use D+SQR+CP ideal decomposition as initialization for short GRAPE refinement in the full dispersive model. This could combine exactness of Strategy B with Hamiltonian-awareness of GRAPE. **UPDATE (Iteration 3)**: Motivation strengthened — ideal decomposition fails at N_cav=8. Warm-started GRAPE from ideal parameters is a promising path but untested.

- **[P2 | LOW] Increase Hilbert space truncation**: N_cav=8 was used. For the 2D computational subspace this is likely sufficient, but convergence should be checked at N_cav=12 and N_cav=15. **UPDATE (Iteration 3)**: RESOLVED — evaluated at N_cav=8, 12, 15. Fidelity converges by N_cav=12 (|F(12)-F(15)| < 0.001). Poor fidelity is physical, not numerical.

- **[P2 | MEDIUM] Include Lindblad dissipators in GRAPE**: Run open-system GRAPE with T1, T2, κ to obtain fidelity estimates that include decoherence directly in the optimization, rather than relying on the first-order coherence estimate.

- **[P2 | LOW] Fine GRAPE resolution near optimum**: Sweep 150–250 ns in 10 ns steps to precisely locate the combined-fidelity maximum (currently estimated at ~200 ns).

## Nice-to-Haves (P3)

- **[P3 | MEDIUM] Full holographic sequential measurement simulation**: The holographic channel and transfer matrix were computed analytically. A full simulation running N sequential qubit-cavity interactions + qubit measurements would demonstrate the protocol end-to-end with cqed_sim dynamics.

- **[P3 | LOW] Waveform-level simulation of GRAPE pulses**: The GRAPE solver finds control amplitudes but we did not propagate these through the full Hamiltonian with realistic pulse shapes. This would verify the GRAPE fidelity is achievable with finite bandwidth.

- **[P3 | LOW] Visualization of optimal GRAPE pulses**: Save and plot the I/Q waveforms for each GRAPE-optimized duration.

## Open Questions

- Why does 150 ns (F=0.956) barely improve over 100 ns (F=0.949)? Is there a control bandwidth limitation near 100–150 ns, or is this purely a GRAPE convergence issue?
- The transfer matrix has eigenvalues {1,0,0,0} implying infinite correlation length. For a finite chain, what is the effective correlation structure for N<10 sites?
- **RESOLVED (Iteration 2)**: Can ConditionalPhaseSQR + Displacement achieve the cluster target with fewer layers than SQR + Displacement? → YES. D+SQR+CP achieves F=1.000 at 2 blocks (9 gates).
- Why does Strategy D (D+R+FE) at 2 blocks saturate at F=0.9999 rather than reaching F=1.000? Is this an inherent limitation of the FreeEvolveCondPhase interaction structure, or a convergence issue of the optimiser?
- What is the minimal number of D+SQR+CP blocks for other entangling target unitaries (e.g., CNOT, iSWAP)?
- Does a full N_cav=12 SQR/CP parameterization with explicit control over n=2,...,11 recover high fidelity, or is the cluster-state transfer unitary fundamentally inefficient in this discrete gate family once higher Fock support is admitted?

## What Was Tried and Did Not Work

- **D-R-SNAP decomposition (v3 script, Phase 3)**: UnitarySynthesizer with Displacement+QubitRotation+SNAP primitives. 2-SNAP result: fidelity stuck at 0.500. 3-SNAP: hung for >10 min. Root cause: SNAP does not entangle qubit and cavity (verified analytically). Any decomposition using only D, R, SNAP will map to (arbitrary qubit rotation) ⊗ (arbitrary cavity displacement+phase), which is always separable.

- **SQR synthesis with Powell optimizer (v4 script)**: 6 layers of Displacement+SQR = 42 optimization parameters. Powell optimizer made ~100 function evaluations in >5 min without reaching F>0.6. The cost landscape is high-dimensional with many local minima. GRAPE or L-BFGS-B with analytic gradients would be more suitable.

- **Holographic state construction from MPS tensors**: Attempted to build the N-qubit cluster state by applying the MPS unitary sequentially to a qubit-cavity register initialized in |0⟩⊗|0⟩, measuring the qubit, and repeating. This produces |+⟩^⊗N (product state), NOT the cluster state, because the bond qubit initialized as |0⟩ is the fixed point of the transfer matrix. The correct approach requires either: (a) post-selection on measurement outcomes, or (b) computing observables via the transfer-matrix formalism without constructing the full state.
- **Iteration 4 bounded-displacement sweep**: Tested the user hypothesis that smaller displacements would keep the dynamics in the logical n=0,1 sector and therefore preserve fidelity at enlarged cavity dimension. Result: not sufficient. For Strategy B, |α|≤0.3 reduces leakage to 0.0399 at N_cav=12 but only achieves F=0.5291. For Strategy D, the best bounded point is |α|≤0.3 with F=0.5465 and leakage 0.2278. Increasing the bound to 0.5 improves the N_cav=2 optimisation metric but immediately drives 48% leakage and collapses F(N_cav=12) back to ≈0.48–0.50. The problem is not merely excessive displacement amplitude; the bounded sequence cannot express the target well enough.
- **Duration-penalized bounded follow-up**: Applied duration weights to both bounded families. For Strategy B, every gate compressed to 20 ns (180 ns total active time) yet F(N_cav=12) changed only from 0.529060 to 0.530455 and leakage from 0.039945 to 0.036845. For Strategy D, F(N_cav=12) changed only from 0.546516 to 0.544929. Time compression is therefore possible, but it does not remove the expressivity limit of the bounded decompositions.

## Compute & Resource Notes

- GRAPE optimization wall times (3 seeds × 300 iterations each, N_cav=8, N_tr=2):
  - 50 ns (10 steps): 75 s
  - 100 ns (20 steps): 76 s
  - 150 ns (30 steps): 100 s
  - 200 ns (40 steps): 130 s
  - 300 ns (60 steps): 214 s
  - 400 ns (80 steps): 275 s
- UnitarySynthesizer with SNAP: 2-layer ~30 s, 3-layer >600 s (did not converge)
- **Iteration 2 — multi-decomposition comparison** (UnitarySynthesizer, ideal mode, N_cav=2):
  - A (D+R+SNAP, 7 gates): 5.6 s (multistart=8, maxiter=400)
  - B (D+SQR+CP, 1 block, 5g): 11.1 s (multistart=8, maxiter=500)
  - B (D+SQR+CP, 2 blocks, 9g): 116.7 s (multistart=8, maxiter=500)
  - B (D+SQR+CP, 3 blocks): 198.3 s
  - B (D+SQR+CP, 4 blocks): 208.0 s
  - C (SQR+CP, no D, 6g): 6.4 s (multistart=12, maxiter=600)
  - D (D+R+FE, 2 blocks, 8g): 105.3 s (multistart=8, maxiter=500)
  - D (D+R+FE, 3 blocks): 166.9 s
  - D (D+R+FE, 4 blocks): 184.5 s
- Python 3.12.10, cqed_sim (local), NumPy, SciPy, matplotlib
- Installed seaborn via `pip install --user seaborn` (used in some figure scripts)
- **Iteration 3 — model-based verification & coherence budget**:
  - iteration3_evaluate.py: ~2 min for all N_cav evaluations (simulate_sequence)
  - iteration3_coherence_budget.py: <5 s (analytical, no simulation)
  - iteration3_figures.py: <3 s
  - iteration3_validation.py: <1 s
  - N_cav=8 re-optimization of Strategy B (multistart=4, maxiter=300): did not complete within 30+ min budget (~90 parameters vs ~18 at N_cav=2)
  - Key bottleneck: UnitarySynthesizer at N_cav=8 with extended SQR/CP parameters (~90 params) is ~5-10x slower per iteration than N_cav=2 (~18 params)
- **Iteration 4 — bounded-displacement sweep and Wigner validation**:
  - cqed_sim import from the Box-mounted repository costs ~93 s before any optimisation starts.
  - B2_amp0.3: 340 s; B2_amp0.5: 205 s; B2_unbounded: 200 s.
  - B2_amp0.3 duration-penalized follow-up: 333 s; all 9 gate durations compressed to 20 ns.
  - D3_amp0.3: 308 s; D3_amp0.5: 170 s; D3_unbounded: 304 s.
  - Duration-penalized follow-up on the best bounded candidate: 186 s.
  - Cross-N_cav embedding (4,6,8,12), figure generation, and Wigner extraction complete within the same v6 run after the six optimisations.
  - Wigner calculation emits a QuTiP singular-matrix warning for some pure-target fidelity evaluations, but figures and NPZ artifacts are still written successfully.

## Resolved

- **Import errors in v1/v2 scripts**: `rotation_xy` renamed to `qubit_rotation_xy` in cqed_sim. Fixed by reading API reference.
- **HolographicChannel import path**: Not in `cqed_sim.holographic_sim` but in `cqed_sim.quantum_algorithms`. Fixed by grep of source code.
- **SQR-based decomposition (Iteration 2)**: Previously failed with 42-parameter Powell optimization. Resolved by using UnitarySynthesizer with GateSequence of D+SQR+CP blocks (multistart=8, maxiter=500). F=1.000 at 2 blocks.
- **GRAPE stochasticity**: Previously P1 due to only 3 seeds. Now contextualized: GRAPE at 400ns achieves F=0.999, but discrete gate strategies B and D both exceed this in ideal mode. The GRAPE stochasticity concern is less critical given the availability of exact decomposition alternatives.
- **ConditionalPhaseSQR question**: Previously open question. D+SQR+CP with 2 blocks achieves F=1.000, confirming CP+SQR+D is sufficient.
- **[P1 | MEDIUM] Ideal-mode verification (Iteration 3)**: RESOLVED by evaluating ideal-mode parameters at N_cav=8. Finding: ideal-mode results are truncation artifacts—Strategy B drops from F=1.000 to F=0.094 (93% leakage). This was the highest-priority improvement and is now fully characterized.
- **[P2 | LOW] Extended GRAPE durations (Iteration 3)**: RESOLVED — 9 durations (50–800 ns) already in results_combined.json. F converges to >0.99999 at ≥500 ns.
- **[P2 | MEDIUM] FE wait-time analysis (Iteration 3)**: RESOLVED — FE0=FE1=177.0 ns = 1.005 × τ_CZ (π/|χ|). Physically consistent.
- **[P2 | LOW] Hilbert space convergence (Iteration 3)**: RESOLVED — evaluated at N_cav=8, 12, 15. Converges by N_cav=12. Poor fidelity is physical.
- **[P2 | MEDIUM] First-order decoherence budget (Iteration 3)**: RESOLVED — analytical coherence budget computed. GRAPE 200ns: F_comb=0.987 (best). Full Lindblad integration remains as future P2.
