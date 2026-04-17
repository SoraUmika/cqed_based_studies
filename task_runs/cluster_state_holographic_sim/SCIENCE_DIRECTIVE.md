# Science Directive — Iteration 3

> Written by Research Loop Orchestrator (acting as Science Director).
> Note: A dedicated high-reasoning model (GPT-5.4 per research_config.json) would produce better scientific planning.

## Objective

Address all remaining P1/P2 improvements from IMPROVEMENTS.md. These are the open items after Iteration 2 that could meaningfully change the study's conclusions.

## Problem Classification

OPT + ANA (Optimization + System Analysis)

## Experiments

### Experiment 1: Model-Based Re-Verification [P1 | MEDIUM] — CRITICAL

**Rationale**: Iteration 2 decomposition results (B: F=1.000, D: F=0.9999) were obtained in ideal mode (χ=K=0, N_cav=2). The physical system has dispersive shift χ/2π = −2.84 MHz, higher-order χ′/2π = −21 kHz, and cavity self-Kerr K/2π = −28 kHz. Gates must be re-optimized in the full dispersive model to confirm viability.

**Method**:
- Set `DriftPhaseModel(chi=CHI, chi2=CHIP, kerr=KERR)` for all SQR, ConditionalPhaseSQR, FreeEvolveCondPhase gates
- Use N_cav=8 (full Hilbert space, 16D total)
- Re-optimize Strategy B (D+SQR+CP, 2 blocks) and Strategy D (D+R+FE, 2 blocks) via UnitarySynthesizer
- Use multistart=12, maxiter=800 (more aggressive than ideal mode to handle the harder landscape)

**Success criterion**: F ≥ 0.99 in model-based mode

### Experiment 2: Extended GRAPE Sweep [P2 | LOW]

**Rationale**: Only 6 duration points (50–400 ns) were tested. Extending to 500, 600, 800 ns maps the asymptotic fidelity limit.

**Method**:
- Run GRAPE at 500, 600, 800 ns using same setup (3 seeds × 300 iterations, N_cav=8, 2 control channels)
- Use `build_control_problem_from_model` with `GrapeSolver(GrapeConfig(...))`

**Success criterion**: Confirm fidelity converges (asymptotic behavior visible)

### Experiment 3: Hilbert Space Convergence [P2 | LOW]

**Rationale**: GRAPE at 400 ns with N_cav=8 gives F=0.999. Convergence should be verified at N_cav=12 and N_cav=15.

**Method**:
- Re-run GRAPE at 400 ns with N_cav=12 and N_cav=15
- Also re-run model-based Strategy B at N_cav=12

**Success criterion**: Fidelity change < 0.1% between N_cav=8 and N_cav=15

### Experiment 4: FE Wait-Time Analysis [P2 | MEDIUM]

**Rationale**: Strategy D uses FreeEvolveCondPhase as the entangling primitive. The optimized wait times should physically correspond to π/|χ| ≈ 176 ns (the CZ-equivalent interaction time for n=1).

**Method**:
- Run Strategy D in model-based mode (Experiment 1 already does this)
- Extract optimized FreeEvolveCondPhase durations from the sequence
- Compare to τ_CZ = π/|χ| and report the ratio

**Success criterion**: Wait times within 2× of τ_CZ

### Experiment 5: Decoherence / Coherence Budget [P2 | MEDIUM]

**Rationale**: All results are unitary. Finite T1, T2 will reduce fidelity for longer gate sequences.

**Method**:
- Analytical coherence budget: F_coh ≈ exp(−t_total/T2) for dephasing, exp(−t_total/(2*T1)) for relaxation
- Compute total gate times for each winning strategy (B-2blocks, D-2blocks, GRAPE at various durations)
- Use T1=30 μs, T2=20 μs (representative cQED values from AGENTS.md default parameters)
- Plot decoherence-limited fidelity vs gate time

**Success criterion**: Identify which strategies remain above 99% fidelity with realistic decoherence

## Ordered Action Items

1. Write `scripts/iteration3_model_based.py` implementing Experiments 1, 3, 4
2. Write `scripts/iteration3_grape_extension.py` implementing Experiment 2
3. Write `scripts/iteration3_coherence_budget.py` implementing Experiment 5
4. Run all scripts, collect results to `data/iteration3_*.json`
5. Generate figures to `figures/`
6. Validate: convergence + sanity checks
7. Update IMPROVEMENTS.md
8. Extend report.tex
9. Compile PDF

## Notes

- The hybrid analytical+GRAPE improvement (P2-MEDIUM) is partially addressed by Experiment 1 (model-based D+SQR+CP), which demonstrates the analytical gate structure works with the full physical Hamiltonian. A full GRAPE warm-start from the D+SQR+CP solution would require translating the gate-based parameterization to control pulses, which is a more complex undertaking best suited for a future iteration.
- Decoherence assessment uses analytical coherence budget rather than full Lindblad simulation, which provides the essential physics (T1/T2 limited fidelity) without requiring the pulse-level bridge from UnitarySynthesizer.
