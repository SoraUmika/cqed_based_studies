# Science Directive: Hybrid Qubit-Cavity Control — Iterations 5–14

## Decision
CONTINUE (10-round auto-research loop requested by user)

## Objective

Extend, reorganize, and strengthen the hybrid_qubit_cavity_control study by:
1. Fixing upstream cqed_sim framework gaps that blocked prior work
2. Implementing ConditionalPhaseSQR waveform bridge using SQR multitone structure
3. Re-verifying all prior results under the fixed framework
4. Running previously-blocked studies (n_tr=3 inner loop, L2d replay, full-sequence GRAPE)
5. Producing a reorganized, comprehensive report

## Assessment of Previous Results

**Strong foundation (iterations 1–4):**
- GRAPE transforms the selective-gate frontier: individual gates reach >99% fidelity
- Hardware constraints validated: all 9 gates pass realizability checks
- Fine duration grid maps the landscape: L1d_S1 peak at χt≈1.0, B_ent_S at χt≈3.61
- SNAP-equivalent pi-SQR confirmed (F=0.5321 = B_ent_S)

**Three active blockers to resolve:**
1. `bloch_xyz_from_qubit_state` rejects n_tr>2 → blocks qutrit-aware inner loop
2. `waveform_sequence_from_gates` missing CPSQR → blocks L2d replay
3. `_expand_target_matrix` mishandles subspace targets → blocks model-based SU(4)

**Key user directives:**
- "For ConditionalPhaseSQR, assume same multitone function as canonical SQR"
- "Try to update upstream cqed_sim for simple implementations"
- "If not possible, build locally"
- "10 iterations of auto-research"

## Iteration Plan

### Iteration 5 (upstream fix 1): bloch_xyz_from_qubit_state for n_tr>2
- Modify `cqed_sim/sim/extractors.py`: truncate to 2-level subspace when n_tr>2
- Modify `conditioned_bloch_xyz` to handle multilevel transmon
- Test with n_tr=3 model
- This unblocks P1.2

### Iteration 6 (upstream fix 2): ConditionalPhaseSQR waveform bridge
- Add `ConditionalPhaseSQRGate` IO class to `cqed_sim/io/gates.py`
- Add `build_cpsqr_multitone_pulse` to `cqed_sim/pulses/builders.py` (reuse SQR multitone)
- Extend `waveform_bridge.py` to handle ConditionalPhaseSQR
- This unblocks P1.4

### Iteration 7 (upstream fix 3): _expand_target_matrix for subspace UnitarySynthesizer
- Fix `_expand_target_matrix` in `targets.py` to handle 4×4 target with Subspace.custom
- Test with model-based SU(4) synthesis
- This unblocks random SU(4) test

### Iteration 8: Re-verify all prior results
- Spot-check all prior data files
- Re-run key simulations with fixed framework
- Verify GRAPE, extension pass 1 and 2 results unchanged

### Iteration 9: n_tr=3 inner-loop study
- Re-run targeted_subspace_multitone with n_tr=3 (now unblocked)
- Compare n_tr=2 vs n_tr=3 fidelities for all 9 refined gates

### Iteration 10: L2d ConditionalPhaseSQR sequence replay
- Replay L2d sequences through new CPSQR waveform bridge
- Update strategy rankings

### Iteration 11: Full-sequence GRAPE optimization
- GRAPE-optimize composites of the full 12-gate L1d sequence
- Measure sequence fidelity

### Iteration 12: Duration-optimized GRAPE
- Sweep GRAPE pulse duration for each bottleneck gate
- Find optimal duration vs. default 352 ns

### Iteration 13: Validation and convergence
- Full validation: sanity, convergence, literature comparison
- Generate validation figures

### Iteration 14: Report reorganization and final compile
- Reorganize report.tex
- Add new sections, update appendices, compile PDF
- Finalize IMPROVEMENTS.md, mark COMPLETE

## Success Criteria
1. All three framework blockers resolved (upstream or local)
2. n_tr=3 inner loop produces valid results
3. L2d replayed through CPSQR bridge
4. Full-sequence GRAPE fidelity measured
5. Duration-optimized GRAPE available
6. All validation checks pass
7. Report compiled, IMPROVEMENTS.md updated
