# Science Directive — Iteration 1
Date: 2026-03-31
Study: studies/fock_resolved_black_box_sqr_inference
Run: task_runs/fock_resolved_black_box_sqr_inference

## Problem Class
ANA | DES

## Core Question
Can ordinary qubit tomography, augmented only by a calibrated cavity displacement
and a known dispersive wait, recover the effective Fock-resolved qubit output
states of a black-box SQR gate on the low-lying cavity sectors?

## Scientific Position Before Numerics
The analytic forward model already rules out the strongest optimistic claim. Under
the allowed operations, the protocol cannot identify `p_n` and `Z_n` separately.
What remains potentially recoverable is the weighted transverse sector content
`p_n X_n` and `p_n Y_n`, plus a combined-protocol residual that can act as a
coherence / model-mismatch witness.

The implementation should therefore pursue two goals simultaneously:
1. Make the strongest positive result that is actually defensible:
   recoverable transverse Fock-resolved information and a useful black-box
   diagnostic.
2. Make the negative result quantitative rather than rhetorical:
   show with analytic rank arguments and numerical MLE restarts that full
   `{p_n, rho_q^(n)}` recovery is non-unique under the stated constraints.

## Hypotheses
1. **Single-qubit baseline**: Cholesky MLE will converge to high fidelity with
   the expected shot-noise scaling and serve as a trustworthy baseline.
2. **Wait-only identifiability**: For Fock-diagonal outputs, wait-only
   tomography is sufficient to recover the weighted transverse sector amplitudes
   `u_n = p_n (X_n + i Y_n)` when the dispersive phase frequencies are distinct.
3. **Displacement-only null result**: Displacement-only tomography produces no
   additional information beyond the ordinary reduced qubit state.
4. **Combined protocol role**: `D(alpha) -> wait` does not remove the
   `p_n / Z_n` nullspace, but it does make cavity coherences visible through
   fit residuals that the diagonal model cannot absorb.
5. **Inference comparison**: On the recoverable transverse subspace, binomial
   MLE should outperform or match least squares at low shot counts, while both
   should agree at moderate to high shots.
6. **Pulse-level realism**: The analytic identifiability conclusions should
   survive when the post-gate states come from actual `cqed_sim` multitone
   waveforms rather than only from ideal block-diagonal operators.

## Required Outputs
1. A single-qubit baseline dataset with fidelity-versus-shot-count curves for
   least squares and MLE.
2. A protocol-comparison dataset for wait-only, displacement-only, and combined
   protocols on controlled Fock-diagonal states.
3. A coherence-witness dataset showing that wait-only is blind to cavity
   coherences while combined `D(alpha) -> wait` reveals them through residuals.
4. Pulse-level black-box case studies:
   - near-ideal multitone SQR-like waveform,
   - nominal/imperfect multitone waveform,
   - at least one noisy replay,
   - at least one leakage-sensitive replay or truncation-sensitive replay.
5. Figures covering the minimum requested set:
   - protocol schematic,
   - single-qubit fidelity scaling,
   - per-sector comparison on recoverable quantities,
   - residual diagnostics,
   - protocol comparison,
   - robustness / noise sensitivity,
   - success and failure examples.
6. A final report whose conclusion explicitly distinguishes:
   - what is recoverable,
   - what is not recoverable,
   - which assumptions are essential,
   - what extra measurement primitive would be needed for full-state recovery.

## Black-Box Case Matrix
### Controlled analytic cases
1. Ideal block-diagonal SQR on a cavity-diagonal input mixture.
2. CPSQR-like conditional phase profile on transverse qubit inputs.
3. Ideal block-diagonal SQR on coherent / superposition cavity inputs to produce
   cavity coherence without pulse imperfections.

### Pulse-level `cqed_sim` cases
1. Optimized multitone SQR-like waveform on a low-dimensional model (`n_tr = 2`)
   as the near-ideal pulse-level case.
2. Nominal or shortened multitone waveform as the imperfect pulse-level case.
3. Noisy replay of the optimized waveform with shot noise, tomography rotation
   error, displacement calibration error, and mild chi / chi-prime mismatch in
   the synthetic data generator.

## Implementation Constraints
1. Use `cqed_sim` for all joint-state generation and all pulse-level black-box
   validation. Do not replace those steps with an ad hoc simulator.
2. Document the local inverse-model layer as a `cqed_sim` gap rather than
   pretending it already exists upstream.
3. Treat every claim about `p_n` or `Z_n` recovery as false unless an explicit
   extra assumption or extra observable removes the analytic nullspace.
4. Keep the pulse-level workload modest. One optimized waveform plus a few
   controlled replays is enough if the analytic case matrix already establishes
   the identifiability boundary.

## Quantitative Success Criteria
1. Single-qubit MLE: mean fidelity should rise monotonically with shot count and
   exceed 0.995 at the high-shot end of the baseline sweep.
2. Fock-diagonal transverse inference: wait-only and/or combined MLE should
   recover the weighted transverse sector amplitudes with small error on the
   diagonal benchmark cases.
3. Displacement-only null result: reconstructed transverse rank or fit quality
   should clearly show non-identifiability.
4. Full-state MLE diagnostic: multiple random restarts should yield materially
   different `{p_n, Z_n}` assignments at nearly equal objective value.
5. Coherence witness: combined-protocol diagonal-model residuals should separate
   coherent from diagonal cases by a clear margin.
6. Pulse-level validation: the near-ideal multitone case should behave
   consistently with the analytic identifiability picture, even if its full
   sector states are not exactly ideal.

## Immediate Execution Plan
1. Build reusable study utilities:
   - runtime/path compatibility,
   - model/constants helpers,
   - Bloch/fidelity utilities,
   - serialization helpers.
2. Implement the single-qubit baseline:
   - forward simulator,
   - Cholesky state parameterization,
   - LS and binomial MLE solvers,
   - shot-count sweep.
3. Implement the black-box inverse-model library:
   - exact wait/displacement kernels for the diagonal model,
   - recoverable transverse LS / MLE,
   - full-state Cholesky-plus-softmax diagnostic fit,
   - coherence-witness residual analysis.
4. Implement case generators:
   - ideal operator cases using `sqr_op`,
   - pulse-level cases using `conditioned_multitone`.
5. Run the study, save machine-readable outputs, generate figures, validate
   convergence, and write the report.
