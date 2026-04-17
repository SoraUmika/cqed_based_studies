# Fock-Resolved Qubit-State Inference for a Black-Box SQR Gate

## Problem Class
ANA | DES

## Motivation
The target question is whether a black-box selective qubit rotation (SQR) can be
validated experimentally without reconstructing the full joint qubit-cavity
process. The allowed analysis operations are intentionally limited to ones that
can plausibly be calibrated independently of the black-box gate: a known
dispersive wait under $(\chi,\chi')$, a calibrated cavity displacement
$D(\alpha)$, and ordinary qubit tomography rotations. The study asks whether
those operations are sufficient to infer the effective Fock-resolved qubit output
states on the low-lying logical sectors, and if not, exactly which information is
recoverable and which information remains hidden.

## Goals
1. Implement and validate a single-qubit forward-simulation and MLE tomography
   baseline, including a fidelity-versus-shot-count scaling curve.
2. Derive the exact forward model for qubit tomography after the analysis
   sequence `displacement -> wait -> qubit pre-rotation -> z readout`, and
   identify the parameter combinations that are and are not observable.
3. Compare three protocol families: wait-only, displacement-only, and combined
   displacement-plus-wait.
4. Test diagonal, coherent, pulse-level, leakage-prone, and noisy black-box SQR
   output cases against solver ground truth.
5. Quantify when constrained least squares and shot-statistics MLE recover the
   correct transverse Fock-resolved information, and when a full
   `{p_n, rho_q^(n)}` fit becomes non-unique.
6. Produce an experimentally honest recommendation that states exactly what this
   protocol can validate, what it cannot validate, and what extra measurement
   primitive would be required to close the gap.

## Methods
- Use `cqed_sim.core.DispersiveTransmonCavityModel` as the physical model for all
  pulse-level black-box validation cases.
- Use `cqed_sim.core.ideal_gates.sqr_op` to construct controlled ideal SQR-like
  reference operators for analytic identifiability tests.
- Use `cqed_sim.calibration.conditioned_multitone` to build realistic multitone
  SQR-like waveforms for near-ideal and intentionally imperfect pulse-level
  black-box cases.
- Use `cqed_sim.sim.prepare_simulation`, `cqed_sim.sim.simulate_sequence`, and
  `cqed_sim.sim.extractors.conditioned_qubit_state` /
  `conditioned_population` to obtain trusted ground-truth post-gate states.
- Use a local inference layer for the inverse problem because the current
  `cqed_sim` public API does not expose the required black-box wait/displacement
  tomography MLE.

## Analytic Preliminary
For the Fock-diagonal model
\[
\rho_{\mathrm{out}}=\sum_{n=0}^{N} p_n |n\rangle\langle n|\otimes \rho_q^{(n)},
\]
the allowed analysis operations are:
\[
\rho \mapsto U_{\mathrm{wait}}(t)\,(I\otimes D(\alpha))\,\rho\,(I\otimes D^\dagger(\alpha))\,U_{\mathrm{wait}}^\dagger(t),
\]
with
\[
U_{\mathrm{wait}}(t)=\sum_{m=0}^{N} |m\rangle\langle m|\otimes
R_z(\phi_m(t)),\quad
\phi_m(t)=[\chi m+\chi' m(m-1)]t.
\]
Tracing out the cavity gives
\[
\rho_q^{\mathrm{meas}}(\alpha,t)=
\sum_{n=0}^{N} p_n \sum_{m=0}^{N}
|\langle m|D(\alpha)|n\rangle|^2\,
R_z(\phi_m)\rho_q^{(n)}R_z^\dagger(\phi_m).
\]
Two immediate consequences follow.

1. `Displacement-only` is exactly uninformative: when `t = 0`,
\[
\rho_q^{\mathrm{meas}}(\alpha,0)=\mathrm{Tr}_c(\rho_{\mathrm{out}}),
\]
because the local cavity unitary disappears under the partial trace.
2. Even the combined `D(\alpha) -> wait` protocol preserves the total qubit
   population observable:
\[
Z(\alpha,t)=\sum_n p_n Z_n,
\]
so the allowed operations cannot identify the Fock populations `p_n` or the
sector-wise longitudinal components `Z_n` separately. What remains observable is
the weighted transverse sector content
\[
u_n \equiv p_n (X_n + iY_n),
\]
because
\[
X(\alpha,t)+iY(\alpha,t)=\sum_n u_n K_n(\alpha,t),
\quad
K_n(\alpha,t)=\sum_m |\langle m|D(\alpha)|n\rangle|^2 e^{i\phi_m(t)}.
\]
The study therefore proceeds with two goals in parallel: recover the observable
transverse sector amplitudes as accurately as possible, and document the
non-identifiability of `{p_n, Z_n}` under the stated experimental constraints.

## cqed_sim Gap Analysis
| Functionality | Needed? | Available in `cqed_sim`? | Plan |
|---|---|---|---|
| Pulse-level dispersive qubit-cavity simulation | Yes | Yes | Use `DispersiveTransmonCavityModel` plus `simulate_sequence` / `prepare_simulation` |
| Exact conditioned qubit-state extraction from joint states | Yes | Yes | Use `conditioned_qubit_state`, `conditioned_population`, and related extractors |
| Ideal SQR-like block-diagonal operators | Yes | Yes | Use `sqr_op` for controlled identifiability cases |
| Multitone SQR-like waveform generation | Yes | Yes | Use `conditioned_multitone` utilities |
| Single-qubit MLE tomography baseline | Yes | No public MLE helper | Implement locally in study scripts |
| Black-box wait/displacement inverse model and MLE | Yes | No | Implement locally and document as a reusable upstream target |
| Full `{p_n, rho_q^(n)}` identifiability diagnostics under wait/displacement-only measurements | Yes | No | Implement locally and use it to expose the null directions explicitly |

## Assumptions
- Physical model: two-mode dispersive transmon-storage system in angular-frequency
  units (rad/s) and seconds unless a `cqed_sim.tomo`-specific helper explicitly
  uses ns units.
- Logical sector window: low-lying Fock sectors `n = 0,1,2,3`, with
  convergence spot checks against larger cavity truncations.
- Analysis order: the informative combined protocol is defined as
  `displacement -> known dispersive wait -> qubit tomography rotation -> z readout`.
- Noise model for synthetic tomography data: shot noise is always included; small
  rotation, displacement, and Hamiltonian-calibration errors are added in the
  robustness study; wait-time decoherence is treated with a local analytic qubit
  channel because the current inverse model is not a first-class `cqed_sim`
  routine.
- Success is judged on recoverable quantities first. Full-state recovery claims
  are allowed only if an added assumption or extra measurement primitive removes
  the analytic nullspace.
- Convergence targets:
  - single-qubit baseline fidelity curves should be stable under repeated Monte
    Carlo averaging and should approach unity with increasing shot count;
  - transverse reconstruction conclusions should be stable under denser wait
    grids, modest alpha-grid changes, and cavity truncation enlargement;
  - pulse-level near-ideal versus imperfect-case rankings should survive a
    larger transmon truncation check.

## Compute & Resource Strategy
- Expected heavy step: one conditioned-multitone optimization to synthesize a
  near-ideal pulse-level black-box SQR case. Scratch testing indicates a
  moderate full-mode optimization run is on the order of one minute.
- Expected light steps: all identifiability sweeps, MLE comparisons, and
  protocol comparisons are small dense linear-algebra problems and should run in
  seconds to a few minutes.
- Acceleration plan: keep the pulse-level workload intentionally small (one
  optimized case, one nominal case, one noisy replay); use vectorized NumPy /
  QuTiP evaluation for the inverse-model sweeps; avoid a large Monte Carlo over
  full solver trajectories.
- Expected bottleneck: multiple random restarts of the non-identifiable
  full-state MLE. Restrict those diagnostics to one representative diagonal case
  because the goal is to expose the nullspace, not to optimize a production
  workflow around a fundamentally underdetermined model.

## Expected Outcomes
- Wait-only tomography will recover the weighted transverse sector amplitudes
  `p_n X_n` and `p_n Y_n` accurately when the true post-gate state is
  Fock-diagonal and the dispersive frequencies are well separated.
- Displacement-only will provide no additional information beyond ordinary qubit
  tomography.
- Combined `D(\alpha) -> wait` will not resolve the `p_n / Z_n` non-identifiability,
  but it will serve as a strong diagnostic for cavity-coherence-induced model
  mismatch because it reintroduces off-diagonal cavity blocks into the reduced
  qubit signal.
- Constrained least squares and binomial MLE will agree closely on recoverable
  transverse quantities at moderate shot counts, with MLE slightly more robust at
  low shot budgets.
- A naive full-state Cholesky-plus-softmax fit will show multiple near-equal-loss
  solutions for `{p_n, Z_n}`, confirming the analytic nullspace.
- The final recommendation will be that the stated protocol is useful as a
  transverse-sector black-box diagnostic, but not sufficient for full
  Fock-resolved qubit-state tomography without an additional cavity-sensitive or
  Fock-selective measurement primitive.

## Known Limitations
- The current measurement budget excludes any direct cavity measurement,
  Fock-selective qubit tag pulse, or cavity parity readout, so the study cannot
  honestly claim full sector-population recovery under the stated assumptions.
- The pulse-level black-box cases will be representative rather than exhaustive:
  one optimized multitone case, one nominal/imperfect case, and one noisy replay.
- The robustness model for the inverse problem uses a compact analytic noise
  surrogate during the wait rather than a full open-system tomography replay for
  every protocol setting.

## Validation
- [x] Sanity checks
- [x] Convergence
- [x] Literature comparison (if applicable)
  Original protocol-design study; literature reproduction comparison not applicable.

## Status
COMPLETE
