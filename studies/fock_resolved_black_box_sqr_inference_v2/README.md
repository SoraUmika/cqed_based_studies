# Fock-Resolved Qubit-State Inference for a Black-Box SQR Gate (v2)

## Problem Class
ANA | DES

## Motivation

A black-box SQR gate acts on a joint qubit+cavity system. We do not assume we
know its waveform-level implementation. We do assume:
- dispersive parameters $\chi$ and $\chi'$ are known,
- cavity displacement $D(\alpha)$ is calibrated,
- qubit tomography rotations are calibrated,
- a numerical solver can serve as a trusted validation backend.

The scientific goal is to determine whether qubit-state tomography, augmented by
known dispersive free evolution and calibrated cavity displacement, can reliably
recover the **effective Fock-resolved qubit output states**
$\rho_q^{(n)} \propto \langle n|\rho_{\rm out}|n\rangle$
for the relevant low-lying Fock sectors — and if so, via which protocol and under
what conditions.

This is v2. v1 (in `fock_resolved_black_box_sqr_inference`) pre-committed to the
recoverable subspace analytically and validated that conclusion numerically. v2
follows the investigative arc specified in the open study plan: it first attempts
full $\{p_n,\rho_q^{(n)}\}$ recovery, documents the non-uniqueness failure mode,
then identifies and validates the recoverable subspace. It also adds per-sector
fidelity and trace-distance metrics, a systematic Model B coherence sweep, the
full planned pulse-level case library, and a per-Fock Bloch-vector comparison
figure that was absent in v1.

## Goals

1. Implement and validate a single-qubit MLE tomography baseline including a
   $\mathcal{F}$ vs. $N$ scaling curve (this engine feeds Part 2).
2. Derive the exact forward model for the analysis sequence
   $D(\alpha) \to \text{wait}(t) \to \text{qubit pre-rotation} \to Z\text{-readout}$
   and identify which parameter combinations are and are not observable.
3. **Attempt full $\{p_n,\rho_q^{(n)}\}$ recovery first** via Cholesky-softmax
   MLE; characterize the non-uniqueness failure mode rigorously.
4. Having demonstrated non-identifiability, characterize the recoverable subspace
   $u_n = p_n(X_n+iY_n)$ and validate its recovery with LS and MLE.
5. Compare wait-only, displacement-only, and combined $D(\alpha)+t$ protocols on
   identifiability, accuracy, and robustness.
6. Characterize **Model B** (general block model with cavity coherences):
   sweep off-diagonal coherence weight systematically and map how accuracy
   degrades.
7. Validate across a full planned black-box case library: near-ideal, distorted,
   CPSQR-like, shorter pulse, longer pulse, leakage, non-ideal spectator.
8. Report per-sector fidelity, trace distance, and Bloch-vector errors for every
   case; produce per-Fock Bloch-vector comparison figures.
9. Answer the six comparison questions from the study plan quantitatively.
10. Produce an experimentally honest recommendation stating exactly what this
    protocol can and cannot validate.

## Methods

- `cqed_sim.core.DispersiveTransmonCavityModel` — physical model for all
  pulse-level black-box validation cases.
- `cqed_sim.core.ideal_gates.sqr_op` — controlled ideal SQR-like reference
  operators for analytic identifiability tests.
- `cqed_sim.calibration.conditioned_multitone` — realistic multitone waveforms
  for near-ideal, distorted, shorter-pulse, and longer-pulse cases.
- `cqed_sim.sim.prepare_simulation`, `cqed_sim.sim.extractors.conditioned_qubit_state`,
  `conditioned_population` — trusted ground-truth post-gate states.
- Local inference layer: Cholesky-softmax full-state MLE and design-matrix
  weighted-transverse LS/MLE (not in public cqed_sim API).

## Analytic Preliminary

For the Fock-diagonal model
$$\rho_{\rm out}=\sum_{n=0}^{N}p_n|n\rangle\langle n|\otimes\rho_q^{(n)},$$
the analysis sequence $D(\alpha)\to U_{\rm wait}(t)$ followed by partial trace gives
$$\rho_q^{\rm meas}(\alpha,t)=\sum_n p_n\sum_m|\langle m|D(\alpha)|n\rangle|^2
R_z(\phi_m)\,\rho_q^{(n)}\,R_z^\dagger(\phi_m),\quad
\phi_m(t)=[\chi m+\chi'm(m-1)]t.$$

Two immediate consequences:

1. **Displacement-only is exactly uninformative** at $t=0$:
   $\rho_q^{\rm meas}(\alpha,0)={\rm Tr}_c(\rho_{\rm out})$, independent of $\alpha$.

2. **The $Z$-projection is invariant** under all allowed operations:
   $Z(\alpha,t)=\sum_n p_n Z_n$,
   so neither $p_n$ nor $Z_n$ individually are observable from these measurements.

The observable is the weighted transverse content:
$$u_n\equiv p_n(X_n+iY_n),\qquad
X(\alpha,t)+iY(\alpha,t)=\sum_n u_n K_n(\alpha,t),$$
$$K_n(\alpha,t)=\sum_m|\langle m|D(\alpha)|n\rangle|^2 e^{-i\phi_m(t)}.$$

The study proceeds by first attempting to recover the full
$\{p_n,\rho_q^{(n)}\}$ via MLE, observing the nullspace empirically, then
validating recovery of the achievable target $\{u_n\}$.

## cqed_sim Gap Analysis

| Functionality | Needed? | In `cqed_sim`? | Plan |
|---|---|---|---|
| Dispersive qubit-cavity simulation | Yes | Yes | `DispersiveTransmonCavityModel` + `prepare_simulation` |
| Conditioned qubit-state extraction | Yes | Yes | `conditioned_qubit_state`, `conditioned_population` |
| Ideal SQR-like operators | Yes | Yes | `sqr_op` |
| Multitone waveform generation | Yes | Yes | `conditioned_multitone` |
| Single-qubit MLE baseline | Yes | No | Implement locally |
| Black-box wait/displacement MLE | Yes | No | Implement locally; document as upstream target |
| Full $\{p_n,\rho_q^{(n)}\}$ identifiability diagnostics | Yes | No | Implement locally |
| Per-sector fidelity / trace-distance reporting | Yes | Partial extractors | Compute locally from extractor output |
| Coherence-sweep degradation characterization | Yes | No | Implement locally |

## Assumptions

- Two-mode dispersive transmon-storage model in rad/s and seconds.
- Logical sector window: $n=0,1,2,3$ ($N_{\rm active}=4$), spot-checked against
  larger truncations.
- Analysis order: $D(\alpha)\to{\rm wait}(t)\to{\rm qubit\;pre\text{-}rotation}
  \to Z{\rm\text{-}readout}$.
- Noise: shot noise always included; rotation, displacement, and Hamiltonian
  calibration errors in robustness sweep; wait-time decoherence modelled with
  analytic qubit channel (not full open-system sim).
- Success judged on recoverable quantities ($u_n$) first; full-state recovery
  claims require additional constraints or measurement primitives.
- Convergence targets: single-qubit baseline $\mathcal{F}\to1$ with shot count;
  transverse inference conclusions stable under denser grids and modest truncation
  enlargement; near-ideal vs. imperfect pulse rankings survive $n_{\rm tr}$ check.

## Compute & Resource Strategy

- Heavy step: two conditioned-multitone optimizations (shorter and longer pulse).
  Each ≈ 1 min. Run sequentially.
- Light steps: identifiability sweeps, LS/MLE sweeps, coherence-degradation sweep
  — all dense linear-algebra, seconds to a few minutes.
- Acceleration: vectorized NumPy/QuTiP evaluation; restrict full-state restart
  diagnostics to one representative case (goal is to expose nullspace, not
  optimize an underdetermined model).
- Expected total runtime: 3–5 min on active workstation with full profile.

## Expected Outcomes

- Full $\{p_n,\rho_q^{(n)}\}$ MLE will show multiple near-equal-loss solutions
  across random restarts, confirming the analytic nullspace.
- Weighted transverse $u_n$ is fully recoverable via wait-only or combined; LS
  and MLE agree closely at moderate shot counts.
- Displacement-only provides no information beyond ordinary qubit tomography.
- Combined protocol residuals serve as a coherence witness: small for diagonal
  states, large when off-diagonal cavity blocks are present.
- Model B coherence sweep will show a threshold below which inference is still
  useful and above which residuals flag the model violation.
- Per-sector fidelity and trace-distance metrics will track closely with weighted
  RMSE; the per-Fock Bloch-vector figure will visually confirm agreement for
  near-ideal cases and failure for leakage/coherent cases.
- Final recommendation: use combined protocol for transverse sector recovery plus
  coherence detection; do not claim $p_n$/$Z_n$ recovery without extra primitive.

## Known Limitations

- No direct cavity measurement, parity readout, or Fock-selective tag pulse; full
  sector-population recovery is analytically impossible under the stated budget.
- Pulse-level cases representative, not exhaustive; CPSQR-like waveform uses
  analytic phase gate rather than a true CPSQR optimization.
- Decoherence model during wait is a compact analytic channel, not a full
  open-system replay.

## Validation

- [ ] Sanity checks
- [ ] Convergence
- [ ] Literature comparison (N/A — original study)

## Status
COMPLETE
