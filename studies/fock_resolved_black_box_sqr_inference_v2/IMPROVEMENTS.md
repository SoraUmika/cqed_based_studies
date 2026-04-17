# Improvement Log: Fock-Resolved SQR Inference v2

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)

- **[P1 | HIGH]** The stated measurement set `{qubit tomography, cavity displacement,
  dispersive wait}` has a provable analytic nullspace: `p_n` and `Z_n` cannot be
  identified separately. Any future extension that claims full Fock-resolved
  qubit-state recovery must add a new measurement primitive (parity readout,
  Fock-selective tag pulse, etc.) or an explicit prior that removes the degeneracy.
- **[P1 | HIGH]** Per-sector oracle metrics (`oracle_fidelity`, `oracle_trace_distance`,
  `oracle_bloch_error`) require the true `p_n` to normalise the inferred `u_n`.
  They are **oracle** quantities — not directly accessible in experiment. This must
  be emphasised when interpreting per-sector figures.
- **[P1 | HIGH]** Combined-protocol inference must not be interpreted literally when
  the fit residual is large. Large residuals signal diagonal-model failure (cavity
  coherences or leakage), not a large `u_n` error.

## Recommended Improvements (P2)

- **[P2 | MEDIUM]** The CPSQR-like pulse-level case uses an analytic operator
  (`cpsqr_like_operator`) rather than a truly optimized CPSQR waveform. Replace
  with an actual CPSQR pulse optimization once `cqed_sim` gains a matching
  calibration target.
- **[P2 | MEDIUM]** The non-ideal spectator case is simulated by offsetting the
  drive frame frequency, not by including a second physical qubit. A more realistic
  model would add a weak coupling to a spectator mode.
- **[P2 | MEDIUM]** Add Bayesian model comparison between the diagonal model and a
  low-rank cavity-coherence surrogate, so combined residuals become a quantitative
  coherence witness rather than only a heuristic flag.
- **[P2 | MEDIUM]** Upstream the recoverable-subspace inference kernel
  (`diagonal_kernel_vector`, `infer_weighted_ls`, `infer_weighted_mle`) into a
  `cqed_sim.tomo` or `cqed_sim.calibration_targets` helper, with a documented
  statement of the recoverable subspace.

## Nice-to-Haves (P3)

- **[P3 | LOW]** Expand the pulse-level case library with more waveform families
  (multi-target SQR, echoed SQR) once the core identifiability conclusion is stable.
- **[P3 | LOW]** Add a benchmark against an oracle that has access to a cavity
  parity readout, clearly labelled as a stronger measurement setting, to quantify
  how much information the extra primitive provides.
- **[P3 | LOW]** Replace the analytic T1/T2 surrogate in the robustness model with
  a Lindblad open-system simulation for one representative point.

## Open Questions

- If independent cavity-population information is supplied from a separate
  calibration measurement, does the remaining wait/displacement protocol become
  sufficient to estimate `Z_n`, or is an additional operation still required?
- What is the minimum alpha-grid density needed for reliable cavity coherence
  detection under realistic calibration drift?
- How much of the pulse-level failure budget in the imperfect cases comes from
  leakage vs. coherent in-subspace error?
- Does the Model-B coherence sweep threshold (below which inference is still
  useful) depend significantly on the specific SQR rotation angles?

## What Was Tried and Did Not Work

- **Displacement-only Fock-sector tagging**: the exact analytic calculation
  confirms that at `t = 0`, every displacement `D(alpha)` maps to the same
  reduced qubit state `Tr_c(rho_out)` regardless of `alpha`. Displacement-only
  is exactly non-informative and was confirmed numerically to rank-2 by SVD.

## Compute & Resource Notes

- Study initialised and implemented 2026-03-31.
- Planned runtime split: two conditioned-multitone optimizations (~1 min each),
  all identifiability / coherence-sweep / robustness sweeps fast (seconds each).
- Expected total runtime with `STUDY_PROFILE=full`: 3–5 min.

## Resolved

- **`cqed_sim` access path**: active installation at
  `C:\Users\jl82323\Box\Shyam Shankar Quantum Circuits Group\Users\Users_JianJun\cQED_simulation`.
- **Wait-phase sign convention**: `e^{-i phi_n(t)}` confirmed correct (matches
  exact joint simulator to machine precision in v1).
- **Full-state non-identifiability**: demonstrated in v1 via random-restart MLE
  and algebraic gauge family; re-demonstrated in v2 as the primary Part 2B result.
