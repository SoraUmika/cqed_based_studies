# Improvement Log: Runtime-Validated Hybrid Unitary Synthesis with Native Entangling Evolution

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)

- **[P1 | HIGH] Shorten or re-optimize the replayable local surrogates**: the replay gap is closed, but the winning runtime candidates are dominated by three `1.28 us` local-surrogate blocks per sequence. `R2_exact_runtime_to_exact_runtime` and `R2_A_runtime_to_A_runtime` both take `4.432 us` total and fall to nominal-noise average probe fidelities `0.614` and `0.580`. The next study should attack duration first, not architecture.
- **[P1 | MEDIUM] Improve surrogate convergence and target fidelity**: all six Phase 5 GRAPE restarts hit the iteration limit (`maxiter = 60`) with `success = false`. The best nominal local-surrogate fidelities are only `0.94288` for `exact_hc` and `0.94762` for `A_local`, so the current runtime ranking is likely limited by optimizer budget rather than a demonstrated local-control ceiling.
- **[P1 | MEDIUM] Enforce a fidelity floor in ranking**: the weighted scalar cost can select the one-wait replay baseline `R1_exact_runtime_to_exact_runtime` whenever the infidelity weight is reduced by `25%`, even though its process fidelity is only `0.2234`. Future ranking logic should use Pareto reporting or a hard minimum fidelity constraint before cost comparisons.

## Recommended Improvements (P2)

- **[P2 | MEDIUM] Expand the noise and robustness model**: Phase 5 only tested one nominal noise point (`T1_q = 30 us`, `T2_q = 20 us`, `T1_c = 250 us`) at `n_cav = 12`, `n_tr = 3`. Add amplitude, duration, `chi`, and thermal-population perturbations before making any experimental recommendation beyond architecture choice.
- **[P2 | MEDIUM] Map the exact-versus-archive-local tradeoff under shorter controls**: `R2_exact_runtime_to_exact_runtime` wins on process and probe fidelity, while `R2_A_runtime_to_A_runtime` has the higher final Wigner overlap (`0.9814` versus `0.9638`). It is not yet clear whether this ordering is intrinsic or just an artifact of the current surrogate family.
- **[P2 | MEDIUM] Benchmark richer local-control ansatz families**: the current surrogates use four rotating-frame channels, 40 slices, and square-envelope exports. Test longer optimizer budgets, smoother parameterizations, or structured local ansatzes before concluding that `~0.90` closed-system process fidelity is the local-control ceiling.
- **[P2 | LOW] Broaden the device-point sweep**: the ranking is stable across `n_cav = 10, 12, 14`, but only at one nominal dispersive operating point. Vary `chi`, `chi'`, and Kerr terms to see whether the runtime winner is architecture-stable or point-specific.

## Nice-to-Haves (P3)

- Compare a calibration-complexity score based on tone crowding, amplitude span, phase spread, and pulse duration rather than tone count alone.
- Add cached surrogate loading to avoid rerunning the same GRAPE solve when only downstream replay analysis changes.
- Extend the same runtime-validation workflow to a small random-SU(4) benchmark set after `U_target` is settled.

## Open Questions

- Why does `R2_A_runtime_to_A_runtime` preserve final Wigner structure better than `R2_exact_runtime_to_exact_runtime` while still underperforming in process fidelity and average probe fidelity?
- Can the two surrogate-backed two-wait families be shortened below `2 us` total local-control time without losing the current `~0.90` closed-system process fidelity level?
- Are the runtime surrogates exploiting off-logical excursions that help the closed-system metric but amplify decoherence sensitivity during the `4.432 us` total sequence?
- Does a direct model-backed implementation of `SNAP` outperform the current GRAPE surrogate route once replay duration is included in the score?

## What Was Tried and Did Not Work

- **Inherited from the earlier hybrid study**: the ideal `L1d` structured winner collapses under full pulse replay (`pulse_fidelity ~ 0.004`) despite `ideal_fidelity ~ 0.875`. Future agents should not assume that a high ideal structured score translates to a viable physical sequence.
- **Inherited from the earlier hybrid study**: CPSQR-based `L2d` replayed substantially better than the SQR-heavy `L1d` route despite similar gate counts. The relevant discriminator is therefore not gate count alone.
- **Inherited from the earlier hybrid study**: sequence-level GRAPE can overfit its internal propagator, producing a large nominal-vs-strict gap on long sequences. Use strict replay as the decision metric if GRAPE refinement is introduced again.
- **Inherited workflow warning**: the local API reference and the older study notes disagree on bridge support for some conditional gates. Treat documentation and implementation as separate claims until both are tested in this study.
- **Phase 3.1 replay-support check**: the shortlisted `N2_*` native-heavy candidates cannot go through `waveform_sequence_from_gates(...)` as currently written because the bridge rejects `PrimitiveGate`. The installed error message confirms support only for `QubitRotation`, `Displacement`, `SQR`, and `ConditionalPhaseSQR`; it explicitly routes `SNAP` and `FreeEvolveCondPhase` sequences toward the model-backed simulation path instead.
- **Phase 5 surrogate prototype with storage-only controls**: a storage-only GRAPE surrogate was not competitive. It plateaued around `0.46`--`0.48` nominal fidelity for the local-H targets, so the replayable surrogate path had to include both storage and qubit I/Q controls.
- **Phase 5 direct replay of the archive `B_local` family**: bridge-compatible direct replay is not sufficient by itself. `R2_B_local_replay` reaches only `0.16665` process fidelity, `0.18729` average probe fidelity, and `0.04885` nominal-noise average probe fidelity while taking `7.912 us` total runtime.
- **Phase 5 one-wait replay baseline**: `R1_exact_runtime_to_exact_runtime` remains a lower bound after replay as well as symbolically. It has process fidelity `0.2234` and average probe fidelity `0.2120`, so no amount of cost-weight tuning should promote it without a fidelity threshold.

## Compute & Resource Notes

- Bootstrap sources already exist in [studies/hybrid_qubit_cavity_control](studies/hybrid_qubit_cavity_control), so the first pass of this study should reuse those artifacts instead of rerunning the full earlier campaign.
- The new high-cost work items are expected to be: depth-resolved replay, Wigner snapshot generation, and any fresh optimal-control refinement.
- The native-entangler bootstrap should be cheap because the first script only aggregates existing JSON outputs into a new weighted frontier.
- The symbolic Phase 4 diagnostics were inexpensive compared with any pulse-level replay. The largest figure set was the 18-panel Wigner grid for `N2_A_local_to_A_local`, which still completed comfortably on the local machine.
- Phase 5 runtime validation is dominated by surrogate optimization, not by the replay sweeps. Each surrogate uses 4 channels, 40 time slices, and 3 GRAPE restarts at `maxiter = 60`; replay over `n_cav = 10, 12, 14` is comparatively cheap once the surrogate pulses exist.
- The Phase 5 GRAPE payload save path in `cqed_sim.optimal_control.result.save()` is not JSON-safe for complex-valued payloads on this machine. Use the study-local JSON serializer instead of the framework convenience save method when storing runtime artifacts.

## Resolved

- **[Resolved] Close the symbolic-versus-pulse validation gap**: `phase5_runtime_validation.py` now constructs replayable runtime candidates, generates surrogate and decomposition artifacts, and reports symbolic and runtime metrics side by side.
- **[Resolved] Replace or bypass the exact local Hadamard references in the top native-heavy candidates**: exact qubit Hadamards are replaced by replayable `QubitRotation` pairs, and the non-replayable `exact_hc` / `A_local` local blocks are replaced by GRAPE-derived replayable surrogates.
- **[Resolved] Decide whether symbolic validation is sufficient for this study milestone**: the study no longer stops at the symbolic milestone. The ranking has been updated using replay-backed closed-system metrics, truncation checks, and a nominal-noise pass.
- **[Resolved] Generate the missing depth diagnostics**: `phase4_depth_diagnostics.py` now produces checkpointed `X/Y/Z` and Wigner-versus-depth outputs for the one-wait baseline, the exact two-wait upper bound, and the fully archive-driven two-wait candidate.
- **[Resolved] Verify current replay support for conditional free-evolution gates**: `phase3_replay_support_check.py` confirmed the installed bridge behavior directly instead of relying on conflicting notes. `ConditionalPhaseSQR` is bridge-representable in the current build, while the shortlisted native-heavy sequences remain blocked by `PrimitiveGate` plus `SNAP` / `FreeEvolveCondPhase` path limitations.