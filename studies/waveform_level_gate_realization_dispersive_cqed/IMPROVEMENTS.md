# Improvement Log: Unconditional Cavity Displacement in Dispersive cQED

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)

- **[P1 | MEDIUM] No open-system benchmarking.** Every reported protocol is evaluated under unitary evolution only. Add `T1`, `Tphi`, and cavity loss `kappa` before treating the reported protocol ranking as hardware-final. This matters most for the longer naive, echoed, and 80 ns two-tone cases.

- **[P1 | MEDIUM] The optimal-control benchmark is intentionally lightweight.** The hardware-aware GRAPE run uses `n_cav = 12`, only two candidate durations, eight time slices, and a modest optimization state set. It is a meaningful constrained benchmark, but not the ceiling of what optimal control could achieve.

- **[P1 | MEDIUM] Two-tone calibration is vacuum-matched, not logical-state-matched.** The best two-tone protocol nearly equalizes the vacuum branches (`delta_alpha = 8.15e-4` at `T = 20 ns`) and gives excellent superposition performance there, but longer-duration two-tone pulses can still lose fidelity on coherent-state tests. If logical-state performance is the goal, the calibration objective must be broadened beyond vacuum branch matching.

## Recommended Improvements (P2)

- **[P2 | MEDIUM] Precompensate the successful short two-tone pulse rather than low-dimensional reshaping.** The explicit four-segment shaped-two-tone optimization tested here converged but still failed badly on the broad state set. Any further shaping effort should therefore start from the already successful `20 ns` two-tone solution and include hardware-aware precompensation directly, not just a low-dimensional segment reweighting.

- **[P2 | HIGH] If multiplex interpretability is still required, move directly to a materially richer hardware-aware parameterization.** The tested full-duration multicarrier, segmented branch-resonant, and low-parameter jointly optimized shaped-two-tone families all remained poor on the broad state set. A useful next attempt would need either more expressive hardware-aware segmentation, more target states in the fit, or acceptance of full sampled optimal control complexity.

- **[P2 | HIGH] Increase the optimal-control training set.** Include larger coherent amplitudes and at least one nonclassical logical state in the optimization objective. The current constrained solve is the strongest bounded sampled-waveform benchmark (`0.958` mean fidelity), but it is no longer the best tested protocol overall after the broad-state two-tone follow-up.

- **[P2 | MEDIUM] Sweep `chi'`, `K`, and filter cutoff more aggressively.** The present hierarchy sweep shows that `chi` dominates at the representative point, but a broader map is still needed before claiming the same conclusion across devices with different higher-order nonlinearities or transfer functions.

- **[P2 | MEDIUM] Design a manifold-uniform echo rather than reusing a vacuum-calibrated `pi` pulse.** The present echo failure is not proof that all echoed displacement ideas fail; it shows that this practical version fails because the inserted qubit inversion is not uniform across the displaced cavity manifolds.

- **[P2 | LOW] Add an experiment-facing calibration notebook.** The report outlines the calibration steps conceptually, but a follow-up notebook should walk through extracting branch frequencies, solving for two-tone weights, and validating the result with branch-conditioned Wigner data and qubit tomography.

## Nice-to-Haves (P3)

- **[P3 | LOW] Add cat-state or binomial-state validation.** The current test set includes low Fock states and a modest coherent state. A bosonic-code-facing study should also include at least one protected-state example.

- **[P3 | LOW] Add a precompensated fast-pulse study.** The low-pass-filter test shows that very short pulses stay nearly unconditional but suffer strong amplitude loss if the filter is not inverted. A lightweight precompensation pass would clarify whether "make it fast" remains practical under hardware filtering.

- **[P3 | LOW] Add branch-resolved phase-space movies.** The saved static figures are enough for the report, but short trajectory animations would make the conditionality mechanism immediately obvious in group discussions.

## Open Questions

- **Why does the short two-tone protocol compete so well with optimal control on vacuum and qubit-superposition tests?** This may indicate that most of the unconditional-displacement problem is captured by branch-frequency compensation alone at the present `chi`, `chi'`, and `K`.

- **Would a materially richer hardware-aware multiplex parameterization ever beat the short two-tone pulse?** The direct full-duration multicarrier fit, the segmented branch-resonant family, and the low-parameter jointly optimized shaped-two-tone family all failed. The remaining open question is therefore much narrower and higher cost than it appeared after the first multiplex follow-up.

- **Would a larger optimal-control basis beat two-tone compensation decisively, or would it mainly learn a slightly shaped two-tone solution?** On the explicit 14-state benchmark used in the follow-up, the bounded 40 ns optimal waveform no longer beats the short 20 ns two-tone pulse. The remaining question is therefore whether a substantially richer optimal-control run can reclaim a real advantage.

- **Can a filter-precompensated short Gaussian or cosine pulse recover the target amplitude without giving back the unconditionality gained from large bandwidth?** The current low-pass test only examined the un-precompensated case.

- **At what coherent amplitude does vacuum-matched two-tone calibration become insufficient even at short duration?** The selected protocol comparison hints at this issue, but the boundary has not been mapped.

## What Was Tried and Did Not Work

- **[NEGATIVE RESULT] A direct full-duration multiplex compression of the best constrained waveform did not work.** We projected the best `40 ns` optimal-control waveform onto `K = 2, 3, 4, 5, 6, 8` explicit carriers over the full pulse duration and re-simulated the resulting multicarrier drives. The best tested case was the `8`-tone approximation, but it still achieved only mean state fidelity `0.524`, minimum fidelity `1.89e-3`, vacuum branch mismatch `0.770`, and vacuum superposition fidelity `0.589`. Extra carriers without temporal segmentation are therefore not a viable substitute for the sampled waveform used by the constrained optimizer.

- **[NEGATIVE RESULT] Structured low-dimensional multiplex refinements also failed.** We tested a segmented branch-resonant family using the two branch carriers repeated over `1, 2, 4, 8` time windows and a four-segment shaped-two-tone family in which only one complex scale per segment was jointly optimized. The best segmented case (`8` segments) still reached only mean fidelity `0.361`, minimum fidelity `8.14e-3`, and vacuum `F_{+x} = 0.657`. The jointly optimized shaped-two-tone case converged successfully after `375` function evaluations, but still reached only mean fidelity `0.158`, minimum fidelity `0.0459`, and vacuum `F_{+x} = 0.228`. The structured-multiplex gap is therefore not just a no-segmentation artifact.

- **[NEGATIVE RESULT] A practical echoed displacement did not solve the unconditional-control problem.** The tested sequence `D(alpha/2) -> X_pi -> D(alpha/2) -> X_pi`, implemented with a vacuum-calibrated Gaussian-DRAG `pi` pulse of `20 ns`, never beat the best simple protocols on the main superposition benchmark. Its best `|+x>` result occurred at total duration `60 ns`, where `delta_alpha = 0.0787` and the vacuum superposition fidelity was only `0.891`, far worse than the short two-tone and fast-pulse options. Root cause: the `pi` pulse itself becomes Fock-manifold dependent once the cavity is populated, so the toggling-frame cancellation picture breaks at the waveform level.

- **[NEGATIVE RESULT] Short pulses plus filtering are not automatically usable.** A `5 ns` Gaussian pulse remains nearly branch-matched after a `40 MHz` low-pass filter (`delta_alpha = 0.0236`), but the filter removes so much pulse area that both branch fidelities collapse to about `0.782`. The lesson is that bandwidth limitation must be accompanied by precompensation or recalibration; otherwise the pulse is "unconditional" only because it no longer performs the intended displacement.

- **[NEGATIVE RESULT] Long naive pulses fail exactly in the way the simple dispersive picture predicts.** For the baseline square pulse at `alpha = 1`, the superposition benchmark degrades from `0.994` at `5 ns` to `0.770` at `80 ns` and `0.563` at `160 ns`, while entanglement grows from `8.5e-3` bits to `0.483` and then `0.792` bits. This confirms that `|chi| T`, not `chi'` or Kerr, is the main driver of unconditionality failure here.

## Compute & Resource Notes

- **Main unconditional-displacement driver:** about `10.1 s` total wall-clock time, including sweeps, constrained optimal control, artifact writes, and figure generation.
- **Single-tone sweep block:** about `3.07 s`.
- **Two-tone sweep block:** about `0.39 s`.
- **Echo sweep block:** about `0.64 s`.
- **Hardware-aware optimal-control block:** about `1.57 s`.
- **Memory footprint:** comfortably below the earlier `~200 MB` convergence-study peak from the broader waveform study because the new workflow reuses the same compact dispersive model and keeps the optimal-control Hilbert space modest.

## Historical Carry-Over From The Earlier Broader Waveform Study

- **[P2 | LOW] Sweep `chi` more broadly across device designs.** The current unconditional-displacement sweep varies `chi` by scale factors, but a future generalized regime map should still cover the full experimentally relevant range.

- **[P2 | MEDIUM] Include thermal cavity occupation.** The present study starts from pure low-Fock and coherent states. A thermal initial state is still relevant for hardware realism.

- **[P2 | MEDIUM] Study multi-pulse sequence accumulation.** The current work characterizes single unconditional-displacement primitives, but practical bosonic protocols chain displacements with other operations and may accumulate additional conditional phases.

- **[P3 | MEDIUM] Add logical-subspace process tomography.** The present study emphasizes explicit state-test fidelity, branch mismatch, and entanglement. Full logical-subspace process reconstruction would be a stronger gate-level diagnostic.

- **Historical open question retained:** the transition from broadband control to intentionally number-selective control remains an interesting follow-up because the same `|chi| T` scaling underlies both regimes.

## Resolved

- **Broad-state protocol ranking on the explicit 14-state set:** Resolved by the multiplex follow-up benchmark. Once the same broad metric is evaluated for the calibrated `20 ns` and `40 ns` two-tone pulses, the `20 ns` two-tone case is the strongest tested protocol (`mean fidelity = 0.9857`, `minimum fidelity = 0.9242`), outperforming the bounded `40 ns` optimal-control benchmark on this specific score.

- **Low-dimensional structured multiplex follow-up:** Resolved by the extended benchmark. The full-duration multicarrier fit, segmented branch-resonant family, and four-segment jointly optimized shaped-two-tone family all remain far below the short `20 ns` two-tone pulse and the bounded sampled waveform on the broad state metric, so the remaining multiplex question is now restricted to materially richer hardware-aware parameterizations.

- **Displacement pulse phase convention:** Fixed by computing the pulse phase from `epsilon = i alpha / T` instead of hard-coding zero phase.

- **Session-variable reuse in an earlier validation script:** Fixed by creating a dedicated representative session rather than reusing the last loop variable.

- **Optimal-control initial-schedule shape mismatch:** Fixed by constructing the held-sample schedule from the actual `sample_period_s` and duration before calling the GRAPE solver.
