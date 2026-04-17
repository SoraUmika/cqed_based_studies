# Improvement Log: Holographic Cluster-State Synthesis Unified Study

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- [P1 | MEDIUM] The next restricted-target search should start directly at the 2- and 3-block frontiers rather than warm-starting from the best full-target representatives. The current ground-sector follow-up already reaches `F_ground,12 > 0.99` at 3 blocks in both families and leaves the 2-block frontier just below threshold (`0.987565` for SQR, `0.988681` for CPSQR), so a fresh low-depth search could still change the minimum-depth conclusion.
- [P1 | HIGH] The new ground-sector winners are still closed-system results only. The current physically plausible winners are very clean under ancilla/support/truncation checks, but no open-system validation has yet been performed on the restricted-target branch, so the shallow 3-block ranking could still reorder under realistic decoherence.

## Recommended Improvements (P2)
- [P2 | MEDIUM] Re-run the strongest restricted-target structures directly at `N_cav = 12` from scratch within each candidate block depth, rather than using only the warm-started refinement from the best full-target representatives, to check whether shallower or shorter solutions exist.
- [P2 | MEDIUM] Repeat the restricted-target local refinement with multiple seeds and a larger optimizer budget. The current ground-sector pass used one seed and `maxiter = 16` per block representative; the numerical outcomes are stable across `N_cav = 10, 12, 14`, but local-minimum risk is not fully ruled out.
- [P2 | MEDIUM] After the minimum block count is fixed, expand the finalist-only level-subset refinement around the winning orderings (`RDS` for SQR, `DRCP` for CPSQR) to reduce `n_active` without losing the `99%` threshold. Tone-count reduction should remain a follow-up compression step, not the primary axis.
- [P2 | MEDIUM] Add a retained-candidate GRAPE Wigner column if a clean waveform-to-state pipeline is exposed locally, so the structured-family Wigner comparison can be anchored against the validated pulse reference.
- [P2 | MEDIUM] Re-run the dedicated SQR ordering sweep with a second refinement budget or seed set to confirm that the fixed-budget ranking `RDS > DSR > SDR > SRD > RSD > DRS` is not a one-seed artifact. The present result strongly disfavors the original `DRS` ordering, but only one local refinement budget has been checked.
- [P2 | MEDIUM] Test replay-aware or noise-aware GRAPE objectives for the full target unitary, preferably with the JAX-enabled frontier runner. The extended `N_cav = 12` frontier through `500, 600, 800, 1000 ns` did not yet reach replay fidelity `0.99` under the explored durations and `maxiter = 250` budget; the best tested point is `800 ns` with replay `0.9760` and open-process fidelity `0.9187`, and every seed stops at the fixed iteration cap.
- [P2 | MEDIUM] Rework the primitive-level GRAPE objective before treating a hybrid structured-plus-GRAPE route as viable. The current representative multiseed runs for the retained `S0` and `CP2` blocks improve the cavity marginals more than the full joint state, but they still stop at `F_sub = 0.6959` (SQR) and `0.4754` (CPSQR), so the present primitive diagnostics are informative rather than deployment-ready.
- [P2 | LOW] Promote the active-subspace diagnostics and the ordered-family builder helpers into reusable `cqed_sim` or study-common utilities; they are now used by the comprehensive design-space runner and no longer belong only to the earlier corrected-scope rerun.

## Nice-to-Haves (P3)
- [P3 | LOW] Test whether CPSQR orderings with the strongest median screen performance (`CPDR`) can be made competitive with the best final `DRCP` solution under a targeted local refinement.

## Open Questions
- Can either family push the 2-block restricted-target frontier above `0.99` once the search is run directly for the reduced objective rather than warm-started from the old full-target winners?
- Is the 3-block restricted-target winner still the best practical option once open-system validation is included, or do the deeper but even cleaner 4- and 5-block solutions win on hardware?
- Why does the restricted objective favor solutions with very poor full-target fidelity (`~0.17` to `0.41`) even after local refinement, rather than preserving partial compatibility with the original joint logical map?
- Why does the SQR family switch from the best median screen ordering `DSR` to the best refined winner `RDS` once the finalists are re-optimized at `N_cav = 12`?
- Why does the best CPSQR solution use only two active tones on levels `{1,2}` even though the strongest coarse screen replay came from four-tone `RDCP` or `DCPR` structures?
- Does the very long but very clean SQR winner remain preferable once realistic decoherence is included, or does the shorter CPSQR solution become the better hardware-level compromise?
- The residual Wigner mismatch does not appear to come from unconstrained spectator action alone: why do spectator-constrained re-fits degrade both families instead of improving phase-space agreement?
- Is the current full-target GRAPE shortfall to replay fidelity `0.99` mainly a consequence of the closed-system objective and fixed `maxiter = 250` seed budget, or does it persist after replay-aware/noise-aware objectives and JAX-accelerated reruns?

## What Was Tried and Did Not Work
- Simply reusing the old full-target winners under the restricted ground-sector metric does not work. The saved `5`-block SQR and CPSQR winners rescored at `N_cav = 12` land near `F_ground ≈ 0.50` because they excite the ancilla almost completely (`ancilla_worst ≈ 0.998` to `0.9999`), so the reduced objective must be optimized explicitly.
- The original unified report inherited mixed-family conclusions from upstream studies. That consolidation is not valid for the corrected scope because it compares families the user explicitly excluded.
- A naive exhaustive screen over all ordered level subsets for every structural setting expanded the search to 810 cases and was too slow to complete reliably. The workable approach was the staged pipeline now in `scripts/run_design_space_study.py`: 216 structural cases, then 30 finalist level-subset refinements, then 12 physical finalists.
- The earlier corrected-scope conclusion that only a low-confidence CPSQR family survived was incomplete. Once `n_active`, block count, and ordering were searched systematically and the strongest structures were re-optimized at `N_cav = 12`, both exclusive families produced retained `>99%` structured solutions.
- Enlarged spectator-constrained Wigner re-fits did not support the hypothesis that residual Wigner mismatch is caused mainly by unconstrained spectator evolution. For `D + R + SQR`, the deeper stable rerun (`logical=4`, `augmented=4`, non-fast path) kept the logical refit excellent (`F12 = 0.999641`, mean Wigner RMS `0.00204`) but the augmented refit still degraded badly (`F12 = 0.887459`, mean Wigner RMS `0.01883`).
- `D + R + CPSQR` is numerically less stable under those deeper basis-extension refits on the current Windows machine. The deepest complete local result uses a mixed backend path: logical refit at `maxiter=3` with `use_fast_path=False`, augmented refit at `maxiter=2` with `use_fast_path=True`. That run still matches the earlier negative conclusion (`logical mean Wigner RMS = 0.00470`, augmented mean Wigner RMS = 0.02824`, interpretation `does_not_support`). Attempts to push the CPSQR logical refit to `maxiter=4` or the spectator-constrained augmented refit beyond `maxiter=2` terminated without a Python traceback.
- The `run_structured_extension_analysis.py` fast path can terminate silently on this Windows setup during deeper logical refits. Adding incremental checkpoints plus a `--disable-fast-path` option allowed the `D + R + SQR` 4/4 basis study to complete and preserved partial progress for `D + R + CPSQR`, but did not eliminate the CPSQR augmented-refit instability.
- Representative primitive-level GRAPE optimization of the retained `S0` and `CP2` blocks (`N_cav: 8 -> 12`, seeds `17/42/73`, `maxiter = 200`, JAX engine) did not produce hybrid-ready primitives. `S0` tops out at `F_sub = 0.6959` with mean full-state / reduced-cavity probe fidelities `0.5793 / 0.8429`; `CP2` tops out at `F_sub = 0.4754` with mean full-state / reduced-cavity probe fidelities `0.2892 / 0.5948`. All six seed runs terminate with the L-BFGS-B `FACTR*EPSMCH` condition at zero recorded optimizer iterations, so simply raising `maxiter` is not enough under the current primitive setup.

## Compute & Resource Notes
- Comprehensive runner: `scripts/run_design_space_study.py`
- Ground-sector follow-up runner: `scripts/run_ground_sector_followup.py`
- Final search schedule: 216 structural screen cases at `N_cav = 4`, 30 finalist level-subset refinements, 12 physical `N_cav = 12` refinements, then validation at `N_cav = 10, 12, 14`.
- Restricted-target follow-up schedule: rescore 9 saved block representatives at `N_cav = 10, 12, 14`, then run 9 warm-started local restricted-target refinements at `N_cav = 12` with `maxiter = 16`, then re-check all refined representatives at `N_cav = 10, 12, 14`.
- The design-space run is only practical with multiprocessing. The finalized runner uses an 8-worker process pool; smaller worker counts make the 4- and 5-block families prohibitively slow.
- Figure-only regeneration is cheap via `python run_design_space_study.py --skip-search`.
- Notebook execution of the full restricted-target follow-up took about `163 s` end-to-end on this machine and wrote `data/ground_sector_followup_summary.json`, `data/ground_sector_followup_candidates.csv`, `artifacts/ground_sector_best_{sqr,cpsqr}.json`, and the two `ground_sector_*` figures.
- Fixed-budget SQR ordering sweep (`data/sqr_ordering_fixed_budget.csv`) took roughly `59` to `92` seconds per ordering at `maxiter=24`, with the best ordering `RDS` reaching `F12 = 0.998389` and the worst `DRS` falling to `F12 = 0.992561` with worst leakage `2.50e-2`.
- Extended full-target GRAPE frontier (`scripts/run_grape_frontier_extension.py`) at `N_cav = 12`, `maxiter = 250`, seeds `{17, 42, 73}`: `500 ns -> best replay 0.962800, open 0.911215`; `600 ns -> 0.964843, open 0.909132`; `800 ns -> 0.976048, open 0.918705`; `1000 ns -> 0.971178, open 0.901765`. Every seed terminated with `STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT`.
- JAX is available on this machine, and the frontier runner now exposes `--engine auto|numpy|jax` plus `--jax-device` so future GRAPE reruns can use the JAX solver path without modifying study code.
- Primitive-level GRAPE diagnostics (`scripts/run_primitive_grape_diagnostics.py`, seeds `17/42/73`, `maxiter = 200`, JAX engine) take about `21` to `22` seconds per full two-primitive pass on this machine and write `data/primitive_grape_diagnostics.json`, `artifacts/primitive_grape_diagnostics_arrays.npz`, `figures/primitive_grape_summary.pdf`, and `figures/primitive_grape_wigner.pdf`.

## Resolved
- Added reusable ground-sector helpers to `scripts/common.py`: restricted fidelity, ancilla/support leakage accounting, and sequence evaluation on the `{|g,0>, |g,1>}` subspace.
- Completed the restricted-target follow-up in `scripts/run_ground_sector_followup.py` and saved machine-readable summaries, figures, and per-family artifacts.
- Established that both exclusive structured families cross the restricted `0.99` threshold at 3 blocks after local ground-sector refinement.
- Established the current best absolute restricted-target candidates:
	- `D + R + SQR`: 5 blocks, `RDS`, `n_active = 4`, `F_ground,12 = 0.9999260`.
	- `D + R + CPSQR`: 4 blocks, `RDCP`, `n_active = 4`, `F_ground,12 = 0.9999284`.
- Re-ran the unified study under exclusive D + R + SQR and D + R + CPSQR families only.
- Enforced physical retention at N_cav = 10, 12, 14 with N_cav = 12 as the default final truncation.
- Added active-subspace tracking from sampled full-sequence evolution and reported the retained/discarded candidates accordingly.
- Moved the appendix Wigner figure to a single local figure directory and regenerated it from corrected-scope outputs.
- Performed the requested design-space sweep over `n_active`, block count, and gate ordering for both exclusive structured families.
- Identified retained `>99%` closed-system structured solutions in both families:
	- `D + R + SQR`: 5 blocks, `n_active = 4`, ordering `RDS`, levels `{0,1,2,3}`, `F_{12} = 0.9996347`.
	- `D + R + CPSQR`: 5 blocks, `n_active = 2`, ordering `DRCP`, levels `{1,2}`, `F_{12} = 0.9974917`.
