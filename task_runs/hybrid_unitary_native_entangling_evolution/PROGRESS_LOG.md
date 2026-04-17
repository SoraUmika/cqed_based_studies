# Progress Log

## 2026-03-24T00:00:00Z - Study initialized and archive bootstrap completed
- Objective: determine the best physically realizable decomposition of the established hybrid U_target when native chi-generated entangling evolution is treated as the dominant expensive resource.
- Created the new study directory at studies/hybrid_unitary_native_entangling_evolution.
- Added README.md and IMPROVEMENTS.md defining the entangling-weighted scope, metrics, outputs, and known risks.
- Added shared utilities and a Windows runtime-compat shim for future cqed_sim work.
- Added scripts/phase1_candidate_bootstrap.py and generated the first machine-readable candidate frontier plus figure from the earlier hybrid study archive.
- Result of the bootstrap: native baselines remain strongest on cost, B_ent remains the best selective physically replayed archive candidate, and L1d remains the key example where ideal structure diverges from physical viability.
- Next: add a Phase 2 native-entangler-biased search over full U_target decompositions using the optimized local and entangling building blocks from the legacy study.

## 2026-03-24T01:00:00Z - Phase 2 native-block search completed
- Added studies/hybrid_unitary_native_entangling_evolution/scripts/phase2_native_block_search.py.
- The new script reuses optimized legacy blocks (`A_local`, `B_local`, `D_ent`) and composes full-U_target candidates with one or two native chi-wait entanglers.
- Generated studies/hybrid_unitary_native_entangling_evolution/data/phase2_native_block_search.json and matching PNG/PDF figures.
- Main scientific result: two native entanglers are enough to recover the target at very high ideal fidelity, while one native entangler is not. The lower-bound candidate `N2_exact_hc_to_exact_hc` reaches fidelity ~1.0 with two waits and 352 ns total entangling time.
- Practical local-control result: replacing the ideal cavity-H with the archive `A_local` block still leaves the two-wait construction near 0.99 ideal fidelity (`N2_exact_hc_to_A_local`), while `B_local`-based constructions are heavily penalized by duration and reduced fidelity.
- Plan change: the study is now ready to shift from Phase 2 implementation into validation, beginning with replay-support verification and depth-indexed diagnostics on the top two native-heavy candidates.

## 2026-03-24T02:00:00Z - Replay-support validation completed
- Added studies/hybrid_unitary_native_entangling_evolution/scripts/phase3_replay_support_check.py and generated studies/hybrid_unitary_native_entangling_evolution/data/phase3_replay_support_check.json.
- The shortlisted native-heavy candidates (`N2_exact_hc_to_exact_hc`, `N2_exact_hc_to_A_local`, `N2_A_local_to_A_local`) all fail direct waveform-bridge conversion in the installed cqed_sim build.
- The blocking reason is explicit: the current candidate constructions use `PrimitiveGate` local Hadamard references, and `waveform_sequence_from_gates(...)` rejects `PrimitiveGate`. The installed bridge reports support for `QubitRotation`, `Displacement`, `SQR`, and `ConditionalPhaseSQR` only, and routes `SNAP` / `FreeEvolveCondPhase` sequences toward the model-backed path.
- Consequence: the next validation step must either replace the exact local-H references with replayable local gates or shift the shortlisted native-heavy candidates onto a model-backed simulation path that does not require waveform-bridge conversion.

## 2026-03-24T03:00:00Z - Phase 4 symbolic depth diagnostics completed
- Added studies/hybrid_unitary_native_entangling_evolution/scripts/phase4_depth_diagnostics.py.
- Generated studies/hybrid_unitary_native_entangling_evolution/data/phase4_depth_diagnostics.json together with Bloch-versus-depth, Wigner-versus-depth, and final-probe-fidelity summary figures.
- The checkpointed cqed_sim gate-model diagnostics cleanly separate the three main architectural regimes:
	- `N1_exact_hc_to_exact_hc` remains a hard one-wait lower bound with probe fidelities pinned at ~0.25.
	- `N2_exact_hc_to_exact_hc` remains the symbolic two-wait upper bound with probe fidelities effectively at 1.0.
	- `N2_A_local_to_A_local` remains the best experimentally grounded symbolic candidate, with average final probe fidelity ~0.969 and visible mid-sequence Bloch/Wigner excursions before recovery at the end.
- The study is no longer missing the requested depth diagnostics. The remaining technical blocker is specifically pulse-backed replay for the native-heavy candidates, not diagnostic coverage.

## 2026-03-24T04:00:00Z - Report drafted and compiled
- Added studies/hybrid_unitary_native_entangling_evolution/report/report.tex and references.bib.
- Compiled studies/hybrid_unitary_native_entangling_evolution/report/report.pdf using the direct `pdflatex` / `bibtex` chain because `latexmk` is unusable on this machine without `perl`.
- The report captures the final symbolic study recommendation: the two-native-wait architecture is the correct family, the exact two-wait candidate is the symbolic upper bound, and `N2_A_local_to_A_local` is the best experimentally grounded symbolic focal candidate.
- The open pulse-backed replay issue is now documented as future work rather than an unrecorded gap.

## 2026-03-25T00:00:00Z - Phase 5 runtime validation completed
- Added studies/hybrid_unitary_native_entangling_evolution/scripts/phase5_runtime_validation.py and patched the shared study utilities so the local `cQED_simulation` checkout is resolved dynamically on this machine.
- Implemented a replayable runtime path for the shortlisted native-heavy families by replacing the exact qubit Hadamard with a `QubitRotation` pair, replaying `FreeEvolveCondPhase` through explicit idle evolution, and replacing non-replayable `exact_hc` / `A_local` locals with GRAPE-derived replayable surrogates.
- Generated machine-readable runtime artifacts in studies/hybrid_unitary_native_entangling_evolution/artifacts together with studies/hybrid_unitary_native_entangling_evolution/data/phase5_runtime_validation.json, comparison/convergence CSV tables, and Phase 5 comparison, convergence, sensitivity, Bloch, and Wigner figures.
- Main runtime result: the symbolic recommendation does not survive unchanged. `N2_A_local_to_A_local` was the best experimentally grounded symbolic candidate, but the best replay-backed closed-system candidate is `R2_exact_runtime_to_exact_runtime` with process fidelity `0.90269`, average probe fidelity `0.83495`, leakage `0.04318`, and nominal-noise average probe fidelity `0.61448` at `n_cav = 12`, `n_tr = 3`.
- The replay-backed archive-local surrogate `R2_A_runtime_to_A_runtime` remains close in process fidelity (`0.90154`) and has the better final Wigner overlap (`0.98142` versus `0.96384`), but it is worse in average probe fidelity (`0.81219`) and nominal-noise average probe fidelity (`0.58042`).
- The direct replayable archive reference `R2_B_local_replay` remains strongly disfavored: process fidelity `0.16665`, average probe fidelity `0.18729`, leakage `0.14324`, runtime duration `7.912 us`, and nominal-noise average probe fidelity `0.04885`.
- Truncation stability checks across `n_cav = 10, 12, 14` show very small metric spans for all runtime candidates, so the ranking shift is not a truncation artifact.
- New caution: the scalar weighted-cost ranking can incorrectly prefer the one-wait replay baseline if the infidelity term is weakened. Raw metrics and a fidelity floor must remain primary in any future recommendation.