# Execution Summary

## Current State

- Study scaffold created at studies/hybrid_unitary_native_entangling_evolution.
- Archive-derived Phase 1 bootstrap completed with JSON and figure outputs.
- Phase 2 native-block composition search completed and written to JSON plus figure outputs.
- Phase 4 symbolic depth diagnostics completed and written to JSON plus figure outputs.
- Phase 5 runtime validation completed with replayable surrogate locals, truncation sweeps at `n_cav = 10, 12, 14`, nominal-noise replay at `n_cav = 12`, and saved runtime artifacts.
- report/report.tex, references.bib, and report/report.pdf now exist for the runtime-validated follow-up milestone.

## Available Artifacts

- studies/hybrid_unitary_native_entangling_evolution/README.md
- studies/hybrid_unitary_native_entangling_evolution/IMPROVEMENTS.md
- studies/hybrid_unitary_native_entangling_evolution/data/phase1_candidate_bootstrap.json
- studies/hybrid_unitary_native_entangling_evolution/figures/phase1_candidate_bootstrap.pdf
- studies/hybrid_unitary_native_entangling_evolution/data/phase2_native_block_search.json
- studies/hybrid_unitary_native_entangling_evolution/figures/phase2_native_block_search.pdf
- studies/hybrid_unitary_native_entangling_evolution/data/phase3_replay_support_check.json
- studies/hybrid_unitary_native_entangling_evolution/data/phase4_depth_diagnostics.json
- studies/hybrid_unitary_native_entangling_evolution/data/phase5_candidate_comparison.csv
- studies/hybrid_unitary_native_entangling_evolution/data/phase5_convergence_table.csv
- studies/hybrid_unitary_native_entangling_evolution/data/phase5_runtime_metrics.csv
- studies/hybrid_unitary_native_entangling_evolution/data/phase5_symbolic_metrics.csv
- studies/hybrid_unitary_native_entangling_evolution/data/phase5_runtime_validation.json
- studies/hybrid_unitary_native_entangling_evolution/data/phase5_weight_sensitivity.csv
- studies/hybrid_unitary_native_entangling_evolution/figures/phase4_probe_fidelity_summary.pdf
- studies/hybrid_unitary_native_entangling_evolution/figures/phase4_bloch_N2_A_local_to_A_local_qx_plus_0.pdf
- studies/hybrid_unitary_native_entangling_evolution/figures/phase4_wigner_N2_A_local_to_A_local_g0.pdf
- studies/hybrid_unitary_native_entangling_evolution/figures/phase5_candidate_comparison.pdf
- studies/hybrid_unitary_native_entangling_evolution/figures/phase5_convergence.pdf
- studies/hybrid_unitary_native_entangling_evolution/figures/phase5_weight_sensitivity.pdf
- studies/hybrid_unitary_native_entangling_evolution/figures/phase5_bloch_R2_exact_runtime_to_exact_runtime_qx_plus_0.pdf
- studies/hybrid_unitary_native_entangling_evolution/figures/phase5_bloch_R2_A_runtime_to_A_runtime_qx_plus_0.pdf
- studies/hybrid_unitary_native_entangling_evolution/figures/phase5_wigner_compare_R2_exact_runtime_to_exact_runtime_g0.pdf
- studies/hybrid_unitary_native_entangling_evolution/figures/phase5_wigner_compare_R2_A_runtime_to_A_runtime_g0.pdf
- studies/hybrid_unitary_native_entangling_evolution/artifacts/phase5_local_surrogate_exact_hc.json
- studies/hybrid_unitary_native_entangling_evolution/artifacts/phase5_local_surrogate_exact_hc_grape.json
- studies/hybrid_unitary_native_entangling_evolution/artifacts/phase5_local_surrogate_A_local.json
- studies/hybrid_unitary_native_entangling_evolution/artifacts/phase5_local_surrogate_A_local_grape.json
- studies/hybrid_unitary_native_entangling_evolution/artifacts/phase5_runtime_candidate_R2_exact_runtime_to_exact_runtime.json
- studies/hybrid_unitary_native_entangling_evolution/artifacts/phase5_runtime_candidate_R2_A_runtime_to_A_runtime.json
- studies/hybrid_unitary_native_entangling_evolution/artifacts/phase5_runtime_candidate_R2_B_local_replay.json
- studies/hybrid_unitary_native_entangling_evolution/report/report.tex
- studies/hybrid_unitary_native_entangling_evolution/report/report.pdf

## Phase 2 Results

- The two-native-wait architecture is the dominant native-heavy family for the full target.
- Best lower-bound candidate: `N2_exact_hc_to_exact_hc`
	- ideal fidelity `~1.000000`
	- entangling gate count `2`
	- total entangling time `352 ns`
	- depth `11`
- Best archive-local candidate: `N2_exact_hc_to_A_local`
	- ideal fidelity `~0.988734`
	- entangling gate count `2`
	- total entangling time `352 ns`
	- depth `13`
- Best fully archive-driven native candidate: `N2_A_local_to_A_local`
	- ideal fidelity `~0.985736`
	- entangling gate count `2`
	- total entangling time `352 ns`
	- depth `17`
- One native entangler is not enough for this target: the best one-wait candidate saturates at fidelity `0.25`.
- `B_local`-based variants are strongly disfavored because they add substantial local complexity and entangling exposure without competitive fidelity.

## Phase 3.1 Replay-Support Result

- Replay-support check completed in studies/hybrid_unitary_native_entangling_evolution/data/phase3_replay_support_check.json.
- All shortlisted two-wait native-heavy candidates fail direct `waveform_sequence_from_gates(...)` conversion.
- Immediate blocker: the current candidate constructions use exact `PrimitiveGate` local Hadamard references, and the installed bridge rejects `PrimitiveGate`.
- Bridge support message from the installed build: `QubitRotation`, `Displacement`, `SQR`, and `ConditionalPhaseSQR` are supported; `SNAP` and `FreeEvolveCondPhase` should go through a model-backed path rather than the waveform bridge.

## Phase 5 Runtime-Validation Results

- The symbolic-only blocker is resolved. The follow-up runtime path uses:
	- a replayable qubit-H replacement: `QubitRotation(theta = pi/2, phi = pi/2)` followed by `QubitRotation(theta = pi, phi = 0)`
	- explicit zero-amplitude idle pulses for `FreeEvolveCondPhase`
	- GRAPE-derived replayable surrogates for the non-replayable `exact_hc` and `A_local` local blocks
- Surrogate quality at reference truncation `n_cav = 12`, `n_tr = 3`:
	- `exact_hc_runtime`: nominal fidelity `0.94288`, total duration `1.28 us`, `160` pulses
	- `A_local_runtime`: nominal fidelity `0.94762`, total duration `1.28 us`, `160` pulses
	- all GRAPE restarts stopped at the iteration limit, so these are usable but not converged local surrogates
- Runtime shortlist at `n_cav = 12`, `n_tr = 3`:
	- `R2_exact_runtime_to_exact_runtime`
		- process fidelity `0.90269`
		- average probe fidelity `0.83495`
		- leakage `0.04318`
		- nominal-noise average probe fidelity `0.61448`
		- Wigner overlap `0.96384`
		- total runtime `4.432 us`
	- `R2_A_runtime_to_A_runtime`
		- process fidelity `0.90154`
		- average probe fidelity `0.81219`
		- leakage `0.04725`
		- nominal-noise average probe fidelity `0.58042`
		- Wigner overlap `0.98142`
		- total runtime `4.432 us`
	- `R2_B_local_replay`
		- process fidelity `0.16665`
		- average probe fidelity `0.18729`
		- leakage `0.14324`
		- nominal-noise average probe fidelity `0.04885`
		- total runtime `7.912 us`
	- `R1_exact_runtime_to_exact_runtime`
		- process fidelity `0.22340`
		- average probe fidelity `0.21197`
		- leakage `0.02860`
		- nominal-noise average probe fidelity `0.04801`
		- total runtime `2.856 us`
- Ranking shift:
	- Best experimentally grounded symbolic candidate: `N2_A_local_to_A_local`
	- Best replay-backed closed-system candidate: `R2_exact_runtime_to_exact_runtime`
- Truncation stability across `n_cav = 10, 12, 14` is strong. The process-fidelity span is `7.67e-06` for `R2_exact_runtime_to_exact_runtime` and `2.79e-05` for `R2_A_runtime_to_A_runtime`, so the ranking shift is not a truncation artifact.
- Cost-function warning: if the infidelity term is weakened by `25%`, the one-wait replay baseline can become the lowest scalar cost despite failing the task. Any future ranking should enforce a fidelity floor or use Pareto reporting.

## Phase 3.2 / 3.3 Depth Diagnostics

- Added studies/hybrid_unitary_native_entangling_evolution/scripts/phase4_depth_diagnostics.py and generated checkpointed cqed_sim depth diagnostics on the symbolic gate models.
- The one-wait baseline `N1_exact_hc_to_exact_hc` remains uniformly bad across the logical probe set:
	- average final probe fidelity `~0.250000`
	- final Bloch-probe target overlap `~0.250000`
	- final Wigner-probe target overlap `~0.250000`
- The exact two-wait upper bound `N2_exact_hc_to_exact_hc` is structurally faithful throughout the probe set:
	- average final probe fidelity `~0.9999997`
	- final Bloch-probe target overlap `~0.9999993`
	- final Wigner-probe target overlap `~0.99999999`
- The fully archive-driven candidate `N2_A_local_to_A_local` is now the best experimentally grounded symbolic recommendation:
	- average final probe fidelity `~0.9692`
	- final Bloch-probe target overlap `~0.9539`
	- final Wigner-probe target overlap `~0.9694`
	- maximum Wigner negativity on the `g0` probe track `~0.0110`
- Interpretation: the `A_local`-based two-wait candidate preserves the native-entangler architecture and recovers the target well at the end, but it shows noticeably larger transient Bloch and Wigner distortions than the exact two-wait upper bound.

## Notes for Resume

- The bootstrap already normalized legacy candidates onto one bookkeeping layer.
- The requested depth diagnostics and the runtime-validation path are both complete.
- The report has been updated to distinguish symbolic upper bounds from replay-backed conclusions and to state explicitly that the follow-up is architecture-validating rather than hardware-ready.
- Any future extension work should start from the current runtime winner `R2_exact_runtime_to_exact_runtime` and focus on reducing local-surrogate duration and improving nominal-noise performance rather than reopening the missing-replay-path question.