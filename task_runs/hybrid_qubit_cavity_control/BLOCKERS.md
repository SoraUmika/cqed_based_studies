# Blockers

## Active Blockers
(none)

## Resolved Blockers
- **Random SU(4) model-based synthesis**: Previously reported as `_expand_target_matrix` shape mismatch. **RESOLVED (Iteration 7)**: Tested model-based `UnitarySynthesizer` with 4×4 target + Subspace.custom(full_dim=16, indices=(0,1,8,9)) — both ideal and model-based synthesis succeed (objective=0.530). The original (8,8) error was likely from an earlier framework version or incorrect target matrix sizing.
- **P1.2 n_tr=3 inner loop**: `targeted_subspace_multitone` raised `ValueError: bloch_xyz_from_qubit_state requires a two-level reduced state` when `n_tr=3`. **RESOLVED (Iteration 5)**: Added `truncate_to_qubit_subspace()` upstream to `cqed_sim/sim/extractors.py`; `bloch_xyz_from_qubit_state` now accepts dim ≥ 2 with automatic projection to {|g>,|e>} subspace.
- **P1.4 ConditionalPhaseSQR waveform bridge**: `waveform_sequence_from_gates` did not support ConditionalPhaseSQR gates. L2d sequences could not be replayed. **RESOLVED (Iteration 6)**: Added `ConditionalPhaseSQRGate` to `cqed_sim/io/gates.py`, added CPSQR handler to `cqed_sim/unitary_synthesis/waveform_bridge.py` using SQR multitone hardware with theta=0 (zero drive; conditional phases from dispersive drift).
- **Fresh noisy SNAP rerun instability**: a direct rerun inside `run_speed_limit_feasibility.py` failed with `ValueError: Function to integrate must not return a tuple.` Resolution: reuse the already validated logical SNAP payload from the earlier follow-up study rather than drop the local primitive from the speed-limit comparison.
- **Targeted multitone qutrit mismatch**: running the targeted-subspace multitone path with `n_tr=3` failed because the reduced-qubit analysis expects a two-level state. Resolution: run the physical frontier at `n_tr=2` and record the missing leakage channel as a study limitation.
- **Optimizer-heavy gate sweep too slow**: a full per-point `optimize_targeted_subspace_multitone(...)` sweep was not tractable for the complete gate set. Resolution: switch to a direct validation scan over duration and Gaussian width, then carry the resulting frontier into the strategy ranking.
- **Missing fair local SQR physical route**: the first pass compared `B_local` only at the ideal level. Resolution: load the serialized optimized `B_local` gate list from the gate-set archive and compile both local selective blocks in the second pass.
- **No direct check of hidden `|f>` leakage**: the first pass ended at `n_tr=2`. Resolution: replay the refined waveforms on an `n_tr=3` model and record the maximum `|f>` population for every refined gate.
- **Parallel duplicate script runs during debugging**: an attached debug run and a background run were accidentally started together. Resolution: kill the background terminal, allow the attached copy to finish, and keep the final artifacts from the single surviving run.
- **GRAPE target matrix shape mismatch**: initial GRAPE implementation used 8×8 target matrix (n_cav dimension) instead of 4×4 (subspace dimension). Resolution: use `Subspace.custom(full_dim=16, indices=(0,1,8,9))` ordering for the 4×4 target.
- **TargetUnitary API**: `TargetUnitary(operator=matrix)` failed — API requires positional arg `TargetUnitary(matrix, ignore_global_phase=True)`.
- **UnitarySynthesizer API**: `UnitarySynthesizer(sequence=, n_cav=)` + `.run()` didn't work — correct pattern is `primitives=`, `subspace=`, `.fit()`.
- **Gate constructor parameter names**: `QubitRotation(angle=, axis_angle=)` → `QubitRotation(theta=, phi=)`. `Displacement(drift_model=)` doesn't accept drift_model.
- **Terminal buffer pollution**: foreground terminal accumulated ~16KB of old traceback output, preventing clean command output. Resolution: used background terminals for clean output.