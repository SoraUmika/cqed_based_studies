# Improvement Log: Cluster-State Unitary Decomposition with Native Free Evolution and SNAP

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- **[P1 | HIGH] Larger-truncation search coverage is still narrow**: The `2026-03-25` improvement pass re-optimized only the two depth-6 winners at `n_cav=6`. That is enough to confirm the native-vs-SNAP ranking under a stronger truncation, but it does not prove that the `n_cav=4` threshold frontiers remain globally optimal once all families and block counts are re-searched at `n_cav>=6`.

## Recommended Improvements (P2)
- **[P2 | HIGH] Replace the local Lindblad surrogate with a direct compiled noisy replay**: The `cqed_sim` waveform bridge on this machine rejects `FreeEvolveCondPhase` and `SNAP`, and the gate-sequence `pulse` backend ignores `NoiseSpec` during state propagation. The current noise result is therefore a documented surrogate rather than a native compiled replay. Upstream support for FE/SNAP waveform compilation, or a direct noisy Hamiltonian-slice runner for gate sequences, would close the remaining realism gap.
- **[P2 | LOW] Expand complexity metrics beyond gate counts**: The current report compares gate depth, total wait time, explicit SNAP count, and parameter count. A follow-up should add calibration burden, number of independently calibrated phase channels, and sensitivity to FE timing errors.
- **[P2 | LOW] Explore lower-depth SNAP ansaetze with stronger multistart coverage at larger truncation**: Interleaved SNAP clearly helps at the `0.95` fidelity threshold, but the current search does not prove that depth `18` is the true minimum, especially once the search space is moved from `n_cav=4` to `n_cav>=6`.

## Nice-to-Haves (P3)
- **[P3 | LOW] Test alternate cluster-target conventions**: The study optimized the default `make_target("cluster", n_match=1)` convention. A small follow-up could repeat the comparison for the equivalent `which="u2"` target to confirm that the phase-attribution story is convention-independent.
- **[P3 | LOW] Add a pulse-calibration-oriented summary table**: The present report focuses on logical performance and phase attribution. A one-page engineering appendix could summarize the number of displacement amplitudes, qubit rotation angles, FE waits, and SNAP phases that must be calibrated in each family.

## Open Questions
- Does interleaved SNAP remain the best family if the full block-count sweep is repeated at `n_cav=6` or `n_cav=8`, rather than only re-optimizing the original depth-6 winners?
- At larger truncation, is the remaining fidelity gap between the best native and best SNAP-assisted sequences still dominated by cavity-only phase correction, or by different leakage behavior?
- How much of the Lindblad-surrogate advantage for interleaved SNAP survives a future direct compiled noisy replay once FE/SNAP waveform support exists?

## What Was Tried and Did Not Work
- **Tail-only SNAP correction**: `snap_tail_b4` reached only `0.8050` fidelity, essentially matching the native four-block baseline. A single terminal SNAP layer does not compensate for an underpowered free-evolution scaffold.
- **Free-evolution removal ablations**: Removing all `FreeEvolveCondPhase` gates collapses the best native sequence to `0.1675` fidelity and the best FE+SNAP sequence to `0.3223`, so FE is the indispensable entangling resource.
- **Fixed-depth duration weighting as a route to lower FE wait**: Duration-weight frontier scans compressed total sequence duration, but they did not materially reduce the entangling wait budget of the best-depth native or interleaved SNAP winners. The practical wait-time improvement comes from switching families at a lower fidelity threshold (`0.95`), not from retuning duration weights at fixed depth.
- **Direct `NoiseSpec` on the gate-sequence `pulse` backend**: Passing `NoiseSpec` into `simulate_sequence(..., backend="pulse", state_inputs=...)` produced output states identical to the noiseless replay because the backend propagates the states gate-by-gate and does not invoke the noisy simulator. The improvement pass therefore switched to a local Lindblad surrogate built from the per-gate pulse unitaries.

## Compute & Resource Notes
- Full `run_study.py` screening + refinement + postprocessing took about `1280 s` wall-clock on the local workstation before a lightweight postprocess rerun completed the artifact serialization fix.
- The dominant cost was the ideal-gate optimization stage: coarse screening across all families plus two-seed refinement of each family winner. The frontier, ablation, and truncation checks were comparatively cheap.
- Screening settings: seed `17`, `20` optimizer iterations. Refinement settings: seeds `17` and `42`, `35` optimizer iterations. Truncation validation replayed the stored winners at `n_cav in {4,6,8}`.
- Best-case resource picture:
  native `native_fe_b6` -> fidelity `0.9639`, depth `20`, FE wait `1018 ns`, `48` continuous parameters.
  interleaved SNAP `snap_interleaved_b6` -> fidelity `0.9881`, depth `26`, FE wait `1027 ns`, `78` continuous parameters including `24` SNAP phases.
- Improvement-pass highlights (`data/improvement_pass_summary.json`):
  baseline replay at `n_cav=6` -> native `0.9177`, interleaved SNAP `0.9473`.
  re-optimized at `n_cav=6` -> native `0.9234`, interleaved SNAP `0.9908`.
  coherent replay at `n_cav=8` for the re-optimized winners -> native `0.9238`, interleaved SNAP `0.9867`.
  Lindblad-surrogate mean probe fidelity at `n_cav=8` -> native `0.8907`, interleaved SNAP `0.9351`.

## Resolved
- **Study initialization placeholders removed**: The improvement log now records the actual search failures, validation shortfalls, and compute budget observed during execution rather than generic startup notes.
- **[Resolved | 2026-03-25] Larger-truncation re-optimization for the two winning depth-6 families**: `run_improvement_pass.py` re-optimized `native_fe_b6` and `snap_interleaved_b6` at `n_cav=6`, lifting the native family from replay fidelity `0.9177` to `0.9234` and the interleaved SNAP family from `0.9473` to `0.9908`.
- **[Resolved | 2026-03-25] Coherent control-surrogate replay**: The same extension replayed the re-optimized winners with the `cqed_sim` pulse-unitary backend. The ranking survives at `n_cav=8`, with coherent replay fidelities `0.9238` (native) and `0.9867` (interleaved SNAP).
- **[Resolved | 2026-03-25] Nominal-noise ranking via local Lindblad surrogate**: Because the direct gate-sequence noisy replay path is not functional for FE/SNAP on this build, the study now includes a documented Lindblad surrogate using the per-gate pulse unitaries and nominal `T1`, `Tphi`, and cavity loss. That surrogate keeps the interleaved SNAP family ahead at both `n_cav=6` and `n_cav=8`.
