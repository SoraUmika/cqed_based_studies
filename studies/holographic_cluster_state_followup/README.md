# Holographic Cluster-State Control in cQED: Corrected Follow-up Study

## Problem Class
OPT | DES | ANA | REP

## Motivation
The prior `cluster_state_holographic_sim` report identified a promising GRAPE route for the cluster-state per-site holographic transfer unitary, but several decomposition-level conclusions were later shown to be artifacts of reduced Hilbert-space optimization and incomplete validation. This follow-up treats that earlier report as an intermediate draft rather than a final result.

The present study re-runs the comparison with stricter physical standards:
- reduced-subspace fidelity is necessary but not sufficient;
- embedded fidelity at larger cavity truncation matters;
- leakage outside the logical sector matters;
- cavity-state Wigner agreement matters;
- waveform-level replay matters;
- open-system performance matters.

The goal is to determine the most physically credible and experimentally relevant control strategy for implementing the cluster-state per-site holographic transfer unitary in the local dispersive transmon-cavity model encoded by `cqed_sim`.

## Goals
1. Recompute the sequential cluster-state transfer-channel spectrum and correct the transfer-matrix / correlation-length discussion so that the abstract, main text, and appendices all use one consistent convention.
2. Reassess the prior decomposition claims under strict validation by reporting ideal reduced fidelity, embedded fidelity, leakage, Wigner agreement, gate count, active time, active-tone count, and truncation stability for each candidate family.
3. Perform an explicit qubit-conditional follow-up for SQR / ConditionalPhaseSQR-like candidate families and determine whether their added expressivity survives embedding and phase-space validation.
4. Upgrade the GRAPE analysis from a promising draft result to a validated control result through a stronger multiseed duration study, independent waveform replay, truncation convergence, and open-system evaluation.
5. Produce a final recommendation that clearly separates validated conclusions, negative results, invalidated intermediate claims, and exploratory future directions.

## Methods
- `cqed_sim.unitary_synthesis.targets.make_target("cluster", n_match=1)` for the per-site target unitary.
- `cqed_sim.quantum_algorithms.HolographicChannel` for Kraus extraction and transfer-channel analysis.
- `cqed_sim.unitary_synthesis` primitives and `UnitarySynthesizer` for ideal-gate decomposition studies.
- `cqed_sim.unitary_synthesis.simulate_sequence` for embedded gate-sequence evaluation at enlarged cavity truncation.
- `cqed_sim.optimal_control.GrapeSolver` and `build_control_problem_from_model(...)` for direct waveform-level optimization.
- `ControlResult.to_pulses()`, `SequenceCompiler`, and `cqed_sim.sim.simulate_sequence(...)` for independent waveform replay.
- `cqed_sim.sim.noise.NoiseSpec` and control-evaluation replay for open-system analysis.
- `cqed_sim.sim.extractors.cavity_wigner` and reduced cavity-state extraction for phase-space validation.

## Expected Outcomes
1. A corrected transfer-channel interpretation with explicit eigenvalues and a clear statement of what the correlation-length formula does and does not mean for this sequential cluster-state channel.
2. A decomposition summary in which any route that fails embedding or Wigner validation is explicitly marked unsuccessful or incomplete.
3. A sharpened answer to whether qubit-conditional SQR-like structure closes the expressivity gap in a physically credible way.
4. A replay-validated GRAPE duration frontier that identifies the best current experimental candidate and quantifies the cost of decoherence.
5. A polished final report with reproducible artifacts, figures, and appendices.

## Known Limitations
- The coarse multiseed GRAPE sweep was performed at `N_cav = 8`; the final recommendation relies on a focused direct `N_cav = 12` rescue pass at `300 ns` and `400 ns` rather than a dense larger-truncation duration grid.
- The installed `cqed_sim` waveform bridge supports `QubitRotation`, `Displacement`, `SQR`, and `ConditionalPhaseSQR`, but not `SNAP` or `FreeEvolveCondPhase`. Those families therefore cannot yet be judged on exactly the same waveform-replay footing as GRAPE.
- Open-system performance was evaluated by replay of the best closed-system pulses under Lindblad dynamics, not by noisy re-optimization.
- Wigner validation was performed on the logical basis inputs for the shortlisted final candidates; a broader superposition-probe suite remains future work.

## Suggested Upstreaming
- Add a first-class pulse-export or waveform-bridge route for `SNAP` and `FreeEvolveCondPhase`.
- Expose a convenience transfer-superoperator helper for `HolographicChannel` diagnostics so channel-spectrum analyses do not need local reconstruction.
- Add a reusable replay helper that evaluates GRAPE pulses across enlarged cavity truncations with the same exported schedule.

## Status
COMPLETE
