# Improvement Log: Literature-Informed Selective Pulse Primitives for Dispersive cQED

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- **Missing upstream pulse builders for SNAP and ConditionalPhaseSQR**: this study can validate study-local pulse definitions, but downstream reuse remains awkward until the builders move into `cqed_sim`.
- **No direct open-system optimal control**: noisy replay can rank literature pulse families, but it cannot yet prove global optimality under dissipation.

## Recommended Improvements (P2)
- **Three-mode replay with readout mode**: repeat the best pulse points in a qubit-storage-readout model to quantify measurement-induced degradation. [MEDIUM difficulty]
- **Multilevel transmon replay**: promote the best two-level pulse definitions to `n_tr = 3` and quantify `|f>` leakage. [MEDIUM]
- **Device-calibrated parameter import**: replace the typical-parameter defaults with measured Hamiltonian and noise parameters from a target device. [LOW]

## Nice-to-Haves (P3)
- **Hardware-distorted replay**: add AWG filtering, DAC quantization, and crosstalk for the optimized primitives. [MEDIUM]
- **Cat and binomial encodings**: retest the same pulse families in non-Fock logical subspaces. [HIGH]

## Open Questions
- How much of the residual SQR strict-fidelity gap is removable by echoed selective pulses versus by a subsequent virtual-Z or SNAP cleanup layer?
- Does the literature SNAP sequence remain competitive once qubit dephasing, not only `T1`, is included?
- At what logical dimension does wait-based conditional phase become too distorted by `chi'` and Kerr to remain the best simple phase primitive?

## What Was Tried and Did Not Work
- **Global process-fidelity metric for SQR conditional phase**: treating the selective qubit rotation as one fixed two-qubit unitary underestimated the usable selective-gate quality because the physically relevant target is only defined up to branch-local qubit-$Z$ phases. The production metric now follows the consolidated SQR study and scores branch-resolved relaxed fidelity instead.
- **Direct Kerr-free correction using only inferred `chi_higher`**: replaying the learner with an updated higher-order dispersive shift but still omitting cavity Kerr did not improve the noisy recommendation set. The dominant residual distortion is therefore not captured by `chi_higher` alone.

## Compute & Resource Notes
- Full literature-family production sweep via `scripts/run_study.py`: about 369 seconds on the current workstation, dominated by noisy pulse replay across the SQR and SNAP duration grids.
- Validation via `scripts/validate_results.py`: about 13 seconds, including time-step and truncation reruns for the best operating points.

## Resolved
- **No reusable literature-anchored selective pulse baseline in the repository**: resolved by adding this study with machine-readable optimization results, figures, validation, and a compiled report for downstream SQR and hybrid-control work.
