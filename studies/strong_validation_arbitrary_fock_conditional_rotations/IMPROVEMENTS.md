# Improvement Log: Strong Validation of SQR / CPSQR for Arbitrary Fock-Conditional Qubit Rotations

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- **No true GRAPE benchmark [P1 | HIGH]**: The higher-expressivity comparison uses a basis-expanded constrained waveform family, not GRAPE. This is enough to show that simply adding envelope freedom does not automatically rescue strict control, but it is not the last word on controllability.
- **No stronger-model replay in the saved final table [P1 | MEDIUM]**: The packaged results remain closed-system and two-level in the main reported dataset. A qutrit or Kerr replay of the strict and relaxed highlight cases is the most valuable next validation step.

## Recommended Improvements (P2)
- **Random ensemble expansion [P2 | MEDIUM]**: Extend the current five-seed `chi + chi'` random ensemble to at least 12 seeds and add a `chi_only` companion slice. The present random result is already suggestive, but a larger sample would make the negative strict-SU(2) claim more statistically robust.
- **Explicit left-vs-right Z-gauge study [P2 | MEDIUM]**: The current CPSQR relaxation is a left-multiplicative per-block `R_z(\delta_n)` gauge. Check whether allowing a right-gauge or symmetric gauge changes the apparent reachability boundary for random targets.

## Nice-to-Haves (P3)
- **Expanded symmetry family [P3 | LOW]**: Add frequency-reflected and time-reversed echoed second-half transforms if the current symmetry subset shows near-ties.
- **Kerr replay [P3 | LOW]**: Turn on cavity Kerr for the strongest highlight cases and check whether the operator conclusions move materially.

## Open Questions
- For structured SU(2) targets, why does direct single-pulse SQR remain surprisingly strong up to `N_active = 3`, while the basis-expanded benchmark never overtakes it in the saved dataset?
- For random SU(2) targets, does segmented relaxed control remain near-perfect under the current CPSQR definition as the ensemble size grows, or is the present `5`-seed result unusually favorable?
- For echoed control, is the large relaxed improvement on the stress target a robust refocusing effect or a fragile consequence of the specific half-target decomposition?

## What Was Tried and Did Not Work
- **Symmetry-constrained echo variants on the stress target**: `echo_identical`, `echo_phaseflip`, and `echo_conjugated` all failed badly on strict joint fidelity and remained clearly worse than independently optimized echo. The phase-flip variant was the worst, with mean strict joint fidelity near zero in the saved stress slice.
- **Basis-expanded benchmark as a universal fix**: The higher-expressivity direct family never produced a positive strict-joint gain relative to the single-pulse Gaussian baseline in the saved paired comparisons. This strongly argues that “just add more envelope parameters” is not enough for the hard cases explored here.

## Compute & Resource Notes
- Convention audit will rely on the patched package tests:
  - `tests/test_25_tensor_product_convention.py`
  - `tests/test_20_gaussian_iq_convention.py`
  - `tests/test_23_sqr_additive_amplitude_correction.py`
  - `tests/test_36_targeted_subspace_multitone.py`
- The saved result table contains `114` family evaluations and `114` corresponding machine-readable artifacts.
- The study was completed in staged resumed passes because a naive single run exceeded a one-hour command limit. Incremental JSON/CSV checkpointing prevented data loss.

## Resolved
- **Study scaffold complete**: README, figures, machine-readable results, summary files, and the reproducibility notebook are all present in the study directory.
