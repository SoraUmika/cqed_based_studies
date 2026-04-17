# Strong Validation of SQR / CPSQR for Arbitrary Fock-Conditional Qubit Rotations

## Executive Summary
- Strict arbitrary blockwise SU(2) success is much rarer than relaxed per-block Z-gauge success once six-state and cross-block validation are enforced.
- Best strict case: single_pulse_gaussian on inplane_axes with strict joint process fidelity 0.9957.
- Best relaxed case: segmented_relaxed on random_su2 with relaxed joint process fidelity 1.0000.
- The higher-expressivity benchmark is used as the main separator between Gaussian-ansatz failure and deeper control difficulty.
- Cross-block superposition tests remain stricter than single-block reduced diagnostics and are part of the reported success criterion.

## Headline Cases
- Best strict: `single_pulse_gaussian` on `chi_plus_chiprime_inplane_axes_na2_chiT5p0` with strict joint `0.9957`.
- Best relaxed: `segmented_relaxed` on `chi_only_random_su2_na3_chiT5p0_seed9100` with relaxed joint `1.0000`.
- Largest benchmark gain: `random_su2` `chi_plus_chiprime` `N_active=3` `|chi|T/2pi=5.0` with strict-joint gain `-0.0436`.

## Family Summary
- `segmented_relaxed`: strict joint `0.3718`, relaxed joint `1.0000`.
- `single_pulse_gaussian`: strict joint `0.8904`, relaxed joint `0.9307`.
- `echo_independent`: strict joint `0.5612`, relaxed joint `0.9264`.
- `benchmark_basis_expanded`: strict joint `0.7329`, relaxed joint `0.7984`.
- `echo_phaseflip`: strict joint `0.0146`, relaxed joint `0.2406`.
- `echo_identical`: strict joint `0.1509`, relaxed joint `0.2349`.
- `echo_conjugated`: strict joint `0.0999`, relaxed joint `0.1779`.
