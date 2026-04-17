# Native / Rich Multitone Feasibility for Ideal SQR and CPSQR

## Executive Summary
- Direct native multitone does not achieve a convincing general strict ideal-SQR claim across the tested family set once quartet validation and full joint metrics are enforced.
- Relaxed CPSQR is materially easier than strict SQR: the best CPSQR joint process fidelity reached 1.0000 for echoed_independent on chi_plus_chiprime_staggered_x_na2_chiT5p0.
- The strongest strict-SQR joint result was 0.9988 for reduced_unitary_direct on chi_only_smooth_x_na2_chiT5p0.
- Echo helps in some cases, but the dominant improvement is more robust for relaxed CPSQR than for full strict-SQR success.
- The reduced qubit-only view is consistently more optimistic than the full joint-unitary view.

## Selected Families
- `reduced_unitary_direct`
- `native_direct_strict`
- `gaussian_seed`
- `echoed_independent`
- `echoed_asymmetric`

## Key Findings
- Selected families after the hard-model screen: reduced_unitary_direct, native_direct_strict, gaussian_seed, echoed_independent, echoed_asymmetric.
- Best strict-SQR family: reduced_unitary_direct with reduced quartet 0.9994 and joint process 0.9988.
- Best CPSQR family: echoed_independent with reduced quartet 1.0000 and joint process 1.0000.
- Legacy negative control: saved strict joint metric 0.4949 became reduced quartet 0.8747 and full quartet 0.8747 under the patched-package replay.

## Negative Control
- `artifact_path`: `C:\Users\jl82323\Box\Shyam Shankar Quantum Circuits Group\Users\Users_JianJun\cqed_based_study\studies\ideal_sqr_direct_vs_echoed_multitone\artifacts\cases\chi_plus_chiprime_smooth_x_na3_chiT3p0_direct_multitone.json`
- `saved_restricted_process_fidelity`: `0.4949180359503573`
- `recomputed_strict_joint_process_fidelity`: `0.8241907925828335`
- `recomputed_strict_reduced_single_ground_mean`: `0.8329805318532896`
- `recomputed_strict_reduced_quartet_mean`: `0.8746835797434201`
- `recomputed_strict_full_quartet_mean`: `0.8746835797434201`
