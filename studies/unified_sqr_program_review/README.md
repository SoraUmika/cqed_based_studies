# Unified SQR Program Review

## Problem Class
AUD | REP | ANA

## Motivation

This study is a **critical audit and unified synthesis** of the nine most recent studies
in this repository. The previous nine studies were run autonomously and contain a mixture
of genuine progress, misleading metrics, unresolved weaknesses, and incomplete validation.

This review:
1. Audits all nine studies for correctness, completeness, and metric consistency.
2. Reruns key scripts to confirm reproducibility.
3. Consolidates all results under a unified notation and metric hierarchy.
4. Produces a single comprehensive report that replaces the fragmented prior reports.

## The Nine Studies Reviewed

| # | Study | Role | Date |
|---|---|---|---|
| 1 | `corrected_sqr_conditioned_rotation_metric` | Establishes correct effective-qubit metric | Mar 27 |
| 2 | `parameterized_waveform_residual_z_cancellation` | Pilot richer-waveform families | Mar 27 |
| 3 | `ideal_sqr_direct_vs_echoed_multitone` | First strict ideal SQR target; baseline | Mar 27 |
| 4 | `native_multitone_sqr_fixed_multiinput_validation` | Multi-input probe ladder validation | Mar 27 |
| 5 | `multitone_sqr_arbitrary_fock_conditional_rotations` | Arbitrary SU(2) targets, broader grid | Mar 28 |
| 6 | `strong_validation_arbitrary_fock_conditional_rotations` | 5-family strict + CPSQR evaluation | Mar 28 |
| 7 | `native_rich_multitone_sqr_cpsqr_feasibility` | Most comprehensive 10-family study | Mar 27 |
| 8 | `the_definitive_ideal_sqr_gate_study` | Aggregation and standalone diagnostics | Mar 28 |
| 9 | `cluster_state_holographic_unified` | Cluster-state ground-sector synthesis | Mar 27 |

## Scientific Questions

1. Can a parameterized multitone qubit-drive waveform implement strict ideal x-axis SQR
   to joint process fidelity > 0.99 in dispersive cQED?

2. Under what conditions (n_active, chi_T, model) does this succeed or fail?

3. Can a structured D+R+SQR or D+R+CPSQR gate sequence implement the restricted
   ground-sector holographic cluster-state transfer to fidelity > 0.99?

## Rerun Results (2026-03-28)

### native_multitone_sqr_fixed_multiinput_validation
Rerun `scripts/run_study.py`. Results match saved data exactly.

| Source study | F_joint | Quartet mean | Quartet min | Assessment |
|---|---|---|---|---|
| Study 2 (parameterized baseline) | 0.0084 | 0.443 | 0.003 | FAIL |
| Study 5 (arbitrary direct) | 0.601 | 0.734 | 0.346 | PARTIAL |
| Study 3 (ideal SQR direct) | 0.824 | 0.875 | 0.720 | PARTIAL |
| Study 1 (corrected unitary-opt) | 0.9986 | 1.000 | 1.000 | PASS |

### the_definitive_ideal_sqr_gate_study
Rerun `scripts/run_full_study.py --full`. Completed in 36.98 s.
All sanity checks pass. Results match saved data exactly.

## Key Conclusions

### Robustly supported
- Strict ideal SQR at n_active=2, chi_T >= 3: F_joint = 0.9990 (reduced_unitary_direct/native_direct_strict)
- CPSQR near-unit fidelity for n_active=2,3: echoed_independent achieves F=1.0000
- Echoed constructions fail for strict SQR due to manifold-dependent X_pi quality
- Multi-input validation is required: single-state metrics are not reliable
- Cluster state: both D+R+SQR and D+R+CPSQR first exceed 0.99 at 3 blocks

### Uncertain
- Whether n_active>=3 failures are ansatz-limited or physics-limited (no GRAPE comparison)
- Robustness to parameter drift
- Hardware relevance of cluster-state results (no open-system validation)

### Not supported
- "Strict SQR is generally achievable": mean F_strict = 0.65-0.72 across all cases
- "Richer envelopes improve strict SQR": complex/basis-expanded families do not help
- "Echoed constructions achieve strict SQR": echoed is consistently worse than direct

## Report

The comprehensive LaTeX report is at `report/report.tex`.

## Status

COMPLETE (audit and synthesis)
