# Execution Summary

Date: 2026-03-31
Study: `studies/fock_resolved_black_box_sqr_inference`
Run: `task_runs/fock_resolved_black_box_sqr_inference`

## Headline Findings
- The allowed measurement set `{qubit tomography, calibrated displacement, known dispersive wait}` does **not** identify `p_n` and sector-resolved `Z_n` separately.
- Wait-only and combined protocols are both full-rank on the recoverable weighted-transverse subspace `u_n = p_n (X_n + i Y_n)`.
- Displacement-only remains rank deficient and empirically non-informative for sector separation.
- The combined protocol is a strong coherence witness: in the exact coherent benchmark its diagonal-model residual rises to `0.1387`, while the wait-only residual remains at machine precision.

## Quantitative Results
- Single-qubit tomography baseline:
  - Best mean MLE fidelity at the high-shot end: `0.9986`
- Recoverable-subspace protocol benchmark:
  - Wait-only mean weighted RMSE at `10^4` shots: `4.14e-3`
  - Combined mean weighted RMSE at `10^4` shots: `2.95e-3`
  - Displacement-only mean weighted RMSE at `10^4` shots: `2.17e-1`
- Pulse-level black-box validation:
  - Best case: combined protocol on `pulse_optimized_mix_g` with weighted RMSE `2.60e-3`
  - Worst case: combined protocol on `pulse_optimized_leakage_replay` with weighted RMSE `2.57e-1`

## Validation Status
- Sanity checks: passed
- Convergence checks: passed
- Literature comparison: not applicable for this original protocol-design study
- Automated validation script: passed all checks

## Artifacts Produced
- Report: `studies/fock_resolved_black_box_sqr_inference/report/report.pdf`
- Reproducibility notebook: `studies/fock_resolved_black_box_sqr_inference/scripts/reproducibility_notebook.ipynb`
- Machine-readable results: `studies/fock_resolved_black_box_sqr_inference/data/study_results.json`
- Validation summary: `studies/fock_resolved_black_box_sqr_inference/data/validation_summary.json`
