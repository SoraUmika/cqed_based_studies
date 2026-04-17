# Review Request
Date: 2026-04-06
Study: `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed`
Run: `task_runs/multitone_sqr_no_go_decoupled_echo_dispersive_cqed`
Status: READY_FOR_REVIEW

## What This Study Claims
1. Under the strict simultaneous shared-line dispersive model with no artificial per-tone detuning and amplitude-plus-azimuth corrections only, exact ideal multitone SQR is generically unavailable.
2. The obstruction is visible analytically as a blockwise second-order `Z` term and numerically in the full exact shared-line propagation.
3. A stronger decoupled-block approximation does allow ideal SQR exactly, but it is not the same physical shared-line problem.
4. The echoed sequence `half-SQR -> pi -> half-SQR -> pi` suppresses some residual `Z` accumulation only approximately and does not rescue the strict gate.

## Evidence Prepared For Review
- Full report with derivation, numerics, limitations, and final verdict:
  - `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/report/report.pdf`
- Machine-readable results:
  - `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/data/study_results.json`
  - `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/data/study_summary.json`
  - `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/data/analytic_summary.json`
  - `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/data/validation_summary.json`
  - `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/data/prior_audit.json`
- Reproducibility notebook:
  - `studies/multitone_sqr_no_go_decoupled_echo_dispersive_cqed/scripts/reproducibility_notebook.ipynb`

## Requested Review Focus
Please review the study as a science-director pass with emphasis on:
1. Whether the two-block and many-block no-go statements are scoped and worded carefully enough.
2. Whether the numerical evidence is sufficient to support the strict shared-line conclusion without overclaiming.
3. Whether the report keeps the decoupled-block success clearly separated from the physical simultaneous shared-line problem.
4. Whether the echo section states precisely what is canceled, what is not, and under which assumptions.
