# Execution Summary

Date: 2026-04-14
Study: `studies/fast_robust_storage_vacuum_reset_comparison`
Run: `task_runs/fast_robust_storage_vacuum_reset_comparison`

## Main Findings
- Best overall physical scheme remains `Continuous Raman-like`.
- Fastest practical physical scheme remains `Pulsed ladder`.
- The report/artifact contradictions from review iteration 2 have been removed by enforcing a single canonical headline metric set based on end-of-run `final_*` values.

## Key Numbers
- Common continuous comparison window: `3999.9 ns`
- Selected pulsed single-photon protocol duration: `1000.0 ns`
- Pulsed single-photon headline metrics: `tau_e = 11.0 ns`, `final_storage_n = 0.00406`
- Bright-state headline metrics: `tau_e = 13.0 ns`, `final_storage_n = 0.00109`
- Raman-like headline metrics: `tau_e = 243.5 ns`, `final_storage_n = 0.000608`
- Main full-run runtime: `2510.6 s`

## Repair Actions Completed
- Reconciled `scheme_summary.csv` and `study_results.json` onto canonical `final_*` headline values while preserving `steady_*` only as diagnostic fields.
- Regenerated the initial-state comparison with matched pulsed ladder depth and the same `0.5 ns` baseline step used by the headline single-photon summary.
- Rewrote the report sections that previously mixed `11 ns`, `995 ns`, and `1990 ns` pulsed timings without defining whether they referred to an e-fold or a full protocol duration.
- Replaced the incorrect `dt <= 1.0 ns` convergence prose with a caveated statement that matches the saved convergence table.

## Remaining Open Items
- The reproducibility notebook has not been re-run after the artifact refresh.
- Continuous-scheme end-of-window residuals remain time-step sensitive, so convergence is still tracked as incomplete in the README.
