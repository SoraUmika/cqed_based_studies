# Progress Log

## 2026-04-13
- Initialized the comparative study scaffold in `studies/fast_robust_storage_vacuum_reset_comparison`.
- Classified the task as `DES`, `ANA`, and `OPT`.
- Locked the scope to four scheme classes: pulsed ladder, continuous resonant readout-assisted cooling, detuned Raman-like virtual-transmon cooling, and an auxiliary autonomous-cooling benchmark tied to the effective `L \propto a_s` limit.
- Recorded the first-principles Hamiltonian picture, the expected elimination formulas, and the explicit `cqed_sim` gap note before any standalone reduced-model code is written.
- Identified the key reuse path: earlier validated pulse winners from `storage_active_cooling_gf_sideband` and `gf_sideband_waveform_optimization` will seed the new comparison.

## 2026-04-14
- Science Director review iteration 1 issued `NEEDS_REWORK`.
- The saved artifact package appears usable, but the report was rejected because it lacks literature grounding, a real methods section, quantitative convergence and uncertainty presentation, a working summary figure, and a complete reproducibility appendix.
- Wrote `REVIEW_DIRECTIVE.md` and `FOLLOWUP_PROMPT.md` to drive the next execution iteration.
- Science Director review iteration 2 issued `NEEDS_REWORK`.
- Iteration 2 fixed most of the earlier manuscript-structure problems, but the revised report is still blocked by quantitative inconsistencies across the abstract, summary table, higher-Fock table, convergence table, and machine-readable artifacts.
- Key mismatches recorded in the review directive include the pulsed-ladder single-photon e-fold time reported as both `11.0 ns` and `995 ns`, the contradiction between `scheme_summary.csv` and `study_results.json` for baseline residual occupations, and the `dt <= 1.0 ns` convergence prose claiming at-most-factor-of-two changes despite the printed convergence table showing much larger bright-state and Raman-like changes.
- Wrote the updated `REVIEW_DIRECTIVE.md` and `FOLLOWUP_PROMPT.md` requesting canonical metric reconciliation before the next submission.
- Execution-engineer repair pass: updated `run_study.py` so the headline scheme summary uses canonical end-of-run `final_*` metrics, saved explicit `metric_definitions` into `study_results.json`, and preserved the tail-averaged `steady_*` values only as diagnostic fields.
- Execution-engineer repair pass: fixed the pulsed initial-state comparison policy so the ladder depth matches the initial state (`n=1` for $|1\rangle$, `n=3` for $|3\rangle$, and full available depth for coherent and thermal states) instead of always running a four-rung sequence.
- Re-ran the full study generator to refresh `scheme_summary.csv`, `convergence_summary.csv`, `study_results.json`, and the affected figures with the corrected metric definitions.
- Found and fixed a second artifact-level mismatch: the initial-state comparison had been using the coarser `1.0 ns` sweep step while the headline summary used the finer `0.5 ns` trajectory step, so a helper regeneration script was added and the initial-state CSV/figure were refreshed at the same baseline time step.
- Rewrote the report abstract, results, validation, discussion, and figure caption text so the pulsed-ladder `11 ns` e-fold is explicitly distinguished from its full `1000 ns` pulse-plus-ringdown protocol duration, the higher-Fock table uses the regenerated values, and the convergence narrative now matches the saved convergence table.
- Updated README and improvement-log bookkeeping to mark the study `ACTIVE`, leave convergence unchecked, and record that the notebook has not yet been re-verified after the artifact refresh.
