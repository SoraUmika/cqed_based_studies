# Science Directive

Date: 2026-04-13
Study: `studies/fast_robust_storage_vacuum_reset_comparison`
Run: `task_runs/fast_robust_storage_vacuum_reset_comparison`

## Problem Class
DES, ANA, OPT

## Central Question
Which active-cooling architecture is best for a storage-transmon-readout cQED system once speed, transmon-decoherence robustness, leakage, detuning/amplitude sensitivity, imperfect reset, and thermal loading are all evaluated together?

## Required Scheme Classes
1. Pulsed transmon-assisted ladder cooling.
2. Continuous resonant readout-assisted sideband cooling.
3. Continuous detuned Raman-like / virtual-transmon cooling.
4. A documented autonomous-cooling benchmark tied to the effective `L \propto a_s` limit.

## Non-Negotiable Modeling Boundaries
- Use `cqed_sim` as the primary simulation engine.
- Reuse earlier validated pulse winners where possible instead of re-optimizing everything.
- Keep the reduced eliminated benchmark clearly labeled as auxiliary because `cqed_sim` does not natively expose it.
- Make the evidence-to-claim chain explicit: every major recommendation must point to saved figures, tables, or artifacts.

## Deliverables
- Executive summary with fastest, most robust, and best overall scheme.
- Scheme-by-scheme comparison tables and robustness heatmaps.
- Time-domain cooling traces for pure, coherent, and mixed initial states.
- Saved machine-readable artifacts, report, PDF, and reproducibility notebook.
