# Execution Summary - Iteration 2
Date: 2026-03-28

## Headline Results
- The definitive-study folder now unifies all three requested prior studies plus the patched native-rich extension in one rebuilt report, figure set, and machine-readable output bundle.
- The original strict ideal-SQR baseline remains a negative control: best direct average gate fidelity is about `0.7245`, and best strict process fidelity is about `0.6557`.
- The patched native-rich extension already contains strict ideal-SQR cases above `0.99`; the best loaded strict joint fidelity is `0.99899`.
- The same extension reaches essentially unit conditional-phase SQR on echoed cases; best loaded CPSQR joint fidelity is `0.999999998`.
- The broad-grid landmark is still not met: the mean best strict fidelity across the structured native-rich case set is `0.8392`, well below the requested `0.95` grid-average target.

## New Work Added In This Iteration
- Rebuilt the definitive-study pipeline so it now saves:
  - `data/new_results/objective_summaries.json`
  - `data/new_results/xpi_characterization.json`
  - `data/new_results/spectral_crowding.json`
  - updated `master_results_table.csv`, `fundamental_limits.json`, and `scaling_analysis.json`
- Added new definitive figures for:
  - standalone refocusing-pulse manifold dependence
  - strict-vs-CPSQR joint comparison
  - one-state-vs-quartet validation
  - addressed-manifold scaling
  - spectral crowding
- Rewrote the definitive report so it is a self-contained synthesis rather than a lightweight artifact aggregator.

## Interpretation
- Strict ideal SQR is achievable in the repository corpus, but only conditionally.
- Direct multitone remains the strongest strict-SQR route.
- Echoed constructions remain strongest for the relaxed CPSQR objective, not for strict joint SQR.
- The remaining hard cases look more like a mix of finite refocusing limitations and spectral crowding than like a simple absence of viable strict-SQR solutions.

## Remaining Gaps
- No full echoed-grid rerun with the compromise `X_pi` pulse.
- No fresh hybrid direct-plus-echo study.
- No finite-difference drift sweeps or open-system replay behind the practical ranking.
- No new unconstrained GRAPE upper bound in this definitive folder.
