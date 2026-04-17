# Execution Summary

## Completed Work
- Audited the local `cqed_sim` API reference and package surface relevant to dispersive hybrid control.
- Consolidated the validated evidence from the displacement, waveform-level, literature-informed selective-pulse, strict-SQR, relaxed-CPSQR, arbitrary-conditional-control, and runtime hybrid-unitary studies.
- Added a first-principles timing and phase-budget analysis.
- Generated:
  - `data/synthesis_summary.json`
  - `artifacts/primitive_verdicts.json`
  - `artifacts/analytic_phase_budget.json`
  - `figures/timescale_hierarchy.{png,pdf}`
  - `figures/phase_budget.{png,pdf}`
  - `scripts/reproducibility_notebook.ipynb`
  - `report/report.pdf`

## Main Scientific Result
- The strict ideal primitive gate set does **not** survive literally once realistic dispersive dynamics are enforced.
- A weaker phase-aware constructive library does survive:
  - short branch-compensated displacement,
  - short spectator-limited qubit pulses,
  - relaxed selective control,
  - explicit gauge cleanup such as SNAP or virtual-Z correction.
- A fully pulse-backed non-GRAPE universal stack is **not yet demonstrated**.

## Key Quantitative Results
- Best simple approximate unconditional displacement: mean fidelity `0.9857` at `20 ns` (`|chi|T/2pi = 0.057`).
- Vacuum-calibrated `X_pi` pulse: `0.99984` in vacuum at `40 ns`, but not unconditional across occupied manifolds.
- Best noisy relaxed selective rotation: `0.9903` at `|chi|T/2pi = 1.0`.
- Best inherited strict full-joint ideal SQR case: `0.9988` at `|chi|T/2pi = 5.0`, but only on easier low-dimensional cases.
- Best relaxed arbitrary conditional-control case: relaxed joint-process fidelity essentially `1.0` on a hard random-target slice.
- Best current replay-backed sequence-level candidate: noisy average probe fidelity `0.61448` at total duration `4.432 us`, still requiring GRAPE-derived local surrogates.
