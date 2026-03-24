# Improvement Log: Procedural Optimization of Dispersive Readout Pulses with `cqed_sim` `TransmonCavity`

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- **No native readout-specific control objective in `cqed_sim`**: full-model optimization
  must currently be wrapped externally around `simulate_sequence(...)`. This is manageable
  for the study but slows iteration and makes cross-study reuse harder.
- **Two-mode replay does not resolve strong-drive non-QND physics**: across all
  representative families the measured QND defect collapsed to the same `T1` floor
  (`Q_QND approx 0.9933` at the representative `240 ns` protocol), with essentially zero
  family-dependent leakage. This means the current model cannot settle whether nonlinear
  backaction kills the procedural-pulse advantage.
- **Current `T1` curve is only a conservative reference**: the survival-weighted
  matched-filter formula used in the bound hierarchy is useful for trend analysis, but it is
  not a rigorous `T1`-aware upper bound because it does not solve the optimal decay
  change-point discrimination problem.

## Recommended Improvements (P2)
- **Add hardware distortion replay for readout pulses**: include finite bandwidth and simple
  mixer errors once the nominal pulse families are identified.
- **Promote the study-local readout optimization helpers into a reusable package**: only
  after the objective definitions stabilize.
- **Replace the representative detector-efficiency parameter with a calibrated device value**:
  needed for quantitative hardware recommendations.
- **Replay the same pulse families in a richer model**: either a three-mode readout-aware
  model or an explicitly beyond-dispersive model is needed to determine whether the
  procedural gains survive once measurement-induced transitions are present.

## Nice-to-Haves (P3)
- **Extend to the three-mode readout model**: compare the two-mode `TransmonCavity`
  conclusions against `DispersiveReadoutTransmonStorageModel`.
- **Add Bayesian or ensemble robust optimization**: useful if the first-pass pulses prove
  sensitive to `chi`, `kappa`, or amplitude calibration.
- **Benchmark against CLEAR-like analytic depletion pulses**: worthwhile if the procedural
  family shows a real emptying advantage.
- **Export a reusable comparison-table CSV**: the report table already exists, but a
  machine-readable CSV would help later metaanalysis across studies.

## Open Questions
- What fraction of the procedural pulse benefit survives once multilevel leakage is included?
- Is residual cavity occupation actually the dominant predictor of repeated-readout
  inconsistency in this model, or do transmon transitions dominate first?
- Does allowing phase variation buy materially better simultaneous emptying of `|g>` and
  `|e>` trajectories than amplitude-only segmentation?
- Can one derive a rigorous `T1`-aware readout upper bound that remains cheap enough to use
  inside the optimizer?

## What Was Tried and Did Not Work
- **Unrestricted linear piecewise reference as a practical pulse**: at `240 ns`, thev   
  8-segment linear reference reached `F_eta = 0.9847` but left `n_res = 5.57` photons. It
  is useful only as an information-theoretic benchmark, not as a practical readout pulse.
- **Treating the survival-weighted matched-filter curve as a hard upper bound**: the
  realistic replay can exceed that curve because the formula is conservative rather than
  rigorous. Future agents should call it a `T1` reference unless they replace it with a
  sharper derivation.
- **Assuming the analytically nulled family is always exactly feasible**: the optimized
  nulling-tail representative at `240 ns` required clipped tail amplitudes
  (`tail_clipped = True`), so exact linear nulling and bounded-hardware feasibility are
  distinct questions.

## Compute & Resource Notes
- **Full optimization sweep**: approximately `2059 s` wall-clock (`34.3 min`) on this
  workstation for the six-duration, multi-family study run through `run_study.py`,
  excluding the follow-up figure-only rerun after a plotting typo fix.
- **Figure regeneration**: approximately `5 s` from the saved `study_summary.json` and
  `representative_traces.npz`.
- **Validation**: approximately `11 s` for `scripts/validate_results.py`.
- **LaTeX compilation**: `latexmk` was unavailable because the local MiKTeX install lacks
  Perl; direct `pdflatex -> bibtex -> pdflatex -> pdflatex` succeeded.

## Resolved
- **No native arbitrary-envelope readout-chain simulator**: resolved at the package level by
  `ReadoutChain.simulate_waveform(...)` and the new waveform-capable measurement helpers in
  `cqed_sim.measurement`.
- **No stochastic continuous-measurement replay**: resolved at the package level by
  `simulate_continuous_readout(...)` plus monitored/unmonitored collapse partitioning via
  `split_collapse_operators(...)`. The procedural study itself has not yet been re-optimized
  around this new path.
