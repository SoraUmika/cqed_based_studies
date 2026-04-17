# Improvement Log: Fast and Robust Active Cooling / Vacuum Reset in a Transmon-Storage-Readout cQED System

> Written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- **[P1 | HIGH]** The continuous schemes are still driven with constant square overlaps only: the present comparison shows the right qualitative speed-versus-robustness ordering, but it does not yet test STIRAP-like timing or open-system optimal control.
- **[P1 | HIGH]** The autonomous `L \propto a_s` result remains an auxiliary benchmark rather than a native `cqed_sim` replay.
- **[P1 | MEDIUM]** The continuous-scheme end-of-window residual occupations remain sensitive to the time step: the manuscript now quotes the final-state numbers from the `0.5 ns` baseline and uses the `1.0 ns` sweeps only for robustness ranking, but a tighter convergence campaign is still needed before claiming fully converged absolute residuals.

## Recommended Improvements (P2)
- **[P2 | MEDIUM]** Add pump-aware Stark-shift and parasitic-channel modeling on top of the present effective sideband layer.
- **[P2 | MEDIUM]** Extend the continuous schemes to multi-tone or shaped photon-number-aware driving so higher-Fock cooling is not limited by the `n=1`-centered line choice.
- **[P2 | MEDIUM]** Re-run the Raman-like protocol with engineered larger readout linewidth, since the current device linewidth is not yet fully in the autonomous bad-cavity regime.

## Nice-to-Haves (P3)
- **[P3 | LOW]** Add explicit measurement-backaction and readout-heating models if those become relevant on hardware.

## Open Questions
- How close can a shaped counter-intuitive two-tone protocol get to the autonomous benchmark while staying faster than the present readout lifetime?
- Does the present best pulsed protocol remain best once the sideband controls are embedded in a microscopic pump model?
- What readout linewidth increase is required before the Raman-like protocol becomes clearly superior on both speed and robustness?

## What Was Tried and Did Not Work
- **Constant resonant continuous driving as a generic multi-photon solution**: it remains fast for `n=1` but cools higher-Fock support much less cleanly because the real transmon path is heavily occupied and the fixed carrier is still centered on the lowest manifold.
- **Assuming the most virtual detuned protocol is automatically best on the current device**: on the present readout linewidth it is more robust to transmon coherence, but still slower and less complete than the pulsed ladder within the same wall-clock window.

## Compute & Resource Notes
- Main comparative run: `2510.6 s`
- Continuous candidate scans reused earlier validated pulse winners instead of re-running a global waveform search.

## Resolved
- **Metric-definition mismatch across summary artifacts**: `scheme_summary.csv`, `study_results.json`, and the report now use end-of-run `final_*` values consistently for headline residuals, while tail-averaged `steady_*` values remain available only as diagnostic fields.
- **Pulsed initial-state ladder-depth mismatch**: the pulsed initial-state comparisons now use a matched ladder depth (`n=1` for $|1\rangle$, `n=3` for $|3\rangle$, and the full available ladder for coherent and thermal states) instead of forcing a four-rung sequence for every input state.
- **Single-photon summary versus initial-state timestep mismatch**: the initial-state comparison artifact now uses the same `0.5 ns` baseline step as the headline single-photon summary, so the single-photon rows agree across the report and saved tables.
