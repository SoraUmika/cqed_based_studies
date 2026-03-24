# Improvement Log: Dispersive Readout Pulse Optimization

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- **Uncalibrated disturbance model**: The phenomenological strong-drive mixing coefficients (ge scale 0.09, ef scale 0.05) are not calibrated to experimental data. The leakage threshold map is therefore useful for design-space exploration but not yet a predictive hardware policy tool. Calibrate against experiment or a microscopic beyond-dispersive model.
- **Stochastic replay instability at high power**: The `simulate_continuous_readout` wrapper becomes numerically unstable at the strongest operating points. Representative histograms at those points use a deterministic final-state mixture fallback. Stabilize the native continuous-readout replay so all operating points use one consistent stochastic path.

## Recommended Improvements (P2)
- **Calibrated amplifier/noise model**: Re-run the mitigation scan with an experimentally calibrated amplifier model so the assignment-fidelity proxy no longer saturates across the entire search window. [MEDIUM difficulty]
- **Higher-ladder continuation presets**: Add calibrated presets or fitting utilities for the reusable higher-ladder continuation scales so high-manifold spreading can be parameterized without per-study retuning. [MEDIUM]
- **Unified cost function**: The four component studies used different objective weightings. A single unified cost function incorporating SNR², residual photons, QND preservation, and leakage would enable direct cross-model comparison. [MEDIUM]
- **Active depletion revisited**: The single-segment heuristic failed in the linear study. Revisit with a properly designed multi-segment depletion protocol informed by the nulling-tail algebraic insight. [LOW]

## Nice-to-Haves (P3)
- **Extended bandwidth scan**: Current bandwidth scan covers only 35–150 MHz low-pass range. Extend to measured transfer functions or more extreme distortion profiles. [LOW]
- **Driven depletion + readout integration**: Combine active cavity emptying with the readout information-extraction phase in a single optimized pulse. [MEDIUM]
- **Temperature sweep**: Map how the readout prescription changes across a realistic fridge temperature range (10–50 mK). [MEDIUM]

## Open Questions
- Does a more microscopic model preserve the same ordering of regimes (QND breakdown first, broad high-manifold spreading later)?
- Is the near-threshold histogram improvement seen in the fallback record model physically real, or is it an artifact of the final-state mixture approximation?
- Should optimal readout always be at the midpoint drive, or can detuned readout improve QND character at high power?
- At what model complexity does GRAPE start outperforming structured families?

## What Was Tried and Did Not Work
- **Single-segment active depletion** (Study 1): Heuristic depletion pulse did not beat passive ring-down in the linear model. The implemented approach used a fixed-phase depletion tone at the estimated cavity frequency; the approach fails because both conditioned amplitudes must be simultaneously nulled.
- **Unbounded GRAPE in the linear model** (Study 1): GRAPE converges to the same square-like solution within numerical precision. Under peak-amplitude constraints, GRAPE adds no meaningful SNR² advantage in the linear dispersive setting.
- **Free procedural segments under hardware/QND model** (Study 3): While excellent for information extraction, free procedural pulses pay higher QND cost than ring-hold due to sharp phase structures and high-occupancy excursions. Not recommended as default under the full model.
- **Bounded Powell refinement** (Study 3): Unstable for procedural segments under the rich model; setting n_local=0 reduced optimizer fragility without changing the frontier conclusion.

## Compute & Resource Notes
- Linear model sweeps: ~seconds per point, full sweep <1 minute.
- Multilevel cqed_sim replay: ~2–5 seconds per pulse evaluation at n_tr=3, n_cav=14.
- Full hardware+nonlinear replay: ~5–10 seconds per point.
- Extended transmon (5-level, 18 Fock) leakage study: ~10–30 seconds per point; full regime map ~hours.
- Convergence checks (increased truncation): factor ~3× slower per dimension increase.

## Resolved
(None yet)
