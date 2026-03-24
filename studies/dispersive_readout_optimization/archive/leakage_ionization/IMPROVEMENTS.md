# Improvement Log: Predictive Modeling of Measurement-Induced Leakage and Transmon Ionization at High Readout Power

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- **Phenomenological strong-readout disturbance model**: the production regime map depends
  on occupancy- and slew-activated auxiliary drives, not on a microscopic beyond-dispersive
  transmon model. This is good enough for workflow development and qualitative threshold
  mapping, but not yet trustworthy enough for hard experimental operating limits. Calibrate
  the disturbance coefficients against fridge data or a validated higher-fidelity model.
- **High-power stochastic replay instability**: the native
  `simulate_continuous_readout(...)` path works in tests and moderate-power examples, but
  the strongest representative study points still return NaNs in the continuous replay
  output. Representative histograms therefore still use a deterministic final-state mixture
  fallback at the affected amplitudes.

## Recommended Improvements (P2)
- **Fit the disturbance model to experiment**: use measured post-readout `P_g`, `P_e`,
  `P_f`, and repeated-readout consistency data to fit the onset ratio, occupancy exponent,
  and effective transition scales. This is the main path from qualitative regime maps to
  quantitative predictive value.
- **Stabilize native stochastic replay for strong-drive points**: inspect solver settings,
  monitored-channel scaling, and record post-processing so the representative histograms can
  come from one consistent SME path across the full operating window.
- **Improve the mitigation metric**: the current assignment-fidelity proxy saturates at low
  power, so the mitigation sweep mostly ranks residual photons. A calibrated amplifier/noise
  model would give a more discriminating tradeoff curve.

## Nice-to-Haves (P3)
- **Extend the bandwidth scan to measured transfer functions**: the tested `35-150 MHz`
  low-pass family only weakly perturbed the threshold. Measured line responses may matter
  more.
- **Add automated repeated-readout summary plots**: the study computes QND consistency, but
  future reports would benefit from a standardized conditional-transition visualization.
- **Generalize the strong-readout disturbance helper to multi-mode devices**: useful if the
  storage-assisted or three-mode readout stack becomes the next target.

## Open Questions
- Does a more microscopic model preserve the same regime ordering seen here, namely QND
  breakdown and nearby-level leakage first, with high-manifold occupation staying much
  smaller over the stable simulation window?
- Is the apparent near-threshold histogram advantage in the deterministic fallback model a
  real physical effect or an artifact of collapsing the full dynamics into final-state
  mixtures?
- How much of the eventual experiment mismatch will come from disturbance-model uncertainty
  versus unknown control-line transfer functions?

## What Was Tried and Did Not Work
- **Direct reliance on the original public `cqed_sim` measurement stack**: before the
  framework changes, the package lacked native continuous-readout replay and a reusable
  strong-readout disturbance interface, which blocked a full experiment-facing study.
- **Using solver-step reduction by recompiling the waveform at smaller `dt` as a convergence
  check**: this falsely changed the pulse discretization and made the convergence test look
  worse than it really was. The corrected validation holds the compiled waveform fixed and
  only tightens the ODE solver maximum step.
- **Native stochastic replay at the representative high-power points**: this produced NaNs
  in the saved measurement records even after disabling the disturbance channels for the
  representative replay pass, so the histogram figure still falls back to the deterministic
  mixture proxy at those points.

## Compute & Resource Notes
- **Production full-study runtime**: the rerun after upstreaming the higher-ladder
  continuation took about `893 s` wall-clock time on the current workstation, roughly
  `14.9 min`.
- **Main bottleneck**: the amplitude-duration regime sweep dominates the runtime because each
  point replays two initial states and computes repeated-readout diagnostics.
- **Convergence checks are inexpensive compared with the full sweep**: the larger-truncation
  and smaller-step validation points added little cost relative to the production scan.

## Resolved
- **Missing package support for continuous monitored replay**: resolved by adding
  `simulate_continuous_readout(...)`, monitored/unmonitored collapse partitioning, and
  waveform-capable readout-chain helpers to `cqed_sim`.
- **No reusable strong-readout disturbance helper**: resolved at the package level for the
  primary `ge` and `ef` auxiliary channels via `build_strong_readout_disturbance(...)`.
- **Study-local higher-ladder continuation**: resolved by extending
  `StrongReadoutMixingSpec` with reusable `higher_ladder_scales` support and rerunning the
  production study against the upstreamed package path.
