# Procedural Optimization of Dispersive Readout Pulses with `cqed_sim` `TransmonCavity`

## Problem Class
ANA | OPT | DES

## Motivation
Dispersive cQED readout is usually tuned with a small collection of heuristic pulse
shapes even though the control objective is intrinsically multi-objective: strong
state discrimination, low residual cavity occupation, short latency, and minimal
QND damage. The existing local study `studies/readout_pulse_optimization` established
useful linear-model bounds, but it stopped short of answering the more important
question for practical use inside `cqed_sim`: whether a low-dimensional,
physics-informed, procedurally parameterized pulse family remains competitive once
we replay it through the multilevel `DispersiveTransmonCavityModel` with cavity
damping and qubit decoherence.

This study treats that question as an open problem. The goal is not to assume that
procedural optimization works, but to determine where it helps, where linear theory
remains predictive, and where nonlinear/QND constraints erase the apparent advantage.

## Goals
1. Build an explicit hierarchy of readout bounds versus total duration `T`:
   ideal linear-dispersive, detection-limited, `T1`-limited, and realistic
   multilevel replay bounds.
2. Define reusable readout objectives and metrics covering distinguishability,
   matched-filtered I/Q separation, residual photons, peak photon safety,
   QND state preservation, repeated-readout consistency, and robustness.
3. Implement study-local reusable readout optimization utilities on top of
   `cqed_sim` for multiple pulse families and constrained optimization.
4. Compare baseline heuristic pulses, procedural segmented pulses, smooth
   basis-expanded pulses, and a higher-dimensional reference parameterization.
5. Quantify the time-to-fidelity, fidelity-to-QND, and fidelity-to-residual-photon
   frontiers for the nominal `TransmonCavity` readout configuration.
6. Test whether active late-time cavity emptying can keep `n_g(T)` and `n_e(T)`
   small without giving up most of the intermediate-time signal separation.
7. Deliver a scientifically honest recommendation on whether procedural readout
   pulse optimization is a promising practical direction in this stack.

## Methods
### `cqed_sim` modules/functions to use
- `cqed_sim.core.DispersiveTransmonCavityModel`
- `cqed_sim.core.FrameSpec`
- `cqed_sim.pulses.Pulse`
- `cqed_sim.sequence.SequenceCompiler`
- `cqed_sim.sim.SimulationConfig`, `prepare_simulation`, `simulate_sequence`
- `cqed_sim.sim.extractors.reduced_qubit_state`
- `cqed_sim.sim.extractors.qubit_conditioned_mode_moments`
- `cqed_sim.sim.extractors.transmon_level_populations`
- `cqed_sim.sim.noise.NoiseSpec`
- `cqed_sim.measurement.ReadoutResonator`
- `cqed_sim.measurement.ReadoutChain`

### Study-local extension plan
- `scripts/readout_opt/pulse_families.py`: baseline, segmented, basis-expanded,
  and higher-dimensional pulse parameterizations.
- `scripts/readout_opt/bounds.py`: linear-dispersive trajectory solver, matched-filter
  SNR bounds, detector-efficiency and `T1` upper bounds.
- `scripts/readout_opt/metrics.py`: assignment fidelity, residual/peak photon metrics,
  leakage, QND preservation, repeated-readout consistency, and robustness summaries.
- `scripts/readout_opt/simulate.py`: `cqed_sim` replay helpers for cavity-drive pulses
  on `DispersiveTransmonCavityModel`.
- `scripts/readout_opt/optimize.py`: constrained SciPy-based optimization wrappers.
- `scripts/readout_opt/experiments.py`: sweeps over duration, robustness, and Pareto
  fronts.
- `scripts/readout_opt/plots.py`: publication-quality figures.

### Required framework-gap documentation
1. `cqed_sim.measurement.ReadoutResonator.response_trace()` currently supports
   constant-amplitude drive response, not arbitrary time-varying readout envelopes.
   We will therefore keep the earlier exact linear ODE/integral machinery for the
   ideal-bound layer and use `cqed_sim` itself for all realistic multilevel replay.
2. `cqed_sim` does not currently provide a built-in readout-specific optimizer with
   objectives such as matched-filter SNR, residual-photon constraints, or QND
   preservation. The optimization layer in this study is therefore study-local glue
   code that evaluates objectives through supported `cqed_sim` simulations.
3. The current simulator provides deterministic Lindblad evolution but not a native
   stochastic continuous-measurement unraveling. Repeated-readout QND metrics will be
   estimated from two-pulse replay, state-preservation probabilities, and readout
   classification consistency rather than a full quantum-trajectory measurement record.

## Assumptions
- Nominal two-mode readout model: `DispersiveTransmonCavityModel` interpreted as
  transmon + readout resonator.
- Initial nominal parameters:
  `omega_q/2pi = 6.150 GHz`, `omega_c/2pi = 8.597 GHz`,
  `alpha/2pi = -255 MHz`, `chi/2pi = -2.84 MHz`,
  `kappa/2pi = 2.4 MHz`, `kerr/2pi = -28 kHz`.
- Initial truncations: `n_tr = 3`, `n_cav = 14`; representative points must be
  checked at larger truncation.
- Nominal decoherence for realistic replay: `T1 = 30 us`, `T2 = 20 us`; cavity loss
  set by `kappa`.
- Ideal bound uses `eta = 1`; detector-limited comparisons additionally report
  `eta = 0.35` as a representative chain efficiency.
- Readout durations will be scanned over approximately `T in [80, 900] ns` unless
  convergence or compute cost forces a narrower window.
- Default cavity-drive replay uses a rotating frame with `omega_q_frame = omega_q`
  and `omega_c_frame = omega_d`, so the readout pulse is represented as a baseband
  envelope on the cavity channel.
- Convergence targets:
  time-step stability better than `5e-3` on the primary metrics and truncation
  stability better than `1e-3` on representative points.

## Expected Outcomes
- A sharper statement of when the linear midpoint-drive bound is still predictive and
  when multilevel/QND effects lower the practical ceiling.
- A quantified answer to whether a 3- to 5-segment procedural pulse can reproduce
  most of the benefit of a more flexible parameterization.
- A measured time overhead for active cavity emptying relative to unconstrained
  high-SNR pulses.
- A recommendation of one pulse family for practical use and one higher-dimensional
  family for benchmarking.

## Known Limitations
- This study inherits the lack of native stochastic measurement trajectories in
  `cqed_sim`, so QND metrics will be replay-based rather than trajectory-based.
- The nominal detector-efficiency value is representative rather than device-calibrated.
- Real hardware distortions (finite AWG bandwidth, mixer skew, amplifier saturation)
  are not yet included unless they become necessary during validation.
- The nominal self-Kerr value is a weak-nonlinearity placeholder that must be tested
  for sensitivity rather than treated as a fully calibrated device constant.
- In the present two-mode dispersive replay, family-dependent QND degradation is almost
  entirely the ordinary `T1` floor; stronger measurement-induced state-mixing channels
  are not resolved without a richer model.
- The plotted `T1` curve is a conservative survival-weighted reference rather than a
  rigorous information-theoretic upper bound.

## Suggested Upstreaming
- Add arbitrary-envelope readout simulation to `cqed_sim.measurement`.
- Add a readout-objective optimization surface that can score matched-filter SNR,
  empty-cavity constraints, and QND penalties.
- Add a helper for two-pulse replay-based repeated-readout diagnostics.

## Status
COMPLETE
