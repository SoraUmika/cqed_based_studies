# Predictive Modeling of Measurement-Induced Leakage and Transmon Ionization at High Readout Power

## Problem Class
ANA | DES | OPT

## Motivation
High-power dispersive readout shortens measurement windows and can improve detector-limited
separation, but the same operating regime can drive the transmon out of the computational
subspace and break the QND assumption. This study extends `cqed_sim` so the framework can
replay hardware-distorted readout envelopes, partition monitored versus unmonitored noise,
and model effective strong-readout disturbance channels inside the existing package
workflow rather than through an ad hoc simulator.

The main experimental question is where readout crosses from safe dispersive operation into
measurable leakage and then, potentially, into broader high-manifold occupation that is
operationally ionization-like. The study is intended as an experiment-facing regime map and
validation plan, not as a literal continuum-ionization calculation.

## Goals
1. Build a reusable `cqed_sim` workflow for strong-readout replay with level-resolved
   populations, leakage metrics, residual-photon diagnostics, and representative readout
   records.
2. Map leakage and repeated-measurement breakdown versus readout amplitude, pulse duration,
   drive detuning, and control-chain bandwidth.
3. Distinguish moderate nearby-level leakage from broader high-manifold population using an
   operational ionization metric.
4. Validate the reported thresholds with truncation and solver-step convergence checks.
5. Produce a report-ready set of figures, quantitative thresholds, and experiment-facing
   validation recommendations.

## Methods
### `cqed_sim` modules/functions used
- `cqed_sim.core.DispersiveTransmonCavityModel`
- `cqed_sim.core.FrameSpec`
- `cqed_sim.pulses.Pulse`
- `cqed_sim.pulses.hardware.HardwareConfig`
- `cqed_sim.sequence.SequenceCompiler`
- `cqed_sim.sim.SimulationConfig`
- `cqed_sim.sim.prepare_simulation`
- `cqed_sim.sim.simulate_sequence`
- `cqed_sim.sim.split_collapse_operators`
- `cqed_sim.sim.extractors.transmon_level_populations`
- `cqed_sim.sim.extractors.reduced_subsystem_state`
- `cqed_sim.measurement.ReadoutResonator`
- `cqed_sim.measurement.ReadoutChain`
- `cqed_sim.measurement.build_strong_readout_disturbance`
- `cqed_sim.measurement.simulate_continuous_readout`

### Framework improvements applied in this study
1. Added reusable readout-chain waveform replay helpers so the measurement stack can accept
   arbitrary compiled envelopes instead of only idealized scalar drives.
2. Added `split_collapse_operators(...)` so monitored continuous-readout channels can be
   separated cleanly from unmonitored Lindblad terms.
3. Added `simulate_continuous_readout(...)` to expose a native stochastic replay path
   through the package measurement API.
4. Added `build_strong_readout_disturbance(...)` to generate occupancy- and slew-activated
   effective disturbance envelopes for strong-readout studies, including reusable
   higher-ladder continuation via `higher_ladder_scales`.

### Production modeling hierarchy
1. Deterministic multilevel dispersive replay with hardware-distorted readout envelopes.
2. Effective strong-readout disturbance channels added as auxiliary transition drives.
3. Representative continuous-readout replay at selected operating points, with a documented
   deterministic mixture fallback when the native stochastic path becomes unstable.

## Expected Outcomes
- Regime maps separating approximately QND operation from moderate leakage and stronger
  non-QND readout.
- Quantitative threshold estimates for the chosen device and pulse family.
- Detector-facing plots that can guide validation experiments and safer high-speed readout
  choices.
- A clear list of what is now reusable in `cqed_sim` and what still needs upstream work.

## Assumptions
- Device parameters were fixed at:
  `omega_q / 2pi = 6.150 GHz`,
  `omega_c / 2pi = 8.597 GHz`,
  `alpha / 2pi = -255 MHz`,
  `chi / 2pi = -2.84 MHz`,
  `K / 2pi = -28 kHz`,
  and `kappa / 2pi = 2.4 MHz`.
- The production scan used `n_tr = 5`, `n_cav = 18`, and `dt = 4 ns`, with convergence
  cross-checks at `n_tr = 7`, `n_cav = 24`, and `max_step = dt / 2` while holding the
  compiled waveform fixed.
- The operational high-manifold metric used `P_ion = sum_{m >= 3} P_m`, so `n_ion = 3`.
- The visible-leakage benchmark used throughout the regime map was `P_leak > 1e-4`.
- The baseline pulse family was a cosine envelope with a `75 MHz` low-pass control-chain
  replay model.
- Strong-readout disturbance strengths were calibrated operationally inside the study using
  occupancy- and slew-activated auxiliary envelopes; they are not yet fitted to fridge
  data.

## Convergence Criteria
- Key observables (`P_leak`, `P_ion`, peak cavity occupancy, and repeated-readout defect)
  must remain stable under larger Hilbert-space truncation and smaller solver step.
- Threshold locations are considered stable if the leakage benchmark crossing does not move
  under the production convergence checks.

## Expected Outcomes vs. Current Results
- The implemented workflow now produces regime maps, detuning and bandwidth scans,
  mitigation comparisons, representative histogram artifacts, and a compiled report.
- In the production run, the representative `240 ns` slice crosses the visible-leakage
  benchmark at `5.5 MHz`.
- The threshold is delayed to `7.5 MHz` at `120 ns` and shifts down to `4.5 MHz` by
  `480 ns`.
- The strongest production point tested (`480 ns`, `7.5 MHz`) reached
  `mean_p_leak = 2.18e-3`,
  `mean_p_ion = 1.25e-6`,
  and `qnd_defect = 0.495`.

## Known Limitations
- The strong-readout disturbance layer is phenomenological. It captures threshold-like
  leakage and QND breakdown trends, but it is not a microscopic beyond-dispersive
  ionization Hamiltonian.
- Native continuous-readout replay now exists in `cqed_sim`, but it becomes numerically
  unstable at the stronger operating points used for the representative histogram study.
- Those histogram figures therefore fall back to a deterministic final-state mixture model
  at the affected amplitudes.
- The operational ionization metric now uses the reusable higher-ladder continuation in
  `cqed_sim.measurement.strong_readout`, but the chosen scale factors remain
  phenomenological and have not yet been fitted to experiment.
- No fridge data were available in this run, so validation is internal (sanity plus
  convergence) rather than simulation-versus-experiment agreement.
- The mitigation scan is detector-limited in the current proxy model: assignment fidelity
  saturates early, so residual cavity population is the main discriminator.

## Suggested Upstreaming
- Harden the continuous-readout replay path so high-power operating points can be simulated
  without falling back to deterministic mixture histograms.
- Add calibrated presets or fitting utilities for the reusable higher-ladder continuation
  scales inside the strong-readout helper API.
- Add experiment-calibrated disturbance presets or fitting utilities for strong-readout
  studies.

## Status
COMPLETE
