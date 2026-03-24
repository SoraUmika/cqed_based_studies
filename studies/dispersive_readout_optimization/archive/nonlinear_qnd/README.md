# Nonlinear-QND Limits and Hardware-Realistic Optimization of Dispersive Readout Pulses in cQED

## Problem Class
OPT | ANA | DES

## Motivation
The completed study [studies/procedural_readout_pulse_sequence_optimization/README.md](studies/procedural_readout_pulse_sequence_optimization/README.md) established that low-dimensional procedural readout pulses can outperform square-pulse baselines in the current two-mode dispersive replay stack. It also identified the critical unresolved gap: repeated-readout degradation in that stack is dominated almost entirely by ordinary $T_1$, so the simulator could not yet determine whether strong-drive, family-dependent non-QND effects erase the procedural-pulse advantage.

This follow-up study addresses that open question directly. The goal is not to re-run the earlier optimizer with a larger budget, but to upgrade the replay stack so it includes both hardware transport constraints and an explicit strong-drive non-QND mechanism that can generate family-dependent readout disturbance.

## Goals
1. Reuse the prior pulse-family and optimization infrastructure where it remains valid, and document the exact parts that need to be upgraded.
2. Introduce a richer readout replay hierarchy that includes multilevel transmon relaxation, readout-aware hardware distortion, and an explicit effective strong-drive mixing mechanism that scales with instantaneous readout occupancy near the dispersive critical-photon scale.
3. Evaluate square, smoothed-square, ring-hold, free procedural segments, analytic nulling-tail, Fourier-basis, and a higher-dimensional piecewise reference family under matched total durations.
4. Track detector-limited discrimination, realistic replay fidelity, residual and peak photon occupancy, repeated-readout consistency, leakage, measurement-induced transition probability, and balanced practical score separately.
5. Produce a sharper benchmark hierarchy for realistic readout performance: linear detector limit, prior $T_1$ reference, high-dimensional physical reference under the richer model, and a QND-constrained practical frontier.
6. Determine whether low-dimensional procedural structure still captures most of the useful gain once nonlinear-QND penalties and hardware transport constraints are included.

## Methods
### `cqed_sim` modules/functions to use
- `cqed_sim.core.DispersiveTransmonCavityModel`
- `cqed_sim.core.DispersiveReadoutTransmonStorageModel`
- `cqed_sim.core.FrameSpec`
- `cqed_sim.pulses.Pulse`
- `cqed_sim.pulses.hardware.HardwareConfig`
- `cqed_sim.sequence.SequenceCompiler`
- `cqed_sim.sim.SimulationConfig`, `prepare_simulation`
- `cqed_sim.sim.noise.NoiseSpec`
- `cqed_sim.sim.extractors.transmon_level_populations`
- `cqed_sim.sim.extractors.qubit_conditioned_mode_moments`
- `cqed_sim.sim.extractors.readout_response_by_qubit_state`

### Study-local extension plan
- `scripts/readout_opt/config.py`: extend the prior nominal configuration with explicit hardware profiles, multilevel lifetimes, and effective-mixing parameters.
- `scripts/readout_opt/simulate.py`: add richer replay regimes combining multilevel Lindblad evolution, `SequenceCompiler` hardware distortion replay, and an effective strong-drive mixing pulse mapped onto explicit qubit transition channels.
- `scripts/readout_opt/metrics.py`: add measurement-induced transition metrics, QND defect decomposition, and hardware robustness summaries.
- `scripts/readout_opt/experiments.py`: run richer-model frontiers, distorted replay, re-optimization after distortion, QND stress sweeps, and benchmark hierarchies.
- `scripts/readout_opt/plots.py`: update figure generation for the new frontiers, distortion visualizations, and comparison tables.

### Required framework-gap documentation
1. `cqed_sim` provides native waveform hardware replay through `SequenceCompiler(hardware=...)`, so bandwidth, quantization, timing quantization, and IQ skew will be handled natively rather than in study-only post-processing.
2. `cqed_sim` does not currently expose a built-in, microscopic strong-drive measurement-induced transition model for arbitrary readout envelopes. This study therefore adds a documented effective mixing layer inside the normal `cqed_sim` pulse-and-simulation workflow by constructing occupancy-dependent ancilla drive terms on top of the supported multilevel Hamiltonian.
3. `cqed_sim` still uses time-independent Lindblad collapse operators in the public replay path. Any explicitly measurement-induced, occupancy-dependent disturbance channel must therefore be represented through effective driven mixing and multilevel ladder loss rather than a native time-dependent collapse operator.
4. The resulting QND-aware benchmark is model-based rather than theorem-tight: it is a reference envelope produced by the richer physical model and a high-dimensional comparison family, not a closed-form proof of optimality.

## Assumptions
- Base frequencies and nonlinearities start from the existing device values: $\omega_q / 2\pi = 6.150$ GHz, $\omega_r / 2\pi = 8.597$ GHz, $\alpha / 2\pi = -255$ MHz, $\chi / 2\pi = -2.84$ MHz, $\kappa_r / 2\pi = 2.4$ MHz, $K_r / 2\pi = -28$ kHz.
- The main optimization layer uses the two-mode `DispersiveTransmonCavityModel` as the readout mode, augmented by multilevel relaxation and effective strong-drive mixing. A three-mode `DispersiveReadoutTransmonStorageModel` validation slice may be used for cross-checks where compute allows.
- Multilevel transmon relaxation will be represented with explicit ladder times `(T1_ge, T1_fe)`, starting from approximately `(30 us, 12 us)` and varied during stress tests.
- Effective strong-drive disturbance will be parameterized by the instantaneous estimated readout occupancy relative to $n_\mathrm{crit}$ and by waveform slew, so that fast, high-occupancy protocols can induce family-dependent excitation, relaxation, and leakage.
- Hardware replay will include finite bandwidth, zero-order hold, amplitude quantization, timing quantization, overall gain error, IQ gain mismatch, quadrature skew, and image leakage through `HardwareConfig`.
- Nominal readout durations will be scanned over approximately $T \in [96, 720]$ ns using the prior duration grid unless convergence or runtime constraints require a reduced subset.
- Validation targets remain: representative truncation stability better than $10^{-3}$ on key metrics where feasible, and practical metric stability better than $3 \times 10^{-2}$ under finer time resolution and larger truncation.

## Expected Outcomes
- A direct answer to whether the earlier procedural-pulse advantage survives after the replay stack includes both hardware transport distortion and a family-dependent non-QND mechanism.
- A quantitative map of when analytic nulling tails remain useful once finite bandwidth and slew constraints distort the programmed waveform.
- A clearer tradeoff surface separating detector-limited discrimination, physically achievable fidelity, repeated-readout consistency, and genuine QND preservation.
- A practical recommendation for which pulse family should be used in short-, intermediate-, and long-duration readout regimes under the richer model.

## Results Summary
- The richer replay hierarchy changes the control conclusion, but it does not erase the value of procedural structure. Across the executed three-duration frontier `(96, 240, 496) ns`, the best rich-model family is `ring_hold`, with matched-filter fidelities `0.6480`, `0.9599`, and `0.9982`, respectively.
- The free `procedural_segments` family still reaches the highest representative-duration information fidelity among the explicitly compared rich-model families at `240 ns` (`F_eta = 0.9547`), but it pays a visibly larger QND cost than smoother structured families: `Q_QND = 0.9886` and measurement-induced transition probability `4.66e-3`.
- Hardware distortion and nonlinear backaction affect different parts of the performance budget. For the representative `procedural_segments` pulse at `240 ns`, hardware replay accounts for most of the nominal fidelity loss (`0.9628 -> 0.9552`), while the nonlinear mixing layer accounts for most of the new QND penalty (`Q_QND: 0.9933 -> 0.9898`).
- The higher-dimensional `piecewise_reference` family remains useful as a benchmark, but not as a practical recommendation. Its information-optimized point reaches `F_eta = 0.9760` at `240 ns`, while leaving `4.31` residual photons and reducing QND preservation to `0.9802`.
- Under amplitude stress (`1.3x` the nominal amplitude probes), `procedural_segments` degrades the fastest among the tested rich-model families (`Q_QND = 0.9618`, induced transition `3.15e-2`), while smoother structured families such as `fourier_basis` fail more gracefully (`Q_QND = 0.9890`, induced transition `4.33e-3`).

## Validation
- Sanity checks passed: the analytic nulling-tail construction drives both conditioned cavity amplitudes to machine zero in the linear model (`|alpha_g(T)| = |alpha_e(T)| = 3.54e-16`), and the zero-drive replay produces zero residual photons with only `1.38e-12` induced transition probability.
- Convergence checks passed at the representative operating point: the fine-time-step replay changes fidelity by `1.43e-3`, and the larger cavity truncation changes fidelity by `4.73e-5`.
- Hardware replay is active rather than nominally bypassed: the distortion-activation diagnostic is `6.33e-1`.
- The richer model generates family-dependent disturbance as intended: the induced-transition spread across families is `2.72e-2` in the validation sweep.
- The validation script `scripts/validate_results.py` passes against the saved study artifacts.

## Known Limitations
- The effective strong-drive mixing layer is intentionally phenomenological; it is designed to expose control consequences of nonlinear-QND breakdown, not to claim a complete microscopic derivation of beyond-RWA readout physics.
- A full stochastic continuous-measurement unraveling is still outside the current public `cqed_sim` workflow, so repeated-readout analysis remains replay-based rather than trajectory-conditioned.
- The three-mode validation layer may remain limited to selected representative pulses if full-duration optimization proves too expensive.
- Any claimed benchmark beyond the linear detector limit should be interpreted as a physically informed reference within the chosen model hierarchy, not as a theorem-tight upper bound.

## Suggested Upstreaming
- Add a reusable occupancy-dependent readout-disturbance interface to `cqed_sim` so strong-drive measurement backaction can be modeled without study-local pulse synthesis glue.
- Expose time-dependent collapse operators or a structured stochastic-measurement replay path for readout studies.
- Add a built-in readout benchmarking utility that reports detector-limited, QND-constrained, and hardware-replayed frontiers from a common API.

## Status
COMPLETE