# Literature-Informed Selective Pulse Primitives for Dispersive cQED

## Problem Class
OPT | DES | ANA

## Motivation
The consolidated SQR and hybrid-control studies rely on SNAP-like number-selective phases, conditional-phase-relaxed selective qubit rotations, and native dispersive phase accumulation, but the repository does not yet contain one study that turns the literature pulse prescriptions for those primitives into a single optimized, noise-aware reference dataset. This study fills that gap using the typical dispersive cQED parameter set already adopted across the repository.

## Goals
1. Extract literature-backed pulse definitions for number-selective qubit control relevant to SNAP, SQR, and chi-wait conditional phase.
2. Implement those pulse families through the normal `cqed_sim` pulse compilation and simulation path.
3. Numerically optimize representative primitive instances under typical dispersive cQED parameters.
4. Replay the optimized primitives under realistic Lindblad noise and identify the best practical operating points.
5. Produce reusable figures, machine-readable data, and a report that downstream studies can cite.

## Methods
- `cqed_sim.DispersiveTransmonCavityModel` for the qubit-storage system.
- `cqed_sim.FrameSpec`, `SequenceCompiler`, `SimulationConfig`, `prepare_simulation`, and `simulate_sequence` for pulse-level replay.
- `cqed_sim.Pulse` and standard envelope helpers for analytic waveform definitions.
- `cqed_sim.NoiseSpec` for qubit relaxation, qubit dephasing, storage loss, and storage thermal occupation.
- `scipy.optimize` for low-dimensional numerical optimization of pulse duration, envelope width, and phase parameters.

## Expected Outcomes
- A literature-backed selective-Gaussian baseline for SQR and SNAP that is directly runnable in this repository.
- A comparison between Gaussian, cosine-squared, and truncated-flat-top selective pulses at the same cQED operating point.
- Noise-aware recommendations for when short pulses outperform more selective pulses because of decoherence.
- A reusable reference showing which primitive should be preferred for:
  1. branch-selective qubit rotation,
  2. single-branch cavity phase shifts,
  3. wait-based conditional phase accumulation.

## Known Limitations
- `cqed_sim` does not currently provide first-class pulse builders or waveform-bridge support for `SNAP` or `ConditionalPhaseSQR`; study-local pulse constructors are therefore required.
- The study focuses on a representative logical window and typical device parameters, not a device-calibrated dataset from one experiment.
- The main optimization loop is low dimensional and physics-informed; it is not a full open-system optimal-control solve.
- The baseline model is qubit plus storage cavity only; readout-assisted effects are treated through noise, not an explicit readout resonator.

## Status
COMPLETE

## Assumptions
- Typical device parameters follow the repository defaults: `omega_q / 2pi = 6.150 GHz`, `omega_c / 2pi = 5.241 GHz`, `alpha / 2pi = -255 MHz`, `chi / 2pi = -2.84 MHz`, `chi' / 2pi = -21 kHz`, and `K / 2pi = -28 kHz`.
- Number selectivity requires pulse bandwidth smaller than the branch spacing set by `|chi|`, so the study optimizes over durations on the order of `1 / |chi|` or longer.
- Representative noise uses physically typical dispersive-cQED values in the same regime as the consolidated studies: qubit `T1` and `Tphi` in the tens of microseconds, storage decay in the few-kHz to tens-of-kHz regime, and nonzero storage thermal occupation.
- Primitive benchmarks are run on a small logical window first and then checked for truncation and time-step convergence.

## Suggested Upstreaming
- Promote the study-local selective SNAP builder into `cqed_sim.pulses.builders`.
- Add waveform-bridge support for `SNAP` and `ConditionalPhaseSQR`.
- Add a standard library helper for literature-style number-selective Gaussian and flat-top Gaussian pulses with explicit sigma and cutoff conventions.
