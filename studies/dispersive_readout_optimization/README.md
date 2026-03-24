# Dispersive Readout Pulse Optimization

## Problem Class
OPT | ANA | DES

## Motivation
Dispersive readout is the canonical measurement mechanism for superconducting qubits in circuit QED. Optimizing the readout pulse shape requires balancing signal-to-noise ratio, cavity emptying speed, QND character, and hardware realizability. This unified study consolidates four progressive investigations into a single authoritative reference covering the full model-refinement chain from linear dispersive theory through hardware-realistic nonlinear QND modeling and leakage/ionization mapping.

## Consolidated Studies
This study merges the following component investigations:

| # | Original Study | Model Level | Key Contribution |
|---|----------------|-------------|------------------|
| 1 | `readout_pulse_optimization` | Linear dispersive | Baseline GRAPE vs. analytical families; midpoint drive optimality |
| 2 | `procedural_readout_pulse_sequence_optimization` | Multilevel dispersive | Structured procedural families; analytic nulling tails |
| 3 | `nonlinear_qnd_hardware_realistic_readout` | Hardware-realistic + QND | Ring-hold as practical default; transport/QND decomposition |
| 4 | `measurement_induced_leakage_ionization_modeling` | Extended transmon (5-level) | Leakage onset mapping; ionization regime characterization |

## Goals
1. Identify the optimal readout pulse family across a progressive model hierarchy (linear → multilevel → hardware-realistic → leakage-aware).
2. Quantify the performance of structured procedural pulses vs. GRAPE and square baselines.
3. Map the hardware-realistic QND boundary and characterize the separation between leakage and ionization regimes.
4. Provide actionable experiment-facing pulse recommendations with quantified tradeoffs.

## Methods
- `cqed_sim.DispersiveTransmonCavityModel` for multilevel replay
- `cqed_sim.SequenceCompiler` with hardware transport model
- `cqed_sim.NoiseSpec` for Lindblad dissipation  
- Study-local: linear dispersive propagation, matched-filter SNR, occupancy-activated disturbance model, pulse family parametrizations

## Expected Outcomes
- Ring-hold family as best practical default at 240 ns under full hardware/QND model
- Leakage onset at ε/(2π) ≈ 5.5 MHz for 240 ns pulses
- Procedural families outperform square by ~10% in fidelity at intermediate durations
- QND breakdown dominates over high-manifold ionization

## Known Limitations
- Phenomenological disturbance model (not calibrated to experiment)
- Stochastic replay numerically unstable at highest drive powers
- No fridge data available for experimental validation
- Detector-limited fidelity proxy saturates, obscuring assignment-level mitigation tradeoffs

## Status
COMPLETE

## Source Data & Scripts
All original simulation scripts, data files, and figures are preserved in the component study directories listed above. This unified report synthesizes their results into a single narrative.
