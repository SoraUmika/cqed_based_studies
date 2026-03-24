# Selective Qubit Rotation (SQR) Gate Design in Dispersive cQED

## Problem Class
OPT | ANA | DES

## Motivation
Selective Qubit Rotation (SQR) is a key primitive for bosonic quantum error correction in circuit QED, enabling Fock-number-conditioned qubit rotations. This unified study consolidates three complementary investigations covering waveform optimization, simultaneous multi-branch control, and open-system performance under realistic noise, Purcell effects, and concurrent readout.

## Consolidated Studies
| # | Original Study | Focus | Key Contribution |
|---|----------------|-------|------------------|
| 1 | `sqr_pulse_waveform_design` | Closed-system waveform comparison and phase compilation | Cosine-squared best envelope; GRAPE ceiling; phase-compilation analysis |
| 2 | `simultaneous_multitone_sqr_design` | Simultaneous multi-branch SQR | Negative result: common Gaussian multitone fails to rotate target branches |
| 3 | `sqr_open_system_deep_dive` | Realistic noise, Purcell, three-mode model | Square-family best under noise; GRAPE retains advantage; concurrent readout penalty |

## Goals
1. Identify the optimal SQR waveform family across closed- and open-system models.
2. Determine whether simultaneous multi-branch SQR is viable with common multitone waveforms.
3. Quantify realistic noise budgets including multilevel relaxation, thermal occupation, and Purcell effects.
4. Compare parametric baselines against GRAPE under realistic noise.
5. Assess concurrent readout penalty on SQR gate fidelity.

## Methods
- `cqed_sim.DispersiveTransmonCavityModel` for two-mode replay
- `cqed_sim.DispersiveReadoutTransmonStorageModel` for three-mode model  
- `cqed_sim.GrapeSolver` and `build_control_problem_from_model` for GRAPE
- `cqed_sim.ReadoutChain`, `ReadoutResonator`, `PurcellFilter` for Purcell analysis
- `cqed_sim.UnitarySynthesizer` for gate decomposition
- `cqed_sim.NoiseSpec` for Lindblad dissipation

## Expected Outcomes
- Cosine-squared best closed-system envelope in practical operating regime (χT/(2π) ∼ 2–3)
- Square-family best under realistic noise (peak F = 0.944 near |χ|T/(2π) = 1)
- Simultaneous multitone SQR not viable (target angle response ∼ 10⁻⁴)
- GRAPE retains noisy-fidelity advantage over parametric baselines

## Known Limitations
- Open-system study uses phenomenological noise model, not calibrated to specific device
- GRAPE controls are archived/replayed, not directly optimized under noise
- Concurrent readout uses conservative χ_sr = 0 assumption
- No experimental validation data available

## Status
COMPLETE

## Source Data & Scripts
All original simulation scripts, data files, and figures are preserved in the component study directories listed above.
