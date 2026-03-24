# Hybrid Qubit-Cavity Control: Benchmark, Follow-Up Optimization, and Implementation Study

## Problem Class
`OPT` | `DES` | `ANA`

## Motivation
This unified study addresses the experiment-facing control question for a dispersive transmon + storage-cavity system: which control library is actually the best practical route to hybrid universal control once fidelity, leakage, duration, calibration burden, and waveform realism are all considered together?

The study now contains two layers:
1. A consolidated gate-library comparison on the logical block $\{|g,0\rangle, |g,1\rangle, |e,0\rangle, |e,1\rangle\}$.
2. A deeper follow-up optimization and implementation pass that exposes the actual optimized decompositions, optimized parameters, pulse waveforms, and optimizer-diagnostic evidence for the most important candidates.

The key follow-up result is that the previous depth-11 structured SQR synthesis was not the final answer. Adding one extra displacement gate lifts the best structured $U_\mathrm{target}$ synthesis from `F=0.9212, L=0.1050` to `F=0.9953, L=0.0065`, materially changing the study's conclusion about what is worth implementing next.

## Goals
1. Compare the leading gate families for local cavity control and hybrid entanglement on a common Fock-encoded 2x2 logical block.
2. Identify the best practical entangler and the best practical local logical primitive for the nominal cQED device point.
3. Synthesize a specific maximally entangling `U_target` with structured gate families and determine whether the best solution is near a meaningful optimum or just a shallow local minimum.
4. Expose the actual optimized solutions: gate sequence, parameter values, pulse durations, waveform structure, and implementation burden.
5. Test whether tuning selective pulse parameterizations and a constructive `selective pi + broadband pi` SNAP-like shortcut produce experimentally useful improvements.

## Methods
### Device parameters
- `omega_q/2pi = 6.150 GHz`
- `omega_c/2pi = 5.241 GHz`
- `alpha/2pi = -255 MHz`
- `chi/2pi = -2.84 MHz`
- `chi'/2pi = -21 kHz`
- `K/2pi = -28 kHz`
- `n_cav = 8`, `n_tr = 2` for synthesis and `n_tr = 3` for pulse-backed replay

### cqed_sim modules used
- `DispersiveTransmonCavityModel`, `FrameSpec`
- `Subspace`, `TargetUnitary`, `UnitarySynthesizer`, `GateSequence`
- `QubitRotation`, `Displacement`, `SNAP`, `SQR`, `ConditionalPhaseSQR`, `FreeEvolveCondPhase`
- `GrapeSolver`, `GrapeConfig`, `build_control_problem_from_model`
- `build_displacement_pulse`, `build_rotation_pulse`, `build_sqr_multitone_pulse`
- `simulate_sequence`, `leakage_metrics`, `subspace_unitary_fidelity`

### Follow-up workflow additions
- Multistart restart sweeps for the key structured candidates
- Local refinement of the best depth-11 solutions
- Depth-plus-one ansatz extension (`L1d`, `L2d`)
- Pulse-parameter grids for logical-window SNAP
- Numerical test of the constructive `selective pi + broadband pi` shortcut
- Automatic figure and LaTeX table generation via `scripts/run_followup_optimization.py`

## Expected Outcomes
- A benchmark-level recommendation for which family to use first in experiment
- Explicit optimized gate decompositions and parameter tables for the top candidates
- Waveform figures for the pulse-backed selective controls
- Stronger evidence about whether the reported solutions are stable optima within the tested ansatze
- A deployability-focused ranking that distinguishes ideal decomposition benchmarks from pulse-backed controls and GRAPE lower bounds

## Known Limitations
- All benchmark conclusions are still specific to the strict Fock `|0>, |1>` encoding
- `ConditionalPhaseSQR`, `ConditionalDisplacement`, and native exchange primitives are not yet fully waveform-backed in this study
- Structured synthesis still uses `n_tr = 2`; direct `|f>` leakage is not present during optimization
- The best structured sequence (`L1d`) is accurate but long (`7.3 us`)
- The follow-up shows depth-11 was not optimal, but it does not prove depth-12 is globally minimal
- Only one representative device parameter point is studied

## Suggested Upstreaming
- Add waveform-bridge coverage for `ConditionalPhaseSQR` and the native effective-interaction primitives so the full comparison can be replayed at equal implementation fidelity.

## Status
COMPLETE
