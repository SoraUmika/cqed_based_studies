# Open-System Deep Dive for Selective Qubit Rotation in Dispersive cQED

## Problem Class

ANA | OPT

## Motivation

The completed `sqr_pulse_waveform_design` study established the closed-system and simple
open-system behavior of selective qubit rotation (SQR) pulses, but its noise treatment was
intentionally lightweight: a single transmon `T1 = T2 = 20 us` model without storage thermal
occupation, readout-induced dephasing during pulsing, or Purcell-limited decay from the
readout resonator. Those omissions are now the dominant gap between the existing SQR report
and a device-relevant recommendation.

This study extends the flagship SQR results into a realistic noise regime using `cqed_sim`
first. The goal is to quantify which noise channels actually limit SQR fidelity on the
nominal device, determine whether the preferred operating window from the completed SQR study
survives under realistic decoherence, and identify when pulse re-optimization is justified.

## Goals

1. Quantify SQR fidelity degradation from explicit multilevel transmon relaxation
   `transmon_t1 = (T1_ge, T1_fe)` relative to the legacy single-`T1` approximation.
2. Map the effect of storage thermal occupation `n_th,storage in [0, 0.1]` on cphase-SQR
   and true-SQR fidelity for the representative pulse families from the completed SQR study.
3. Include Purcell decay inferred from the readout chain via `ReadoutChain.purcell_rate()`
   and `ReadoutChain.purcell_limited_t1()` for the nominal readout resonator
   `omega_r/2pi = 8.597 GHz`, `kappa_r/2pi = 2.4 MHz`.
4. Simulate the full three-mode system with `DispersiveReadoutTransmonStorageModel` to
   quantify readout-induced storage dephasing and spectator effects during SQR control.
5. Determine the gate-duration window that best balances selectivity and decoherence for
   each dominant noise channel.
6. Compare the best calibrated parametric SQR pulses against GRAPE-based references replayed
   through the noisy simulator path.

## Methods

### cqed_sim modules and functions

- `DispersiveReadoutTransmonStorageModel` for the three-mode qubit + storage + readout model.
- `FrameSpec` for lab-frame and rotating-frame consistency across the SQR and readout modes.
- `SimulationConfig`, `simulate_sequence`, `prepare_simulation`, and `simulate_batch` for
  single-run and sweep execution.
- `NoiseSpec(t1=..., transmon_t1=(...), tphi=..., tphi_storage=..., tphi_readout=...,
  kappa_storage=..., kappa_readout=..., nth_storage=..., nth_readout=...)` for Lindblad noise.
- `ReadoutChain`, `ReadoutResonator`, `PurcellFilter`, `ReadoutChain.purcell_rate()`, and
  `ReadoutChain.purcell_limited_t1()` for measurement-induced dephasing and Purcell analysis.
- `calibrate_sqr_gate()`, `load_or_calibrate_sqr_gate()`, `evaluate_sqr_gate_levels()`, and
  `extract_effective_qubit_unitary()` for reuse of the SQR calibration workflow.
- `build_control_problem_from_model()`, `GrapeSolver`, and
  `evaluate_control_with_simulator()` for closed-system control optimization followed by
  simulator-backed noisy replay of the archived control schedule.

### Reference inputs from existing studies

- Reuse the best closed-system pulse families and calibration settings from
  `studies/sqr_pulse_waveform_design` as the baseline controls to replay under realistic noise.
- Reuse the nominal device constants from `AGENTS.md` and the readout resonator parameters
  already used elsewhere in this repository.

### Coverage assessment

- Fully supported in `cqed_sim`: three-mode Hamiltonian construction, Lindblad noise replay,
  storage/readout pure dephasing, explicit multilevel transmon `T1`, Purcell-rate analysis,
  SQR calibration, and simulator-backed replay of optimized controls.
- Partially supported: direct optimization of an explicitly open-system, measurement-backaction-
  weighted SQR objective is not a single built-in `cqed_sim` objective. The implementation
  plan is therefore to optimize closed-system or parametric controls with `cqed_sim` and then
  replay them with `NoiseSpec` plus readout-chain-derived Purcell and dephasing terms. Any
  glue code needed for this comparison will be documented as a study-local extension.

## Assumptions

- Device constants start from the repository defaults: `omega_q/2pi = 6.150 GHz`,
  `omega_s/2pi = 5.241 GHz`, `omega_r/2pi = 8.597 GHz`, `chi_s/2pi = -2.84 MHz`,
  `chi_prime/2pi = -21 kHz`, `K_s/2pi = -28 kHz`.
- Readout resonator linewidth for Purcell analysis: `kappa_r/2pi = 2.4 MHz`.
- Representative storage decay for thermal-occupation sweeps: `kappa_storage/2pi = 10 kHz`.
- Transmon relaxation sweep: `T1_ge in [20, 50] us`, `T1_fe in [5, 15] us`.
- Storage thermal occupation sweep: `nth_storage in {0.00, 0.01, 0.02, 0.05, 0.10}`.
- Initial three-mode truncation: `n_storage = 8`, `n_readout = 6`, `n_tr = 3`; representative
  points must be checked for convergence by increasing storage and readout truncations.
- Initial three-mode implementation defaults to `chi_r = chi_s` and `chi_sr = 0` unless a
  more specific device-calibrated scenario is supplied during execution.
- Baseline pulse families and target metrics are inherited from the completed SQR study:
  cphase-SQR and true-SQR remain the primary comparison targets.
- Readout-induced dephasing is modeled through `cqed_sim` replay and readout-chain analysis,
  not through a custom stochastic-measurement solver unless the README is updated to document
  a concrete framework gap.
- Convergence target: representative fidelity metrics must remain stable to `5e-4` under
  increased Hilbert-space truncation and finer replay settings.
- Windows runs must reuse the repository runtime-compatibility shim before importing
  `qutip` or `cqed_sim`.

## Expected Outcomes

- A ranked breakdown of which open-system channels dominate SQR fidelity loss on the nominal
  device, with explicit comparison to the legacy single-`T1` approximation.
- A quantitative threshold for storage thermal occupation above which the preferred SQR pulse
  family or gate-duration recommendation changes materially.
- A Purcell-limited `T1` estimate for the nominal readout chain and a determination of whether
  it is negligible, comparable, or dominant relative to the intrinsic transmon `T1` budget.
- A recommended open-system operating window for SQR pulses that preserves at least
  `F_block >= 0.99` on representative targets if such a window exists; otherwise a documented
  blocker showing which physical limit prevents that regime.
- A clear answer on whether GRAPE replay under realistic noise meaningfully outperforms the
  calibrated parametric baselines once decoherence is included.

## Validation

- [x] Sanity checks - passed for the archived sweeps: legacy fidelity-to-ideal traces are monotone non-increasing in `chiT`, Purcell filtering improves the Purcell-limited `T1` across the detuning sweep, thermal-occupation scans are monotone non-increasing within tolerance, and mean three-mode backaction worsens monotonically with readout amplitude.
- [x] Two-mode convergence - passed at the representative multilevel point with `delta_t = 4.70e-4` and truncation delta `4.27e-9`.
- [x] Three-mode convergence - passed for the representative reduced-fidelity metric with `delta_t = 3.12e-4` and truncation delta `3.67e-6`.
- [x] GRAPE noisy replay convergence - passed after switching the archived-control replay to
  `evaluate_control_with_simulator(...)`. For the representative `|chi_s| T / 2pi = 2`
  control, the noisy target fidelity is `0.958318` at `512` replay substeps per slice and
  `0.958358` at `1024`, giving a convergence delta of `3.97e-5`.

## Status

COMPLETE
