# Gray-Box Adaptive Control for cQED Systems

## Problem Class
OPT | ANA

## Motivation
Gate fidelity in dispersive cQED can degrade sharply when the learner model uses an
incorrect dispersive shift `chi`. This study tests whether a gray-box workflow,
consisting of a targeted multi-Fock Ramsey probe followed by model-based GRAPE on a
corrected learner Hamiltonian, is a practical alternative to both fixed-model control
and model-free black-box optimization.

The key question is whether this approach remains useful under realistic cQED
imperfections: finite shot budgets, readout confusion, slow parameter drift, and
Hamiltonian omissions such as `chi_higher` and cavity Kerr.

## Goals
1. Compare nominal, gray-box, perfect-knowledge, and black-box control strategies
   across chi mismatch levels from 0% to 40%.
2. Quantify how much of the perfect-model fidelity gray-box correction recovers at the
   30% mismatch operating point.
3. Test robustness to decoherence, readout confusion, probe-budget reduction, slow chi
   drift, and omitted higher-order Hamiltonian terms.
4. Determine whether the resulting workflow is scientifically meaningful and
   experimentally relevant for dispersive qubit-cavity control.
5. Produce a validated report and reproducibility checks for the consolidated study.

## Methods
- **Framework**: `cqed_sim` provides the dispersive transmon-cavity model, GRAPE
  optimization, and cross-model evaluation.
- **Models**: `DispersiveTransmonCavityModel` with `n_cav = 4` and `n_tr = 3` for the
  production study. Validation spot checks also test `n_cav = 5` and `n_tr = 4`.
- **Target**: simultaneous qubit X gate on the logical block
  `{|g,0>, |g,1>, |g,2>, |g,3>, |e,0>, |e,1>, |e,2>, |e,3>}`.
- **Optimization**: `GrapeSolver` with `UnitaryObjective` and `LeakagePenalty`, using
  16 slices of 10 ns each, qubit I/Q control, and multistart seeds `{2, 9, 14}`.
- **Gray-box probe**: analytical multi-Fock Ramsey model implemented in
  `scripts/probe_library.py`, because this specific calibration target is not yet in
  `cqed_sim.calibration_targets`.
- **Validation**: archived data checks plus fresh convergence spot checks in
  `scripts/validate_results.py`, with results written to `data/validation_summary.json`.

### Assumptions
- The control Hamiltonian is modeled in the dispersive rotating-wave regime.
- The gray-box control update feeds back the inferred `chi`; inferred `chi_higher` is
  retained for diagnosis but is not yet used in the production learner model.
- The truth model includes `chi_higher` and cavity Kerr, while the production learner
  omits Kerr and, by default, keeps `chi_higher = 0`.
- Noise is represented by aggregate `T1`, `Tphi`, and cavity loss rather than a
  hardware-calibrated noise model.

### Convergence Criteria
- Production results use multistart GRAPE and the best training-model fidelity.
- Validation requires stable conclusions under expanded multistart seeds and under
  enlarged Hilbert-space truncations.
- Probe-budget validation checks the expected `1/sqrt(N)` uncertainty scaling.

## Expected Outcomes
- At 0% mismatch, nominal, gray-box, and perfect-knowledge control should agree within
  numerical tolerance.
- At 30% mismatch, gray-box correction should recover nearly all of the perfect-model
  advantage over the nominal strategy.
- Readout confusion and finite shot counts should perturb the inferred `chi` only weakly
  over the tested ranges.
- Periodic recalibration should improve long-horizon performance under slow chi drift.

## Known Limitations
- The multi-Fock Ramsey probe is implemented analytically rather than as a first-class
  `cqed_sim` calibration target.
- The production gray-box loop corrects only `chi`; a fully consistent update of
  `chi_higher` and Kerr remains future work.
- The noise model is phenomenological rather than device-calibrated.
- The study targets a single logical task and one device parameter set.
- No one-to-one experimental benchmark is included in this repository pass.

## Status
COMPLETE

## Known API Gaps
1. `cqed_sim.calibration_targets` does not yet expose a multi-Fock dispersive Ramsey
   probe or matching inference helper.
2. The production gray-box workflow does not yet have a first-class reusable config
   object that encapsulates probe, inference, and re-optimization.

## Suggested Upstreaming
1. Add `run_chi_ramsey_dispersive` to `cqed_sim.calibration_targets`.
2. Add `infer_chi_dispersive` with uncertainty reporting.
3. Add a reusable `GrayBoxAdaptiveConfig`-style interface for probe -> infer -> update
   -> GRAPE loops.
