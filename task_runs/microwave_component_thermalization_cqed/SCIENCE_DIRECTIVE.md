# Science Directive

Date: 2026-04-02
Study: `studies/microwave_component_thermalization_cqed`
Run: `task_runs/microwave_component_thermalization_cqed`

## Problem Class
ANA, DES

## Central Question
How does a hot microwave component, modeled as an effective thermal bath or hot auxiliary mode, map onto experimentally measurable qubit and cavity observables in a dispersive cQED system, and which observables most reliably diagnose dangerous thermal back-action before qubit coherence or readout performance becomes unacceptable?

## Hypotheses
1. The cavity thermal occupation is the most linear thermometer for the hottest regimes, but the qubit excited-state population provides the highest relative sensitivity at low effective occupation because it up-converts small thermal leakage into a rare-event signal.
2. Thermal-photon-induced dephasing becomes the dominant coherence limit before substantial qubit heating once the readout-mode occupation exceeds the perturbative low-occupation regime.
3. A weakly coupled hot auxiliary mode can degrade the qubit-cavity system well before it is directly populated if the auxiliary linewidth and detuning place it near a hybridization-enhanced thermal leakage window.
4. The intrinsic cQED response to a step in bath occupation is set by quantum dissipation rates such as `kappa` and auxiliary-mode linewidths; experimentally slower VTS traces must therefore come from external thermalization of the hardware component rather than the internal quantum subsystem alone.

## Model Hierarchy
### Thrust A: Minimal thermometer model
- Three-level transmon plus one readout/storage mode.
- Dispersive Hamiltonian with `omega_q`, `omega_r`, `alpha`, and `chi`.
- Hot bosonic bath through `NoiseSpec(kappa, nth)`, qubit bath kept cold.

### Thrust B: Thermal dephasing model
- Same qubit-readout model.
- Ramsey-like free evolution initialized in a qubit superposition.
- Thermal cavity bath retained during the idle period.

### Thrust C: Multimode hot-component model
- `UniversalCQEDModel` with one transmon and two bosonic modes: local readout mode plus hot auxiliary mode.
- Include qubit dispersive shifts to both modes, cross-Kerr if needed, and an exchange term between the bosonic modes to mimic cable-mediated hybridization or leakage.

### Thrust D: Transient thermal step
- Time-dependent sequence of steady-state or finite-time evolutions after a step in effective bath occupation.
- Compare instantaneous bath-step predictions to slower externally imposed thermal ramps as a future-model interface.

## Primary Observables
- Qubit excited-state population and non-computational transmon population.
- Cavity and auxiliary-mode occupations.
- Qubit coherence magnitude `|rho_ge|` during Ramsey-like evolution.
- Effective dephasing rate extracted from envelope fits.
- Dispersive readout proxies: cavity occupation, steady-state response contrast, and inferred qubit-state separation.
- Sensitivity metrics `d observable / dT` and `d observable / dn_th`.

## Quantitative Success Criteria
1. Produce monotonic calibration curves from temperature or `n_th` to at least three experimentally meaningful observables.
2. Identify at least one low-occupation and one high-occupation regime where the preferred thermometer changes.
3. Extract a dephasing-versus-occupation curve and compare it against a simple dispersive analytic estimate with explicit residual error.
4. Produce at least one multimode heat map showing a safe-versus-dangerous coupling or detuning boundary.
5. Demonstrate from simulation that the internal quantum response time is much shorter than plausible hardware thermalization times, supporting a modular quantum-plus-thermal modeling program.

## Validation Plan
- Zero-temperature limit: all thermal occupations and induced qubit heating must vanish numerically within solver tolerance.
- Weak-coupling limit: multimode back-action must disappear continuously as the auxiliary coupling approaches zero.
- Truncation convergence: repeat representative points with larger transmon and bosonic dimensions.
- Analytic checks: compare cavity occupation to Bose-Einstein expectation and compare low-occupation dephasing trends to the perturbative dispersive formula.

## Deliverables
- Shared study library plus execution scripts for all four thrusts.
- Machine-readable JSON artifacts for observables, fitted rates, and validation summaries.
- Figure set in `.png` and `.pdf`.
- Markdown memo for internal experimental use.
- LaTeX report and compiled PDF for the repo workflow.
- Reproducibility notebook under `scripts/reproducibility_notebook.ipynb`.
