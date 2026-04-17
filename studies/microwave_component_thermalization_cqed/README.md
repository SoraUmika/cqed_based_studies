# Thermalization of Microwave Components in cQED

## Problem Class
ANA, DES

## Motivation
Hot microwave components can inject thermal photons and excess dissipation into a superconducting cQED stack even when the mixing chamber remains cold. This study uses a dispersive transmon-plus-readout architecture as a quantum probe of effective bath occupation to determine which observables are most sensitive to hot attenuators, cables, dielectrics, filters, and auxiliary modes, and to identify experimentally safe versus dangerous operating regimes for coherence and readout.

## Goals
1. Build a quantitatively explicit map from effective hot-component occupation to measurable qubit and cavity observables.
2. Identify which observable is the most sensitive thermometer across realistic cQED parameter regimes.
3. Quantify thermal-photon-induced dephasing, population heating, and readout degradation as functions of bath occupation and component temperature.
4. Determine how a hot auxiliary mode or cable-like resonance back-acts on a nominally cold qubit-cavity system.
5. Separate intrinsic quantum response times from the slower macroscopic thermalization of real hardware components and state how future thermal-transport models should be coupled to the cQED layer.

## Methods
- `cqed_sim.core.DispersiveTransmonCavityModel` for the minimal qubit-plus-readout model.
- `cqed_sim.core.UniversalCQEDModel` and multimode coupling specs for hot auxiliary-mode studies.
- `cqed_sim.sequence.SequenceCompiler` with idle compiled sequences to propagate free evolution under the drift Hamiltonian.
- `cqed_sim.sim.{SimulationConfig, NoiseSpec, simulate_sequence}` for Lindblad dynamics with thermal bosonic baths.
- `cqed_sim.sim.noise.collapse_operators` for validated thermal collapse channels.
- `cqed_sim.observables` plus reduced-state post-processing for qubit populations, cavity occupation, coherence, and Fock distributions.
- A small local QuTiP `steadystate` bridge built on top of `cqed_sim` Hamiltonians and collapse operators, because the public package does not yet expose a dedicated steady-state helper.

## Analytic Preliminary
The first-pass analytic picture is the standard dispersive thermal-photon model. For a bosonic mode of angular frequency `omega` coupled to a bath at temperature `T`, the effective occupation is
\[
\bar{n}_{\mathrm{th}}(T,\omega)=\frac{1}{e^{\hbar \omega /(k_B T)}-1}.
\]
In the weak-excitation dispersive regime, the transmon sees a photon-number-dependent transition
\[
\omega_{ge}(n)=\omega_{ge}(0)+\chi n,
\]
so cavity photon-number fluctuations convert into qubit dephasing. In the simplest bad-cavity thermal limit, one expects
\[
\Gamma_{\varphi,\mathrm{th}} \sim \frac{4 \chi^2}{\kappa}\bar{n}_{\mathrm{th}}(\bar{n}_{\mathrm{th}}+1)
\]
when `|chi| << kappa`, with crossover corrections once `|chi|/kappa` is not perturbatively small. The qubit excited-state population should remain exponentially suppressed if its own bath is cold, but can rise indirectly through residual dressing or multimode leakage from a hot resonator. These expectations motivate the numerical thrusts and provide limiting-case checks, but they do not capture multimode leakage, finite truncation, or explicit readout observables.

## cqed_sim Gap Analysis
| Functionality | Needed? | Available in cqed_sim? | Plan |
|---|---|---|---|
| Dispersive qubit-cavity Hamiltonian | Yes | Yes | Use `DispersiveTransmonCavityModel` directly |
| Multimode dispersive plus exchange model | Yes | Yes | Use `UniversalCQEDModel` and `DispersiveReadoutTransmonStorageModel` |
| Thermal Lindblad bosonic bath with `nth` | Yes | Yes | Use `NoiseSpec` and `collapse_operators` |
| Idle/free evolution under drift | Yes | Yes | Use empty compiled sequences with explicit `t_end` |
| Direct steady-state solver helper | Yes | Partial | Use `cqed_sim` Hamiltonian and collapse operators, then call QuTiP `steadystate` locally |
| Time-dependent hardware thermal transport | Yes | No | Represent only the induced quantum bath occupation versus time, not macroscopic heat diffusion |

## Assumptions
- The central cQED subsystem remains in the dispersive regime throughout the study: `|chi|`, auxiliary couplings, and linewidths are small compared with qubit-cavity detunings.
- A three-level transmon is sufficient for the baseline thermometer and coherence studies; truncation convergence checks test whether larger cutoffs materially change the main observables.
- The hot component is modeled as an effective bosonic bath occupation or hot auxiliary mode rather than a microscopic phonon or boundary-resistance model.
- Lindblad dynamics are adequate for the slow thermal channels considered here; non-Markovian cable or dielectric physics is out of scope for this first pass.
- Exact hardware parameters from the project slides are incomplete in the present workspace, so the baseline study uses experimentally reasonable sweeps anchored to standard cQED scales and states clearly where conclusions depend on those choices.

## Compute & Resource Strategy
The study fits within a single-machine CPU workflow. The heaviest step is the steady-state temperature and multimode sweep, so the implementation uses the smallest Hilbert truncations consistent with convergence and calls QuTiP `steadystate` directly instead of integrating to very long times.

Realized runtimes on this workstation were:
- Thermometer sweep: about 11 s
- Dephasing sweep: about 20 s
- Multimode heating maps: about 49 s
- Transient temperature-step sweep: about 3 s
- Full reproducibility notebook execution: about 112 s

`joblib` threaded parallelism was used for the multimode grid so the notebook remains stable on Windows/Jupyter without `loky` resource-tracker warnings. GPU acceleration was not necessary for the present sweep sizes.

## Expected Outcomes
- Calibration-style curves from bath temperature or occupation to qubit excited-state population, cavity thermal occupation, spectroscopy shift, and readout contrast.
- Quantified dephasing and coherence-time limits versus thermal occupation and component temperature.
- Multimode regime maps that identify safe coupling and detuning windows for hot auxiliary structures.
- Transient response estimates that cleanly separate quantum equilibration from real hardware thermalization.
- Actionable recommendations for which next experiments best constrain attenuator, cable, Teflon, and multimode thermalization.

## Known Limitations
- The study does not model microscopic thermal transport, Kapitza resistance, or distributed temperature gradients inside a component.
- Auxiliary-mode parameter ranges are only qualitatively motivated because the exact slide values are not available in the current workspace.
- The readout analysis uses dispersive-response surrogates rather than a full amplifier-chain optimization study.
- Any inference from bath temperature to actual component temperature depends on an assumed coupling between the component and the relevant electromagnetic mode.

## Validation
- [x] Sanity checks
- [x] Convergence
- [x] Literature comparison (if applicable)

## Suggested Upstreaming
- Add a public `cqed_sim` steady-state helper that accepts a model and `NoiseSpec`.
- Add a documented idle-evolution utility for thermalization and coherence studies.
- Add calibration helpers for cavity-thermal Ramsey and hot-bath qubit thermometry.

## Status
COMPLETE
