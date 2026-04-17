---
name: cqed-sim-lookup
description: "Look up cqed_sim API before writing simulation code. Use when: planning a simulation, checking if functionality exists in cqed_sim, finding the right model/function/class, writing pulse sequences, configuring simulation parameters, or mapping a first-principles analytic model onto cqed_sim. Provides the API reference and usage guidance."
argument-hint: "What you need to simulate, e.g. 'dispersive readout with noise'"
---

# cqed_sim API Lookup

## When to Use

- Before writing **any** simulation code (mandatory per AGENTS.md Critical Rule #1)
- Checking whether a feature exists in cqed_sim before writing standalone code
- Finding the correct class, function, or parameter name
- Understanding pulse construction, simulation config, or analysis tools

## API Reference Location

The canonical API reference is maintained alongside the cqed_sim source code:

- **Local file:** `C:\Users\dazzl\Box\Shyam Shankar Quantum Circuits Group\Users\Users_JianJun\cQED_simulation\API_REFERENCE.md`
- **GitHub:** [API_REFERENCE.md](https://github.com/SoraUmika/qubox_cQEDsim/blob/main/API_REFERENCE.md)

**Always read the local file** — it reflects the installed version and avoids network dependencies.

## Procedure

### 1. Sketch the First-Principles Model

Before searching the API, write down the minimal first-principles model you are trying to simulate:

- What Hamiltonian, equations of motion, or conservation laws define the problem?
- Which controlled approximations are you planning to use (RWA, dispersive, perturbative, adiabatic elimination, etc.)?
- What validity conditions must hold for those approximations?

This step tells you what capabilities the API must support and prevents API-driven rather than physics-driven study design.

### 2. Read the API Reference

Read the local `API_REFERENCE.md` file. If looking for specific functionality, search for keywords related to:

- **Models:** `DispersiveTransmonCavityModel`, `DispersiveReadoutTransmonStorageModel`, `UniversalCQEDModel`
- **Pulses:** `Pulse`, envelopes (`GaussianEnvelope`, `FlatTopEnvelope`, etc.), `PulseBuilder`
- **Simulation:** `simulate_sequence`, `SimulationConfig`, `SimulationResult`, `SimulationSession`
- **Noise:** `NoiseSpec`, Lindblad operators, T1/T2 parameters
- **Analysis:** `cqed_sim.analysis` module
- **Gates:** `cqed_sim.gates`, ideal gate operators
- **Tomography:** `cqed_sim.tomo`
- **Optimization:** `UnitarySynthesizer`, optimal control
- **Calibration:** `cqed_sim.calibration`, `cqed_sim.calibration_targets`
- **Plotting:** `cqed_sim.plotting`

### 3. Assess Coverage

After reading, classify the required functionality:

| Coverage | Action |
|----------|--------|
| **Fully supported** | Use cqed_sim directly. Do not write standalone code. |
| **Partially supported** | Use cqed_sim as the foundation and extend. Document the gap in the study README. |
| **Not supported** | Write standalone code. Document the gap and add a `## Suggested Upstreaming` section to the README. |

### 4. Document the Plan

In the study README `## Methods` and `## Analytic Preliminary` sections, list:
- The first-principles model being implemented
- Any controlled approximations and their validity conditions
- Exact classes and functions to be used
- Parameter configurations
- Any gaps or extensions needed

## Physics Conventions

For sign definitions, Hamiltonian algebra, and rotating-frame conventions, consult:

- **Local:** `C:\Users\dazzl\Box\Shyam Shankar Quantum Circuits Group\Users\Users_JianJun\cQED_simulation\physics_and_conventions\physics_conventions_report.tex`
- **GitHub:** [physics_conventions_report.tex](https://github.com/SoraUmika/qubox_cQEDsim/blob/main/physics_and_conventions/physics_conventions_report.tex)

## Key Reminders

- **Units:** cqed_sim is unit-coherent. The recommended convention is **rad/s** and **seconds**. Constructors with `_hz` or `_ns` suffixes accept those specific units.
- **Dependencies:** NumPy ≥ 1.24, SciPy ≥ 1.10, QuTiP ≥ 5.0. Optional: JAX for dense-matrix backend.
- **Backends:** QuTiP (default) or JAX. Select via `SimulationConfig.backend`.
