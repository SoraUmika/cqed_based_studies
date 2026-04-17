# Improvement Log: Thermalization of Microwave Components in cQED

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- **[P1 | HIGH] No native steady-state helper in the public `cqed_sim` workflow**: The package exposes thermal collapse operators and transient Lindblad simulation, but not a dedicated steady-state solver. This study therefore uses `cqed_sim` Hamiltonians and noise operators with a small QuTiP `steadystate` wrapper. Upstreaming a first-class helper would remove this local bridge.
- **[P1 | HIGH] Macroscopic thermal transport is not captured**: The quantum model only sees an effective bath occupation or auxiliary-mode temperature. It cannot predict real component temperature profiles, boundary resistances, or VTS heat-flow delays without a separate thermal model.

## Recommended Improvements (P2)
- **[P2 | MEDIUM] Replace qualitative auxiliary-mode scales with measured mode spectra**: Once the exact SNAILmon-cable or device-mode data are available, rerun the multimode maps with those frequencies, linewidths, and couplings.
- **[P2 | MEDIUM] Add a full readout-chain performance layer**: Convert the current dispersive-response proxies into integrated assignment-fidelity forecasts including amplifier noise and finite integration time.

## Nice-to-Haves (P3)
- **[P3 | LOW] Add non-Markovian colored-noise surrogates**: A finite-bandwidth or telegraph-noise bath model could test whether simple Lindblad predictions understate dephasing in some cable or dielectric scenarios.
- **[P3 | LOW] Couple the transient study to a lumped thermal RC model**: This would allow direct fits of VTS heating traces using a small set of thermal time constants.

## Open Questions
- How distinguishable are hot-readout-induced qubit heating and direct qubit-bath heating using only standard spectroscopy and thermometry traces?
- Does a weakly coupled auxiliary mode create a unique fingerprint in the joint trend of qubit excitation and Ramsey decay that could separate cable-mode leakage from simple broadband thermalization?

## What Was Tried and Did Not Work
- **Process-based `joblib` parallelism inside Jupyter on Windows**: The original multimode map used the default `loky` process backend. The study still ran, but notebook execution emitted noisy resource-tracker warnings even though the calculations succeeded. Threaded parallelism was adopted instead because the wall-clock time stayed acceptable and the notebook became cleanly reproducible.

## Compute & Resource Notes
- Installed `joblib` with `python -m pip install --user joblib` to parallelize the multimode parameter grid.
- Thermometer sweep runtime: about 11 s.
- Dephasing sweep runtime: about 20 s.
- Multimode heating-map runtime: about 49 s using threaded `joblib` parallelism.
- Transient temperature-step runtime: about 3 s.
- Full reproducibility notebook execution time: about 112 s.
- No GPU acceleration was needed for the present grid sizes.

## Resolved
- **Windows notebook parallel-backend noise**: Switched the multimode sweep from `joblib`'s default process backend to threaded parallelism, which preserved acceptable runtime while avoiding notebook-time `loky` resource-tracker warnings.
