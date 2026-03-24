# Improvement Log: Thermal Noise Cavity Sensing in a cQED System

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- **cqed_sim lacks multi-bath and steady-state support**: The study uses QuTiP directly for multi-thermal-bath Lindblad dynamics and `qt.steadystate()`. These are not available in cqed_sim. Upstreaming `MultiThermalBath` and a steady-state solver would enable future sensing studies to stay within the framework.

## Recommended Improvements (P2)
- **Photon-number-resolved detection efficiency**: The current model assumes perfect state-resolved measurement. Including realistic detection efficiency and SPAM errors would give more experimentally relevant sensitivity estimates.
- **Bayesian inference for simultaneous κ_tot and n̄_ss extraction**: Currently the minimal observable set {n̄_ss, κ_tot} is identified analytically. A Bayesian time-resolved Ramsey approach could extract both simultaneously from a single measurement protocol.
- **Higher-order dispersive terms**: The model uses linear dispersive Hamiltonian H/ℏ = χ a†a |e⟩⟨e|. Including χ' and K corrections would quantify their impact on sensing accuracy at high photon numbers.
- **Non-Markovian bath corrections**: The model assumes flat-spectrum (Markov) baths. For structured environments (e.g., frequency-dependent attenuation), non-Markovian corrections could be important.

## Nice-to-Haves (P3)
- Cross-validation against experimental ring-down data from the cQED setup.
- Extension to multi-mode sensing (using storage and readout cavities simultaneously).
- Temperature-resolved sensitivity maps for practical thermometry applications.

## Open Questions
- At what bath coupling fraction κ_target/κ_tot does the sensing precision become limited by ancilla back-action rather than photon statistics?
- Can the Ramsey coherence generating function be inverted numerically for non-thermal (e.g., squeezed) bath states?
- How does finite qubit population (thermal excitation of the ancilla) modify the identifiability structure?

## What Was Tried and Did Not Work
- Nothing recorded as failed—all 51 validation tests pass and the study is marked COMPLETE. The approach is comprehensive within its stated assumptions.

## Compute & Resource Notes
- 51 validation tests all pass. Report compiled to 11-page PDF.
- Fock-space truncation: N=30 default. The geometric P_n distribution ensures rapid convergence.
- Steady-state solving with QuTiP is fast for single-mode problems at N=30.
