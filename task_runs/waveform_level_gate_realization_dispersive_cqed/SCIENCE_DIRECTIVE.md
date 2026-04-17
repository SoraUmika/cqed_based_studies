# Science Directive

Date: 2026-04-02
Study: `studies/waveform_level_gate_realization_dispersive_cqed`
Run: `task_runs/waveform_level_gate_realization_dispersive_cqed`

## Problem Class
ANA, DES, OPT

## Central Question
Which control strategy most accurately approximates an unconditional cavity displacement
\[
U_{\mathrm{target}} = I_q \otimes D(\alpha_{\mathrm{target}})
\]
in a realistic dispersive qubit-cavity system once the qubit may begin in `|g>`, `|e>`, or a superposition?

## Hypotheses
1. A naive single-tone cavity pulse calibrated in the `|g>` branch will fail once `|chi| T` becomes appreciable, producing both branch-dependent displacements and residual qubit-cavity entanglement for superposition inputs.
2. Very short broadband single-tone pulses will reduce branch mismatch but will become vulnerable to hardware-filter distortion unless the command waveform is precompensated.
3. A physically interpretable two-tone pulse addressing both branch resonances should strongly suppress vacuum-branch mismatch and outperform naive single-tone control at matched duration.
4. A practical echoed displacement using vacuum-calibrated `pi` pulses will not fully realize the ideal toggling-frame cancellation picture because the inserted qubit pulse itself becomes manifold dependent.
5. Hardware-aware optimal control will provide the best mean fidelity on the chosen state-test set, but may sacrifice interpretability relative to the two-tone solution.

## Experimental Design
### Hamiltonian hierarchy
- Minimal dispersive model with `chi`
- Higher-order model adding `chi'`
- Full model adding self-Kerr `K`

### Protocol families
- Naive single-tone square pulse
- Fast single-tone square, Gaussian, and cosine pulses
- Two-tone branch-compensated displacement
- Echoed displacement `D(alpha/2) -> X_pi -> D(alpha/2) -> X_pi`
- Hardware-aware GRAPE waveform

### Required metrics
- Branch displacement mismatch `delta_alpha = |alpha_g - alpha_e|`
- Superposition-state entanglement entropy and branch overlap
- Fidelity to the ideal displaced state on the explicit state-test set
- Wigner-function comparison for representative failing and successful cases
- Protocol complexity and experimental practicality

### Validation states
- Qubit: `|g>`, `|e>`, `(|g> + |e>)/sqrt(2)`, `(|g> + i|e>)/sqrt(2)`
- Cavity: `|0>`, `|1>`, `|2>`, `|3>`, and a modest coherent state

## Quantitative Success Criteria
1. Baseline failure must be visible and numerically explicit for at least one representative long-duration naive pulse.
2. At least one simple interpretable protocol must reduce `delta_alpha` by more than an order of magnitude relative to the naive `80 ns`, `alpha = 1` square pulse.
3. The final study must report a protocol ranking using the same state-test fidelity metric for every family.
4. The report must distinguish vacuum-branch matching from broader state-set fidelity rather than conflating them.
5. The final recommendation must state both the best overall protocol and the best experimentally interpretable protocol.

## Deliverables
- Updated study README and improvement log
- Machine-readable unconditional-displacement artifacts
- Publication-style report and compiled PDF
- Reproducibility notebook using the new `unconditional_*` artifacts
- Execution summary and review handoff files
