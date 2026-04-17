# Unconditional Cavity Displacement in Dispersive cQED

## Problem Class
ANA, DES, OPT

## Motivation
An ideal cavity drive is often modeled as the unconditional gate
\[
U_{\mathrm{target}} = I_q \otimes D(\alpha),
\]
but this approximation breaks down in the dispersive regime once the qubit is not pinned in `|g>`. The cavity resonance shifts by the qubit-state-dependent pull `chi`, and realistic control is further perturbed by `chi'`, cavity self-Kerr, finite-duration envelopes, bandwidth limits, and hardware transfer effects. The central question in this study is therefore not whether a cavity pulse displaces the oscillator in one branch, but whether it can realize nearly the same displacement in both qubit branches while avoiding residual qubit-cavity entanglement for superposition inputs.

This study benchmarks five protocol families under a common `cqed_sim` workflow:

1. Naive single-tone displacement
2. Fast broadband single-tone displacement
3. Two-tone branch-compensated displacement, interpreted as the minimal multiplex drive for the two-branch dispersive problem
4. Echoed displacement with inserted qubit `pi` pulses
5. Hardware-aware optimal-control synthesis

The goal is a quantitative recommendation for when unconditional displacement is reliable, when it fails, and which control strategy offers the best balance of fidelity, simplicity, and experimental interpretability.

## Goals
1. Demonstrate quantitatively why a cavity pulse calibrated in the `|g>` branch fails to realize `I_q \otimes D(alpha)` once the qubit starts in `|e>` or a superposition.
2. Benchmark naive, fast, two-tone, echoed, and numerically optimized protocols against a common target displacement amplitude and a common test-state set, while clarifying when the two-tone solution should be regarded as a simple multiplex drive and when richer multitone structure may be needed.
3. Report branch mismatch `delta_alpha`, superposition-state entanglement, state fidelities on a representative cavity basis, and Wigner-function agreement for successful and failing cases.
4. Separate what is already true in the minimal dispersive model from what changes after adding `chi'`, self-Kerr, and hardware filtering.
5. Produce a clear final recommendation for practical unconditional displacement in this parameter regime.

## Methods
- `DispersiveTransmonCavityModel` for the dispersive transmon-cavity Hamiltonian with `chi`, `chi_higher`, and `kerr`
- `FrameSpec` for the bare-frequency rotating frame used throughout the study
- `Pulse`, shaped displacement envelopes, and `SequenceCompiler` for waveform construction and propagation
- `SimulationConfig`, `prepare_simulation`, and shared simulation-session helpers in [`common.py`](scripts/common.py)
- `displacement` / `displacement_op` for ideal cavity translations
- Branch-resolved cavity-frequency extraction through `cavity_branch_transition_frequency`
- Two-tone calibration through a linear response solve on the simulated branch displacements, interpreted experimentally as a two-component multiplex cavity drive
- A multiplex follow-up using Fourier-compressed multicarrier approximations, segmented branch-resonant fits, and a jointly optimized shaped-two-tone family, all benchmarked on the same broad state-test set
- Hardware-aware GRAPE using `build_control_problem_from_model`, `HeldSampleParameterization`, `GrapeSolver`, `FirstOrderLowPassHardwareMap`, `SmoothIQRadiusLimitHardwareMap`, and `BoundaryWindowHardwareMap`
- QuTiP diagnostics for cavity reductions, trace distance, entanglement entropy, and Wigner functions
- Main driver script: [`unconditional_displacement_study.py`](scripts/unconditional_displacement_study.py)
- Multiplex follow-up script: [`multiplex_displacement_followup.py`](scripts/multiplex_displacement_followup.py)

## Analytic Preliminary
In the minimal dispersive model,
\[
H/\hbar = \Delta_c a^\dagger a + \frac{\chi}{2} a^\dagger a \sigma_z + \epsilon(t) a^\dagger + \epsilon^*(t)a,
\]
the cavity sees different effective detunings in the two qubit manifolds. For a qubit eigenstate `|b>` with branch detuning `Delta_b`, the final coherent displacement is
\[
\alpha_b(T) = -i \int_0^T \epsilon(t) e^{-i \Delta_b (T-t)} dt.
\]
For a square pulse with constant complex amplitude `epsilon`, this reduces to
\[
\alpha_b(T) = \epsilon \frac{1 - e^{-i \Delta_b T}}{\Delta_b},
\]
with the unconditional limit `alpha_b(T) -> -i epsilon T` recovered only when `|Delta_b| T << 1`. This predicts the basic failure mode of a naive pulse: one qubit branch is close to resonance while the other acquires a rotated and reduced displacement.

For a qubit superposition, the final state takes the form
\[
\frac{|g\rangle |\psi_g\rangle + |e\rangle |\psi_e\rangle}{\sqrt{2}},
\]
so the residual entanglement is controlled by the branch overlap `|<psi_g|psi_e>|`. A protocol can therefore look acceptable on `|g>` and `|e>` separately while still failing badly on `(|g> + |e>)/sqrt(2)`.

The echoed sequence `D_1 -> X_pi -> D_2 -> X_pi` cancels the first-order dispersive phase only if the inserted `pi` pulses are effectively instantaneous and uniform across the populated cavity manifolds. In a realistic dispersive cavity, however, the qubit `pi` pulse inherits number-dependent detuning during the sequence, so first-order echo cancellation is not guaranteed once the cavity has been displaced.

The two-tone strategy is analytically motivated by driving both branch resonances and solving for complex tone weights `w_g` and `w_e` such that the branch response matrix satisfies
\[
R \begin{bmatrix} w_g \\ w_e \end{bmatrix}
= \alpha_{\mathrm{target}} \begin{bmatrix} 1 \\ 1 \end{bmatrix}.
\]
This makes the compensation mechanism physically interpretable even before the numerical sweep. In the present two-branch problem, the same construction can be viewed as the minimal multiplex drive: two complex spectral weights are sufficient to equalize two branch responses. A richer multitone or frequency-comb drive would only become necessary if the control objective must also track higher occupied cavity manifolds, coherent-state support, or hardware-distorted spectral structure.

## cqed_sim Gap Analysis
| Functionality | Needed? | Available in cqed_sim? | Plan |
|---|---|---|---|
| Dispersive Hamiltonian with `chi`, `chi'`, and self-Kerr | Yes | Yes | Use `DispersiveTransmonCavityModel` directly |
| Shaped cavity displacement pulses | Yes | Yes | Build with shared pulse helpers in the study |
| Branch-resolved cavity frequencies | Yes | Partial | Compute using existing transition-frequency helpers plus a small shared helper |
| Time-domain simulation with stored final states | Yes | Yes | Use the shared simulation-session wrapper |
| Two-tone calibrated displacement | Yes | Partial | Reuse `cqed_sim` propagation and solve the 2x2 complex response system locally |
| Hardware-aware optimal control | Yes | Yes | Use GRAPE plus the available hardware maps |
| Wigner and reduced-state diagnostics | Yes | Yes | Use QuTiP-based diagnostics already present in the repo |

No blocking simulator gap was found for the requested protocol comparison. The only local additions were reusable analysis helpers and pulse-construction utilities layered on top of `cqed_sim`.

## Suggested Upstreaming
- The shaped-displacement helpers, branch-frequency helper, and state-distance utilities added to [`common.py`](scripts/common.py) are reusable across future bosonic-control studies.
- The hardware-aware waveform export path used for the optimal-control benchmark could be factored into a more general reusable study utility for constrained oscillator-control design.

## Assumptions
- Bare-frequency rotating frame at the qubit and cavity frequencies
- Closed-system unitary evolution; no `T1`, `Tphi`, or cavity-loss channels are included
- Default device parameters:
  - `omega_q / 2pi = 6.150 GHz`
  - `omega_c / 2pi = 5.241 GHz`
  - `alpha / 2pi = -255 MHz`
  - `chi / 2pi = -2.84 MHz`
  - `chi' / 2pi = -21 kHz`
  - `K / 2pi = -28 kHz`
- Default truncation:
  - transmon dimension `n_tr = 3`
  - cavity dimension `n_cav = 15`
  - optimal-control benchmark uses `n_cav = 12` to keep the constrained solve lightweight
- Default propagation step `dt = 0.5 ns`
- Logical cavity benchmark states include `|0>`, `|1>`, `|2>`, `|3>`, and a modest coherent state
- Superposition validation uses `(|g> + |e>)/sqrt(2)` and `(|g> + i|e>)/sqrt(2)`
- The main summary metric is the mean state fidelity on the explicit state test set saved in the artifacts

## Compute & Resource Strategy
- Upfront expectation: this study is dominated by many short closed-system propagations rather than by one long simulation.
- Chosen strategy: keep the state-space modest, reuse compiled sessions aggressively, and reserve the optimal-control benchmark for only two candidate durations.
- No GPU backend was required for the present parameter sizes.
- Realized wall-clock time for the main unconditional-displacement driver: about `10.1 s` on one CPU core.
- The most expensive step is the hardware-aware optimal-control solve, but it remained fast enough that additional parallelization was unnecessary for this study size.

## Expected Outcomes
- A quantitative demonstration that the naive single-tone pulse is not an unconditional displacement once `|chi| T` becomes appreciable
- A protocol-by-protocol comparison table with fidelity, branch mismatch, entanglement, duration, and implementation complexity
- A recommendation for the best experimentally interpretable protocol and the best overall protocol under the stated constraints
- A regime statement identifying when an experimentalist should stop modeling a cavity pulse as `I_q \otimes D(alpha)` and instead treat it as a conditional hybrid operation, plus guidance on when the minimal two-tone multiplex picture is likely insufficient

## Known Limitations
- Decoherence is omitted, so all reported fidelities are optimistic upper bounds for hardware.
- The constrained optimal-control benchmark is intentionally lightweight. A larger control basis, longer solve, or larger optimization state set could improve its performance further.
- The two-tone calibration is performed for the vacuum branch-matching objective; longer-duration coherent-state performance can still degrade even when the vacuum branches are nearly equalized.
- The study now includes a focused multiplex follow-up with full-duration multicarrier fits, segmented branch-resonant fits, and a low-parameter jointly optimized shaped-two-tone family. All three structured families remain substantially worse than the short `20 ns` two-tone pulse and the bounded sampled waveform on the broad state metric, so any further multiplex improvement likely requires a materially richer hardware-aware parameterization.
- The echo protocol uses vacuum-calibrated qubit `pi` pulses. It is therefore a test of practical echoing, not of an idealized manifold-uniform qubit inversion.
- The study focuses on one representative dispersive operating point. The `chi` sweep clarifies scaling, but a full multi-parameter device-design map remains future work.

## Validation
- [x] Sanity checks
  - The minimal analytic picture is reproduced numerically: short pulses behave nearly unconditionally, while branch mismatch grows rapidly with `|chi| T`.
  - The hierarchy sweep at `alpha = 1`, `T = 80 ns` shows that adding `chi'` and self-Kerr only perturbs the naive result slightly (`delta_alpha = 0.682` minimal, `0.686` with `chi'`, `0.686` full), confirming that `chi` is the dominant mechanism in this regime.
  - The best two-tone case (`T = 20 ns`) equalizes the vacuum-branch displacements to `delta_alpha = 8.15e-4`, matching the intended compensation picture.
- [x] Convergence
  - A direct unconditional-displacement spot check for the representative naive square pulse (`alpha = 1`, `T = 80 ns`) gave:
    - `n_cav = 12 -> 15 -> 18`: `delta_alpha = 0.685867, 0.685870, 0.685871`
    - `n_tr = 3 -> 4`: `delta_alpha = 0.685870 -> 0.685866`
    - `dt = 0.25 -> 0.5 -> 1.0 ns`: `delta_alpha = 0.681879, 0.685870, 0.693775`
  - The study also reuses the same propagation stack and truncation settings already validated earlier in this study family, so the numerical conclusions are stable at the chosen defaults.
- [x] Literature comparison
  - The naive single-tone crossover is consistent with the dispersive criterion `|chi| T << 1` and with the number-splitting intuition from the dispersive cQED literature.
  - The extracted guideline `1 / |chi| = 56.0 ns` matches the observed transition from near-unconditional behavior at `5-10 ns` to strong conditionality by `80-160 ns`.

## Status
COMPLETE
