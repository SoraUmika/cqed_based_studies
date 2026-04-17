# Multitone SQR No-Go, Decoupled-Block Limit, and Echo Alternative in Dispersive cQED

## Problem Class
ANA | REP | DES

## Motivation
This study tests a narrow but important claim about shared-line multitone selective qubit rotations in dispersive cQED. The question is not whether richer waveform engineering can realize high-fidelity conditional gates in general, but whether an ideal simultaneous multitone SQR can be realized exactly when each addressed Fock block is driven by a tone placed at its own dispersive qubit transition, no artificial per-tone detuning is allowed, and the only free corrections are per-tone amplitude and azimuth.

Several earlier studies in this repository are related but do not answer that exact question directly. The closest prior ideal-SQR and arbitrary-block studies optimized `d_omega` frequency corrections or explored broader composite families, so they cannot by themselves establish either the strict no-detuning no-go statement or its limits. This study therefore re-derives the result from first principles, then tries to falsify it numerically using the same `cqed_sim` stack.

All Hilbert-space objects follow the `cqed_sim` qubit-first convention: qubit tensor cavity, with logical basis ordered as `(|g,0>, |e,0>, |g,1>, |e,1>, ...)`.

## Goals
1. Prove a formal no-go statement for exact ideal simultaneous shared-line multitone SQR with amplitude and azimuth corrections only and no artificial per-tone detuning.
2. Make the obstruction explicit through a two-block proof, then generalize to multiple addressed Fock blocks.
3. Numerically attempt to falsify the no-go claim with exact shared-line `cqed_sim` simulations and extensive amplitude/phase optimization.
4. Show analytically and numerically that a stronger decoupled-block approximation does allow ideal blockwise SQR.
5. Study the echoed sequence `half-SQR -> pi -> half-SQR -> pi`, identify what it cancels, what it preserves, and where it fails.
6. Correct the record where earlier repository studies made stronger claims than their assumptions justify.

## Methods
- `cqed_sim.core.DispersiveTransmonCavityModel` and `FrameSpec` for the shared-line dispersive Hamiltonian.
- `cqed_sim.core.frequencies.manifold_transition_frequency` and `carrier_for_transition_frequency` for exact block frequencies and replay carriers.
- `cqed_sim.calibration.conditioned_multitone` for multitone tone construction and shared-line waveform compilation.
- `cqed_sim.calibration.targeted_subspace_multitone` for exact logical-subspace operator extraction and fidelity diagnostics, with optimization restricted to `d_lambda` and `d_alpha` only.
- `cqed_sim.sequence.SequenceCompiler` and `cqed_sim.sim.prepare_simulation` for full pulse-sequence replay of the echoed construction.
- Local study helpers only where `cqed_sim` has no public API for the exact task:
  - a decoupled-block reduced model that drops spectator tones by construction,
  - an exact reduced blockwise replay of the compiled shared-line waveform,
  - blockwise generator and axis-angle diagnostics for residual `X`, `Y`, and `Z` behavior.

## Analytic Preliminary
Start from the minimal block-resolved dispersive Hamiltonian
```text
H_0 = sum_n (omega_n / 2) sigma_z tensor |n><n|,
omega_n = omega_q + chi n + chi' n(n-1) + ...
```
driven through one shared qubit-control line by
```text
H_d(t) = sum_{m in S} (Omega_m f(t) / 2)
  [exp(-i(omega_m t - phi_m)) sigma_+ + exp(+i(omega_m t - phi_m)) sigma_-]
```
with each tone placed exactly at its selected block frequency and no independent `Z`-compensation term.

In the interaction frame of `H_0`, block `n` sees one resonant transverse term from tone `m=n` plus all off-resonant spectator tones with beat notes `Delta_nm = omega_m - omega_n`. The first Magnus term already distorts the intended transverse axis and angle, while the second Magnus term produces blockwise effective `Z` generators. For the square-pulse two-block problem, the second-order coefficients are
```text
zeta_0 = -lambda_1^2 K(Delta,T) + lambda_0 lambda_1 L(Delta,T,delta)
zeta_1 = +lambda_0^2 K(Delta,T) - lambda_0 lambda_1 L(Delta,T,delta)
```
with
```text
K(Delta,T) = [Delta T - sin(Delta T)] / (Delta^2 T) > 0
```
for every `Delta > 0`, and `delta = phi_0 - phi_1`. Setting both `zeta_0` and `zeta_1` to zero forces `lambda_0 = lambda_1` and `L = K`, so after the transverse target has already fixed the available amplitude and azimuth knobs, exact cancellation survives only on a lower-dimensional tuned set. That is the controlled no-go statement used in the report.

The report defines "generically impossible" precisely: the exact-cancellation set has no open interior in the space of nontrivial targets and durations. Accidental special cases may exist, but they are fine tuned rather than robust.

## cqed_sim Gap Analysis
| Functionality | Needed? | Available in cqed_sim? | Plan |
|---|---|---|---|
| Full shared-line dispersive multitone propagation | Yes | Yes | Use `DispersiveTransmonCavityModel` plus compiled multitone waveforms |
| No-detuning amplitude/azimuth-only optimization | Yes | Yes | Reuse targeted-subspace optimizer with `parameters=("d_lambda", "d_alpha")` only |
| Logical-subspace unitary and leakage diagnostics | Yes | Yes | Use targeted-subspace diagnostics from `cqed_sim` |
| Requested `half-SQR -> pi -> half-SQR -> pi` replay | Yes | Partial | Use `SequenceCompiler` and local sequence wrapper |
| Decoupled-block model with spectator tones dropped exactly per block | Yes | No public helper | Implemented locally and validated against analytic one-tone limits |
| Exact reduced blockwise replay of the compiled shared-line waveform | Yes | No public helper | Implemented locally and validated against the full strict model |
| Blockwise effective-generator decomposition into `X`, `Y`, `Z` components | Yes | Partial | Implemented local diagnostics layered on the extracted restricted operator |

## Assumptions
- Dispersive regime with a block-diagonal Fock-resolved qubit Hamiltonian as the analytic baseline.
- Shared qubit control line carrying the full simultaneous multitone waveform.
- Baseline no-go proof neglects inter-Fock leakage and treats each cavity Fock block as a driven qubit sector.
- Main numerical optimization keeps the qubit as a two-level system (`n_tr = 2`) unless convergence checks require a higher truncation.
- Primary strict target: ideal blockwise `XY` rotations with no residual block-dependent `Z` phase.
- Echo analysis distinguishes ideal instantaneous `pi` pulses from practical finite-duration Gaussian `pi` pulses.
- "Generic impossibility" means impossibility away from accidental parameter relations of codimension at least one, not impossibility for every specially tuned case.

## Compute & Resource Strategy
- Expected bottleneck: repeated short unitary propagations for multi-start amplitude/phase optimization across several durations and active-window sizes.
- Planned acceleration: keep the active window modest (`N_active = 2, 3, 4`), reuse compiled waveforms where possible, and parallelize if the wall-clock time grew beyond a few minutes.
- Realized cost: the full production sweep completed in about `529.4 s` on CPU. The workflow remained small enough that additional package installation, GPU backends, or multiprocessing were not required for this study.

## Expected Outcomes
- A defensible analytical no-go statement for the strict simultaneous shared-line no-detuning model.
- Numerical evidence that amplitude/azimuth optimization alone cannot recover the ideal SQR operator under the strict shared-line model.
- A clear counterpoint showing that ideal SQR becomes realizable under the stronger decoupled-block approximation.
- A careful statement of when the echoed construction helps, when it fails, and whether its success is exact, approximate, or restricted to special aligned-axis regimes.

## Known Limitations
- The analytical no-go is a controlled result inside the dispersive block-resolved model, not an all-regime nonperturbative theorem for arbitrary strong driving.
- The numerical study samples a broad but finite family of target rotations rather than every conceivable target.
- The decoupled-block model is deliberately stronger than the physical shared-line problem and must not be conflated with it.
- The finite echo study uses vacuum-calibrated Gaussian `pi` pulses; more manifold-aware refocusing pulses were outside scope.

## Validation
- [x] Sanity checks
  - The exact reduced blockwise replay matches the full strict shared-line result with mean and minimum comparison fidelity `1.0`.
  - The decoupled-block construction reproduces the ideal target with fidelity `1.0` in every tested case.
- [x] Convergence
  - Representative strict case (`chi_plus_chiprime`, structured `XY`, `N_active = 3`, `|chi|T/2pi = 3`) stayed at `0.696774` restricted average gate fidelity under a larger optimization budget.
  - Changing `dt` from `2 ns` to `1 ns` shifted the same case only to `0.697201`, and increasing the transmon truncation from `2` to `3` levels shifted it to `0.696958`.
- [x] Literature comparison (if applicable)
  - The analytical structure is consistent with standard dispersive and echo-frame reasoning.
  - The repository-scope audit confirmed that earlier nearby "positive" studies relied on extra resources such as `d_omega` or richer waveform families, so they are not direct counterexamples to the strict claim.

## Suggested Upstreaming
- Add a public `cqed_sim` helper for exact reduced blockwise replay of a compiled shared-line waveform.
- Add public blockwise residual-generator diagnostics (`X`, `Y`, `Z` components and best-fit block gauges).
- Add public support for targeted-subspace optimization over multi-segment echoed sequences.

## Status
COMPLETE
