# Waveform Optimization for the `|f,n-1> <-> |g,n>` Sideband Interaction in cQED

## Problem Class
OPT, ANA, DES

## Motivation
This study asks which control-waveform family is best for driving the effective
\[
|f,n-1\rangle \leftrightarrow |g,n\rangle
\]
sideband interaction in the local simulator, with separate conclusions for storage-mode and readout-mode sidebands. The goal is to produce experimentally useful guidance, not only a closed-system ranking, so the final workflow now includes explicit device-parameter provenance, open-system replay, and a transmon-decoherence sensitivity extension.
sideband interaction in the local simulator, with separate conclusions for storage-mode and readout-mode sidebands. The goal is to produce experimentally useful guidance, not only a closed-system ranking, so the final workflow now includes explicit device-parameter provenance, open-system replay, a transmon-decoherence sensitivity extension, and a simultaneous two-tone extension that asks whether a continuous storage-sideband plus readout-sideband drive can coherently move a single storage photon into the readout mode quickly enough to matter on this device.

## Goals
1. Determine the dressed storage-sideband and readout-sideband resonance frequencies for the `|f,n-1> <-> |g,n>` ladder using the actual local device parameters.
2. Compare square, Gaussian, cosine, flat-top cosine, flat-top Gaussian, smooth compact-support, and Blackman envelopes with family-specific parameter variation.
3. Define explicit operational metrics for selective control, fast/unselective control, leakage, neighboring-manifold selectivity, and gate-oriented projected-subspace quality.
4. Identify the best selective and fastest/unselective waveform families for the storage and readout sidebands in the closed-system truth model.
5. Re-evaluate the recommendations once mode dissipation and transmon decoherence are included.
6. Export reproducible scripts, figures, machine-readable artifacts, a markdown summary, a LaTeX report/PDF, and a reproducibility notebook.
7. Evaluate a simultaneous two-tone storage-to-readout transfer scheme, identify the physical mechanism that enables it, and determine the achievable transfer times in the local cQED parameter regime.

## Methods
- `cqed_sim.core.DispersiveReadoutTransmonStorageModel` for the full multilevel transmon-storage-readout drift Hamiltonian.
- `cqed_sim.core.drive_targets.SidebandDriveSpec` for effective storage-sideband and readout-sideband activation.
- `cqed_sim.core.frequencies.sideband_transition_frequency` and `basis_energy` for exact model-consistent spectroscopy.
- `cqed_sim.sequence.scheduler.SequenceCompiler` and `cqed_sim.sim.runner.simulate_sequence` for time-domain closed- and open-system replay.
- `cqed_sim.sim.noise.NoiseSpec` for storage, readout, and transmon decoherence channels.
- Study-local helpers for waveform construction, amplitude calibration, projected two-state diagnostics, neighboring-manifold scans, robustness maps, transmon-noise scenario management, and artifact/report generation.
- A study-local simultaneous-drive extension using overlapping storage-sideband and readout-sideband channels, with comparison against the minimal reduced ladder Hamiltonian for the single-photon case.

## Analytic Preliminary
### Basis, Units, and Tensor Ordering
The simulator tensor order is `(transmon, storage, readout)`. Human-facing state labels are reordered into the requested notation `|q,n_r,n_s>`. Internal frequencies are `rad/s`, times are `s`, and all reported values are converted into `Hz`, `MHz`, `GHz`, and `ns`.

### Drift Hamiltonian
The truth-model drift Hamiltonian is
\[
\frac{H_0}{\hbar} =
\omega_q \hat n_q
+ \frac{\alpha}{2}\hat n_q(\hat n_q-1)
+ \omega_s \hat n_s
+ \omega_r \hat n_r
+ \chi_s \hat n_q \hat n_s
+ \chi_r \hat n_q \hat n_r
+ \chi_{sr} \hat n_s \hat n_r
+ \frac{K_s}{2}\hat n_s(\hat n_s-1)
+ \frac{K_r}{2}\hat n_r(\hat n_r-1).
\]

The key first-principles expectation is that the red-sideband matrix element grows like `sqrt(n)` while spectral selectivity is limited by the same-mode line spacing `|2 chi - K|`. This predicts:
- faster transfer at larger `n` for fixed amplitude,
- broader spectral spillover for shorter pulses,
- smoother edges reducing neighboring-manifold excitation,
- and a tension between selective control and practical duration.

### Effective Sideband Picture
For a selected mode `m in {storage, readout}`, the simulator applies an effective interaction of the form
\[
H_{\mathrm{sb}}(t)/\hbar
\sim \Omega(t)\left(
|f,n_m-1\rangle\langle g,n_m|
+ |g,n_m\rangle\langle f,n_m-1|
\right)
\]
inside the intended two-state manifold while still evolving the full multilevel Hilbert space. The study therefore uses the full simulator for ranking and only uses the reduced picture for intuition.

### Working Definitions
- `Selective pulse`: `P_target >= 0.99`, `P_leak <= 0.02`, and `P_neighbor^max <= 0.01`.
- `Fast unselective pulse`: `P_target >= 0.985` and `P_leak <= 0.03`.
- `Gate-oriented diagnostic`: the fast-transfer floor plus `P_neighbor^max <= 0.02`, then rank by projected `2 x 2` SWAP fidelity.

### Simultaneous Two-Tone Transfer Extension
For the single-photon case `|g,0,1> -> |g,1,0>`, the simultaneous drive adds a storage red sideband and a readout red sideband at the same time, producing the minimal three-state ladder
\[
|g,0,1\rangle \leftrightarrow |f,0,0\rangle \leftrightarrow |g,1,0\rangle.
\]
Within that reduced basis, the rotating-frame model is
\[
\frac{H_{\mathrm{2tone}}}{\hbar} =
\begin{pmatrix}
0 & g_s & 0 \\
g_s & \Delta & g_r \\
0 & g_r & 0
\end{pmatrix},
\]
where `g_s` and `g_r` are the effective storage-leg and readout-leg couplings and `Delta` is the common one-photon detuning of the intermediate `|f>` state.

This gives two controlled limits.
- `Resonant bright-state transfer`: `Delta = 0`. For equal couplings `g_s = g_r = g`, the first perfect transfer occurs at `T_res = pi / (sqrt(2) g)` and the intermediate state reaches order-unity occupation.
- `Detuned Raman-like transfer`: `|Delta| >> g`. The effective storage-readout exchange rate becomes `g_eff ~ g_s g_r / Delta`, so the transfer slows to `T_Raman ~ pi Delta / (2 g_s g_r)` while the peak `|f>` occupation is reduced to order `(g / Delta)^2`.

The extension also reveals a crucial many-photon caveat. The `n = 1` case is an isolated three-state ladder, but `n > 1` opens a longer chain under a constant simultaneous drive, so the clean near-unit transfer result is specific to the single-photon protocol rather than to arbitrary storage Fock states.

## cqed_sim Gap Analysis
| Functionality | Needed? | Available in cqed_sim? | Plan |
|---|---|---|---|
| Full multilevel dispersive truth model for transmon + storage + readout | Yes | Yes | Use the native model directly |
| Effective storage and readout `g-f` red-sideband controls | Yes | Yes | Use `SidebandDriveSpec` with `lower_level=0`, `upper_level=2` |
| Exact model-consistent sideband spectroscopy | Yes | Yes | Use `sideband_transition_frequency` and direct basis-energy checks |
| Simultaneous storage-sideband and readout-sideband replay on separate channels | Yes | Yes | Use overlapping pulses compiled together through `SequenceCompiler` |
| Family-by-family arbitrary-envelope replay | Yes | Partial | Use study-local envelope builders on top of the native pulse stack |
| Transmon `T1/Tphi` replay | Yes | Yes | Use `NoiseSpec(t1=..., tphi=...)` and document parameter provenance |
| Direct built-in optimizer over all requested analytic families | Helpful | No | Use principled family-specific sweeps seeded by the analytic `pi`-area estimate |
| Microscopic derivation of the pump-induced sideband rate and Stark shift | Yes | No | Keep the effective-control boundary explicit and record it as a remaining limitation |

## Assumptions
- The local sideband-reset example supplies the authoritative Hamiltonian and mode-noise parameters.
- The local tomography workflow supplies a reasonable transmon-decoherence sensitivity anchor because the sideband-reset example itself does not export transmon `T1/T2`.
- The effective sideband-control abstraction is acceptable for waveform ranking within the current environment.
- Low-lying manifolds `n = 1,2,3` are sufficient for the main comparison.
- The working truncations are adequate because the finalist ranking is explicitly convergence-checked.
- The single-photon simultaneous-transfer case is the main experiment-facing target; higher storage Fock states are treated as a diagnostic extension because a constant two-tone drive opens a longer conversion chain.

## Compute & Resource Strategy
The dominant cost is the waveform-family sweep across mode, `n`, duration, amplitude, shape parameter, and open-system follow-up. The implementation used:

- analytic `pi`-area seeding to avoid a broad amplitude search at each duration,
- exact dressed frequencies instead of repeated frequency searches,
- single-process execution because import overhead is large on this Windows machine,
- finalist-only robustness maps,
- and a transmon-reference reranking pass over the shortlist rather than over the full amplitude grid.

Realized runtimes:
- main waveform study plus transmon-reference reranking: `534.8 s`
- simultaneous two-tone extension: `21.1 s`
- validation pass: `3.6 s`
- notebook build + execution after the extension: `~6.8 s`
- report compilation after the extension (three `pdflatex` passes): `~6 s`

## Expected Outcomes
The completed study now supports three layers of recommendation:
- closed-system ranking:
  - best selective family: `gaussian (sigma_fraction = 0.24)` for both modes
  - best fast family: `square` for both modes
- mode-noise-only open-system replay:
  - storage recommendations survive
  - selective readout control fails
- transmon-inclusive replay using the matched local reference:
  - no family remains threshold-valid in the selective regime
  - only the storage fast-square pulse remains threshold-valid overall
- simultaneous two-tone storage-to-readout transfer:
  - the single-photon resonant bright-state protocol reaches `0.999998` closed-system target probability in `29.5 ns` at equal `12 MHz` leg couplings
  - the best detuned Raman-like single-photon case lowers the peak `|f>` population to `0.121` but stretches the closed-system transfer to `348 ns`
  - under readout decay and the matched transmon reference, the resonant single-photon case still peaks near `0.953`, whereas the detuned case drops to about `0.556`
  - for `n > 1`, a constant simultaneous drive leaks into a longer ladder and no longer realizes the clean single-transfer picture

## Known Limitations
- The sideband is still an effective control operator rather than a microscopic pump Hamiltonian.
- The sideband-reset example still does not export device-matched transmon `T1/T2`, so the transmon-noise extension is a sensitivity study rather than a full device-verification pass.
- The gate-oriented extension adds a ranking diagnostic, but not a direct unitary-level optimization.
- The winner pulses remain state-transfer primitives, not phase-clean coherent SWAP gates.
- The simultaneous two-tone extension only tests constant square drives; it does not yet include counter-intuitive STIRAP-like timing, shaped overlapping pulses, or direct open-system optimal control.

## Validation
- [x] Sanity checks
- [x] Convergence
- [x] Literature comparison (if applicable)

Validation notes:
- exact rotating-frame sideband frequencies agree with the analytic dispersive formulas to numerical precision for both modes and all `n = 1,2,3`
- finalist baselines use `0.25 ns`; repeating them at `0.5 ns` changes target transfer by only `1.8e-4` to `1.4e-3`
- increasing the truncation from `(n_tr, n_s, n_r) = (4,5,5)` to `(5,6,6)` changes finalist target transfer by at most `5.7e-5`
- the closed-system control ranking matches the usual sideband expectation that smooth pulses win in the selective regime and square pulses win in the fastest regime
- the extended open-system pass shows that practical threshold-valid control is much more restrictive once the matched transmon reference is included
- the simultaneous two-tone single-photon ladder matches the reduced three-state model to within `4.1e-7` in peak target probability and exactly in peak time for the resonant selected case
- the two-tone validation pass shows `<= 2.5e-3` peak-target variation when the single-photon representative cases are rerun at `1.0 ns` instead of `0.25 ns`, and `< 2.3e-6` peak-target variation under the larger `(5,6,6)` truncation for the selected single-photon cases

## Status
COMPLETE
