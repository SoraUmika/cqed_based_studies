# Improvement Log: Sequential Active Cooling via the `|g,0_r,n_s> <-> |f,0_r,n_s-1>` Ladder

> Written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- **[P1 | HIGH] No microscopic pump-generation model inside the native study workflow**: `cqed_sim` exposes the storage and readout sidebands as effective control operators. That is good enough for spectroscopy, pulse screening, and cooling-cycle simulation, but it does not derive the actual coupling rate, Stark shifts, or parasitic channels from the circuit nonlinearity and drive routing. Any future claim about hardware power calibration or microscopic selection rules needs a stronger model.
- **[P1 | HIGH] Effective amplitudes cannot yet be translated into instrument settings**: The recommended `6-14 MHz` Step A amplitudes and `12 MHz` Step B amplitudes are internal effective rates, not AWG voltages or room-temperature source powers. Experimental transfer will require an explicit calibration bridge.

## Recommended Improvements (P2)
- **[P2 | MEDIUM] Add a microscopic or semi-microscopic sideband calibration layer**: The next valuable extension is a pump-aware Hamiltonian or calibration wrapper that predicts how the effective sideband rate depends on the experimental drive amplitude, detuning, and hardware chain.
- **[P2 | MEDIUM] Extend the strong-drive analysis beyond the target doublet**: The Floquet calculation reproduces the expected quasienergy splitting, but every recommended pulse also triggers a truncation-boundary warning. A higher-confidence strong-drive study should increase the accessible Hilbert space or add targeted reduced models to separate real leakage from truncation artifacts.
- **[P2 | MEDIUM] Add explicit readout-reset optimization**: The passive `4/kappa_r` wait leaves residual readout occupation around `1.7e-2` after one cycle. Repeated-cycle hardware experiments may need either a longer wait or an active readout-reset primitive.
- **[P2 | MEDIUM] Validate the qutrit-readout versus dump-witness strategy on hardware**: The spectroscopy follow-up indicates that Step A should be found through `P_f` first, with a Step-B-assisted readout witness as the best fallback when `f` discrimination is weak. A hardware study should compare those two observables directly and quantify which one gives the best line contrast and calibration repeatability.

## Nice-to-Haves (P3)
- **[P3 | LOW] Extend the ladder beyond `n_s = 4`**: The present study establishes the low-photon regime. Scaling to larger Fock support would test whether the `5.680842 MHz` line spacing becomes prohibitive.
- **[P3 | LOW] Add calibration-aware robust-control optimization**: The present sweep already includes local detuning and amplitude sensitivity, but a future pass could optimize directly against those uncertainties instead of evaluating them afterward.
- **[P3 | LOW] Add finite storage thermal population or explicit heating channels**: The current model starts from idealized initial states and does not include an active heating bath.

## Open Questions
- The old locally stored nominal `n_s = 1` storage-sideband frequency is lower than the static model prediction by `581.406 kHz`, while the new effective-model detuning scan finds its optimum at `0.0 MHz`. That mismatch could reflect an older calibration point, extra dressing absent from the effective model, or a frame-convention issue; it still needs an experiment-facing explanation.
- The repeated coherent-state ladder leaves `4.17%` transmon excitation after the first `n_s = 4` cycle. It is still unclear whether that is mostly unavoidable multi-manifold crowding or whether a more selective robust pulse could reduce it.
- Step A is best served by the paper-inspired `bump` family for `n_s = 1, 3, 4`, but `n_s = 2` prefers a short square pulse. It would be useful to understand whether that is a genuine even-manifold resonance effect or simply a grid-resolution artifact of the present scan.

## What Was Tried and Did Not Work
- **Direct transmon `g-f` carrier as the cooling step**: A carrier pulse populates `|f,0_r,n_s>` with probability above `0.99999994` while producing zero population in `|g,0_r,n_s-1>`. It is useful as a control diagnostic only and should not be described as cooling.
- **DRAG-like Gaussian derivative correction as a default improvement**: The `gaussian_drag` family repeatedly triggered an ODE solver stiffness failure (`Excess work done on this call... increase nsteps`) in both Step A and Step B scans. In this effective model it is screened out as numerically unreliable rather than promoted as a recommendation.
- **Phase-modulated bump as a proxy for stronger parametric richness**: The `phase_modulated_bump` family underperformed badly, reaching only about `0.677` transfer in the best Step B cases. The added phase modulation appears to smear the desired two-state swap more than it helps in the present effective model.
- **Assuming strong-drive Stark shifts would automatically move the best frequency away from the static dressed line**: The `+-2 MHz` calibration scans around every recommended pulse found the best point at `0.0 MHz` detuning. The sideband paper motivated the check, but this effective model does not show a resolved Stark-shifted optimum for the recommended operating points.

## Compute & Resource Notes
- Full study rerun (`scripts/run_study.py`): `157.23 s` wall clock on CPU.
- Targeted convergence pass (`scripts/convergence_checks.py`): `13.80 s`.
- Reproducibility notebook execution (`python -m nbconvert --to notebook --execute --inplace reproducibility_notebook.ipynb`): `10.44 s`.
- Floquet analysis was restricted to the recommended pulses rather than the full scan grid to keep runtime moderate.
- No additional package installs were required; the study used the existing system Python 3.12 environment.

## Resolved
- **[Resolved] User-facing basis ordering now matches the study prompt**: All exported tables, figures, and labels now report states in the `|q,n_r,n_s>` convention even though the internal tensor ordering remains `(transmon, storage, readout)`.
- **[Resolved] Paper-guided pulse-family comparison was added explicitly**: The final artifact set now compares square, Gaussian, DRAG-like, bump, phase-modulated bump, cosine-squared, and composite-pulse-inspired envelopes, with clear winners and failure cases.
- **[Resolved] Strong-drive checks were added instead of guessed**: The final study now includes a direct frequency-calibration scan and a Floquet doublet-splitting summary to test whether the paper's Stark-shift and Floquet lessons actually appear in the local effective model.
- **[Resolved] Convergence was updated for the new baseline truncation**: The final `n_tr = 4`, `n_storage = 7`, `n_readout = 3` baseline was revalidated, with only `O(10^-4)` changes under nearby truncation and timestep variations.
- **[Resolved] Experiment-facing spectroscopy guidance now exists**: The study now includes a dedicated follow-up showing that the correct first-pass Step A measurement is `P_f`, that direct `g-f` carrier control is useful for Step B calibration but not for cooling, and that a dump-assisted witness is the right fallback when direct `f` readout is weak.
