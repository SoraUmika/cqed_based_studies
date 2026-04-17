# Improvement Log: Can an Ideal Multitone SQR Realize Arbitrary Fock-Conditional Qubit Rotations?

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- **[P1 | MEDIUM] Test independently parameterized echoed halves**: The present echoed extension uses the exact repeated-half schedule requested by the user. That isolates refocusing, but it does not test whether a composite sequence with independently optimized first and second SQR halves can recover the structured-family fidelity that the repeated-half echo destroys.
- **[P1 | MEDIUM] Distinguish finite-pulse versus ideal-refocusing effects**: The inserted `X_pi` pulses are finite Gaussian pulses resonant with the `n = 0` transition. A follow-up should compare against an idealized instantaneous refocusing pulse or a more broadband unconditional `X_pi` construction to separate echoed-sequence limitations from refocusing-pulse imperfections.

## Recommended Improvements (P2)
- **[P2 | MEDIUM] Replay the best echoed cases on a qutrit model**: The echoed comparison remains in the strict `n_tr = 2` workflow. Replaying the strongest random-case improvement and the strongest structured-family failure at `n_tr = 3` would test whether the echoed sequence introduces additional transmon leakage or whether the dominant defect remains coherent block mismatch.
- **[P2 | MEDIUM] Add multistart or a second optimizer for echoed cases**: The explicit echoed comparison already shows strong target-class asymmetry, but the best random-case improvements should still be checked against optimizer stagnation with a modest multistart budget.
- **[P2 | LOW] Add one or two intermediate durations near `|chi| T / 2 pi = 2`**: The current echoed extension uses `1, 3, 5`. A moderate-duration interpolation would test whether the weak random-target improvement of `echoed_fixed_active` is smooth or whether there is a narrow timing sweet spot.

## Nice-to-Haves (P3)
- **[P3 | LOW] Add Bloch-sphere visual summaries for the best random improvement and the strongest structured failure**: The saved residual-Z/transverse block decomposition already shows the main effect, but a Bloch-rotation visual could communicate the axis mismatch more directly.
- **[P3 | LOW] Add robustness sweeps in pi-pulse duration and phase**: The current extension fixes `tau_pi = 40 ns` and `phase = 0`. A small sweep would show whether the modest random-target benefits are robust or fine-tuned.

## Open Questions
- Why does the repeated-half echoed sequence improve a substantial subset of random family D cases while catastrophically degrading family C, even though family C is the strongest structured baseline family?
- For the random ensemble, is the observed benefit of `echoed_fixed_active` mainly a weak reduction of matched-case residual-Z error, or is it really selecting a different subset of targets that happen to align better with the repeated-half echoed ansatz?
- Would an independently parameterized second half recover the structured-family frontier without losing the random-case benefits, or is the repeated-half refocusing idea itself too restrictive for arbitrary conditional block control?

## What Was Tried and Did Not Work
- **Original single-pulse-only report**: The first report tested only the single Gaussian multitone SQR ansatz. It did not include an echoed or composite refocusing sequence, so its negative conclusion applied narrowly to the single-pulse ansatz.
- **Explicit repeated-half echoed extension**: The exact time-ordered schedule `half SQR -> X_pi -> half SQR -> X_pi` was implemented and validated under two fairness conventions (`echoed_fixed_total`, `echoed_fixed_active`). For the structured family C subset, both echoed branches failed badly: the single-pulse median fidelity was `0.7768`, while the medians dropped to `0.2109` (`echoed_fixed_total`) and `0.1718` (`echoed_fixed_active`). Improved-fidelity count for family C was `0/18` in both echoed branches.
- **Residual-Z-only explanation for family C**: The structured-family failure is not just a block-phase problem. On family C, the echoed branches increased the matched-case median transverse coherent error by `+0.2475 rad` (`echoed_fixed_total`) and `+0.3012 rad` (`echoed_fixed_active`) while also worsening fidelity by about `-0.56` to `-0.59` in median.
- **Uniform echoed success hypothesis**: The echoed construction is not uniformly useless either. On the random family D subset, `echoed_fixed_active` improved fidelity in `38/72` matched cases and reduced residual-Z error in `38/72`. The strongest random-case improvement came from `echoed_fixed_total` on `chi_only_na3_chiT5p0_familyD_seed317160`, where fidelity increased by `+0.428582` to `0.606132`.
- **Fixed-total echoed refocusing as a general solution**: `echoed_fixed_total` did not show a general random-target benefit. On family D its matched-case median fidelity change was slightly negative (`-0.002212`) and its matched-case median residual-Z change was slightly positive (`+0.026359 rad`).

## Compute & Resource Notes
- The original 192-case single-pulse sweep remains the heavier baseline study and has already been validated in `data/validation_summary.json`.
- The echoed extension reused the baseline single-pulse subset and added `180` echoed optimizations across `90` matched operating points. Terminal timings for individual echoed cases ranged from about `1.5 s` to about `26.9 s`, with the longest cases appearing in the `N_active = 4`, long-duration random subset.
- Echo-extension validation replays for the best echoed case (`chi_only_na3_chiT5p0_familyD_seed317160`, `echoed_fixed_total`) changed average gate fidelity by only `-8.05e-4` under a finer `2 ns` step and by only `+6.26e-7` under two extra cavity levels.

## Resolved
- **Explicit echoed sequence implemented**: The study now explicitly tests `half SQR -> X_pi -> half SQR -> X_pi` and records the exact sequence order in `data/echo_comparison_validation.json`.
- **Fair duration conventions added**: The extension now compares `single_pulse`, `echoed_fixed_total`, and `echoed_fixed_active` on matched operating points.
- **Residual-Z versus transverse-error separation**: The comparison dataset and figures now separate blockwise residual-Z error from transverse coherent error so the report can distinguish refocusable phase error from deeper controllability limits.
- **Best-case artifacts saved**: Dedicated highlight artifacts now exist for the best single-pulse and best echoed cases under `artifacts/echo_comparison/highlights/`.# Improvement Log: Can an Ideal Multitone SQR Realize Arbitrary Fock-Conditional Qubit Rotations?

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- **[P1 | MEDIUM] Verify optimizer robustness beyond local minima**: The standard multitone fit will be judged on a strict block-diagonal target, so local minima can masquerade as controllability limits. If final random-ensemble tails are broad, add multistart or a second optimizer before drawing strong no-go conclusions.

## Recommended Improvements (P2)
- **[P2 | MEDIUM] Add qutrit replay follow-up**: The main strict-logical workflow uses `n_tr = 2`. Replay the best and worst closed-system waveforms at `n_tr = 3` to distinguish coherent block mismatch from transmon leakage.
- **[P2 | LOW] Expand duration sweep in the moderate regime**: If the fidelity-duration frontier has a sharp knee, add extra points between `chi T / 2 pi = 2` and `3`.
- **[P2 | MEDIUM] Compare Gaussian and segmented multitone ansatze**: This study intentionally fixes the standard Gaussian multitone SQR family. If the Gaussian family fails, the next question is whether the failure is fundamental or only ansatz-limited.

## Nice-to-Haves (P3)
- **[P3 | LOW] Add robustness sweeps in `chi` and calibration offsets**: Once the ideal closed-system picture is established, perturb `chi`, `chi'`, duration, and per-tone frequencies to estimate calibration sensitivity.
- **[P3 | LOW] Add Bloch-sphere visual summaries per block**: A compact visual summary of target versus realized Bloch rotations would help presentation but is not required for the core conclusion.

## Open Questions
- Does the `d_omega` correction parameter genuinely unlock arbitrary blockwise Z control, or does the strict fidelity remain limited because the Gaussian ansatz couples phase and amplitude too strongly across manifolds?
- When random targets fail, is the dominant defect spectral crowding, insufficient control parameters, or a hidden block-phase structure that the raw waveform cannot realize directly?
- Why does `chi'` slightly improve the sampled random-target medians at `N_active = 3, 4, 5` in this sweep even though it reduces structured-family performance at larger active windows?

## What Was Tried and Did Not Work
- **Initial strict-SU(2) pilot, `N_active = 2`, family B, `chi`-only, `|chi| T / 2 pi = 3`**: The optimized Gaussian multitone waveform reached average gate fidelity only `~0.45` with process fidelity `~0.31`, while the single-tone crosstalk proxy stayed at `~5e-8`. This indicates that low fidelity can arise even when off-manifold excitation is negligible; the defect is not purely spectator crosstalk and is consistent with insufficient controllability of the standard Gaussian ansatz for the requested blockwise SU(2) map.
- **`N_active = 1` sanity pilot**: Different arbitrary SU(2) targets are not equally easy even without spectral crowding. Representative average gate fidelities at `|chi| T / 2 pi = 5` were `~0.47` (family A), `~0.54` (family B), and `~0.95` (family C). This suggests that the standard waveform can realize some isolated arbitrary blocks well, but not arbitrary SU(2) targets uniformly.
- **Full 192-case sweep**: No case exceeded average gate fidelity `0.8737`, and the random-family median over the entire study was only `0.1917` with a best random case of `0.5022`. All 192 cases remained in the study's failure tier, so the strict arbitrary-control claim is not supported anywhere on the explored grid.
- **Spectator-excitation hypothesis**: Median crosstalk maxima across families A/B/C/D were only about `5e-9`, `3e-8`, `7e-9`, and `1e-8`, respectively. The best random case still had crosstalk off-diagonal max `~2.5e-10` while its average gate fidelity was only `~0.50`, so suppressing spectator response is not sufficient to realize the target block operator.
- **Monotonic-`chi'`-penalty hypothesis**: The random-target median at `N_active = 4`, `|chi| T / 2 pi = 3` shifted from `0.1326` (`chi` only) to `0.1510` (`chi + chi'`). The higher-order shift therefore did not act as a simple monotonic penalty in the sampled random ensemble.

## Compute & Resource Notes
- Early pilot optimizations on the local system Python 3.12.10 environment took about `7 s` per `N_active = 2` case with the current two-stage Powell + L-BFGS-B configuration.
- The completed 192-case sweep accumulated about `2959 s` of optimizer wall time in the saved case summaries (`~49.3 min` total, `15.4 s` mean per case, `72.8 s` max case runtime).