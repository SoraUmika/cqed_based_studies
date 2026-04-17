# Improvement Log: Multitone SQR No-Go, Decoupled-Block Limit, and Echo Alternative in Dispersive cQED

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- [P1][MEDIUM] **The formal no-go is controlled, not all-regime**: The analytical proof uses the dispersive block-resolved model and Magnus expansion. If a future study wants a theorem beyond that regime, it needs a different argument rather than a stronger rhetorical version of the present one.
- [P1][MEDIUM] **Finite echoed rescue is still not ruled out universally**: The extension pass broadened the test set to include a symmetry-aligned manifold-aware multitone refocusing pulse and found only partial improvement, not rescue. A future study that wants a general impossibility claim must test a broader refocusing family rather than cite only the Gaussian and current manifold-aware square-like constructions.

## Recommended Improvements (P2)
- [P2][LOW] **Upstream the exact reduced blockwise replay helper**: The equality between the full strict model and the reduced blockwise replay was one of the most important sanity checks in this study. That helper should live in `cqed_sim` rather than being reimplemented locally.
- [P2][LOW] **Upstream blockwise residual-generator diagnostics**: Per-block `X`, `Y`, `Z`, and best-fit block-gauge diagnostics would make future SQR studies much easier to audit.
- [P2][MEDIUM] **Broaden the echoed refocusing family beyond the tested symmetry-aligned square design**: The extension pass found that a manifold-aware multitone refocusing pulse can help on some aligned-`x` checkpoints, but it still falls well short of ideal SQR. Segmented or jointly optimized refocusing pulses are the next meaningful test.

## Nice-to-Haves (P3)
- [P3][LOW] **Extend the target family beyond square shared-line bursts**: The strict square-envelope ansatz answered the present prompt cleanly, but a future study could chart the minimum extra control structure needed to escape the obstruction.
- [P3][LOW] **Compare against explicit per-tone detuning compensation on the same metrics**: The local audit established that earlier studies used `d_omega`, but this study did not rescore those richer controls side-by-side on the exact same blockwise diagnostics.

## Open Questions
- What is the smallest additional control resource beyond amplitude and azimuth that genuinely removes the obstruction: per-tone detuning, explicit `Z` compensation, or segmented control?
- Why does the exact shared-line `off-plus` checkpoint outperform the analytically tuned equal-angle root in the full model even though both remain far from ideal?
- For aligned-`x` targets, how much of the remaining error after manifold-aware refocusing is due to higher-order noncommutation and how much is due to the restricted square-like echo ansatz?
- Can a broader manifold-aware or segmented refocusing family turn the present `partial rescue` into a near-ideal gate, or is the remaining error structural?

## What Was Tried and Did Not Work
- **Strict amplitude/azimuth optimization as a falsification attempt**: Even after multi-start Powell optimization over amplitudes and azimuths only, the best strict shared-line case reached restricted average gate fidelity `0.8058`, the mean was `0.6094`, and the worst case was `0.3011`. The optimizer could not recover the ideal blockwise unitary.
- **Best-fit block-gauge rescue**: Allowing a posteriori blockwise gauge fitting improved the mean strict process fidelity by only about `0.0023`. Hidden removable block phases were not the main source of failure.
- **Ideal instantaneous echo as a rescue mechanism**: The ideal echo reduced the matched-set mean maximum residual-`Z` error from `0.0786 rad` to `0.0135 rad`, but the mean fidelity collapsed from `0.7133` to `0.2018`. It canceled part of the wrong error without restoring the gate.
- **Finite Gaussian echo**: Two `40 ns` Gaussian `pi` pulses made the result worse in every matched case. The mean maximum residual-`Z` error increased to `0.4146 rad`, the mean transverse error increased to `1.9640 rad`, and no finite-echo case outperformed the corresponding plain strict pulse.
- **Second-order tuned cancellation mapping plus exact checkpointing**: The equal-angle aligned-`x` tuned root at `|chi| T / (2 pi) = 0.7151483266` did not produce a hidden exact gate. The tuned shared-line checkpoint still reached only fidelity `0.4270` with maximum residual `Z` `1.3824 rad`, while `off-minus` and `off-plus` stayed at `0.3466` and `0.5176` fidelity respectively.
- **Including `chi'` in the tuned checkpoint**: Repeating the tuned case with `chi + chi'` in the aligned-block spacing changed the result negligibly (`0.4270` fidelity, `1.3822 rad` maximum residual `Z`). The accidental second-order set is not stabilized by this modest nonlinear correction.
- **Finite manifold-aware multitone echo**: A symmetry-aligned manifold-aware refocusing pulse improved both fidelity and residual `Z` relative to the plain pulse on the tuned and `off-plus` checkpoints, but the best cases still stalled at `0.4799` fidelity with `0.7172 rad` residual `Z` (tuned) and `0.5374` fidelity with `0.6908 rad` residual `Z` (`off-plus`). It is a partial help, not a rescue.

## Compute & Resource Notes
- Full production sweep (`run_study.py --n-starts 1 --maxiter 5`): about `529.4 s` wall clock on CPU.
- Extension sweep (`run_extension_study.py`): about `47.6 s` wall clock on CPU.
- Representative convergence and sanity reruns were inexpensive compared with the main sweep and did not require GPU or multiprocessing support.
- No additional packages were installed for this study.

## Resolved
- **Strict-model audit completed**: Earlier local "positive" ideal-SQR studies were traced to extra resources such as `d_omega` or broader waveform families, so they are no longer ambiguous evidence against the strict no-detuning claim.
- **Decoupled-block helper validated**: The stronger approximation was checked against analytic one-tone limits and matched the ideal target with unit fidelity in every tested case.
- **Exact reduced blockwise replay validated**: The reduced replay of the compiled shared-line waveform matched the full strict propagation to machine precision, confirming that the failure is already present in the block-resolved dynamics.
- **Accidental tuned-set map generated and checked**: The extension pass built the explicit two-block tuned-set map and selected the first aligned-`x` equal-angle root at `|chi| T / (2 pi) = 0.7151483266`. Exact shared-line checkpoints showed that this accidental second-order condition does not overturn the no-go.
- **Finite-echo family broadened beyond Gaussian refocusing**: The extension pass added a manifold-aware multitone finite-refocusing test. It improved the tuned and `off-plus` aligned-`x` checkpoints over the plain pulse on both fidelity and residual `Z`, but still failed to approach ideal SQR.
