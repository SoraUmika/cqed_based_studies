# Executive Summary

The previous report tested only a single Gaussian multitone SQR pulse. It did **not** test any echoed, composite, or refocused SQR sequence, so its negative conclusion applied narrowly to that single-pulse ansatz.

This extension adds the explicit echoed schedule

\[
\text{half SQR} \rightarrow X_\pi \rightarrow \text{half SQR} \rightarrow X_\pi
\]

under two fair timing conventions:

- `echoed_fixed_total`: same total gate duration as the single pulse.
- `echoed_fixed_active`: same total active SQR time as the single pulse, with the two finite `X_pi` pulses added on top.

What changed technically:

- Added `scripts/run_echo_comparison.py` for the explicit echoed comparison.
- Added `scripts/validate_echo_comparison.py` for artifact and convergence validation of the echoed extension.
- Saved dedicated highlight artifacts for the best single-pulse and best echoed cases under `artifacts/echo_comparison/highlights/`.
- Generated new comparison figures: branch means, duration scans, echo-minus-single tradeoff scatter, best-waveform plots, and block-error breakdowns.

Main verdict:

- The exact repeated-half echoed sequence does **not** rescue the structured family C cases. It is dramatically worse than the baseline single pulse there.
- The echoed sequence is **not uniformly harmful** either. On the random family D subset, `echoed_fixed_active` improves fidelity in `38/72` matched cases and reduces residual-Z error in `38/72` matched cases.
- The strongest echoed improvement reaches fidelity `0.606132` on `chi_only_na3_chiT5p0_familyD_seed317160` (`echoed_fixed_total`), improving that specific random case by `+0.428582` over its single-pulse baseline.
- Even so, the best echoed result remains well below the best single-pulse structured result (`0.873625`), and the echoed branch means are lower than the single-pulse mean.

Interpretation:

- The old report should be read as a failure of the **single-pulse Gaussian multitone SQR ansatz**, not as a general impossibility theorem.
- The newly tested explicit echoed construction provides only **partial, target-dependent mitigation**. It sometimes helps random targets and sometimes reduces residual-Z error, but it does not solve the broader arbitrary-control problem and it badly degrades the strongest structured family.
- The best current practical verdict is therefore: the single-pulse Gaussian ansatz is inadequate for strict arbitrary block control, and the exact repeated-half echoed extension only partially mitigates that limitation rather than overturning it.