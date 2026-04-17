# Improvement Log: Holographic Cluster-State Control in cQED

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- **[P1 | HIGH] Dense larger-truncation GRAPE sweep is still missing**: The final recommendation is based on direct `N_cav = 12` rescue runs at `300 ns` and `400 ns`. This was enough to overturn the misleading `N_cav = 8` replay picture and identify a credible `400 ns` pulse, but it does not yet prove that `400 ns` is globally optimal over a dense `N_cav >= 12` duration grid. Re-run the direct larger-truncation sweep at intermediate and slightly longer durations.
- **[P1 | MEDIUM] Noisy optimization has not yet been done**: Open-system performance was evaluated by replaying the best closed-system pulses with `T1`, `Tphi`, and cavity `kappa`, not by optimizing in the noisy model. The current best `400 ns` pulse still wins after noisy replay, but a true Lindblad optimization could shift the preferred duration or waveform.
- **[P1 | HIGH] `SNAP` and `FreeEvolveCondPhase` still lack waveform-export support**: This is now the main simulator limitation blocking a fair replay-level comparison between structured decompositions and GRAPE. Until those primitives can be exported and replayed through the same pulse path, decomposition families remain validation-asymmetric.

## Recommended Improvements (P2)
- **[P2 | MEDIUM] Add robustness objectives to larger-truncation GRAPE**: The current best pulse is optimized for nominal Hamiltonian parameters only. Penalizing sensitivity to small `chi`, drive-scale, and detuning errors would make the experiment-facing recommendation more credible.
- **[P2 | MEDIUM] Expand the Wigner probe set beyond logical basis inputs**: The logical basis already separates the successful GRAPE candidate from the unsuccessful structured routes, but superposition probes would provide a stricter operator-level phase-space test.
- **[P2 | LOW] Add automatic truncation-adaptive replay into the optimization workflow**: The `N_cav = 8` sweep looked excellent before replay at `N_cav = 10` and `12` exposed the truncation failure. The next pass should fail fast when replay fidelity collapses at larger truncation.

## Nice-to-Haves (P3)
- **[P3 | LOW] Add a long-chain observable follow-up from the corrected channel picture**: The transfer-channel inconsistency is fixed, but the present report still focuses on the per-site control step. A later study could connect the corrected `xi = 0` channel statement to explicit long-chain cluster-state observables and sampling consequences.
- **[P3 | LOW] Include hardware-calibrated AWG / transfer-function effects**: Current replay uses the simulator control stack directly rather than a hardware-calibrated pulse-distortion model.

## Open Questions
- Can a direct `N_cav = 12` GRAPE pulse shorter than `400 ns` match the current `0.956` replay fidelity and `0.901` open-system process fidelity once optimized with comparable care?
- Does robustness-aware or noise-aware GRAPE materially change the preferred duration?
- Is there a structured qubit-conditional family that becomes competitive only after the missing pulse-export paths are added?

## What Was Tried and Did Not Work
- **Plain `D-R-SNAP` remained non-expressive**: The best reduced-subspace fidelity was only `0.5000`, and the bounded version stayed at `0.4422` despite negligible embedded leakage (`6e-4`). This confirms that cavity-only `SNAP` in the installed framework cannot supply the needed entangling structure.
- **Exact reduced-subspace `D-SQR-CPSQR` agreement did not survive physical validation**: The unbounded candidate reached ideal fidelity `1.0000` at `N_cav = 2`, but embedded fidelity fell to `0.5977` at `N_cav = 12` with `48.8%` average leakage. This invalidates the earlier reduced-model overclaim.
- **Bounded SQR-like routes stayed physically weak after replay**: The bounded `D-SQR-CPSQR` candidate had modest embedded fidelity (`0.5288`) and low replay leakage (`3.23%` average), but pulse replay still reached only `0.3495`. The bounded `D-R-SQR-CPSQR` variant replayed even worse (`0.2043`).
- **The coarse `N_cav = 8` GRAPE frontier was not truncation-converged**: The best `400 ns` pulse replayed at `0.9793` in the optimization truncation, but collapsed to `0.1559` at `N_cav = 10` and `0.1228` at `N_cav = 12`, with leakage above `80%`. Do not reuse the coarse-sweep recommendation without larger-truncation replay.
- **Not all larger-truncation seeds were good just because nominal fidelity was high**: In the direct `N_cav = 12` rescue, the `400 ns` seed `42` reached nominal fidelity `0.9883` yet replayed at only `0.4059` with `56.6%` leakage. The replay check is essential, not optional.

## Compute & Resource Notes
- `run_followup_study.py --stages transfer decomp grape wigner summary` generated the core data products: corrected channel diagnostics, structured-candidate metrics, the coarse multiseed GRAPE sweep, Wigner grids, and the final comparison summary.
- Direct `N_cav = 12` rescue optimization used three seeds (`17`, `42`, `73`) per duration with 400 iterations each. The `300 ns` seeds took about `19-21 s` per seed (`~61 s` total), while the `400 ns` seeds took about `28-29 s` per seed (`~85 s` total) on this workstation.
- Wigner-grid generation and figure rendering were light compared with control optimization; the dominant cost in this study was repeated GRAPE solve + replay at larger truncation.

## Resolved
- **Transfer-channel inconsistency removed**: The follow-up now uses the installed `HolographicChannel.from_unitary(...)` convention consistently. The transfer spectrum is `{1, 0, 0, 0}` and the ordinary channel correlation length is `xi = 0`.
- **Independent waveform replay is now part of the final ranking**: The best candidates are no longer ranked by optimizer-reported or reduced-subspace fidelity alone.
- **Wigner validation was added for the final shortlist**: The best decomposition, best SQR-like route, and best GRAPE candidate were all checked on logical basis outputs in phase space.
- **Higher-truncation GRAPE rescue identified a credible final candidate**: Direct re-optimization at `N_cav = 12` recovered a strong `400 ns` pulse with replay fidelity `0.9561`, open-system process fidelity `0.9009`, and Wigner cavity-state fidelities from `0.8829` to `0.9705`.
