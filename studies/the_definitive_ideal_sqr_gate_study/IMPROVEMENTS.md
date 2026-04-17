# Improvement Log: The Definitive Ideal SQR Gate Study

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- **[P1 | HIGH] No full echoed-grid rerun with the compromise `X_pi` pulse yet**: The definitive pass now includes a standalone manifold-resolved `X_pi` audit and a lightweight compromise pulse, but it does not yet replay the full echoed construction grid with that updated pulse.
- **[P1 | HIGH] No new sensitivity sweep behind the practical ranking**: `F_practical` is currently a decoherence-only proxy. A real finite-difference drift analysis for the best direct and echoed constructions is still required.
- **[P1 | MEDIUM] No fresh GRAPE or unconstrained upper bound in this study folder**: The native-rich source study already gives strong direct positive cases, but the present definitive pass does not yet establish whether harder failures are ansatz-limited or physics-limited.

## Recommended Improvements (P2)
- **[P2 | MEDIUM] Extend the aggregation to `n_tr = 3` follow-ups**: The strongest strict-SQR and CPSQR cases should be replayed with a third transmon level before the conclusion is frozen.
- **[P2 | MEDIUM] Add explicit scaling follow-up for `N_active = 4, 5`**: The current report shows the qualitative breakdown of broad strict-SQR success, but not a dedicated fresh scaling campaign under one consistent direct ansatz.
- **[P2 | LOW] Improve the practical-ranking model**: Replace the current decoherence-only score with the full prompt definition once parameter-drift sweeps are available.

## Nice-to-Haves (P3)
- **[P3 | LOW] Add a figure for baseline-versus-native-rich waveform examples**: The current report emphasizes fidelity and error generators more than waveform shape comparisons.
- **[P3 | LOW] Add explicit external literature citations for selective conditioned rotations**: This pass focuses on repository-internal study provenance.

## Open Questions
- Are the remaining hard strict-SQR failures dominated by residual conditional phase, transverse under-rotation, or optimization-budget limits once a robust `X_pi` is available?
- Can the best direct strict-SQR construction remain above 0.99 under realistic drift in `chi`, qubit frequency, and drive amplitude?
- Does a truly unconstrained time-domain drive outperform the current direct multitone ansatz materially on the hard `chi_plus_chiprime`, `N_active = 3` cases?

## What Was Tried and Did Not Work
- The first pass at the definitive study attempted to use `tqdm` unconditionally in the new parallel helper. The environment did not have `tqdm` installed, so the helper was changed to degrade gracefully without it.
- The first generated master ranking table duplicated native-rich cases because the same case appeared in multiple source-study stages. The runner now deduplicates by `(study, case_id, construction)` before ranking.

## Compute & Resource Notes
- The artifact-driven definitive study rebuild (`python studies/the_definitive_ideal_sqr_gate_study/scripts/run_full_study.py --full`) completed in about 30-43 seconds on March 28, 2026, including the standalone refocusing-pulse scan, spectral-crowding analysis, and LaTeX rebuild.
- The current pass still does not launch a fresh long multitone optimization campaign; the heavy numerical work remains anchored in previously saved source-study artifacts.

## Resolved
- **Missing progress-bar dependency**: `parallel_utils.py` now falls back cleanly when `tqdm` is unavailable.
- **Duplicate native-rich rows in the ranking table**: `run_full_study.py` now deduplicates repeated `(study, case_id, construction)` rows before writing the definitive outputs.
