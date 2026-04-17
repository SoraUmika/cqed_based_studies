# Review Request
Date: 2026-04-06
Study: `studies/multitone_sqr_echo_rigorous_followup`
Run: `task_runs/multitone_sqr_echo_rigorous_followup`
Status: READY_FOR_REVIEW

## What Changed Relative To The Earlier Strict Study
1. The echo section is no longer inherited half-target replay only.
2. The new study optimizes the echoed ansatz directly, with independent segment corrections and a toggling-consistent second-half seed.
3. The metric suite is now explicit about phase-sensitive and probe-state diagnostics.
4. A total-duration-matched direct comparator and a representative manifold-aware refocusing benchmark were added.

## What This Follow-Up Now Claims
1. The strict simultaneous shared-line no-detuning problem still does not realize the ideal gate exactly.
2. The decoupled-block success remains exact but is still a different physical problem.
3. Replayed ideal echo can look falsely good on residual-`Z` alone, so that metric must not be used by itself.
4. Jointly optimized ideal instantaneous echo is only an upper bound: it can help some long-duration special cases, but it never realizes the ideal gate and does not win overall.
5. Physical finite echoed constructions, including the manifold-aware refocusing benchmark, still do not rescue the gate.

## Requested Review Focus
Please review the package with particular attention to:
1. Whether the report is careful enough in distinguishing ideal instantaneous echo from the physically implemented finite echoes.
2. Whether the metric discussion is clear that residual-`Z` alone can hide failure.
3. Whether the finite Gaussian and multitone-refocus data support the practical negative verdict without overstating the idealized upper-bound result.
4. Whether the decoupled-block exact success is kept clearly separate from the physical shared-line problem.
