# Science Directive - Iteration 1
Date: 2026-04-06
Study: `studies/multitone_sqr_echo_rigorous_followup`
Run: `task_runs/multitone_sqr_echo_rigorous_followup`

## Classification
- Primary class: `ANA`
- Secondary classes: `REP`, `DES`

## Core Question
If the echoed multitone sequence is treated as an actual optimizable ansatz rather than a replay of a separately optimized half target, does it rescue the strict simultaneous shared-line no-detuning SQR problem, or does it still fail once judged on phase-sensitive blockwise metrics?

## Hypotheses
1. Jointly optimizing the two echoed halves will improve over inherited-half replay, especially for aligned-`x` targets, but will still not recover the ideal gate generically.
2. Ideal instantaneous refocusing is an optimistic upper bound; finite refocusing will reintroduce manifold-dependent errors that spoil exact cancellation.
3. A manifold-aware shared-line refocusing pulse may reduce the gap relative to vacuum-calibrated Gaussian refocusing, but it will not change the exact strict verdict.

## Ordered Tasks
1. Reproduce the strict direct baseline and decoupled-block comparator.
2. Implement explicit metric definitions and machine-readable metric payloads.
3. Implement and optimize echoed ansatz variants:
   - ideal instantaneous `X_pi`,
   - finite Gaussian `X_pi`,
   - manifold-aware shared-line multitone `X_pi` if tractable.
4. Add a duration-matched direct baseline to remove unfair timing comparisons.
5. Validate the most favorable aligned-`x` cases carefully before writing conclusions.
6. Regenerate the report, notebook, execution summary, and review request.

## Acceptance Criteria
1. The report must say exactly which echo variants were optimized and which were only replayed.
2. The metrics section must define every scalar shown in plots and tables.
3. Any negative claim about the echo must be supported by optimized echo data, not only inherited-half replay.
4. The final verdict must separate exact impossibility, practical failure under tested ansatzes, and remaining unresolved possibilities.
