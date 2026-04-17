# Science Directive — Iteration 1
Date: 2026-03-28

## Problem Class
OPT | REP | ANA

## Scope
This iteration converts the initialized study folder into a functioning definitive-study aggregation pass. The work is intentionally anchored on the repository's existing SQR numerical corpus rather than silently launching a new multi-hour optimization campaign.

## Scientific Questions
1. Do the existing repository studies already contain strict ideal-SQR cases above 0.99?
2. If so, are those successes broad across the tested grid or restricted to selected addressed manifolds?
3. What error channels dominate the saved ideal-SQR failures when examined block by block?
4. Do echoed constructions help strict SQR, or mainly the relaxed CPSQR objective?

## Hypotheses
1. The original direct-vs-echoed ideal-SQR baseline remains a low-fidelity negative result.
2. The later native-rich follow-up already contains strict ideal-SQR cases above 0.99 on restricted addressed windows.
3. Echoed constructions remain much more compelling for relaxed CPSQR than for strict ideal SQR.
4. The dominant hard-case error is mixed transverse plus conditional-phase error, not a pure residual-Z term.

## Execution Plan
1. Normalize all four relevant SQR studies into one cross-study table.
2. Recompute per-manifold error-generator components for every saved ideal-SQR artifact that supports blockwise analysis.
3. Generate a unified prior-work figure set, case-construction heatmap, duration trends, and parameter-count comparison.
4. Build a definitive study report and reproducibility notebook from the aggregated results.
5. Record what remains missing for a true fresh optimization campaign: robust refocusing, sensitivity sweeps, and GRAPE bounds.

## Success Criteria
1. The study folder contains machine-readable prior-study snapshots, new aggregated results, validation JSON, figures, a compiled PDF, and a reproducibility notebook.
2. The report states clearly whether strict ideal SQR already exceeds 0.99 anywhere in the saved corpus.
3. The report distinguishes strict-SQR success from CPSQR success and backs that distinction with figures.
4. Remaining gaps are logged explicitly rather than implied away.

## Compute Plan
- Use artifact-level parallelism only.
- No fresh long `cqed_sim` optimization runs in this iteration.
- Expected runtime: under one minute end to end.
