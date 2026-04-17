# Improvement Log: Realistic Universal Hybrid Control in the Dispersive cQED Regime

> Written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- **[P1 | HIGH] No fully pulse-backed non-GRAPE universal stack yet**: The synthesis supports a phase-aware primitive library, but not a complete pulse-backed universal architecture. What to do: build and validate one explicit composite stack using branch-compensated displacement, spectator-limited qubit pulses, relaxed selective control, and gauge cleanup without hiding the local layers inside unconstrained GRAPE.
- **[P1 | MEDIUM] Strict-joint composite optimization is still missing**: Echoed and segmented families look excellent under relaxed metrics but fail the strict full-joint metric on hard cases. What to do: optimize directly against the full strict joint operator including inter-manifold phase structure.
- **[P1 | MEDIUM] Higher-level and noisy validation are not yet universalized**: Many inherited strict-control studies remain closed-system and optimize with `n_tr = 2`. What to do: replay the strongest strict and relaxed cases with `n_tr = 3`, Kerr, and noise inside the optimization loop.

## Recommended Improvements (P2)
- **[P2 | MEDIUM] Compare SNAP cleanup against virtual-Z cleanup**: The synthesis argues for a phase-aware library, but the cheapest cleanup layer is not yet settled. What to do: quantify the cost and fidelity tradeoff between explicit SNAP cleanup and virtual-Z style cleanup after the best relaxed selective pulses.
- **[P2 | MEDIUM] Extend the synthesis to non-Fock logical encodings**: The current verdict is based on low-dimensional Fock-window evidence. What to do: repeat the same claim-boundary synthesis for cat and binomial encodings.
- **[P2 | LOW] Replace narrative source pulls with more structured machine-readable summaries**: A few inherited conclusions still come from markdown summaries rather than dedicated JSON schemas. What to do: upstream compact summary JSON outputs to the inherited studies when they are next touched.

## Nice-to-Haves (P3)
- **[P3 | LOW] Broaden the analytic phase-budget figure**: The current plot already captures the main timing conflict, but it could include branch-resolved cavity displacement mismatch more explicitly.
- **[P3 | LOW] Add a report-side appendix table of inherited study paths and exact artifact roles**: The current appendix is concise. A fuller audit table would help future meta-studies.

## Open Questions
- Is the present strict-control failure fundamentally tied to conditional-Z structure, or could a richer but still interpretable composite family restore the missing inter-manifold phase relations?
- How far can the surviving phase-aware primitive library be pushed before the cleanup overhead erases its practical advantage over a short GRAPE-assisted sequence?
- Does the balance between branch-compensated displacement and relaxed selective control change qualitatively in larger bosonic encodings?

## What Was Tried and Did Not Work
- **Literal inheritance of the ideal primitive story**: The synthesis did not support the claim that the three abstract primitives survive unchanged once realistic dispersive dynamics are enforced. Why it failed: the same `chi` that enables selectivity also obstructs unconditional control.
- **Naive echo as a universal fix**: The inherited evidence shows that vacuum-calibrated echo constructions can improve relaxed conditional-phase metrics without restoring strict full-joint control. Why it failed: the inserted refocusing pulse is itself manifold dependent.
- **Treating unconstrained GRAPE as the conceptual answer**: This was intentionally not adopted because it would obscure the physical claim boundary the user asked for. Why it was rejected: the study objective was to understand the structured control problem itself, not to replace it by a black-box optimum.

## Resolved
- **Resolved source integration**: The synthesis now consolidates the inherited displacement, waveform-level, selective-pulse, strict-SQR, relaxed-CPSQR, arbitrary-control, and runtime sequence evidence into one machine-readable summary plus report.

## Compute & Resource Notes
- No new packages were installed.
- The synthesis dataset build and figure generation completed quickly on CPU because the workflow reused validated inherited artifacts rather than rerunning the full repository.
- The report was built with explicit `pdflatex -> bibtex -> pdflatex -> pdflatex` passes.
