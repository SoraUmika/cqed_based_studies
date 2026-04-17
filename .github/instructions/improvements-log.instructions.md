---
description: "IMPROVEMENTS.md structure enforcement. Ensures priority/difficulty tags, never-delete rule, and proper section organization."
applyTo: "studies/**/IMPROVEMENTS.md"
---

# Improvement Log Conventions

## Required Sections
Every IMPROVEMENTS.md must contain these sections:

1. `## Critical Gaps (P1)` — Things that could make results qualitatively wrong
2. `## Recommended Improvements (P2)` — Meaningful accuracy or scope improvements
3. `## Nice-to-Haves (P3)` — Lower-priority enhancements
4. `## Open Questions` — Unresolved physics or numerical observations
5. `## What Was Tried and Did Not Work` — Failed approaches with enough detail to understand WHY
6. `## Compute & Resource Notes` — Wall-clock times, memory usage, bottlenecks

## Tagging Rules
- Every improvement item must have a priority tag: **P1** (critical), **P2** (meaningful), **P3** (nice-to-have)
- Every improvement item must have a difficulty tag: **LOW**, **MEDIUM**, **HIGH**
- Format: `- **[P2 | MEDIUM]** Description of the improvement.`

## Never-Delete Rule
- Never delete entries from IMPROVEMENTS.md.
- Move resolved items to a `## Resolved` section at the bottom with a note on how they were fixed.

## Specificity Rule
- Be specific about failures. Not: "optimization didn't converge". Instead: "Nelder-Mead on 12-parameter pulse stalled at fidelity 0.987 after 500 iterations; cost landscape appears flat near the minimum; GRAPE or gradient-based method likely needed."

## Real-Time Updates
- Start this file in Step 1 (Initialize), even with placeholder headings.
- Update in real time during Steps 3-4 (Implement & Validate).
- Log limitations immediately when encountered — do not defer to the report phase.
