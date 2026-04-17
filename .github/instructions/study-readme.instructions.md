---
description: "Study README structure enforcement. Ensures all mandatory sections from AGENTS.md are present and properly formatted."
applyTo: "studies/**/README.md"
---

# Study README Conventions

## Required Sections (all must be present)
Every study README.md must contain these sections in order:

1. `# <Study Title>` — Top-level heading
2. `## Problem Class` — One or more of: OPT, REP, DES, ANA
3. `## Motivation` — Why this study matters; link to paper if REP class
4. `## Goals` — Numbered, concrete, falsifiable goals
5. `## Methods` — Which cqed_sim modules/functions will be used
6. `## Analytic Preliminary` — First-principles model and closed-form or limiting-case reasoning attempted before numerics; list controlled approximations and why they are valid (or explain why no useful analytic foothold exists)
7. `## cqed_sim Gap Analysis` — Table: Functionality | Needed? | Available? | Plan
8. `## Assumptions` — Physical assumptions, parameter ranges, convergence criteria
9. `## Compute & Resource Strategy` — Cost estimates, acceleration plans, bottlenecks
10. `## Expected Outcomes` — Quantitative success criteria
11. `## Known Limitations` — Updated throughout the study
12. `## Validation` — Three checkboxes that must be kept current:
    - `- [ ] Sanity checks`
    - `- [ ] Convergence`
    - `- [ ] Literature comparison (if applicable)`
13. `## Status` — One of: ACTIVE | COMPLETE | BLOCKED

## Status Updates
- Update `## Status` as work progresses.
- Mark validation checkboxes as `[x]` only when the check has actually passed with documented evidence.
- Do not mark COMPLETE unless: all validation checks pass, report appendices are present, and at least one machine-readable artifact exists per headline result.

## Study Composition (optional)
- If the study consolidates or extends earlier work, include `## Study Composition` mapping each inherited component to its original study, key result, and role.

## Analytic Preliminary Quality Bar
- Start from the first-principles model whenever feasible, not just a heuristic summary.
- State each controlled approximation explicitly and give its validity condition or regime.
- If the problem truly has no useful analytic foothold, say why before proceeding numerically.
