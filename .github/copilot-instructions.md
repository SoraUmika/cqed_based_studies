---
description: "Global Copilot instructions for the cQED Autonomous Research Platform. Applied to every agent invocation. Provides baseline context: repo purpose, key paths, conventions, and non-negotiable rules."
---

# cQED Autonomous Research Platform — Global Instructions

## Repo Purpose

This repository hosts autonomous cQED (circuit quantum electrodynamics) simulation studies using the `cqed_sim` framework. Every study follows a structured lifecycle: Initialize → Plan → Implement → Validate → Report, executed by a two-model agent loop (Opus 4.6 as Execution Engineer, Codex 5.4 xHigh as Science Director).

## Key Paths

| What | Where |
|------|-------|
| Framework spec | `AGENTS.md` |
| Research loop architecture | `RESEARCH_LOOP.md` |
| Loop configuration | `research_config.json` |
| Study root | `studies/<study_name>/` |
| Task run state | `task_runs/<study_name>/` |
| Orchestration tools | `tools/*.ps1` |
| cqed_sim API Reference | [API_REFERENCE.md](https://github.com/SoraUmika/qubox_cQEDsim/blob/main/API_REFERENCE.md) |
| cqed_sim local copy | `C:\Users\dazzl\Box\Shyam Shankar Quantum Circuits Group\Users\Users_JianJun\cQED_simulation` |

## Environment

- **Python**: 3.12.10 (system — no venvs, no conda)
- **OS**: Windows 11 with PowerShell 5.1
- **Package install**: `pip install <package> --user` (log in IMPROVEMENTS.md)
- **Core packages**: numpy, scipy, matplotlib, qutip 5.x, lmfit, seaborn

## Non-Negotiable Rules

1. Always use `cqed_sim` — no ad-hoc simulation code unless a gap is documented.
2. Never skip workflow steps (Initialize → Plan → Implement → Validate → Report).
3. Stop immediately on any `cqed_sim` inconsistency — log in BLOCKERS.md.
4. No virtual environments — use system Python directly.
5. Every study needs `scripts/reproducibility_notebook.ipynb`.
6. Consult the API Reference before writing any simulation code.
7. Start from first principles whenever feasible. For physics, math, or study questions, first formulate the minimal analytic model, then introduce only controlled approximations with stated validity before moving to numerics or code.

## Default Reasoning Pattern

For study work and technical questions, prefer this order unless the user explicitly wants something else:

1. First-principles model or governing equations
2. Minimal controlled approximations, with validity conditions stated
3. Numerical simulation, code, or implementation details

## Naming Conventions

- Study folders: `studies/<lowercase_descriptive_name>/`
- Scripts: `snake_case.py`
- Figures: `<figure_description>.{png,pdf}` (save both formats, 300 dpi PNG)
- Artifacts: JSON for parameters/metadata, NPZ for arrays, CSV for tabular data

## Report Standards

- Use `revtex4-2` two-column format with `siunitx` for units
- No filenames or `snake_case` identifiers in the main text prose
- Every claim backed by a figure, table, or artifact
- Mandatory sections: Abstract, Introduction, System & Methods, Results, Validation, Discussion, Conclusion, Limitations & Future Work, Appendix (Detailed Results + Reproducibility)

## PowerShell 5.1 Quirks

- `Join-Path` only takes 2 args — nest calls: `Join-Path (Join-Path $a $b) $c`
- Keep `.ps1` files ASCII-only (no Unicode arrows or box-drawing characters)
- Use `@()` wrapper for `.Count` on potentially null values

## Context Management

Long agent sessions degrade in quality as the context window fills. Prefer **context resets** over compaction.

**Checkpoint-and-reset protocol:**
1. When you notice coherence degrading (repeating yourself, losing track of state, or exceeding ~60% of context) — checkpoint.
2. Save all progress to persistent state files: update `PROGRESS_LOG.md` (append-only), `TASK_CHECKLIST.md` (mark completed items), `study_state.json`, and any `BLOCKERS.md`.
3. Write `RESUME_PROMPT.md` via `tools/research_loop.ps1 -Action recover` (or manually).
4. Signal the user to start a fresh session with the recovery prompt.

**What to persist at checkpoints:**
- Completed task list with quantitative results achieved
- Active blockers with diagnostic detail
- Key numerical results (fidelities, parameters, timings)
- Current phase and next action

**What NOT to carry forward:**
- Full file contents already saved to disk
- Intermediate debugging output
- Exploratory dead ends (summarize in PROGRESS_LOG instead)

## Cross-Study Learning

Before starting any new study, read `LESSONS_LEARNED.md` at the repo root. It contains cross-study insights about failed approaches, cqed_sim quirks, parameter ranges, and compute timing data that prevent repeating past mistakes.
