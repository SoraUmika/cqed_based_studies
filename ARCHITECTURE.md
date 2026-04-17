# Architecture Overview

> Quick-parse entry point for agents. For full protocol details see `AGENTS.md`.
> For loop state management see `RESEARCH_LOOP.md`.

---

## Repository Layout

```
cqed_based_study/
├── AGENTS.md                 # Full agent protocol (600+ lines) — the source of truth
├── ARCHITECTURE.md           # ← YOU ARE HERE — fast orientation for agents
├── RESEARCH_LOOP.md          # Four-stage loop state management details
├── LESSONS_LEARNED.md        # Cross-study lessons — READ BEFORE starting any new study
├── research_config.json      # Loop configuration (models, iterations, retry, report settings)
│
├── studies/                  # One subdirectory per study
│   └── <study_name>/
│       ├── README.md         # Study scope, goals, methods, status
│       ├── IMPROVEMENTS.md   # Limitation log (never delete entries)
│       ├── study_state.json  # Machine-readable loop state
│       ├── scripts/          # Python simulation scripts + reproducibility_notebook.ipynb
│       ├── data/             # Raw simulation output (NPZ, CSV)
│       ├── artifacts/        # Processed results (JSON metadata, NPZ arrays)
│       ├── figures/          # Plots in both .png (300 dpi) and .pdf (vector)
│       └── report/           # report.tex, references.bib, report.pdf
│
├── task_runs/                # Per-study orchestration state
│   └── <study_name>/
│       ├── SCIENCE_DIRECTIVE.md    # Research plan (Opus writes)
│       ├── EXECUTION_SUMMARY.md    # Results digest (Opus writes)
│       ├── REVIEW_REQUEST.md       # Signals readiness for review
│       ├── REVIEW_DIRECTIVE.md     # Reviewer verdict + required actions (Codex writes)
│       ├── FOLLOWUP_PROMPT.md      # Revision instructions (Codex writes)
│       ├── TASK_CHECKLIST.md       # Resumable task tracking
│       ├── PROGRESS_LOG.md         # Append-only log
│       ├── BLOCKERS.md             # Active and resolved blockers
│       └── RESUME_PROMPT.md        # Auto-generated recovery prompt
│
├── tools/                    # PowerShell orchestration
│   ├── research_loop.ps1     # State machine driver (init, status, execute, review, recover, quickstart)
│   ├── auto_loop.ps1         # Background watcher — auto-fires agents on signal file changes
│   ├── copilot_task_run.ps1  # Lightweight task-run pattern for ad-hoc work
│   └── validate_study.ps1    # Automated pre-review structural validator
│
├── .github/
│   ├── agents/               # Agent definitions (.agent.md)
│   │   ├── execution-engineer.agent.md    # Opus 4.6 — full study lifecycle
│   │   ├── science-director.agent.md      # Codex 5.4 xHigh — critical reviewer
│   │   ├── research-loop.agent.md         # Single-agent fallback
│   │   └── autonomous-*.agent.md          # Resumable task agents
│   ├── instructions/         # File-pattern-scoped instructions (.instructions.md)
│   │   ├── python-study-code.instructions.md   # Python script conventions
│   │   ├── latex-report.instructions.md        # LaTeX report rules
│   │   ├── study-readme.instructions.md        # README section enforcement
│   │   ├── improvements-log.instructions.md    # IMPROVEMENTS.md structure
│   │   └── task-run-state.instructions.md      # State file conventions
│   ├── prompts/              # Reusable prompt templates (.prompt.md)
│   ├── skills/               # Domain-specific skill packages
│   │   ├── cqed-sim-lookup/           # API reference lookup before simulation code
│   │   ├── study-init/                # Study folder scaffolding
│   │   ├── study-validator/           # Pre-review structural completeness check
│   │   ├── red-green-validation/      # Test-first validation (Red/Green TDD)
│   │   ├── parallel-sweep/            # Parallel subagent orchestration for sweeps
│   │   ├── validate-results/          # Three-check validation gate
│   │   ├── latex-report/              # Report generation and compilation
│   │   ├── report-preflight/          # Pre-compilation lint scan
│   │   ├── report-review/             # Critical review protocol
│   │   ├── publication-figures/       # Matplotlib style and palettes
│   │   ├── reproducibility-notebook/  # End-to-end notebook creation
│   │   ├── cross-study-analysis/      # Multi-study comparison
│   │   └── ...                        # See .github/skills/ for full list
│   └── copilot-instructions.md  # Global instructions applied to every agent
│
└── schemas/
    └── study_state.schema.json  # JSON Schema for study_state.json
```

---

## Two-Model Agent Loop

```
  Opus 4.6 (Claude Code)              Codex 5.4 xHigh (GitHub Copilot)
  ─────────────────────               ──────────────────────────────────
  PLAN → IMPLEMENT → VALIDATE         
           → REPORT →                    REVIEW (APPROVE / REVISE / NEEDS_REWORK)
                                           │
                       ┌───────────────────┘
                       ↓
           REFINE (if REVISE)  ──→       RE-REVIEW
                                           │
                       ┌───────────────────┘
                       ↓
           POLISH (after APPROVE)  →     COMPLETE
```

Signal files drive transitions. See `RESEARCH_LOOP.md` for the full state machine.

---

## cqed_sim Integration

| Resource | Location |
|----------|----------|
| Source (GitHub) | [SoraUmika/qubox_cQEDsim](https://github.com/SoraUmika/qubox_cQEDsim) |
| API Reference | [API_REFERENCE.md](https://github.com/SoraUmika/qubox_cQEDsim/blob/main/API_REFERENCE.md) |
| Local copy | `C:\Users\dazzl\Box\...\cQED_simulation` |

**Usage rule:** Always prefer `cqed_sim`. Document any gap in the study README before writing standalone code. See AGENTS.md §7 for Path A (upstream) vs Path B (local) gap handling.

---

## Environment

- **Python** 3.12.10 (system — no venvs)
- **OS** Windows 11, PowerShell 5.1
- **Core packages** numpy, scipy, matplotlib, qutip 5.x, lmfit, seaborn
- **Install** `pip install <pkg> --user` — log in IMPROVEMENTS.md

---

## Quick Reference for New Agents

1. Read `ARCHITECTURE.md` (this file) for orientation
2. Read `LESSONS_LEARNED.md` for cross-study insights
3. Read `AGENTS.md` for the full protocol
4. Read the study `README.md` for scope and goals
5. Read `research_config.json` for loop settings
6. Use the `cqed-sim-lookup` skill before writing simulation code
