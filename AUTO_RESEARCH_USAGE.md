# Auto-Research Workflow — Usage Guide

This guide explains how to run, configure, and recover the automated two-model
research loop for cQED studies.

---

## The single command

There is one command that covers all three cases — new study, resume, or
extend a completed study:

```powershell
.\tools\research_loop.ps1 -Action quickstart `
    -StudyName "my_study" `
    -StudyGoal "Your research question here"
```

It does three things automatically:

1. Creates the study directories and state files (or detects an existing/completed study)
2. Builds a full-context `@research-loop` prompt for Copilot Chat
3. **Copies the prompt to your clipboard**

Then: **Ctrl+V into Copilot Chat** → Enter. That is everything.

### The three cases

| Situation | What quickstart does |
|-----------|---------------------|
| Study does not exist | Creates all directories and files, builds a new-study prompt |
| Study in-progress | Detects the current phase, builds a recovery/resume prompt |
| Study is COMPLETE | Marks it for extension, appends your new goal, builds an extension prompt |

You can also run it from VS Code: `Ctrl+Shift+P` → **Tasks: Run Task** →
**Research: Quickstart (start or extend a study)**.

---

## Table of Contents

1. [Architecture — two models, two roles](#1-architecture--two-models-two-roles)
2. [One-time setup — configure models](#2-one-time-setup--configure-models)
3. [Starting a new study](#3-starting-a-new-study)
4. [Running the loop — step by step](#4-running-the-loop--step-by-step)
5. [VS Code tasks — point and click](#5-vs-code-tasks--point-and-click)
6. [Recovery — when an agent stops unexpectedly](#6-recovery--when-an-agent-stops-unexpectedly)
7. [Report preservation — how follow-up runs extend the report](#7-report-preservation--how-follow-up-runs-extend-the-report)
8. [Configuration reference](#8-configuration-reference)
9. [File reference](#9-file-reference)
10. [Troubleshooting](#10-troubleshooting)

---

## 1. Architecture — two models, two roles

The loop uses **two models in two distinct roles**:

| Role | Agent | Model | Responsible for |
|------|-------|-------|----------------|
| Science Director | `@science-director` | High-reasoning (e.g. o3, GPT-5.4 codex) | Physics reasoning, experiment design, result review, CONTINUE/VALIDATE/STOP decisions |
| Execution Engineer | `@execution-engineer` | Claude Opus 4.6 | Code, simulation, debugging, figures, report writing |

A third agent acts as a single-model fallback:

| Role | Agent | Model | When to use |
|------|-------|-------|-------------|
| Research Loop (unified) | `@research-loop` | Claude Opus 4.6 | Full autonomous loop in one session, or recovery |

The orchestration script (`tools/research_loop.ps1`) manages state and prints the
exact Copilot invocation to run at each step. You paste that invocation into
Copilot Chat.

```text
PLAN (Science Director)
  → IMPLEMENT (Execution Engineer)
  → REVIEW (Science Director)
     ├── CONTINUE / REVISE → back to IMPLEMENT
     └── VALIDATE → VALIDATE (Execution Engineer)
                 → REPORT  (Execution Engineer)
                 → COMPLETE
```

---

## 2. One-time setup — configure models

Open [research_config.json](research_config.json) and set the model IDs to match
the model names shown in your Copilot model picker:

```json
"models": {
    "planning":   { "model_id": "GPT-5.4"             },
    "execution":  { "model_id": "claude-opus-4-6" },
    "documentation": { "model_id": "claude-opus-4-6" }
}
```

> **Important:** these are documentation values, not automatic routing.
> You must manually select the correct model in the Copilot Chat model picker
> before invoking each agent. The config is what reminds you which one to pick.

To see the currently configured models:

```powershell
.\tools\research_loop.ps1 -Action config -StudyName <any-study-name>
```

Or use the VS Code task **Research: Show Loop Config**.

---

## 3. Starting a new study

### Step 1 — run quickstart

```powershell
.\tools\research_loop.ps1 -Action quickstart `
    -StudyName "my_new_study" `
    -StudyGoal "Optimize dispersive readout fidelity to >99.5%"
```

The script creates:

```text
studies/my_new_study/
    README.md
    IMPROVEMENTS.md
    study_state.json
    scripts/   data/   figures/   report/

task_runs/my_new_study/
    PROGRESS_LOG.md   BLOCKERS.md   TASK_CHECKLIST.md
    RESUME_PROMPT.md        ← pre-filled recovery prompt
```

It then prints and copies the Copilot invocation to your clipboard.

### Step 2 — open Copilot Chat

In VS Code: `Ctrl+Alt+I` (or the Copilot Chat icon in the sidebar).

Select the model, paste, press Enter.

---

## 4. Running the loop — step by step

Each step follows the same pattern:

1. Run the orchestration script to transition state and get the invocation
2. In Copilot Chat, **select the correct model**
3. Paste the invocation and press Enter
4. The agent writes its output file and updates state
5. Move to the next step

---

### Step A — PLAN (Science Director)

```powershell
.\tools\research_loop.ps1 -Action plan -StudyName my_new_study
```

In Copilot Chat — select **GPT-5.4** (or your reasoning model), then paste:

```text
@science-director study=studies/my_new_study run=task_runs/my_new_study phase=plan
```

The Science Director will:

- Read the study README and cqed_sim API
- Design experiments with specific parameters
- Write `task_runs/my_new_study/SCIENCE_DIRECTIVE.md`

---

### Step B — IMPLEMENT (Execution Engineer)

```powershell
.\tools\research_loop.ps1 -Action execute -StudyName my_new_study
```

In Copilot Chat — select **Claude Opus 4.6**, then paste:

```text
@execution-engineer study=studies/my_new_study run=task_runs/my_new_study phase=implement
```

The Execution Engineer will:

- Read `SCIENCE_DIRECTIVE.md`
- Write simulation scripts to `scripts/`, run them, save data and figures
- Debug failures (up to `max_retries_per_phase` attempts per task)
- Write `task_runs/my_new_study/EXECUTION_SUMMARY.md`

---

### Step C — REVIEW (Science Director)

```powershell
.\tools\research_loop.ps1 -Action review -StudyName my_new_study
```

In Copilot Chat — select **GPT-5.4** (or your reasoning model), then paste:

```text
@science-director study=studies/my_new_study run=task_runs/my_new_study phase=review
```

The Science Director reads `EXECUTION_SUMMARY.md` and makes one decision:

| Decision | Meaning | Next step |
|----------|---------|-----------|
| **CONTINUE** | Results on track but need more data | New directive → back to Step B |
| **REVISE** | Approach is wrong, need a different strategy | New directive → back to Step B |
| **VALIDATE** | Results look publication-quality | Proceed to Step D |
| **STOP** | Blocked on something only you can resolve | Loop stops |

The loop runs B → C until the decision is VALIDATE (or until `max_iterations`
is reached, after which VALIDATE is forced).

---

### Step D — VALIDATE (Execution Engineer)

```powershell
.\tools\research_loop.ps1 -Action validate -StudyName my_new_study
```

In Copilot Chat — select **Claude Opus 4.6**, then paste the invocation.

The agent runs sanity checks, convergence tests, and literature comparisons.

---

### Step E — REPORT (Execution Engineer)

```powershell
.\tools\research_loop.ps1 -Action report -StudyName my_new_study
```

In Copilot Chat — select **Claude Opus 4.6**, then paste the invocation.

The agent writes (or extends) `report/report.tex`, compiles the PDF, and sets
`study_state.json` status to `COMPLETE`.

---

### Checking status at any time

```powershell
.\tools\research_loop.ps1 -Action status -StudyName my_new_study
```

Shows: current status, iteration count, open tasks, active blockers, and the
exact Copilot invocation to run next.

---

### Single-agent shortcut — no model switching

To run the complete loop with one model and no manual switching:

```text
@research-loop study=studies/my_study goal="Your research question here"
```

Or to resume an existing study:

```text
@research-loop study=studies/my_study resume
```

This is what `quickstart` uses. It is less optimal for planning quality (the
single model does both physics reasoning and implementation) but is fully
autonomous.

---

## 5. VS Code tasks — point and click

Open the task runner: `Ctrl+Shift+P` → **Tasks: Run Task**, then select from:

| Task | What it does |
|------|-------------|
| **Research: Quickstart (start or extend a study)** | The single command — init or detect extension, copy prompt to clipboard |
| **Research: New Study** | Runs `init` only — prompts for study name and goal |
| **Research: Study Status** | Runs `status` — shows current state |
| **Research: Resume Study** | Runs `resume` — prints the next Copilot invocation |
| **Research: Generate Recovery Prompt** | Runs `recover` — writes `RESUME_PROMPT.md`, copies invocation to clipboard |
| **Research: Show Loop Config** | Runs `config` — prints all active settings |
| **Research: Open Loop Config** | Opens `research_config.json` in the editor |
| **Research: Run Loop Action** | Runs any action — shows a pick list |

---

## 6. Recovery — when an agent stops unexpectedly

Copilot agents can stop for many reasons: context limit, network error, VS Code
crash, session timeout. The loop is designed to survive all of these.

### Option A — VS Code task (easiest)

1. Run the task **Research: Generate Recovery Prompt**
2. Enter the study name when prompted
3. The agent invocation is now in your clipboard
4. Paste it into Copilot Chat and press Enter

### Option B — terminal

```powershell
.\tools\research_loop.ps1 -Action recover -StudyName my_new_study
```

This writes `task_runs/my_new_study/RESUME_PROMPT.md` with:

- The current loop state (status, iteration, objective)
- The exact `@agent` invocation to paste
- A list of all context files the agent must read
- Active blockers and config settings

The agent invocation is also copied to your clipboard.

### Option C — quickstart (also works for recovery)

```powershell
.\tools\research_loop.ps1 -Action quickstart -StudyName my_new_study
```

Detects that the study is in-progress and builds a recovery prompt automatically.

### What the agent does on recovery

When any agent sees `RESUME_PROMPT.md` in the run directory it:

1. Reads all state files (study_state.json, SCIENCE_DIRECTIVE.md, EXECUTION_SUMMARY.md,
   TASK_CHECKLIST.md, BLOCKERS.md)
2. Announces: "RESUMING from [phase] — iteration N/max. Completed tasks: X. Open tasks: Y."
3. Continues from the current phase
4. Skips all tasks already marked `[x]` in TASK_CHECKLIST.md

### What NOT to do on recovery

- Do NOT run `init` again — it will overwrite your state
- Do NOT restart from scratch
- Do NOT manually edit `study_state.json` unless you know exactly what you are doing

### If the study is BLOCKED

Check blockers with:

```powershell
.\tools\research_loop.ps1 -Action status -StudyName my_new_study
```

Most blockers are framework limitations documented in `BLOCKERS.md`. You can:

- Resolve the blocker and then run `recover`
- Set `blocked_phase_policy = "continue_with_partial"` in `research_config.json`
  (this is the default) to skip past blocked tasks and continue
- Manually set `study_state.json → status` back to `IMPLEMENTING` and run `recover`

---

## 7. Report preservation — how follow-up runs extend the report

By default (`report.preserve_existing_report = true`), follow-up runs **never
overwrite** `report.tex`. Instead the agent:

1. Reads the existing `report.tex` in full
2. Creates a backup: `report.tex.bak`
3. Inserts a new section **before** `\end{document}`:

   ```latex
   % ===== Research Extension — Iteration N =====
   \clearpage
   \section{Extension: <Title> (Iteration N)}
   ...new results, figures, analysis...
   ```

4. Appends new bibliography entries to `references.bib` (no duplicates)
5. Leaves the abstract unchanged after the first follow-up

Each follow-up run adds a new chapter to the same document.

To disable and allow overwriting:

```json
"report": { "preserve_existing_report": false }
```

---

## 8. Configuration reference

All settings live in [research_config.json](research_config.json).
Edit it with the task **Research: Open Loop Config**.

```json
{
  "models": {
    "planning":      { "model_id": "o3"              },
    "execution":     { "model_id": "claude-opus-4-6" },
    "documentation": { "model_id": "claude-opus-4-6" }
  },

  "loop": {
    "max_iterations": 10,
    "min_iterations": 1
  },

  "retry": {
    "auto_recover_on_failure": true,
    "max_retries_per_phase": 3,
    "blocked_phase_policy": "continue_with_partial"
  },

  "report": {
    "preserve_existing_report": true,
    "extension_mode": "append_iteration_section",
    "backup_before_write": true
  },

  "logging": {
    "generate_recovery_prompt": true,
    "recovery_prompt_file": "RESUME_PROMPT.md",
    "max_execution_summary_lines": 500
  }
}
```

| Key | Default | Meaning |
|-----|---------|---------|
| `models.planning.model_id` | `o3` | Select this in Copilot before invoking `@science-director` |
| `models.execution.model_id` | `claude-opus-4-6` | Select this before invoking `@execution-engineer` |
| `loop.max_iterations` | `10` | Hard cap on PLAN→IMPLEMENT→REVIEW cycles; forces VALIDATE at the limit |
| `retry.max_retries_per_phase` | `3` | Per-task debug attempts inside the Execution Engineer |
| `retry.blocked_phase_policy` | `continue_with_partial` | `continue_with_partial` skips blocked tasks; `stop_and_report` halts |
| `report.preserve_existing_report` | `true` | Never overwrite report.tex; always extend |
| `report.backup_before_write` | `true` | Saves report.tex.bak before any write |

---

## 9. File reference

### Per-study files

| File | Location | Purpose |
|------|----------|---------|
| `study_state.json` | `studies/<name>/` | Machine-readable state: status, iteration, tasks, results |
| `README.md` | `studies/<name>/` | Human-readable study overview |
| `IMPROVEMENTS.md` | `studies/<name>/` | Gaps, limitations, open questions — drives extension passes |
| `report/report.tex` | `studies/<name>/` | LaTeX report (extended, never overwritten by default) |

### Per-run coordination files

| File | Location | Written by | Read by |
|------|----------|-----------|---------|
| `SCIENCE_DIRECTIVE.md` | `task_runs/<name>/` | Science Director | Execution Engineer |
| `EXECUTION_SUMMARY.md` | `task_runs/<name>/` | Execution Engineer | Science Director |
| `TASK_CHECKLIST.md` | `task_runs/<name>/` | Execution Engineer | Both agents |
| `PROGRESS_LOG.md` | `task_runs/<name>/` | Both agents | Recovery |
| `BLOCKERS.md` | `task_runs/<name>/` | Execution Engineer | Both agents |
| `RESUME_PROMPT.md` | `task_runs/<name>/` | PS1 script | Both agents on recovery |

---

## 10. Troubleshooting

### "The property 'Count' cannot be found"

Stale PowerShell strict-mode bug. Pull the latest `tools/research_loop.ps1`.

### Agent says "I don't see any SCIENCE_DIRECTIVE.md"

The plan phase was not completed. Run:

```powershell
.\tools\research_loop.ps1 -Action plan -StudyName <name>
```

Then invoke `@science-director`.

### Agent rewrote my report.tex

Check `research_config.json`: `report.preserve_existing_report` must be `true`.
Restore from `report.tex.bak` if needed.

### Loop stuck — Science Director keeps saying CONTINUE

Lower `loop.max_iterations` in `research_config.json` (e.g. to 5) so the loop
is forced to VALIDATE sooner. You can review the study state and increase it
again if more iterations are genuinely needed.

### Study is COMPLETE but I want to run a new extension

Do NOT run `init` again. Instead just run quickstart with a new goal:

```powershell
.\tools\research_loop.ps1 -Action quickstart `
    -StudyName "hybrid_qubit_cavity_control" `
    -StudyGoal "Investigate thermal photon effects on the best SNAP and SQR gates"
```

quickstart detects the COMPLETE status, marks the study for extension, appends
the new goal, and copies the Copilot prompt to your clipboard. The report phase
will add a new section to `report.tex` rather than overwriting it.

### I want one model for everything — no switching

Set both `models.planning.model_id` and `models.execution.model_id` to the same
model. Or just use `@research-loop` directly (what `quickstart` does).

### Extending the iteration limit mid-study

Edit `loop.max_iterations` in `research_config.json` and run:

```powershell
.\tools\research_loop.ps1 -Action recover -StudyName <name>
```

The next invocation will use the updated limit.
