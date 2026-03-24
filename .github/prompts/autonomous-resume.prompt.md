---
description: "Resume a previously interrupted Copilot task from file-based state without redoing completed work."
name: "Autonomous Resume"
argument-hint: "task=RESEARCH_PLAN.md run=task_runs/research_plan"
agent: "autonomous-resume"
---
Resume a stopped run from its state files.

Inputs expected in the chat request:
- `task=<path to the source task file>`
- `run=<path to the run directory>`

Required behavior:
- Read the source task plus `EXECUTION_PLAN.md`, `TASK_CHECKLIST.md`, `PROGRESS_LOG.md`, and `BLOCKERS.md`.
- Reconstruct what is complete, what remains, and what was last attempted.
- Continue from the next unfinished checkpoint.
- Avoid redoing checked items.
- Stop cleanly on real blockers and document them.