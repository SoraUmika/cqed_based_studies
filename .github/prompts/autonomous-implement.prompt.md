---
description: "Execute the next unfinished slice of a planned Copilot task and update file-based checkpoints."
name: "Autonomous Implement"
argument-hint: "task=RESEARCH_PLAN.md run=task_runs/research_plan"
agent: "autonomous-implementer"
---
Execute unfinished checklist items from a task run directory.

Inputs expected in the chat request:
- `task=<path to the source task file>`
- `run=<path to the run directory>`

Required behavior:
- Re-read the source task and all run-state files before making changes.
- Complete one checkpoint or a small cluster of related checklist items.
- Update `PROGRESS_LOG.md`, `TASK_CHECKLIST.md`, and `BLOCKERS.md` before stopping.
- Write `DONE.md` only after explicit final verification.