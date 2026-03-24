---
description: "Create or refresh a resumable execution plan and checklist for a large task document such as RESEARCH_PLAN.md."
name: "Autonomous Plan"
argument-hint: "task=RESEARCH_PLAN.md run=task_runs/research_plan"
agent: "autonomous-planner"
---
Plan a long-running task using file-based state instead of chat-only state.

Inputs expected in the chat request:
- `task=<path to the source task file>`
- `run=<path to the run directory>`

Required behavior:
- Read the source task document first.
- Create or refresh `EXECUTION_PLAN.md` and `TASK_CHECKLIST.md` in the run directory.
- Ensure `PROGRESS_LOG.md` and `BLOCKERS.md` exist.
- Keep planning separate from implementation.
- Define phases, deliverables, validation gates, checkpoints, and completion criteria.
- Preserve existing progress when rerun.