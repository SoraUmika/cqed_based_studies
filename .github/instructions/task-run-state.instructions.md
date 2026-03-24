---
description: "Use when creating or updating autonomous task-run state files such as EXECUTION_PLAN.md, TASK_CHECKLIST.md, PROGRESS_LOG.md, BLOCKERS.md, or DONE.md. Covers resumable checkpointing, append-only progress logging, blocker handling, and completion verification."
applyTo: "task_runs/**/*.md"
---
# Task-Run State File Rules

- Preserve the existing section headings when updating a state file unless a structural change is necessary.
- Treat `TASK_CHECKLIST.md` as the source of truth for completion state.
- Keep checklist task IDs stable once assigned.
- Append new entries to `PROGRESS_LOG.md`; do not rewrite prior checkpoints unless they are factually wrong.
- Keep `BLOCKERS.md` split into active and resolved blockers.
- Only create or update `DONE.md` after explicit final verification.
- When a plan changes, explain why in `PROGRESS_LOG.md`.
- If a task is partially complete but not verified, leave the checkbox unchecked and record details in `PROGRESS_LOG.md`.