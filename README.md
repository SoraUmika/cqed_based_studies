# Resumable Copilot Workflow

This workspace includes a VS Code-native Copilot workflow for long, multi-step tasks that need explicit state, resumability, and bounded continuation.

## Core Pieces
- Custom agents in `.github/agents/` for planning, implementation, and resume.
- Slash prompts in `.github/prompts/` for repeatable invocation from the chat UI.
- State files under `task_runs/<run-name>/` so the workflow survives interrupted sessions.
- A PowerShell helper in `tools/copilot_task_run.ps1` for bootstrapping and status checks.
- VS Code tasks in `.vscode/tasks.json` for launching the helper from the Command Palette.

## Typical Flow
1. Run `Tasks: Run Task` and choose `Copilot: Init Task Run`.
2. In chat, run `/Autonomous Plan` with `task=RESEARCH_PLAN.md run=task_runs/research_plan`.
3. Run `/Autonomous Implement` to execute the next checkpoint.
4. If the run stops midway, run `/Autonomous Resume` with the same inputs.
5. Treat `DONE.md` as the verified completion marker.

## State Files
- `EXECUTION_PLAN.md`: phases, deliverables, success criteria, and checkpoint boundaries.
- `TASK_CHECKLIST.md`: the execution source of truth.
- `PROGRESS_LOG.md`: append-only checkpoints and validation notes.
- `BLOCKERS.md`: active and resolved blockers.
- `DONE.md`: final verification report created only when the run is actually complete.

## Limitations
- Current VS Code Copilot workflows are not a reliable infinite unattended loop.
- The robust pattern is bounded execution slices plus explicit resume from files.
- Tasks can scaffold and inspect run state, but the actual Copilot invocation still happens through chat or agent selection.