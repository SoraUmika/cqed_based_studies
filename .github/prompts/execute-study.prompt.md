---
description: "Launch the Execution Engineer to run a study phase (plan, implement, validate, report, or polish). Select Opus 4.6 in Claude Code BEFORE using this prompt."
name: "Execute Study Phase"
argument-hint: "study=studies/<name> run=task_runs/<slug> phase=plan|implement|validate|report|polish"
agent: "execution-engineer"
---
Execute a research study phase using the full cQED research loop protocol.

**Model:** Select **Opus 4.6** in Claude Code before invoking.

Inputs expected in the chat request:
- `study=studies/<name>` — path to the study directory
- `run=task_runs/<slug>` — path to the task run state directory
- `phase=plan|implement|validate|report|polish` — which phase to execute

Required behavior:
1. Read `research_config.json` to load loop configuration.
2. Read the study `README.md` and any existing state files (`SCIENCE_DIRECTIVE.md`, `EXECUTION_SUMMARY.md`, `REVIEW_DIRECTIVE.md`, `FOLLOWUP_PROMPT.md`).
3. Execute the requested phase following AGENTS.md protocol:
   - **plan**: Self-generate `SCIENCE_DIRECTIVE.md` with classified problem type, a first-principles analytic preliminary, explicit controlled approximations, hypotheses, experiment design, and success criteria.
   - **implement**: Write scripts, run simulations, save data/figures/artifacts; update `IMPROVEMENTS.md` in real time.
   - **validate**: Run all three validation checks (sanity, convergence, literature comparison).
   - **report**: Write `report/report.tex` with mandatory appendices; compile to PDF.
   - **polish**: Final readability pass after Science Director approval; write `POLISH_COMPLETE.md`.
4. After completing the report phase, write `REVIEW_REQUEST.md` and `EXECUTION_SUMMARY.md` in the run directory.
5. Update `study_state.json` and `PROGRESS_LOG.md` throughout.
