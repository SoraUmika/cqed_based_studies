---
description: "Launch the Research Loop orchestrator to run a complete study autonomously in single-agent mode. Works on whichever model is currently selected."
name: "Run Full Research Loop"
argument-hint: "study=studies/<name> goal='Research question here' OR study=studies/<name> resume"
agent: "research-loop"
---
Run the full 4-stage cQED research loop in single-agent mode.

Inputs expected in the chat request:
- `study=studies/<name>` — path to the study directory
- `goal='<research question>'` — for new studies
- `resume` — to continue an interrupted study

Required behavior for new studies:
1. Initialize the study folder structure (README, IMPROVEMENTS.md, scripts/, data/, figures/, artifacts/, report/).
2. Self-generate the research plan (SCIENCE_DIRECTIVE.md) with classified problem type, a first-principles analytic model, and explicit controlled approximations before numerics.
3. Execute: write scripts, run simulations, save data/figures/artifacts.
4. Validate: all three checks (sanity, convergence, literature).
5. Write the LaTeX report with mandatory appendices; compile to PDF.
6. Switch to reviewer hat: critically evaluate the report across all four dimensions.
7. If deficiencies found: write follow-up tasks and iterate (up to `max_iterations` from research_config.json).
8. Once approved: perform final readability polish.
9. Create the reproducibility notebook.
10. Write POLISH_COMPLETE.md and mark study COMPLETE.

Required behavior for resume:
1. Read all state files (study_state.json, PROGRESS_LOG.md, TASK_CHECKLIST.md, BLOCKERS.md).
2. Determine the current phase and last completed checkpoint.
3. Continue from where the study was interrupted.
