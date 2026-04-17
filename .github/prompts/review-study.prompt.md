---
description: "Launch the Science Director to critically review a completed study report. Select Codex 5.4 xHigh in GitHub Copilot BEFORE using this prompt."
name: "Review Study Report"
argument-hint: "study=studies/<name> run=task_runs/<slug> phase=review"
agent: "science-director"
---
Critically review a completed cQED study report as Science Director.

**Model:** Select **Codex 5.4 xHigh** in GitHub Copilot before invoking.

Inputs expected in the chat request:
- `study=studies/<name>` — path to the study directory
- `run=task_runs/<slug>` — path to the task run state directory
- `phase=review`

Required behavior:
1. Read `research_config.json` to check iteration limits and review configuration.
2. Read `REVIEW_REQUEST.md` and `EXECUTION_SUMMARY.md` from the run directory.
3. Read the full `report/report.tex` (not just summaries).
4. Evaluate across all four review dimensions:
   - **A. Writing quality and readability** — scientific prose, no code identifiers, self-contained captions
   - **B. Logical flow and structure** — section transitions, argument order, abstract-conclusion consistency
   - **C. Evidence-claim mapping** — every claim backed by figure/table/artifact, axes labeled with units
   - **D. Physics and methodology** — parameter values reasonable, approximations valid, convergence documented
5. Write `REVIEW_DIRECTIVE.md` with decision: APPROVE, REVISE, or NEEDS_REWORK.
6. If REVISE or NEEDS_REWORK: also write `FOLLOWUP_PROMPT.md` with ordered, actionable tasks for the next Opus iteration.
7. If APPROVE: signal the Execution Engineer to proceed to the polish phase.
