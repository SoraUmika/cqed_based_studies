---
name: cross-study-analysis
description: "Compare results, parameters, and findings across multiple completed cQED studies. Use when: looking for patterns across studies, comparing fidelities/parameters, synthesizing lessons learned, planning new studies based on prior work, or building a meta-analysis."
argument-hint: "Goal description, e.g. 'compare fidelity results across all OPT studies' or 'what chi values have been explored?'"
---

# Cross-Study Analysis

## When to Use

- Comparing results across multiple completed studies
- Finding parameter ranges that have been explored
- Identifying which approaches worked best for a given problem type
- Planning a new study that builds on prior work
- Synthesizing lessons learned across the research program

## Data Sources

For each study in `studies/`, read:

1. `study_state.json` — status, problem class, completed/failed/pending tasks
2. `README.md` — goals, methods, assumptions, validation status, known limitations
3. `IMPROVEMENTS.md` — what was tried, what failed, suggested improvements
4. `EXECUTION_SUMMARY.md` (in `task_runs/<study>/`) — quantitative findings digest

## Analysis Modes

### Mode 1: Parameter Survey

Collect and tabulate parameter values across studies:

```markdown
## Parameter Survey: <parameter_name>

| Study | Value | Unit | Context | Result |
|-------|-------|------|---------|--------|
| study_a | -2.84 | MHz | dispersive readout | F=0.995 |
| study_b | -3.10 | MHz | gate optimization | F=0.987 |
```

### Mode 2: Fidelity Comparison

Compare optimization results across studies:

```markdown
## Fidelity Comparison

| Study | Problem Class | Target | Best Fidelity | Method | Duration | Notes |
|-------|--------------|--------|---------------|--------|----------|-------|
| study_a | OPT | CZ gate | 0.9995 | GRAPE | 200 ns | converged |
| study_b | DES | SNAP | 0.987 | Nelder-Mead | 400 ns | local min suspected |
```

### Mode 3: Lessons Learned Synthesis

Aggregate insights from IMPROVEMENTS.md across studies:

```markdown
## Lessons Learned

### What Consistently Works
- ...

### Common Failure Modes
- ...

### Unresolved Questions (across studies)
- ...

### Suggested Follow-Up Studies
- ...
```

### Mode 4: Gap Analysis

Identify what has NOT been explored:

```markdown
## Research Gaps

| Topic | Explored? | Studies | Missing |
|-------|----------|---------|---------|
| Decoherence effects | Partial | study_a (T1 only) | T_phi, thermal photons |
| Multi-qubit gates | No | — | Entire topic |
```

## Output Format

Always produce:
1. A summary table with the key comparison metric
2. A narrative synthesis interpreting the patterns
3. Recommendations for the user (next study to run, parameters to explore, approaches to try)

## Procedure

1. List all study directories in `studies/`.
2. Read `study_state.json` from each — filter to relevant studies based on the user's query.
3. For each relevant study, read `README.md` and `IMPROVEMENTS.md`.
4. If `task_runs/<study>/EXECUTION_SUMMARY.md` exists, read it for quantitative data.
5. Compile the comparison tables.
6. Write the narrative synthesis.
7. Make recommendations.
