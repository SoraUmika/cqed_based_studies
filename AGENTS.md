# cQED Autonomous Research Platform

Agents use the **cqed_sim** framework to simulate, optimize, analyze, and validate circuit quantum electrodynamics (cQED) problems — then produce publication-quality reports.

---

## Navigation Index

| § | Section | Quick-find |
|---|---------|-----------|
| 1 | [Non-Negotiable Rules](#1-non-negotiable-rules) | Hard rules, consequences, cross-references |
| 2 | [Problem Classes](#2-problem-classes) | OPT / REP / DES / ANA |
| 3 | [Key Paths](#3-key-paths) | Study folder layout, file locations |
| 4 | [Agent Invocation Guide](#4-agent-invocation-guide) | Per-platform commands, state files, VS Code tasks |
| 5 | [Canonical Research Workflow](#5-canonical-research-workflow) | Steps 1–7 (Opus) + Stages 2–4 (Codex review loop) |
| 6 | [Report Format](#6-report-format) | LaTeX template, section rules, prose standards, OPT/DES requirements |
| 7 | [cqed_sim Framework](#7-cqed_sim-framework) | Paths, usage policy, gap handling Path A/B |
| 8 | [Decision Trees](#8-decision-trees) | Standalone code? Install package? Study complete? |
| 9 | [Conventions](#9-conventions) | Naming, code quality, figures |
| 10 | [cQED System Parameters](#10-typical-cqed-system-parameters) | Default parameter table |

---

## §1 Non-Negotiable Rules

> Violating any rule is a **hard failure**. There are no exceptions.

| # | Rule | Consequence of Violation | See Also |
|---|------|--------------------------|---------|
| 1 | Always use `cqed_sim` — document any gap before writing standalone code | Results not reproducible; study invalid | §7, §8 |
| 2 | Never skip workflow steps: Initialize → Plan → Implement → Validate → Report | Incomplete studies; missing validation or documentation | §5 |
| 3 | Stop **immediately** on any `cqed_sim` inconsistency — log in `BLOCKERS.md`, report to user, do not continue | Continuing with potentially invalid results is a hard failure | §5.3, §7 |
| 4 | Install packages as needed (`pip install --user`); log every install in `IMPROVEMENTS.md § Compute & Resource Notes` | Undocumented environment drift | §5.2 |
| 5 | No virtual environments — use system Python 3.12.10 directly | Breaks reproducibility across the platform | — |
| 6 | Every study must include `scripts/reproducibility_notebook.ipynb` | Users and future agents cannot reproduce results interactively | §5.6 |
| 7 | Start from first principles — derive the minimal analytic model first; state every controlled approximation and its validity | Numerical work lacks physical grounding; conclusions may be misinterpreted | §5.2 |

---

## §2 Problem Classes

Classify every task into one or more classes **before** starting work.

| ID | Class | Description | Typical Deliverables |
|----|-------|-------------|---------------------|
| `OPT` | Parameter Optimization | Optimize control parameters for a target unitary or state transfer | Converged parameters, fidelity plots, landscape scans |
| `REP` | Result Reproduction | Reproduce published results — spectra, dynamics, benchmarks | Comparison plots, quantitative agreement metrics |
| `DES` | Experiment Design | Design state preparation, gate implementation, or measurement protocols | Pulse sequences, protocol specs, simulated outcomes |
| `ANA` | System Analysis | Extract physical insights, identify optimal operating points | Parameter sweeps, phase diagrams, operating-point recommendations |

---

## §3 Key Paths

| What | Where |
|------|-------|
| Study root | `studies/<study_name>/` |
| Scripts | `studies/<study_name>/scripts/` |
| Data | `studies/<study_name>/data/` |
| Figures | `studies/<study_name>/figures/` |
| Artifacts | `studies/<study_name>/artifacts/` |
| Report | `studies/<study_name>/report/report.tex` |
| Improvement log | `studies/<study_name>/IMPROVEMENTS.md` |
| Reproducibility notebook | `studies/<study_name>/scripts/reproducibility_notebook.ipynb` |
| API Reference | [API_REFERENCE.md](https://github.com/SoraUmika/qubox_cQEDsim/blob/main/API_REFERENCE.md) |

---

## §4 Agent Invocation Guide

### §4.1 Agent Roles

| Agent | Platform | Model | Responsibilities |
|-------|----------|-------|-----------------|
| **Execution Engineer** | Claude Code | Opus 4.6 | Full study execution: planning, code, simulation, figures, validation, report writing |
| **Critical Reviewer** | GitHub Copilot | Codex 5.4 xHigh | Report quality review, evidence-claim verification, follow-up prompt generation |
| **Research Loop** | Either | Combined | Single-agent click-and-research mode; switches between both roles automatically |

The Execution Engineer is **not** a copy-editor. The Critical Reviewer is **not** a rubber stamp.

### §4.2 Claude Code — Opus 4.6 (Execution Engineer)

| Phase | Invocation |
|-------|-----------|
| Plan (new study) | `@execution-engineer study=studies/<name> run=task_runs/<slug> phase=plan` |
| Implement / refine | `@execution-engineer study=studies/<name> run=task_runs/<slug> phase=implement` |
| Validate | `@execution-engineer study=studies/<name> run=task_runs/<slug> phase=validate` |
| Report | `@execution-engineer study=studies/<name> run=task_runs/<slug> phase=report` |
| Polish (post-APPROVE) | `@execution-engineer study=studies/<name> run=task_runs/<slug> phase=polish` |
| Click-and-research | `@research-loop study=studies/<name> goal='<goal>'` |
| Resume interrupted | `@research-loop study=studies/<name> resume` |

### §4.3 GitHub Copilot — Codex 5.4 xHigh (Science Director)

> **Model selection is critical.** Select **Codex 5.4 xHigh** explicitly — not the Copilot default.

| Phase | Invocation |
|-------|-----------|
| Review | Open GitHub Copilot Chat → Select **Codex 5.4 xHigh** → `@science-director study=studies/<name> run=task_runs/<slug> phase=review` |

Output is written to `task_runs/<slug>/REVIEW_DIRECTIVE.md`.

**Watcher automation:** Run `tools/auto_loop.ps1` to monitor state files and auto-generate prompts. The watcher copies each prompt to clipboard and displays model selection instructions. Paste each prompt into GitHub Copilot Chat selecting the model shown.

### §4.4 State Files

| File | Location | Purpose |
|------|----------|---------|
| `study_state.json` | `studies/<name>/` | Machine-readable study state (single source of truth) |
| `SCIENCE_DIRECTIVE.md` | `task_runs/<name>/` | Science Director → Execution Engineer instructions |
| `EXECUTION_SUMMARY.md` | `task_runs/<name>/` | Execution Engineer → Science Director results |
| `TASK_CHECKLIST.md` | `task_runs/<name>/` | Task completion tracking |
| `PROGRESS_LOG.md` | `task_runs/<name>/` | Append-only log of what happened |
| `BLOCKERS.md` | `task_runs/<name>/` | Active and resolved blockers |

### §4.5 VS Code Tasks

Available from **Terminal → Run Task**:
- **Research: New Study** — Initialize a new study with goal
- **Research: Study Status** — Show current loop state
- **Research: Resume Study** — Detect phase and continue
- **Research: Run Loop Action** — Pick any phase to run

---

## §5 Canonical Research Workflow

> This is the **single, authoritative** workflow. It replaces all other workflow descriptions.

```
[Opus 4.6 — Claude Code — Execution Engineer]
  Step 1: Initialize  →  Create study folder + README + IMPROVEMENTS.md
  Step 2: Plan        →  Write SCIENCE_DIRECTIVE.md; analytic model first; compute cost estimate
  Step 3: Implement   →  Write scripts, run simulations, save data/figures, update IMPROVEMENTS.md
  Step 4: Validate    →  Sanity checks ✓  Convergence ✓  Literature comparison ✓
  Step 5: Report      →  Write report.tex (mandatory appendices) → compile PDF → self-review
  Step 6: Notebook    →  Create scripts/reproducibility_notebook.ipynb
  Step 7: Handoff     →  Write REVIEW_REQUEST.md → signal reviewer

[Codex 5.4 xHigh — GitHub Copilot — Science Director]
  Stage 2: Review     →  Critical review → REVIEW_DIRECTIVE.md (APPROVE / REVISE / NEEDS_REWORK)
  Stage 3: Refine     →  If REVISE/NEEDS_REWORK: write FOLLOWUP_PROMPT.md → Opus repeats Step 3–7
  Stage 4: Polish     →  On APPROVE: signal Opus for final polish → POLISH_COMPLETE.md → COMPLETE
```

### DO / DON'T

**DO:**
- ✅ Start from a first-principles analytic model whenever feasible, then introduce the minimal controlled approximations needed before running numerics — record the argument and approximation-validity conditions in the README
- ✅ Consult the API Reference before writing any simulation code
- ✅ Classify the problem (OPT/REP/DES/ANA) before starting
- ✅ Update `IMPROVEMENTS.md` in real time during implementation
- ✅ Put findings and discussion in the **main text**; put raw data, optimal parameters, and detailed plots in **appendices**
- ✅ Validate results before writing the report (all 3 checks)
- ✅ Tag every suggested improvement with priority (P1/P2/P3) and difficulty (LOW/MEDIUM/HIGH)
- ✅ Save figures in both `.png` (300 dpi) and `.pdf` (vector) formats
- ✅ Verify the compiled PDF has no overlapping text or layout errors; reformat overflowing equations using `multline`, `split`, or `align`
- ✅ Before running any long simulation, estimate its cost; if it may take more than a few minutes, apply GPU backends (CuPy, JAX, QuTiP with GPU) or parallelization (multiprocessing, joblib) before starting the run
- ✅ Provide a reproducibility notebook (`scripts/reproducibility_notebook.ipynb`) that reproduces the main results end-to-end
- ✅ Expose all tunable parameters in a single early cell so users can modify and re-run without rewriting the workflow
- ✅ State all mathematical equations explicitly — objective function, variables, constraints, physical model, and waveform parameterization
- ✅ For OPT/DES studies: report actual optimized parameter values (not just final fidelity); include time-domain and frequency-domain waveform visualizations
- ✅ For OPT/DES studies: test multiple pulse durations and report the duration–fidelity tradeoff

**DON'T:**
- ❌ Write ad-hoc simulation code when `cqed_sim` already has the functionality
- ❌ Skip the appendix — it is **required**, not optional
- ❌ Put raw pulse shapes, full parameter tables, or sweep data dumps in the main Results section
- ❌ Reference filenames, script names, or code identifiers in the main text — these belong in the appendix
- ❌ Use code-style identifiers with underscores or camelCase in main text prose
- ❌ Submit `report.tex` without a self-review pass
- ❌ Delete entries from `IMPROVEMENTS.md` — move resolved items to `## Resolved`
- ❌ Write vague limitation entries ("optimization didn't converge" → say **why** and **what to try next**)
- ❌ Leave layout errors in the compiled PDF
- ❌ Continue a run after detecting any inconsistency in `cqed_sim` — stop, log in `BLOCKERS.md`, report to user

---

### §5.1 Step 1 — Initialize Study

Create the study folder:

```
studies/<descriptive_name>/
├── README.md
├── IMPROVEMENTS.md
├── scripts/
│   ├── *.py
│   └── reproducibility_notebook.ipynb    ← REQUIRED
├── data/
├── figures/
├── artifacts/
└── report/
    ├── report.tex
    ├── references.bib
    └── report.pdf
```

**README.md must contain these sections:**

```markdown
# <Study Title>

## Problem Class
<!-- OPT | REP | DES | ANA -->

## Motivation
<!-- Why this study matters. Link to paper if REP class. -->

## Goals
<!-- Numbered, concrete, falsifiable goals. -->

## Methods
<!-- Which cqed_sim modules/functions will be used. -->

## Analytic Preliminary
<!-- First-principles model, closed-form derivation, or limiting-case reasoning attempted before numerics.
  List every controlled approximation and why it is valid.
  If no useful analytic result exists, explain why. -->

## cqed_sim Gap Analysis
<!-- Table format:
| Functionality | Needed? | Available in cqed_sim? | Plan |
|---|---|---|---|
-->

## Assumptions
<!-- Physical assumptions, parameter ranges, convergence criteria. -->

## Compute & Resource Strategy
<!-- Upfront cost estimate, planned acceleration, expected bottlenecks.
  Update with realized wall-clock times. -->

## Expected Outcomes
<!-- What success looks like — quantitative where possible. -->

## Known Limitations
<!-- Updated throughout the study. Approximations made, compute/framework constraints.
     Feeds directly into the report and IMPROVEMENTS.md. -->

## Validation
<!-- Keep current throughout the study:
- [ ] Sanity checks
- [ ] Convergence
- [ ] Literature comparison (if applicable)
-->

## Status
<!-- ACTIVE | COMPLETE | BLOCKED -->
```

If the study consolidates earlier work, add `## Study Composition` mapping each inherited component to its original study, key result, and role.

**IMPROVEMENTS.md template:**

```markdown
# Improvement Log: <Study Title>

> Written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
<!-- Things that could make results qualitatively wrong. -->
<!-- Format: - **<What>**: <Why it matters>. <What to do>. -->

## Recommended Improvements (P2)
<!-- Things that would meaningfully improve accuracy or scope. -->

## Nice-to-Haves (P3)

## Open Questions
<!-- Unresolved physics or numerical observations worth investigating. -->

## What Was Tried and Did Not Work
<!-- Failed approaches, dead-end parameter ranges, diverged algorithms.
     Include enough detail to understand WHY it failed. -->

## Compute & Resource Notes
<!-- Wall-clock times for key simulations. Memory usage. Bottlenecks. -->
```

**Rules for IMPROVEMENTS.md:**
1. **Start it in Step 1**, even with placeholder headings.
2. **Update in real time** during Steps 3–4. Log limitations immediately, not at the report phase.
3. **Never delete entries.** Move resolved items to `## Resolved` with a note on how they were fixed.
4. **Be specific about failures.** "Optimization didn't converge" is useless. State what was tried, at what fidelity it stalled, and what to try next.

---

### §5.2 Step 2 — Plan & Validate Approach

Before writing any simulation code:

1. **Attempt an analytic answer first.** Write down the minimal first-principles model (Hamiltonian, equations of motion, symmetry argument, conservation law) and determine whether the central question can be answered in closed form in a simplified or limiting case. Introduce only controlled approximations justified in the stated regime; record the approximation and its validity condition. Record in the README under `## Analytic Preliminary`. Only proceed to numerics after establishing the analytic picture.

2. **Estimate compute cost and plan acceleration.** Before writing simulation code, estimate whether any part (parameter sweeps, GRAPE, time evolution over many states, large Hilbert spaces) will exceed a few minutes of wall-clock time. If so, decide upfront:
   - **GPU backend:** prefer CuPy-backed or JAX-backed solvers; check whether `cqed_sim` or QuTiP supports GPU for the required operation.
   - **Parallelization:** use `multiprocessing.Pool`, `joblib.Parallel`, or similar for embarrassingly parallel workloads.
   - **Batching:** group independent runs to amortize Python overhead.
   Record the chosen strategy in the README under `## Compute & Resource Strategy` and log expected vs. actual wall-clock time in `IMPROVEMENTS.md § Compute & Resource Notes`.

3. **Check the [API Reference](https://github.com/SoraUmika/qubox_cQEDsim/blob/main/API_REFERENCE.md)** — confirm required functionality exists in `cqed_sim`.

4. **Identify gaps** — if `cqed_sim` cannot handle the task, document the gap in the README under `## cqed_sim Gap Analysis` and explain why standalone code is needed.

5. **State assumptions** — list all physical assumptions, parameter ranges, and convergence criteria in the README under `## Assumptions`.

> If new reusable functionality is developed, add `## Suggested Upstreaming` to the README.

---

### §5.3 Step 3 — Implement & Execute

- Write scripts in `scripts/`, save data to `data/`, generate figures to `figures/`.
- Update `IMPROVEMENTS.md` in real time — log limitations, failed approaches, and compute notes as they arise.
- Apply GPU/parallelization **before** starting any slow run:
  - Wrap parameter sweeps: `joblib.Parallel(n_jobs=-1)(delayed(fn)(p) for p in params)`
  - Pass `options={"gpu": True}` or equivalent when the solver supports it
  - Use numpy vectorization instead of Python loops wherever possible
  - Log actual wall-clock time in `IMPROVEMENTS.md § Compute & Resource Notes`
- **Stop immediately on any `cqed_sim` inconsistency** — unit/scale mismatch, sign/phase convention conflict, unexpected output dimensions, result contradicting a known analytic limit, or ambiguous API conventions. Log in `BLOCKERS.md`: (a) what was observed, (b) what was expected, (c) the specific API call or line, (d) suspected cause. Report to user and wait for resolution. See Rule 3.

---

### §5.4 Step 4 — Validate Results

**All three checks must pass before reporting:**

- [ ] **Sanity checks** — Verify limiting cases, conservation laws, or known analytic results.
- [ ] **Convergence** — Confirm results are stable with respect to Hilbert space truncation, time steps, and optimization iterations. Report what was varied and by how much.
- [ ] **Literature comparison** (if applicable) — Quantitatively compare to published benchmarks; report percent error or fidelity.

After validation, **finalize `IMPROVEMENTS.md`**: review every limitation, ensure suggested improvements have priority and difficulty tags, record open questions. Update the README `## Validation` section to reflect actual status.

---

### §5.5 Step 5 — Report

Write `report/report.tex` following the Scientific Review Paper Format in §6. The report **must** include:
1. A `Limitations and Future Work` section drawn from `IMPROVEMENTS.md`.
2. A **mandatory Appendix** with detailed results and data; the main text presents findings, the appendix presents supporting data.

Then compile to PDF and **perform a mandatory self-review** (read the full document as a researcher encountering it for the first time):

- [ ] Every sentence in the main text reads as standard scientific writing (Physical Review style).
- [ ] No filenames or script names appear in the main text.
- [ ] No code-style identifiers appear in the prose.
- [ ] Every equation is numbered and referenced in the text.
- [ ] Every figure and table is referenced in the text.
- [ ] The document compiles without layout errors.
- [ ] The README `## Validation` section is fully updated.
- [ ] Every major claim in Results/Discussion points to a figure, table, or saved artifact.
- [ ] `Saved Artifacts` subsection inventories every machine-readable output.

---

### §5.6 Step 6 — Reproducibility Notebook

Create `scripts/reproducibility_notebook.ipynb`. **Mandatory for every study.**

**Notebook structure:**

| Section | Content |
|---------|---------|
| 1. Title & Overview | Study name, goals, one-paragraph summary |
| 2. Environment Setup | Imports, path configuration, shared-module imports |
| 3. User-Tunable Parameters | **Single labeled cell** exposing every adjustable knob (Hilbert dims, optimizer settings, noise model, cost weights, probe states, convergence sweep ranges). Summary print at end. |
| 4. Derived Objects | Builds all simulation objects from Section 3 parameters. Re-running Sections 3→4 propagates any parameter change downstream. |
| 5. Step-by-step Reproduction | One section per major result: markdown cell (what/why/how/which parameters affect it) + "Load saved results" cell (default, fast) + commented-out "Re-run with current parameters" cell. |
| 6. Validation | Reproduces sanity checks and convergence analysis. |
| 7. Key Figures | Re-generates publication figures from the report. |
| 8. Summary | Final markdown cell summarizing what was reproduced and caveats. Table mapping each tunable parameter to its default value and effect on results. |

**Rules:**
- All tunable parameters in one place (Section 3). No magic numbers buried in later cells.
- Derived Objects cell (Section 4) constructs every simulation object from Section 3; re-running Sections 3→4 is sufficient to propagate parameter changes.
- Each reproduction step has (a) a default load-saved-data cell and (b) a commented-out re-run cell, clearly marked with `# --- Load saved results (default) ---` and `# --- Re-run with current parameters ---`.
- Load saved artifacts by default; include commented-out "full run" cells for expensive computations.
- Every code cell preceded by a markdown cell explaining what it does, why, and which parameters affect it.
- Self-contained: a user should be able to run top-to-bottom without reading other files first.
- Include expected output values or plots so the user can verify correctness.

Then update README status to `COMPLETE`.

> Do not mark `COMPLETE` unless the README `## Validation` section is checked, report appendices are present, and at least one machine-readable artifact exists for each headline result.

---

### §5.7 Step 7 — Handoff

Before signaling the reviewer, the following files must exist:

```
task_runs/<study>/
├── SCIENCE_DIRECTIVE.md       ← self-generated research plan
├── EXECUTION_SUMMARY.md       ← quantitative findings digest
├── REVIEW_REQUEST.md          ← signals readiness for review
├── TASK_CHECKLIST.md
└── PROGRESS_LOG.md
studies/<study>/
├── report/report.tex + report.pdf
├── figures/*.{png,pdf}
├── data/
└── artifacts/
```

**Quality bar before signaling review:**
- Every claim in the report is backed by a figure, table, or saved artifact.
- All validation checks are marked `[x]` in the README.
- The compiled PDF is free of layout errors.
- The reproducibility notebook exists and runs.

---

### §5.8 Stage 2 — Independent Critical Review (Codex 5.4 xHigh)

Codex reviews as a **referee for a high-impact physics journal** (Physical Review Letters / Nature Physics standard). **All five dimensions must pass. Default stance is REVISE.**

**APPROVE requires:** The reviewer can honestly say they would recommend this work for publication in a top-tier journal without major revisions. If any doubt exists, the correct decision is REVISE.

**Hard-reject triggers — any one → NEEDS_REWORK:**
- A central claim in the abstract or conclusion has no supporting figure, table, or artifact.
- An approximation is applied outside its stated validity regime without acknowledgment.
- An optimization is called "optimal" without multiple random restarts or a global-optimum argument.
- Convergence is declared but not demonstrated with quantitative data.
- The study adds no new insight beyond what is already known from prior work.

**A. Writing quality and readability**
- Is the abstract self-contained? Does it state the problem, method, main result, and conclusion without requiring the reader to have read the paper?
- Does the introduction survey relevant prior work and state the specific gap this study fills?
- Is every symbol defined at first use, with units, before it is used in an equation?
- Are code-style identifiers (`snake_case`, backtick symbols) present in prose? (Must be removed.)
- Can each figure caption be understood without reading the body?

**B. Evidence-claim mapping**
- For every non-trivial claim in Results, Discussion, and Conclusion: identify the claim, identify the supporting evidence (figure/table/artifact), and evaluate whether the evidence actually supports the claim.
- Common failure modes: quantitative claims backed only by qualitative figures; "significantly improved" with no baseline; convergence claimed but no convergence plot; "X works well" with no failure cases; approximation validity asserted without verification.
- Build an explicit audit table (required in the output).

**C. Physics and methodology**
- Are approximation validity bounds stated quantitatively (not just declared)?
- Are convergence results reported with numbers — what was varied, by how much, resulting change in key observable?
- Are uncertainty estimates or error bars provided for all key results? For OPT/DES: were multiple random restarts performed? Is there evidence or argument that the result is near-global?
- Has parameter sensitivity (±10% or ±1σ) been assessed and reported?
- Has the failure regime been characterized — the conditions under which the approach degrades?
- Is the literature comparison quantitative (percent errors, not "consistent with")?
- Are alternative physical mechanisms acknowledged?

**D. Completeness**
- Is the reproducibility appendix complete (optimized parameters, waveforms, assumptions, step-by-step procedure, saved artifacts)?
- Are all claims in the abstract supported in the body?
- Does the reproducibility notebook exist and appear correct?

**E. Novelty and scientific significance**
- What new insight does this study provide that was not derivable from prior work?
- Are the key metrics competitive with published state-of-the-art?
- Is the contribution clearly distinguished from prior work in the Introduction and Conclusion?
- Is the scope accurately stated (system-specific vs. general)?

**Review output — `REVIEW_DIRECTIVE.md`:**

```markdown
# Review Directive — Iteration {N}
Date: {ISO date}

## Decision
{APPROVE | REVISE | NEEDS_REWORK}

## Journal Review Score
| Dimension | Score (1–5) | Blocking issue? |
|-----------|------------|----------------|
| E. Novelty & scientific significance | | {Yes/No} |
| C. Technical soundness & methodology | | {Yes/No} |
| A. Clarity & presentation | | {Yes/No} |
| D. Reproducibility & completeness | | {Yes/No} |

Scores: 5=publication-ready, 4=minor fixes, 3=significant improvement needed, 2=major rework, 1=fundamental flaw
Equivalent journal verdict: {Accept | Minor Revision | Major Revision | Reject}

## Summary Verdict
{3–5 sentences: what the study accomplishes, what new insight it provides, whether it is convincing,
and what the primary gap is. Be specific — name the main result and the main weakness.}

## A. Writing Quality Assessment
### Strengths
- {specific strengths with location}
### Required Fixes (blocking)
- {section/paragraph}: {specific issue}
### Suggestions (non-blocking)
- {list}

## B. Evidence-Claim Audit
| Claim (exact text, section) | Supporting evidence | Verdict | Action |
|-----------------------------|---------------------|---------|--------|
| {exact claim} | {figure/table/artifact or "none"} | SUPPORTED / WEAK / UNSUPPORTED | {none / qualify / add evidence} |

## C. Physics and Methodology Assessment
### What is correct
- {list with physics reasoning}
### Required Fixes (blocking)
- {issue with location}: {suggested remedy with enough detail to act on}

### Convergence and Uncertainty Audit
- Hilbert space convergence: {what was reported / what is missing}
- Optimizer convergence: {what was reported / what is missing}
- Uncertainty/error bars: {present / absent — if absent, is an argument made?}
- Multiple restarts / global optimum evidence (OPT/DES): {yes / no / N/A}
- Parameter sensitivity: {reported / not reported / N/A}
- Approximation validity bounds: {stated quantitatively / asserted without verification}

## D. Completeness Check
- Reproducibility appendix: {Complete / Incomplete — list what is missing}
- Saved artifacts in artifacts/: {Present and documented / Missing}
- IMPROVEMENTS.md: {Current / Outdated — what's missing}
- Notebook runs end-to-end: {Verified / Not verified}
- All abstract claims supported in body: {Yes / No — list unsupported claims}

## E. Novelty and Significance Assessment
- New insight delivered: {state specifically what is new}
- Competitive with state-of-the-art: {yes / no / not applicable — cite comparison}
- Contribution delineated from prior work: {yes / no}
- Missing prior work: {list or "none"}

## Required Actions for Next Iteration
{Ordered by priority. Each item must be specific enough to execute without physics judgment.}
1. **[ACTION_TYPE]** {specific task}
   - What: {exact change needed}
   - Where: {section/script/figure}
   - Success criterion: {how to verify completion}

## Open Concerns (non-blocking)
{Issues for the record that do not block approval in this iteration.}
```

---

### §5.9 Stage 3 — Iterative Refinement

If Codex issues **REVISE** or **NEEDS_REWORK**, it writes `FOLLOWUP_PROMPT.md` — a complete, self-contained prompt ready to paste into a new Opus invocation. Opus reads the prior `EXECUTION_SUMMARY.md` and `REVIEW_DIRECTIVE.md`, then executes another full Steps 3–7 cycle. The loop continues until Codex issues **APPROVE**.

**Decision thresholds:**

| Decision | Meaning | Opus action |
|----------|---------|-------------|
| **APPROVE** | Study is technically convincing and well written | Proceed to Stage 4 final polish |
| **REVISE** | Core content is sound; specific targeted improvements needed | Execute follow-up tasks; extend report |
| **NEEDS_REWORK** | Fundamental issues: wrong physics, missing controls, unsupported conclusions | Revisit approach; may require redesigning experiments |

**Iteration limit:** If `research_config.json → loop.max_iterations` is reached before APPROVE, Codex must issue APPROVE or document why the study cannot be completed; it must not continue requesting more iterations.

**Follow-up prompt format (`FOLLOWUP_PROMPT.md`):**

```markdown
# Follow-Up Research Prompt — Iteration {N}
Generated by: Codex 5.4 xHigh Science Director
Study: studies/<name>
Run: task_runs/<slug>

## Context
This is iteration {N} of the research loop. The previous report was reviewed and found to have
the following deficiencies. This prompt directs the next Opus invocation to address them.

## Prior Report Status
{One paragraph summarizing what was achieved in the prior iteration and what remains unsatisfactory.}

## Required Actions (ordered by priority)
1. {Specific, actionable task — what to do, where, expected output, success criterion}
2. {Next task}
...

## What Must NOT Change
{List any correct results or well-written sections that must be preserved.}

## Definition of Acceptance for This Iteration
{Explicit criteria: what Codex will check in the next review. Opus should verify these before
signalling review-ready.}
```

---

### §5.10 Stage 4 — Final Polishing (Opus 4.6)

Triggered only after Codex issues **APPROVE** on technical content. Opus performs a dedicated final pass focused exclusively on presentation quality:

- Read the full report as a human researcher would.
- Improve sentence-level clarity, precision, and flow without changing scientific content.
- Ensure section-to-section transitions are logical.
- Verify abstract, introduction, and conclusion are mutually consistent and accurate.
- Check that all figure captions are self-contained and interpretable without reading the body.
- Remove any remaining code-style identifiers from prose.
- Verify the reference list is complete and consistently formatted.
- Confirm the Limitations and Future Work section is honest and specific.

Opus writes the revised `report.tex` (back up the prior version first) and compiles the final PDF. Then writes `POLISH_COMPLETE.md` in the run directory.

---

### §5.11 Quality Standards (Non-Negotiable)

Applied to every study. The reviewer enforces all eight. The bar is **publication in a high-impact physics journal**.

1. **Every important claim must be backed by evidence.** A claim without a supporting figure, table, or artifact is unsupported. Unsupported claims must be flagged as UNSUPPORTED in the evidence audit and resolved before APPROVE.
2. **Plots must be directly interpretable.** Every axis must have a label and unit. Every curve must be identified. A reviewer unfamiliar with the code must be able to understand the figure from the caption alone.
3. **Missing definitions are not acceptable.** Every symbol, approximation, and technical term must be defined at first use with units. Undefined notation is a required-fix issue.
4. **Convergence must be demonstrated, not declared.** "Convergence was verified" without numbers is insufficient and is a hard-reject trigger. State what was varied, by how much, and what change in the key observable was observed. Include a convergence plot or table.
5. **Uncertainty must be quantified.** Key numerical results must carry error bars, confidence intervals, or an explicit argument that uncertainty is negligible. For optimization results, multiple random restarts or a global-optimum argument is required before claiming "optimal." Sensitivity to ±10% parameter variation must be reported or explicitly waived with justification.
6. **Every study must demonstrate new scientific insight.** A technically correct study that adds nothing new to the field does not merit APPROVE. The Introduction and Conclusion must explicitly state what is new.
7. **Weak or incomplete studies must not be treated as final.** If evidence does not support the conclusion, the reviewer must issue REVISE or NEEDS_REWORK. Sunk cost (number of iterations) is never a reason to approve weak work.
8. **APPROVE means the reviewer would stake their scientific reputation on it.** Specifically: the reviewer must be able to say they would recommend this work for publication in a top-tier journal (Physical Review Letters / npj Quantum Information equivalent) without major revisions.

---

## §6 Report Format

### §6.1 Required LaTeX Template

```latex
\documentclass[aps,pra,twocolumn,reprint,amsmath,amssymb]{revtex4-2}

% ── Packages ──────────────────────────────────────────────
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{siunitx}
\usepackage[colorlinks=true,linkcolor=blue,citecolor=blue,urlcolor=blue]{hyperref}
\usepackage{float}
\usepackage{xcolor}

\begin{document}

% ── Metadata ──────────────────────────────────────────────
\title{<Study Title>}
\author{<Author(s)>}
\affiliation{<Affiliation>}
\date{\today}

% ══════════════════════════════════════════════════════════
% 1. ABSTRACT  (150–300 words, self-contained, no citations/equations)
% ══════════════════════════════════════════════════════════
\begin{abstract}
\end{abstract}

\maketitle

% ══════════════════════════════════════════════════════════
% 2. INTRODUCTION
% Context, motivation, prior work (citations), objectives, section roadmap.
% ══════════════════════════════════════════════════════════
\section{Introduction}

% ══════════════════════════════════════════════════════════
% 3. SYSTEM AND METHODS
% ══════════════════════════════════════════════════════════
\section{System and Methods}

\subsection{Hamiltonian}
% Full Hamiltonian in equation form. Define every symbol; give numerical value with units.

\subsection{Simulation Parameters}
% Tabulate all parameters with units (use siunitx):
% \begin{table*}[t]
%   \centering
%   \caption{System parameters.}
%   \begin{tabular}{@{} l S[table-format=+1.3] l @{}}
%     \toprule
%     Parameter & {Value} & Unit \\
%     \midrule
%     ...
%     \bottomrule
%   \end{tabular}
% \end{table*}

\subsection{Analytic Preliminary}
% First-principles model, perturbative/closed-form argument, controlled approximations
% with validity conditions. If no analytic result exists, state that explicitly.

\subsection{Computational Approach}
% cqed_sim modules used (with version), Hilbert space truncation, time-stepping,
% optimizer and convergence criteria (OPT), any standalone code with justification.

% ══════════════════════════════════════════════════════════
% 4. RESULTS
% Findings and interpretation. One subsection per study goal.
% Detailed raw data (pulses, parameter tables, sweeps) → APPENDIX.
% Every figure/table must be captioned, have labeled axes with units,
% and be referenced in the text.
% ══════════════════════════════════════════════════════════
\section{Results}

% ══════════════════════════════════════════════════════════
% 5. VALIDATION
% Report all three required checks with numbers, not declarations.
% ══════════════════════════════════════════════════════════
\section{Validation}
% \subsection{Sanity Checks}
% \subsection{Convergence Analysis}
% \subsection{Literature Comparison}  % if applicable

% ══════════════════════════════════════════════════════════
% 6. DISCUSSION
% Physical interpretation, quantitative comparison to literature,
% limitations, unexpected findings, open questions.
% ══════════════════════════════════════════════════════════
\section{Discussion}

% ══════════════════════════════════════════════════════════
% 7. CONCLUSION
% Concise summary. State whether each study goal was met.
% ══════════════════════════════════════════════════════════
\section{Conclusion}

% ══════════════════════════════════════════════════════════
% 8. LIMITATIONS AND FUTURE WORK  (MANDATORY — agent-to-agent handoff)
% Written for future agents. Be specific, honest, actionable.
% ══════════════════════════════════════════════════════════
\section{Limitations and Future Work}

\subsection{Known Limitations}
% For each limitation: (a) what it is, (b) why it exists, (c) how it affects reliability.
% Categories: optimization (local minima, single algorithm), truncation,
% physics omitted (no decoherence, RWA), compute budget, framework gaps.

\subsection{Suggested Improvements}
% For each: (a) what to do, (b) why it helps, (c) difficulty LOW/MEDIUM/HIGH, (d) priority P1/P2/P3.
% \begin{itemize}
%   \item \textbf{[P1 | HIGH]} Re-run with GRAPE to check for local minima.
%   \item \textbf{[P2 | LOW]} Increase Fock truncation from $N=8$ to $N=15$.
% \end{itemize}

\subsection{Open Questions}
% Physics questions or unexpected observations that emerged but were not resolved.

% ══════════════════════════════════════════════════════════
% 9. REFERENCES
% ══════════════════════════════════════════════════════════
\bibliographystyle{apsrev4-2}
\bibliography{references}

% ══════════════════════════════════════════════════════════
% APPENDICES (REQUIRED for ALL studies)
% Main body: findings, interpretation, discussion.
% Appendices: raw/detailed results that support, reproduce, or extend the work.
% ══════════════════════════════════════════════════════════
\appendix

\section{Detailed Results and Data}
% OPT: optimal pulse shapes (I/Q vs time), full parameter tables, convergence traces,
%      cost landscape cross-sections.
% REP: side-by-side comparison plots, full numerical comparison tables with % errors.
% DES: complete pulse sequence diagrams with timing, protocol parameter tables,
%      state/process tomography matrices.
% ANA: full parameter sweep heatmaps or line cuts, phase diagram data, fit params with uncertainties.

\section{Reproducibility}
% MANDATORY. Everything a future agent or user needs to reproduce results from artifacts alone.
%
% \subsection{Optimized Parameters}
% Full table of every optimized/fitted/tuned parameter: gate angles, pulse amplitudes,
% phases, durations, fitted coefficients. If swept, report value at optimum.
%
% \subsection{Waveform and Pulse Information}
% Per pulse: number of time slices, dt, amplitude bounds, quadrature labels,
% reference to artifact file (preferred) or waveform data in table. Include I/Q vs time plots.
%
% \subsection{Gate Sequence and Decomposition}
% Exact ordered gate sequence to realize the target unitary (type, parameters, duration per gate).
% For GRAPE: total duration and number of control channels. For multi-decomposition comparisons:
% include a comparison table.
%
% \subsection{Modeling and Simulation Assumptions}
% Hilbert space dimensions, solver tolerances, rotating-frame definitions,
% approximations (RWA, dispersive), decoherence model or lack thereof.
% Enough detail to configure an identical simulation.
%
% \subsection{Reproduction Procedure}
% Step-by-step: (1) which scripts to run and in what order,
% (2) expected output files and locations, (3) how to verify agreement.
%
% \subsection{Saved Artifacts}
% List every machine-readable artifact in artifacts/ and data/:
% filename, format (JSON/NPZ/CSV), contents, example load code.

% \section{Supplementary Derivations}  % Optional
% \section{Simulation Configuration}   % Optional

\end{document}
```

---

### §6.2 Report Section Rules

| Section | Required? | Key Constraint |
|---------|-----------|----------------|
| Abstract | **Yes** | 150–300 words; self-contained; no citations, equations, or undefined acronyms |
| Introduction | **Yes** | Establish context, cite prior work, state objectives; end with one-sentence roadmap |
| System and Methods | **Yes** | Full Hamiltonian with all terms defined; parameters table with units (`siunitx`); identify `cqed_sim` modules |
| Results | **Yes** | One subsection per study goal; all figures captioned with labeled axes and units; reference every figure and table |
| Validation | **Yes** | Report all three checks with numbers; include convergence plots where applicable |
| Discussion | **Yes** | Interpret results physically; compare quantitatively to literature; state limitations honestly |
| Conclusion | **Yes** | Summarize findings; state whether each goal was met |
| Limitations & Future Work | **Yes** | **Agent-to-agent handoff section.** Known Limitations (cause + impact) + Suggested Improvements (P1–P3 + difficulty) + Open Questions |
| References | **Yes** | BibTeX (`references.bib`); cite all referenced papers, `cqed_sim`, and foundational references; use `apsrev4-2` style |
| Appendices | **Yes** | **Required for all studies.** Detailed Results + Reproducibility subsections (see §6.1 template) |

---

### §6.3 Main Text Prose Standards

#### No Filenames or Script References

Filenames, script names, and data-file paths must not appear in the main text. They belong exclusively in the appendix (Reproducibility section) and in `IMPROVEMENTS.md`.

| Incorrect (main text) | Correct (main text) |
|-----------------------|---------------------|
| "Results are stored in `grape_200ns.npz`." | "Optimized pulse parameters are provided in Appendix~\ref{app:reproducibility}." |
| "Running `run_sweep.py` produces Fig.~3." | "The parameter sweep (Fig.~3) was computed over the range..." |
| "See `optimal_unitary.json` for full data." | "Full unitary data are archived as machine-readable artifacts (Appendix~\ref{app:reproducibility})." |

**Quick self-check:** search the pre-appendix portion of `report.tex` for `.json`, `.npz`, `.csv`, `.py`, `data/`, and `artifacts/`. Any hit not inside a code listing is a prose violation.

#### No Code-Style Identifiers

The main text must not contain `snake_case`, `camelCase`, or any programming-style identifier.

| Incorrect | Correct |
|-----------|---------|
| `` `reduced_unitary_direct` `` | "the reduced unitary on the logical subspace" |
| `` `chi_shift` `` | "the dispersive shift $\chi$" |
| `` `grape_fidelity_200ns` `` | "the gate fidelity at $T = 200\,\mathrm{ns}$" |
| `` `run_opt_loop` `` | "the optimization loop" |
| `` `TargetStateMapping` `` | "the target-state mapping objective" |
| `` `objective_value` `` | "the reported objective value" |

Code identifiers may appear **only** in code listings (`verbatim`/`lstlisting`) and in the appendix when referring to artifact files or function calls.

#### Self-Review Pass Is Mandatory

After finishing `report.tex`, re-read the entire document as a researcher encountering it for the first time:

- [ ] Every sentence in the main text reads as standard scientific writing (Physical Review style).
- [ ] No filenames or script names appear in the main text.
- [ ] No code-style identifiers appear in the prose.
- [ ] Every equation is numbered and referenced in the text.
- [ ] Every figure and table is referenced in the text.
- [ ] The document compiles without layout errors.

---

### §6.4 OPT/DES — Optimization Report Requirements

Applies to all `OPT` and `DES` studies, and any study involving numerical optimization. **Mandatory.**

| Requirement | What to Include |
|-------------|----------------|
| **State equations explicitly** | Objective function (e.g., $\mathcal{F} = |\mathrm{Tr}(U_\mathrm{target}^\dagger U)|^2/d^2$), variables being optimized, constraints (e.g., $|\Omega(t)| \le \Omega_\mathrm{max}$), full Hamiltonian/Lindblad with every symbol defined and its numerical value with units |
| **Show waveform parameterization** | Amplitude (piecewise-constant, Gaussian, DRAG, Fourier); phase (IQ convention $\Omega(t)=I(t)\cos\omega_d t - Q(t)\sin\omega_d t$); detuning/rotating-frame offset; basis functions; hardware transfer-function assumptions |
| **Report optimized parameters** | Full converged parameter set; fidelity/cost at the converged solution; artifact file (JSON/NPZ in `artifacts/`) with example load snippet |
| **Visualize the optimized object** | Waveform I/Q or amplitude/phase vs. time; cost landscape or parameter sweep if relevant; full detail in appendix; summary figure may appear in main text |
| **Show time and frequency domains** | Time domain: I/Q or amplitude/phase with labeled axes and units. Frequency domain: $|\tilde{\Omega}(f)|$ vs. frequency in MHz/GHz. Both panels in appendix. |
| **Treat duration as an optimization axis** | Test multiple pulse durations; save separate artifact per duration (e.g., `grape_optimal_200ns.json`); report duration–fidelity tradeoff including: summary figure/table, shortest duration achieving acceptable performance, commentary on whether there is a sharp threshold and its physical reason |

---

### §6.5 Main Text vs. Appendix — Content Split

| Content Type | Goes In | Example |
|--------------|---------|---------|
| Key result (fidelity, key parameter) | **Main text** (Results) | "The optimized gate achieves $\mathcal{F} = 0.995$." |
| Physical interpretation | **Main text** (Discussion) | "The fidelity plateau is caused by leakage to $|2\rangle$." |
| Optimal pulse shape (I/Q vs. time) | **Appendix** | Figure: `optimal_pulse_iq.pdf` |
| Full converged parameter table | **Appendix** | Table with all 12 optimized pulse coefficients |
| Summary parameter table (key values only) | **Main text** (Results) | Table with gate time, fidelity, leakage |
| Convergence trace | **Validation** or **Appendix** | If brief: Validation; if many traces: Appendix |
| Parameter sweep heatmap (full) | **Appendix** | 2D grid of $\chi$ vs. $\kappa$ showing fidelity |
| Selected sweep cross-sections | **Main text** (Results) | 1D line cuts at optimal operating point |
| Cost landscape | **Appendix** | Contour plots of cost function |
| Side-by-side literature comparison | **Main text** (Validation, summary) + **Appendix** (full tables) | — |
| Pulse sequence timing diagrams | **Appendix** | Detailed protocol spec |

**Appendix content by problem class:**

| Problem Class | Required Appendix Content |
|---------------|---------------------------|
| **OPT** | Optimal pulse shapes (I/Q plots), full parameter tables, convergence traces, cost landscape cross-sections |
| **REP** | Full numerical comparison tables (simulation vs. published), overlay plots, percent-error breakdowns |
| **DES** | Complete pulse sequence diagrams, protocol parameter tables, simulated measurement outcomes, state/process tomography data |
| **ANA** | Full sweep heatmaps, phase diagram data, extracted fit parameters with uncertainties, all line cuts |

> **Rule of thumb:** If a figure/table shows *what the answer is* (detailed data) → appendix. If it shows *what the answer means* (insight, comparison, trend) → main text. When in doubt, put a summary in the main text and full version in the appendix.

---

### §6.6 Reproducibility Requirements

Every study report **must** include a `\section{Reproducibility}` appendix. This is the critical bridge between a completed study and anyone who wants to verify, extend, or build upon the work.

| Subsection | Content |
|------------|---------|
| **Optimized Parameters** | Full table of every optimized/fitted/tuned parameter. If swept, report value at optimum. |
| **Waveform and Pulse Information** | Per pulse: time slices, `dt`, amplitude bounds, quadrature labels, artifact reference or waveform table; I/Q vs time plots for key pulses. |
| **Gate Sequence and Decomposition** | Exact ordered gate sequence (type, parameters, duration). For GRAPE: total duration, number of control channels. For multi-decomposition comparisons: include comparison table. |
| **Modeling and Simulation Assumptions** | Hilbert space dimensions, solver tolerances, rotating-frame definitions, approximations (RWA, dispersive), decoherence model or lack thereof. Enough detail to configure an identical simulation. |
| **Reproduction Procedure** | Step-by-step: (1) which scripts to run and in what order, (2) expected output files and locations, (3) how to verify agreement (fidelity threshold, figure match). |
| **Saved Artifacts** | List every machine-readable artifact in `artifacts/` and `data/`: filename, format (JSON/NPZ/CSV), contents, example load code. Future agents must be able to load the solution programmatically. |

**Artifact formats and structure:**

| Format | Use Case |
|--------|----------|
| **JSON** | Parameters, metadata, small numerical results, gate sequences |
| **NPZ** | Large arrays (waveforms, sweep grids, density matrices) |
| **CSV** | Tabular data that should be human-readable |

Each artifact file must include metadata fields: `study_name`, `date_created`, `description`, `parameters`, and `load_instructions`.

```
artifacts/
├── grape_optimal_200ns.json     ← GRAPE schedule + fidelity for 200 ns
├── grape_optimal_400ns.json     ← GRAPE schedule + fidelity for 400 ns
├── decomposition_best.json      ← Best gate sequence
└── target_unitary.npz           ← Target unitary matrix
```

---

### §6.7 Limitations and Future Work — Guidance

Write this section for future agents picking up where the current agent left off. Be specific, honest, and actionable.

**Known Limitations categories:**

| Category | Examples |
|----------|----------|
| Optimization | Local minima not ruled out; only one algorithm tested; landscape not fully explored |
| Truncation | Fock basis cutoff not converged; Hilbert space dimension too low for high-photon states |
| Physics omitted | No decoherence channels; no thermal photon population; RWA or dispersive approximation may break down at range edges |
| Compute budget | Parameter sweep was coarse; longer pulse durations not tested; wall-clock limit reached |
| Framework gaps | `cqed_sim` does not support feature X; workaround Y was used instead |

**Priority and difficulty tags:**

| Tag | Meaning |
|-----|---------|
| **P1** | Critical gap — results may be qualitatively wrong without this fix |
| **P2** | Meaningful improvement — results qualitatively correct but quantitatively limited |
| **P3** | Nice to have — would strengthen the study but not essential |
| **LOW** | Change a parameter or re-run a script |
| **MEDIUM** | Write new code or modify the simulation setup |
| **HIGH** | Requires new `cqed_sim` functionality, significant compute, or new physics |

**Open Questions** should capture surprising or unresolved observations — anomalous data points, unexpected parameter sensitivities, or phenomena that deserve their own study.

---

### §6.8 Citation and Reference Guidelines

- **Bibliography file:** Every study must include `references.bib` in `report/`.
- **Citation style:** Numbered references in order of appearance (`apsrev4-2`).
- **Minimum citations:** (a) original papers for any reproduced results, (b) the `cqed_sim` framework, (c) foundational references (transmon, Jaynes–Cummings model).
- **BibTeX format:** Use `@article` for journal papers, `@misc` for preprints/software. Always include `doi` when available.

```bibtex
@article{koch2007transmon,
  author  = {Koch, Jens and Yu, Terri M. and Gambetta, Jay and others},
  title   = {Charge-insensitive qubit design derived from the {Cooper} pair box},
  journal = {Physical Review A},
  volume  = {76},
  pages   = {042319},
  year    = {2007},
  doi     = {10.1103/PhysRevA.76.042319}
}

@misc{cqed_sim,
  author = {<Author(s)>},
  title  = {cqed\_sim: Circuit QED Simulation Framework},
  year   = {2025},
  url    = {https://github.com/SoraUmika/qubox_cQEDsim}
}
```

---

### §6.9 Figure Standards

- Save both `.png` (300 dpi, for quick inspection) and `.pdf` (vector, for LaTeX inclusion).
- Use `\includegraphics` with relative paths from the `report/` directory.
- Every figure must have: descriptive caption, labeled axes with units, legible font size (≥ 8 pt in print), and a colorblind-friendly palette.
- Prefer vector formats (`.pdf`) in the compiled report.

---

### §6.10 Document Layout and Typesetting

The compiled PDF must be free of layout errors. **This is a hard requirement.**

- After compiling, visually inspect for overlapping text, figures running into body text, or captions colliding with other elements. Fix before finalizing.
- Long equations must fit within column/page margins. Never allow overflow.

| Situation | Recommended environment |
|-----------|------------------------|
| Single long equation needing one line break | `multline` |
| Equation with multiple aligned steps | `align` or `align*` |
| Sub-expression too long for one line | `split` inside `equation` |
| Several short equations to align at `=` | `align` |

**Compilation log — warnings that must be fixed:**
- `Overfull \hbox` — text or math overflowing margin; **must fix** (break equation at binary operator, align with `&`)
- `Overfull \vbox` — content overflowing a page or float; investigate and fix
- `Underfull \hbox` with badness ≥ 10000 — severely loose spacing; consider rewording or adjusting hyphenation

Warnings about undefined references (after the first pass) are expected and do not require action beyond ensuring the `.bib` file is complete.

---

### §6.11 Compilation

```bash
cd studies/<study_name>/report/
pdflatex report.tex
bibtex report
pdflatex report.tex
pdflatex report.tex
```

Or equivalently: `latexmk -pdf report.tex`.

---

## §7 cqed_sim Framework

### §7.1 Framework Resources

| Resource | Location |
|----------|----------|
| Source code (GitHub) | [SoraUmika/qubox_cQEDsim](https://github.com/SoraUmika/qubox_cQEDsim) |
| API Reference | [API_REFERENCE.md](https://github.com/SoraUmika/qubox_cQEDsim/blob/main/API_REFERENCE.md) |
| Physics Conventions | [physics_conventions_report.tex](https://github.com/SoraUmika/qubox_cQEDsim/blob/main/physics_and_conventions/physics_conventions_report.tex) |
| Local copy | `C:\Users\dazzl\Box\Shyam Shankar Quantum Circuits Group\Users\Users_JianJun\cQED_simulation` (most up-to-date) |

**Usage policy:** Always prefer `cqed_sim`. Do not duplicate functionality that already exists. Document any deviation in the study README.

### §7.2 Handling cqed_sim Gaps

When `cqed_sim` does not support a required experiment, **never skip it** — apply one of two paths. **Both paths are mandatory** — the study cannot be marked COMPLETE if a required experiment was simply skipped.

#### Path A — Contribute Upstream (preferred when simple)

Use Path A when **all** of the following are true:
1. Self-contained (≤ ~100 lines, no new dependencies).
2. Well-established math (no novel algorithm).
3. Likely reused across studies.

Steps:
1. Locate the relevant module in `cQED_simulation/`.
2. Add the function following existing code style.
3. Write a minimal docstring + one-line usage example.
4. Add a trivial sanity-check assert.
5. Add `### Upstream Additions` to the study README: function name, module, rationale.
6. Proceed with the experiment using the new function.

#### Path B — Local Implementation (when upstream is complex)

Use Path B when the feature requires a research-grade algorithm, a new dependency, significant architectural change to `cqed_sim`, or more than ~1 hour to add safely.

Steps:
1. Create `scripts/local_<feature_name>.py`.
2. Use NumPy / SciPy / QuTiP only unless justified.
3. Add a module docstring:
   ```python
   """
   Local implementation of <feature>.
   Reason not in cqed_sim: <one sentence>.
   Upstreaming priority: HIGH / MEDIUM / LOW.
   """
   ```
4. Log in `BLOCKERS.md`: `cqed_sim gap: <feature>. Local impl: scripts/local_<feature>.py`
5. Log in `IMPROVEMENTS.md` with a priority tag.

---

## §8 Decision Trees

### Can I use standalone simulation code?

```
Is the required functionality in cqed_sim?
├─ YES → Use cqed_sim. Do NOT write standalone code.
├─ PARTIALLY → Extend using cqed_sim as the foundation.
│              Document what is missing in the README.
└─ NO → Write standalone code.
        Document the gap in the README.
        Add a "Suggested Upstreaming" section if reusable.
```

### Should I install a package?

```
Is the package already installed?
├─ YES → Use it.
└─ NO → Would it genuinely improve the analysis or study?
         ├─ NO → Do not install.
         └─ YES → pip install --user <package>
                  Log the install in IMPROVEMENTS.md (Compute & Resource Notes).
```

### Is the study complete?

```
Are all goals in the README met?
├─ NO → Continue work or document blockers.
└─ YES → Have results been validated (sanity, convergence, literature)?
          ├─ NO → Complete validation first.
          └─ YES → Is IMPROVEMENTS.md populated with limitations, suggestions, and failed approaches?
                    ├─ NO → Complete the improvement log first.
                    └─ YES → Are key results saved as artifacts (JSON/NPZ in artifacts/)?
                              ├─ NO → Save optimized parameters, waveforms, and gate sequences.
                              └─ YES → Write report.tex:
                                        ─ Main text: findings, interpretation, discussion
                                        ─ Appendix: detailed data (pulses, parameters, sweeps)
                                        ─ Reproducibility appendix (parameters, waveforms,
                                          gate sequences, assumptions, procedure, artifacts)
                                        ─ Limitations & Future Work section
                                       → Compile PDF
                                       → Is reproducibility notebook in scripts/?
                                          ├─ NO → Create scripts/reproducibility_notebook.ipynb
                                          └─ YES → README status = COMPLETE
```

---

## §9 Conventions

### Naming

- Study folders: `studies/<lowercase_descriptive_name>/` (e.g., `studies/transmon_chi_shift_optimization/`)
- Scripts: `snake_case.py`
- Figures: `<figure_description>.png` or `.pdf`

### Code Quality

- All scripts must be self-contained and runnable from the study folder.
- Include docstrings explaining the physical setup and parameters.
- Pin key simulation parameters (Hilbert space dims, time steps, etc.) as named constants at the top of each script.

### Figures

- Use matplotlib with consistent styling (labeled axes, units, legends).
- Save both `.png` (for README / quick inspection) and `.pdf` (for LaTeX report).
- Use colorblind-friendly palettes where possible (e.g., `tab10`, `colorblind` from `seaborn`).

---

## §10 Typical cQED System Parameters

Use these as default starting values unless the study specifies otherwise.

| Parameter | Symbol | Value | Unit |
|-----------|--------|-------|------|
| Dispersive shift | χ | −2.84 | MHz |
| Second-order dispersive shift | χ′ | −21 | kHz |
| Cavity self-Kerr | K | −28 | kHz |
| Qubit anharmonicity | α | −255 | MHz |
| Qubit frequency | ω_q | 6.150 | GHz |
| Cavity frequency | ω_c | 5.241 | GHz |
| Cavity decay rate | T1_cavity | 4 | kHz |
| Readout resonator frequency | ω_r | 8.597 | GHz |
| Readout resonator linewidth | κ_r | 1.8 | MHz |

