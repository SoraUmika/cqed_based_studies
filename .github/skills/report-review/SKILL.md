---
name: report-review
description: "Critically review a completed cQED study report. Use when: acting as the Critical Reviewer (Codex 5.4 xHigh), evaluating a report after Opus completes Stage 1, deciding APPROVE/REVISE/NEEDS_REWORK, writing REVIEW_DIRECTIVE.md, and generating FOLLOWUP_PROMPT.md for the next Opus iteration."
argument-hint: "study=studies/<name> run=task_runs/<slug>"
---

# Critical Review of a cQED Study Report

## When to Use

- After the Execution Engineer (Opus 4.6) completes a study and writes `REVIEW_REQUEST.md`
- When acting as the Critical Reviewer (`@science-director` agent, Codex 5.4 xHigh)
- Before deciding whether a study is ready for final polish or needs another iteration

## Prerequisites

Before beginning the review, confirm:

1. `REVIEW_REQUEST.md` exists in the run directory — the Execution Engineer has signalled readiness
2. `report/report.pdf` or `report/report.tex` exists in the study directory
3. `EXECUTION_SUMMARY.md` exists in the run directory
4. Key figures are present in `studies/<name>/figures/`

Do not skip reading the full report. The execution summary is self-reported and may omit weaknesses.

---

## Review Procedure

Work through all five checks in order. Document every issue with a specific location (section, paragraph, or figure number).

---

### Check 1 — Writing Quality and Readability

The goal: can a human researcher read this paper fluently and understand it without access to the code?

#### 1a. Abstract
- [ ] States the problem and motivation in the first 1–2 sentences
- [ ] Describes the method briefly
- [ ] States the main quantitative result (not just "we found X" — give the number)
- [ ] States the significance or conclusion
- [ ] All quantities mentioned in the abstract are defined in the abstract (no undefined symbols)
- [ ] Does not require reading the body to make sense

#### 1b. Introduction
- [ ] Provides context: what is the physical system and why does it matter?
- [ ] Surveys relevant prior work with citations
- [ ] States the specific gap this study fills
- [ ] States the study objective clearly
- [ ] Provides a roadmap for the rest of the paper ("In Section II, we...")

#### 1c. Notation and Definitions
- [ ] Every symbol is defined at first use
- [ ] Units are stated the first time each physical quantity appears
- [ ] Acronyms are spelled out at first use
- [ ] No undefined technical terms

**Red flags:**
- A symbol that appears in the results before being defined in the methods
- A frequency stated without specifying whether it is $f$ (Hz) or $\omega$ (rad/s) or $\omega/2\pi$ (Hz)
- "SNR", "RWA", "DRAG" used without prior definition

#### 1d. Prose Quality
- [ ] No `snake_case`, `camelCase`, or backtick-quoted identifiers in the main text
- [ ] No references to filenames, script names, or code paths in the main text
- [ ] Sentences are grammatically correct
- [ ] Each paragraph has a clear topic and does not mix multiple ideas
- [ ] Interpretation is present: results sections should not merely list numbers without explaining what they mean

#### 1e. Figures and Captions
- [ ] Every figure has a caption
- [ ] Every caption is self-contained: the reader can understand the figure from the caption alone
- [ ] Every axis has a label with units
- [ ] Every curve/dataset is identified (in legend or caption)
- [ ] For 2D colormaps: colorbar is labeled with units
- [ ] No figure is missing from the text (every figure is referenced with `\ref`)
- [ ] No figure is orphaned (no figure exists that is not referenced)

#### 1f. Structural Flow
- [ ] Each section follows logically from the previous
- [ ] Section transitions exist (not just abrupt section breaks)
- [ ] The Discussion interprets the results — it does not merely restate them
- [ ] The Conclusion is consistent with the Abstract and Introduction

**Scoring:**
Record a verdict for each subsection: `PASS`, `MINOR ISSUES`, or `FAILS`. Any `FAILS` requires a required fix. `MINOR ISSUES` should be included as suggestions.

---

### Check 2 — Evidence-Claim Mapping

The goal: every non-trivial claim is directly supported by evidence in the report.

#### Procedure

1. Read through the Results and Discussion sections.
2. For every claim that is not a direct restatement of a number from a figure, ask: "What is the evidence?"
3. Fill in the audit table:

| Claim (section, ¶) | Type | Cited evidence | Verdict |
|--------------------|------|---------------|---------|
| {exact claim} | Quantitative / Qualitative / Comparative | {Fig. X / Table Y / Artifact Z} | SUPPORTED / WEAK / UNSUPPORTED |

**Claim types and expected evidence:**

| Claim type | Required evidence |
|------------|-----------------|
| "Fidelity of X% was achieved" | Figure showing the fidelity value; convergence data confirming it is not a local artifact |
| "Method A outperforms method B" | Side-by-side comparison with the same conditions; statistical support if stochastic |
| "The approximation is valid" | Parameter check (e.g., g/Δ ≪ 1 shown explicitly); limiting-case comparison |
| "The result is near-optimal" | Literature comparison or systematic sweep showing no better solution was found |
| "Convergence was verified" | Numbers: what was varied, by how much, what changed in the observable |
| "Results agree with [paper]" | Quantitative comparison table with percent error |

**Verdicts:**
- **SUPPORTED** — the cited evidence directly and specifically supports the claim
- **WEAK** — some evidence exists but it is indirect, too qualitative, or insufficient for the strength of the claim
- **UNSUPPORTED** — no evidence is cited or the cited evidence does not support the claim

Any UNSUPPORTED claim is a required fix before APPROVE. WEAK claims must be either strengthened or qualified.

---

### Check 3 — Physics and Methodology

The goal: the physics is correct, the approximations are justified, and the numerical methods are sound.

#### 3a. Parameter Sanity

For the stated system (transmon, fluxonium, etc.), check:

| Parameter | Typical range | Reported value | Flag? |
|-----------|--------------|----------------|-------|
| Qubit frequency ωq/2π | 4–8 GHz | {value} | {yes/no} |
| Anharmonicity α/2π | −200 to −100 MHz | {value} | {yes/no} |
| Coupling g/2π | 50–200 MHz | {value} | {yes/no} |
| Cavity frequency ωr/2π | 5–10 GHz | {value} | {yes/no} |
| Dispersive shift χ/2π | 0.1–10 MHz | {value} | {yes/no} |
| T1 (qubit) | 10–1000 μs | {value} | {yes/no} |
| κ/2π (cavity) | 0.1–10 MHz | {value} | {yes/no} |

Parameters outside typical ranges must be explained with physical reasoning, not just asserted.

#### 3b. Approximation Validity

| Approximation | Validity condition | Verified in report? |
|--------------|-------------------|-------------------|
| Rotating wave approximation | Drive amplitude ≪ transition frequency | {yes/no} |
| Dispersive approximation | g/Δ ≪ 1 (Δ = qubit-cavity detuning) | {yes/no} |
| Two-level (qubit) approximation | Drive frequency far from higher transitions | {yes/no} |
| Markov (Lindblad) | System correlation time ≪ bath correlation time | {yes/no — often assumed} |

#### 3c. Convergence Quality

Convergence checks must report numbers. "Convergence was verified" is insufficient.

Acceptable: "Doubling N_storage from 20 to 40 changed the fidelity by 1.2×10⁻⁵, confirming convergence at the 10⁻⁴ level."

Unacceptable: "We verified convergence by doubling the Hilbert space dimension."

Check:
- [ ] Hilbert space dimension sweep reported with numbers
- [ ] Time step sweep reported with numbers (if time-domain simulation)
- [ ] Optimization convergence shown (cost function vs. iteration)
- [ ] Results stable across at least two levels of refinement

#### 3d. Sanity Check Quality

Sanity checks must report results. "Sanity checks passed" is insufficient.

Acceptable: "With drive amplitude set to zero, the system remained in the ground state throughout the simulation (fidelity = 1.000 ± 10⁻¹²)."

Unacceptable: "We verified zero-drive and unitarity."

#### 3e. Literature Comparison

- For REP-class studies: quantitative comparison table required
- For other classes: comparison to at least one published or analytically known limiting case
- Percent errors or relative differences must be stated explicitly

---

### Check 4 — Completeness

| Item | Present? | Notes |
|------|---------|-------|
| Reproducibility appendix | {yes/no} | Must include: parameters, artifacts, reproduction steps |
| All optimized parameter values tabulated | {yes/no} | Final result parameters, not just settings |
| Artifacts in `artifacts/` directory | {yes/no} | JSON/NPZ/CSV with metadata |
| `scripts/reproducibility_notebook.ipynb` | {yes/no} | Must run end-to-end |
| IMPROVEMENTS.md updated | {yes/no} | Must reflect actual study limitations |
| All questions from Introduction answered | {yes/no} | |
| Limitations section specific | {yes/no} | Generic boilerplate does not count |

---

### Check 5 — Reviewer Pre-Check (for revision iterations)

If this is not the first iteration, check the `EXECUTION_SUMMARY.md → Reviewer Pre-Check` table:

- Did the Execution Engineer acknowledge all required actions from the prior review?
- For each required action: is the stated evidence actually present in the updated report?

Do not accept "addressed" without verification. Check the specific sections, figures, and data files.

---

## Writing the Review Output

### `REVIEW_DIRECTIVE.md`

After completing all five checks, write a complete `REVIEW_DIRECTIVE.md` in the run directory. Use the template in `science-director.agent.md`.

**Required sections:**
1. Decision (APPROVE / REVISE / NEEDS_REWORK)
2. Summary verdict (2–4 sentences)
3. Writing quality assessment — strengths, required fixes, suggestions
4. Evidence-claim audit table (complete, covering all non-trivial claims)
5. Physics and methodology assessment — what is correct, what is problematic
6. Completeness check table
7. Required actions (if REVISE or NEEDS_REWORK) — ordered, specific, with success criteria
8. Open concerns (non-blocking)

### `FOLLOWUP_PROMPT.md` (if REVISE or NEEDS_REWORK)

Write a complete, self-contained prompt for the next Opus invocation. Use the template in `science-director.agent.md`.

**Required sections:**
1. Context (iteration number, study name, run path)
2. Prior report status (what was done, what was correct, what must be fixed)
3. Required actions — ordered by priority, each specific enough to execute
4. What to preserve (explicit list of correct content not to change)
5. Definition of acceptance for this iteration (criteria Codex will check)
6. Notes to the Execution Engineer

**Quality bar for FOLLOWUP_PROMPT.md:**
A well-written follow-up prompt gives the Execution Engineer everything needed to act without asking clarifying questions. Each required action specifies: what to do, where, expected output, and how to know it is complete.

---

## Decision Thresholds

| Decision | Use when |
|----------|---------|
| **APPROVE** | All checks pass; study is technically sound and clearly written; a human researcher could use this as a credible reference |
| **REVISE** | Core content is sound and physics is correct; targeted improvements needed (fix a figure, strengthen a claim, improve a section) without requiring new experiments |
| **NEEDS_REWORK** | Fundamental issues: wrong physics, missing critical controls, unsupported core conclusions, structural problems with the argument; Opus must revisit the experimental design |

**Bias toward honesty.** APPROVE means you would be comfortable with a human researcher citing this report. If that is not true, issue REVISE.

---

## Common Failure Patterns to Watch For

| Pattern | How to detect | Required action |
|---------|--------------|-----------------|
| Fidelity claimed without convergence | No convergence table with numbers | Add convergence sweep with numbers |
| "Near-optimal" without evidence | No literature comparison or systematic sweep | Add comparison or qualify the claim |
| Figure caption does not match figure | Read caption, then look at the figure | Rewrite caption or regenerate figure |
| Approximation invoked without verification | No parameter check (e.g., g/Δ value) | Add explicit validity check |
| Conclusion overstates the result | Conclusion uses stronger language than the Results data supports | Qualify the language |
| Results section lacks interpretation | Results describe what was computed but not what it means | Add interpretation sentences |
| Sanity check declared passed without numbers | "Sanity checks passed" with no details | Add specific results for each check |
| Abstract does not match body | Abstract mentions a result not in the body, or vice versa | Align abstract with body |
| Missing alternative explanations | Only one interpretation of data given when multiple are plausible | Add discussion of alternatives |
| IMPROVEMENTS.md not updated | Generic/empty limitations section | Update with study-specific limitations |
