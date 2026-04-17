---
name: study-validator
description: "Validate that a study folder is structurally complete before marking COMPLETE. Checks README sections, IMPROVEMENTS.md tags, artifacts, report sections, figures, and reproducibility notebook. Use when: finishing a study, before writing REVIEW_REQUEST.md, or auditing an existing study."
argument-hint: "Path to study folder, e.g. studies/transmon_chi_shift_optimization"
---

# Study Structure Validator

## When to Use

- Before marking a study as COMPLETE
- Before writing REVIEW_REQUEST.md
- When auditing an existing study for completeness
- After the report is compiled but before the review phase

## Validation Checklist

Run through every check below. Report each as PASS or FAIL with a brief explanation.

### 1. README.md Structure

Verify the study `README.md` contains all mandatory sections:

- [ ] `## Problem Class` — contains one or more of: OPT, REP, DES, ANA
- [ ] `## Motivation` — non-empty
- [ ] `## Goals` — numbered list, at least one goal
- [ ] `## Methods` — references cqed_sim modules/functions
- [ ] `## Analytic Preliminary` — present, starts from the first-principles model when feasible, and states controlled approximations (even if it explains why no useful analytic result exists)
- [ ] `## cqed_sim Gap Analysis` — present
- [ ] `## Assumptions` — non-empty
- [ ] `## Compute & Resource Strategy` — present
- [ ] `## Expected Outcomes` — non-empty
- [ ] `## Known Limitations` — present
- [ ] `## Validation` — contains three checkboxes (sanity, convergence, literature)
- [ ] `## Status` — one of ACTIVE, COMPLETE, BLOCKED

### 2. IMPROVEMENTS.md Structure

- [ ] File exists
- [ ] Contains `## Critical Gaps (P1)` section
- [ ] Contains `## Recommended Improvements (P2)` section
- [ ] Contains `## Nice-to-Haves (P3)` section
- [ ] Contains `## What Was Tried and Did Not Work` section
- [ ] Contains `## Compute & Resource Notes` section
- [ ] Every improvement item has a priority tag (P1/P2/P3) and difficulty tag (LOW/MEDIUM/HIGH)

### 3. Artifacts and Data

- [ ] `artifacts/` directory exists and is non-empty
- [ ] At least one machine-readable artifact (JSON, NPZ, or CSV) per headline result
- [ ] JSON artifacts include metadata fields: `study_name`, `date_created`, `description`
- [ ] `data/` directory exists

### 4. Figures

- [ ] `figures/` directory exists and is non-empty
- [ ] Every figure exists in both `.png` and `.pdf` formats
- [ ] Figure filenames are descriptive (not numbered like `fig1.png`)

### 5. Report

- [ ] `report/report.tex` exists
- [ ] `report/references.bib` exists
- [ ] `report/report.pdf` exists and is non-zero size
- [ ] Report contains Abstract section
- [ ] Report contains Introduction section
- [ ] Report contains System and Methods section
- [ ] Report contains Results section
- [ ] Report contains Validation section
- [ ] Report contains Discussion section
- [ ] Report contains Conclusion section
- [ ] Report contains Limitations and Future Work section
- [ ] Report contains appendix: Detailed Results and Data
- [ ] Report contains appendix: Reproducibility
- [ ] No filenames or script references in main text (before `\appendix`)
- [ ] No `snake_case` or `camelCase` identifiers in main text prose

### 6. Reproducibility Notebook

- [ ] `scripts/reproducibility_notebook.ipynb` exists
- [ ] Notebook contains a title and overview cell
- [ ] Notebook contains a user-tunable parameters cell
- [ ] Notebook contains load-saved-results cells
- [ ] Notebook contains figure reproduction cells

### 7. Validation Status

- [ ] README `## Validation` section has at least two checks marked `[x]`
- [ ] Sanity checks are documented with specific evidence
- [ ] Convergence analysis is documented with specific evidence

## Output Format

```markdown
## Study Validation Report: <study_name>

**Overall:** PASS / FAIL (N of M checks passed)

### Failures
| Check | Section | Issue |
|-------|---------|-------|
| ... | ... | ... |

### Warnings
- ...

### Recommendation
READY FOR REVIEW / NEEDS WORK (list specific items to fix)
```

## Action on Failure

If any check fails:
1. List all failures in a clear table.
2. Suggest the specific fix for each failure.
3. Do NOT mark the study as COMPLETE.
4. Do NOT write REVIEW_REQUEST.md until all failures are resolved.
