---
name: study-init
description: "Initialize a new cQED study. Use when: starting a new study, creating study folder, scaffolding simulation project, setting up study README, or framing an analytic-first study plan from first principles. Creates the full folder structure and README with mandatory sections."
argument-hint: "Descriptive study name, e.g. transmon_chi_shift_optimization"
---

# Initialize a New cQED Study

## When to Use

- Starting any new simulation, optimization, analysis, or reproduction task
- User says "new study", "set up a study", "initialize", or names a new research topic
- Beginning Step 1 of the AGENTS.md workflow

## Procedure

### 1. Determine the Study Name

Ask the user for a descriptive name if not provided. Must be `lowercase_with_underscores`.

Examples: `transmon_chi_shift_optimization`, `kerr_cat_state_prep`, `snail_readout_dispersive_shift`

### 2. Consult Cross-Study Lessons

Read `LESSONS_LEARNED.md` at the repo root. Check for:
- Failed approaches relevant to the new study's problem domain
- cqed_sim quirks that apply to the planned methods
- Parameter ranges that have been validated in prior studies
- Compute timing data to inform the resource estimate

Incorporate relevant lessons into the README sections (Assumptions, cqed_sim Gap Analysis, Compute & Resource Strategy).

### 3. Create the Folder Structure

Create the following directory tree under `studies/`:

```
studies/<study_name>/
├── README.md
├── IMPROVEMENTS.md
├── scripts/
│   └── reproducibility_notebook.ipynb  ← REQUIRED (created in Step 6)
├── data/
├── figures/
└── report/
```

### 3. Generate the README

Use the [README template](./assets/README_TEMPLATE.md) as the starting point. Fill in:

- **Study Title** — from the user's description
- **Problem Class** — classify as `OPT`, `REP`, `DES`, `ANA`, or a combination
- **Motivation** — why this study matters; include paper references for `REP` class
- **Goals** — numbered, concrete, falsifiable goals
- **Methods** — which `cqed_sim` modules and functions will be used (consult the API Reference)
- **Analytic Preliminary** — the first-principles model plus the closed-form or limiting-case reasoning attempted before numerics, including any controlled approximations and why they are valid
- **cqed_sim Gap Analysis** — what is needed, what exists upstream, and what must be extended or implemented locally
- **Assumptions** — physical assumptions, parameter ranges, and convergence criteria
- **Compute & Resource Strategy** — cost estimate, planned acceleration, and expected bottlenecks
- **Expected Outcomes** — quantitative success criteria where possible
- **Validation** — initialize the three gate checks as unchecked items
- **Status** — set to `ACTIVE`

If the study consolidates or extends earlier studies, add a `## Study Composition` section that maps inherited components to their original source and role.

### 4. Confirm with the User

Display the generated README and ask whether the scope, goals, and methods look correct before proceeding to implementation.

## Rules

- Every study **must** have a README with all mandatory sections before any code is written.
- The study name must be lowercase with underscores — no spaces, hyphens, or capitals.
- Do not create placeholder scripts or empty Python files. Scripts are written in Step 3 (Implement).
- The README must anchor the planning work: first-principles analytic reasoning, controlled approximations, cqed_sim gap analysis, assumptions, compute strategy, and validation status all belong there, not only in chat.
- Prefer a first-principles Hamiltonian or governing-equation view before simplifying the problem. Any approximation used to make the study tractable must be stated explicitly with its validity condition.
- Every study **must** include `scripts/reproducibility_notebook.ipynb` before it can be marked COMPLETE. The notebook is created in Step 6 (see AGENTS.md).
- The reproducibility notebook **must** expose all tunable parameters (Hilbert space dims, optimizer settings, noise model, cost weights, probe states, sweep ranges) in a single early cell so users can modify and re-run the workflow without rewriting code. See AGENTS.md Step 6 for the full specification.
