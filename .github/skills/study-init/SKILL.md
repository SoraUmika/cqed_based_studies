---
name: study-init
description: "Initialize a new cQED study. Use when: starting a new study, creating study folder, scaffolding simulation project, setting up study README. Creates the full folder structure and README with mandatory sections."
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

### 2. Create the Folder Structure

Create the following directory tree under `studies/`:

```
studies/<study_name>/
├── README.md
├── scripts/
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
- **Expected Outcomes** — quantitative success criteria where possible
- **Status** — set to `ACTIVE`

### 4. Confirm with the User

Display the generated README and ask whether the scope, goals, and methods look correct before proceeding to implementation.

## Rules

- Every study **must** have a README with all mandatory sections before any code is written.
- The study name must be lowercase with underscores — no spaces, hyphens, or capitals.
- Do not create placeholder scripts or empty Python files. Scripts are written in Step 3 (Implement).
