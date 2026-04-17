---
name: reproducibility-notebook
description: "Create the mandatory end-to-end reproducibility notebook for a cQED study. Use when: creating scripts/reproducibility_notebook.ipynb, restructuring a study notebook, exposing tunable parameters, adding load-vs-rerun paths, or making a study reproducible for future agents."
argument-hint: "Path to study folder, e.g. studies/transmon_chi_shift_optimization"
---

# Build a Reproducibility Notebook

## When to Use

- Creating `scripts/reproducibility_notebook.ipynb`
- Repairing a notebook that lacks a clear parameter cell or rerun path
- Converting a one-off analysis notebook into a reusable study handoff

## Required Notebook Structure

### 1. Title and Overview

- State the study title, problem class, and headline result.
- Give one paragraph that explains what the notebook reproduces.

### 2. Environment Setup

- Import required packages.
- Set paths relative to the study folder.
- Import shared helpers such as `common.py` when available.

### 3. User-Tunable Parameters

- Create one clearly labeled markdown cell and one code cell.
- Expose every adjustable knob here: Hilbert-space dimensions, optimizer settings, noise parameters, cost weights, convergence settings, diagnostic toggles, and figure controls.
- Print a short configuration summary at the end of the code cell.

### 4. Derived Objects

- Build the model, subspace, target, noise spec, and any derived grids from the tunable parameters.
- Re-running Sections 3 and 4 must be sufficient to propagate a parameter change downstream.

### 5. One Section per Major Result

For each result:

- Add a markdown cell that explains what the step does, why it matters, and which tunable parameters affect it.
- Add a default code cell labeled `# --- Load saved results (default) ---`.
- Add a commented rerun cell labeled `# --- Re-run with current parameters ---`.
- Show expected outputs or plots so a user can verify success.

### 6. Validation

- Reproduce sanity checks and convergence checks.
- State what should agree numerically and what tolerance is acceptable.

### 7. Key Figures

- Re-generate or display the publication figures from saved results.

### 8. Summary

- End with a markdown summary of what was reproduced.
- Include a table mapping each tunable parameter to its default value and effect on the results.

## Acceptance Checklist

- Every code cell is preceded by a markdown explanation.
- There is exactly one early parameter cell where users can change all main knobs.
- There is a derived-objects cell immediately after the parameter cell.
- Every major result has both a fast load path and a rerun path.
- The notebook can be read top-to-bottom without opening other files first.
- The notebook favors loading saved artifacts over expensive recomputation.

## Common Failure Modes

- Parameters scattered across many later cells.
- Hard-coded paths or hidden magic numbers.
- Rerun cells that do not actually use the tunable parameters.
- Notebook ends without telling the user what outputs to expect.