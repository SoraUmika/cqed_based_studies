---
name: validate-results
description: "Run the 3-check validation gate before reporting. Use when: finishing a simulation study, checking convergence, verifying sanity checks, comparing to literature, preparing to write the final report. Enforces the AGENTS.md validation checklist."
argument-hint: "Path to study folder, e.g. studies/transmon_chi_shift_optimization"
---

# Validate Simulation Results

## When to Use

- After completing simulations (AGENTS.md Step 4)
- Before writing the final report
- When the user asks to "validate", "check convergence", or "compare to literature"

## Procedure

All three checks must pass before reporting. Work through them in order.

### Check 1 — Sanity Checks

Verify the simulation against known limiting cases, conservation laws, or analytic results.

**Common sanity checks for cQED simulations:**

| Check | How |
|-------|-----|
| Zero drive → no state change | Run simulation with all drives set to zero amplitude; verify initial state is preserved |
| Known analytic limit | Compare to textbook results in appropriate limits (e.g., Jaynes-Cummings for weak coupling) |
| Unitarity | For closed-system sims, verify `norm(final_state) ≈ 1` or `det(propagator) ≈ 1` |
| Energy conservation | For time-independent Hamiltonians, check `⟨H⟩` is constant |
| Symmetry | If the system has a known symmetry, verify the simulation respects it |

**Action:** Write a sanity-check script in `studies/<study_name>/scripts/` that runs these tests and prints PASS/FAIL for each.

### Check 2 — Convergence

Confirm results are stable with respect to numerical parameters.

**Parameters to sweep:**

| Parameter | Method |
|-----------|--------|
| Hilbert space dimension | Double `N_transmon` and `N_storage`; check fidelity change < threshold |
| Time steps (`nsteps`) | Double step count; check fidelity change < threshold |
| Optimization iterations | Run additional iterations; confirm cost function is plateaued |

**Procedure:**

1. Identify the key numerical parameters in the simulation scripts.
2. Run the [convergence check script](./scripts/convergence_check.py) or write a study-specific version.
3. Record results in a convergence table:

```markdown
| Parameter        | Baseline | Doubled | Δ Fidelity | Converged? |
|------------------|----------|---------|------------|------------|
| N_storage        | 20       | 40      | 1.2e-5     | ✓          |
| N_transmon       | 4        | 8       | 3.1e-6     | ✓          |
| nsteps           | 1000     | 2000    | 8.7e-7     | ✓          |
```

4. Save the convergence data to `studies/<study_name>/data/convergence/`.
5. Generate a convergence plot and save to `studies/<study_name>/figures/`.
6. If an analytic preliminary exists, explicitly compare the numerical result to that limiting-case expectation and record any discrepancy.

### Check 3 — Literature Comparison (if applicable)

For `REP`-class studies or any study with published benchmarks:

1. Extract quantitative values from the reference paper.
2. Run the simulation with matching parameters.
3. Compute agreement metrics:
   - **Fidelity** for state/gate comparisons
   - **Percent error** for scalar quantities (frequencies, decay rates)
   - **Visual overlay** for spectral or time-domain data
4. Report results in a comparison table:

```markdown
| Quantity | Published | Simulated | Error  |
|----------|-----------|-----------|--------|
| χ/2π     | -1.5 MHz  | -1.48 MHz | 1.3%   |
| T_gate   | 250 ns    | 248 ns    | 0.8%   |
```

## Completion Criteria

Update the study README with the validation status:

```markdown
## Validation
- [x] Sanity checks — PASSED (see scripts/sanity_checks.py)
- [x] Convergence — PASSED (see data/convergence/, figures/convergence_*.png)
- [x] Literature comparison — PASSED (see figures/literature_comparison.png)
```

### Check 4 — Reproducibility Artifacts

Before proceeding to the report, verify that key results are saved as machine-readable artifacts:

1. **`artifacts/` directory exists** with at least one file (JSON, NPZ, or CSV).
2. Each artifact contains enough information to reproduce the result without re-running the optimization (optimized parameters, waveforms, gate sequences).
3. Artifact files include metadata: `study_name`, `date_created`, `description`.
4. The report's future `Saved Artifacts` subsection can inventory the files cleanly. If a file cannot be described in one sentence, it is probably not ready for handoff.

### Check 5 — Reproducibility Notebook

Before marking the study COMPLETE, verify:

1. **`scripts/reproducibility_notebook.ipynb` exists** and is a valid Jupyter notebook.
2. The notebook loads saved data/artifacts and verifies key results without requiring expensive re-runs.
3. Every code cell is preceded by a markdown cell explaining what the step does and why.
4. The notebook is self-contained: a user can run it top-to-bottom without reading other files first.
5. Key figures from the report are re-generated or displayed for visual verification.
6. **A dedicated "User-Tunable Parameters" cell** exists early in the notebook (before any simulation code) that exposes all adjustable knobs: Hilbert space dimensions, logical subspace indices, optimizer settings (GRAPE steps, dt, amplitude bounds, restarts, iterations), noise model (T₁, T₂, cavity T₁), cost-function weights, probe states, convergence sweep ranges, and diagnostic settings (Wigner grid, time step). No magic numbers in later cells.
7. **A "Derived Objects" cell** follows the parameters cell and builds all simulation objects (model, subspace, target unitary, noise spec) from those parameters. Re-running these two cells propagates any user change to all downstream code.
8. **Each reproduction step has dual-path cells**: a default "Load saved results" cell and a commented-out "Re-run with current parameters" cell that uses the tunable parameters. Users can uncomment to re-execute with their choices.
9. The summary cell includes a **parameter-effect table** mapping each tunable parameter to its default value and its impact on results.

This check is required by AGENTS.md. See the Reproducibility Requirements section for full details.

Only proceed to report writing (Step 5) when all applicable checks show `[x]` and artifacts are saved.

Do not mark a study `COMPLETE` if the README `## Validation` section is missing, unchecked, or disconnected from the actual evidence in `data/`, `figures/`, and `artifacts/`.
