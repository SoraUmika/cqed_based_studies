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

Only proceed to report writing (Step 5) when all applicable checks show `[x]`.
