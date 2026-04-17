---
name: red-green-validation
description: "Test-first validation for cQED simulation studies. Write validation tests that encode known analytic limits, conservation laws, and sanity checks BEFORE running numerics. Confirm they fail against placeholder results, then implement until they pass. Inspired by Red/Green TDD adapted for scientific simulation."
argument-hint: "study=studies/<name> — provide the study path to generate validation tests for"
---

# Red/Green Validation — Test-First Scientific Simulation

## When to Use

- **Before** running any simulation in a new study (Step 2–3 of AGENTS.md workflow)
- When adding a new physics module or extending an existing simulation
- When the reviewer flags missing sanity checks or convergence verification
- At any point you need rigorous, executable validation of simulation results

## Core Idea

Adapted from Simon Willison's Red/Green TDD pattern for coding agents:

1. **RED**: Write executable validation tests that encode what you *know* must be true — analytic limits, conservation laws, symmetry constraints, literature benchmarks. Run them against placeholder/zero results. **They must fail.**
2. **GREEN**: Implement the simulation. Run the same tests. **They must pass.**
3. **REFACTOR**: If tests pass but results look physically unreasonable, add more tests.

This forces the agent to articulate expectations *before* seeing results, preventing post-hoc rationalization of wrong answers.

## Procedure

### Step 1: Identify Testable Constraints

Before writing simulation code, list every constraint you can derive from first principles:

| Category | Examples |
|----------|---------|
| **Limiting cases** | Zero drive → no population transfer; infinite detuning → decoupled; χ → 0 → no dispersive shift |
| **Conservation laws** | Trace of density matrix = 1; unitarity of evolution; total excitation number (if conserved) |
| **Symmetry** | Fidelity invariant under global phase; spectrum symmetric under parity (if applicable) |
| **Analytic results** | Rabi oscillation frequency matches $\Omega_R = \sqrt{\Omega^2 + \Delta^2}$; dispersive shift matches perturbation theory |
| **Known benchmarks** | Published fidelities, spectral lines, or dynamics from literature |
| **Convergence** | Increasing Hilbert space dimension changes result by < tolerance; halving time step changes result by < tolerance |

### Step 2: Write the Validation Script

Create `scripts/test_validation.py` (or use the template from `tools/templates/test_validation_template.py`).

Structure:

```python
"""
Red/Green validation tests for <study_name>.

These tests encode known physical constraints that any correct simulation
must satisfy. They are written BEFORE the simulation code and must FAIL
against placeholder results (RED phase), then PASS after implementation
(GREEN phase).

Usage:
    python test_validation.py          # Run all tests
    python test_validation.py --red    # Verify tests fail against placeholders
"""

import sys
import numpy as np

# ── Tolerances ───────────────────────────────────────────
ATOL = 1e-8   # Absolute tolerance for conservation laws
RTOL = 1e-4   # Relative tolerance for physical quantities
FIDELITY_TOL = 1e-3  # Tolerance for fidelity comparisons

# ── Test Results Tracker ─────────────────────────────────
_results = []

def check(name, condition, detail=""):
    """Record a test result."""
    status = "PASS" if condition else "FAIL"
    _results.append((name, status, detail))
    marker = "  [PASS]" if condition else "  [FAIL]"
    print(f"{marker} {name}" + (f" — {detail}" if detail else ""))
    return condition


# ══════════════════════════════════════════════════════════
# LIMITING CASES
# ══════════════════════════════════════════════════════════

def test_zero_drive_no_transfer():
    """With zero drive amplitude, initial state must be preserved."""
    # Load or compute result for zero drive
    # result = load_result("zero_drive")  # uncomment after implementation
    result = None  # RED: placeholder

    if result is None:
        check("zero_drive_no_transfer", False, "Not yet implemented")
        return

    # Check that initial state fidelity is 1.0
    check(
        "zero_drive_no_transfer",
        abs(result["fidelity"] - 1.0) < ATOL,
        f"fidelity = {result['fidelity']:.10f}, expected 1.0"
    )


# ══════════════════════════════════════════════════════════
# CONSERVATION LAWS
# ══════════════════════════════════════════════════════════

def test_density_matrix_trace():
    """Density matrix trace must equal 1 at all times."""
    # rho_t = load_result("density_matrices")  # uncomment after implementation
    rho_t = None  # RED: placeholder

    if rho_t is None:
        check("density_matrix_trace", False, "Not yet implemented")
        return

    for i, rho in enumerate(rho_t):
        tr = np.real(np.trace(rho))
        if abs(tr - 1.0) > ATOL:
            check("density_matrix_trace", False, f"Tr(rho[{i}]) = {tr}")
            return
    check("density_matrix_trace", True)


def test_unitarity():
    """Evolution operator must be unitary: U^dag U = I."""
    # U = load_result("unitary")  # uncomment after implementation
    U = None  # RED: placeholder

    if U is None:
        check("unitarity", False, "Not yet implemented")
        return

    product = U.conj().T @ U
    identity = np.eye(product.shape[0])
    check(
        "unitarity",
        np.allclose(product, identity, atol=ATOL),
        f"max |U^dag U - I| = {np.max(np.abs(product - identity)):.2e}"
    )


# ══════════════════════════════════════════════════════════
# ANALYTIC BENCHMARKS
# ══════════════════════════════════════════════════════════

def test_analytic_benchmark():
    """Compare key result to analytic prediction."""
    # result = load_result("key_observable")  # uncomment after implementation
    result = None  # RED: placeholder

    if result is None:
        check("analytic_benchmark", False, "Not yet implemented")
        return

    # Example: compare dispersive shift to perturbation theory prediction
    # expected = g**2 / delta  # analytic value
    # check("analytic_benchmark", abs(result - expected) / abs(expected) < RTOL,
    #        f"simulation={result:.6f}, analytic={expected:.6f}")


# ══════════════════════════════════════════════════════════
# CONVERGENCE
# ══════════════════════════════════════════════════════════

def test_convergence_hilbert_space():
    """Result must be stable when increasing Hilbert space dimension."""
    # results = load_result("convergence_sweep")  # uncomment after implementation
    results = None  # RED: placeholder

    if results is None:
        check("convergence_hilbert_space", False, "Not yet implemented")
        return

    # Check that the last two dimension values give results within tolerance
    # vals = [r["observable"] for r in results]
    # check("convergence_hilbert_space",
    #        abs(vals[-1] - vals[-2]) / abs(vals[-1]) < RTOL,
    #        f"delta = {abs(vals[-1] - vals[-2]):.2e}")


# ══════════════════════════════════════════════════════════
# RUNNER
# ══════════════════════════════════════════════════════════

def main():
    red_mode = "--red" in sys.argv

    if red_mode:
        print("=" * 60)
        print("RED PHASE: All tests should FAIL (not yet implemented)")
        print("=" * 60)
    else:
        print("=" * 60)
        print("GREEN PHASE: All tests should PASS")
        print("=" * 60)

    print()

    # Run all test functions
    test_zero_drive_no_transfer()
    test_density_matrix_trace()
    test_unitarity()
    test_analytic_benchmark()
    test_convergence_hilbert_space()

    # Summary
    print()
    print("=" * 60)
    n_pass = sum(1 for _, s, _ in _results if s == "PASS")
    n_fail = sum(1 for _, s, _ in _results if s == "FAIL")
    print(f"Results: {n_pass} PASS, {n_fail} FAIL out of {len(_results)} tests")

    if red_mode:
        if n_fail == len(_results):
            print("RED phase confirmed: all tests fail as expected.")
        else:
            print("WARNING: Some tests passed in RED phase — check placeholders.")
    else:
        if n_fail > 0:
            print("VALIDATION INCOMPLETE: Fix failing tests before proceeding to report.")
            sys.exit(1)
        else:
            print("All validations passed.")

    print("=" * 60)


if __name__ == "__main__":
    main()
```

### Step 3: Run RED Phase

```bash
cd studies/<study_name>
python scripts/test_validation.py --red
```

**All tests must fail.** If any test passes against placeholders, the test is vacuous — fix it.

### Step 4: Implement the Simulation

Write the simulation scripts in `scripts/`. Save results to `data/` and `artifacts/`.

### Step 5: Run GREEN Phase

```bash
python scripts/test_validation.py
```

**All tests must pass.** If any test fails:
- Debug the simulation code (not the test) first
- Only modify a test if the test itself encoded a wrong expectation — document why

### Step 6: Add to Report

In the Validation section of the report, reference the test results:

> "Limiting-case tests (zero drive, trace preservation, unitarity) pass with absolute
> tolerance $10^{-8}$. The dispersive shift agrees with the perturbative prediction
> to within 0.02\%. Hilbert space convergence is confirmed: increasing the Fock
> truncation from $N=8$ to $N=15$ changes the gate fidelity by less than $10^{-5}$."

## Integration with Existing Workflow

This skill fits into Step 2–3 of the AGENTS.md canonical workflow:

```
Step 2: Plan        → Write validation tests (RED phase)
Step 3: Implement   → Run simulations → validation tests (GREEN phase)
Step 4: Validate    → Tests already passed; document in report
```

The validation script becomes a machine-checkable artifact that future agents and the reviewer can re-run.

## Anti-Patterns to Avoid

| Anti-Pattern | Why It's Bad | Do This Instead |
|-------------|-------------|-----------------|
| Writing tests *after* seeing results | Post-hoc rationalization; tests fitted to output | Write tests from analytic limits before implementation |
| Tests with huge tolerances (> 10%) | Vacuous — hides real errors | Use physically motivated tolerances |
| Only testing the happy path | Misses failure regimes | Include off-resonance, strong-drive, and edge cases |
| Modifying tests to match wrong results | Defeats the purpose | Fix the simulation, not the test |
| Skipping RED phase | Can't verify tests are actually constraining | Always confirm tests fail first |
