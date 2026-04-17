"""
Red/Green validation test template for cQED simulation studies.

Copy this file into your study's scripts/ directory and customize
the test functions for your specific physics problem.

These tests encode known physical constraints that any correct simulation
must satisfy. They are written BEFORE the simulation code and must FAIL
against placeholder results (RED phase), then PASS after implementation
(GREEN phase).

Usage:
    python test_validation.py          # Run all tests (GREEN phase)
    python test_validation.py --red    # Verify tests fail against placeholders (RED phase)

See .github/skills/red-green-validation/SKILL.md for full instructions.
"""

import sys
import numpy as np
from pathlib import Path

# ── Paths ────────────────────────────────────────────────
STUDY_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = STUDY_DIR / "data"
ARTIFACTS_DIR = STUDY_DIR / "artifacts"

# ── Tolerances ───────────────────────────────────────────
ATOL = 1e-8   # Absolute tolerance for conservation laws
RTOL = 1e-4   # Relative tolerance for physical quantities
FIDELITY_TOL = 1e-3  # Tolerance for fidelity comparisons

# ── Test Results Tracker ─────────────────────────────────
_results = []


def check(name, condition, detail=""):
    """Record and print a test result."""
    status = "PASS" if condition else "FAIL"
    _results.append((name, status, detail))
    marker = "  [PASS]" if condition else "  [FAIL]"
    print(f"{marker} {name}" + (f" -- {detail}" if detail else ""))
    return condition


def load_artifact(filename):
    """Load a JSON artifact. Returns None if file doesn't exist."""
    import json
    path = ARTIFACTS_DIR / filename
    if not path.exists():
        return None
    with open(path, encoding="utf-8-sig") as f:
        return json.load(f)


def load_npz(filename):
    """Load an NPZ data file. Returns None if file doesn't exist."""
    path = DATA_DIR / filename
    if not path.exists():
        path = ARTIFACTS_DIR / filename
    if not path.exists():
        return None
    return np.load(path, allow_pickle=True)


# ══════════════════════════════════════════════════════════
# LIMITING CASES
# Customize these for your specific physics problem.
# ══════════════════════════════════════════════════════════

def test_zero_drive_no_transfer():
    """With zero drive amplitude, initial state must be preserved."""
    # TODO: Load or compute result for zero drive case
    result = None  # RED: placeholder — replace with actual load

    if result is None:
        check("zero_drive_no_transfer", False, "Not yet implemented")
        return

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
    # TODO: Load density matrices from simulation
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
    # TODO: Load unitary from simulation
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
    """Compare key result to analytic prediction.

    Customize: replace with your specific analytic formula and
    the corresponding simulation observable.
    """
    # TODO: Load result and compare to analytic prediction
    result = None  # RED: placeholder

    if result is None:
        check("analytic_benchmark", False, "Not yet implemented")
        return

    # Example:
    # expected = g**2 / delta
    # check("analytic_benchmark",
    #        abs(result - expected) / abs(expected) < RTOL,
    #        f"sim={result:.6f}, analytic={expected:.6f}, err={abs(result-expected)/abs(expected):.2e}")


# ══════════════════════════════════════════════════════════
# CONVERGENCE
# ══════════════════════════════════════════════════════════

def test_convergence_hilbert_space():
    """Result must be stable when increasing Hilbert space dimension."""
    # TODO: Load convergence sweep results
    results = None  # RED: placeholder

    if results is None:
        check("convergence_hilbert_space", False, "Not yet implemented")
        return

    # Example: check last two dimensions give consistent results
    # vals = [r["observable"] for r in results]
    # rel_change = abs(vals[-1] - vals[-2]) / abs(vals[-1])
    # check("convergence_hilbert_space", rel_change < RTOL,
    #        f"relative change = {rel_change:.2e}")


# ══════════════════════════════════════════════════════════
# RUNNER
# ══════════════════════════════════════════════════════════

ALL_TESTS = [
    test_zero_drive_no_transfer,
    test_density_matrix_trace,
    test_unitarity,
    test_analytic_benchmark,
    test_convergence_hilbert_space,
]


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

    for test_fn in ALL_TESTS:
        test_fn()

    # Summary
    print()
    print("=" * 60)
    n_pass = sum(1 for _, s, _ in _results if s == "PASS")
    n_fail = sum(1 for _, s, _ in _results if s == "FAIL")
    total = len(_results)
    print(f"Results: {n_pass} PASS, {n_fail} FAIL out of {total} tests")

    if red_mode:
        if n_fail == total:
            print("RED phase confirmed: all tests fail as expected.")
        else:
            print(f"WARNING: {n_pass} tests passed in RED phase -- check placeholders.")
    else:
        if n_fail > 0:
            print("VALIDATION INCOMPLETE: Fix failing tests before reporting.")
            sys.exit(1)
        else:
            print("All validations passed.")

    print("=" * 60)


if __name__ == "__main__":
    main()
