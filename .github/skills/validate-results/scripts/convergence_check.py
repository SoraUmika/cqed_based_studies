"""
Convergence check utility for cQED simulation studies.

This script provides a reusable framework to verify that simulation results
are stable with respect to key numerical parameters (Hilbert space dimension,
time steps, optimization iterations).

Usage:
    Adapt the `CONVERGENCE_PARAMS` dict and `run_simulation` function to your
    specific study, then run:

        python convergence_check.py

    Results are printed as a table and optionally saved to a CSV file.
"""

import numpy as np
from pathlib import Path

# ---------------------------------------------------------------------------
# Configuration — adapt these to your study
# ---------------------------------------------------------------------------

# Each entry: parameter_name -> list of values to test (ascending).
# The first value is the baseline; subsequent values check convergence.
CONVERGENCE_PARAMS = {
    "N_storage": [10, 20, 40],
    "N_transmon": [3, 4, 8],
    "nsteps": [500, 1000, 2000],
}

# Convergence threshold: maximum acceptable change in the metric between
# successive parameter values.
CONVERGENCE_THRESHOLD = 1e-4

# Output directory for convergence data (relative to study root).
OUTPUT_DIR = Path("data/convergence")


def run_simulation(**params):
    """Run the study simulation with the given parameters.

    Returns
    -------
    float
        The metric to track for convergence (e.g., fidelity, energy, cost).

    NOTE: Replace this stub with your actual simulation call.
    """
    raise NotImplementedError(
        "Replace this function with your study's simulation call. "
        "It should accept keyword arguments matching CONVERGENCE_PARAMS "
        "and return a scalar metric (e.g., fidelity)."
    )


def check_convergence():
    """Sweep each parameter and report convergence."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    results = []
    header = f"{'Parameter':<20} {'Baseline':>10} {'Increased':>10} {'Metric_base':>12} {'Metric_inc':>12} {'Delta':>12} {'Converged':>10}"
    print(header)
    print("-" * len(header))

    for param_name, values in CONVERGENCE_PARAMS.items():
        for i in range(len(values) - 1):
            baseline_val = values[i]
            increased_val = values[i + 1]

            # Run with baseline
            metric_base = run_simulation(**{param_name: baseline_val})

            # Run with increased parameter
            metric_inc = run_simulation(**{param_name: increased_val})

            delta = abs(metric_inc - metric_base)
            converged = delta < CONVERGENCE_THRESHOLD

            row = {
                "parameter": param_name,
                "baseline": baseline_val,
                "increased": increased_val,
                "metric_base": metric_base,
                "metric_inc": metric_inc,
                "delta": delta,
                "converged": converged,
            }
            results.append(row)

            status = "PASS" if converged else "FAIL"
            print(
                f"{param_name:<20} {baseline_val:>10} {increased_val:>10} "
                f"{metric_base:>12.6e} {metric_inc:>12.6e} {delta:>12.6e} {status:>10}"
            )

    # Save results to CSV
    csv_path = OUTPUT_DIR / "convergence_results.csv"
    with open(csv_path, "w") as f:
        f.write("parameter,baseline,increased,metric_base,metric_inc,delta,converged\n")
        for r in results:
            f.write(
                f"{r['parameter']},{r['baseline']},{r['increased']},"
                f"{r['metric_base']:.8e},{r['metric_inc']:.8e},"
                f"{r['delta']:.8e},{r['converged']}\n"
            )
    print(f"\nResults saved to {csv_path}")

    # Summary
    all_passed = all(r["converged"] for r in results)
    if all_passed:
        print("\nAll convergence checks PASSED.")
    else:
        failed = [r for r in results if not r["converged"]]
        print(f"\n{len(failed)} convergence check(s) FAILED:")
        for r in failed:
            print(f"  - {r['parameter']}: {r['baseline']} -> {r['increased']}, delta = {r['delta']:.2e}")

    return results


if __name__ == "__main__":
    check_convergence()
