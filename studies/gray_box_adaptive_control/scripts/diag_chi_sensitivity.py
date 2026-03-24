"""
diag_chi_sensitivity.py — Quick diagnostic for chi sensitivity of GRAPE.

Tests whether the unitary fidelity on the truth model is meaningfully
sensitive to chi mismatch using UnitaryObjective (current control.py settings).

Runs nominal (chi=-2.0 MHz) and perfect (chi=-2.84 MHz) GRAPE at T=160ns
and evaluates both on the truth model.
"""

import sys
import numpy as np
from pathlib import Path

STUDY_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(STUDY_DIR / "scripts"))

from models import (
    make_truth_model, make_learner_model, make_frame,
    make_grape_subspace, make_target_matrix,
    CHI_TRUE, CHI_PRIOR, TRUTH_NOISE,
)
from control import run_grape, eval_on_model, compute_subspace_fidelity_from_result

print("=" * 60)
print("Chi sensitivity diagnostic")
print("=" * 60)

truth_model = make_truth_model()
truth_frame = make_frame(truth_model)
subspace = make_grape_subspace()
target = make_target_matrix()

# Test two chi values: 30% wrong (nominal) and exact (perfect)
chi_tests = [
    ("Nominal (30% wrong)", CHI_PRIOR),
    ("Perfect (oracle)", CHI_TRUE),
]

for label, chi_val in chi_tests:
    print(f"\n--- {label}: chi = {chi_val/(2*np.pi)/1e6:.3f} MHz ---")
    model = make_learner_model(chi=chi_val)
    frame = make_frame(model)

    result, problem, f_learner = run_grape(
        model, frame, subspace, target,
        n_steps=16, dt_s=10e-9, maxiter=200, seed=42
    )
    f_exact = compute_subspace_fidelity_from_result(result)

    print(f"  Training model fidelity (extract): {f_learner:.6f}")
    print(f"  Exact unitary fidelity (system_metrics): {f_exact:.6f}")
    print(f"  GRAPE success: {result.success}")
    print(f"  Objective value: {result.objective_value:.6f}")

    # Print all available metrics
    print(f"  Metrics keys: {list(result.metrics.keys())}")
    for k, v in result.metrics.items():
        print(f"    {k}: {v}")

    # Evaluate on truth model (noiseless)
    avg_f, per_f = eval_on_model(result, problem, truth_model, truth_frame, eval_noise=None)
    print(f"  Avg fidelity on truth (noiseless): {avg_f:.6f}")
    print(f"  Per-state fidelities: {[f'{f:.4f}' for f in per_f]}")

    # Evaluate with noise
    avg_f_noisy, _ = eval_on_model(result, problem, truth_model, truth_frame, eval_noise=TRUTH_NOISE)
    print(f"  Avg fidelity on truth (noisy): {avg_f_noisy:.6f}")

print("\n--- Summary ---")
print("Compare nominal vs perfect to verify chi sensitivity direction.")
