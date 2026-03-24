"""
diag_multiseed.py — Test multiple GRAPE seeds to find best convergence.

For both nominal (chi=-2.0 MHz) and perfect (chi=-2.84 MHz) models,
run GRAPE with seeds 0-14 and maxiter=400 to find the best achievable fidelity.

This identifies:
1. Whether local minima are the issue (sensitivity to seed)
2. The achievable chi sensitivity gap with better convergence
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

print("=" * 65)
print("Multi-seed GRAPE convergence test (maxiter=400, seeds 0-14)")
print("=" * 65)

truth_model = make_truth_model()
truth_frame = make_frame(truth_model)
subspace = make_grape_subspace()
target = make_target_matrix()
SEEDS = list(range(15))
MAXITER = 400

best = {}
for label, chi_val in [("Nominal", CHI_PRIOR), ("Perfect", CHI_TRUE)]:
    print(f"\n{'='*40}")
    print(f"{label}: chi = {chi_val/(2*np.pi)/1e6:.3f} MHz")
    print(f"{'='*40}")
    print(f"{'seed':>5}  {'train_F':>10}  {'truth_F':>10}  {'obj_val':>10}")
    print("-" * 42)

    model = make_learner_model(chi=chi_val)
    frame = make_frame(model)
    best_truth_f = 0.0
    best_seed = -1

    for seed in SEEDS:
        result, problem, f_train = run_grape(
            model, frame, subspace, target,
            n_steps=16, dt_s=10e-9, maxiter=MAXITER, seed=seed
        )
        avg_f, _ = eval_on_model(result, problem, truth_model, truth_frame)
        obj_val = result.objective_value
        print(f"  {seed:3d}  {f_train:10.6f}  {avg_f:10.6f}  {obj_val:10.6f}")

        if avg_f > best_truth_f:
            best_truth_f = avg_f
            best_seed = seed
            best[label] = {
                "seed": seed, "train_F": f_train, "truth_F": avg_f,
                "result": result, "problem": problem
            }

    print(f"\nBest seed for {label}: seed={best_seed}, truth_F={best_truth_f:.6f}")

print("\n" + "="*65)
print("Summary: Best results per model")
print("="*65)
for label in ["Nominal", "Perfect"]:
    if label in best:
        d = best[label]
        print(f"{label:10s}: seed={d['seed']}, train_F={d['train_F']:.6f}, truth_F={d['truth_F']:.6f}")
if "Nominal" in best and "Perfect" in best:
    gap = best["Perfect"]["truth_F"] - best["Nominal"]["truth_F"]
    print(f"\nGap (perfect - nominal) on truth model: {gap:.6f}")
    print(f"Relative infidelity reduction: {gap/(1-best['Nominal']['truth_F'])*100:.1f}%")
