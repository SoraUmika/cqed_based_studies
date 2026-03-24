"""
diag_phase4_quick.py — Quick validation of the phase4 study logic.

Runs only the 0% and 30% mismatch cases (no black-box) to verify:
1. Multi-start GRAPE works correctly
2. Direction is correct (perfect >= gray-box >= nominal)
3. Study is ready for full run

Usage: python scripts/diag_phase4_quick.py
"""

import sys
import numpy as np
from pathlib import Path

STUDY_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(STUDY_DIR / "scripts"))

from models import (
    make_truth_model, make_learner_model, make_frame,
    make_grape_subspace, make_target_matrix,
    CHI_TRUE, TRUTH_NOISE, N_SHOTS_PROBE, CONFUSION_MATRIX,
)
from probe_library import run_chi_ramsey_probe, infer_chi_from_probe, make_ramsey_delays, t2_star_from_noise
from control import run_grape_multistart, eval_on_model

SEEDS = [2, 9, 14]
MISMATCH_FRACTIONS = [0.0, 0.30]
N_STEPS = 16
DT_S = 10e-9
MAXITER = 200

truth_model = make_truth_model()
truth_frame = make_frame(truth_model)
subspace = make_grape_subspace()
target = make_target_matrix()
rng = np.random.default_rng(12345)
t2_star = t2_star_from_noise(TRUTH_NOISE)

print("=" * 65)
print("Phase 4 quick validation (0% and 30% mismatch)")
print(f"Multi-start seeds: {SEEDS}, maxiter={MAXITER}")
print("=" * 65)

for mismatch in MISMATCH_FRACTIONS:
    chi_prior = CHI_TRUE * (1.0 - mismatch)
    print(f"\n{'='*55}")
    print(f"Mismatch = {mismatch*100:.0f}%  chi_prior = {chi_prior/(2*np.pi)/1e6:.3f} MHz")
    print(f"{'='*55}")

    # ---- NOMINAL ----
    print("[NOMINAL] Multi-start GRAPE on prior model...")
    model_nom = make_learner_model(chi=chi_prior)
    frame_nom = make_frame(model_nom)
    gr_nom, prob_nom, f_nom_train, seed_nom = run_grape_multistart(
        model_nom, frame_nom, subspace, target,
        seeds=SEEDS, n_steps=N_STEPS, dt_s=DT_S, maxiter=MAXITER,
    )
    nom_f, nom_per_f = eval_on_model(gr_nom, prob_nom, truth_model, truth_frame)
    print(f"  Train F = {f_nom_train:.4f} [seed={seed_nom}], Truth F = {nom_f:.4f}")

    # ---- GRAY-BOX ----
    if mismatch > 0:
        print("[GRAY-BOX] Probe -> infer -> GRAPE...")
        delays = make_ramsey_delays(chi_prior, n_periods=3.0, n_points=80)
        probe_data = run_chi_ramsey_probe(
            chi_true=truth_model.chi,
            chi_higher_true=float(truth_model.chi_higher[0]) if truth_model.chi_higher else 0.0,
            t2_star=t2_star, confusion_matrix=CONFUSION_MATRIX,
            n_shots=N_SHOTS_PROBE, delays_s=delays, fock_levels=[1, 2, 3], rng=rng,
        )
        infer_result = infer_chi_from_probe(
            probe_data, CONFUSION_MATRIX, N_SHOTS_PROBE, chi_initial=chi_prior,
        )
        chi_hat = infer_result["chi"]
        print(f"  chi_hat = {chi_hat/(2*np.pi)/1e6:.4f} MHz (true: {CHI_TRUE/(2*np.pi)/1e6:.4f} MHz)")
        model_gb = make_learner_model(chi=chi_hat)
        frame_gb = make_frame(model_gb)
        gr_gb, prob_gb, f_gb_train, seed_gb = run_grape_multistart(
            model_gb, frame_gb, subspace, target,
            seeds=SEEDS, n_steps=N_STEPS, dt_s=DT_S, maxiter=MAXITER,
        )
        gb_f, _ = eval_on_model(gr_gb, prob_gb, truth_model, truth_frame)
        print(f"  Train F = {f_gb_train:.4f} [seed={seed_gb}], Truth F = {gb_f:.4f}")
    else:
        gb_f = nom_f
        print("[GRAY-BOX] 0% mismatch => same as nominal")

    # ---- PERFECT ----
    print("[PERFECT] Multi-start GRAPE on true chi model...")
    model_perf = make_learner_model(chi=CHI_TRUE)
    frame_perf = make_frame(model_perf)
    gr_perf, prob_perf, f_perf_train, seed_perf = run_grape_multistart(
        model_perf, frame_perf, subspace, target,
        seeds=SEEDS, n_steps=N_STEPS, dt_s=DT_S, maxiter=MAXITER,
    )
    perf_f, perf_per_f = eval_on_model(gr_perf, prob_perf, truth_model, truth_frame)
    print(f"  Train F = {f_perf_train:.4f} [seed={seed_perf}], Truth F = {perf_f:.4f}")

    print(f"\n  Summary: Nominal={nom_f:.4f}  GrayBox={gb_f:.4f}  Perfect={perf_f:.4f}")
    print(f"  Gap (gray-box vs nominal): {gb_f - nom_f:+.4f}")
    print(f"  Gap (perfect vs nominal):  {perf_f - nom_f:+.4f}")
    print(f"  Per-state fidelities (perfect): {[f'{v:.3f}' for v in perf_per_f]}")
