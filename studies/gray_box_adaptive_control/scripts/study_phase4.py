"""
study_phase4.py — Main 3-way comparison: nominal, gray-box, perfect, black-box.

This script is the primary comparison study. For each chi mismatch level:
  1. NOMINAL: GRAPE on learner(chi=chi_prior) → evaluate on truth model
  2. GRAY-BOX: Probe truth → infer chi → GRAPE on learner(chi=chi_hat) → eval on truth
  3. PERFECT: GRAPE on learner(chi=chi_true) → evaluate on truth model
  4. BLACK-BOX: Direct optimizer on truth model (at 30% mismatch only, most expensive)

All evaluations use noiseless (unitary) simulation on the truth model.
Additionally, the 30% mismatch case is run with truth noise to show noisy degradation.

Data is saved to: data/phase4_results.npz

Usage (from study directory):
    python scripts/study_phase4.py

Runtime estimate: ~5-30 min depending on system speed (5 mismatch levels * 3 GRAPE runs).
"""

from __future__ import annotations

import sys
import os
import time
from pathlib import Path

import numpy as np

# Ensure scripts/ is on the path for local imports
STUDY_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(STUDY_DIR / "scripts"))

from models import (
    make_truth_model,
    make_learner_model,
    make_frame,
    make_grape_subspace,
    make_target_matrix,
    CHI_TRUE,
    CHI_PRIOR,
    N_SHOTS_PROBE,
    CONFUSION_MATRIX,
    TRUTH_NOISE,
)
from probe_library import (
    run_chi_ramsey_probe,
    infer_chi_from_probe,
    make_ramsey_delays,
    t2_star_from_noise,
)
from control import run_grape, run_grape_multistart, eval_on_model
from black_box import run_black_box


# ---------------------------------------------------------------------------
# Study parameters
# ---------------------------------------------------------------------------

# Chi mismatch fractions: chi_prior = chi_true * (1 - frac)
# 0% mismatch = perfectly known chi
# 30% mismatch = chi_prior = 0.7 * chi_true (approximately CHI_PRIOR)
MISMATCH_FRACTIONS = [0.0, 0.10, 0.20, 0.30, 0.40]

# GRAPE settings (reduce for faster runs; increase for final results)
GRAPE_MAXITER = 200
GRAPE_N_STEPS = 16
GRAPE_DT_S = 10e-9
GRAPE_SEED_BASE = 42
# Multi-start seeds: run GRAPE with each seed and take best training-model result.
# Seeds 2, 9, 14 consistently give good solutions for both nominal and perfect models.
GRAPE_MULTISTART_SEEDS = [2, 9, 14]

# Black-box settings (run only at 30% mismatch)
# Using small budget (popsize*n_params=3*16=48 population, 10 generations)
# to keep computation tractable. Black-box is a qualitative comparison only.
BB_MISMATCH = 0.30
BB_N_STEPS = 8
BB_SEED = 42
BB_MAXITER_DE = 10
BB_POPSIZE = 3

# Probe settings for gray-box
PROBE_N_DELAYS = 80
PROBE_N_PERIODS = 3.0

# Output directory
DATA_DIR = STUDY_DIR / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Helper: run a full gray-box pipeline
# ---------------------------------------------------------------------------


def run_gray_box_pipeline(
    truth_model,
    truth_frame,
    subspace,
    target,
    chi_initial: float,
    rng: np.random.Generator,
    verbose: bool = True,
) -> dict:
    """
    Full gray-box pipeline: probe → infer → GRAPE on corrected model.

    Parameters
    ----------
    truth_model : DispersiveTransmonCavityModel
        Truth model (probed for chi measurement).
    truth_frame : FrameSpec
        Frame for truth model evaluation.
    subspace : Subspace
        GRAPE subspace.
    target : np.ndarray
        Target matrix.
    chi_initial : float
        Initial chi estimate (learner's prior before probe).
    rng : np.random.Generator
        Random generator for probe shot noise.
    verbose : bool
        Print progress messages.

    Returns
    -------
    dict with keys:
        'chi_hat', 'chi_higher_hat', 'chi_error_frac', 'grape_result',
        'grape_problem', 'fidelity_on_learner', 'avg_fidelity_on_truth',
        'per_state_fidelities', 'probe_data', 'infer_result',
        'avg_fidelity_noisy' (with TRUTH_NOISE)
    """
    t2_star = t2_star_from_noise(TRUTH_NOISE)

    # Step 1: run chi Ramsey probe on truth model
    delays = make_ramsey_delays(chi_initial, n_periods=PROBE_N_PERIODS, n_points=PROBE_N_DELAYS)
    probe_data = run_chi_ramsey_probe(
        chi_true=truth_model.chi,
        chi_higher_true=float(truth_model.chi_higher[0]) if len(truth_model.chi_higher) > 0 else 0.0,
        t2_star=t2_star,
        confusion_matrix=CONFUSION_MATRIX,
        n_shots=N_SHOTS_PROBE,
        delays_s=delays,
        fock_levels=[1, 2, 3],
        rng=rng,
    )
    if verbose:
        print(f"  [probe] delays: 0 to {delays[-1]*1e6:.1f} us, {len(delays)} points")

    # Step 2: infer chi from probe data
    infer_result = infer_chi_from_probe(
        probe_data=probe_data,
        confusion_matrix=CONFUSION_MATRIX,
        n_shots=N_SHOTS_PROBE,
        chi_initial=chi_initial,
        chi_higher_initial=0.0,
    )
    chi_hat = infer_result["chi"]
    chi_higher_hat = infer_result["chi_higher"]
    chi_error_frac = abs(chi_hat - truth_model.chi) / abs(truth_model.chi)
    if verbose:
        print(f"  [infer] chi_hat = {chi_hat/(2*np.pi)/1e6:.4f} MHz "
              f"(true: {truth_model.chi/(2*np.pi)/1e6:.4f} MHz, error: {chi_error_frac*100:.2f}%)")

    # Step 3: GRAPE on corrected learner model
    learner_corrected = make_learner_model(chi=chi_hat, chi_higher_val=0.0, kerr_val=0.0)
    learner_frame = make_frame(learner_corrected)
    if verbose:
        print(f"  [GRAPE] Running on corrected model (chi={chi_hat/(2*np.pi)/1e6:.4f} MHz)...")

    grape_result, problem, fidelity_on_learner, gb_best_seed = run_grape_multistart(
        learner_corrected, learner_frame, subspace, target,
        seeds=GRAPE_MULTISTART_SEEDS,
        n_steps=GRAPE_N_STEPS, dt_s=GRAPE_DT_S,
        maxiter=GRAPE_MAXITER,
    )
    if verbose:
        print(f"  [GRAPE] Fidelity on corrected model: {fidelity_on_learner:.6f} [best seed={gb_best_seed}]")

    # Step 4: evaluate on truth model (noiseless)
    avg_f, per_f = eval_on_model(grape_result, problem, truth_model, truth_frame, eval_noise=None)
    if verbose:
        print(f"  [eval] Avg fidelity on truth (noiseless): {avg_f:.6f}")

    # Step 5: evaluate on truth model (noisy)
    avg_f_noisy, _ = eval_on_model(grape_result, problem, truth_model, truth_frame, eval_noise=TRUTH_NOISE)
    if verbose:
        print(f"  [eval] Avg fidelity on truth (noisy): {avg_f_noisy:.6f}")

    return {
        "chi_hat": chi_hat,
        "chi_higher_hat": chi_higher_hat,
        "chi_error_frac": chi_error_frac,
        "grape_result": grape_result,
        "grape_problem": problem,
        "fidelity_on_learner": fidelity_on_learner,
        "avg_fidelity_on_truth": avg_f,
        "per_state_fidelities": per_f,
        "avg_fidelity_noisy": avg_f_noisy,
        "probe_data": probe_data,
        "infer_result": infer_result,
    }


# ---------------------------------------------------------------------------
# Main comparison loop
# ---------------------------------------------------------------------------


def run_phase4(verbose: bool = True):
    """
    Run the main 3-way + black-box comparison.

    Iterates over mismatch fractions, runs nominal/gray-box/perfect GRAPE,
    collects results, and saves to data/phase4_results.npz.
    """
    t0_total = time.perf_counter()

    truth_model = make_truth_model()
    truth_frame = make_frame(truth_model)
    subspace = make_grape_subspace()
    target = make_target_matrix()
    rng = np.random.default_rng(12345)

    n_mismatch = len(MISMATCH_FRACTIONS)

    # Result arrays
    nominal_fidelities = np.zeros(n_mismatch)
    gray_box_fidelities = np.zeros(n_mismatch)
    perfect_fidelities = np.zeros(n_mismatch)
    nominal_fidelities_noisy = np.zeros(n_mismatch)
    gray_box_fidelities_noisy = np.zeros(n_mismatch)
    perfect_fidelities_noisy = np.zeros(n_mismatch)

    chi_hat_values = np.zeros(n_mismatch)
    chi_error_fracs = np.zeros(n_mismatch)
    chi_prior_values = np.zeros(n_mismatch)

    # Per-state fidelity arrays (for 30% mismatch case, stored separately)
    per_state_fidelities_30pct = {
        "nominal": None,
        "gray_box": None,
        "perfect": None,
    }

    # Black-box results (only at 30% mismatch)
    bb_fidelity = np.nan
    bb_fidelity_history = []
    bb_n_evaluations = 0

    for i, mismatch in enumerate(MISMATCH_FRACTIONS):
        if verbose:
            print(f"\n{'='*60}")
            print(f"Mismatch {mismatch*100:.0f}% (step {i+1}/{n_mismatch})")
            print(f"{'='*60}")

        # Chi prior for this mismatch level
        chi_prior_i = CHI_TRUE * (1.0 - mismatch)
        chi_prior_values[i] = chi_prior_i

        if verbose:
            print(f"chi_true  = {CHI_TRUE/(2*np.pi)/1e6:.4f} MHz")
            print(f"chi_prior = {chi_prior_i/(2*np.pi)/1e6:.4f} MHz  ({mismatch*100:.0f}% below true)")

        # ---- NOMINAL ----
        if verbose:
            print("\n[NOMINAL] GRAPE on prior model -> evaluate on truth")
        learner_prior = make_learner_model(chi=chi_prior_i)
        learner_prior_frame = make_frame(learner_prior)
        grape_nominal, prob_nominal, f_nominal_learner, nom_best_seed = run_grape_multistart(
            learner_prior, learner_prior_frame, subspace, target,
            seeds=GRAPE_MULTISTART_SEEDS,
            n_steps=GRAPE_N_STEPS, dt_s=GRAPE_DT_S,
            maxiter=GRAPE_MAXITER,
        )
        if verbose:
            print(f"  Fidelity on learner (prior): {f_nominal_learner:.6f} [best seed={nom_best_seed}]")

        nom_f, nom_per_f = eval_on_model(grape_nominal, prob_nominal, truth_model, truth_frame, eval_noise=None)
        nom_f_noisy, _ = eval_on_model(grape_nominal, prob_nominal, truth_model, truth_frame, eval_noise=TRUTH_NOISE)
        nominal_fidelities[i] = nom_f
        nominal_fidelities_noisy[i] = nom_f_noisy
        if verbose:
            print(f"  Avg fidelity on truth (noiseless): {nom_f:.6f}")
            print(f"  Avg fidelity on truth (noisy):    {nom_f_noisy:.6f}")

        if abs(mismatch - 0.30) < 0.001:
            per_state_fidelities_30pct["nominal"] = nom_per_f

        # ---- GRAY-BOX ----
        if verbose:
            print("\n[GRAY-BOX] Probe truth -> infer chi -> GRAPE on corrected model")
        gb_result = run_gray_box_pipeline(
            truth_model, truth_frame, subspace, target,
            chi_initial=chi_prior_i, rng=rng, verbose=verbose,
        )
        gray_box_fidelities[i] = gb_result["avg_fidelity_on_truth"]
        gray_box_fidelities_noisy[i] = gb_result["avg_fidelity_noisy"]
        chi_hat_values[i] = gb_result["chi_hat"]
        chi_error_fracs[i] = gb_result["chi_error_frac"]

        if abs(mismatch - 0.30) < 0.001:
            per_state_fidelities_30pct["gray_box"] = gb_result["per_state_fidelities"]

        # ---- PERFECT ----
        if verbose:
            print("\n[PERFECT] GRAPE on true chi model -> evaluate on truth")
        learner_perfect = make_learner_model(chi=CHI_TRUE)
        learner_perfect_frame = make_frame(learner_perfect)
        grape_perfect, prob_perfect, f_perfect_learner, perf_best_seed = run_grape_multistart(
            learner_perfect, learner_perfect_frame, subspace, target,
            seeds=GRAPE_MULTISTART_SEEDS,
            n_steps=GRAPE_N_STEPS, dt_s=GRAPE_DT_S,
            maxiter=GRAPE_MAXITER,
        )
        perf_f, perf_per_f = eval_on_model(grape_perfect, prob_perfect, truth_model, truth_frame, eval_noise=None)
        perf_f_noisy, _ = eval_on_model(grape_perfect, prob_perfect, truth_model, truth_frame, eval_noise=TRUTH_NOISE)
        perfect_fidelities[i] = perf_f
        perfect_fidelities_noisy[i] = perf_f_noisy
        if verbose:
            print(f"  Fidelity on perfect learner: {f_perfect_learner:.6f} [best seed={perf_best_seed}]")
            print(f"  Avg fidelity on truth (noiseless): {perf_f:.6f}")
            print(f"  Avg fidelity on truth (noisy):    {perf_f_noisy:.6f}")

        if abs(mismatch - 0.30) < 0.001:
            per_state_fidelities_30pct["perfect"] = perf_per_f

        # ---- BLACK-BOX (only at 30%) ----
        if abs(mismatch - BB_MISMATCH) < 0.001:
            if verbose:
                print("\n[BLACK-BOX] Direct optimization on truth model (no model learning)...")
            bb_result = run_black_box(
                truth_model=truth_model,
                frame=truth_frame,
                grape_subspace=subspace,
                target_matrix=target,
                seed=BB_SEED,
                n_steps=BB_N_STEPS,
                dt_s=GRAPE_DT_S,
                maxiter_de=BB_MAXITER_DE,
                popsize=BB_POPSIZE,
            )
            bb_fidelity = bb_result["best_fidelity"]
            bb_fidelity_history = bb_result["fidelity_history"]
            bb_n_evaluations = bb_result["n_evaluations_used"]
            if verbose:
                print(f"  BB best fidelity: {bb_fidelity:.6f}")
                print(f"  BB evaluations used: {bb_n_evaluations}")

    total_time = time.perf_counter() - t0_total
    if verbose:
        print(f"\n{'='*60}")
        print("Phase 4 complete")
        print(f"Total wall time: {total_time:.1f} s")

    # ---- Summary table ----
    if verbose:
        print(f"\n{'Mismatch':>10} {'Nominal':>10} {'GrayBox':>10} {'Perfect':>10}")
        print("-" * 45)
        for i, frac in enumerate(MISMATCH_FRACTIONS):
            print(f"{frac*100:>9.0f}% "
                  f"{nominal_fidelities[i]:>10.6f} "
                  f"{gray_box_fidelities[i]:>10.6f} "
                  f"{perfect_fidelities[i]:>10.6f}")
        if not np.isnan(bb_fidelity):
            print(f"\nBlack-box at {BB_MISMATCH*100:.0f}%: {bb_fidelity:.6f} ({bb_n_evaluations} evals)")

    # ---- Save results ----
    np.savez(
        DATA_DIR / "phase4_results.npz",
        mismatch_fractions=np.array(MISMATCH_FRACTIONS),
        chi_true=np.array(CHI_TRUE),
        chi_prior_values=chi_prior_values,
        chi_hat_values=chi_hat_values,
        chi_error_fracs=chi_error_fracs,
        nominal_fidelities=nominal_fidelities,
        gray_box_fidelities=gray_box_fidelities,
        perfect_fidelities=perfect_fidelities,
        nominal_fidelities_noisy=nominal_fidelities_noisy,
        gray_box_fidelities_noisy=gray_box_fidelities_noisy,
        perfect_fidelities_noisy=perfect_fidelities_noisy,
        bb_fidelity=np.array(bb_fidelity),
        bb_fidelity_history=np.array(bb_fidelity_history) if bb_fidelity_history else np.array([np.nan]),
        bb_n_evaluations=np.array(bb_n_evaluations),
        bb_mismatch=np.array(BB_MISMATCH),
        per_state_nominal_30=np.array(per_state_fidelities_30pct["nominal"] or [np.nan]),
        per_state_gray_box_30=np.array(per_state_fidelities_30pct["gray_box"] or [np.nan]),
        per_state_perfect_30=np.array(per_state_fidelities_30pct["perfect"] or [np.nan]),
        total_wall_time_s=np.array(total_time),
        grape_n_steps=np.array(GRAPE_N_STEPS),
        grape_dt_ns=np.array(GRAPE_DT_S * 1e9),
        grape_maxiter=np.array(GRAPE_MAXITER),
    )
    if verbose:
        print(f"\nSaved to: {DATA_DIR / 'phase4_results.npz'}")

    return {
        "mismatch_fractions": MISMATCH_FRACTIONS,
        "nominal_fidelities": nominal_fidelities,
        "gray_box_fidelities": gray_box_fidelities,
        "perfect_fidelities": perfect_fidelities,
        "nominal_fidelities_noisy": nominal_fidelities_noisy,
        "gray_box_fidelities_noisy": gray_box_fidelities_noisy,
        "perfect_fidelities_noisy": perfect_fidelities_noisy,
        "chi_hat_values": chi_hat_values,
        "chi_error_fracs": chi_error_fracs,
        "bb_fidelity": bb_fidelity,
        "bb_fidelity_history": bb_fidelity_history,
        "total_wall_time_s": total_time,
    }


if __name__ == "__main__":
    results = run_phase4(verbose=True)
