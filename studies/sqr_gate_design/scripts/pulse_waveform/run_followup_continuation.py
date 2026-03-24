"""
Continuation script: assemble Phase B + C results from terminal output,
then run remaining Phases D, E, and C-EXT with proper memory management.

Fixes:
  - gc.collect() + del between optimizations to prevent memory leak
  - Incremental saving after each sub-phase
  - Reduced maxiter for expensive optimizers

Usage:
    python scripts/run_followup_continuation.py
"""
from __future__ import annotations

import gc
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.linalg import block_diag

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import runtime_compat  # noqa: F401

from common import (
    CHI, N_CAV, N_TR, DT,
    build_frame, build_model, duration_from_chi_t,
    target_qubit_unitary, z_corrected_target_fidelity,
    spectator_z_fidelity, spectator_transverse_error,
    identity_fidelity_with_z, extract_branch_unitaries, extract_leakage,
)
from cqed_sim.core.frequencies import (
    carrier_for_transition_frequency,
    manifold_transition_frequency,
)
from cqed_sim.sequence import SequenceCompiler
from cqed_sim.sim import SimulationConfig, prepare_simulation

# Import builders and helpers from the main script
from run_followup_multitone import (
    CHI_T_VALUES, REP_N0, REP_THETA, REP_PHI, REP_LOGICAL_N,
    SIGMA_FRACTION, OPT_METHOD,
    TARGET_BRANCHES, THETA_VALUES, PHI_VALUES,
    build_single_tone_gaussian,
    build_cosine_squared_pulse,
    build_independent_tone_multitone,
    build_segmented_multitone,
    build_smooth_basis_multitone,
    simulate_and_extract,
    compute_all_metrics,
    normalized_gaussian,
    MultitoneTone,
)

DATA_DIR = SCRIPT_DIR.parent / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)
OUTPATH = DATA_DIR / "followup_multitone_results.npz"

# Lower maxiter for continuation phases to control memory / time
CONT_MAXITER = 150


# ===================================================================
# Hardcoded Phase B + C results from terminal output
# ===================================================================
def assemble_existing_results():
    """Build results dict from Phase B + C terminal data."""
    results = {}
    results["chi_t_values"] = CHI_T_VALUES
    results["target_branches"] = np.array(TARGET_BRANCHES)
    results["theta_values"] = np.array(THETA_VALUES)
    results["phi_values"] = np.array(PHI_VALUES)
    results["logical_n_values"] = np.array([3, 4])
    n_chi = len(CHI_T_VALUES)

    # Phase B baselines: shape (3, n_chi)
    families_B = ["single_tone_gaussian", "cosine_squared", "multitone_baseline"]
    results["baseline_families"] = np.array(families_B, dtype=object)

    # Data from terminal: [F_strict, F_block, F_cphase] per family × chi_t
    B_data = {
        "single_tone_gaussian": [
            (0.0492, 0.0789, 0.5548),  # 0.5
            (0.5463, 0.5497, 0.8527),  # 1.0
            (0.1913, 0.2794, 0.9718),  # 1.5
            (0.8691, 0.8743, 0.9978),  # 2.0
            (0.9243, 0.9297, 1.0000),  # 3.0
            (0.9544, 0.9597, 1.0000),  # 5.0
        ],
        "cosine_squared": [
            (0.0556, 0.0926, 0.6101),
            (0.6226, 0.6264, 0.9049),
            (0.2147, 0.3047, 0.9963),
            (0.8920, 0.8973, 0.9981),
            (0.9338, 0.9393, 1.0000),
            (0.9584, 0.9637, 1.0000),
        ],
        "multitone_baseline": [
            (0.0492, 0.0789, 0.5548),
            (0.5463, 0.5497, 0.8527),
            (0.1913, 0.2794, 0.9718),
            (0.8691, 0.8743, 0.9978),
            (0.9243, 0.9297, 1.0000),
            (0.9544, 0.9597, 1.0000),
        ],
    }

    # Initialize baseline metric arrays — we only had strict, block, cphase printed
    # For other metrics we need to re-compute a few key ones (or set NaN)
    for metric in ["strict_logical_fid", "block_phase_relaxed_fid",
                    "branch_cphase_mean", "target_branch_fid", "branch_true_mean",
                    "same_block_pop_mean", "leakage_max",
                    "spectator_phase_spread", "spectator_max_transverse"]:
        results[f"baseline_{metric}"] = np.full((3, n_chi), np.nan)

    for fi, fam in enumerate(families_B):
        for ci in range(n_chi):
            strict, block, cphase = B_data[fam][ci]
            results["baseline_strict_logical_fid"][fi, ci] = strict
            results["baseline_block_phase_relaxed_fid"][fi, ci] = block
            results["baseline_branch_cphase_mean"][fi, ci] = cphase

    # Phase C data from terminal: [F_block, F_strict] per optimizer × chi_t
    C_data = {
        "opt_indep": [
            (0.4419, 0.2017),  # 0.5
            (0.7918, 0.7868),  # 1.0
            (0.3442, 0.1911),  # 1.5
            (0.8997, 0.8943),  # 2.0
            (0.9299, 0.9245),  # 3.0
            (0.9597, 0.9544),  # 5.0
        ],
        "opt_detuned": [
            (0.4457, 0.1944),
            (0.8063, 0.8012),
            (0.3423, 0.1898),
            (0.8985, 0.8931),
            (0.9311, 0.9257),
            (np.nan, np.nan),  # χT/2π=5.0 not completed
        ],
        "opt_smooth": [
            (0.3034, 0.2558),
            (0.7303, 0.7255),
            (0.3100, 0.2467),
            (0.9009, 0.8954),
            (0.9379, 0.9324),
            (np.nan, np.nan),  # χT/2π=5.0 not completed
        ],
    }

    for prefix in ["opt_indep", "opt_detuned", "opt_smooth"]:
        for metric in ["target_branch_fid", "branch_true_mean", "branch_cphase_mean",
                        "strict_logical_fid", "block_phase_relaxed_fid",
                        "same_block_pop_mean", "leakage_max",
                        "spectator_phase_spread", "spectator_max_transverse"]:
            results[f"{prefix}_{metric}"] = np.full(n_chi, np.nan)
        results[f"{prefix}_opt_cost"] = np.full(n_chi, np.nan)

        for ci in range(n_chi):
            block, strict = C_data[prefix][ci]
            results[f"{prefix}_block_phase_relaxed_fid"][ci] = block
            results[f"{prefix}_strict_logical_fid"][ci] = strict

    return results


# ===================================================================
# Re-compute detailed metrics for baselines (quick, no optimization)
# ===================================================================
def fill_baseline_details(results, model, frame):
    """Re-compute full baseline metrics (fast, no optimization)."""
    print("Recomputing baseline detailed metrics...")
    families = ["single_tone_gaussian", "cosine_squared", "multitone_baseline"]
    builders = {
        "single_tone_gaussian": lambda T: build_single_tone_gaussian(
            model, frame, REP_N0, REP_THETA, REP_PHI, T),
        "cosine_squared": lambda T: build_cosine_squared_pulse(
            model, frame, REP_N0, REP_THETA, REP_PHI, T),
        "multitone_baseline": lambda T: build_independent_tone_multitone(
            model, frame, REP_LOGICAL_N, REP_N0, REP_THETA, REP_PHI, T),
    }

    for fi, fam in enumerate(families):
        for ci, chi_t in enumerate(CHI_T_VALUES):
            T = duration_from_chi_t(chi_t)
            pulses, dops, T_tot = builders[fam](T)
            full_op, states = simulate_and_extract(
                model, frame, pulses, dops, REP_LOGICAL_N, T_tot)
            m = compute_all_metrics(
                full_op, states, model, REP_LOGICAL_N, REP_N0, REP_THETA, REP_PHI)
            for metric in ["target_branch_fid", "branch_true_mean", "branch_cphase_mean",
                           "strict_logical_fid", "block_phase_relaxed_fid",
                           "same_block_pop_mean", "leakage_max",
                           "spectator_phase_spread", "spectator_max_transverse"]:
                results[f"baseline_{metric}"][fi, ci] = m[metric]
            del full_op, states, m
        gc.collect()
        print(f"  {fam} done")

    return results


# ===================================================================
# Memory-safe optimization wrappers (with gc.collect)
# ===================================================================
def optimize_independent_tone_safe(model, frame, logical_n, n0, theta, phi, T,
                                    objective="block_phase_relaxed_fid",
                                    maxiter=CONT_MAXITER):
    n_free = int(logical_n) - 1
    non_target = [n for n in range(int(logical_n)) if n != int(n0)]
    x0 = np.zeros(2 * n_free)

    eval_count = [0]

    def cost(x):
        eval_count[0] += 1
        amp_ratios = np.zeros(int(logical_n))
        phases_arr = np.zeros(int(logical_n))
        amp_ratios[int(n0)] = 1.0
        phases_arr[int(n0)] = float(phi)
        for i, n in enumerate(non_target):
            amp_ratios[n] = x[i]
            phases_arr[n] = x[n_free + i]
        try:
            pulses, drive_ops, T_tot = build_independent_tone_multitone(
                model, frame, logical_n, n0, theta, phi, T,
                amp_ratios=amp_ratios, phases=phases_arr
            )
            full_op, final_states = simulate_and_extract(
                model, frame, pulses, drive_ops, logical_n, T_tot
            )
            metrics = compute_all_metrics(
                full_op, final_states, model, logical_n, n0, theta, phi
            )
            val = -metrics[objective]
            del full_op, final_states, metrics, pulses
            if eval_count[0] % 50 == 0:
                gc.collect()
            return val
        except Exception:
            return 0.0

    result = minimize(cost, x0, method=OPT_METHOD,
                      options={"maxiter": maxiter, "xatol": 1e-5, "fatol": 1e-6})
    gc.collect()

    amp_ratios = np.zeros(int(logical_n))
    phases_arr = np.zeros(int(logical_n))
    amp_ratios[int(n0)] = 1.0
    phases_arr[int(n0)] = float(phi)
    for i, n in enumerate(non_target):
        amp_ratios[n] = result.x[i]
        phases_arr[n] = result.x[n_free + i]
    return amp_ratios, phases_arr, result


def optimize_2seg_safe(model, frame, logical_n, n0, theta, phi, T,
                        objective="block_phase_relaxed_fid",
                        maxiter=CONT_MAXITER):
    n_tones = int(logical_n)
    non_target = [n for n in range(n_tones) if n != int(n0)]
    n_free = len(non_target)
    n_per_seg = n_free + n_tones
    x0 = np.zeros(2 * n_per_seg)
    x0[n_free + int(n0)] = float(phi)
    x0[n_per_seg + n_free + int(n0)] = -float(phi)

    eval_count = [0]

    def cost(x):
        eval_count[0] += 1
        seg_params = []
        for seg in range(2):
            offset = seg * n_per_seg
            amp_ratios = np.zeros(n_tones)
            phases_arr = np.zeros(n_tones)
            amp_ratios[int(n0)] = 1.0
            for i, n in enumerate(non_target):
                amp_ratios[n] = x[offset + i]
            for n in range(n_tones):
                phases_arr[n] = x[offset + n_free + n]
            seg_params.append({
                "amp_ratios": amp_ratios,
                "phases": phases_arr,
                "detunings": None,
            })
        try:
            pulses, drive_ops, T_tot = build_segmented_multitone(
                model, frame, logical_n, n0, theta, phi, T,
                n_segments=2, segment_params_list=seg_params
            )
            full_op, final_states = simulate_and_extract(
                model, frame, pulses, drive_ops, logical_n, T_tot
            )
            metrics = compute_all_metrics(
                full_op, final_states, model, logical_n, n0, theta, phi
            )
            val = -metrics[objective]
            del full_op, final_states, metrics, pulses
            if eval_count[0] % 50 == 0:
                gc.collect()
            return val
        except Exception:
            return 0.0

    result = minimize(cost, x0, method=OPT_METHOD,
                      options={"maxiter": maxiter, "xatol": 1e-5, "fatol": 1e-6})
    gc.collect()
    return result


# ===================================================================
# GRAPE benchmark (reuse from main script, add gc)
# ===================================================================
def run_grape_benchmark(model, frame, logical_n, n0, theta, phi, T, cphase=True):
    from cqed_sim import (
        GrapeConfig, GrapeSolver, ModelControlChannelSpec,
        PiecewiseConstantTimeGrid, UnitaryObjective,
        build_control_problem_from_model,
    )
    from cqed_sim.unitary_synthesis import Subspace

    n_slices = 48
    dt_grape = T / n_slices
    amp_bound = 2 * np.pi * 50e6

    indices = []
    labels = []
    for n in range(int(logical_n)):
        indices.append(0 * int(model.n_cav) + n)
        indices.append(1 * int(model.n_cav) + n)
        labels.append(f"|g,{n}>")
        labels.append(f"|e,{n}>")
    sub = Subspace(
        full_dim=int(model.n_tr) * int(model.n_cav),
        indices=tuple(indices), labels=tuple(labels)
    )

    I2 = np.eye(2, dtype=np.complex128)
    R = target_qubit_unitary(theta, phi)
    tgt_blocks = [R if n == int(n0) else I2 for n in range(int(logical_n))]
    target_sub = block_diag(*tgt_blocks)

    obj_kwargs = dict(
        target_operator=target_sub, subspace=sub,
        ignore_global_phase=True,
    )
    if cphase:
        phase_blocks = tuple((2 * n, 2 * n + 1) for n in range(int(logical_n)))
        obj_kwargs["phase_blocks"] = phase_blocks

    problem = build_control_problem_from_model(
        model, frame=frame,
        time_grid=PiecewiseConstantTimeGrid.uniform(steps=n_slices, dt_s=dt_grape),
        channel_specs=(
            ModelControlChannelSpec(
                name="qubit_I", target="qubit", quadratures=("I",),
                amplitude_bounds=(-amp_bound, amp_bound),
            ),
            ModelControlChannelSpec(
                name="qubit_Q", target="qubit", quadratures=("Q",),
                amplitude_bounds=(-amp_bound, amp_bound),
            ),
        ),
        objectives=(UnitaryObjective(**obj_kwargs),),
    )

    config = GrapeConfig(maxiter=300, seed=42)
    result = GrapeSolver(config).solve(problem)

    if "nominal_fidelity" in result.metrics:
        fid = result.metrics["nominal_fidelity"]
    elif "fidelity" in result.metrics:
        fid = result.metrics["fidelity"]
    else:
        fid = 1.0 - result.objective_value if result.objective_value <= 1.0 else 0.0

    del problem, result
    gc.collect()
    return {"fidelity": float(fid)}


# ===================================================================
# Incremental save
# ===================================================================
def save_results(results):
    np.savez(str(OUTPATH), **results)
    print(f"  [Saved to {OUTPATH}]")


# ===================================================================
# Main continuation
# ===================================================================
def main():
    t_start = time.time()

    # Step 1: Assemble existing data
    print("Assembling existing Phase B + C results from terminal output...")
    results = assemble_existing_results()

    # Step 2: Build model
    print("Building model...")
    model = build_model()
    frame = build_frame(model)

    # Step 3: Fill detailed baseline metrics
    results = fill_baseline_details(results, model, frame)
    save_results(results)
    print(f"  Baselines done ({time.time() - t_start:.1f}s)\n")

    # Step 4: Complete Phase C missing points (detuned + smooth at χT/2π=5.0)
    print("=" * 60)
    print("Completing Phase C: χT/2π = 5.0 (detuned, smooth-basis)")
    print("=" * 60)
    T5 = duration_from_chi_t(5.0)

    # Detuned at χT/2π=5.0
    print("  Detuned multitone at χT/2π=5.0...", end=" ", flush=True)
    t0 = time.time()
    from run_followup_multitone import optimize_detuned_multitone
    amp_r2, ph_r2, det_r, res2 = optimize_detuned_multitone(
        model, frame, REP_LOGICAL_N, REP_N0, REP_THETA, REP_PHI, T5,
        objective="block_phase_relaxed_fid"
    )
    pulses, dops, T_tot = build_independent_tone_multitone(
        model, frame, REP_LOGICAL_N, REP_N0, REP_THETA, REP_PHI, T5,
        amp_ratios=amp_r2, phases=ph_r2, detunings=det_r
    )
    full_op, states = simulate_and_extract(model, frame, pulses, dops, REP_LOGICAL_N, T_tot)
    m = compute_all_metrics(full_op, states, model, REP_LOGICAL_N, REP_N0, REP_THETA, REP_PHI)
    results["opt_detuned_block_phase_relaxed_fid"][5] = m["block_phase_relaxed_fid"]
    results["opt_detuned_strict_logical_fid"][5] = m["strict_logical_fid"]
    for metric in ["target_branch_fid", "branch_true_mean", "branch_cphase_mean",
                    "same_block_pop_mean", "leakage_max",
                    "spectator_phase_spread", "spectator_max_transverse"]:
        results[f"opt_detuned_{metric}"][5] = m[metric]
    del full_op, states, m
    gc.collect()
    print(f"F_block={results['opt_detuned_block_phase_relaxed_fid'][5]:.4f} ({time.time()-t0:.1f}s)")

    # Smooth-basis at χT/2π=5.0
    print("  Smooth-basis at χT/2π=5.0...", end=" ", flush=True)
    t0 = time.time()
    from run_followup_multitone import optimize_smooth_basis
    sb_coeffs, sb_phases, res3 = optimize_smooth_basis(
        model, frame, REP_LOGICAL_N, REP_N0, REP_THETA, REP_PHI, T5,
        n_basis=3, objective="block_phase_relaxed_fid"
    )
    pulses, dops, T_tot = build_smooth_basis_multitone(
        model, frame, REP_LOGICAL_N, REP_N0, REP_THETA, REP_PHI, T5,
        n_basis=3, all_coeffs=sb_coeffs, all_phases=sb_phases
    )
    full_op, states = simulate_and_extract(model, frame, pulses, dops, REP_LOGICAL_N, T_tot)
    m = compute_all_metrics(full_op, states, model, REP_LOGICAL_N, REP_N0, REP_THETA, REP_PHI)
    results["opt_smooth_block_phase_relaxed_fid"][5] = m["block_phase_relaxed_fid"]
    results["opt_smooth_strict_logical_fid"][5] = m["strict_logical_fid"]
    for metric in ["target_branch_fid", "branch_true_mean", "branch_cphase_mean",
                    "same_block_pop_mean", "leakage_max",
                    "spectator_phase_spread", "spectator_max_transverse"]:
        results[f"opt_smooth_{metric}"][5] = m[metric]
    del full_op, states, m
    gc.collect()
    print(f"F_block={results['opt_smooth_block_phase_relaxed_fid'][5]:.4f} ({time.time()-t0:.1f}s)")
    save_results(results)

    # Step 5: Fill missing opt_indep detailed metrics (re-optimize at each χT)
    print("\n" + "=" * 60)
    print("Filling opt_indep detailed metrics")
    print("=" * 60)
    for ci, chi_t in enumerate(CHI_T_VALUES):
        T = duration_from_chi_t(chi_t)
        print(f"  χT/2π={chi_t:.1f}...", end=" ", flush=True)
        t0 = time.time()
        amp_r, ph_r, _ = optimize_independent_tone_safe(
            model, frame, REP_LOGICAL_N, REP_N0, REP_THETA, REP_PHI, T)
        pulses, dops, T_tot = build_independent_tone_multitone(
            model, frame, REP_LOGICAL_N, REP_N0, REP_THETA, REP_PHI, T,
            amp_ratios=amp_r, phases=ph_r)
        full_op, states = simulate_and_extract(
            model, frame, pulses, dops, REP_LOGICAL_N, T_tot)
        m = compute_all_metrics(
            full_op, states, model, REP_LOGICAL_N, REP_N0, REP_THETA, REP_PHI)
        for metric in ["target_branch_fid", "branch_true_mean", "branch_cphase_mean",
                        "strict_logical_fid", "block_phase_relaxed_fid",
                        "same_block_pop_mean", "leakage_max",
                        "spectator_phase_spread", "spectator_max_transverse"]:
            results[f"opt_indep_{metric}"][ci] = m[metric]
        del full_op, states, m
        gc.collect()
        print(f"F_block={results['opt_indep_block_phase_relaxed_fid'][ci]:.4f} ({time.time()-t0:.1f}s)")
    save_results(results)

    # Step 6: Phase D — 2-segment multitone
    print("\n" + "=" * 60)
    print("PHASE D: Segmented multitone optimization")
    print("=" * 60)
    n_chi = len(CHI_T_VALUES)
    for metric in ["target_branch_fid", "branch_true_mean", "branch_cphase_mean",
                    "strict_logical_fid", "block_phase_relaxed_fid",
                    "same_block_pop_mean", "leakage_max",
                    "spectator_phase_spread", "spectator_max_transverse"]:
        results[f"opt_2seg_{metric}"] = np.full(n_chi, np.nan)
    results["opt_2seg_opt_cost"] = np.full(n_chi, np.nan)

    for ci, chi_t in enumerate(CHI_T_VALUES):
        T = duration_from_chi_t(chi_t)
        print(f"  χT/2π={chi_t:.1f} 2-segment...", end=" ", flush=True)
        t0 = time.time()
        res4 = optimize_2seg_safe(
            model, frame, REP_LOGICAL_N, REP_N0, REP_THETA, REP_PHI, T)

        # Reconstruct and evaluate
        n_tones = REP_LOGICAL_N
        non_target = [n for n in range(n_tones) if n != REP_N0]
        n_free = len(non_target)
        n_per_seg = n_free + n_tones
        seg_params = []
        for seg in range(2):
            offset = seg * n_per_seg
            amp_ratios = np.zeros(n_tones)
            phases_arr = np.zeros(n_tones)
            amp_ratios[REP_N0] = 1.0
            for i, n in enumerate(non_target):
                amp_ratios[n] = res4.x[offset + i]
            for n in range(n_tones):
                phases_arr[n] = res4.x[offset + n_free + n]
            seg_params.append({
                "amp_ratios": amp_ratios, "phases": phases_arr, "detunings": None})
        pulses, dops, T_tot = build_segmented_multitone(
            model, frame, REP_LOGICAL_N, REP_N0, REP_THETA, REP_PHI, T,
            n_segments=2, segment_params_list=seg_params)
        full_op, states = simulate_and_extract(
            model, frame, pulses, dops, REP_LOGICAL_N, T_tot)
        m = compute_all_metrics(
            full_op, states, model, REP_LOGICAL_N, REP_N0, REP_THETA, REP_PHI)
        for metric in ["target_branch_fid", "branch_true_mean", "branch_cphase_mean",
                        "strict_logical_fid", "block_phase_relaxed_fid",
                        "same_block_pop_mean", "leakage_max",
                        "spectator_phase_spread", "spectator_max_transverse"]:
            results[f"opt_2seg_{metric}"][ci] = m[metric]
        results["opt_2seg_opt_cost"][ci] = res4.fun
        del full_op, states, m, res4
        gc.collect()
        print(f"F_block={results['opt_2seg_block_phase_relaxed_fid'][ci]:.4f} "
              f"F_strict={results['opt_2seg_strict_logical_fid'][ci]:.4f} ({time.time()-t0:.1f}s)")
    save_results(results)
    print(f"Phase D complete ({time.time() - t_start:.1f}s)\n")

    # Step 7: Phase E — GRAPE comparison
    print("=" * 60)
    print("PHASE E: GRAPE comparison")
    print("=" * 60)
    grape_chi_t = np.array([1.0, 2.0, 3.0, 5.0])
    results["grape_chi_t"] = grape_chi_t
    results["grape_cphase_fid"] = np.zeros(len(grape_chi_t))
    results["grape_true_fid"] = np.zeros(len(grape_chi_t))

    for ci, chi_t in enumerate(grape_chi_t):
        T = duration_from_chi_t(chi_t)
        print(f"  χT/2π={chi_t:.1f}", end=" ", flush=True)
        t0 = time.time()

        r_cp = run_grape_benchmark(model, frame, REP_LOGICAL_N, REP_N0,
                                    REP_THETA, REP_PHI, T, cphase=True)
        results["grape_cphase_fid"][ci] = r_cp["fidelity"]
        print(f"cphase={r_cp['fidelity']:.6f}", end=" ", flush=True)

        r_tr = run_grape_benchmark(model, frame, REP_LOGICAL_N, REP_N0,
                                    REP_THETA, REP_PHI, T, cphase=False)
        results["grape_true_fid"][ci] = r_tr["fidelity"]
        print(f"true={r_tr['fidelity']:.6f} ({time.time()-t0:.1f}s)")
        save_results(results)

    print(f"Phase E complete ({time.time() - t_start:.1f}s)\n")

    # Step 8: Phase C-EXT — Branch scan
    print("=" * 60)
    print("PHASE C-EXT: Branch and angle scans")
    print("=" * 60)
    scan_branches = TARGET_BRANCHES
    scan_chi_t = np.array([1.0, 2.0, 3.0, 5.0])
    results["scan_branches"] = np.array(scan_branches)
    results["scan_chi_t"] = scan_chi_t

    for prefix in ["baseline_gauss", "opt_indep"]:
        for metric in ["strict_logical_fid", "block_phase_relaxed_fid",
                        "branch_cphase_mean", "leakage_max"]:
            results[f"branch_scan_{prefix}_{metric}"] = np.zeros(
                (len(scan_branches), len(scan_chi_t)))

    for bi, n0 in enumerate(scan_branches):
        for ci, chi_t in enumerate(scan_chi_t):
            T = duration_from_chi_t(chi_t)
            ln = min(REP_LOGICAL_N, max(n0 + 2, 3))

            # Baseline
            pulses, dops, T_tot = build_single_tone_gaussian(
                model, frame, n0, REP_THETA, REP_PHI, T)
            full_op, states = simulate_and_extract(
                model, frame, pulses, dops, ln, T_tot)
            m = compute_all_metrics(full_op, states, model, ln, n0, REP_THETA, REP_PHI)
            for metric in ["strict_logical_fid", "block_phase_relaxed_fid",
                            "branch_cphase_mean", "leakage_max"]:
                results[f"branch_scan_baseline_gauss_{metric}"][bi, ci] = m[metric]

            # Optimized
            amp_r, ph_r, _ = optimize_independent_tone_safe(
                model, frame, ln, n0, REP_THETA, REP_PHI, T, maxiter=100)
            pulses, dops, T_tot = build_independent_tone_multitone(
                model, frame, ln, n0, REP_THETA, REP_PHI, T,
                amp_ratios=amp_r, phases=ph_r)
            full_op, states = simulate_and_extract(
                model, frame, pulses, dops, ln, T_tot)
            m = compute_all_metrics(full_op, states, model, ln, n0, REP_THETA, REP_PHI)
            for metric in ["strict_logical_fid", "block_phase_relaxed_fid",
                            "branch_cphase_mean", "leakage_max"]:
                results[f"branch_scan_opt_indep_{metric}"][bi, ci] = m[metric]

            del full_op, states, m
            gc.collect()
            print(f"  n0={n0} χT/2π={chi_t:.1f}  "
                  f"gauss_F={results[f'branch_scan_baseline_gauss_block_phase_relaxed_fid'][bi,ci]:.3f}  "
                  f"opt_F={results[f'branch_scan_opt_indep_block_phase_relaxed_fid'][bi,ci]:.3f}")
        save_results(results)

    # Angle scan
    scan_thetas = THETA_VALUES
    scan_phis = PHI_VALUES
    results["scan_thetas"] = np.array(scan_thetas)
    results["scan_phis"] = np.array(scan_phis)

    for prefix in ["baseline_gauss", "opt_indep"]:
        for metric in ["strict_logical_fid", "block_phase_relaxed_fid"]:
            results[f"angle_scan_{prefix}_{metric}"] = np.zeros(
                (len(scan_thetas), len(scan_phis), len(scan_chi_t)))

    for ti, theta in enumerate(scan_thetas):
        for pi_idx, phi in enumerate(scan_phis):
            for ci, chi_t in enumerate(scan_chi_t):
                T = duration_from_chi_t(chi_t)
                # Baseline
                pulses, dops, T_tot = build_single_tone_gaussian(
                    model, frame, REP_N0, theta, phi, T)
                full_op, states = simulate_and_extract(
                    model, frame, pulses, dops, REP_LOGICAL_N, T_tot)
                m = compute_all_metrics(
                    full_op, states, model, REP_LOGICAL_N, REP_N0, theta, phi)
                for metric in ["strict_logical_fid", "block_phase_relaxed_fid"]:
                    results[f"angle_scan_baseline_gauss_{metric}"][ti, pi_idx, ci] = m[metric]

                # Optimized
                amp_r, ph_r, _ = optimize_independent_tone_safe(
                    model, frame, REP_LOGICAL_N, REP_N0, theta, phi, T, maxiter=100)
                pulses, dops, T_tot = build_independent_tone_multitone(
                    model, frame, REP_LOGICAL_N, REP_N0, theta, phi, T,
                    amp_ratios=amp_r, phases=ph_r)
                full_op, states = simulate_and_extract(
                    model, frame, pulses, dops, REP_LOGICAL_N, T_tot)
                m = compute_all_metrics(
                    full_op, states, model, REP_LOGICAL_N, REP_N0, theta, phi)
                for metric in ["strict_logical_fid", "block_phase_relaxed_fid"]:
                    results[f"angle_scan_opt_indep_{metric}"][ti, pi_idx, ci] = m[metric]

                del full_op, states, m
                gc.collect()
                print(f"  θ={theta/np.pi:.2f}π φ={phi/np.pi:.2f}π χT/2π={chi_t:.1f}  "
                      f"gauss={results[f'angle_scan_baseline_gauss_block_phase_relaxed_fid'][ti,pi_idx,ci]:.3f}  "
                      f"opt={results[f'angle_scan_opt_indep_block_phase_relaxed_fid'][ti,pi_idx,ci]:.3f}")
        save_results(results)

    print(f"\nAll phases complete. Total time: {time.time() - t_start:.1f}s")
    print(f"Results saved to {OUTPATH}")


if __name__ == "__main__":
    main()
