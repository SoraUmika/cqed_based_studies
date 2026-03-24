"""
Validation harness for the gray-box adaptive control study.

This script validates that the consolidated study is reproducible, internally
consistent, and numerically stable at the key 30% mismatch operating point.
It performs two kinds of checks:

1. Static checks on the archived NPZ data and generated figures.
2. Fresh cqed_sim spot checks for convergence with respect to truncation and
   multistart optimization choices.

The summary is written to data/validation_summary.json.

Usage (from the study directory):
    python scripts/validate_results.py
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np

STUDY_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = STUDY_DIR / "data"
FIGURES_DIR = STUDY_DIR / "figures"
OUTPUT_PATH = DATA_DIR / "validation_summary.json"

sys.path.insert(0, str(STUDY_DIR / "scripts"))

from cqed_sim import DispersiveTransmonCavityModel, FrameSpec  # noqa: E402
from cqed_sim.unitary_synthesis import Subspace  # noqa: E402
from control import eval_on_model, run_grape_multistart  # noqa: E402
from models import (  # noqa: E402
    ALPHA,
    CHI_HIGHER_TRUE,
    CHI_TRUE,
    KERR_TRUE,
    OMEGA_C,
    OMEGA_Q,
    TRUTH_NOISE,
    CONFUSION_MATRIX,
    N_SHOTS_PROBE,
    make_frame,
    make_target_matrix,
    make_truth_model,
)
from probe_library import (  # noqa: E402
    infer_chi_from_probe,
    make_ramsey_delays,
    run_chi_ramsey_probe,
    t2_star_from_noise,
)

REQUIRED_DATA_FILES = [
    "phase4_results.npz",
    "phase5_1_chi_mismatch.npz",
    "phase5_2_noise_sweep.npz",
    "phase5_3_readout_sweep.npz",
    "phase5_4_probe_budget.npz",
    "phase5_5_drift.npz",
    "phase5_6_omission.npz",
]

REQUIRED_FIGURE_FILES = [
    "phase4_main_comparison.pdf",
    "phase4_per_fock.pdf",
    "phase5_1_chi_mismatch.pdf",
    "phase5_2_noise.pdf",
    "phase5_3_readout.pdf",
    "phase5_4_probe_budget.pdf",
    "phase5_5_drift.pdf",
    "phase5_6_omission.pdf",
    "viability_summary.pdf",
]

THRESHOLDS = {
    "phase4_gap_close_30pct": 2.0e-4,
    "phase4_gray_advantage_30pct": 5.0e-2,
    "phase5_1_gray_range": 5.0e-4,
    "phase5_2_min_advantage": 5.0e-2,
    "phase5_3_gray_range": 5.0e-4,
    "phase5_3_max_chi_error": 3.0e-4,
    "phase5_4_sigma_scaling_rel_spread": 2.0e-2,
    "phase5_5_improvement_0p1pct": 3.0e-2,
    "phase5_5_improvement_1pct": 1.0e-2,
    "phase5_6_max_penalty_with_correct_chi": 3.0e-3,
    "convergence_truncation_delta": 5.0e-4,
    "convergence_multistart_delta": 1.5e-3,
}


def _configure_utf8_output() -> None:
    for stream_name in ("stdout", "stderr"):
        stream = getattr(sys, stream_name, None)
        if stream is not None and hasattr(stream, "reconfigure"):
            stream.reconfigure(encoding="utf-8")


def _load_npz(name: str) -> dict[str, Any]:
    path = DATA_DIR / name
    with np.load(path, allow_pickle=True) as data:
        return {key: data[key] for key in data.files}


def _to_jsonable(value: Any) -> Any:
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {k: _to_jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_to_jsonable(v) for v in value]
    return value


def _record_check(
    checks: list[dict[str, Any]],
    name: str,
    passed: bool,
    details: dict[str, Any],
) -> None:
    checks.append(
        {
            "name": name,
            "passed": bool(passed),
            "details": _to_jsonable(details),
        }
    )


def _make_truth_model_for_dims(n_cav: int, n_tr: int) -> DispersiveTransmonCavityModel:
    return DispersiveTransmonCavityModel(
        omega_c=OMEGA_C,
        omega_q=OMEGA_Q,
        alpha=ALPHA,
        chi=CHI_TRUE,
        chi_higher=(CHI_HIGHER_TRUE,),
        kerr=KERR_TRUE,
        n_cav=n_cav,
        n_tr=n_tr,
    )


def _make_learner_model_for_dims(n_cav: int, n_tr: int, chi: float) -> DispersiveTransmonCavityModel:
    return DispersiveTransmonCavityModel(
        omega_c=OMEGA_C,
        omega_q=OMEGA_Q,
        alpha=ALPHA,
        chi=chi,
        chi_higher=(),
        kerr=0.0,
        n_cav=n_cav,
        n_tr=n_tr,
    )


def _make_logical_subspace(n_cav: int, n_tr: int) -> Subspace:
    full_dim = n_cav * n_tr
    indices = list(range(4)) + list(range(n_cav, n_cav + 4))
    labels = [f"|g,{n}>" for n in range(4)] + [f"|e,{n}>" for n in range(4)]
    return Subspace.custom(full_dim=full_dim, indices=indices, labels=labels)


def _run_perfect_spot_check(
    n_cav: int,
    n_tr: int,
    seeds: list[int],
    maxiter: int,
) -> dict[str, float]:
    learner = _make_learner_model_for_dims(n_cav=n_cav, n_tr=n_tr, chi=CHI_TRUE)
    truth = _make_truth_model_for_dims(n_cav=n_cav, n_tr=n_tr)
    frame = FrameSpec(
        omega_c_frame=float(learner.omega_c),
        omega_q_frame=float(learner.omega_q),
    )
    subspace = _make_logical_subspace(n_cav=n_cav, n_tr=n_tr)
    target = make_target_matrix()
    result, problem, train_fidelity, best_seed = run_grape_multistart(
        learner,
        frame,
        subspace,
        target,
        seeds=seeds,
        n_steps=16,
        dt_s=10e-9,
        maxiter=maxiter,
    )
    eval_fidelity, _ = eval_on_model(result, problem, truth, frame, eval_noise=None)
    return {
        "train_fidelity": float(train_fidelity),
        "eval_fidelity": float(eval_fidelity),
        "best_seed": int(best_seed),
    }


def _run_static_checks(checks: list[dict[str, Any]]) -> None:
    missing_data = [name for name in REQUIRED_DATA_FILES if not (DATA_DIR / name).exists()]
    missing_figures = [name for name in REQUIRED_FIGURE_FILES if not (FIGURES_DIR / name).exists()]
    _record_check(
        checks,
        "required_artifacts_present",
        not missing_data and not missing_figures,
        {
            "missing_data_files": missing_data,
            "missing_figure_files": missing_figures,
        },
    )

    phase4 = _load_npz("phase4_results.npz")
    idx_30 = int(np.where(np.isclose(phase4["mismatch_fractions"], 0.30))[0][0])
    gray_gap = abs(
        float(phase4["gray_box_fidelities"][idx_30])
        - float(phase4["perfect_fidelities"][idx_30])
    )
    gray_advantage = float(phase4["gray_box_fidelities"][idx_30] - phase4["nominal_fidelities"][idx_30])
    _record_check(
        checks,
        "phase4_gray_box_closes_perfect_gap_at_30pct",
        gray_gap < THRESHOLDS["phase4_gap_close_30pct"]
        and gray_advantage > THRESHOLDS["phase4_gray_advantage_30pct"],
        {
            "gray_box_30pct": float(phase4["gray_box_fidelities"][idx_30]),
            "perfect_30pct": float(phase4["perfect_fidelities"][idx_30]),
            "nominal_30pct": float(phase4["nominal_fidelities"][idx_30]),
            "gap_to_perfect": gray_gap,
            "advantage_over_nominal": gray_advantage,
        },
    )

    bb_fidelity = float(phase4["bb_fidelity"])
    _record_check(
        checks,
        "phase4_black_box_underperforms_gray_box",
        bb_fidelity < float(phase4["gray_box_fidelities"][idx_30])
        and bb_fidelity < float(phase4["nominal_fidelities"][idx_30]),
        {
            "black_box_fidelity": bb_fidelity,
            "gray_box_30pct": float(phase4["gray_box_fidelities"][idx_30]),
            "nominal_30pct": float(phase4["nominal_fidelities"][idx_30]),
            "black_box_evaluations": int(phase4["bb_n_evaluations"]),
        },
    )

    phase51 = _load_npz("phase5_1_chi_mismatch.npz")
    gray_range = float(np.max(phase51["gray_box_fidelities"]) - np.min(phase51["gray_box_fidelities"]))
    _record_check(
        checks,
        "phase5_1_gray_box_remains_flat_across_large_mismatch_range",
        gray_range < THRESHOLDS["phase5_1_gray_range"],
        {
            "gray_box_range": gray_range,
            "gray_box_min": float(np.min(phase51["gray_box_fidelities"])),
            "gray_box_max": float(np.max(phase51["gray_box_fidelities"])),
        },
    )

    phase52 = _load_npz("phase5_2_noise_sweep.npz")
    gray_advantages = phase52["gray_box_fidelities"] - phase52["nominal_fidelities"]
    _record_check(
        checks,
        "phase5_2_gray_box_advantage_survives_noise_sweep",
        bool(np.all(gray_advantages > THRESHOLDS["phase5_2_min_advantage"])),
        {
            "gray_advantages": gray_advantages,
            "min_advantage": float(np.min(gray_advantages)),
        },
    )

    phase53 = _load_npz("phase5_3_readout_sweep.npz")
    gray_range_readout = float(np.max(phase53["gray_box_fidelities"]) - np.min(phase53["gray_box_fidelities"]))
    max_chi_error = float(np.max(phase53["chi_errors"]))
    _record_check(
        checks,
        "phase5_3_readout_robustness",
        gray_range_readout < THRESHOLDS["phase5_3_gray_range"]
        and max_chi_error < THRESHOLDS["phase5_3_max_chi_error"],
        {
            "gray_box_fidelity_range": gray_range_readout,
            "max_fractional_chi_error": max_chi_error,
        },
    )

    phase54 = _load_npz("phase5_4_probe_budget.npz")
    scaled_sigmas = phase54["chi_uncertainties"] * np.sqrt(phase54["n_shots_values"])
    rel_spread = float(np.std(scaled_sigmas) / np.mean(scaled_sigmas))
    _record_check(
        checks,
        "phase5_4_shot_noise_scaling",
        rel_spread < THRESHOLDS["phase5_4_sigma_scaling_rel_spread"],
        {
            "scaled_sigma_values": scaled_sigmas,
            "relative_spread": rel_spread,
        },
    )

    phase55 = _load_npz("phase5_5_drift.npz")
    drift_rates = phase55["drift_rates"]
    recal_labels = [str(label) for label in phase55["recal_labels"]]
    never_idx = recal_labels.index("never")
    every5_idx = recal_labels.index("every_5")
    idx_0p1 = int(np.where(np.isclose(drift_rates, 0.001))[0][0])
    idx_1p0 = int(np.where(np.isclose(drift_rates, 0.01))[0][0])
    final_never_0p1 = float(phase55["fidelity_traces"][idx_0p1, never_idx, -1])
    final_every5_0p1 = float(phase55["fidelity_traces"][idx_0p1, every5_idx, -1])
    final_never_1p0 = float(phase55["fidelity_traces"][idx_1p0, never_idx, -1])
    final_every5_1p0 = float(phase55["fidelity_traces"][idx_1p0, every5_idx, -1])
    _record_check(
        checks,
        "phase5_5_recalibration_improves_drifted_operation",
        (final_every5_0p1 - final_never_0p1) > THRESHOLDS["phase5_5_improvement_0p1pct"]
        and (final_every5_1p0 - final_never_1p0) > THRESHOLDS["phase5_5_improvement_1pct"],
        {
            "final_0p1pct_never": final_never_0p1,
            "final_0p1pct_every5": final_every5_0p1,
            "improvement_0p1pct": final_every5_0p1 - final_never_0p1,
            "final_1pct_never": final_never_1p0,
            "final_1pct_every5": final_every5_1p0,
            "improvement_1pct": final_every5_1p0 - final_never_1p0,
        },
    )

    phase56 = _load_npz("phase5_6_omission.npz")
    true_idx = int(np.where(np.isclose(phase56["chi_higher_multipliers"], 1.0))[0][0])
    correct_penalty = abs(
        float(phase56["fidelities_chi_correct"][0])
        - float(phase56["fidelities_chi_correct"][true_idx])
    )
    chi_penalty = float(
        phase56["fidelities_chi_correct"][true_idx]
        - phase56["fidelities_chi_wrong"][true_idx]
    )
    _record_check(
        checks,
        "phase5_6_wrong_chi_dominates_wrong_chi_higher",
        correct_penalty < THRESHOLDS["phase5_6_max_penalty_with_correct_chi"]
        and chi_penalty > THRESHOLDS["phase4_gray_advantage_30pct"],
        {
            "penalty_from_omitting_chi_higher_with_correct_chi": correct_penalty,
            "penalty_from_wrong_chi_even_with_correct_chi_higher": chi_penalty,
        },
    )


def _run_convergence_checks(checks: list[dict[str, Any]]) -> None:
    baseline = _run_perfect_spot_check(n_cav=4, n_tr=3, seeds=[2, 9], maxiter=120)
    expanded_cavity = _run_perfect_spot_check(n_cav=5, n_tr=3, seeds=[2, 9], maxiter=120)
    expanded_transmon = _run_perfect_spot_check(n_cav=4, n_tr=4, seeds=[2, 9], maxiter=120)

    cavity_delta = abs(expanded_cavity["eval_fidelity"] - baseline["eval_fidelity"])
    transmon_delta = abs(expanded_transmon["eval_fidelity"] - baseline["eval_fidelity"])
    _record_check(
        checks,
        "convergence_truncation_stable_at_30pct_operating_point",
        cavity_delta < THRESHOLDS["convergence_truncation_delta"]
        and transmon_delta < THRESHOLDS["convergence_truncation_delta"],
        {
            "baseline_eval_fidelity": baseline["eval_fidelity"],
            "expanded_cavity_eval_fidelity": expanded_cavity["eval_fidelity"],
            "expanded_transmon_eval_fidelity": expanded_transmon["eval_fidelity"],
            "delta_n_cav": cavity_delta,
            "delta_n_tr": transmon_delta,
        },
    )

    default_multistart = _run_perfect_spot_check(n_cav=4, n_tr=3, seeds=[2, 9, 14], maxiter=200)
    expanded_multistart = _run_perfect_spot_check(
        n_cav=4,
        n_tr=3,
        seeds=[2, 9, 14, 21, 33],
        maxiter=200,
    )
    multistart_delta = abs(expanded_multistart["eval_fidelity"] - default_multistart["eval_fidelity"])
    _record_check(
        checks,
        "convergence_multistart_seed_sensitivity_small",
        multistart_delta < THRESHOLDS["convergence_multistart_delta"],
        {
            "default_multistart_eval_fidelity": default_multistart["eval_fidelity"],
            "expanded_multistart_eval_fidelity": expanded_multistart["eval_fidelity"],
            "delta": multistart_delta,
            "default_best_seed": default_multistart["best_seed"],
            "expanded_best_seed": expanded_multistart["best_seed"],
        },
    )


def _run_probe_feedback_spot_check(checks: list[dict[str, Any]]) -> None:
    truth = make_truth_model()
    truth_frame = make_frame(truth)
    subspace = _make_logical_subspace(n_cav=4, n_tr=3)
    target = make_target_matrix()
    chi_prior = CHI_TRUE * 0.7
    rng = np.random.default_rng(12345)
    t2_star = t2_star_from_noise(TRUTH_NOISE)
    delays = make_ramsey_delays(chi_prior, n_periods=3.0, n_points=80)
    probe_data = run_chi_ramsey_probe(
        chi_true=truth.chi,
        chi_higher_true=float(truth.chi_higher[0]),
        t2_star=t2_star,
        confusion_matrix=CONFUSION_MATRIX,
        n_shots=N_SHOTS_PROBE,
        delays_s=delays,
        fock_levels=[1, 2, 3],
        rng=rng,
    )
    infer_result = infer_chi_from_probe(
        probe_data=probe_data,
        confusion_matrix=CONFUSION_MATRIX,
        n_shots=N_SHOTS_PROBE,
        chi_initial=chi_prior,
        chi_higher_initial=0.0,
    )

    def evaluate_feedback(chi_higher_value: float) -> float:
        learner = DispersiveTransmonCavityModel(
            omega_c=OMEGA_C,
            omega_q=OMEGA_Q,
            alpha=ALPHA,
            chi=float(infer_result["chi"]),
            chi_higher=(float(chi_higher_value),) if chi_higher_value != 0.0 else (),
            kerr=0.0,
            n_cav=4,
            n_tr=3,
        )
        frame = FrameSpec(
            omega_c_frame=float(learner.omega_c),
            omega_q_frame=float(learner.omega_q),
        )
        result, problem, _, _ = run_grape_multistart(
            learner,
            frame,
            subspace,
            target,
            seeds=[2, 9, 14],
            n_steps=16,
            dt_s=10e-9,
            maxiter=200,
        )
        fidelity, _ = eval_on_model(result, problem, truth, truth_frame, eval_noise=None)
        return fidelity

    fidelity_chi_only = evaluate_feedback(0.0)
    fidelity_with_inferred_chi_higher = evaluate_feedback(float(infer_result["chi_higher"]))
    _record_check(
        checks,
        "probe_feedback_spot_check_partial_chi_higher_feedback_not_adopted",
        fidelity_chi_only >= fidelity_with_inferred_chi_higher,
        {
            "fidelity_chi_only_feedback": fidelity_chi_only,
            "fidelity_with_inferred_chi_higher_feedback": fidelity_with_inferred_chi_higher,
            "inferred_chi_higher_khz": float(infer_result["chi_higher"] / (2.0 * np.pi) / 1.0e3),
        },
    )


def main() -> None:
    _configure_utf8_output()
    start = time.perf_counter()
    checks: list[dict[str, Any]] = []

    print("Running archived-data checks...")
    _run_static_checks(checks)

    print("Running fresh convergence checks...")
    _run_convergence_checks(checks)

    print("Running probe-feedback spot check...")
    _run_probe_feedback_spot_check(checks)

    passed = sum(1 for check in checks if check["passed"])
    failed = len(checks) - passed
    summary = {
        "study": "gray_box_adaptive_control",
        "generated_at": time.strftime("%Y-%m-%d %H:%M:%S"),
        "wall_time_s": time.perf_counter() - start,
        "all_passed": failed == 0,
        "passed_checks": passed,
        "failed_checks": failed,
        "checks": checks,
    }

    OUTPUT_PATH.write_text(json.dumps(_to_jsonable(summary), indent=2), encoding="utf-8")

    print()
    for check in checks:
        status = "PASS" if check["passed"] else "FAIL"
        print(f"[{status}] {check['name']}")
    print()
    print(f"Wrote {OUTPUT_PATH}")
    print(f"Summary: {passed} passed, {failed} failed")

    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
