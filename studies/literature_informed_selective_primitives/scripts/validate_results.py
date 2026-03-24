"""Validation for the literature-informed selective-primitives study."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from common import (
    DATA_DIR,
    DT,
    FIGURES_DIR,
    LOGICAL_N,
    SNAP_PHASE,
    SQR_PHI,
    SQR_THETA,
    TARGET_BRANCH,
    build_frame,
    build_model,
    build_noise_spec,
    build_session,
    build_snap_pulses,
    build_sqr_pulses,
    extract_restricted_operator,
    logical_indices,
    sqr_strict_target_operator,
    snap_target_operator,
    state_fidelity_to_pure_target,
)
from run_study import evaluate_snap_candidate, evaluate_sqr_candidate


RESULTS_PATH = DATA_DIR / "study_results.json"
VALIDATION_PATH = DATA_DIR / "validation_summary.json"


def _strict_process_fidelity(target: np.ndarray, actual: np.ndarray) -> float:
    dim = float(target.shape[0])
    return float(np.clip(abs(np.trace(target.conj().T @ actual)) ** 2 / (dim * dim), 0.0, 1.0))


def _load_results() -> dict[str, object]:
    return json.loads(RESULTS_PATH.read_text(encoding="utf-8"))


def _record(checks: list[dict[str, object]], name: str, passed: bool, details: dict[str, object]) -> None:
    checks.append({"name": name, "passed": bool(passed), "details": details})
    status = "PASS" if passed else "FAIL"
    print(f"{status}: {name}")
    for key, value in details.items():
        print(f"  {key}: {value}")


def _reevaluate_best_sqr(best: dict[str, object], *, n_cav: int | None = None, dt: float = DT) -> dict[str, object]:
    model = build_model(n_cav=(n_cav if n_cav is not None else build_model().n_cav))
    frame = build_frame(model)
    noise = build_noise_spec()
    return evaluate_sqr_candidate(
        model,
        frame,
        noise,
        family=str(best["family"]),
        chi_t=float(best["chi_t"]),
        shape_parameter=float(best["shape_parameter"]),
        amplitude_scale=float(best["amplitude_scale"]),
        dt=dt,
    )


def _reevaluate_best_snap(best: dict[str, object], *, n_cav: int | None = None, dt: float = DT) -> dict[str, object]:
    model = build_model(n_cav=(n_cav if n_cav is not None else build_model().n_cav))
    frame = build_frame(model)
    noise = build_noise_spec()
    return evaluate_snap_candidate(
        model,
        frame,
        noise,
        family=str(best["family"]),
        chi_t=float(best["chi_t"]),
        shape_parameter=float(best["shape_parameter"]),
        amplitude_scale=float(best["amplitude_scale"]),
        dt=dt,
    )


def _zero_angle_sqr_basis_preservation(best: dict[str, object]) -> float:
    model = build_model()
    frame = build_frame(model)
    pulses, drive_ops, total_duration = build_sqr_pulses(
        model,
        frame,
        family=str(best["family"]),
        branch=TARGET_BRANCH,
        theta=0.0,
        phi=SQR_PHI,
        duration_s=float(best["duration_us"]) * 1.0e-6,
        shape_parameter=float(best["shape_parameter"]),
        amplitude_scale=float(best["amplitude_scale"]),
    )
    session = build_session(model, frame, pulses, drive_ops, total_duration_s=total_duration, noise=None)
    basis_states = [
        model.basis_state(qubit_level, storage_level)
        for storage_level in range(LOGICAL_N)
        for qubit_level in (0, 1)
    ]
    finals = [item.final_state for item in session.run_many(basis_states)]
    return float(np.mean([state_fidelity_to_pure_target(psi0, psi1) for psi0, psi1 in zip(basis_states, finals, strict=True)]))


def _zero_phase_snap_basis_preservation(best: dict[str, object]) -> float:
    model = build_model()
    frame = build_frame(model)
    pulses, drive_ops, total_duration = build_snap_pulses(
        model,
        frame,
        family=str(best["family"]),
        branch=TARGET_BRANCH,
        phase_angle=0.0,
        duration_s=float(best["pi_pulse_duration_us"]) * 1.0e-6,
        shape_parameter=float(best["shape_parameter"]),
        amplitude_scale=float(best["amplitude_scale"]),
    )
    session = build_session(model, frame, pulses, drive_ops, total_duration_s=total_duration, noise=None)
    basis_states = [model.basis_state(0, storage_level) for storage_level in range(LOGICAL_N)]
    finals = [item.final_state for item in session.run_many(basis_states)]
    return float(np.mean([state_fidelity_to_pure_target(psi0, psi1) for psi0, psi1 in zip(basis_states, finals, strict=True)]))


def main() -> None:
    results = _load_results()
    checks: list[dict[str, object]] = []

    sqr_best = dict(results["sqr"]["best_noisy_overall"])
    snap_best = dict(results["snap"]["best_noisy_overall"])

    required_artifacts = [
        RESULTS_PATH,
        FIGURES_DIR / "sqr_family_scan.png",
        FIGURES_DIR / "sqr_family_scan.pdf",
        FIGURES_DIR / "snap_family_scan.png",
        FIGURES_DIR / "snap_family_scan.pdf",
        FIGURES_DIR / "best_waveforms_and_summary.png",
        FIGURES_DIR / "best_waveforms_and_summary.pdf",
    ]
    artifact_ok = all(path.exists() and path.stat().st_size > 0 for path in required_artifacts)
    _record(
        checks,
        "artifacts_present",
        artifact_ok,
        {path.name: (path.exists() and path.stat().st_size > 0) for path in required_artifacts},
    )

    _record(
        checks,
        "performance_thresholds",
        (
            float(sqr_best["closed_relaxed_branch_mean"]) > 0.99
            and float(sqr_best["noisy_relaxed_avg_state_fidelity"]) > 0.98
            and float(snap_best["closed_process_fidelity"]) > 0.90
            and float(snap_best["noisy_avg_state_fidelity"]) > 0.92
        ),
        {
            "sqr_closed_relaxed": sqr_best["closed_relaxed_branch_mean"],
            "sqr_noisy_avg": sqr_best["noisy_relaxed_avg_state_fidelity"],
            "snap_closed_process": snap_best["closed_process_fidelity"],
            "snap_noisy_avg": snap_best["noisy_avg_state_fidelity"],
        },
    )

    sqr_dt_fine = _reevaluate_best_sqr(sqr_best, dt=1.0e-9)
    snap_dt_fine = _reevaluate_best_snap(snap_best, dt=1.0e-9)
    _record(
        checks,
        "time_step_convergence",
        (
            abs(float(sqr_dt_fine["closed_relaxed_branch_mean"]) - float(sqr_best["closed_relaxed_branch_mean"])) < 5.0e-3
            and abs(float(sqr_dt_fine["noisy_relaxed_avg_state_fidelity"]) - float(sqr_best["noisy_relaxed_avg_state_fidelity"])) < 5.0e-3
            and abs(float(snap_dt_fine["closed_process_fidelity"]) - float(snap_best["closed_process_fidelity"])) < 1.5e-2
            and abs(float(snap_dt_fine["noisy_avg_state_fidelity"]) - float(snap_best["noisy_avg_state_fidelity"])) < 1.5e-2
        ),
        {
            "sqr_closed_delta": float(sqr_dt_fine["closed_relaxed_branch_mean"]) - float(sqr_best["closed_relaxed_branch_mean"]),
            "sqr_noisy_delta": float(sqr_dt_fine["noisy_relaxed_avg_state_fidelity"]) - float(sqr_best["noisy_relaxed_avg_state_fidelity"]),
            "snap_closed_delta": float(snap_dt_fine["closed_process_fidelity"]) - float(snap_best["closed_process_fidelity"]),
            "snap_noisy_delta": float(snap_dt_fine["noisy_avg_state_fidelity"]) - float(snap_best["noisy_avg_state_fidelity"]),
        },
    )

    sqr_big = _reevaluate_best_sqr(sqr_best, n_cav=build_model().n_cav + 1)
    snap_big = _reevaluate_best_snap(snap_best, n_cav=build_model().n_cav + 1)
    _record(
        checks,
        "truncation_convergence",
        (
            abs(float(sqr_big["closed_relaxed_branch_mean"]) - float(sqr_best["closed_relaxed_branch_mean"])) < 2.0e-2
            and abs(float(sqr_big["noisy_relaxed_avg_state_fidelity"]) - float(sqr_best["noisy_relaxed_avg_state_fidelity"])) < 2.0e-2
            and abs(float(snap_big["closed_process_fidelity"]) - float(snap_best["closed_process_fidelity"])) < 3.0e-2
            and abs(float(snap_big["noisy_avg_state_fidelity"]) - float(snap_best["noisy_avg_state_fidelity"])) < 3.0e-2
        ),
        {
            "sqr_closed_delta": float(sqr_big["closed_relaxed_branch_mean"]) - float(sqr_best["closed_relaxed_branch_mean"]),
            "sqr_noisy_delta": float(sqr_big["noisy_relaxed_avg_state_fidelity"]) - float(sqr_best["noisy_relaxed_avg_state_fidelity"]),
            "snap_closed_delta": float(snap_big["closed_process_fidelity"]) - float(snap_best["closed_process_fidelity"]),
            "snap_noisy_delta": float(snap_big["noisy_avg_state_fidelity"]) - float(snap_best["noisy_avg_state_fidelity"]),
        },
    )

    sqr_identity = _zero_angle_sqr_basis_preservation(sqr_best)
    snap_identity = _zero_phase_snap_basis_preservation(snap_best)
    _record(
        checks,
        "basis_preservation_limits",
        sqr_identity > 0.999 and snap_identity > 0.995,
        {
            "sqr_zero_angle_basis_fidelity": sqr_identity,
            "snap_zero_phase_basis_fidelity": snap_identity,
        },
    )

    noisy_shorter_ok = float(sqr_best["duration_us"]) <= max(
        float(entry["duration_us"]) for entry in results["sqr"]["best_closed_by_family"].values()
    )
    _record(
        checks,
        "noise_prefers_shorter_sqr_point",
        noisy_shorter_ok,
        {
            "best_noisy_duration_us": sqr_best["duration_us"],
            "max_best_closed_duration_us": max(float(entry["duration_us"]) for entry in results["sqr"]["best_closed_by_family"].values()),
        },
    )

    passed = all(item["passed"] for item in checks)
    summary = {
        "all_passed": passed,
        "checks": checks,
        "sqr_best_noisy": sqr_best,
        "snap_best_noisy": snap_best,
    }
    VALIDATION_PATH.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    print(f"Saved {VALIDATION_PATH}")
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
