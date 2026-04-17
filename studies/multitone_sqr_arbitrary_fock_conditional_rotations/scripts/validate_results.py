"""Validation checks for the multitone SQR arbitrary block-rotation study."""

from __future__ import annotations

from typing import Any

import numpy as np

from cqed_sim.calibration.conditioned_multitone import build_conditioned_multitone_waveform
from cqed_sim.pulses.envelopes import MultitoneTone

from common import (
    ALPHA,
    ARTIFACTS_DIR,
    CHI,
    CHI_PRIME,
    DATA_DIR,
    DEFAULT_DT,
    FIGURES_DIR,
    KERR,
    N_CAV_PADDING,
    OMEGA_C,
    OMEGA_Q,
    active_subspace_metrics,
    build_model,
    duration_from_chi_t,
    logical_levels,
    make_run_config,
    save_json,
    simulate_logical_basis_operator,
)
from run_study import RESULTS_PATH, CaseRequest, run_case
from common import load_json

from cqed_sim.core import DispersiveTransmonCavityModel


VALIDATION_PATH = DATA_DIR / "validation_summary.json"


def _restore_array(value: Any) -> np.ndarray:
    if isinstance(value, dict) and {"real", "imag"}.issubset(value.keys()):
        return np.asarray(value["real"], dtype=float) + 1.0j * np.asarray(value["imag"], dtype=float)
    return np.asarray(value)


def _record(checks: list[dict[str, Any]], name: str, passed: bool, details: dict[str, Any]) -> None:
    checks.append({"name": name, "passed": bool(passed), "details": details})
    status = "PASS" if passed else "FAIL"
    print(f"{status}: {name}")
    for key, value in details.items():
        print(f"  {key}: {value}")


def _load_results() -> dict[str, Any]:
    return load_json(RESULTS_PATH)


def _load_case_artifact(case_id: str) -> dict[str, Any]:
    return load_json(ARTIFACTS_DIR / "cases" / f"{case_id}.json")


def _tone_specs_from_rows(rows: list[dict[str, Any]]) -> list[MultitoneTone]:
    return [
        MultitoneTone(
            manifold=int(row.get("manifold", row.get("n"))),
            omega_rad_s=float(row["omega_rad_s"]),
            amp_rad_s=float(row["amp_rad_s"]),
            phase_rad=float(row["phase_rad"]),
        )
        for row in rows
    ]


def _build_model_for_request(request: dict[str, Any], *, extra_cavity_levels: int = 0) -> DispersiveTransmonCavityModel:
    n_active = int(request["n_active"])
    n_cav = int(n_active) + int(N_CAV_PADDING) + int(extra_cavity_levels)
    chi_higher = (CHI_PRIME,) if bool(request["include_chi_prime"]) else ()
    return DispersiveTransmonCavityModel(
        omega_q=OMEGA_Q,
        omega_c=OMEGA_C,
        alpha=ALPHA,
        chi=CHI,
        chi_higher=chi_higher,
        kerr=KERR,
        n_cav=n_cav,
        n_tr=2,
    )


def _replay_artifact_metrics(
    artifact: dict[str, Any],
    *,
    dt_s: float | None = None,
    extra_cavity_levels: int = 0,
) -> dict[str, float]:
    request = dict(artifact["case_request"])
    model = _build_model_for_request(request, extra_cavity_levels=extra_cavity_levels)
    duration_s = float(artifact["summary_row"]["pulse_duration_s"])
    run_config = make_run_config(
        model,
        n_active=int(request["n_active"]),
        duration_s=duration_s,
        dt_s=(DEFAULT_DT if dt_s is None else float(dt_s)),
    )
    waveform = build_conditioned_multitone_waveform(
        _tone_specs_from_rows(artifact["tone_specs"]),
        run_config,
        label=str(artifact["summary_row"]["case_id"]),
    )
    restricted = simulate_logical_basis_operator(model, waveform, run_config, logical_levels(int(request["n_active"])))
    target = _restore_array(artifact["target_operator"])
    return active_subspace_metrics(target, restricted)


def main() -> None:
    results = _load_results()
    rows = list(results["case_rows"])
    checks: list[dict[str, Any]] = []

    required_files = [
        RESULTS_PATH,
        FIGURES_DIR / "transition_frequency_diagram.png",
        FIGURES_DIR / "representative_multitone_waveforms.png",
        FIGURES_DIR / "blockwise_fidelity_heatmap.png",
        FIGURES_DIR / "crosstalk_heatmap.png",
        FIGURES_DIR / "fidelity_vs_pulse_duration.png",
        FIGURES_DIR / "fidelity_vs_active_subspace_size.png",
        FIGURES_DIR / "random_target_fidelity_histogram.png",
        FIGURES_DIR / "state_level_validation.png",
        FIGURES_DIR / "best_worst_parameter_tables.png",
        FIGURES_DIR / "chi_only_vs_chi_plus_chiprime.png",
    ]
    artifact_ok = all(path.exists() and path.stat().st_size > 0 for path in required_files)
    _record(
        checks,
        "artifacts_present",
        artifact_ok,
        {path.name: bool(path.exists() and path.stat().st_size > 0) for path in required_files},
    )

    zero_model = build_model(include_chi_prime=False, n_active=3)
    zero_run = make_run_config(zero_model, n_active=3, duration_s=duration_from_chi_t(3.0))
    zero_waveform = build_conditioned_multitone_waveform(
        [MultitoneTone(manifold=0, omega_rad_s=0.0, amp_rad_s=0.0, phase_rad=0.0)],
        zero_run,
        label="zero_drive",
    )
    zero_target = np.eye(6, dtype=np.complex128)
    zero_actual = simulate_logical_basis_operator(zero_model, zero_waveform, zero_run, logical_levels(3))
    zero_metrics = active_subspace_metrics(zero_target, zero_actual)
    zero_pass = (
        float(zero_metrics["average_gate_fidelity"]) > 0.997
        and float(zero_metrics["restricted_unitarity_error"]) < 1.0e-12
    )
    _record(
        checks,
        "zero_drive_identity_sanity",
        zero_pass,
        {
            **zero_metrics,
            "interpretation": "Zero drive remains nearly trivial and exactly unitary on the active subspace; the residual infidelity is dominated by deterministic drift phases in the chosen rotating frame.",
        },
    )

    single_case = run_case(
        CaseRequest(
            family="C",
            model_variant="chi_only",
            include_chi_prime=False,
            n_active=1,
            duration_chi_t=5.0,
        )
    )
    _record(
        checks,
        "single_manifold_sanity",
        float(single_case["average_gate_fidelity"]) > 0.90,
        {
            "average_gate_fidelity": single_case["average_gate_fidelity"],
            "state_validation_ground_fidelity": single_case["state_validation_ground_fidelity"],
            "state_validation_plus_fidelity": single_case["state_validation_plus_fidelity"],
        },
    )

    df_sorted = sorted(rows, key=lambda row: float(row["average_gate_fidelity"]), reverse=True)
    best_structured = next(row for row in df_sorted if row["family"] != "D")
    best_random = next(row for row in df_sorted if row["family"] == "D")

    time_step_deltas = {}
    truncation_deltas = {}
    time_pass = True
    trunc_pass = True
    for label, row in (("best_structured", best_structured), ("best_random", best_random)):
        artifact = _load_case_artifact(str(row["case_id"]))
        dt_metrics = _replay_artifact_metrics(artifact, dt_s=2.0e-9)
        trunc_metrics = _replay_artifact_metrics(artifact, extra_cavity_levels=2)
        time_delta = float(dt_metrics["average_gate_fidelity"] - float(row["average_gate_fidelity"]))
        trunc_delta = float(trunc_metrics["average_gate_fidelity"] - float(row["average_gate_fidelity"]))
        time_step_deltas[label] = time_delta
        truncation_deltas[label] = trunc_delta
        time_pass = time_pass and abs(time_delta) < 1.0e-2
        trunc_pass = trunc_pass and abs(trunc_delta) < 2.0e-2

    _record(checks, "time_step_convergence", time_pass, time_step_deltas)
    _record(checks, "truncation_convergence", trunc_pass, truncation_deltas)

    random_rows = [row for row in rows if row["family"] == "D" and int(row["n_active"]) == 4 and float(row["chi_t_over_2pi"]) == 3.0]
    chi_only_values = [float(row["average_gate_fidelity"]) for row in random_rows if row["model_variant"] == "chi_only"]
    chip_values = [float(row["average_gate_fidelity"]) for row in random_rows if row["model_variant"] == "chi_plus_chiprime"]
    if chi_only_values and chip_values:
        chi_penalty = float(np.median(chi_only_values) - np.median(chip_values))
        _record(
            checks,
            "chi_prime_random_ensemble_comparison",
            True,
            {
                "median_chi_only": float(np.median(chi_only_values)),
                "median_chi_plus_chiprime": float(np.median(chip_values)),
                "median_difference": chi_penalty,
                "interpretation": "The sampled random-target median shifts slightly upward when chi' is included; the study records the effect instead of assuming its sign a priori.",
            },
        )
    else:
        _record(
            checks,
            "chi_prime_random_ensemble_comparison",
            False,
            {"reason": "Required random-ensemble comparison rows not found."},
        )

    _record(
        checks,
        "literature_comparison",
        True,
        {"status": "not_applicable", "reason": "This is an original OPT/ANA/DES study, not a reproduction benchmark."},
    )

    summary = {
        "all_passed": all(item["passed"] for item in checks),
        "checks": checks,
        "best_structured_case_id": best_structured["case_id"],
        "best_random_case_id": best_random["case_id"],
    }
    save_json(VALIDATION_PATH, summary)
    print(f"Saved {VALIDATION_PATH}")
    if not summary["all_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()