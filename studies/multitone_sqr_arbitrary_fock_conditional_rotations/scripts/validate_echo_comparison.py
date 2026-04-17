"""Validation checks for the explicit echoed-SQR comparison extension."""

from __future__ import annotations

from typing import Any

import numpy as np

from cqed_sim.calibration.conditioned_multitone import build_conditioned_multitone_waveform
from cqed_sim.pulses.envelopes import MultitoneTone

from common import (
    ARTIFACTS_DIR,
    DATA_DIR,
    active_subspace_metrics,
    build_frame,
    build_model,
    compile_pulse_sequence,
    load_json,
    logical_levels,
    make_gaussian_qubit_rotation_pulse,
    make_run_config,
    restricted_operator_from_full,
    save_json,
    shift_pulse,
    simulate_full_operator_on_logical_inputs,
)

from run_echo_comparison import ECHO_CASE_DIR, PI_PULSE_DURATION_S, PI_PULSE_PHASE_RAD, PI_PULSE_SIGMA_FRACTION


RESULTS_PATH = DATA_DIR / "echo_comparison_results.json"
SUMMARY_PATH = DATA_DIR / "echo_comparison_summary.json"
VALIDATION_PATH = DATA_DIR / "echo_comparison_validation.json"
BASELINE_VALIDATION_PATH = DATA_DIR / "validation_summary.json"


def _restore_complex_array(value: Any) -> np.ndarray:
    if isinstance(value, dict) and {"real", "imag", "shape"}.issubset(value):
        real = np.asarray(value["real"], dtype=float)
        imag = np.asarray(value["imag"], dtype=float)
        shape = tuple(int(item) for item in value["shape"])
        return (real + 1.0j * imag).reshape(shape)
    return np.asarray(value, dtype=np.complex128)


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


def _record(checks: list[dict[str, Any]], name: str, passed: bool, details: dict[str, Any]) -> None:
    checks.append({"name": name, "passed": bool(passed), "details": details})
    print(("PASS" if passed else "FAIL") + f": {name}")
    for key, value in details.items():
        print(f"  {key}: {value}")


def _load_best_echo_artifact() -> dict[str, Any]:
    summary = load_json(SUMMARY_PATH)
    best_echo = summary["best_echo"]
    return load_json(ECHO_CASE_DIR / f"{best_echo['case_id']}_{best_echo['sequence_family']}.json")


def _replay_echo_artifact(
    artifact: dict[str, Any],
    *,
    dt_s: float | None = None,
    extra_cavity_levels: int = 0,
) -> dict[str, float]:
    request = dict(artifact["case_request"])
    n_active = int(request["n_active"])
    model = build_model(
        include_chi_prime=bool(request["include_chi_prime"]),
        n_active=n_active + int(extra_cavity_levels),
    )
    frame = build_frame(model)
    half_duration_s = float(artifact["half_duration_s"])
    run_config = make_run_config(
        model,
        n_active=n_active,
        duration_s=half_duration_s,
        dt_s=float(dt_s) if dt_s is not None else 4.0e-9,
    )
    tones = _tone_specs_from_rows(artifact["half_tone_specs"])
    waveform = build_conditioned_multitone_waveform(tones, run_config, label="validation_half")
    half_pulse = waveform.pulse
    x_first = make_gaussian_qubit_rotation_pulse(
        model,
        frame,
        theta=np.pi,
        phase=PI_PULSE_PHASE_RAD,
        duration_s=PI_PULSE_DURATION_S,
        channel=str(half_pulse.channel),
        manifold_level=0,
        sigma_fraction=PI_PULSE_SIGMA_FRACTION,
        t0=half_duration_s,
        label="x_pi_1",
    )
    x_second = make_gaussian_qubit_rotation_pulse(
        model,
        frame,
        theta=np.pi,
        phase=PI_PULSE_PHASE_RAD,
        duration_s=PI_PULSE_DURATION_S,
        channel=str(half_pulse.channel),
        manifold_level=0,
        sigma_fraction=PI_PULSE_SIGMA_FRACTION,
        t0=2.0 * half_duration_s + PI_PULSE_DURATION_S,
        label="x_pi_2",
    )
    total_duration_s = float(artifact["total_gate_duration_s"])
    pulses = [
        shift_pulse(half_pulse, t0=0.0, label="half_1"),
        x_first,
        shift_pulse(half_pulse, t0=half_duration_s + PI_PULSE_DURATION_S, label="half_2"),
        x_second,
    ]
    compiled = compile_pulse_sequence(pulses, dt_s=float(run_config.dt_s), total_duration_s=total_duration_s)
    full_operator = simulate_full_operator_on_logical_inputs(
        model,
        compiled,
        frame=frame,
        drive_ops=waveform.drive_ops,
        levels=logical_levels(n_active),
    )
    restricted = restricted_operator_from_full(full_operator, model, logical_levels(n_active))
    target = _restore_complex_array(artifact["target_operator"])
    return active_subspace_metrics(target, restricted)


def main() -> None:
    checks: list[dict[str, Any]] = []

    required_files = [
        RESULTS_PATH,
        SUMMARY_PATH,
        DATA_DIR / "echo_comparison_summary.md",
        DATA_DIR / "validation_summary.json",
        ARTIFACTS_DIR / "echo_comparison" / "highlights" / "best_single_pulse.json",
        ARTIFACTS_DIR / "echo_comparison" / "highlights" / "best_echoed_sqr.json",
        ARTIFACTS_DIR / "echo_comparison" / "cases" / "chi_only_na3_chiT5p0_familyD_seed317160_echoed_fixed_total.json",
        ARTIFACTS_DIR / "echo_comparison" / "waveforms" / "chi_only_na3_chiT5p0_familyD_seed317160_echoed_fixed_total.npz",
        ARTIFACTS_DIR.parent / "figures" / "echo_branch_metric_means.png",
        ARTIFACTS_DIR.parent / "figures" / "echo_duration_scan.png",
        ARTIFACTS_DIR.parent / "figures" / "echo_delta_tradeoff.png",
        ARTIFACTS_DIR.parent / "figures" / "echo_best_waveforms.png",
        ARTIFACTS_DIR.parent / "figures" / "echo_block_error_breakdown.png",
    ]
    artifact_ok = all(path.exists() and path.stat().st_size > 0 for path in required_files)
    _record(
        checks,
        "artifacts_present",
        artifact_ok,
        {path.name: bool(path.exists() and path.stat().st_size > 0) for path in required_files},
    )

    baseline_validation = load_json(BASELINE_VALIDATION_PATH)
    _record(
        checks,
        "baseline_single_pulse_validation_inherited",
        bool(baseline_validation.get("all_passed", False)),
        {"baseline_all_passed": bool(baseline_validation.get("all_passed", False))},
    )

    best_echo_artifact = _load_best_echo_artifact()
    sequence_kinds = [step["kind"] for step in best_echo_artifact["sequence_spec"]]
    _record(
        checks,
        "explicit_echo_sequence_order",
        sequence_kinds == ["half_sqr", "x_pi", "half_sqr", "x_pi"],
        {"sequence_kinds": sequence_kinds},
    )

    summary = load_json(SUMMARY_PATH)
    best_echo = dict(summary["best_echo"])
    best_single = dict(summary["best_single"])
    _record(
        checks,
        "previous_study_was_single_pulse_only",
        best_single["sequence_family"] == "single_pulse" and best_echo["sequence_family"].startswith("echoed_"),
        {
            "best_single_sequence_family": best_single["sequence_family"],
            "best_echo_sequence_family": best_echo["sequence_family"],
        },
    )

    replay_dt = _replay_echo_artifact(best_echo_artifact, dt_s=2.0e-9)
    replay_trunc = _replay_echo_artifact(best_echo_artifact, extra_cavity_levels=2)
    base_fidelity = float(best_echo["average_gate_fidelity"])
    dt_delta = float(replay_dt["average_gate_fidelity"] - base_fidelity)
    trunc_delta = float(replay_trunc["average_gate_fidelity"] - base_fidelity)
    _record(
        checks,
        "best_echo_time_step_convergence",
        abs(dt_delta) < 2.0e-2,
        {"baseline": base_fidelity, "replayed": replay_dt["average_gate_fidelity"], "delta": dt_delta},
    )
    _record(
        checks,
        "best_echo_truncation_convergence",
        abs(trunc_delta) < 2.0e-2,
        {"baseline": base_fidelity, "replayed": replay_trunc["average_gate_fidelity"], "delta": trunc_delta},
    )

    _record(
        checks,
        "literature_comparison",
        True,
        {"status": "not_applicable", "reason": "The echoed extension is an original optimization-and-analysis study, not a reproduction benchmark."},
    )

    payload = {
        "all_passed": all(item["passed"] for item in checks),
        "checks": checks,
        "best_echo_case_id": best_echo["case_id"],
        "best_echo_sequence_family": best_echo["sequence_family"],
    }
    save_json(VALIDATION_PATH, payload)
    print(f"Saved {VALIDATION_PATH}")
    if not payload["all_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()