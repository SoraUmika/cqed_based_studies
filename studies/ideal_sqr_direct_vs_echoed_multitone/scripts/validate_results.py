"""Validation checks for the ideal-SQR direct-vs-echoed follow-up study."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from cqed_sim.calibration.conditioned_multitone import ConditionedMultitoneCorrections, compile_conditioned_multitone_waveform

from common import (
    DATA_DIR,
    PI_PULSE_DURATION_S,
    PI_PULSE_SIGMA_FRACTION,
    average_gate_fidelity,
    build_frame,
    build_model,
    build_multitone_waveform_from_corrections,
    build_target_operator,
    corrections_to_dict,
    duration_from_chi_t,
    load_json,
    logical_levels,
    make_gaussian_qubit_rotation_pulse,
    make_run_config,
    restricted_operator_from_full,
    save_json,
    shift_pulse,
    simulate_full_operator_on_logical_inputs,
    target_spec,
)
from run_study import CONSTRUCTIONS, CaseRequest, build_case_context, evaluate_direct_case, evaluate_symmetric_echo_case


RESULTS_PATH = DATA_DIR / "study_results.json"
SUMMARY_PATH = DATA_DIR / "study_summary.json"
VALIDATION_PATH = DATA_DIR / "validation_summary.json"


def _corrections_from_dict(payload: dict[str, Any]) -> ConditionedMultitoneCorrections:
    return ConditionedMultitoneCorrections(
        d_lambda=tuple(float(x) for x in payload["d_lambda"]),
        d_alpha=tuple(float(x) for x in payload["d_alpha"]),
        d_omega_rad_s=tuple(float(x) for x in payload["d_omega_rad_s"]),
    )


def _rebuild_direct(artifact: dict[str, Any], *, dt_s: float, n_cav_padding: int) -> float:
    request = artifact["case_request"]
    spec = target_spec(str(request["target_family"]), int(request["n_active"]))
    model = build_model(
        include_chi_prime=bool(request["include_chi_prime"]),
        n_active=int(request["n_active"]),
        n_cav_padding=int(n_cav_padding),
    )
    run_config = make_run_config(
        model,
        n_active=int(request["n_active"]),
        duration_s=duration_from_chi_t(float(request["chi_t_over_2pi"])),
        dt_s=float(dt_s),
    )
    levels = logical_levels(int(request["n_active"]))
    target_operator = build_target_operator(spec, levels)
    corrections = _corrections_from_dict(artifact["metadata"]["corrections_segment_1"])
    waveform, _ = build_multitone_waveform_from_corrections(model, spec, run_config, corrections=corrections, label="rebuild_direct")
    compiled = compile_conditioned_multitone_waveform(waveform, run_config)
    full_operator = simulate_full_operator_on_logical_inputs(model, compiled, frame=build_frame(model), drive_ops=waveform.drive_ops, levels=levels)
    restricted = restricted_operator_from_full(full_operator, model, levels)
    return float(average_gate_fidelity(target_operator, restricted))


def _rebuild_echo(artifact: dict[str, Any], *, dt_s: float, n_cav_padding: int) -> float:
    request = artifact["case_request"]
    spec = target_spec(str(request["target_family"]), int(request["n_active"]))
    half_spec = type(spec)(
        family=spec.family,
        theta_values=tuple(float(value / 2.0) for value in spec.theta_values),
        phi_values=spec.phi_values,
        metadata=dict(spec.metadata),
    )
    model = build_model(
        include_chi_prime=bool(request["include_chi_prime"]),
        n_active=int(request["n_active"]),
        n_cav_padding=int(n_cav_padding),
    )
    frame = build_frame(model)
    levels = logical_levels(int(request["n_active"]))
    target_operator = build_target_operator(spec, levels)
    meta = artifact["metadata"]
    corr_1 = _corrections_from_dict(meta["corrections_segment_1"])
    corr_2 = _corrections_from_dict(meta["corrections_segment_2"])

    run_1 = make_run_config(model, n_active=int(request["n_active"]), duration_s=float(meta["segment_1_duration_s"]), dt_s=float(dt_s))
    run_2 = make_run_config(model, n_active=int(request["n_active"]), duration_s=float(meta["segment_2_duration_s"]), dt_s=float(dt_s))
    waveform_1, _ = build_multitone_waveform_from_corrections(model, half_spec, run_1, corrections=corr_1, label="rebuild_seg1")
    waveform_2, _ = build_multitone_waveform_from_corrections(model, half_spec, run_2, corrections=corr_2, label="rebuild_seg2")

    x_first = make_gaussian_qubit_rotation_pulse(
        model,
        frame,
        theta=np.pi,
        phase=0.0,
        duration_s=float(meta["pi_pulse_duration_s"]),
        channel=str(waveform_1.pulse.channel),
        manifold_level=0,
        sigma_fraction=PI_PULSE_SIGMA_FRACTION,
        t0=float(meta["segment_1_duration_s"]),
        label="rebuild_x1",
    )
    x_second = make_gaussian_qubit_rotation_pulse(
        model,
        frame,
        theta=np.pi,
        phase=0.0,
        duration_s=float(meta["pi_pulse_duration_s"]),
        channel=str(waveform_1.pulse.channel),
        manifold_level=0,
        sigma_fraction=PI_PULSE_SIGMA_FRACTION,
        t0=float(meta["segment_1_duration_s"] + meta["pi_pulse_duration_s"] + meta["segment_2_duration_s"]),
        label="rebuild_x2",
    )
    from common import compile_pulse_sequence

    pulses = [
        shift_pulse(waveform_1.pulse, t0=0.0, label="rebuild_seg1"),
        x_first,
        shift_pulse(waveform_2.pulse, t0=float(meta["segment_1_duration_s"] + meta["pi_pulse_duration_s"]), label="rebuild_seg2"),
        x_second,
    ]
    compiled = compile_pulse_sequence(pulses, dt_s=float(dt_s), total_duration_s=float(meta["total_gate_duration_s"]))
    full_operator = simulate_full_operator_on_logical_inputs(model, compiled, frame=frame, drive_ops=waveform_1.drive_ops, levels=levels)
    restricted = restricted_operator_from_full(full_operator, model, levels)
    return float(average_gate_fidelity(target_operator, restricted))


def _artifact_for(case_id: str, construction: str) -> dict[str, Any]:
    path = Path(__file__).resolve().parent.parent / "artifacts" / "cases" / f"{case_id}_{construction}.json"
    return load_json(path)


def main() -> None:
    results = load_json(RESULTS_PATH)
    summary = load_json(SUMMARY_PATH)
    rows = list(results["case_rows"])
    n_cases = len({str(row["case_id"]) for row in rows})
    expected_rows = n_cases * len(CONSTRUCTIONS)

    sanity_request = CaseRequest(
        model_variant="chi_only",
        include_chi_prime=False,
        target_family="smooth_x",
        n_active=1,
        chi_t_over_2pi=5.0,
    )
    sanity_context = build_case_context(sanity_request)
    sanity_direct = evaluate_direct_case(sanity_context)[0]
    sanity_symmetric = evaluate_symmetric_echo_case(sanity_context)[0]

    best_direct = next(item for item in summary["best_by_construction"] if item["construction"] == "direct_multitone")
    best_echo = max(
        (item for item in summary["best_by_construction"] if item["construction"] != "direct_multitone"),
        key=lambda item: float(item["average_gate_fidelity"]),
    )

    best_direct_artifact = _artifact_for(str(best_direct["case_id"]), "direct_multitone")
    best_echo_artifact = _artifact_for(str(best_echo["case_id"]), str(best_echo["construction"]))

    convergence_rows = []
    for label, baseline, artifact, construction in (
        ("best_direct", best_direct, best_direct_artifact, "direct_multitone"),
        ("best_echo", best_echo, best_echo_artifact, str(best_echo["construction"])),
    ):
        base_fid = float(baseline["average_gate_fidelity"])
        if construction == "direct_multitone":
            dt2_fid = _rebuild_direct(artifact, dt_s=2.0e-9, n_cav_padding=2)
            pad3_fid = _rebuild_direct(artifact, dt_s=4.0e-9, n_cav_padding=3)
        else:
            dt2_fid = _rebuild_echo(artifact, dt_s=2.0e-9, n_cav_padding=2)
            pad3_fid = _rebuild_echo(artifact, dt_s=4.0e-9, n_cav_padding=3)
        convergence_rows.append(
            {
                "label": label,
                "case_id": str(baseline["case_id"]),
                "construction": construction,
                "baseline_fidelity": base_fid,
                "dt2_fidelity": float(dt2_fid),
                "padding3_fidelity": float(pad3_fid),
                "delta_dt2": float(dt2_fid - base_fid),
                "delta_padding3": float(pad3_fid - base_fid),
            }
        )

    payload = {
        "n_rows": int(len(rows)),
        "expected_rows": int(expected_rows),
        "row_count_ok": bool(len(rows) == expected_rows),
        "construction_set_ok": sorted({str(row["construction"]) for row in rows}) == sorted(CONSTRUCTIONS),
        "sanity_case": {
            "request": sanity_request.__dict__,
            "direct_fidelity": float(sanity_direct["average_gate_fidelity"]),
            "symmetric_echo_fidelity": float(sanity_symmetric["average_gate_fidelity"]),
            "delta": float(sanity_symmetric["average_gate_fidelity"] - sanity_direct["average_gate_fidelity"]),
            "pass": bool(
                sanity_direct["average_gate_fidelity"] > 0.80
                and abs(sanity_symmetric["average_gate_fidelity"] - sanity_direct["average_gate_fidelity"]) < 0.02
            ),
        },
        "convergence_rows": convergence_rows,
        "convergence_pass": bool(
            max(abs(item["delta_dt2"]) for item in convergence_rows) < 0.03
            and max(abs(item["delta_padding3"]) for item in convergence_rows) < 0.03
        ),
        "literature_comparison": "not_applicable_for_this_followup",
    }
    save_json(VALIDATION_PATH, payload)


if __name__ == "__main__":
    main()
