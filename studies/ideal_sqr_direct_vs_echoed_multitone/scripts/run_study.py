"""Run the ideal-SQR direct-vs-echoed multitone follow-up study."""

from __future__ import annotations

import argparse
import time
from dataclasses import dataclass
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import Bounds, minimize

from cqed_sim.calibration.conditioned_multitone import (
    ConditionedMultitoneCorrections,
    ConditionedOptimizationConfig,
    compile_conditioned_multitone_waveform,
)
from cqed_sim.calibration.targeted_subspace_multitone import (
    TargetedSubspaceObjectiveWeights,
    TargetedSubspaceOptimizationConfig,
    analyze_targeted_subspace_operator,
    build_spanning_state_transfer_set,
    optimize_targeted_subspace_multitone,
)

from common import (
    ARTIFACTS_DIR,
    DATA_DIR,
    FIGURES_DIR,
    IDEAL_X_PI,
    PI_PULSE_DURATION_S,
    PI_PULSE_SIGMA_FRACTION,
    STUDY_DIR,
    TWO_PI,
    apply_plot_style,
    average_gate_fidelity,
    block_rotation_metrics,
    build_frame,
    build_model,
    build_multitone_waveform_from_corrections,
    build_target_operator,
    channel_waveform_samples,
    compile_pulse_sequence,
    conditioned_targets_from_target_spec,
    corrections_from_vector,
    corrections_to_dict,
    corrections_to_vector,
    duration_from_chi_t,
    frobenius_error,
    json_ready,
    load_json,
    logical_levels,
    make_gaussian_qubit_rotation_pulse,
    make_run_config,
    restricted_blocks,
    restricted_operator_from_full,
    save_json,
    save_waveform_npz,
    shift_pulse,
    simulate_full_operator_on_logical_inputs,
    state_validation_summary_for_compiled,
    target_spec,
)


RESULTS_PATH = DATA_DIR / "study_results.json"
SUMMARY_PATH = DATA_DIR / "study_summary.json"
CSV_PATH = DATA_DIR / "study_results.csv"
MARKDOWN_PATH = DATA_DIR / "study_summary.md"
ANALYTIC_PATH = DATA_DIR / "analytic_summary.json"

CASE_DIR = ARTIFACTS_DIR / "cases"
WAVEFORM_DIR = ARTIFACTS_DIR / "waveforms"

MODEL_VARIANTS = (
    ("chi_only", False),
    ("chi_plus_chiprime", True),
)
TARGET_FAMILIES = ("smooth_x", "staggered_x")
ACTIVE_GRID = (2, 3)
DURATION_GRID = (3.0, 5.0)
CONSTRUCTIONS = (
    "direct_multitone",
    "echoed_symmetric",
    "echoed_independent",
    "echoed_asymmetric",
)

OBJECTIVE_WEIGHTS = TargetedSubspaceObjectiveWeights(
    qubit_weight=0.15,
    subspace_weight=1.0,
    preservation_weight=0.35,
    leakage_weight=0.35,
)

CONSTRUCTION_COLORS = {
    "direct_multitone": "#4477AA",
    "echoed_symmetric": "#EE6677",
    "echoed_independent": "#228833",
    "echoed_asymmetric": "#CCBB44",
}


@dataclass(frozen=True)
class CaseRequest:
    model_variant: str
    include_chi_prime: bool
    target_family: str
    n_active: int
    chi_t_over_2pi: float

    @property
    def case_id(self) -> str:
        duration_label = str(self.chi_t_over_2pi).replace(".", "p")
        return f"{self.model_variant}_{self.target_family}_na{int(self.n_active)}_chiT{duration_label}"


@dataclass
class CaseContext:
    request: CaseRequest
    model: Any
    frame: Any
    levels: tuple[int, ...]
    duration_s: float
    run_config: Any
    spec: Any
    targets: Any
    target_operator: np.ndarray
    transfer_set: Any
    direct_optimization: Any
    half_spec: Any
    half_targets: Any
    half_target_operator: np.ndarray
    half_transfer_set: Any
    half_optimization: Any


def tone_rows(tone_specs: list[Any] | tuple[Any, ...]) -> list[dict[str, float]]:
    rows = []
    for tone in tone_specs:
        rows.append(
            {
                "manifold": int(tone.manifold),
                "omega_rad_s": float(tone.omega_rad_s),
                "omega_hz": float(tone.omega_rad_s / TWO_PI),
                "amp_rad_s": float(tone.amp_rad_s),
                "phase_rad": float(tone.phase_rad),
            }
        )
    return rows


def optimization_config(n_active: int, *, stage1: int = 10, stage2: int = 16) -> TargetedSubspaceOptimizationConfig:
    conditioned = ConditionedOptimizationConfig(
        active_levels=tuple(range(int(n_active))),
        parameters=("d_lambda", "d_alpha", "d_omega"),
        maxiter_stage1=int(stage1 + 2 * n_active),
        maxiter_stage2=int(stage2 + 2 * n_active),
        d_lambda_bounds=(-0.75, 0.75),
        d_alpha_bounds=(-np.pi, np.pi),
        d_omega_hz_bounds=(-3.0e6, 3.0e6),
        regularization_lambda=5.0e-4,
        regularization_alpha=5.0e-4,
        regularization_omega=5.0e-4,
    )
    return TargetedSubspaceOptimizationConfig(conditioned=conditioned, include_block_phase=False)


def case_requests() -> list[CaseRequest]:
    rows: list[CaseRequest] = []
    for model_variant, include_chi_prime in MODEL_VARIANTS:
        for target_family in TARGET_FAMILIES:
            for n_active in ACTIVE_GRID:
                for chi_t in DURATION_GRID:
                    rows.append(
                        CaseRequest(
                            model_variant=model_variant,
                            include_chi_prime=include_chi_prime,
                            target_family=target_family,
                            n_active=n_active,
                            chi_t_over_2pi=chi_t,
                        )
                    )
    return rows


def build_case_context(request: CaseRequest) -> CaseContext:
    spec = target_spec(request.target_family, request.n_active)
    targets = conditioned_targets_from_target_spec(spec)
    model = build_model(include_chi_prime=request.include_chi_prime, n_active=request.n_active)
    frame = build_frame(model)
    duration_s = duration_from_chi_t(request.chi_t_over_2pi)
    run_config = make_run_config(model, n_active=request.n_active, duration_s=duration_s)
    levels = logical_levels(request.n_active)
    target_operator = build_target_operator(spec, levels)
    transfer_set = build_spanning_state_transfer_set(target_operator)

    direct = optimize_targeted_subspace_multitone(
        model,
        targets,
        run_config,
        logical_levels=levels,
        optimization_config=optimization_config(request.n_active),
        objective_weights=OBJECTIVE_WEIGHTS,
        target_operator=target_operator,
        transfer_set=transfer_set,
        label=f"direct_{request.case_id}",
    )

    half_spec = type(spec)(
        family=spec.family,
        theta_values=tuple(float(value / 2.0) for value in spec.theta_values),
        phi_values=spec.phi_values,
        metadata={
            **dict(spec.metadata),
            "description": f"Half-angle target derived from {spec.family}",
            "theta_values_rad": [float(value / 2.0) for value in spec.theta_values],
            "parent_family": spec.family,
        },
    )
    half_targets = conditioned_targets_from_target_spec(half_spec)
    half_duration_s = 0.5 * duration_s
    half_run_config = make_run_config(model, n_active=request.n_active, duration_s=half_duration_s)
    half_target_operator = build_target_operator(half_spec, levels)
    half_transfer_set = build_spanning_state_transfer_set(half_target_operator)

    half = optimize_targeted_subspace_multitone(
        model,
        half_targets,
        half_run_config,
        logical_levels=levels,
        optimization_config=optimization_config(request.n_active),
        objective_weights=OBJECTIVE_WEIGHTS,
        target_operator=half_target_operator,
        transfer_set=half_transfer_set,
        label=f"half_{request.case_id}",
    )

    return CaseContext(
        request=request,
        model=model,
        frame=frame,
        levels=levels,
        duration_s=duration_s,
        run_config=run_config,
        spec=spec,
        targets=targets,
        target_operator=target_operator,
        transfer_set=transfer_set,
        direct_optimization=direct,
        half_spec=half_spec,
        half_targets=half_targets,
        half_target_operator=half_target_operator,
        half_transfer_set=half_transfer_set,
        half_optimization=half,
    )


def half_run_config_for_duration(context: CaseContext, duration_s: float):
    return make_run_config(context.model, n_active=context.request.n_active, duration_s=float(duration_s), dt_s=float(context.run_config.dt_s))


def build_echo_compiled(
    context: CaseContext,
    *,
    corrections_1: ConditionedMultitoneCorrections,
    corrections_2: ConditionedMultitoneCorrections,
    segment_1_duration_s: float,
    segment_2_duration_s: float,
    construction: str,
):
    run_1 = half_run_config_for_duration(context, segment_1_duration_s)
    run_2 = half_run_config_for_duration(context, segment_2_duration_s)
    waveform_1, tones_1 = build_multitone_waveform_from_corrections(
        context.model,
        context.half_spec,
        run_1,
        corrections=corrections_1,
        label=f"{construction}_seg1",
    )
    waveform_2, tones_2 = build_multitone_waveform_from_corrections(
        context.model,
        context.half_spec,
        run_2,
        corrections=corrections_2,
        label=f"{construction}_seg2",
    )

    channel = str(waveform_1.pulse.channel)
    x_first = make_gaussian_qubit_rotation_pulse(
        context.model,
        context.frame,
        theta=np.pi,
        phase=0.0,
        duration_s=PI_PULSE_DURATION_S,
        channel=channel,
        manifold_level=0,
        sigma_fraction=PI_PULSE_SIGMA_FRACTION,
        t0=segment_1_duration_s,
        label=f"{construction}_xpi_1",
    )
    x_second = make_gaussian_qubit_rotation_pulse(
        context.model,
        context.frame,
        theta=np.pi,
        phase=0.0,
        duration_s=PI_PULSE_DURATION_S,
        channel=channel,
        manifold_level=0,
        sigma_fraction=PI_PULSE_SIGMA_FRACTION,
        t0=segment_1_duration_s + PI_PULSE_DURATION_S + segment_2_duration_s,
        label=f"{construction}_xpi_2",
    )
    pulses = [
        shift_pulse(waveform_1.pulse, t0=0.0, label=f"{construction}_seg1"),
        x_first,
        shift_pulse(waveform_2.pulse, t0=segment_1_duration_s + PI_PULSE_DURATION_S, label=f"{construction}_seg2"),
        x_second,
    ]
    total_gate_duration_s = segment_1_duration_s + segment_2_duration_s + 2.0 * PI_PULSE_DURATION_S
    compiled = compile_pulse_sequence(pulses, dt_s=float(context.run_config.dt_s), total_duration_s=total_gate_duration_s)
    metadata = {
        "time_order": "half_sqr -> x_pi -> half_sqr -> x_pi",
        "construction": str(construction),
        "segment_1_duration_s": float(segment_1_duration_s),
        "segment_2_duration_s": float(segment_2_duration_s),
        "total_gate_duration_s": float(total_gate_duration_s),
        "active_sqr_duration_s": float(segment_1_duration_s + segment_2_duration_s),
        "duration_asymmetry_eta": float((segment_1_duration_s - segment_2_duration_s) / (segment_1_duration_s + segment_2_duration_s)),
        "pi_pulse_duration_s": float(PI_PULSE_DURATION_S),
        "pi_pulse_phase_rad": 0.0,
        "pi_pulse_axis": "x",
        "pi_pulse_type": "finite_gaussian",
        "half_waveforms_identical": bool(
            np.allclose(corrections_to_vector(corrections_1), corrections_to_vector(corrections_2))
            and np.isclose(segment_1_duration_s, segment_2_duration_s)
        ),
        "corrections_segment_1": corrections_to_dict(corrections_1),
        "corrections_segment_2": corrections_to_dict(corrections_2),
        "tone_rows_segment_1": tone_rows(tones_1),
        "tone_rows_segment_2": tone_rows(tones_2),
        "sequence_spec": [
            {"kind": "half_sqr", "label": "segment_1", "t0_s": 0.0, "duration_s": float(segment_1_duration_s)},
            {"kind": "x_pi", "label": "x_pi_1", "t0_s": float(segment_1_duration_s), "duration_s": float(PI_PULSE_DURATION_S)},
            {"kind": "half_sqr", "label": "segment_2", "t0_s": float(segment_1_duration_s + PI_PULSE_DURATION_S), "duration_s": float(segment_2_duration_s)},
            {"kind": "x_pi", "label": "x_pi_2", "t0_s": float(segment_1_duration_s + PI_PULSE_DURATION_S + segment_2_duration_s), "duration_s": float(PI_PULSE_DURATION_S)},
        ],
    }
    return compiled, waveform_1.drive_ops, metadata


def evaluate_compiled_candidate(
    context: CaseContext,
    *,
    construction: str,
    compiled: Any,
    drive_ops: dict[str, str],
    metadata: dict[str, Any],
    optimizer_payload: dict[str, Any] | None,
    target_rule: str,
) -> tuple[dict[str, Any], dict[str, Any]]:
    full_operator = simulate_full_operator_on_logical_inputs(
        context.model,
        compiled,
        frame=context.frame,
        drive_ops=drive_ops,
        levels=context.levels,
    )
    restricted_operator = restricted_operator_from_full(full_operator, context.model, context.levels)
    validation = analyze_targeted_subspace_operator(
        full_operator,
        context.model,
        context.targets,
        logical_levels=context.levels,
        objective_weights=OBJECTIVE_WEIGHTS,
        target_operator=context.target_operator,
        transfer_set=context.transfer_set,
        metadata={"construction": str(construction), "target_rule": str(target_rule), **dict(metadata)},
    )
    per_block = [
        {"level": int(level), **block_rotation_metrics(target_block, actual_block)}
        for level, target_block, actual_block in zip(
            context.levels,
            restricted_blocks(context.target_operator),
            restricted_blocks(restricted_operator),
            strict=True,
        )
    ]
    residual_z = np.asarray([float(item["residual_z_error_rad"]) for item in per_block], dtype=float)
    transverse = np.asarray([float(item["transverse_error_rad"]) for item in per_block], dtype=float)
    state_summary = state_validation_summary_for_compiled(
        context.model,
        compiled,
        frame=context.frame,
        drive_ops=drive_ops,
        levels=context.levels,
        target_operator=context.target_operator,
    )
    row = {
        "case_id": context.request.case_id,
        "construction": str(construction),
        "model_variant": context.request.model_variant,
        "include_chi_prime": bool(context.request.include_chi_prime),
        "target_family": context.request.target_family,
        "n_active": int(context.request.n_active),
        "chi_t_over_2pi": float(context.request.chi_t_over_2pi),
        "theta_values_rad": [float(x) for x in context.spec.theta_values],
        "phi_values_rad": [float(x) for x in context.spec.phi_values],
        "pulse_duration_ns": float(context.duration_s * 1.0e9),
        "active_sqr_duration_ns": float(metadata["active_sqr_duration_s"] * 1.0e9),
        "total_gate_duration_ns": float(metadata["total_gate_duration_s"] * 1.0e9),
        "segment_1_duration_ns": float(metadata["segment_1_duration_s"] * 1.0e9),
        "segment_2_duration_ns": float(metadata["segment_2_duration_s"] * 1.0e9),
        "duration_asymmetry_eta": float(metadata["duration_asymmetry_eta"]),
        "pi_pulse_duration_ns": float(metadata["pi_pulse_duration_s"] * 1.0e9),
        "half_waveforms_identical": bool(metadata["half_waveforms_identical"]),
        "average_gate_fidelity": float(average_gate_fidelity(context.target_operator, restricted_operator)),
        "restricted_process_fidelity": float(validation.restricted_process_fidelity),
        "restricted_fro_error": float(frobenius_error(context.target_operator, restricted_operator)),
        "restricted_unitarity_error": float(validation.restricted_unitarity_error),
        "state_transfer_fidelity_mean": float(validation.state_transfer_fidelity_mean),
        "state_transfer_fidelity_min": float(validation.state_transfer_fidelity_min),
        "same_block_population_mean": float(validation.same_block_population_mean),
        "same_block_population_min": float(validation.same_block_population_min),
        "other_target_population_mean": float(validation.other_target_population_mean),
        "other_target_population_max": float(validation.other_target_population_max),
        "leakage_outside_target_mean": float(validation.leakage_outside_target_mean),
        "leakage_outside_target_max": float(validation.leakage_outside_target_max),
        "weighted_loss": float(validation.weighted_loss),
        "mean_residual_z_error_rad": float(np.mean(residual_z)),
        "max_residual_z_error_rad": float(np.max(residual_z)),
        "mean_transverse_error_rad": float(np.mean(transverse)),
        "max_transverse_error_rad": float(np.max(transverse)),
        "per_block_average_gate_fidelities": [float(item["average_gate_fidelity"]) for item in per_block],
        "per_block_process_fidelities": [float(item["process_fidelity"]) for item in per_block],
        "per_block_rotation_angle_errors_rad": [float(item["rotation_angle_error_rad"]) for item in per_block],
        "per_block_rotation_axis_errors_rad": [float(item["rotation_axis_error_rad"]) for item in per_block],
        "per_block_residual_z_errors_rad": [float(item["residual_z_error_rad"]) for item in per_block],
        "per_block_transverse_errors_rad": [float(item["transverse_error_rad"]) for item in per_block],
        "state_validation_ground_fidelity": float(state_summary["states"][0]["state_fidelity"]),
        "state_validation_plus_fidelity": float(state_summary["states"][1]["state_fidelity"]),
        "optimizer_kind": None if optimizer_payload is None else optimizer_payload.get("kind"),
        "optimizer_success": None if optimizer_payload is None else optimizer_payload.get("success"),
        "optimizer_nfev": None if optimizer_payload is None else optimizer_payload.get("nfev"),
        "optimizer_runtime_s": None if optimizer_payload is None else optimizer_payload.get("runtime_s"),
        "target_operator_rule": str(target_rule),
    }
    waveform_samples = channel_waveform_samples(compiled)
    artifact = {
        "study_name": STUDY_DIR.name,
        "date_created": time.strftime("%Y-%m-%d"),
        "description": "Ideal-SQR feasibility study comparing direct and echoed multitone waveform parameterizations.",
        "case_request": json_ready(context.request.__dict__),
        "construction": str(construction),
        "target_spec": context.spec.metadata,
        "target_operator": context.target_operator,
        "restricted_operator": restricted_operator,
        "full_operator_columns_on_logical_inputs": full_operator,
        "ideal_x_pi": IDEAL_X_PI,
        "metadata": metadata,
        "per_block_metrics": per_block,
        "validation": validation.as_dict(),
        "state_validation": state_summary,
        "waveform_samples": waveform_samples,
        "optimizer": optimizer_payload,
        "summary_row": row,
        "load_instructions": "Load this JSON and inspect `metadata`, `optimizer`, `restricted_operator`, `per_block_metrics`, and `waveform_samples` to reproduce the saved direct or echoed case.",
    }
    return row, artifact


def save_case_artifact(case_id: str, construction: str, artifact: dict[str, Any]) -> None:
    stem = f"{case_id}_{construction}"
    npz_path = WAVEFORM_DIR / f"{stem}.npz"
    save_waveform_npz(npz_path, artifact["waveform_samples"])
    artifact["waveform_npz"] = str(npz_path.relative_to(STUDY_DIR))
    save_json(CASE_DIR / f"{stem}.json", artifact)


def evaluate_direct_case(context: CaseContext) -> tuple[dict[str, Any], dict[str, Any]]:
    compiled = compile_conditioned_multitone_waveform(context.direct_optimization.optimized_result.waveform, context.run_config)
    metadata = {
        "time_order": "single_direct_multitone",
        "construction": "direct_multitone",
        "segment_1_duration_s": float(context.duration_s),
        "segment_2_duration_s": 0.0,
        "total_gate_duration_s": float(context.duration_s),
        "active_sqr_duration_s": float(context.duration_s),
        "duration_asymmetry_eta": 0.0,
        "pi_pulse_duration_s": 0.0,
        "pi_pulse_phase_rad": None,
        "pi_pulse_axis": None,
        "pi_pulse_type": None,
        "half_waveforms_identical": True,
        "corrections_segment_1": corrections_to_dict(context.direct_optimization.optimized_corrections),
        "corrections_segment_2": None,
        "tone_rows_segment_1": tone_rows(context.direct_optimization.optimized_result.waveform.tone_specs),
        "tone_rows_segment_2": None,
        "sequence_spec": [{"kind": "direct_multitone", "label": "single_direct", "t0_s": 0.0, "duration_s": float(context.duration_s)}],
    }
    optimizer_payload = {
        "kind": "cqed_sim_targeted_subspace",
        "success": bool(context.direct_optimization.success_stage1 and context.direct_optimization.success_stage2),
        "nfev": int(len(context.direct_optimization.history)),
        "runtime_s": float(sum(item.get("runtime_s", 0.0) for item in context.direct_optimization.history if isinstance(item, dict))),
        "optimized_corrections": corrections_to_dict(context.direct_optimization.optimized_corrections),
        "history": context.direct_optimization.history,
        "improvement_summary": context.direct_optimization.improvement_summary(),
        "message_stage1": str(context.direct_optimization.message_stage1),
        "message_stage2": str(context.direct_optimization.message_stage2),
    }
    return evaluate_compiled_candidate(
        context,
        construction="direct_multitone",
        compiled=compiled,
        drive_ops=context.direct_optimization.optimized_result.waveform.drive_ops,
        metadata=metadata,
        optimizer_payload=optimizer_payload,
        target_rule="direct ideal-SQR target U = sum_n |n><n| tensor R_x(theta_n)",
    )


def evaluate_symmetric_echo_case(context: CaseContext) -> tuple[dict[str, Any], dict[str, Any]]:
    half_corrections = context.half_optimization.optimized_corrections
    compiled, drive_ops, metadata = build_echo_compiled(
        context,
        corrections_1=half_corrections,
        corrections_2=half_corrections,
        segment_1_duration_s=0.5 * context.duration_s,
        segment_2_duration_s=0.5 * context.duration_s,
        construction="echoed_symmetric",
    )
    optimizer_payload = {
        "kind": "cqed_sim_half_target_seed",
        "success": bool(context.half_optimization.success_stage1 and context.half_optimization.success_stage2),
        "nfev": int(len(context.half_optimization.history)),
        "runtime_s": float(sum(item.get("runtime_s", 0.0) for item in context.half_optimization.history if isinstance(item, dict))),
        "optimized_corrections": corrections_to_dict(half_corrections),
        "history": context.half_optimization.history,
        "improvement_summary": context.half_optimization.improvement_summary(),
        "message_stage1": str(context.half_optimization.message_stage1),
        "message_stage2": str(context.half_optimization.message_stage2),
    }
    return evaluate_compiled_candidate(
        context,
        construction="echoed_symmetric",
        compiled=compiled,
        drive_ops=drive_ops,
        metadata=metadata,
        optimizer_payload=optimizer_payload,
        target_rule="ideal echoed algebra with identical half targets R_x(theta_n/2) and two finite X_pi refocusing pulses",
    )


def echo_objective(
    context: CaseContext,
    vector: np.ndarray,
    *,
    asymmetric: bool,
    reference_vector: np.ndarray,
) -> float:
    n_active = int(context.request.n_active)
    arr = np.asarray(vector, dtype=float).reshape(-1)
    eta = 0.0 if not asymmetric else float(arr[-1])
    payload = arr if not asymmetric else arr[:-1]
    corr_1 = corrections_from_vector(payload[: 3 * n_active], n_active)
    corr_2 = corrections_from_vector(payload[3 * n_active :], n_active)
    if asymmetric:
        seg_1 = 0.5 * context.duration_s * (1.0 + eta)
        seg_2 = 0.5 * context.duration_s * (1.0 - eta)
    else:
        seg_1 = seg_2 = 0.5 * context.duration_s

    compiled, drive_ops, _ = build_echo_compiled(
        context,
        corrections_1=corr_1,
        corrections_2=corr_2,
        segment_1_duration_s=seg_1,
        segment_2_duration_s=seg_2,
        construction="echoed_asymmetric" if asymmetric else "echoed_independent",
    )
    full_operator = simulate_full_operator_on_logical_inputs(
        context.model,
        compiled,
        frame=context.frame,
        drive_ops=drive_ops,
        levels=context.levels,
    )
    restricted_operator = restricted_operator_from_full(full_operator, context.model, context.levels)
    validation = analyze_targeted_subspace_operator(
        full_operator,
        context.model,
        context.targets,
        logical_levels=context.levels,
        objective_weights=OBJECTIVE_WEIGHTS,
        target_operator=context.target_operator,
        transfer_set=context.transfer_set,
        metadata={"objective_mode": "composite_echo_local"},
    )
    per_block = [
        block_rotation_metrics(target_block, actual_block)
        for target_block, actual_block in zip(
            restricted_blocks(context.target_operator),
            restricted_blocks(restricted_operator),
            strict=True,
        )
    ]
    mean_residual_z = float(np.mean([item["residual_z_error_rad"] for item in per_block]))
    drift_penalty = 1.0e-5 * float(np.linalg.norm(arr[: reference_vector.size] - reference_vector) ** 2)
    asym_penalty = 0.01 * abs(eta)
    return float(validation.weighted_loss) + 0.05 * mean_residual_z + drift_penalty + asym_penalty


def optimize_echo_variant(context: CaseContext, *, asymmetric: bool) -> tuple[dict[str, Any], dict[str, Any]]:
    n_active = int(context.request.n_active)
    seed_corr = context.half_optimization.optimized_corrections
    reference_vector = np.concatenate([corrections_to_vector(seed_corr), corrections_to_vector(seed_corr)])
    x0 = reference_vector.copy()
    lower = np.concatenate(
        [
            np.full(n_active, -0.75),
            np.full(n_active, -np.pi),
            np.full(n_active, -3.0e6 * TWO_PI),
            np.full(n_active, -0.75),
            np.full(n_active, -np.pi),
            np.full(n_active, -3.0e6 * TWO_PI),
        ]
    )
    upper = np.concatenate(
        [
            np.full(n_active, 0.75),
            np.full(n_active, np.pi),
            np.full(n_active, 3.0e6 * TWO_PI),
            np.full(n_active, 0.75),
            np.full(n_active, np.pi),
            np.full(n_active, 3.0e6 * TWO_PI),
        ]
    )
    if asymmetric:
        x0 = np.concatenate([x0, np.asarray([0.0], dtype=float)])
        lower = np.concatenate([lower, np.asarray([-0.15], dtype=float)])
        upper = np.concatenate([upper, np.asarray([0.15], dtype=float)])

    start = time.perf_counter()
    result = minimize(
        lambda x: echo_objective(context, x, asymmetric=asymmetric, reference_vector=reference_vector),
        x0,
        method="Powell",
        bounds=Bounds(lower, upper),
        options={"maxiter": 5, "maxfev": 80 if asymmetric else 70, "xtol": 5.0e-3, "ftol": 5.0e-3},
    )
    runtime_s = time.perf_counter() - start

    eta = 0.0 if not asymmetric else float(result.x[-1])
    payload = np.asarray(result.x[:-1] if asymmetric else result.x, dtype=float)
    corr_1 = corrections_from_vector(payload[: 3 * n_active], n_active)
    corr_2 = corrections_from_vector(payload[3 * n_active :], n_active)
    seg_1 = 0.5 * context.duration_s * (1.0 + eta)
    seg_2 = 0.5 * context.duration_s * (1.0 - eta)
    construction = "echoed_asymmetric" if asymmetric else "echoed_independent"

    compiled, drive_ops, metadata = build_echo_compiled(
        context,
        corrections_1=corr_1,
        corrections_2=corr_2,
        segment_1_duration_s=seg_1,
        segment_2_duration_s=seg_2,
        construction=construction,
    )
    optimizer_payload = {
        "kind": "powell_composite_echo",
        "success": bool(result.success),
        "message": str(result.message),
        "nfev": int(getattr(result, "nfev", -1)),
        "nit": int(getattr(result, "nit", -1)) if getattr(result, "nit", None) is not None else -1,
        "fun": float(result.fun),
        "runtime_s": float(runtime_s),
        "x": np.asarray(result.x, dtype=float).tolist(),
        "optimized_corrections_segment_1": corrections_to_dict(corr_1),
        "optimized_corrections_segment_2": corrections_to_dict(corr_2),
    }
    target_rule = (
        "ideal echoed algebra with independent half corrections and weak duration asymmetry"
        if asymmetric
        else "ideal echoed algebra with independent half corrections and equal durations"
    )
    return evaluate_compiled_candidate(
        context,
        construction=construction,
        compiled=compiled,
        drive_ops=drive_ops,
        metadata=metadata,
        optimizer_payload=optimizer_payload,
        target_rule=target_rule,
    )


def run_case(context: CaseContext) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []

    direct_row, direct_artifact = evaluate_direct_case(context)
    save_case_artifact(context.request.case_id, "direct_multitone", direct_artifact)
    rows.append(direct_row)

    symmetric_row, symmetric_artifact = evaluate_symmetric_echo_case(context)
    save_case_artifact(context.request.case_id, "echoed_symmetric", symmetric_artifact)
    rows.append(symmetric_row)

    independent_row, independent_artifact = optimize_echo_variant(context, asymmetric=False)
    save_case_artifact(context.request.case_id, "echoed_independent", independent_artifact)
    rows.append(independent_row)

    asymmetric_row, asymmetric_artifact = optimize_echo_variant(context, asymmetric=True)
    save_case_artifact(context.request.case_id, "echoed_asymmetric", asymmetric_artifact)
    rows.append(asymmetric_row)
    return rows


def analytic_summary() -> dict[str, Any]:
    notes = [
        {
            "id": "audit_target_mismatch",
            "statement": "The earlier baseline arbitrary-rotation study optimized arbitrary SU(2) block targets rather than an ideal x-axis SQR gate.",
        },
        {
            "id": "audit_echo_mismatch",
            "statement": "The earlier residual-Z follow-up used only one mid-sequence X_pi pulse, so it did not test the requested half-SQR -> pi -> half-SQR -> pi construction.",
        },
        {
            "id": "direct_first_order",
            "equation": "U_n^dir approx exp[-i (Theta_n X + Phi_n Z)/2], with Theta_n = int Omega_n dt and Phi_n = int delta_n dt.",
            "conclusion": "A direct ideal SQR requires Theta_n = theta_n and Phi_n = 0 for every active manifold.",
        },
        {
            "id": "echo_first_order",
            "equation": "U_n^echo approx exp[-i ((Theta_{n,1}+Theta_{n,2}) X + (Phi_{n,1}-Phi_{n,2}) Z)/2] when X_pi is ideal and the target axis is X.",
            "conclusion": "A symmetric echoed sequence cancels first-order Z-type phase if the two halves accumulate the same Phi_n.",
        },
        {
            "id": "finite_pi_obstruction",
            "statement": "Once the inserted X_pi pulses become manifold dependent, the clean toggling-frame cancellation is spoiled and the echoed protocol can fail even when the ideal first-order algebra is favorable.",
        },
    ]
    save_json(ANALYTIC_PATH, {"notes": notes})
    return {"notes": notes}


def build_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    df = pd.DataFrame(rows)
    direct = df[df["construction"] == "direct_multitone"].copy()
    echoes = df[df["construction"] != "direct_multitone"].copy()
    deltas = echoes.merge(
        direct[[
            "case_id",
            "average_gate_fidelity",
            "mean_residual_z_error_rad",
            "mean_transverse_error_rad",
            "total_gate_duration_ns",
        ]],
        on="case_id",
        suffixes=("", "_direct"),
    )
    deltas["delta_fidelity"] = deltas["average_gate_fidelity"] - deltas["average_gate_fidelity_direct"]
    deltas["delta_residual_z"] = deltas["mean_residual_z_error_rad"] - deltas["mean_residual_z_error_rad_direct"]
    deltas["delta_transverse"] = deltas["mean_transverse_error_rad"] - deltas["mean_transverse_error_rad_direct"]
    deltas["delta_total_duration_ns"] = deltas["total_gate_duration_ns"] - deltas["total_gate_duration_ns_direct"]

    construction_summary = (
        df.groupby("construction", as_index=False)
        .agg(
            avg_fidelity_mean=("average_gate_fidelity", "mean"),
            avg_fidelity_best=("average_gate_fidelity", "max"),
            residual_z_mean=("mean_residual_z_error_rad", "mean"),
            transverse_mean=("mean_transverse_error_rad", "mean"),
            leakage_mean=("leakage_outside_target_mean", "mean"),
            count=("case_id", "count"),
        )
        .sort_values("avg_fidelity_mean", ascending=False)
    )
    target_construction = (
        df.groupby(["target_family", "construction"], as_index=False)
        .agg(
            avg_fidelity_mean=("average_gate_fidelity", "mean"),
            avg_fidelity_median=("average_gate_fidelity", "median"),
            residual_z_mean=("mean_residual_z_error_rad", "mean"),
            transverse_mean=("mean_transverse_error_rad", "mean"),
        )
        .sort_values(["target_family", "avg_fidelity_mean"], ascending=[True, False])
    )
    best_by_construction = (
        df.sort_values("average_gate_fidelity", ascending=False)
        .groupby("construction", as_index=False)
        .first()
        .to_dict(orient="records")
    )
    best_overall = df.sort_values("average_gate_fidelity", ascending=False).iloc[0].to_dict()
    largest_echo_gain = None if deltas.empty else deltas.sort_values("delta_fidelity", ascending=False).iloc[0].to_dict()
    largest_residual_reduction = None if deltas.empty else deltas.sort_values("delta_residual_z", ascending=True).iloc[0].to_dict()
    return {
        "study": STUDY_DIR.name,
        "n_rows": int(len(rows)),
        "n_cases": int(df["case_id"].nunique()),
        "constructions": list(CONSTRUCTIONS),
        "audit_findings": analytic_summary()["notes"],
        "best_overall": best_overall,
        "best_by_construction": best_by_construction,
        "construction_summary": construction_summary.to_dict(orient="records"),
        "target_construction_summary": target_construction.to_dict(orient="records"),
        "largest_echo_gain": largest_echo_gain,
        "largest_residual_reduction": largest_residual_reduction,
        "rows": rows,
        "deltas": deltas.to_dict(orient="records"),
    }


def save_figure(fig, stem: str) -> None:
    fig.savefig(FIGURES_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_figures(df: pd.DataFrame) -> None:
    apply_plot_style()
    construction_order = df.groupby("construction")["average_gate_fidelity"].mean().sort_values(ascending=False).index.tolist()

    metrics = (
        df.groupby("construction")
        .agg(
            avg_fidelity=("average_gate_fidelity", "mean"),
            residual_z=("mean_residual_z_error_rad", "mean"),
            transverse=("mean_transverse_error_rad", "mean"),
        )
        .loc[construction_order]
    )
    colors = [CONSTRUCTION_COLORS[item] for item in metrics.index]

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.2))
    axes[0].bar(metrics.index, metrics["avg_fidelity"], color=colors)
    axes[0].set_ylabel("Mean average gate fidelity")
    axes[1].bar(metrics.index, metrics["residual_z"], color=colors)
    axes[1].set_ylabel("Mean residual Z error (rad)")
    axes[2].bar(metrics.index, metrics["transverse"], color=colors)
    axes[2].set_ylabel("Mean transverse error (rad)")
    for axis in axes:
        axis.set_xlabel("Construction")
        axis.tick_params(axis="x", rotation=25)
    save_figure(fig, "construction_metric_means")

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.4), sharey=True)
    for idx, target_family in enumerate(TARGET_FAMILIES):
        subset = df[df["target_family"] == target_family]
        ax = axes[idx]
        for construction in construction_order:
            branch = subset[subset["construction"] == construction]
            grouped = branch.groupby("chi_t_over_2pi", as_index=False).median(numeric_only=True)
            ax.plot(
                grouped["chi_t_over_2pi"],
                grouped["average_gate_fidelity"],
                marker="o",
                label=construction,
                color=CONSTRUCTION_COLORS[construction],
            )
        ax.set_title(target_family.replace("_", " "))
        ax.set_xlabel(r"$|\chi|T/2\pi$")
        ax.set_ylabel("Median average fidelity")
    axes[0].legend(frameon=False)
    save_figure(fig, "duration_fidelity_tradeoff")

    pivot = df.pivot(index="case_id", columns="construction", values="average_gate_fidelity").reindex(columns=construction_order)
    case_labels = []
    for case_id in pivot.index:
        row = df[df["case_id"] == case_id].iloc[0]
        case_labels.append(
            f"{row['model_variant']} | {row['target_family']} | na={int(row['n_active'])} | chiT={row['chi_t_over_2pi']:g}"
        )
    fig_height = max(5.5, 0.38 * len(case_labels))
    fig, ax = plt.subplots(figsize=(9.4, fig_height))
    image = ax.imshow(pivot.to_numpy(dtype=float), aspect="auto", cmap="viridis", vmin=0.0, vmax=1.0)
    ax.set_xticks(np.arange(len(pivot.columns)), pivot.columns, rotation=25, ha="right")
    ax.set_yticks(np.arange(len(case_labels)), case_labels)
    ax.set_xlabel("Construction")
    ax.set_ylabel("Study case")
    fig.colorbar(image, ax=ax, label="Average gate fidelity")
    save_figure(fig, "case_construction_heatmap")

    direct = df[df["construction"] == "direct_multitone"].copy()
    echoes = df[df["construction"] != "direct_multitone"].copy()
    deltas = echoes.merge(
        direct[["case_id", "average_gate_fidelity", "mean_residual_z_error_rad"]],
        on="case_id",
        suffixes=("", "_direct"),
    )
    deltas["delta_fidelity"] = deltas["average_gate_fidelity"] - deltas["average_gate_fidelity_direct"]
    deltas["delta_residual_z"] = deltas["mean_residual_z_error_rad"] - deltas["mean_residual_z_error_rad_direct"]
    fig, ax = plt.subplots(figsize=(6.8, 5.0))
    for construction in ("echoed_symmetric", "echoed_independent", "echoed_asymmetric"):
        branch = deltas[deltas["construction"] == construction]
        ax.scatter(
            branch["delta_residual_z"],
            branch["delta_fidelity"],
            s=54,
            alpha=0.85,
            color=CONSTRUCTION_COLORS[construction],
            label=construction,
        )
    ax.axhline(0.0, color="black", linewidth=0.8, linestyle="--")
    ax.axvline(0.0, color="black", linewidth=0.8, linestyle="--")
    ax.set_xlabel(r"$\Delta$ mean residual Z error (rad) vs direct")
    ax.set_ylabel(r"$\Delta$ average fidelity vs direct")
    ax.legend(frameon=False)
    save_figure(fig, "echo_delta_tradeoff")

    focus = df[
        (df["model_variant"] == "chi_plus_chiprime")
        & (df["target_family"] == "staggered_x")
        & (df["n_active"] == max(ACTIVE_GRID))
        & (df["chi_t_over_2pi"] == max(DURATION_GRID))
    ]
    if focus.empty:
        focus = df[df["construction"] == "direct_multitone"].nlargest(1, "average_gate_fidelity")
    case_id = str(focus.iloc[0]["case_id"])
    fig, axes = plt.subplots(len(construction_order), 1, figsize=(8.6, 2.4 * len(construction_order)), sharex=True)
    if len(construction_order) == 1:
        axes = [axes]
    for axis, construction in zip(axes, construction_order, strict=True):
        artifact = load_json(CASE_DIR / f"{case_id}_{construction}.json")
        samples = artifact["waveform_samples"]
        time_ns = 1.0e9 * np.asarray(samples["time_s"], dtype=float)
        signal = np.asarray(samples["baseband_real"], dtype=float) + 1.0j * np.asarray(samples["baseband_imag"], dtype=float)
        axis.plot(time_ns, np.real(signal), label="I", color="#4477AA")
        axis.plot(time_ns, np.imag(signal), label="Q", color="#EE6677")
        axis.set_ylabel(construction)
    axes[0].legend(frameon=False, ncol=2)
    axes[-1].set_xlabel("Time (ns)")
    save_figure(fig, "representative_waveforms")


def write_markdown_summary(summary: dict[str, Any]) -> None:
    lines = [
        f"# Summary: {STUDY_DIR.name}",
        "",
        "## Audit Findings",
    ]
    for item in summary["audit_findings"]:
        text = item.get("statement") or item.get("conclusion") or item.get("equation")
        lines.append(f"- {text}")
    lines.extend(["", "## Best Overall"])
    lines.append(
        (
            "- {construction} on {case_id} with average gate fidelity {average_gate_fidelity:.6f}, "
            "mean residual-Z {mean_residual_z_error_rad:.6f} rad, mean transverse {mean_transverse_error_rad:.6f} rad."
        ).format(**summary["best_overall"])
    )
    lines.extend(["", "## Construction Means"])
    for row in summary["construction_summary"]:
        lines.append(
            "- {construction}: mean fidelity {avg_fidelity_mean:.6f}, best fidelity {avg_fidelity_best:.6f}, mean residual-Z {residual_z_mean:.6f} rad, mean transverse {transverse_mean:.6f} rad.".format(
                **row
            )
        )
    if summary["largest_echo_gain"] is not None:
        lines.extend(["", "## Largest Echo Gain"])
        lines.append(
            (
                "- {construction} on {case_id}: delta fidelity {delta_fidelity:+.6f}, "
                "delta residual-Z {delta_residual_z:+.6f} rad."
            ).format(**summary["largest_echo_gain"])
        )
    if summary["largest_residual_reduction"] is not None:
        heading = "## Largest Residual-Z Reduction" if float(summary["largest_residual_reduction"]["delta_residual_z"]) < 0.0 else "## Smallest Residual-Z Penalty"
        lines.extend(["", heading])
        lines.append(
            (
                "- {construction} on {case_id}: delta residual-Z {delta_residual_z:+.6f} rad, "
                "delta fidelity {delta_fidelity:+.6f}."
            ).format(**summary["largest_residual_reduction"])
        )
    MARKDOWN_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=str, default="", help="Optional case_id filter for a single study case.")
    args = parser.parse_args()

    rows: list[dict[str, Any]] = []
    requests = case_requests()
    if args.case:
        requests = [item for item in requests if item.case_id == args.case]
        if not requests:
            raise SystemExit(f"No case matched '{args.case}'.")

    for request in requests:
        print(f"[case] {request.case_id}")
        context = build_case_context(request)
        case_rows = run_case(context)
        rows.extend(case_rows)
        for row in case_rows:
            print(
                "  {construction}: fidelity={average_gate_fidelity:.6f} resZ={mean_residual_z_error_rad:.6f} "
                "trans={mean_transverse_error_rad:.6f}".format(**row)
            )

    payload = {"study": STUDY_DIR.name, "case_rows": rows}
    save_json(RESULTS_PATH, payload)
    pd.DataFrame(rows).to_csv(CSV_PATH, index=False)

    summary = build_summary(rows)
    save_json(SUMMARY_PATH, summary)
    write_markdown_summary(summary)
    plot_figures(pd.DataFrame(rows))


if __name__ == "__main__":
    main()
