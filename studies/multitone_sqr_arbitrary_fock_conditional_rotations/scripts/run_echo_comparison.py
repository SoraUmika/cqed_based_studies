"""Run the explicit echoed-SQR comparison study for arbitrary conditional rotations."""

from __future__ import annotations

import argparse
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cqed_sim.calibration.conditioned_multitone import ConditionedOptimizationConfig
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
    STUDY_DIR,
    active_subspace_metrics,
    apply_plot_style,
    average_gate_fidelity,
    block_rotation_metrics,
    build_frame,
    build_model,
    channel_waveform_samples,
    compile_pulse_sequence,
    conditioned_targets_from_blocks,
    echoed_half_target_blocks,
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
)


BASELINE_RESULTS_PATH = DATA_DIR / "study_results.json"
BASELINE_ARTIFACTS_DIR = ARTIFACTS_DIR / "cases"

RESULTS_PATH = DATA_DIR / "echo_comparison_results.json"
SUMMARY_PATH = DATA_DIR / "echo_comparison_summary.json"
CSV_PATH = DATA_DIR / "echo_comparison_results.csv"
MARKDOWN_PATH = DATA_DIR / "echo_comparison_summary.md"

ECHO_ARTIFACTS_DIR = ARTIFACTS_DIR / "echo_comparison"
ECHO_CASE_DIR = ECHO_ARTIFACTS_DIR / "cases"
ECHO_WAVEFORM_DIR = ECHO_ARTIFACTS_DIR / "waveforms"
HIGHLIGHTS_DIR = ECHO_ARTIFACTS_DIR / "highlights"

TARGET_FAMILIES = ("C", "D")
ACTIVE_GRID = (2, 3, 4)
DURATION_GRID = (1.0, 3.0, 5.0)
PI_PULSE_DURATION_S = 40.0e-9
PI_PULSE_SIGMA_FRACTION = 0.25
PI_PULSE_PHASE_RAD = 0.0

OBJECTIVE_WEIGHTS = TargetedSubspaceObjectiveWeights(
    qubit_weight=0.15,
    subspace_weight=1.0,
    preservation_weight=0.35,
    leakage_weight=0.35,
)

SEQUENCE_PALETTE = {
    "single_pulse": "#4477AA",
    "echoed_fixed_total": "#EE6677",
    "echoed_fixed_active": "#228833",
}


@dataclass(frozen=True)
class EchoMode:
    name: str
    fixed_total_duration: bool


ECHO_MODES = (
    EchoMode(name="echoed_fixed_total", fixed_total_duration=True),
    EchoMode(name="echoed_fixed_active", fixed_total_duration=False),
)


@dataclass(frozen=True)
class BaselineCase:
    case_id: str
    family: str
    model_variant: str
    include_chi_prime: bool
    n_active: int
    chi_t_over_2pi: float
    pulse_duration_s: float
    pulse_duration_ns: float
    random_seed: int | None


@dataclass
class EchoContext:
    case: BaselineCase
    model: Any
    frame: Any
    levels: tuple[int, ...]
    target_blocks: tuple[np.ndarray, ...]
    target_operator: np.ndarray
    targets: Any
    transfer_set: Any


def optimization_config(n_active: int) -> TargetedSubspaceOptimizationConfig:
    conditioned = ConditionedOptimizationConfig(
        active_levels=tuple(range(int(n_active))),
        parameters=("d_lambda", "d_alpha", "d_omega"),
        maxiter_stage1=18,
        maxiter_stage2=30,
        d_lambda_bounds=(-0.75, 0.75),
        d_alpha_bounds=(-np.pi, np.pi),
        d_omega_hz_bounds=(-3.0e6, 3.0e6),
        regularization_lambda=5.0e-4,
        regularization_alpha=5.0e-4,
        regularization_omega=5.0e-4,
    )
    return TargetedSubspaceOptimizationConfig(conditioned=conditioned, include_block_phase=False)


def _restore_complex_array(value: Any) -> np.ndarray:
    if isinstance(value, dict) and {"real", "imag", "shape"}.issubset(value):
        real = np.asarray(value["real"], dtype=float)
        imag = np.asarray(value["imag"], dtype=float)
        shape = tuple(int(item) for item in value["shape"])
        return (real + 1.0j * imag).reshape(shape)
    return np.asarray(value, dtype=np.complex128)


def _nan_to_none(value: Any) -> Any:
    if isinstance(value, float) and np.isnan(value):
        return None
    return value


def selected_baseline_cases() -> list[BaselineCase]:
    payload = load_json(BASELINE_RESULTS_PATH)
    df = pd.DataFrame(payload["case_rows"])
    subset = df[
        df["family"].isin(TARGET_FAMILIES)
        & df["n_active"].isin(ACTIVE_GRID)
        & df["chi_t_over_2pi"].isin(DURATION_GRID)
    ].copy()
    subset = subset.sort_values(["model_variant", "n_active", "chi_t_over_2pi", "family", "random_seed"])
    cases: list[BaselineCase] = []
    for row in subset.itertuples(index=False):
        cases.append(
            BaselineCase(
                case_id=str(row.case_id),
                family=str(row.family),
                model_variant=str(row.model_variant),
                include_chi_prime=bool(row.include_chi_prime),
                n_active=int(row.n_active),
                chi_t_over_2pi=float(row.chi_t_over_2pi),
                pulse_duration_s=float(row.pulse_duration_s),
                pulse_duration_ns=float(row.pulse_duration_ns),
                random_seed=None if _nan_to_none(row.random_seed) is None else int(row.random_seed),
            )
        )
    return cases


def load_baseline_artifact(case_id: str) -> dict[str, Any]:
    return load_json(BASELINE_ARTIFACTS_DIR / f"{case_id}.json")


def load_baseline_summary_row(case_id: str) -> dict[str, Any]:
    payload = load_json(BASELINE_RESULTS_PATH)
    for row in payload["case_rows"]:
        if str(row["case_id"]) == str(case_id):
            return dict(row)
    raise KeyError(f"Baseline case_id '{case_id}' not found in {BASELINE_RESULTS_PATH}.")


def build_context(case: BaselineCase) -> EchoContext:
    artifact = load_baseline_artifact(case.case_id)
    model = build_model(include_chi_prime=case.include_chi_prime, n_active=case.n_active)
    frame = build_frame(model)
    levels = logical_levels(case.n_active)
    target_operator = _restore_complex_array(artifact["target_operator"])
    target_blocks = restricted_blocks(target_operator)
    targets = conditioned_targets_from_blocks(target_blocks)
    transfer_set = build_spanning_state_transfer_set(target_operator)
    return EchoContext(
        case=case,
        model=model,
        frame=frame,
        levels=levels,
        target_blocks=target_blocks,
        target_operator=target_operator,
        targets=targets,
        transfer_set=transfer_set,
    )


def make_row_from_baseline(case: BaselineCase, artifact: dict[str, Any]) -> dict[str, Any]:
    row = load_baseline_summary_row(case.case_id)
    row.update(
        {
            "sequence_family": "single_pulse",
            "comparison_group": case.case_id,
            "active_sqr_duration_ns": float(case.pulse_duration_ns),
            "total_gate_duration_ns": float(case.pulse_duration_ns),
            "half_sqr_duration_ns": float(case.pulse_duration_ns),
            "pi_pulse_duration_ns": 0.0,
            "half_waveforms_identical": True,
            "second_half_transform": "not_applicable",
            "x_pi_pulse_type": "none",
            "x_pi_axis": "none",
            "x_pi_phase_rad": None,
            "timing_mode": "single",
            "worst_block_average_gate_fidelity": float(np.min(np.asarray(row["per_block_average_gate_fidelities"], dtype=float))),
            "max_rotation_angle_error_rad": float(np.max(np.asarray(row["per_block_rotation_angle_errors_rad"], dtype=float))),
            "max_rotation_axis_error_rad": float(np.nanmax(np.asarray(row["per_block_rotation_axis_errors_rad"], dtype=float))),
            "target_operator_rule": "direct single-pulse target",
            "baseline_artifact": str((BASELINE_ARTIFACTS_DIR / f"{case.case_id}.json").relative_to(STUDY_DIR)),
        }
    )
    return row


def _run_config_for_duration(context: EchoContext, duration_s: float) -> Any:
    return make_run_config(context.model, n_active=context.case.n_active, duration_s=float(duration_s))


def optimize_half_waveform(
    context: EchoContext,
    half_duration_s: float,
    mode_name: str,
) -> tuple[Any, tuple[np.ndarray, ...], np.ndarray, Any]:
    half_blocks = echoed_half_target_blocks(context.target_blocks)
    half_target_operator = np.zeros_like(context.target_operator)
    for index, block in enumerate(half_blocks):
        half_target_operator[2 * index : 2 * index + 2, 2 * index : 2 * index + 2] = block
    half_targets = conditioned_targets_from_blocks(half_blocks)
    half_transfer_set = build_spanning_state_transfer_set(half_target_operator)
    run_config = _run_config_for_duration(context, half_duration_s)
    optimization = optimize_targeted_subspace_multitone(
        context.model,
        half_targets,
        run_config,
        logical_levels=context.levels,
        optimization_config=optimization_config(context.case.n_active),
        objective_weights=OBJECTIVE_WEIGHTS,
        target_operator=half_target_operator,
        transfer_set=half_transfer_set,
        label=f"{mode_name}_{context.case.case_id}",
    )
    return optimization, half_blocks, half_target_operator, run_config


def evaluate_echo_sequence(context: EchoContext, mode: EchoMode) -> tuple[dict[str, Any], dict[str, Any]]:
    total_reference_s = float(context.case.pulse_duration_s)
    if mode.fixed_total_duration:
        active_total_s = total_reference_s - 2.0 * PI_PULSE_DURATION_S
        total_gate_s = total_reference_s
        timing_note = "fixed_total_duration"
    else:
        active_total_s = total_reference_s
        total_gate_s = total_reference_s + 2.0 * PI_PULSE_DURATION_S
        timing_note = "fixed_active_sqr_duration"
    if active_total_s <= 0.0:
        raise ValueError(f"Echo active time is non-positive for {context.case.case_id} in mode {mode.name}.")
    half_duration_s = 0.5 * active_total_s
    optimization, half_blocks, half_target_operator, half_run_config = optimize_half_waveform(context, half_duration_s, mode.name)
    half_waveform = optimization.optimized_result.waveform
    half_pulse = half_waveform.pulse
    x_first = make_gaussian_qubit_rotation_pulse(
        context.model,
        context.frame,
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
        context.model,
        context.frame,
        theta=np.pi,
        phase=PI_PULSE_PHASE_RAD,
        duration_s=PI_PULSE_DURATION_S,
        channel=str(half_pulse.channel),
        manifold_level=0,
        sigma_fraction=PI_PULSE_SIGMA_FRACTION,
        t0=2.0 * half_duration_s + PI_PULSE_DURATION_S,
        label="x_pi_2",
    )
    pulses = [
        shift_pulse(half_pulse, t0=0.0, label="sqr_half_1"),
        x_first,
        shift_pulse(half_pulse, t0=half_duration_s + PI_PULSE_DURATION_S, label="sqr_half_2"),
        x_second,
    ]
    compiled = compile_pulse_sequence(pulses, dt_s=float(half_run_config.dt_s), total_duration_s=total_gate_s)
    full_operator = simulate_full_operator_on_logical_inputs(
        context.model,
        compiled,
        frame=context.frame,
        drive_ops=half_waveform.drive_ops,
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
        metadata={
            "sequence_family": mode.name,
            "timing_mode": timing_note,
            "time_order": "half_sqr -> x_pi -> half_sqr -> x_pi",
            "half_target_rule": "W_n = X_pi^dagger sqrt(R_n)",
        },
    )
    per_block = [
        {"level": int(level), **block_rotation_metrics(target_block, actual_block)}
        for level, target_block, actual_block in zip(
            context.levels,
            context.target_blocks,
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
        drive_ops=half_waveform.drive_ops,
        levels=context.levels,
        target_operator=context.target_operator,
    )
    row = {
        "case_id": context.case.case_id,
        "comparison_group": context.case.case_id,
        "family": context.case.family,
        "model_variant": context.case.model_variant,
        "include_chi_prime": bool(context.case.include_chi_prime),
        "n_active": int(context.case.n_active),
        "chi_t_over_2pi": float(context.case.chi_t_over_2pi),
        "random_seed": context.case.random_seed,
        "sequence_family": mode.name,
        "timing_mode": timing_note,
        "pulse_duration_ns": float(context.case.pulse_duration_ns),
        "active_sqr_duration_ns": float(active_total_s * 1.0e9),
        "total_gate_duration_ns": float(total_gate_s * 1.0e9),
        "half_sqr_duration_ns": float(half_duration_s * 1.0e9),
        "pi_pulse_duration_ns": float(PI_PULSE_DURATION_S * 1.0e9),
        "half_waveforms_identical": True,
        "second_half_transform": "identical_repeat",
        "x_pi_pulse_type": "finite_gaussian",
        "x_pi_axis": "x",
        "x_pi_phase_rad": float(PI_PULSE_PHASE_RAD),
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
        "worst_block_average_gate_fidelity": float(np.min([item["average_gate_fidelity"] for item in per_block])),
        "max_rotation_angle_error_rad": float(np.max([item["rotation_angle_error_rad"] for item in per_block])),
        "max_rotation_axis_error_rad": float(np.nanmax(np.asarray([item["rotation_axis_error_rad"] for item in per_block], dtype=float))),
        "state_validation_ground_fidelity": float(state_summary["states"][0]["state_fidelity"]),
        "state_validation_plus_fidelity": float(state_summary["states"][1]["state_fidelity"]),
        "optimizer_success": bool(getattr(optimization, "success", True)),
        "optimizer_runtime_s": float(sum(item.get("runtime_s", 0.0) for item in optimization.history if isinstance(item, dict))),
        "target_operator_rule": "half_target W_n = X_pi^dagger sqrt(R_n), repeated identically in both halves",
    }
    waveform_samples = channel_waveform_samples(compiled)
    artifact = {
        "study_name": STUDY_DIR.name,
        "date_created": time.strftime("%Y-%m-%d"),
        "description": "Explicit echoed-SQR comparison using the time-ordered pulse schedule half SQR -> X_pi -> half SQR -> X_pi.",
        "sequence_family": mode.name,
        "time_order": "half_sqr -> x_pi -> half_sqr -> x_pi",
        "timing_mode": timing_note,
        "comparison_group": context.case.case_id,
        "single_case_id": context.case.case_id,
        "case_request": json_ready(context.case.__dict__),
        "target_operator": context.target_operator,
        "target_blocks": context.target_blocks,
        "half_target_operator": half_target_operator,
        "half_target_blocks": half_blocks,
        "ideal_x_pi": IDEAL_X_PI,
        "half_optimization": optimization.improvement_summary(),
        "half_optimization_history": optimization.history,
        "half_optimized_corrections": {
            "d_lambda": list(optimization.optimized_corrections.d_lambda),
            "d_alpha": list(optimization.optimized_corrections.d_alpha),
            "d_omega_rad_s": list(optimization.optimized_corrections.d_omega_rad_s),
            "d_omega_hz": [float(value / (2.0 * np.pi)) for value in optimization.optimized_corrections.d_omega_rad_s],
        },
        "half_tone_specs": half_waveform.tone_rows(),
        "half_duration_s": float(half_duration_s),
        "active_sqr_duration_s": float(active_total_s),
        "total_gate_duration_s": float(total_gate_s),
        "pi_pulse": {
            "duration_s": float(PI_PULSE_DURATION_S),
            "sigma_fraction": float(PI_PULSE_SIGMA_FRACTION),
            "phase_rad": float(PI_PULSE_PHASE_RAD),
            "axis": "x",
            "pulse_type": "finite_gaussian",
            "channel": str(half_pulse.channel),
            "manifold_level": 0,
        },
        "sequence_spec": [
            {"kind": "half_sqr", "t0_s": 0.0, "duration_s": float(half_duration_s)},
            {"kind": "x_pi", "t0_s": float(half_duration_s), "duration_s": float(PI_PULSE_DURATION_S)},
            {"kind": "half_sqr", "t0_s": float(half_duration_s + PI_PULSE_DURATION_S), "duration_s": float(half_duration_s)},
            {"kind": "x_pi", "t0_s": float(2.0 * half_duration_s + PI_PULSE_DURATION_S), "duration_s": float(PI_PULSE_DURATION_S)},
        ],
        "restricted_operator": restricted_operator,
        "full_operator_columns_on_logical_inputs": full_operator,
        "active_subspace_metrics": active_subspace_metrics(context.target_operator, restricted_operator),
        "per_block_metrics": per_block,
        "validation": validation.as_dict(),
        "state_validation": state_summary,
        "waveform_samples": waveform_samples,
        "summary_row": row,
        "load_instructions": "Load this JSON and inspect `half_tone_specs`, `pi_pulse`, `sequence_spec`, `restricted_operator`, and `waveform_samples` to reproduce the explicit echoed comparison case.",
    }
    return row, artifact


def save_echo_artifact(mode: EchoMode, case_id: str, artifact: dict[str, Any]) -> Path:
    stem = f"{case_id}_{mode.name}"
    npz_path = ECHO_WAVEFORM_DIR / f"{stem}.npz"
    save_waveform_npz(npz_path, artifact["waveform_samples"])
    artifact["waveform_npz"] = str(npz_path.relative_to(STUDY_DIR))
    path = ECHO_CASE_DIR / f"{stem}.json"
    save_json(path, artifact)
    return path


def _merge_deltas(df: pd.DataFrame) -> pd.DataFrame:
    single = df[df["sequence_family"] == "single_pulse"].copy()
    echoes = df[df["sequence_family"] != "single_pulse"].copy()
    merged = echoes.merge(
        single[[
            "comparison_group",
            "average_gate_fidelity",
            "worst_block_average_gate_fidelity",
            "mean_residual_z_error_rad",
            "mean_transverse_error_rad",
            "state_validation_ground_fidelity",
            "state_validation_plus_fidelity",
        ]],
        on="comparison_group",
        suffixes=("", "_single"),
    )
    merged["delta_fidelity"] = merged["average_gate_fidelity"] - merged["average_gate_fidelity_single"]
    merged["delta_worst_block"] = merged["worst_block_average_gate_fidelity"] - merged["worst_block_average_gate_fidelity_single"]
    merged["delta_residual_z"] = merged["mean_residual_z_error_rad"] - merged["mean_residual_z_error_rad_single"]
    merged["delta_transverse"] = merged["mean_transverse_error_rad"] - merged["mean_transverse_error_rad_single"]
    merged["delta_ground_state"] = merged["state_validation_ground_fidelity"] - merged["state_validation_ground_fidelity_single"]
    merged["delta_plus_state"] = merged["state_validation_plus_fidelity"] - merged["state_validation_plus_fidelity_single"]
    return merged


def build_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    df = pd.DataFrame(rows)
    deltas = _merge_deltas(df)
    branch_summary = (
        df.groupby("sequence_family", as_index=False)
        .agg(
            average_gate_fidelity_mean=("average_gate_fidelity", "mean"),
            average_gate_fidelity_best=("average_gate_fidelity", "max"),
            worst_block_mean=("worst_block_average_gate_fidelity", "mean"),
            residual_z_mean=("mean_residual_z_error_rad", "mean"),
            transverse_mean=("mean_transverse_error_rad", "mean"),
            plus_state_mean=("state_validation_plus_fidelity", "mean"),
        )
        .sort_values("average_gate_fidelity_mean", ascending=False)
    )
    family_branch = (
        df.groupby(["family", "sequence_family"], as_index=False)
        .agg(
            average_gate_fidelity_mean=("average_gate_fidelity", "mean"),
            average_gate_fidelity_median=("average_gate_fidelity", "median"),
            worst_block_median=("worst_block_average_gate_fidelity", "median"),
            residual_z_median=("mean_residual_z_error_rad", "median"),
            transverse_median=("mean_transverse_error_rad", "median"),
            count=("case_id", "count"),
        )
        .sort_values(["family", "average_gate_fidelity_median"], ascending=[True, False])
    )
    duration_summary = (
        df.groupby(["family", "sequence_family", "chi_t_over_2pi"], as_index=False)
        .agg(
            average_gate_fidelity_median=("average_gate_fidelity", "median"),
            residual_z_median=("mean_residual_z_error_rad", "median"),
            transverse_median=("mean_transverse_error_rad", "median"),
        )
        .sort_values(["family", "sequence_family", "chi_t_over_2pi"])
    )
    delta_summary = (
        deltas.groupby(["family", "sequence_family"], as_index=False)
        .agg(
            delta_fidelity_mean=("delta_fidelity", "mean"),
            delta_fidelity_median=("delta_fidelity", "median"),
            delta_worst_block_median=("delta_worst_block", "median"),
            delta_residual_z_median=("delta_residual_z", "median"),
            delta_transverse_median=("delta_transverse", "median"),
            improved_fidelity_count=("delta_fidelity", lambda s: int(np.sum(np.asarray(s, dtype=float) > 0.0))),
            improved_residual_z_count=("delta_residual_z", lambda s: int(np.sum(np.asarray(s, dtype=float) < 0.0))),
            count=("delta_fidelity", "count"),
        )
        .sort_values(["family", "sequence_family"])
    )
    best_single = df[df["sequence_family"] == "single_pulse"].sort_values("average_gate_fidelity", ascending=False).iloc[0].to_dict()
    best_echo = df[df["sequence_family"] != "single_pulse"].sort_values("average_gate_fidelity", ascending=False).iloc[0].to_dict()
    random_deltas = deltas[deltas["family"] == "D"].sort_values("delta_fidelity", ascending=False)
    strongest_random_improvement = None if random_deltas.empty else random_deltas.iloc[0].to_dict()
    strongest_phase_cancellation = None if deltas.empty else deltas.sort_values("delta_residual_z", ascending=True).iloc[0].to_dict()
    return {
        "study": STUDY_DIR.name,
        "n_rows": int(len(rows)),
        "n_single_rows": int(np.sum(df["sequence_family"] == "single_pulse")),
        "n_echo_rows": int(np.sum(df["sequence_family"] != "single_pulse")),
        "pi_pulse_duration_ns": float(PI_PULSE_DURATION_S * 1.0e9),
        "best_single": best_single,
        "best_echo": best_echo,
        "strongest_random_improvement": strongest_random_improvement,
        "strongest_phase_cancellation": strongest_phase_cancellation,
        "branch_summary": branch_summary.to_dict(orient="records"),
        "family_branch_summary": family_branch.to_dict(orient="records"),
        "duration_summary": duration_summary.to_dict(orient="records"),
        "delta_summary": delta_summary.to_dict(orient="records"),
        "rows": rows,
    }


def save_figure(fig, stem: str) -> None:
    fig.savefig(FIGURES_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def _extract_waveform_components(samples: dict[str, Any]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    time_s = np.asarray(samples["time_s"], dtype=float)
    if "baseband_real" in samples:
        real = np.asarray(samples["baseband_real"], dtype=float)
        imag = np.asarray(samples["baseband_imag"], dtype=float)
    else:
        real = np.asarray(samples["real"], dtype=float)
        imag = np.asarray(samples["imag"], dtype=float)
    return time_s, real, imag


def _numeric_vector(value: Any) -> np.ndarray:
    arr = np.asarray(value, dtype=object)
    if arr.shape == ():
        inner = arr.item()
        return np.asarray(inner, dtype=float).reshape(-1)
    if arr.shape == (1,) and isinstance(arr[0], (list, tuple, np.ndarray)):
        return np.asarray(arr[0], dtype=float).reshape(-1)
    return np.asarray(value, dtype=float).reshape(-1)


def plot_figures(df: pd.DataFrame, summary: dict[str, Any]) -> None:
    apply_plot_style()

    branch_order = ["single_pulse", "echoed_fixed_total", "echoed_fixed_active"]
    branch_metrics = (
        df.groupby("sequence_family")
        .agg(
            avg_fidelity=("average_gate_fidelity", "mean"),
            worst_block=("worst_block_average_gate_fidelity", "mean"),
            residual_z=("mean_residual_z_error_rad", "mean"),
            transverse=("mean_transverse_error_rad", "mean"),
        )
        .reindex(branch_order)
    )
    fig, axes = plt.subplots(1, 4, figsize=(16.0, 4.0))
    for axis, column, ylabel in zip(
        axes,
        ("avg_fidelity", "worst_block", "residual_z", "transverse"),
        ("Mean average fidelity", "Mean worst-block fidelity", "Mean residual Z error (rad)", "Mean transverse error (rad)"),
        strict=True,
    ):
        axis.bar(branch_metrics.index, branch_metrics[column], color=[SEQUENCE_PALETTE[item] for item in branch_metrics.index])
        axis.set_ylabel(ylabel)
        axis.tick_params(axis="x", rotation=25)
    save_figure(fig, "echo_branch_metric_means")

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.5), sharex="col")
    for row_index, family in enumerate(("C", "D")):
        subset = df[df["family"] == family]
        for branch in branch_order:
            branch_df = subset[subset["sequence_family"] == branch]
            grouped = branch_df.groupby("chi_t_over_2pi", as_index=False).median(numeric_only=True)
            axes[row_index, 0].plot(
                grouped["chi_t_over_2pi"],
                grouped["average_gate_fidelity"],
                marker="o",
                label=branch,
                color=SEQUENCE_PALETTE[branch],
            )
            axes[row_index, 1].plot(
                grouped["chi_t_over_2pi"],
                grouped["mean_residual_z_error_rad"],
                marker="o",
                label=branch,
                color=SEQUENCE_PALETTE[branch],
            )
        axes[row_index, 0].set_ylabel(f"Family {family} median fidelity")
        axes[row_index, 1].set_ylabel(f"Family {family} median residual Z (rad)")
    axes[0, 0].legend(frameon=False)
    axes[1, 0].set_xlabel(r"$|\chi| T / 2\pi$")
    axes[1, 1].set_xlabel(r"$|\chi| T / 2\pi$")
    save_figure(fig, "echo_duration_scan")

    deltas = _merge_deltas(df)
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.5), sharey=True)
    markers = {"C": "o", "D": "s"}
    for family in ("C", "D"):
        family_df = deltas[deltas["family"] == family]
        for branch in ("echoed_fixed_total", "echoed_fixed_active"):
            branch_df = family_df[family_df["sequence_family"] == branch]
            axes[0].scatter(
                branch_df["delta_residual_z"],
                branch_df["delta_fidelity"],
                color=SEQUENCE_PALETTE[branch],
                marker=markers[family],
                alpha=0.8,
                label=f"{branch} | family {family}",
            )
            axes[1].scatter(
                branch_df["delta_transverse"],
                branch_df["delta_fidelity"],
                color=SEQUENCE_PALETTE[branch],
                marker=markers[family],
                alpha=0.8,
                label=f"{branch} | family {family}",
            )
    axes[0].axvline(0.0, color="#666666", linewidth=1.0, linestyle="--")
    axes[1].axvline(0.0, color="#666666", linewidth=1.0, linestyle="--")
    axes[0].axhline(0.0, color="#666666", linewidth=1.0, linestyle="--")
    axes[1].axhline(0.0, color="#666666", linewidth=1.0, linestyle="--")
    axes[0].set_xlabel(r"$\Delta$ mean residual Z error (rad)")
    axes[1].set_xlabel(r"$\Delta$ mean transverse error (rad)")
    axes[0].set_ylabel(r"$\Delta$ average gate fidelity")
    axes[0].legend(frameon=False, fontsize=8)
    save_figure(fig, "echo_delta_tradeoff")

    best_single_artifact = load_baseline_artifact(str(summary["best_single"]["case_id"]))
    best_echo_artifact = load_json(ECHO_CASE_DIR / f"{summary['best_echo']['case_id']}_{summary['best_echo']['sequence_family']}.json")
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 5.8), sharex="col")
    for column, (title, artifact) in enumerate((("Best single pulse", best_single_artifact), ("Best echoed sequence", best_echo_artifact))):
        time_s, real, imag = _extract_waveform_components(artifact["waveform_samples"])
        time_ns = 1.0e9 * time_s
        axes[0, column].plot(time_ns, real, label="I", color="#4477AA")
        axes[0, column].plot(time_ns, imag, label="Q", color="#EE6677")
        axes[0, column].set_title(title)
        axes[0, column].set_ylabel("Drive amplitude (rad/s)")
        axes[0, column].legend(frameon=False)
        signal = real + 1.0j * imag
        dt_s = float(np.mean(np.diff(time_s)))
        freqs_hz = np.fft.fftfreq(signal.size, d=dt_s)
        positive = freqs_hz >= 0.0
        spectrum = np.abs(np.fft.fft(signal))[positive]
        if np.max(spectrum) > 0.0:
            spectrum = spectrum / np.max(spectrum)
        axes[1, column].plot(freqs_hz[positive] / 1.0e6, spectrum, color="#228833")
        axes[1, column].set_xlabel("Baseband frequency (MHz)")
        axes[1, column].set_ylabel("Normalized spectrum")
        axes[1, column].set_xlim(0.0, 125.0)
    save_figure(fig, "echo_best_waveforms")

    focus = _merge_deltas(df).sort_values("delta_fidelity", ascending=False).iloc[0]
    single_focus_row = df[(df["comparison_group"] == focus["comparison_group"]) & (df["sequence_family"] == "single_pulse")].iloc[0]
    echo_focus_artifact = load_json(ECHO_CASE_DIR / f"{focus['comparison_group']}_{focus['sequence_family']}.json")
    echo_block_metrics = echo_focus_artifact["per_block_metrics"]
    single_residual = _numeric_vector(single_focus_row["per_block_residual_z_errors_rad"])
    single_transverse = _numeric_vector(single_focus_row["per_block_transverse_errors_rad"])
    levels = np.arange(single_residual.size)
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3), sharex=True)
    echo_residual = np.asarray([float(row["residual_z_error_rad"]) for row in echo_block_metrics], dtype=float)
    echo_transverse = np.asarray([float(row["transverse_error_rad"]) for row in echo_block_metrics], dtype=float)
    axes[0].plot(levels, single_residual, marker="o", label="single pulse", color=SEQUENCE_PALETTE["single_pulse"])
    axes[0].plot(levels, echo_residual, marker="s", label=str(focus["sequence_family"]), color=SEQUENCE_PALETTE[str(focus["sequence_family"])])
    axes[0].set_xlabel("Fock block n")
    axes[0].set_ylabel("Residual Z error (rad)")
    axes[0].legend(frameon=False)
    axes[1].plot(levels, single_transverse, marker="o", label="single pulse", color=SEQUENCE_PALETTE["single_pulse"])
    axes[1].plot(levels, echo_transverse, marker="s", label=str(focus["sequence_family"]), color=SEQUENCE_PALETTE[str(focus["sequence_family"])])
    axes[1].set_xlabel("Fock block n")
    axes[1].set_ylabel("Transverse coherent error (rad)")
    axes[1].legend(frameon=False)
    save_figure(fig, "echo_block_error_breakdown")


def write_markdown_summary(summary: dict[str, Any]) -> None:
    lines = [
        f"# Echo Comparison Summary: {summary['study']}",
        "",
        "## What changed relative to the previous report",
        "",
        "- The previous report tested only the single Gaussian multitone SQR ansatz.",
        "- This extension adds an explicit echoed sequence with the time-ordered schedule: half SQR -> X_pi -> half SQR -> X_pi.",
        f"- The inserted X_pi pulses are finite Gaussian pulses of duration {summary['pi_pulse_duration_ns']:.1f} ns about the x axis.",
        "- Two fairness conventions are included: fixed total gate duration and fixed active SQR duration.",
        "",
        "## Headline results",
        "",
        f"- Best single pulse: {summary['best_single']['case_id']} with fidelity {summary['best_single']['average_gate_fidelity']:.6f}.",
        f"- Best echoed case: {summary['best_echo']['sequence_family']} on {summary['best_echo']['case_id']} with fidelity {summary['best_echo']['average_gate_fidelity']:.6f}.",
        (
            f"- Strongest random-target fidelity gain: {summary['strongest_random_improvement']['sequence_family']} on {summary['strongest_random_improvement']['comparison_group']} with delta fidelity {summary['strongest_random_improvement']['delta_fidelity']:+.6f}."
            if summary["strongest_random_improvement"] is not None
            else "- Strongest random-target fidelity gain: not available in this filtered run."
        ),
        (
            f"- Strongest residual-Z reduction: {summary['strongest_phase_cancellation']['sequence_family']} on {summary['strongest_phase_cancellation']['comparison_group']} with delta residual Z {summary['strongest_phase_cancellation']['delta_residual_z']:+.6f} rad."
            if summary["strongest_phase_cancellation"] is not None
            else "- Strongest residual-Z reduction: not available in this filtered run."
        ),
        "",
        "## Branch means",
        "",
    ]
    for row in summary["branch_summary"]:
        lines.append(
            "- {sequence_family}: mean fidelity {average_gate_fidelity_mean:.6f}, best fidelity {average_gate_fidelity_best:.6f}, mean worst-block fidelity {worst_block_mean:.6f}, mean residual Z {residual_z_mean:.6f} rad, mean transverse {transverse_mean:.6f} rad".format(
                **row,
            )
        )
    lines.append("")
    lines.append("## Echo minus single medians")
    lines.append("")
    for row in summary["delta_summary"]:
        lines.append(
            "- Family {family}, {sequence_family}: delta fidelity median {delta_fidelity_median:+.6f}, delta worst-block median {delta_worst_block_median:+.6f}, delta residual Z median {delta_residual_z_median:+.6f} rad, delta transverse median {delta_transverse_median:+.6f} rad, improved-fidelity count {improved_fidelity_count}/{count}, reduced-residual-Z count {improved_residual_z_count}/{count}".format(
                **row,
            )
        )
    MARKDOWN_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def save_highlights(summary: dict[str, Any]) -> None:
    best_single_case_id = str(summary["best_single"]["case_id"])
    best_echo_case_id = str(summary["best_echo"]["case_id"])
    best_echo_branch = str(summary["best_echo"]["sequence_family"])
    best_single_artifact = load_baseline_artifact(best_single_case_id)
    best_echo_artifact = load_json(ECHO_CASE_DIR / f"{best_echo_case_id}_{best_echo_branch}.json")
    save_json(HIGHLIGHTS_DIR / "best_single_pulse.json", best_single_artifact)
    save_json(HIGHLIGHTS_DIR / "best_echoed_sqr.json", best_echo_artifact)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=str, default="", help="Optional baseline case_id to run.")
    args = parser.parse_args()

    rows: list[dict[str, Any]] = []
    cases = selected_baseline_cases()
    if args.case:
        cases = [case for case in cases if case.case_id == args.case]
        if not cases:
            raise SystemExit(f"No selected baseline case matched '{args.case}'.")

    for index, case in enumerate(cases, start=1):
        print(f"[{index:03d}/{len(cases):03d}] Preparing {case.case_id}", flush=True)
        baseline_artifact = load_baseline_artifact(case.case_id)
        rows.append(make_row_from_baseline(case, baseline_artifact))
        context = build_context(case)
        for mode in ECHO_MODES:
            start = time.perf_counter()
            row, artifact = evaluate_echo_sequence(context, mode)
            save_echo_artifact(mode, case.case_id, artifact)
            rows.append(row)
            print(
                f"     {mode.name}: fidelity={row['average_gate_fidelity']:.6f} "
                f"resZ={row['mean_residual_z_error_rad']:.6f} "
                f"trans={row['mean_transverse_error_rad']:.6f} "
                f"elapsed={time.perf_counter() - start:.1f}s",
                flush=True,
            )

    save_json(RESULTS_PATH, {"study": STUDY_DIR.name, "rows": rows})
    df = pd.DataFrame(rows)
    df.to_csv(CSV_PATH, index=False)
    summary = build_summary(rows)
    save_json(SUMMARY_PATH, summary)
    write_markdown_summary(summary)
    plot_figures(df, summary)
    save_highlights(summary)
    print(f"Saved {RESULTS_PATH}")
    print(f"Saved {SUMMARY_PATH}")


if __name__ == "__main__":
    main()