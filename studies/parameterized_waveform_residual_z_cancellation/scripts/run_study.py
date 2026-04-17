"""Pilot study runner for residual-Z cancellation with richer multitone families."""

from __future__ import annotations

import argparse
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import Bounds, minimize

from cqed_sim.calibration.conditioned_multitone import ConditionedOptimizationConfig, build_conditioned_multitone_waveform
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
    STUDY_DIR,
    apply_plot_style,
    average_gate_fidelity,
    block_diag_target,
    block_rotation_metrics,
    build_frame,
    build_model,
    channel_waveform_samples,
    compile_pulse_sequence,
    conditioned_targets_from_blocks,
    duration_from_chi_t,
    frobenius_error,
    gaussian_samples,
    json_ready,
    load_json,
    logical_levels,
    make_gaussian_qubit_rotation_pulse,
    make_run_config,
    make_target_blocks,
    normalize_complex_samples,
    orthogonal_basis,
    restricted_blocks,
    restricted_operator_from_full,
    save_json,
    save_waveform_npz,
    shift_pulse,
    simulate_full_operator_on_logical_inputs,
    single_pulse_sequence,
    state_validation_summary_for_sequence,
)


RESULTS_PATH = DATA_DIR / "study_results.json"
SUMMARY_PATH = DATA_DIR / "study_summary.json"
CSV_PATH = DATA_DIR / "study_results.csv"
MARKDOWN_SUMMARY_PATH = DATA_DIR / "study_summary.md"

MODEL_VARIANTS = (
    ("chi_only", False),
    ("chi_plus_chiprime", True),
)
STUDY_DURATIONS = (3.0, 5.0)
STUDY_ACTIVE_GRID = (2, 3, 4)
STRUCTURED_TARGETS = ("C",)
RANDOM_ENSEMBLE_SIZE = 3
BASE_RANDOM_SEED = 90210

WAVEFORM_FAMILIES = (
    "baseline_multitone",
    "symmetric_two_segment",
    "echoed_multitone",
    "complex_envelope",
    "basis_expanded",
)

OBJECTIVE_WEIGHTS = TargetedSubspaceObjectiveWeights(
    qubit_weight=0.15,
    subspace_weight=1.0,
    preservation_weight=0.35,
    leakage_weight=0.35,
)


@dataclass(frozen=True)
class CaseRequest:
    model_variant: str
    include_chi_prime: bool
    n_active: int
    chi_t_over_2pi: float
    target_family: str
    random_seed: int | None = None

    @property
    def case_id(self) -> str:
        duration_label = str(self.chi_t_over_2pi).replace(".", "p")
        seed_part = "" if self.random_seed is None else f"_seed{int(self.random_seed)}"
        return (
            f"{self.model_variant}_na{int(self.n_active)}_chiT{duration_label}_target{self.target_family}{seed_part}"
        )


@dataclass
class CaseContext:
    request: CaseRequest
    model: Any
    frame: Any
    levels: tuple[int, ...]
    duration_s: float
    run_config: Any
    targets: Any
    blocks: tuple[np.ndarray, ...]
    target_operator: np.ndarray
    transfer_set: Any
    baseline_optimization: Any


def baseline_optimization_config(n_active: int) -> TargetedSubspaceOptimizationConfig:
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


def case_requests() -> list[CaseRequest]:
    requests: list[CaseRequest] = []
    for model_variant, include_chi_prime in MODEL_VARIANTS:
        for n_active in STUDY_ACTIVE_GRID:
            for chi_t in STUDY_DURATIONS:
                for target_family in STRUCTURED_TARGETS:
                    requests.append(
                        CaseRequest(
                            model_variant=model_variant,
                            include_chi_prime=include_chi_prime,
                            n_active=n_active,
                            chi_t_over_2pi=chi_t,
                            target_family=target_family,
                            random_seed=None,
                        )
                    )
                for seed_offset in range(RANDOM_ENSEMBLE_SIZE):
                    requests.append(
                        CaseRequest(
                            model_variant=model_variant,
                            include_chi_prime=include_chi_prime,
                            n_active=n_active,
                            chi_t_over_2pi=chi_t,
                            target_family="D",
                            random_seed=(
                                BASE_RANDOM_SEED
                                + 1000 * int(n_active)
                                + 100 * int(include_chi_prime)
                                + 10 * int(round(chi_t))
                                + seed_offset
                            ),
                        )
                    )
    return requests


def build_case_context(request: CaseRequest) -> CaseContext:
    rng = None if request.random_seed is None else np.random.default_rng(int(request.random_seed))
    blocks, _ = make_target_blocks(request.target_family, request.n_active, rng=rng)
    target_operator = block_diag_target(blocks)
    targets = conditioned_targets_from_blocks(blocks)
    model = build_model(include_chi_prime=request.include_chi_prime, n_active=request.n_active)
    frame = build_frame(model)
    duration_s = duration_from_chi_t(request.chi_t_over_2pi)
    run_config = make_run_config(model, n_active=request.n_active, duration_s=duration_s)
    transfer_set = build_spanning_state_transfer_set(target_operator)
    baseline = optimize_targeted_subspace_multitone(
        model,
        targets,
        run_config,
        logical_levels=logical_levels(request.n_active),
        optimization_config=baseline_optimization_config(request.n_active),
        objective_weights=OBJECTIVE_WEIGHTS,
        target_operator=target_operator,
        transfer_set=transfer_set,
        label=f"baseline_{request.case_id}",
    )
    return CaseContext(
        request=request,
        model=model,
        frame=frame,
        levels=logical_levels(request.n_active),
        duration_s=duration_s,
        run_config=run_config,
        targets=targets,
        blocks=blocks,
        target_operator=target_operator,
        transfer_set=transfer_set,
        baseline_optimization=baseline,
    )


def family_bounds(family: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if family == "symmetric_two_segment":
        x0 = np.asarray([0.0, 0.0, 0.0, 0.0, 1.0, 0.0], dtype=float)
        lower = np.asarray([-0.9, -0.9, -1.2, -1.2, 0.2, -np.pi], dtype=float)
        upper = np.asarray([0.9, 0.9, 1.2, 1.2, 1.6, np.pi], dtype=float)
        return x0, lower, upper
    if family == "echoed_multitone":
        x0 = np.asarray([0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0], dtype=float)
        lower = np.asarray([-0.8, -0.8, -1.2, -1.2, 0.4, -np.pi, -np.pi, 0.6], dtype=float)
        upper = np.asarray([0.8, 0.8, 1.2, 1.2, 1.6, np.pi, np.pi, 1.4], dtype=float)
        return x0, lower, upper
    if family == "complex_envelope":
        x0 = np.zeros(6, dtype=float)
        lower = np.asarray([-0.9, -0.9, -0.9, -1.4, -1.4, -1.4], dtype=float)
        upper = np.asarray([0.9, 0.9, 0.9, 1.4, 1.4, 1.4], dtype=float)
        return x0, lower, upper
    if family == "basis_expanded":
        x0 = np.zeros(8, dtype=float)
        lower = np.asarray([-0.8] * 4 + [-1.0] * 4, dtype=float)
        upper = np.asarray([0.8] * 4 + [1.0] * 4, dtype=float)
        return x0, lower, upper
    raise ValueError(f"Unsupported family '{family}'.")


def _single_segment_run_config(run_config: Any, duration_s: float) -> Any:
    return type(run_config)(
        frame=run_config.frame,
        duration_s=float(duration_s),
        dt_s=float(run_config.dt_s),
        sigma_fraction=float(run_config.sigma_fraction),
        tone_cutoff=float(run_config.tone_cutoff),
        include_all_levels=bool(run_config.include_all_levels),
        max_step_s=run_config.max_step_s,
        fock_fqs_hz=run_config.fock_fqs_hz,
    )


def build_symmetric_two_segment_sequence(context: CaseContext, params: np.ndarray) -> tuple[list[Any], dict[str, str], dict[str, Any]]:
    half_duration = 0.5 * float(context.duration_s)
    count = max(16, int(round(half_duration / float(context.run_config.dt_s))))
    base = gaussian_samples(count, sigma_fraction=float(context.run_config.sigma_fraction))
    cos_basis = orthogonal_basis(count, 2, kind="cos")
    amp_profile = 1.0 + float(params[0]) * cos_basis[0] + float(params[1]) * cos_basis[1]
    phase_profile = float(params[2]) * cos_basis[0] + float(params[3]) * cos_basis[1]
    samples_1 = normalize_complex_samples(base * amp_profile * np.exp(1.0j * phase_profile))
    samples_2 = normalize_complex_samples(float(params[4]) * samples_1[::-1].conj() * np.exp(1.0j * float(params[5])))
    run_half = _single_segment_run_config(context.run_config, half_duration)
    tone_specs = context.baseline_optimization.optimized_result.waveform.tone_specs
    wave_1 = build_conditioned_multitone_waveform(
        tone_specs,
        run_half,
        base_samples=samples_1,
        sample_rate=1.0 / float(run_half.dt_s),
        label="symmetric_seg1",
    )
    wave_2 = build_conditioned_multitone_waveform(
        tone_specs,
        run_half,
        base_samples=samples_2,
        sample_rate=1.0 / float(run_half.dt_s),
        label="symmetric_seg2",
    )
    pulses = [
        shift_pulse(wave_1.pulse, t0=0.0, label="symmetric_seg1"),
        shift_pulse(wave_2.pulse, t0=half_duration, label="symmetric_seg2"),
    ]
    return pulses, wave_1.drive_ops, {
        "segment_count": 2,
        "segment_duration_s": half_duration,
        "parameters": params.tolist(),
    }


def build_echoed_multitone_sequence(context: CaseContext, params: np.ndarray) -> tuple[list[Any], dict[str, str], dict[str, Any]]:
    echo_duration_s = float(np.clip(40.0e-9 * float(params[7]), 24.0e-9, min(80.0e-9, 0.2 * float(context.duration_s))))
    segment_duration_s = 0.5 * float(context.duration_s - echo_duration_s)
    count = max(16, int(round(segment_duration_s / float(context.run_config.dt_s))))
    base = gaussian_samples(count, sigma_fraction=float(context.run_config.sigma_fraction))
    cos_basis = orthogonal_basis(count, 2, kind="cos")
    amp_profile = 1.0 + float(params[0]) * cos_basis[0] + float(params[1]) * cos_basis[1]
    phase_profile = float(params[2]) * cos_basis[0] + float(params[3]) * cos_basis[1]
    samples_1 = normalize_complex_samples(base * amp_profile * np.exp(1.0j * phase_profile))
    samples_2 = normalize_complex_samples(float(params[4]) * samples_1[::-1].conj() * np.exp(1.0j * float(params[5])))
    run_half = _single_segment_run_config(context.run_config, segment_duration_s)
    tone_specs = context.baseline_optimization.optimized_result.waveform.tone_specs
    wave_1 = build_conditioned_multitone_waveform(
        tone_specs,
        run_half,
        base_samples=samples_1,
        sample_rate=1.0 / float(run_half.dt_s),
        label="echoed_seg1",
    )
    wave_2 = build_conditioned_multitone_waveform(
        tone_specs,
        run_half,
        base_samples=samples_2,
        sample_rate=1.0 / float(run_half.dt_s),
        label="echoed_seg2",
    )
    echo_pulse = make_gaussian_qubit_rotation_pulse(
        context.model,
        context.frame,
        theta=np.pi,
        phase=float(params[6]),
        duration_s=echo_duration_s,
        channel=str(wave_1.pulse.channel),
        manifold_level=0,
        sigma_fraction=0.25,
        t0=segment_duration_s,
        label="echo_pi",
    )
    pulses = [
        shift_pulse(wave_1.pulse, t0=0.0, label="echoed_seg1"),
        echo_pulse,
        shift_pulse(wave_2.pulse, t0=segment_duration_s + echo_duration_s, label="echoed_seg2"),
    ]
    return pulses, wave_1.drive_ops, {
        "segment_count": 3,
        "segment_duration_s": segment_duration_s,
        "echo_duration_s": echo_duration_s,
        "parameters": params.tolist(),
    }


def build_complex_envelope_sequence(context: CaseContext, params: np.ndarray) -> tuple[list[Any], dict[str, str], dict[str, Any]]:
    count = max(32, int(round(float(context.duration_s) / float(context.run_config.dt_s))))
    base = gaussian_samples(count, sigma_fraction=float(context.run_config.sigma_fraction))
    real_basis = orthogonal_basis(count, 2, kind="cos") + orthogonal_basis(count, 1, kind="sin")
    phase_basis = orthogonal_basis(count, 1, kind="cos") + orthogonal_basis(count, 2, kind="sin")
    amp_profile = np.ones(count, dtype=float)
    phase_profile = np.zeros(count, dtype=float)
    for coeff, basis_vec in zip(params[:3], real_basis, strict=True):
        amp_profile = amp_profile + float(coeff) * basis_vec
    for coeff, basis_vec in zip(params[3:], phase_basis, strict=True):
        phase_profile = phase_profile + float(coeff) * basis_vec
    samples = normalize_complex_samples(base * amp_profile * np.exp(1.0j * phase_profile))
    pulses, drive_ops = single_pulse_sequence(
        context.baseline_optimization.optimized_result.waveform.tone_specs,
        context.run_config,
        base_samples=samples,
        label="complex_envelope",
    )
    return pulses, drive_ops, {"segment_count": 1, "parameters": params.tolist()}


def build_basis_expanded_sequence(context: CaseContext, params: np.ndarray) -> tuple[list[Any], dict[str, str], dict[str, Any]]:
    count = max(32, int(round(float(context.duration_s) / float(context.run_config.dt_s))))
    base = gaussian_samples(count, sigma_fraction=float(context.run_config.sigma_fraction))
    cos_basis = orthogonal_basis(count, 4, kind="cos")
    amp_profile = np.ones(count, dtype=float)
    phase_profile = np.zeros(count, dtype=float)
    for coeff, basis_vec in zip(params[:4], cos_basis, strict=True):
        amp_profile = amp_profile + float(coeff) * basis_vec
    for coeff, basis_vec in zip(params[4:], cos_basis, strict=True):
        phase_profile = phase_profile + float(coeff) * basis_vec
    samples = normalize_complex_samples(base * amp_profile * np.exp(1.0j * phase_profile))
    pulses, drive_ops = single_pulse_sequence(
        context.baseline_optimization.optimized_result.waveform.tone_specs,
        context.run_config,
        base_samples=samples,
        label="basis_expanded",
    )
    return pulses, drive_ops, {"segment_count": 1, "parameters": params.tolist()}


FAMILY_BUILDERS: dict[str, Callable[[CaseContext, np.ndarray], tuple[list[Any], dict[str, str], dict[str, Any]]]] = {
    "symmetric_two_segment": build_symmetric_two_segment_sequence,
    "echoed_multitone": build_echoed_multitone_sequence,
    "complex_envelope": build_complex_envelope_sequence,
    "basis_expanded": build_basis_expanded_sequence,
}


def evaluate_sequence(
    context: CaseContext,
    waveform_family: str,
    pulses: Sequence[Any],
    drive_ops: dict[str, str],
    metadata: dict[str, Any],
) -> dict[str, Any]:
    total_duration_s = float(max(float(pulse.t1) for pulse in pulses))
    compiled = compile_pulse_sequence(pulses, dt_s=float(context.run_config.dt_s), total_duration_s=total_duration_s)
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
        metadata=metadata,
    )
    per_block = []
    for level, target_block, actual_block in zip(
        context.levels,
        context.blocks,
        restricted_blocks(restricted_operator),
        strict=True,
    ):
        row = {"level": int(level), **block_rotation_metrics(target_block, actual_block)}
        per_block.append(row)
    residual_z = np.asarray([float(row["residual_z_error_rad"]) for row in per_block], dtype=float)
    transverse = np.asarray([float(row["transverse_error_rad"]) for row in per_block], dtype=float)
    state_summary = state_validation_summary_for_sequence(
        context.model,
        compiled,
        frame=context.frame,
        drive_ops=drive_ops,
        levels=context.levels,
        target_operator=context.target_operator,
    )
    return {
        "compiled": compiled,
        "full_operator": full_operator,
        "restricted_operator": restricted_operator,
        "validation": validation,
        "per_block": per_block,
        "mean_residual_z_error_rad": float(np.mean(residual_z)),
        "max_residual_z_error_rad": float(np.max(residual_z)),
        "mean_transverse_error_rad": float(np.mean(transverse)),
        "max_transverse_error_rad": float(np.max(transverse)),
        "state_summary": state_summary,
        "waveform_samples": channel_waveform_samples(compiled),
        "waveform_family": waveform_family,
    }


def family_objective(context: CaseContext, family: str, params: np.ndarray) -> float:
    pulses, drive_ops, metadata = FAMILY_BUILDERS[family](context, np.asarray(params, dtype=float))
    evaluation = evaluate_sequence(context, family, pulses, drive_ops, metadata)
    validation = evaluation["validation"]
    regularization = 1.0e-3 * float(np.linalg.norm(np.asarray(params, dtype=float)) ** 2)
    return float(validation.weighted_loss) + 0.20 * float(evaluation["mean_residual_z_error_rad"]) + regularization


def optimize_custom_family(context: CaseContext, family: str) -> tuple[dict[str, Any], Any]:
    x0, lower, upper = family_bounds(family)
    start = time.perf_counter()
    result = minimize(
        lambda x: family_objective(context, family, x),
        x0,
        method="Powell",
        bounds=Bounds(lower, upper),
        options={"maxiter": 10, "maxfev": 120, "xtol": 5.0e-3, "ftol": 5.0e-3},
    )
    runtime_s = time.perf_counter() - start
    pulses, drive_ops, metadata = FAMILY_BUILDERS[family](context, np.asarray(result.x, dtype=float))
    evaluation = evaluate_sequence(context, family, pulses, drive_ops, metadata)
    evaluation["optimizer"] = {
        "success": bool(result.success),
        "message": str(result.message),
        "nfev": int(getattr(result, "nfev", -1)),
        "nit": int(getattr(result, "nit", -1)) if getattr(result, "nit", None) is not None else -1,
        "fun": float(result.fun),
        "x": np.asarray(result.x, dtype=float).tolist(),
        "runtime_s": float(runtime_s),
    }
    return evaluation, result


def baseline_evaluation(context: CaseContext) -> dict[str, Any]:
    baseline_waveform = context.baseline_optimization.optimized_result.waveform
    pulses = [baseline_waveform.pulse]
    return evaluate_sequence(
        context,
        "baseline_multitone",
        pulses,
        baseline_waveform.drive_ops,
        {
            "segment_count": 1,
            "tone_specs": baseline_waveform.tone_rows(),
            "optimization": context.baseline_optimization.improvement_summary(),
        },
    )


def row_from_evaluation(context: CaseContext, evaluation: dict[str, Any]) -> dict[str, Any]:
    validation = evaluation["validation"]
    state_summary = evaluation["state_summary"]
    optimizer = evaluation.get("optimizer", {})
    return {
        "case_id": context.request.case_id,
        "model_variant": context.request.model_variant,
        "include_chi_prime": bool(context.request.include_chi_prime),
        "n_active": int(context.request.n_active),
        "chi_t_over_2pi": float(context.request.chi_t_over_2pi),
        "target_family": context.request.target_family,
        "random_seed": context.request.random_seed,
        "waveform_family": str(evaluation["waveform_family"]),
        "pulse_duration_ns": float(context.duration_s * 1.0e9),
        "average_gate_fidelity": float(average_gate_fidelity(context.target_operator, evaluation["restricted_operator"])),
        "restricted_process_fidelity": float(validation.restricted_process_fidelity),
        "restricted_fro_error": float(frobenius_error(context.target_operator, evaluation["restricted_operator"])),
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
        "mean_residual_z_error_rad": float(evaluation["mean_residual_z_error_rad"]),
        "max_residual_z_error_rad": float(evaluation["max_residual_z_error_rad"]),
        "mean_transverse_error_rad": float(evaluation["mean_transverse_error_rad"]),
        "max_transverse_error_rad": float(evaluation["max_transverse_error_rad"]),
        "state_validation_ground_fidelity": float(state_summary["states"][0]["state_fidelity"]),
        "state_validation_plus_fidelity": float(state_summary["states"][1]["state_fidelity"]),
        "optimizer_success": None if not optimizer else bool(optimizer.get("success", False)),
        "optimizer_nfev": None if not optimizer else int(optimizer.get("nfev", -1)),
        "optimizer_runtime_s": None if not optimizer else float(optimizer.get("runtime_s", 0.0)),
    }


def save_case_artifact(context: CaseContext, evaluation: dict[str, Any], row: dict[str, Any]) -> None:
    stem = f"{context.request.case_id}_{evaluation['waveform_family']}"
    npz_path = ARTIFACTS_DIR / "waveforms" / f"{stem}.npz"
    save_waveform_npz(npz_path, evaluation["waveform_samples"])
    artifact_payload = {
        "study_name": STUDY_DIR.name,
        "date_created": time.strftime("%Y-%m-%d"),
        "description": "Residual-Z cancellation comparison for richer multitone waveform families across structured and random conditional-rotation targets.",
        "case_request": json_ready(context.request.__dict__),
        "waveform_family": evaluation["waveform_family"],
        "summary_row": row,
        "target_operator": context.target_operator,
        "restricted_operator": evaluation["restricted_operator"],
        "full_operator_columns_on_logical_inputs": evaluation["full_operator"],
        "per_block_metrics": evaluation["per_block"],
        "validation": evaluation["validation"].as_dict(),
        "state_validation": evaluation["state_summary"],
        "waveform_samples": evaluation["waveform_samples"],
        "waveform_npz": str(npz_path.relative_to(STUDY_DIR)),
        "optimizer": evaluation.get("optimizer"),
        "load_instructions": "Load this JSON and inspect `waveform_samples`, `waveform_npz`, `per_block_metrics`, and `validation` to reproduce the reported comparison case.",
    }
    save_json(ARTIFACTS_DIR / "cases" / f"{stem}.json", artifact_payload)


def run_case(context: CaseContext) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    baseline_eval = baseline_evaluation(context)
    baseline_row = row_from_evaluation(context, baseline_eval)
    save_case_artifact(context, baseline_eval, baseline_row)
    rows.append(baseline_row)
    for family in WAVEFORM_FAMILIES[1:]:
        evaluation, _ = optimize_custom_family(context, family)
        row = row_from_evaluation(context, evaluation)
        save_case_artifact(context, evaluation, row)
        rows.append(row)
    return rows


def build_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    df = pd.DataFrame(rows)
    family_group = (
        df.groupby("waveform_family", as_index=False)
        .agg(
            avg_fidelity_mean=("average_gate_fidelity", "mean"),
            avg_fidelity_best=("average_gate_fidelity", "max"),
            residual_z_mean=("mean_residual_z_error_rad", "mean"),
            residual_z_best=("mean_residual_z_error_rad", "min"),
            transverse_mean=("mean_transverse_error_rad", "mean"),
            transverse_best=("mean_transverse_error_rad", "min"),
            leakage_mean=("leakage_outside_target_mean", "mean"),
        )
        .sort_values("avg_fidelity_mean", ascending=False)
    )
    family_target_group = (
        df.groupby(["waveform_family", "target_family"], as_index=False)
        .agg(
            avg_fidelity_mean=("average_gate_fidelity", "mean"),
            residual_z_mean=("mean_residual_z_error_rad", "mean"),
            transverse_mean=("mean_transverse_error_rad", "mean"),
            count=("case_id", "count"),
        )
        .sort_values(["target_family", "avg_fidelity_mean"], ascending=[True, False])
    )
    best_overall = df.sort_values("average_gate_fidelity", ascending=False).iloc[0].to_dict()
    lowest_residual_z = df.sort_values("mean_residual_z_error_rad", ascending=True).iloc[0].to_dict()
    lowest_transverse = df.sort_values("mean_transverse_error_rad", ascending=True).iloc[0].to_dict()
    per_case_best = (
        df.sort_values("average_gate_fidelity", ascending=False)
        .groupby("case_id", as_index=False)
        .first()
        .to_dict(orient="records")
    )
    return {
        "study": STUDY_DIR.name,
        "n_rows": int(len(rows)),
        "best_overall": best_overall,
        "lowest_residual_z": lowest_residual_z,
        "lowest_transverse": lowest_transverse,
        "family_summary": family_group.to_dict(orient="records"),
        "family_target_summary": family_target_group.to_dict(orient="records"),
        "best_by_case": per_case_best,
        "rows": rows,
    }


def save_figure(fig, stem: str) -> None:
    fig.savefig(FIGURES_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_summary_figures(df: pd.DataFrame) -> None:
    apply_plot_style()

    family_order = list(df.groupby("waveform_family")["average_gate_fidelity"].mean().sort_values(ascending=False).index)
    family_metrics = df.groupby("waveform_family").agg(
        avg_fidelity=("average_gate_fidelity", "mean"),
        avg_residual_z=("mean_residual_z_error_rad", "mean"),
        avg_transverse=("mean_transverse_error_rad", "mean"),
    )
    family_metrics = family_metrics.loc[family_order]

    fig, axes = plt.subplots(1, 3, figsize=(15.2, 4.0))
    axes[0].bar(family_metrics.index, family_metrics["avg_fidelity"], color="#4477AA")
    axes[0].set_ylabel("Average gate fidelity")
    axes[0].set_xlabel("Waveform family")
    axes[0].tick_params(axis="x", rotation=25)
    axes[1].bar(family_metrics.index, family_metrics["avg_residual_z"], color="#EE6677")
    axes[1].set_ylabel("Mean residual Z error (rad)")
    axes[1].set_xlabel("Waveform family")
    axes[1].tick_params(axis="x", rotation=25)
    axes[2].bar(family_metrics.index, family_metrics["avg_transverse"], color="#228833")
    axes[2].set_ylabel("Mean transverse error (rad)")
    axes[2].set_xlabel("Waveform family")
    axes[2].tick_params(axis="x", rotation=25)
    save_figure(fig, "family_metric_means")

    case_labels = [
        (
            f"{row.model_variant} | na={row.n_active} | chiT={row.chi_t_over_2pi:g}"
            f" | T{row.target_family}{'' if pd.isna(row.random_seed) else f' | s{int(row.random_seed) % 1000:03d}'}"
        )
        for row in df[["case_id", "model_variant", "n_active", "chi_t_over_2pi", "target_family", "random_seed"]]
        .drop_duplicates()
        .itertuples(index=False)
    ]
    pivot = df.pivot(index="case_id", columns="waveform_family", values="average_gate_fidelity")
    pivot = pivot.reindex(columns=family_order)
    fig_height = max(5.2, 0.28 * len(pivot.index))
    fig, ax = plt.subplots(figsize=(9.0, fig_height))
    image = ax.imshow(pivot.to_numpy(dtype=float), aspect="auto", cmap="viridis", vmin=0.0, vmax=1.0)
    ax.set_xticks(np.arange(len(pivot.columns)), pivot.columns, rotation=25, ha="right")
    ax.set_yticks(np.arange(len(pivot.index)), case_labels)
    ax.set_xlabel("Waveform family")
    ax.set_ylabel("Study case")
    fig.colorbar(image, ax=ax, label="Average gate fidelity")
    save_figure(fig, "case_family_fidelity_heatmap")

    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.5))
    palette = {
        "baseline_multitone": "#4477AA",
        "symmetric_two_segment": "#EE6677",
        "echoed_multitone": "#AA3377",
        "complex_envelope": "#228833",
        "basis_expanded": "#CCBB44",
    }
    for family in family_order:
        subset = df[df["waveform_family"] == family]
        axes[0].scatter(
            subset["mean_residual_z_error_rad"],
            subset["average_gate_fidelity"],
            s=42,
            label=family,
            color=palette.get(family),
            alpha=0.85,
        )
        axes[1].scatter(
            subset["mean_transverse_error_rad"],
            subset["average_gate_fidelity"],
            s=42,
            label=family,
            color=palette.get(family),
            alpha=0.85,
        )
    axes[0].set_xlabel("Mean residual Z error (rad)")
    axes[0].set_ylabel("Average gate fidelity")
    axes[1].set_xlabel("Mean transverse error (rad)")
    axes[1].set_ylabel("Average gate fidelity")
    axes[0].legend(frameon=False)
    save_figure(fig, "fidelity_tradeoff_planes")

    fig, ax = plt.subplots(figsize=(6.2, 4.8))
    for family in family_order:
        subset = df[df["waveform_family"] == family]
        ax.scatter(
            subset["mean_residual_z_error_rad"],
            subset["mean_transverse_error_rad"],
            s=46,
            label=family,
            color=palette.get(family),
            alpha=0.85,
        )
    ax.set_xlabel("Mean residual Z error (rad)")
    ax.set_ylabel("Mean transverse error (rad)")
    ax.legend(frameon=False)
    save_figure(fig, "residual_z_vs_transverse")

    representative = df[
        (df["target_family"] == "D")
        & (df["model_variant"] == "chi_plus_chiprime")
        & (df["n_active"] == df["n_active"].max())
        & (df["chi_t_over_2pi"] == max(STUDY_DURATIONS))
    ].copy()
    if representative.empty:
        representative = df.sort_values("average_gate_fidelity", ascending=False).head(len(family_order)).copy()
    else:
        representative = representative.sort_values(["random_seed", "average_gate_fidelity"], ascending=[True, False])
    rep_case_id = str(representative.iloc[0]["case_id"])
    fig, ax = plt.subplots(figsize=(6.8, 4.8))
    for family in family_order:
        artifact = load_json(ARTIFACTS_DIR / "cases" / f"{rep_case_id}_{family}.json")
        samples = artifact["waveform_samples"]
        signal = np.asarray(samples["baseband_real"], dtype=float) + 1.0j * np.asarray(samples["baseband_imag"], dtype=float)
        signal = signal - np.mean(signal)
        dt_s = float(np.mean(np.diff(np.asarray(samples["time_s"], dtype=float))))
        freqs_hz = np.fft.fftfreq(signal.size, d=dt_s)
        positive = freqs_hz >= 0.0
        freqs_mhz = freqs_hz[positive] / 1.0e6
        spectrum = np.abs(np.fft.fft(signal))[positive]
        if np.max(spectrum) > 0.0:
            spectrum = spectrum / np.max(spectrum)
        ax.plot(freqs_mhz, spectrum, linewidth=1.8, label=family, color=palette.get(family))
    ax.set_xlim(0.0, 125.0)
    ax.set_xlabel("Baseband frequency (MHz)")
    ax.set_ylabel("Normalized spectrum")
    ax.legend(frameon=False)
    save_figure(fig, "representative_waveform_spectra")


def write_markdown_summary(summary: dict[str, Any]) -> None:
    lines = [
        f"# Study Summary: {summary['study']}",
        "",
        f"- Rows: {summary['n_rows']}",
        f"- Best overall: {summary['best_overall']['waveform_family']} on {summary['best_overall']['case_id']} with fidelity {summary['best_overall']['average_gate_fidelity']:.6f}",
        f"- Lowest residual Z: {summary['lowest_residual_z']['waveform_family']} on {summary['lowest_residual_z']['case_id']} with mean residual Z {summary['lowest_residual_z']['mean_residual_z_error_rad']:.6f} rad",
        f"- Lowest transverse error: {summary['lowest_transverse']['waveform_family']} on {summary['lowest_transverse']['case_id']} with mean transverse error {summary['lowest_transverse']['mean_transverse_error_rad']:.6f} rad",
        "",
        "## Family summary",
        "",
    ]
    for row in summary["family_summary"]:
        lines.append(
            "- {waveform_family}: mean fidelity {avg_fidelity_mean:.6f}, best fidelity {avg_fidelity_best:.6f}, mean residual Z {residual_z_mean:.6f} rad, mean transverse {transverse_mean:.6f} rad".format(
                **row,
            )
        )
    MARKDOWN_SUMMARY_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=str, default="", help="Optional single case id to run.")
    args = parser.parse_args()

    requests = case_requests()
    if args.case:
        requests = [request for request in requests if request.case_id == args.case]
        if not requests:
            raise SystemExit(f"No study case matched '{args.case}'.")

    all_rows: list[dict[str, Any]] = []
    for index, request in enumerate(requests, start=1):
        print(f"[{index:02d}/{len(requests):02d}] Running {request.case_id}")
        case_start = time.perf_counter()
        context = build_case_context(request)
        rows = run_case(context)
        all_rows.extend(rows)
        print(f"     wrote {len(rows)} waveform-family rows in {time.perf_counter() - case_start:.1f} s")

    payload = {"study": STUDY_DIR.name, "rows": all_rows}
    save_json(RESULTS_PATH, payload)
    summary = build_summary(all_rows)
    save_json(SUMMARY_PATH, summary)
    write_markdown_summary(summary)
    df = pd.DataFrame(all_rows)
    df.to_csv(CSV_PATH, index=False)
    plot_summary_figures(df)
    print(f"Saved {RESULTS_PATH}")
    print(f"Saved {SUMMARY_PATH}")


if __name__ == "__main__":
    main()