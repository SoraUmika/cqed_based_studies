from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np

import runtime_compat  # noqa: F401

from cqed_sim.calibration.conditioned_multitone import ConditionedOptimizationConfig
from cqed_sim.calibration.targeted_subspace_multitone import (
    TargetedSubspaceObjectiveWeights,
    TargetedSubspaceOptimizationConfig,
    analyze_targeted_subspace_operator,
    build_spanning_state_transfer_set,
    optimize_targeted_subspace_multitone,
)

from common import (
    TWO_PI,
    block_rotation_metrics,
    build_model,
    build_target_operator,
    channel_waveform_samples,
    compile_pulse_sequence,
    conditioned_targets_from_target_spec,
    duration_from_chi_t,
    logical_levels,
    make_gaussian_qubit_rotation_pulse,
    make_run_config,
    restricted_blocks,
    restricted_operator_from_full,
    save_json,
    save_waveform_npz,
    simulate_full_operator_on_logical_inputs,
    state_validation_summary_for_compiled,
    target_spec,
)
from metrics import (
    PROBE_QUBIT_STATES,
    addressed_indices,
    average_gate_fidelity,
    build_cpsqr_joint_target,
    channel_process_fidelity_to_unitary,
    coherent_error_decomposition,
    fit_cpsqr_block,
    fit_cpsqr_channel,
    full_state_fidelity,
    leakage_outside_indices,
    nearest_unitary,
    operator_2norm_error,
    probe_tier_summary,
    process_fidelity,
    qubit_channel_kraus_from_full,
    qubit_probe_fidelity_rows_for_channel,
    qubit_rx,
    same_manifold_block,
    unitary_rotation_parameters,
)


OBJECTIVE_WEIGHTS = TargetedSubspaceObjectiveWeights(
    qubit_weight=0.15,
    subspace_weight=1.0,
    preservation_weight=0.35,
    leakage_weight=0.35,
)


@dataclass(frozen=True)
class CaseRequest:
    stage: str
    model_variant: str
    include_chi_prime: bool
    target_family: str
    n_active: int
    chi_t_over_2pi: float
    random_seed: int | None = None

    @property
    def case_id(self) -> str:
        duration_label = str(self.chi_t_over_2pi).replace(".", "p")
        seed_part = "" if self.random_seed is None else f"_seed{int(self.random_seed)}"
        return f"{self.model_variant}_{self.target_family}_na{int(self.n_active)}_chiT{duration_label}{seed_part}"


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
    half_target_operator: np.ndarray
    half_optimization: Any


def correction_bounds(n_active: int) -> tuple[np.ndarray, np.ndarray]:
    n = int(n_active)
    lower = np.concatenate(
        [
            np.full(n, -0.75, dtype=float),
            np.full(n, -np.pi, dtype=float),
            np.full(n, -2.0 * np.pi * 3.0e6, dtype=float),
        ]
    )
    upper = np.concatenate(
        [
            np.full(n, 0.75, dtype=float),
            np.full(n, np.pi, dtype=float),
            np.full(n, 2.0 * np.pi * 3.0e6, dtype=float),
        ]
    )
    return lower, upper


def strict_optimization_config(n_active: int, *, scale: float = 1.0) -> TargetedSubspaceOptimizationConfig:
    conditioned = ConditionedOptimizationConfig(
        active_levels=tuple(range(int(n_active))),
        parameters=("d_lambda", "d_alpha", "d_omega"),
        maxiter_stage1=max(10, int(round(scale * (8 + 2 * n_active)))),
        maxiter_stage2=max(14, int(round(scale * (10 + 2 * n_active)))),
        d_lambda_bounds=(-0.75, 0.75),
        d_alpha_bounds=(-np.pi, np.pi),
        d_omega_hz_bounds=(-3.0e6, 3.0e6),
        regularization_lambda=5.0e-4,
        regularization_alpha=5.0e-4,
        regularization_omega=5.0e-4,
    )
    return TargetedSubspaceOptimizationConfig(conditioned=conditioned, include_block_phase=False)


def build_case_context(request: CaseRequest) -> CaseContext:
    spec = target_spec(request.target_family, request.n_active, seed=request.random_seed)
    targets = conditioned_targets_from_target_spec(spec)
    model = build_model(include_chi_prime=request.include_chi_prime, n_active=request.n_active)
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
        optimization_config=strict_optimization_config(request.n_active, scale=0.9 if request.stage == "screen" else 1.0),
        objective_weights=OBJECTIVE_WEIGHTS,
        target_operator=target_operator,
        transfer_set=transfer_set,
        label=f"strict_{request.case_id}",
    )
    half_spec = type(spec)(
        family=spec.family,
        theta_values=tuple(float(value / 2.0) for value in spec.theta_values),
        phi_values=spec.phi_values,
        metadata={
            **dict(spec.metadata),
            "derived_from": str(spec.family),
            "description": f"Half-angle target derived from {spec.family}",
            "theta_values_rad": [float(value / 2.0) for value in spec.theta_values],
        },
    )
    half_target_operator = build_target_operator(half_spec, levels)
    half_targets = conditioned_targets_from_target_spec(half_spec)
    half_run_config = make_run_config(model, n_active=request.n_active, duration_s=0.5 * duration_s)
    half = optimize_targeted_subspace_multitone(
        model,
        half_targets,
        half_run_config,
        logical_levels=levels,
        optimization_config=strict_optimization_config(request.n_active, scale=0.9 if request.stage == "screen" else 1.0),
        objective_weights=OBJECTIVE_WEIGHTS,
        target_operator=half_target_operator,
        transfer_set=build_spanning_state_transfer_set(half_target_operator),
        label=f"half_{request.case_id}",
    )
    return CaseContext(
        request=request,
        model=model,
        frame=run_config.frame,
        levels=levels,
        duration_s=duration_s,
        run_config=run_config,
        spec=spec,
        targets=targets,
        target_operator=target_operator,
        transfer_set=transfer_set,
        direct_optimization=direct,
        half_spec=half_spec,
        half_target_operator=half_target_operator,
        half_optimization=half,
    )


def quick_candidate_metrics(context: CaseContext, full_operator: np.ndarray) -> dict[str, Any]:
    restricted_operator = restricted_operator_from_full(full_operator, context.model, context.levels)
    strict_validation = analyze_targeted_subspace_operator(
        full_operator,
        context.model,
        context.targets,
        logical_levels=context.levels,
        objective_weights=OBJECTIVE_WEIGHTS,
        target_operator=context.target_operator,
        transfer_set=context.transfer_set,
        metadata={"stage": context.request.stage},
    )
    cpsqr_target, cpsqr_rows = build_cpsqr_joint_target(restricted_operator, context.spec.theta_values)
    strict_values = []
    cpsqr_values = []
    for level, theta_target in zip(context.levels, context.spec.theta_values, strict=True):
        kraus_ops = qubit_channel_kraus_from_full(full_operator, int(context.model.n_cav), int(level))
        strict_values.append(channel_process_fidelity_to_unitary(kraus_ops, qubit_rx(theta_target)))
        cpsqr_values.append(fit_cpsqr_channel(kraus_ops, theta_target).process_fidelity)
    strict_reduced_mean = float(np.mean(strict_values))
    cpsqr_reduced_mean = float(np.mean(cpsqr_values))
    return {
        "restricted_operator": restricted_operator,
        "strict_validation": strict_validation,
        "cpsqr_target": cpsqr_target,
        "cpsqr_rows": cpsqr_rows,
        "strict_reduced_mean": strict_reduced_mean,
        "cpsqr_reduced_mean": cpsqr_reduced_mean,
        "strict_joint": float(strict_validation.restricted_process_fidelity),
        "cpsqr_joint": float(process_fidelity(cpsqr_target, restricted_operator)),
        "strict_objective": float(strict_validation.weighted_loss + 0.10 * (1.0 - strict_reduced_mean)),
        "cpsqr_objective": float(
            0.15 * (1.0 - cpsqr_reduced_mean)
            + 1.0 * (1.0 - process_fidelity(cpsqr_target, restricted_operator))
            + 0.35 * float(strict_validation.other_target_population_mean)
            + 0.35 * float(strict_validation.leakage_outside_target_mean)
        ),
    }


def evaluate_sequence_fast(
    context: CaseContext,
    pulses: Sequence[Any],
    drive_ops: dict[str, str],
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
    metrics = quick_candidate_metrics(context, full_operator)
    metrics["compiled"] = compiled
    metrics["full_operator"] = full_operator
    return metrics


def full_probe_rows(context: CaseContext, full_operator: np.ndarray, cpsqr_target: np.ndarray) -> list[dict[str, Any]]:
    full = np.asarray(full_operator, dtype=np.complex128)
    n_cav = int(context.model.n_cav)
    keep_indices = addressed_indices(n_cav, context.levels)
    rows: list[dict[str, Any]] = []
    for row_index, level in enumerate(context.levels):
        strict_block = context.target_operator[2 * row_index : 2 * row_index + 2, 2 * row_index : 2 * row_index + 2]
        cpsqr_block_local = cpsqr_target[2 * row_index : 2 * row_index + 2, 2 * row_index : 2 * row_index + 2]
        for probe_name, qubit_state in PROBE_QUBIT_STATES.items():
            psi0 = np.zeros(2 * n_cav, dtype=np.complex128)
            psi0[int(level)] = complex(qubit_state[0])
            psi0[n_cav + int(level)] = complex(qubit_state[1])
            actual = full @ psi0
            strict_target = np.zeros_like(actual)
            strict_target[int(level)] = strict_block[0, 0] * qubit_state[0] + strict_block[0, 1] * qubit_state[1]
            strict_target[n_cav + int(level)] = strict_block[1, 0] * qubit_state[0] + strict_block[1, 1] * qubit_state[1]
            cpsqr_target_state = np.zeros_like(actual)
            cpsqr_target_state[int(level)] = cpsqr_block_local[0, 0] * qubit_state[0] + cpsqr_block_local[0, 1] * qubit_state[1]
            cpsqr_target_state[n_cav + int(level)] = cpsqr_block_local[1, 0] * qubit_state[0] + cpsqr_block_local[1, 1] * qubit_state[1]
            same_manifold = float(abs(actual[int(level)]) ** 2 + abs(actual[n_cav + int(level)]) ** 2)
            rows.append(
                {
                    "level": int(level),
                    "probe_label": str(probe_name),
                    "strict_fidelity": full_state_fidelity(actual, strict_target),
                    "cpsqr_fidelity": full_state_fidelity(actual, cpsqr_target_state),
                    "same_manifold_population": float(same_manifold),
                    "addressed_window_population": float(np.sum(np.abs(actual[keep_indices]) ** 2)),
                    "leakage_outside_addressed": leakage_outside_indices(actual, keep_indices),
                }
            )
    return rows


def classify_result(row: dict[str, Any]) -> tuple[int, str]:
    single = float(row["strict_reduced_single_ground_mean"])
    strict_reduced = float(row["strict_reduced_quartet_mean"])
    strict_full = float(row["strict_full_quartet_mean"])
    strict_joint = float(row["strict_joint_process_fidelity"])
    cpsqr_reduced = float(row["cpsqr_reduced_quartet_mean"])
    cpsqr_full = float(row["cpsqr_full_quartet_mean"])
    cpsqr_joint = float(row["cpsqr_joint_process_fidelity"])
    if strict_reduced >= 0.99 and strict_full >= 0.99 and strict_joint >= 0.99:
        return 5, "full_ideal_sqr_success"
    if strict_reduced >= 0.99:
        return 4, "reduced_ideal_sqr_success"
    if cpsqr_reduced >= 0.99 and cpsqr_full >= 0.99 and cpsqr_joint >= 0.99:
        return 3, "cpsqr_success"
    if max(strict_reduced, cpsqr_reduced) >= 0.99:
        return 2, "reduced_state_success_only"
    if single >= 0.99 and max(strict_reduced, cpsqr_reduced) < 0.95:
        return 1, "single_input_success_only"
    return 0, "no_useful_control"


def evaluate_candidate_full(
    context: CaseContext,
    family_name: str,
    *,
    pulses: Sequence[Any],
    drive_ops: dict[str, str],
    metadata: dict[str, Any],
    optimizer_payload: dict[str, Any] | None,
    objective_mode: str,
) -> tuple[dict[str, Any], dict[str, Any]]:
    fast = evaluate_sequence_fast(context, pulses, drive_ops)
    compiled = fast["compiled"]
    full_operator = fast["full_operator"]
    restricted_operator = fast["restricted_operator"]
    strict_validation = fast["strict_validation"]
    cpsqr_target = fast["cpsqr_target"]
    cpsqr_rows = fast["cpsqr_rows"]
    reduced_level_rows: list[dict[str, Any]] = []
    reduced_probe_rows: list[dict[str, Any]] = []
    for level, theta_target in zip(context.levels, context.spec.theta_values, strict=True):
        kraus_ops = qubit_channel_kraus_from_full(full_operator, int(context.model.n_cav), int(level))
        same_block = same_manifold_block(full_operator, int(context.model.n_cav), int(level))
        dominant_unitary = nearest_unitary(same_block)
        achieved_theta, achieved_phi, achieved_axis_z = unitary_rotation_parameters(dominant_unitary)
        strict_target_block = qubit_rx(theta_target)
        cpsqr_fit_channel = fit_cpsqr_channel(kraus_ops, theta_target)
        cpsqr_fit_block = fit_cpsqr_block(dominant_unitary, theta_target)
        probe_rows = qubit_probe_fidelity_rows_for_channel(
            kraus_ops,
            strict_unitary=strict_target_block,
            cpsqr_unitary=cpsqr_fit_channel.target_block,
        )
        for row in probe_rows:
            row["level"] = int(level)
        reduced_probe_rows.extend(probe_rows)
        reduced_level_rows.append(
            {
                "level": int(level),
                "target_theta_rad": float(theta_target),
                "achieved_theta_rad": float(achieved_theta),
                "achieved_phi_rad": float(achieved_phi),
                "achieved_axis_z": float(achieved_axis_z),
                "strict_process_fidelity": float(process_fidelity(strict_target_block, dominant_unitary)),
                "cpsqr_process_fidelity": float(cpsqr_fit_channel.process_fidelity),
                "cpsqr_delta_rad": float(cpsqr_fit_channel.delta_rad),
                **{f"strict_{key}": float(value) for key, value in coherent_error_decomposition(dominant_unitary, strict_target_block).items()},
                **{f"cpsqr_{key}": float(value) for key, value in coherent_error_decomposition(dominant_unitary, cpsqr_fit_block.target_block).items()},
            }
        )
    reduced_strict_tiers = probe_tier_summary(reduced_probe_rows, "strict_fidelity")
    reduced_cpsqr_tiers = probe_tier_summary(reduced_probe_rows, "cpsqr_fidelity")
    full_probe = full_probe_rows(context, full_operator, cpsqr_target)
    full_strict_tiers = probe_tier_summary(full_probe, "strict_fidelity")
    full_cpsqr_tiers = probe_tier_summary(full_probe, "cpsqr_fidelity")
    strict_per_block = [
        {"level": int(level), **block_rotation_metrics(target_block, actual_block)}
        for level, target_block, actual_block in zip(
            context.levels,
            restricted_blocks(context.target_operator),
            restricted_blocks(restricted_operator),
            strict=True,
        )
    ]
    cpsqr_per_block = [
        {"level": int(idx), **block_rotation_metrics(cpsqr_target[2 * idx : 2 * idx + 2, 2 * idx : 2 * idx + 2], restricted_operator[2 * idx : 2 * idx + 2, 2 * idx : 2 * idx + 2])}
        for idx in range(len(context.levels))
    ]
    row = {
        "stage": str(context.request.stage),
        "case_id": context.request.case_id,
        "family_name": str(family_name),
        "objective_mode": str(objective_mode),
        "model_variant": str(context.request.model_variant),
        "include_chi_prime": bool(context.request.include_chi_prime),
        "target_family": str(context.request.target_family),
        "random_seed": None if context.request.random_seed is None else int(context.request.random_seed),
        "n_active": int(context.request.n_active),
        "chi_t_over_2pi": float(context.request.chi_t_over_2pi),
        "duration_ns": float(context.duration_s * 1.0e9),
        "active_duration_ns": float(metadata.get("active_duration_s", context.duration_s) * 1.0e9),
        "total_gate_duration_ns": float(metadata.get("total_gate_duration_s", context.duration_s) * 1.0e9),
        "fairness_mode": str(metadata.get("fairness_mode", "fixed_total_duration")),
        "strict_reduced_process_mean": float(np.mean([item["strict_process_fidelity"] for item in reduced_level_rows])),
        "cpsqr_reduced_process_mean": float(np.mean([item["cpsqr_process_fidelity"] for item in reduced_level_rows])),
        "strict_reduced_single_ground_mean": float(reduced_strict_tiers["single_ground"]["mean_fidelity"]),
        "strict_reduced_pair_mean": float(reduced_strict_tiers["selected_pair"]["mean_fidelity"]),
        "strict_reduced_quartet_mean": float(reduced_strict_tiers["spanning_quartet"]["mean_fidelity"]),
        "strict_reduced_quartet_min": float(reduced_strict_tiers["spanning_quartet"]["min_fidelity"]),
        "cpsqr_reduced_single_ground_mean": float(reduced_cpsqr_tiers["single_ground"]["mean_fidelity"]),
        "cpsqr_reduced_pair_mean": float(reduced_cpsqr_tiers["selected_pair"]["mean_fidelity"]),
        "cpsqr_reduced_quartet_mean": float(reduced_cpsqr_tiers["spanning_quartet"]["mean_fidelity"]),
        "cpsqr_reduced_quartet_min": float(reduced_cpsqr_tiers["spanning_quartet"]["min_fidelity"]),
        "strict_full_quartet_mean": float(full_strict_tiers["spanning_quartet"]["mean_fidelity"]),
        "strict_full_quartet_min": float(full_strict_tiers["spanning_quartet"]["min_fidelity"]),
        "cpsqr_full_quartet_mean": float(full_cpsqr_tiers["spanning_quartet"]["mean_fidelity"]),
        "cpsqr_full_quartet_min": float(full_cpsqr_tiers["spanning_quartet"]["min_fidelity"]),
        "strict_joint_process_fidelity": float(strict_validation.restricted_process_fidelity),
        "strict_joint_average_gate_fidelity": float(average_gate_fidelity(context.target_operator, restricted_operator)),
        "strict_joint_operator_2norm": float(operator_2norm_error(context.target_operator, restricted_operator)),
        "cpsqr_joint_process_fidelity": float(process_fidelity(cpsqr_target, restricted_operator)),
        "cpsqr_joint_average_gate_fidelity": float(average_gate_fidelity(cpsqr_target, restricted_operator)),
        "cpsqr_joint_operator_2norm": float(operator_2norm_error(cpsqr_target, restricted_operator)),
        "same_block_population_mean": float(strict_validation.same_block_population_mean),
        "same_block_population_min": float(strict_validation.same_block_population_min),
        "other_target_population_mean": float(strict_validation.other_target_population_mean),
        "leakage_outside_target_mean": float(strict_validation.leakage_outside_target_mean),
        "strict_mean_residual_z_error_rad": float(np.mean([item["residual_z_error_rad"] for item in strict_per_block])),
        "strict_mean_transverse_error_rad": float(np.mean([item["transverse_error_rad"] for item in strict_per_block])),
        "cpsqr_mean_residual_z_error_rad": float(np.mean([item["residual_z_error_rad"] for item in cpsqr_per_block])),
        "cpsqr_mean_transverse_error_rad": float(np.mean([item["transverse_error_rad"] for item in cpsqr_per_block])),
    }
    classification_level, classification_label = classify_result(row)
    row["classification_level"] = int(classification_level)
    row["classification_label"] = str(classification_label)
    artifact = {
        "summary_row": row,
        "metadata": metadata,
        "optimizer": optimizer_payload,
        "target_spec": {
            "family": str(context.spec.family),
            "theta_values_rad": [float(x) for x in context.spec.theta_values],
            "phi_values_rad": [float(x) for x in context.spec.phi_values],
            "metadata": dict(context.spec.metadata),
        },
        "strict_validation": strict_validation.as_dict(),
        "cpsqr_fit_rows": cpsqr_rows,
        "reduced_level_rows": reduced_level_rows,
        "reduced_probe_rows": reduced_probe_rows,
        "full_probe_rows": full_probe,
        "strict_per_block": strict_per_block,
        "cpsqr_per_block": cpsqr_per_block,
        "restricted_operator": restricted_operator,
        "strict_target_operator": context.target_operator,
        "cpsqr_target_operator": cpsqr_target,
        "full_operator_columns_on_logical_inputs": full_operator,
        "waveform_samples": channel_waveform_samples(compiled, channel="qubit"),
        "state_validation_summary": state_validation_summary_for_compiled(
            context.model,
            compiled,
            frame=context.frame,
            drive_ops=drive_ops,
            levels=context.levels,
            target_operator=context.target_operator,
        ),
    }
    return row, artifact


def save_case_artifact(case_id: str, family_name: str, artifact: dict[str, Any], waveform_samples: dict[str, Any] | None = None) -> None:
    from common import ARTIFACTS_DIR  # local import to keep this module light

    json_file = ARTIFACTS_DIR / "cases" / f"{case_id}_{family_name}.json"
    save_json(json_file, artifact)
    if waveform_samples is not None:
        npz_file = ARTIFACTS_DIR / "waveforms" / f"{case_id}_{family_name}.npz"
        save_waveform_npz(npz_file, waveform_samples)
