"""Run the rigorous echoed-ansatz follow-up for strict no-detuning multitone SQR."""

from __future__ import annotations

import argparse
import time
from dataclasses import dataclass
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import Bounds, minimize

from cqed_sim.calibration.targeted_subspace_multitone import analyze_targeted_subspace_operator, build_spanning_state_transfer_set

from common import (
    ARTIFACTS_DIR,
    DATA_DIR,
    DEFAULT_DT,
    FIGURES_DIR,
    IDEAL_X_PI,
    OBJECTIVE_WEIGHTS,
    PI_PULSE_DURATION_S,
    STUDY_DIR,
    CaseRequest,
    OptimizationResult,
    analyze_full_operator,
    apply_plot_style,
    average_gate_fidelity,
    block_rotation_metrics,
    build_frame,
    build_model,
    build_square_multitone_compiled,
    build_target_operator,
    channel_waveform_samples,
    compile_pulse_sequence,
    conditioned_targets_from_target_spec,
    corrections_from_vector,
    corrections_to_dict,
    corrections_to_vector,
    decoupled_block_operator,
    duration_from_chi_t,
    embed_restricted_operator_in_full,
    embed_qubit_operator,
    json_ready,
    logical_levels,
    magnus_effective_blocks,
    make_gaussian_qubit_rotation_pulse,
    make_run_config,
    optimize_square_multitone,
    probe_state_metrics_from_full_operator,
    process_fidelity,
    reduced_blockwise_operator,
    restricted_blocks,
    save_json,
    save_waveform_npz,
    scaled_target_spec,
    shift_pulse,
    simulate_full_operator_on_logical_inputs,
    target_spec,
)


RESULTS_PATH = DATA_DIR / "study_results.json"
SUMMARY_PATH = DATA_DIR / "study_summary.json"
CSV_PATH = DATA_DIR / "study_results.csv"
MARKDOWN_PATH = DATA_DIR / "study_summary.md"
VALIDATION_PATH = DATA_DIR / "validation_summary.json"
METRICS_PATH = DATA_DIR / "metric_definitions.json"

CASE_DIR = ARTIFACTS_DIR / "cases"
WAVEFORM_DIR = ARTIFACTS_DIR / "waveforms"

MODEL_VARIANTS = (
    ("chi_only", False),
    ("chi_plus_chiprime", True),
)
BASE_FAMILIES = ("aligned_x", "structured_xy")
ACTIVE_GRID = (2, 3)
DURATION_GRID = (3.0, 5.0)


@dataclass(frozen=True)
class EchoSolveResult:
    corrections_1: Any
    corrections_2: Any
    eta: float
    analysis: Any
    full_operator: np.ndarray
    metadata: dict[str, Any]
    optimizer_payload: dict[str, Any]
    segment_artifacts: dict[str, Any]


@dataclass(frozen=True)
class RefocusPulseResult:
    spec: Any
    optimization: OptimizationResult
    validation: Any
    target_operator: np.ndarray


@dataclass(frozen=True)
class CaseContext:
    request: CaseRequest
    model: Any
    frame: Any
    levels: tuple[int, ...]
    spec: Any
    targets: Any
    target_operator: np.ndarray
    transfer_set: Any
    duration_s: float
    run_config: Any
    matched_total_duration_s: float
    matched_run_config: Any


def stable_seed(case_id: str, offset: int = 0) -> int:
    total = 0
    for index, char in enumerate(case_id):
        total += (index + 1) * ord(char)
    return int(total + 97 * offset)


def case_requests() -> list[CaseRequest]:
    rows: list[CaseRequest] = []
    for model_variant, include_chi_prime in MODEL_VARIANTS:
        for family in BASE_FAMILIES:
            for n_active in ACTIVE_GRID:
                for chi_t in DURATION_GRID:
                    rows.append(
                        CaseRequest(
                            model_variant=model_variant,
                            include_chi_prime=include_chi_prime,
                            family=family,
                            n_active=n_active,
                            chi_t_over_2pi=chi_t,
                        )
                    )
    return rows


def build_context(request: CaseRequest) -> CaseContext:
    spec = target_spec(request.family, request.n_active, seed=request.seed)
    targets = conditioned_targets_from_target_spec(spec)
    model = build_model(include_chi_prime=request.include_chi_prime, n_active=request.n_active)
    frame = build_frame(model)
    levels = logical_levels(request.n_active)
    duration_s = duration_from_chi_t(request.chi_t_over_2pi)
    run_config = make_run_config(model, n_active=request.n_active, duration_s=duration_s)
    matched_total_duration_s = duration_s + 2.0 * PI_PULSE_DURATION_S
    matched_run_config = make_run_config(model, n_active=request.n_active, duration_s=matched_total_duration_s)
    target_operator = build_target_operator(spec, levels)
    transfer_set = build_spanning_state_transfer_set(target_operator)
    return CaseContext(
        request=request,
        model=model,
        frame=frame,
        levels=levels,
        spec=spec,
        targets=targets,
        target_operator=target_operator,
        transfer_set=transfer_set,
        duration_s=duration_s,
        run_config=run_config,
        matched_total_duration_s=matched_total_duration_s,
        matched_run_config=matched_run_config,
    )


def use_multitone_pi_benchmark(request: CaseRequest) -> bool:
    return bool(request.include_chi_prime and request.chi_t_over_2pi == max(DURATION_GRID))


def metric_definitions() -> dict[str, Any]:
    return {
        "restricted_process_fidelity": "Process fidelity between the ideal restricted operator and the actual restricted operator, including block-dependent phase errors.",
        "restricted_average_gate_fidelity": "Average gate fidelity computed from the restricted process fidelity.",
        "best_fit_restricted_process_fidelity": "Restricted process fidelity after the framework's best-fit logical block-phase correction.",
        "state_transfer_fidelity_mean": "Mean fidelity over the framework spanning transfer set.",
        "same_block_population_mean": "Mean retained population inside the correct logical block over basis probes.",
        "other_target_population_mean": "Mean population transferred into the wrong addressed blocks.",
        "leakage_outside_target_mean": "Mean population leaked outside the addressed logical subspace.",
        "mean_residual_z_error_rad": "Mean absolute blockwise residual-Z component extracted from the SU(2) error generator relative to the target block rotation.",
        "max_residual_z_error_rad": "Maximum absolute blockwise residual-Z component across addressed manifolds.",
        "mean_transverse_error_rad": "Mean norm of the blockwise transverse error-generator component.",
        "probe_fidelity_mean": "Mean state fidelity over a fixed explicit logical probe set that includes cavity superpositions and fixed random logical states.",
        "probe_fidelity_min": "Worst probe-state fidelity over that explicit logical probe set.",
        "refocus_target": "For the manifold-aware refocusing pulse, the target is blockwise Rx(pi) on every addressed manifold rather than the SQR target.",
    }


def analysis_against_target(context: CaseContext, full_operator: np.ndarray, *, metadata: dict[str, Any]) -> Any:
    return analyze_targeted_subspace_operator(
        np.asarray(full_operator, dtype=np.complex128),
        context.model,
        context.targets,
        logical_levels=context.levels,
        objective_weights=OBJECTIVE_WEIGHTS,
        target_operator=context.target_operator,
        transfer_set=context.transfer_set,
        metadata=dict(metadata),
    )


def block_rows_from_analysis(target_operator: np.ndarray, actual_operator: np.ndarray) -> list[dict[str, float]]:
    return [
        block_rotation_metrics(target_block, actual_block)
        for target_block, actual_block in zip(restricted_blocks(target_operator), restricted_blocks(actual_operator), strict=True)
    ]


def row_from_analysis(
    context: CaseContext,
    *,
    construction: str,
    comparison_spec: Any,
    comparison_target_operator: np.ndarray,
    analysis: Any,
    full_operator: np.ndarray,
    extra: dict[str, Any] | None = None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    restricted = np.asarray(analysis.restricted_operator, dtype=np.complex128)
    block_rows = block_rows_from_analysis(comparison_target_operator, restricted)
    probe_metrics = probe_state_metrics_from_full_operator(
        context.model,
        context.levels,
        full_operator=np.asarray(full_operator, dtype=np.complex128),
        target_operator=np.asarray(comparison_target_operator, dtype=np.complex128),
    )
    row = {
        "case_id": str(context.request.case_id),
        "construction": str(construction),
        "model_variant": str(context.request.model_variant),
        "include_chi_prime": bool(context.request.include_chi_prime),
        "family": str(context.request.family),
        "seed": None if context.request.seed is None else int(context.request.seed),
        "n_active": int(context.request.n_active),
        "chi_t_over_2pi": float(context.request.chi_t_over_2pi),
        "active_duration_ns": float(context.duration_s * 1.0e9),
        "evaluation_target_family": str(comparison_spec.family),
        "target_theta_values_rad": [float(x) for x in comparison_spec.theta_values],
        "target_phi_values_rad": [float(x) for x in comparison_spec.phi_values],
        "restricted_process_fidelity": float(analysis.restricted_process_fidelity),
        "restricted_average_gate_fidelity": float(average_gate_fidelity(comparison_target_operator, restricted)),
        "best_fit_restricted_process_fidelity": float(analysis.best_fit_restricted_process_fidelity),
        "uncorrected_restricted_process_fidelity": float(analysis.uncorrected_restricted_process_fidelity),
        "restricted_fro_error": float(analysis.restricted_fro_error),
        "restricted_unitarity_error": float(analysis.restricted_unitarity_error),
        "state_transfer_fidelity_mean": float(analysis.state_transfer_fidelity_mean),
        "state_transfer_fidelity_min": float(analysis.state_transfer_fidelity_min),
        "same_block_population_mean": float(analysis.same_block_population_mean),
        "same_block_population_min": float(analysis.same_block_population_min),
        "other_target_population_mean": float(analysis.other_target_population_mean),
        "other_target_population_max": float(analysis.other_target_population_max),
        "leakage_outside_target_mean": float(analysis.leakage_outside_target_mean),
        "leakage_outside_target_max": float(analysis.leakage_outside_target_max),
        "weighted_loss": float(analysis.weighted_loss),
        "best_fit_block_phase_rms_rad": float(
            float("nan") if analysis.block_phase_diagnostics is None else analysis.block_phase_diagnostics.rms_block_phase_error_rad
        ),
        "probe_fidelity_mean": float(probe_metrics["probe_fidelity_mean"]),
        "probe_fidelity_min": float(probe_metrics["probe_fidelity_min"]),
        "per_block_process_fidelities": [float(item["process_fidelity"]) for item in block_rows],
        "per_block_average_gate_fidelities": [float(item["average_gate_fidelity"]) for item in block_rows],
        "per_block_rotation_angle_errors_rad": [float(item["rotation_angle_error_rad"]) for item in block_rows],
        "per_block_rotation_axis_errors_rad": [float(item["rotation_axis_error_rad"]) for item in block_rows],
        "per_block_residual_z_error_rad": [float(item["residual_z_error_rad"]) for item in block_rows],
        "per_block_transverse_error_rad": [float(item["transverse_error_rad"]) for item in block_rows],
        "mean_block_average_gate_fidelity": float(np.mean([item["average_gate_fidelity"] for item in block_rows])),
        "worst_block_average_gate_fidelity": float(np.min([item["average_gate_fidelity"] for item in block_rows])),
        "mean_residual_z_error_rad": float(np.mean([item["residual_z_error_rad"] for item in block_rows])),
        "max_residual_z_error_rad": float(np.max([item["residual_z_error_rad"] for item in block_rows])),
        "mean_transverse_error_rad": float(np.mean([item["transverse_error_rad"] for item in block_rows])),
        "max_transverse_error_rad": float(np.max([item["transverse_error_rad"] for item in block_rows])),
    }
    if extra:
        row.update(extra)
    artifact = {
        "study_name": STUDY_DIR.name,
        "case_request": json_ready(context.request.__dict__),
        "construction": str(construction),
        "comparison_spec": json_ready(comparison_spec.metadata),
        "restricted_operator": restricted,
        "full_operator_columns_on_logical_inputs": np.asarray(full_operator, dtype=np.complex128),
        "validation": analysis.as_dict(),
        "probe_metrics": probe_metrics,
        "summary_row": row,
    }
    return row, artifact


def save_case_artifact(case_id: str, construction: str, artifact: dict[str, Any]) -> None:
    stem = f"{case_id}_{construction}"
    waveform_samples = artifact.get("waveform_samples")
    if waveform_samples is not None:
        npz_path = WAVEFORM_DIR / f"{stem}.npz"
        save_waveform_npz(npz_path, waveform_samples)
        artifact["waveform_npz"] = str(npz_path.relative_to(STUDY_DIR))
    save_json(CASE_DIR / f"{stem}.json", artifact, description=f"Artifact for {construction} on {case_id}.")


def echo_bounds(n_active: int) -> Bounds:
    lower = np.concatenate(
        [
            np.full(n_active, -0.75),
            np.full(n_active, -np.pi),
            np.full(n_active, -0.75),
            np.full(n_active, -np.pi),
            np.asarray([-0.30], dtype=float),
        ]
    )
    upper = np.concatenate(
        [
            np.full(n_active, 0.75),
            np.full(n_active, np.pi),
            np.full(n_active, 0.75),
            np.full(n_active, np.pi),
            np.asarray([0.30], dtype=float),
        ]
    )
    return Bounds(lower, upper)


def parse_echo_vector(vector: np.ndarray, n_active: int) -> tuple[Any, Any, float]:
    arr = np.asarray(vector, dtype=float).reshape(-1)
    first = corrections_from_vector(arr[: 2 * n_active], n_active=n_active)
    second = corrections_from_vector(arr[2 * n_active : 4 * n_active], n_active=n_active)
    eta = float(arr[-1])
    return first, second, eta


def echo_segment_specs(spec: Any, eta: float) -> tuple[Any, Any]:
    theta_1 = 0.5 * (1.0 + float(eta))
    theta_2 = 0.5 * (1.0 - float(eta))
    spec_1 = scaled_target_spec(
        spec,
        theta_scale=theta_1,
        phi_sign=1.0,
        family_suffix="_echo_seg1",
        description="First echoed active segment with original azimuths.",
    )
    spec_2 = scaled_target_spec(
        spec,
        theta_scale=theta_2,
        phi_sign=-1.0,
        family_suffix="_echo_seg2",
        description="Second echoed active segment with toggling-consistent conjugated azimuths.",
    )
    return spec_1, spec_2


def run_config_for_duration(context: CaseContext, duration_s: float):
    return make_run_config(
        context.model,
        n_active=context.request.n_active,
        duration_s=float(duration_s),
        dt_s=float(context.run_config.dt_s),
    )


def square_payload(context: CaseContext, spec: Any, duration_s: float, corrections: Any, *, label: str) -> dict[str, Any]:
    run_config = run_config_for_duration(context, duration_s)
    waveform, tone_specs, compiled = build_square_multitone_compiled(
        context.model,
        spec,
        run_config,
        corrections=corrections,
        label=label,
    )
    full_operator = simulate_full_operator_on_logical_inputs(
        context.model,
        compiled,
        frame=context.frame,
        drive_ops=waveform.drive_ops,
        levels=context.levels,
    )
    return {
        "spec": spec,
        "run_config": run_config,
        "waveform": waveform,
        "tone_specs": tone_specs,
        "compiled": compiled,
        "full_operator": np.asarray(full_operator, dtype=np.complex128),
    }


def build_echo_full_operator(
    context: CaseContext,
    *,
    corrections_1: Any,
    corrections_2: Any,
    eta: float,
    refocus_mode: str,
    refocus_waveform: Any | None = None,
    label: str,
) -> dict[str, Any]:
    duration_1 = 0.5 * context.duration_s * (1.0 + float(eta))
    duration_2 = 0.5 * context.duration_s * (1.0 - float(eta))
    spec_1, spec_2 = echo_segment_specs(context.spec, eta)
    first = square_payload(context, spec_1, duration_1, corrections_1, label=f"{label}_seg1")
    second = square_payload(context, spec_2, duration_2, corrections_2, label=f"{label}_seg2")
    metadata = {
        "eta": float(eta),
        "segment_1_duration_s": float(duration_1),
        "segment_2_duration_s": float(duration_2),
        "active_sqr_duration_s": float(duration_1 + duration_2),
        "segment_1_spec": json_ready(spec_1.metadata),
        "segment_2_spec": json_ready(spec_2.metadata),
        "corrections_segment_1": corrections_to_dict(corrections_1),
        "corrections_segment_2": corrections_to_dict(corrections_2),
        "refocus_mode": str(refocus_mode),
    }
    if refocus_mode == "ideal_instantaneous":
        x_full = embed_qubit_operator(IDEAL_X_PI, n_cav=int(context.model.n_cav))
        full_operator = x_full @ second["full_operator"] @ x_full @ first["full_operator"]
        metadata["total_gate_duration_s"] = float(context.duration_s)
        return {
            "full_operator": np.asarray(full_operator, dtype=np.complex128),
            "metadata": metadata,
            "segment_artifacts": {
                "segment_1_waveform_samples": channel_waveform_samples(first["compiled"]),
                "segment_2_waveform_samples": channel_waveform_samples(second["compiled"]),
            },
        }

    if refocus_mode == "gaussian_pi":
        pi_1 = make_gaussian_qubit_rotation_pulse(
            context.model,
            context.frame,
            theta=np.pi,
            phase=0.0,
            duration_s=PI_PULSE_DURATION_S,
            channel=str(first["waveform"].pulse.channel),
            manifold_level=0,
            t0=float(duration_1),
            label=f"{label}_xpi_1",
        )
        pi_2 = make_gaussian_qubit_rotation_pulse(
            context.model,
            context.frame,
            theta=np.pi,
            phase=0.0,
            duration_s=PI_PULSE_DURATION_S,
            channel=str(first["waveform"].pulse.channel),
            manifold_level=0,
            t0=float(duration_1 + PI_PULSE_DURATION_S + duration_2),
            label=f"{label}_xpi_2",
        )
        pulses = [
            shift_pulse(first["waveform"].pulse, t0=0.0, label=f"{label}_seg1"),
            pi_1,
            shift_pulse(second["waveform"].pulse, t0=float(duration_1 + PI_PULSE_DURATION_S), label=f"{label}_seg2"),
            pi_2,
        ]
        total_duration = float(duration_1 + duration_2 + 2.0 * PI_PULSE_DURATION_S)
        compiled = compile_pulse_sequence(pulses, dt_s=float(context.run_config.dt_s), total_duration_s=total_duration)
        full_operator = simulate_full_operator_on_logical_inputs(
            context.model,
            compiled,
            frame=context.frame,
            drive_ops=first["waveform"].drive_ops,
            levels=context.levels,
        )
        metadata["total_gate_duration_s"] = total_duration
        return {
            "full_operator": np.asarray(full_operator, dtype=np.complex128),
            "compiled": compiled,
            "metadata": metadata,
            "segment_artifacts": {
                "sequence_waveform_samples": channel_waveform_samples(compiled),
                "segment_1_waveform_samples": channel_waveform_samples(first["compiled"]),
                "segment_2_waveform_samples": channel_waveform_samples(second["compiled"]),
            },
        }

    if refocus_mode == "multitone_pi":
        if refocus_waveform is None:
            raise ValueError("multitone_pi refocus requested without a refocus waveform.")
        total_duration = float(duration_1 + duration_2 + 2.0 * PI_PULSE_DURATION_S)
        pulses = [
            shift_pulse(first["waveform"].pulse, t0=0.0, label=f"{label}_seg1"),
            shift_pulse(refocus_waveform.pulse, t0=float(duration_1), label=f"{label}_refocus_1"),
            shift_pulse(second["waveform"].pulse, t0=float(duration_1 + PI_PULSE_DURATION_S), label=f"{label}_seg2"),
            shift_pulse(refocus_waveform.pulse, t0=float(duration_1 + PI_PULSE_DURATION_S + duration_2), label=f"{label}_refocus_2"),
        ]
        compiled = compile_pulse_sequence(pulses, dt_s=float(context.run_config.dt_s), total_duration_s=total_duration)
        full_operator = simulate_full_operator_on_logical_inputs(
            context.model,
            compiled,
            frame=context.frame,
            drive_ops=first["waveform"].drive_ops,
            levels=context.levels,
        )
        metadata["total_gate_duration_s"] = total_duration
        return {
            "full_operator": np.asarray(full_operator, dtype=np.complex128),
            "compiled": compiled,
            "metadata": metadata,
            "segment_artifacts": {
                "sequence_waveform_samples": channel_waveform_samples(compiled),
                "segment_1_waveform_samples": channel_waveform_samples(first["compiled"]),
                "segment_2_waveform_samples": channel_waveform_samples(second["compiled"]),
            },
        }
    raise ValueError(f"Unsupported refocus_mode '{refocus_mode}'.")


def optimize_echo(
    context: CaseContext,
    *,
    construction: str,
    refocus_mode: str,
    x0: np.ndarray,
    maxiter: int,
    maxfev: int,
    refocus_waveform: Any | None = None,
) -> EchoSolveResult:
    n_active = context.request.n_active
    history: list[dict[str, float]] = []

    def objective(vector: np.ndarray) -> float:
        corr_1, corr_2, eta = parse_echo_vector(vector, n_active)
        payload = build_echo_full_operator(
            context,
            corrections_1=corr_1,
            corrections_2=corr_2,
            eta=eta,
            refocus_mode=refocus_mode,
            refocus_waveform=refocus_waveform,
            label=f"{construction}_eval",
        )
        analysis = analysis_against_target(context, payload["full_operator"], metadata={"construction": construction, **payload["metadata"]})
        value = float(analysis.weighted_loss) + 1.0e-4 * float(np.mean(np.asarray(vector, dtype=float) ** 2))
        history.append(
            {
                "objective": float(value),
                "restricted_process_fidelity": float(analysis.restricted_process_fidelity),
                "best_fit_restricted_process_fidelity": float(analysis.best_fit_restricted_process_fidelity),
                "leakage_outside_target_mean": float(analysis.leakage_outside_target_mean),
                "eta": float(eta),
            }
        )
        return value

    start_time = time.perf_counter()
    result = minimize(
        objective,
        np.asarray(x0, dtype=float),
        method="Powell",
        bounds=echo_bounds(n_active),
        options={"maxiter": int(maxiter), "maxfev": int(maxfev), "xtol": 5.0e-3, "ftol": 5.0e-3},
    )
    runtime_s = time.perf_counter() - start_time
    corr_1, corr_2, eta = parse_echo_vector(np.asarray(result.x, dtype=float), n_active)
    payload = build_echo_full_operator(
        context,
        corrections_1=corr_1,
        corrections_2=corr_2,
        eta=eta,
        refocus_mode=refocus_mode,
        refocus_waveform=refocus_waveform,
        label=f"{construction}_best",
    )
    analysis = analysis_against_target(context, payload["full_operator"], metadata={"construction": construction, **payload["metadata"]})
    optimizer_payload = {
        "kind": str(refocus_mode),
        "success": bool(result.success),
        "message": str(result.message),
        "nfev": int(getattr(result, "nfev", -1)),
        "nit": int(getattr(result, "nit", -1)) if getattr(result, "nit", None) is not None else -1,
        "runtime_s": float(runtime_s),
        "objective": float(result.fun),
        "history_tail": history[-80:],
    }
    return EchoSolveResult(
        corrections_1=corr_1,
        corrections_2=corr_2,
        eta=float(eta),
        analysis=analysis,
        full_operator=np.asarray(payload["full_operator"], dtype=np.complex128),
        metadata=dict(payload["metadata"]),
        optimizer_payload=optimizer_payload,
        segment_artifacts=dict(payload["segment_artifacts"]),
    )


def refocus_pi_spec(n_active: int):
    spec_type = type(target_spec("aligned_x", n_active))
    return spec_type(
        family="refocus_x_pi",
        theta_values=tuple(float(np.pi) for _ in range(n_active)),
        phi_values=tuple(0.0 for _ in range(n_active)),
        metadata={
            "family": "refocus_x_pi",
            "description": "Manifold-aware shared-line multitone refocusing target with blockwise Rx(pi).",
            "theta_values_rad": [float(np.pi) for _ in range(n_active)],
            "phi_values_rad": [0.0 for _ in range(n_active)],
            "n_active": int(n_active),
        },
    )


def optimize_refocus_pulse(context: CaseContext, *, maxiter: int, n_starts: int) -> RefocusPulseResult:
    spec = refocus_pi_spec(context.request.n_active)
    run_config = make_run_config(
        context.model,
        n_active=context.request.n_active,
        duration_s=PI_PULSE_DURATION_S,
        dt_s=float(context.run_config.dt_s),
    )
    optimization = optimize_square_multitone(
        context.model,
        spec,
        run_config,
        n_starts=n_starts,
        maxiter=maxiter,
        random_seed=stable_seed(context.request.case_id, 700),
        label_prefix=f"{context.request.case_id}_refocus",
    )
    target_operator = build_target_operator(spec, context.levels)
    validation = analyze_full_operator(
        np.asarray(optimization.validation.full_operator, dtype=np.complex128),
        context.model,
        spec,
        levels=context.levels,
        metadata={"construction": "refocus_multitone_pi"},
    )
    return RefocusPulseResult(
        spec=spec,
        optimization=optimization,
        validation=validation,
        target_operator=target_operator,
    )


def figure_duration_tradeoff(df: pd.DataFrame) -> None:
    apply_plot_style()
    subset = df[df["construction"].isin(["full_shared_line", "full_shared_line_total_matched", "echo_opt_ideal", "echo_opt_gaussian"])].copy()
    if subset.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2), sharey=True)
    palette = {
        "full_shared_line": "#4477AA",
        "full_shared_line_total_matched": "#66CCEE",
        "echo_opt_ideal": "#228833",
        "echo_opt_gaussian": "#EE6677",
    }
    for axis, family in zip(axes, BASE_FAMILIES, strict=True):
        family_df = subset[subset["family"] == family]
        grouped = family_df.groupby(["construction", "chi_t_over_2pi"], as_index=False)["restricted_average_gate_fidelity"].median()
        for construction in ("full_shared_line", "full_shared_line_total_matched", "echo_opt_ideal", "echo_opt_gaussian"):
            rows = grouped[grouped["construction"] == construction]
            if rows.empty:
                continue
            axis.plot(
                rows["chi_t_over_2pi"],
                rows["restricted_average_gate_fidelity"],
                marker="o",
                linewidth=2.0,
                color=palette[construction],
                label=construction.replace("_", " "),
            )
        axis.set_title(family.replace("_", " "))
        axis.set_xlabel(r"$|\chi|T/2\pi$")
        axis.grid(True, alpha=0.25)
        axis.set_ylim(0.0, 1.02)
    axes[0].set_ylabel("Median restricted average gate fidelity")
    axes[1].legend(frameon=False, loc="lower left")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "duration_fidelity_tradeoff.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "duration_fidelity_tradeoff.pdf", bbox_inches="tight")
    plt.close(fig)


def figure_residual_z(df: pd.DataFrame) -> None:
    apply_plot_style()
    subset = df[df["construction"].isin(["full_shared_line", "echo_replay_ideal", "echo_opt_ideal", "echo_opt_gaussian"])].copy()
    if subset.empty:
        return
    grouped = subset.groupby(["construction", "chi_t_over_2pi"], as_index=False)["max_residual_z_error_rad"].mean()
    palette = {
        "full_shared_line": "#4477AA",
        "echo_replay_ideal": "#CCBB44",
        "echo_opt_ideal": "#228833",
        "echo_opt_gaussian": "#EE6677",
    }
    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    for construction in ("full_shared_line", "echo_replay_ideal", "echo_opt_ideal", "echo_opt_gaussian"):
        rows = grouped[grouped["construction"] == construction]
        if rows.empty:
            continue
        ax.plot(
            rows["chi_t_over_2pi"],
            rows["max_residual_z_error_rad"],
            marker="o",
            linewidth=2.0,
            color=palette[construction],
            label=construction.replace("_", " "),
        )
    ax.set_xlabel(r"$|\chi|T/2\pi$")
    ax.set_ylabel("Mean max block residual-Z error [rad]")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "blockwise_residual_z_vs_duration.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "blockwise_residual_z_vs_duration.pdf", bbox_inches="tight")
    plt.close(fig)


def figure_echo_comparison(df: pd.DataFrame) -> None:
    apply_plot_style()
    order = ["full_shared_line", "full_shared_line_total_matched", "echo_replay_ideal", "echo_opt_ideal", "echo_opt_gaussian"]
    labels = ["direct", "direct matched", "echo replay", "echo ideal opt", "echo gauss opt"]
    subset = df[df["construction"].isin(order)].copy()
    if subset.empty:
        return
    grouped = subset.groupby(["family", "construction"], as_index=False)["restricted_average_gate_fidelity"].mean()
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.0), sharey=True)
    palette = ["#4477AA", "#66CCEE", "#CCBB44", "#228833", "#EE6677"]
    for axis, family in zip(axes, BASE_FAMILIES, strict=True):
        family_df = grouped[grouped["family"] == family].set_index("construction")
        values = [float(family_df.loc[key, "restricted_average_gate_fidelity"]) if key in family_df.index else np.nan for key in order]
        axis.bar(labels, values, color=palette)
        axis.set_title(family.replace("_", " "))
        axis.set_ylim(0.0, 1.02)
        axis.grid(True, axis="y", alpha=0.25)
        axis.tick_params(axis="x", rotation=20)
    axes[0].set_ylabel("Mean restricted average gate fidelity")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "plain_vs_echo_comparison.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "plain_vs_echo_comparison.pdf", bbox_inches="tight")
    plt.close(fig)


def figure_echo_tradeoff(df: pd.DataFrame) -> None:
    apply_plot_style()
    direct = df[df["construction"] == "full_shared_line_total_matched"].copy()
    echoes = df[df["construction"].isin(["echo_opt_gaussian", "echo_opt_multitone_pi"])].copy()
    if direct.empty or echoes.empty:
        return
    merged = echoes.merge(
        direct[["case_id", "restricted_average_gate_fidelity", "max_residual_z_error_rad", "probe_fidelity_mean"]],
        on="case_id",
        suffixes=("", "_direct"),
    )
    merged["delta_fidelity"] = merged["restricted_average_gate_fidelity"] - merged["restricted_average_gate_fidelity_direct"]
    merged["delta_residual_z"] = merged["max_residual_z_error_rad"] - merged["max_residual_z_error_rad_direct"]
    palette = {"echo_opt_gaussian": "#EE6677", "echo_opt_multitone_pi": "#228833"}
    fig, ax = plt.subplots(figsize=(6.8, 5.0))
    for construction in ("echo_opt_gaussian", "echo_opt_multitone_pi"):
        rows = merged[merged["construction"] == construction]
        if rows.empty:
            continue
        ax.scatter(
            rows["delta_residual_z"],
            rows["delta_fidelity"],
            s=58,
            alpha=0.85,
            color=palette[construction],
            label=construction.replace("_", " "),
        )
    ax.axhline(0.0, color="black", linewidth=0.8, linestyle="--")
    ax.axvline(0.0, color="black", linewidth=0.8, linestyle="--")
    ax.set_xlabel(r"$\Delta$ max residual-Z error [rad] vs direct matched")
    ax.set_ylabel(r"$\Delta$ restricted average fidelity vs direct matched")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "echo_delta_tradeoff.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "echo_delta_tradeoff.pdf", bbox_inches="tight")
    plt.close(fig)


def figure_refocus_benchmark(df: pd.DataFrame) -> None:
    apply_plot_style()
    subset = df[df["construction"] == "refocus_multitone_pi"].copy()
    if subset.empty:
        return
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    x = np.arange(len(subset))
    ax.bar(x, subset["restricted_average_gate_fidelity"], color="#228833")
    ax.set_xticks(x, [f"{row.family}\nN={int(row.n_active)}" for row in subset.itertuples()], rotation=0)
    ax.set_ylabel(r"Restricted average fidelity to blockwise $R_x(\pi)$")
    ax.set_xlabel("Representative hard cases")
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "refocus_multitone_pi_benchmark.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "refocus_multitone_pi_benchmark.pdf", bbox_inches="tight")
    plt.close(fig)


def build_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    df = pd.DataFrame(rows)
    summary_rows = (
        df.groupby("construction", as_index=False)
        .agg(
            mean_fidelity=("restricted_average_gate_fidelity", "mean"),
            best_fidelity=("restricted_average_gate_fidelity", "max"),
            mean_residual_z=("max_residual_z_error_rad", "mean"),
            mean_probe=("probe_fidelity_mean", "mean"),
            mean_leakage=("leakage_outside_target_mean", "mean"),
            count=("case_id", "count"),
        )
        .sort_values("mean_fidelity", ascending=False)
    )
    best_overall = None if df.empty else df.sort_values("restricted_average_gate_fidelity", ascending=False).iloc[0].to_dict()
    return {
        "study": STUDY_DIR.name,
        "n_rows": int(len(rows)),
        "n_cases": int(df["case_id"].nunique()),
        "metric_definitions": metric_definitions(),
        "construction_summary": summary_rows.to_dict(orient="records"),
        "best_overall": best_overall,
    }


def validation_payload(
    rows: list[dict[str, Any]],
    *,
    representative_direct: dict[str, Any] | None,
    representative_echo: dict[str, Any] | None,
) -> dict[str, Any]:
    frame = pd.DataFrame(rows)
    return {
        "sanity_checks": {
            "decoupled_block_min_fidelity": float(frame.query("construction == 'decoupled_block'")["restricted_average_gate_fidelity"].min()),
            "reduced_replay_min_match": float(frame.query("construction == 'blockwise_exact_reduced'")["restricted_average_gate_fidelity"].min()),
        },
        "convergence_representative_direct": representative_direct,
        "convergence_representative_echo": representative_echo,
        "notes": {
            "dt_default_s": float(DEFAULT_DT),
            "gaussian_pi_duration_s": float(PI_PULSE_DURATION_S),
        },
    }


def representative_convergence(context: CaseContext, *, direct_budget: int, echo_budget: int) -> tuple[dict[str, Any], dict[str, Any]]:
    direct = optimize_square_multitone(
        context.model,
        context.spec,
        context.run_config,
        n_starts=2,
        maxiter=direct_budget,
        random_seed=stable_seed(context.request.case_id, 901),
        label_prefix=f"{context.request.case_id}_conv_direct",
    )
    direct_row = {
        "restricted_average_gate_fidelity": float(average_gate_fidelity(context.target_operator, np.asarray(direct.validation.restricted_operator))),
        "restricted_process_fidelity": float(direct.validation.restricted_process_fidelity),
        "best_fit_restricted_process_fidelity": float(direct.validation.best_fit_restricted_process_fidelity),
    }

    spec_1, spec_2 = echo_segment_specs(context.spec, 0.0)
    seed_1 = optimize_square_multitone(
        context.model,
        spec_1,
        run_config_for_duration(context, 0.5 * context.duration_s),
        n_starts=1,
        maxiter=max(8, echo_budget // 2),
        random_seed=stable_seed(context.request.case_id, 902),
        label_prefix=f"{context.request.case_id}_conv_seg1",
    )
    seed_2 = optimize_square_multitone(
        context.model,
        spec_2,
        run_config_for_duration(context, 0.5 * context.duration_s),
        n_starts=1,
        maxiter=max(8, echo_budget // 2),
        random_seed=stable_seed(context.request.case_id, 903),
        label_prefix=f"{context.request.case_id}_conv_seg2",
    )
    echo_x0 = np.concatenate(
        [
            corrections_to_vector(seed_1.corrections, n_active=context.request.n_active),
            corrections_to_vector(seed_2.corrections, n_active=context.request.n_active),
            np.asarray([0.0], dtype=float),
        ]
    )
    echo = optimize_echo(
        context,
        construction="echo_conv_gaussian",
        refocus_mode="gaussian_pi",
        x0=echo_x0,
        maxiter=echo_budget,
        maxfev=180,
    )
    echo_rows = block_rows_from_analysis(context.target_operator, np.asarray(echo.analysis.restricted_operator))
    echo_row = {
        "restricted_average_gate_fidelity": float(average_gate_fidelity(context.target_operator, np.asarray(echo.analysis.restricted_operator))),
        "restricted_process_fidelity": float(echo.analysis.restricted_process_fidelity),
        "best_fit_restricted_process_fidelity": float(echo.analysis.best_fit_restricted_process_fidelity),
        "max_residual_z_error_rad": float(np.max([item["residual_z_error_rad"] for item in echo_rows])),
    }
    return direct_row, echo_row


def run_case(
    context: CaseContext,
    *,
    direct_starts: int,
    direct_maxiter: int,
    segment_maxiter: int,
    echo_ideal_maxiter: int,
    echo_gaussian_maxiter: int,
    refocus_maxiter: int,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []

    direct = optimize_square_multitone(
        context.model,
        context.spec,
        context.run_config,
        n_starts=direct_starts,
        maxiter=direct_maxiter,
        random_seed=stable_seed(context.request.case_id, 1),
        label_prefix=f"{context.request.case_id}_direct",
    )
    direct_row, direct_artifact = row_from_analysis(
        context,
        construction="full_shared_line",
        comparison_spec=context.spec,
        comparison_target_operator=context.target_operator,
        analysis=direct.validation,
        full_operator=np.asarray(direct.validation.full_operator, dtype=np.complex128),
        extra={
            "total_gate_duration_ns": float(context.duration_s * 1.0e9),
            "optimizer_success": bool(direct.success),
            "optimizer_message": str(direct.message),
            "corrections_vector": corrections_to_vector(direct.corrections, n_active=context.request.n_active).tolist(),
        },
    )
    direct_artifact.update(
        {
            "target_spec": json_ready(context.spec.metadata),
            "tone_rows": [tone.as_dict() for tone in direct.tone_specs],
            "waveform_samples": channel_waveform_samples(direct.validation.compiled),
        }
    )
    save_case_artifact(context.request.case_id, "full_shared_line", direct_artifact)
    rows.append(direct_row)

    reduced_operator = reduced_blockwise_operator(
        context.model,
        direct.validation.compiled,
        direct.waveform,
        context.run_config,
        levels=context.levels,
    )
    reduced_full = embed_restricted_operator_in_full(reduced_operator, context.model, context.levels)
    reduced_analysis = analysis_against_target(context, reduced_full, metadata={"construction": "blockwise_exact_reduced"})
    reduced_row, reduced_artifact = row_from_analysis(
        context,
        construction="blockwise_exact_reduced",
        comparison_spec=context.spec,
        comparison_target_operator=context.target_operator,
        analysis=reduced_analysis,
        full_operator=reduced_full,
        extra={
            "total_gate_duration_ns": float(context.duration_s * 1.0e9),
            "reduced_vs_full_restricted_process_fidelity": float(process_fidelity(np.asarray(direct.validation.restricted_operator), reduced_operator)),
            "reduced_vs_full_restricted_average_gate_fidelity": float(average_gate_fidelity(np.asarray(direct.validation.restricted_operator), reduced_operator)),
        },
    )
    save_case_artifact(context.request.case_id, "blockwise_exact_reduced", reduced_artifact)
    rows.append(reduced_row)

    magnus_operator = np.zeros_like(context.target_operator)
    for block_index, block in enumerate(
        magnus_effective_blocks(direct.tone_specs, model=context.model, run_config=context.run_config, levels=context.levels)
    ):
        magnus_operator[2 * block_index : 2 * block_index + 2, 2 * block_index : 2 * block_index + 2] = block
    magnus_full = embed_restricted_operator_in_full(magnus_operator, context.model, context.levels)
    magnus_analysis = analysis_against_target(context, magnus_full, metadata={"construction": "magnus_effective"})
    magnus_row, magnus_artifact = row_from_analysis(
        context,
        construction="magnus_effective",
        comparison_spec=context.spec,
        comparison_target_operator=context.target_operator,
        analysis=magnus_analysis,
        full_operator=magnus_full,
        extra={
            "total_gate_duration_ns": float(context.duration_s * 1.0e9),
            "magnus_vs_full_restricted_process_fidelity": float(process_fidelity(np.asarray(direct.validation.restricted_operator), magnus_operator)),
        },
    )
    save_case_artifact(context.request.case_id, "magnus_effective", magnus_artifact)
    rows.append(magnus_row)

    _waveform, ideal_tones, _compiled = build_square_multitone_compiled(
        context.model,
        context.spec,
        context.run_config,
        corrections=None,
        label=f"{context.request.case_id}_ideal",
    )
    decoupled_operator = decoupled_block_operator(ideal_tones, levels=context.levels, duration_s=context.duration_s)
    decoupled_full = embed_restricted_operator_in_full(decoupled_operator, context.model, context.levels)
    decoupled_analysis = analysis_against_target(context, decoupled_full, metadata={"construction": "decoupled_block"})
    decoupled_row, decoupled_artifact = row_from_analysis(
        context,
        construction="decoupled_block",
        comparison_spec=context.spec,
        comparison_target_operator=context.target_operator,
        analysis=decoupled_analysis,
        full_operator=decoupled_full,
        extra={
            "total_gate_duration_ns": float(context.duration_s * 1.0e9),
            "same_block_population_mean": 1.0,
            "leakage_outside_target_mean": 0.0,
        },
    )
    save_case_artifact(context.request.case_id, "decoupled_block", decoupled_artifact)
    rows.append(decoupled_row)

    direct_total = optimize_square_multitone(
        context.model,
        context.spec,
        context.matched_run_config,
        n_starts=max(1, direct_starts - 1),
        maxiter=direct_maxiter,
        random_seed=stable_seed(context.request.case_id, 2),
        label_prefix=f"{context.request.case_id}_direct_total",
    )
    direct_total_row, direct_total_artifact = row_from_analysis(
        context,
        construction="full_shared_line_total_matched",
        comparison_spec=context.spec,
        comparison_target_operator=context.target_operator,
        analysis=direct_total.validation,
        full_operator=np.asarray(direct_total.validation.full_operator, dtype=np.complex128),
        extra={
            "total_gate_duration_ns": float(context.matched_total_duration_s * 1.0e9),
            "optimizer_success": bool(direct_total.success),
            "optimizer_message": str(direct_total.message),
        },
    )
    direct_total_artifact.update({"waveform_samples": channel_waveform_samples(direct_total.validation.compiled)})
    save_case_artifact(context.request.case_id, "full_shared_line_total_matched", direct_total_artifact)
    rows.append(direct_total_row)

    spec_1, spec_2 = echo_segment_specs(context.spec, 0.0)
    seed_1 = optimize_square_multitone(
        context.model,
        spec_1,
        run_config_for_duration(context, 0.5 * context.duration_s),
        n_starts=1,
        maxiter=segment_maxiter,
        random_seed=stable_seed(context.request.case_id, 10),
        label_prefix=f"{context.request.case_id}_echo_seed_1",
    )
    seed_2 = optimize_square_multitone(
        context.model,
        spec_2,
        run_config_for_duration(context, 0.5 * context.duration_s),
        n_starts=1,
        maxiter=segment_maxiter,
        random_seed=stable_seed(context.request.case_id, 11),
        label_prefix=f"{context.request.case_id}_echo_seed_2",
    )

    replay_payload = build_echo_full_operator(
        context,
        corrections_1=seed_1.corrections,
        corrections_2=seed_2.corrections,
        eta=0.0,
        refocus_mode="ideal_instantaneous",
        label=f"{context.request.case_id}_echo_replay",
    )
    replay_analysis = analysis_against_target(context, replay_payload["full_operator"], metadata={"construction": "echo_replay_ideal", **replay_payload["metadata"]})
    replay_row, replay_artifact = row_from_analysis(
        context,
        construction="echo_replay_ideal",
        comparison_spec=context.spec,
        comparison_target_operator=context.target_operator,
        analysis=replay_analysis,
        full_operator=np.asarray(replay_payload["full_operator"], dtype=np.complex128),
        extra={
            "total_gate_duration_ns": float(replay_payload["metadata"]["total_gate_duration_s"] * 1.0e9),
            "echo_eta": 0.0,
            "optimizer_kind": "segmentwise_replay",
        },
    )
    replay_artifact.update(replay_payload["segment_artifacts"])
    save_case_artifact(context.request.case_id, "echo_replay_ideal", replay_artifact)
    rows.append(replay_row)

    echo_x0 = np.concatenate(
        [
            corrections_to_vector(seed_1.corrections, n_active=context.request.n_active),
            corrections_to_vector(seed_2.corrections, n_active=context.request.n_active),
            np.asarray([0.0], dtype=float),
        ]
    )
    echo_ideal = optimize_echo(
        context,
        construction="echo_opt_ideal",
        refocus_mode="ideal_instantaneous",
        x0=echo_x0,
        maxiter=echo_ideal_maxiter,
        maxfev=180,
    )
    echo_ideal_row, echo_ideal_artifact = row_from_analysis(
        context,
        construction="echo_opt_ideal",
        comparison_spec=context.spec,
        comparison_target_operator=context.target_operator,
        analysis=echo_ideal.analysis,
        full_operator=echo_ideal.full_operator,
        extra={
            "total_gate_duration_ns": float(echo_ideal.metadata["total_gate_duration_s"] * 1.0e9),
            "echo_eta": float(echo_ideal.eta),
            **echo_ideal.optimizer_payload,
        },
    )
    echo_ideal_artifact.update(echo_ideal.segment_artifacts)
    save_case_artifact(context.request.case_id, "echo_opt_ideal", echo_ideal_artifact)
    rows.append(echo_ideal_row)

    gaussian_x0 = np.concatenate(
        [
            corrections_to_vector(echo_ideal.corrections_1, n_active=context.request.n_active),
            corrections_to_vector(echo_ideal.corrections_2, n_active=context.request.n_active),
            np.asarray([echo_ideal.eta], dtype=float),
        ]
    )
    echo_gaussian = optimize_echo(
        context,
        construction="echo_opt_gaussian",
        refocus_mode="gaussian_pi",
        x0=np.asarray(gaussian_x0, dtype=float),
        maxiter=echo_gaussian_maxiter,
        maxfev=110,
    )
    echo_gaussian_row, echo_gaussian_artifact = row_from_analysis(
        context,
        construction="echo_opt_gaussian",
        comparison_spec=context.spec,
        comparison_target_operator=context.target_operator,
        analysis=echo_gaussian.analysis,
        full_operator=echo_gaussian.full_operator,
        extra={
            "total_gate_duration_ns": float(echo_gaussian.metadata["total_gate_duration_s"] * 1.0e9),
            "echo_eta": float(echo_gaussian.eta),
            **echo_gaussian.optimizer_payload,
        },
    )
    echo_gaussian_artifact.update(echo_gaussian.segment_artifacts)
    echo_gaussian_artifact["waveform_samples"] = echo_gaussian.segment_artifacts.get("sequence_waveform_samples")
    save_case_artifact(context.request.case_id, "echo_opt_gaussian", echo_gaussian_artifact)
    rows.append(echo_gaussian_row)

    if use_multitone_pi_benchmark(context.request):
        refocus = optimize_refocus_pulse(context, maxiter=refocus_maxiter, n_starts=1)
        refocus_row, refocus_artifact = row_from_analysis(
            context,
            construction="refocus_multitone_pi",
            comparison_spec=refocus.spec,
            comparison_target_operator=refocus.target_operator,
            analysis=refocus.validation,
            full_operator=np.asarray(refocus.optimization.validation.full_operator, dtype=np.complex128),
            extra={
                "total_gate_duration_ns": float(PI_PULSE_DURATION_S * 1.0e9),
                "optimizer_success": bool(refocus.optimization.success),
                "optimizer_message": str(refocus.optimization.message),
            },
        )
        refocus_artifact.update(
            {
                "waveform_samples": channel_waveform_samples(refocus.optimization.validation.compiled),
                "tone_rows": [tone.as_dict() for tone in refocus.optimization.tone_specs],
            }
        )
        save_case_artifact(context.request.case_id, "refocus_multitone_pi", refocus_artifact)
        rows.append(refocus_row)

        echo_multitone = optimize_echo(
            context,
            construction="echo_opt_multitone_pi",
            refocus_mode="multitone_pi",
            refocus_waveform=refocus.optimization.waveform,
            x0=np.asarray(gaussian_x0, dtype=float),
            maxiter=max(6, echo_gaussian_maxiter - 1),
            maxfev=90,
        )
        echo_multitone_row, echo_multitone_artifact = row_from_analysis(
            context,
            construction="echo_opt_multitone_pi",
            comparison_spec=context.spec,
            comparison_target_operator=context.target_operator,
            analysis=echo_multitone.analysis,
            full_operator=echo_multitone.full_operator,
            extra={
                "total_gate_duration_ns": float(echo_multitone.metadata["total_gate_duration_s"] * 1.0e9),
                "echo_eta": float(echo_multitone.eta),
                **echo_multitone.optimizer_payload,
            },
        )
        echo_multitone_artifact.update(echo_multitone.segment_artifacts)
        echo_multitone_artifact["waveform_samples"] = echo_multitone.segment_artifacts.get("sequence_waveform_samples")
        save_case_artifact(context.request.case_id, "echo_opt_multitone_pi", echo_multitone_artifact)
        rows.append(echo_multitone_row)

    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=str, default="", help="Optional case_id filter.")
    parser.add_argument("--direct-starts", type=int, default=2)
    parser.add_argument("--direct-maxiter", type=int, default=18)
    parser.add_argument("--segment-maxiter", type=int, default=10)
    parser.add_argument("--echo-ideal-maxiter", type=int, default=14)
    parser.add_argument("--echo-gaussian-maxiter", type=int, default=8)
    parser.add_argument("--refocus-maxiter", type=int, default=12)
    parser.add_argument("--skip-plots", action="store_true")
    args = parser.parse_args()

    CASE_DIR.mkdir(parents=True, exist_ok=True)
    WAVEFORM_DIR.mkdir(parents=True, exist_ok=True)

    requests = case_requests()
    if args.case:
        requests = [request for request in requests if request.case_id == args.case]
        if not requests:
            raise SystemExit(f"No case matched '{args.case}'.")

    save_json(METRICS_PATH, {"metric_definitions": metric_definitions()}, description="Metric definitions for the rigorous echo follow-up study.")

    start_time = time.perf_counter()
    all_rows: list[dict[str, Any]] = []
    for request in requests:
        print(f"[case] {request.case_id}")
        context = build_context(request)
        case_rows = run_case(
            context,
            direct_starts=args.direct_starts,
            direct_maxiter=args.direct_maxiter,
            segment_maxiter=args.segment_maxiter,
            echo_ideal_maxiter=args.echo_ideal_maxiter,
            echo_gaussian_maxiter=args.echo_gaussian_maxiter,
            refocus_maxiter=args.refocus_maxiter,
        )
        all_rows.extend(case_rows)
        for row in case_rows:
            print(
                "  {construction}: fid={restricted_average_gate_fidelity:.6f} pf={restricted_process_fidelity:.6f} "
                "resZ={max_residual_z_error_rad:.6f} probe={probe_fidelity_mean:.6f}".format(**row)
            )

    df = pd.DataFrame(all_rows)
    df.to_csv(CSV_PATH, index=False)

    if not args.skip_plots and not df.empty:
        figure_duration_tradeoff(df)
        figure_residual_z(df)
        figure_echo_comparison(df)
        figure_echo_tradeoff(df)
        figure_refocus_benchmark(df)

    representative_case = next(
        (
            item
            for item in requests
            if item.model_variant == "chi_plus_chiprime" and item.family == "aligned_x" and item.n_active == 3 and item.chi_t_over_2pi == 5.0
        ),
        requests[0] if requests else None,
    )
    rep_direct = None
    rep_echo = None
    if representative_case is not None and len(requests) > 1:
        rep_direct, rep_echo = representative_convergence(
            build_context(representative_case),
            direct_budget=max(24, args.direct_maxiter + 6),
            echo_budget=max(12, args.echo_gaussian_maxiter + 4),
        )

    summary = build_summary(all_rows)
    validation = validation_payload(all_rows, representative_direct=rep_direct, representative_echo=rep_echo)
    save_json(
        RESULTS_PATH,
        {"case_rows": all_rows, "runtime_s": float(time.perf_counter() - start_time)},
        description="Machine-readable result table for the rigorous echoed-ansatz follow-up study.",
    )
    save_json(SUMMARY_PATH, summary, description="Headline summary for the rigorous echoed-ansatz follow-up study.")
    save_json(VALIDATION_PATH, validation, description="Validation payload for the rigorous echoed-ansatz follow-up study.")

    lines = [
        f"# Summary: {STUDY_DIR.name}",
        "",
        "## Construction Means",
    ]
    for row in summary["construction_summary"]:
        lines.append(
            "- {construction}: mean fidelity {mean_fidelity:.6f}, best fidelity {best_fidelity:.6f}, mean max residual-Z {mean_residual_z:.6f} rad, mean probe fidelity {mean_probe:.6f}.".format(
                **row
            )
        )
    if summary["best_overall"] is not None:
        lines.extend(["", "## Best Overall"])
        lines.append(
            "- {construction} on {case_id}: restricted average fidelity {restricted_average_gate_fidelity:.6f}, restricted process fidelity {restricted_process_fidelity:.6f}, max residual-Z {max_residual_z_error_rad:.6f} rad.".format(
                **summary["best_overall"]
            )
        )
    MARKDOWN_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
