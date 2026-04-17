"""Extension pass for tuned-set mapping and echoed refocusing follow-up."""

from __future__ import annotations

import argparse
import time
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import brentq

from common import (
    ARTIFACTS_DIR,
    DATA_DIR,
    DEFAULT_DT,
    FIGURES_DIR,
    IDEAL_X_PI,
    PI_PULSE_DURATION_S,
    TWO_PI,
    CaseRequest,
    TargetSpec,
    apply_plot_style,
    average_gate_fidelity,
    build_frame,
    build_model,
    build_square_multitone_waveform,
    build_target_operator,
    channel_waveform_samples,
    compile_pulse_sequence,
    corrections_to_vector,
    duration_from_chi_t,
    embed_qubit_operator,
    logical_levels,
    make_gaussian_qubit_rotation_pulse,
    make_run_config,
    manifold_transition_frequencies_rad_s,
    optimize_square_multitone,
    process_fidelity,
    reduced_blockwise_operator,
    restricted_operator_from_full,
    save_json,
    save_waveform_npz,
    shift_pulse,
    simulate_full_operator_on_logical_inputs,
)
from run_study import row_from_operator, save_case_artifact


CASE_DIR = ARTIFACTS_DIR / "cases"
WAVEFORM_DIR = ARTIFACTS_DIR / "waveforms"

TUNED_MAP_PATH = DATA_DIR / "extension_tuned_set_map.json"
CHECKPOINT_PATH = DATA_DIR / "extension_checkpoint_summary.json"
ECHO_PATH = DATA_DIR / "extension_echo_summary.json"
RESULTS_PATH = DATA_DIR / "extension_results.json"
SUMMARY_PATH = DATA_DIR / "extension_summary.json"

TUNED_MAP_FIGURE_BASENAME = "extension_tuned_set_map"
CHECKPOINT_FIGURE_BASENAME = "extension_exact_checkpoint_comparison"
ECHO_FIGURE_BASENAME = "extension_echo_followup_comparison"

EQUAL_ANGLE_X_THETA_RAD = 0.50 * np.pi
ROOT_SEARCH_X_MIN = 0.25
ROOT_SEARCH_X_MAX = 12.0
ROOT_SEARCH_SAMPLES = 8000
TUNED_MAP_X_MAX = 10.0
TUNED_MAP_X_POINTS = 451
TUNED_MAP_DELTA_POINTS = 361
OFF_TUNED_OFFSET_CHI_T = 0.10
CHECKPOINT_N_STARTS = 3
CHECKPOINT_MAXITER = 60
ECHO_HALF_N_STARTS = 3
ECHO_HALF_MAXITER = 60
REFOCUS_N_STARTS = 3
REFOCUS_MAXITER = 45
MANIFOLD_AWARE_PI_DURATION_S = 80.0e-9


def equal_angle_aligned_x_target(
    *,
    theta_rad: float,
    n_active: int = 2,
    family: str,
    description: str,
) -> TargetSpec:
    theta_values = tuple(float(theta_rad) for _ in range(int(n_active)))
    phi_values = tuple(0.0 for _ in range(int(n_active)))
    metadata = {
        "family": str(family),
        "description": str(description),
        "n_active": int(n_active),
        "theta_values_rad": list(theta_values),
        "phi_values_rad": list(phi_values),
        "symmetry": "aligned_x_equal_angle",
    }
    return TargetSpec(
        family=str(family),
        theta_values=theta_values,
        phi_values=phi_values,
        metadata=metadata,
    )


def half_target(spec: TargetSpec) -> TargetSpec:
    theta_values = tuple(float(value / 2.0) for value in spec.theta_values)
    metadata = {
        **dict(spec.metadata),
        "description": f"Half-angle target derived from {spec.family}.",
        "parent_family": str(spec.family),
        "theta_values_rad": list(theta_values),
    }
    return TargetSpec(
        family=f"{spec.family}_half",
        theta_values=theta_values,
        phi_values=tuple(float(value) for value in spec.phi_values),
        metadata=metadata,
    )


def k_kernel_dimensionless(x: np.ndarray | float) -> np.ndarray:
    arr = np.asarray(x, dtype=float)
    out = np.empty_like(arr)
    small = np.abs(arr) <= 1.0e-8
    out[small] = arr[small] / 6.0
    xx = arr[~small]
    out[~small] = (xx - np.sin(xx)) / (xx * xx)
    return out


def l_kernel_dimensionless(x: np.ndarray | float, delta: np.ndarray | float) -> np.ndarray:
    x_arr = np.asarray(x, dtype=float)
    delta_arr = np.asarray(delta, dtype=float)
    out = np.empty(np.broadcast_shapes(x_arr.shape, delta_arr.shape), dtype=float)
    xx = np.broadcast_to(x_arr, out.shape)
    dd = np.broadcast_to(delta_arr, out.shape)
    small = np.abs(xx) <= 1.0e-8
    out[small] = np.cos(dd[small]) * xx[small] / 6.0
    large_x = xx[~small]
    large_d = dd[~small]
    out[~small] = (
        large_x * (np.cos(large_d) + np.cos(large_x + large_d))
        + 2.0 * np.sin(large_d)
        - 2.0 * np.sin(large_x + large_d)
    ) / (large_x * large_x)
    return out


def equal_amplitude_root_function(x: float) -> float:
    value = float(x)
    return value * np.cos(value) - np.sin(value)


def find_equal_amplitude_roots(*, x_min: float, x_max: float, n_samples: int) -> list[float]:
    xs = np.linspace(float(x_min), float(x_max), int(n_samples))
    values = xs * np.cos(xs) - np.sin(xs)
    roots: list[float] = []
    for left, right, left_value, right_value in zip(
        xs[:-1],
        xs[1:],
        values[:-1],
        values[1:],
        strict=True,
    ):
        if left_value == 0.0:
            candidate = float(left)
        elif left_value * right_value > 0.0:
            continue
        else:
            candidate = float(brentq(equal_amplitude_root_function, float(left), float(right)))
        if candidate <= 1.0e-6:
            continue
        if not roots or abs(candidate - roots[-1]) > 1.0e-5:
            roots.append(candidate)
    return roots


def build_tuned_map() -> tuple[dict[str, Any], np.ndarray, np.ndarray, np.ndarray, list[float]]:
    x_values = np.linspace(ROOT_SEARCH_X_MIN, TUNED_MAP_X_MAX, TUNED_MAP_X_POINTS)
    delta_values = np.linspace(-np.pi, np.pi, TUNED_MAP_DELTA_POINTS)
    xx, dd = np.meshgrid(x_values, delta_values)
    k_bar = k_kernel_dimensionless(xx)
    l_bar = l_kernel_dimensionless(xx, dd)
    equal_amplitude_residual = np.abs(k_bar - l_bar)
    required_ratio = np.divide(
        k_bar,
        l_bar,
        out=np.full_like(k_bar, np.nan),
        where=np.abs(l_bar) > 1.0e-10,
    )
    roots = find_equal_amplitude_roots(
        x_min=ROOT_SEARCH_X_MIN,
        x_max=ROOT_SEARCH_X_MAX,
        n_samples=ROOT_SEARCH_SAMPLES,
    )
    payload = {
        "x_values": x_values.tolist(),
        "delta_values_rad": delta_values.tolist(),
        "chi_t_over_2pi_values": (x_values / TWO_PI).tolist(),
        "equal_amplitude_residual": equal_amplitude_residual.tolist(),
        "required_amplitude_ratio": required_ratio.tolist(),
        "equal_amplitude_roots": [
            {
                "x_root": float(root),
                "chi_t_over_2pi": float(root / TWO_PI),
                "description": "Nontrivial equal-angle aligned-x tuned root solving x cos(x) - sin(x) = 0.",
            }
            for root in roots
        ],
        "notes": {
            "condition": "Equal-amplitude aligned-x tuned loci satisfy K(Delta,T) = L(Delta,T,0), equivalently x cos(x) - sin(x) = 0 with x = |Delta| T.",
            "generic_statement": "Away from these lower-dimensional tuned loci, the equal-amplitude residual remains nonzero, and unequal-angle targets miss the tuned set because lambda_1/lambda_0 != 1.",
        },
    }
    return payload, x_values, delta_values, equal_amplitude_residual, roots


def figure_tuned_map(
    x_values: np.ndarray,
    delta_values: np.ndarray,
    residual: np.ndarray,
    roots: list[float],
) -> None:
    apply_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.4))

    image = axes[0].imshow(
        np.log10(np.maximum(residual, 1.0e-8)),
        origin="lower",
        aspect="auto",
        extent=(x_values[0] / TWO_PI, x_values[-1] / TWO_PI, delta_values[0] / np.pi, delta_values[-1] / np.pi),
        cmap="cividis",
    )
    axes[0].contour(
        x_values / TWO_PI,
        delta_values / np.pi,
        residual,
        levels=[5.0e-3],
        colors=["white"],
        linewidths=1.2,
    )
    axes[0].scatter(
        np.asarray(roots, dtype=float) / TWO_PI,
        np.zeros(len(roots), dtype=float),
        color="#EE6677",
        s=20.0,
        label="Aligned-x tuned roots",
    )
    axes[0].set_xlabel(r"$|\chi| T / 2\pi$")
    axes[0].set_ylabel(r"$\delta / \pi$")
    axes[0].set_title("Equal-amplitude tuned-set residual")
    axes[0].legend(frameon=False, loc="upper right")
    colorbar = fig.colorbar(image, ax=axes[0])
    colorbar.set_label(r"$\log_{10}|K-L|$")

    aligned_slice = residual[np.argmin(np.abs(delta_values)), :]
    axes[1].plot(
        x_values / TWO_PI,
        aligned_slice,
        color="#4477AA",
        linewidth=2.0,
        label=r"Aligned $x$ slice ($\delta = 0$)",
    )
    axes[1].scatter(
        np.asarray(roots, dtype=float) / TWO_PI,
        np.full(len(roots), 1.0e-8, dtype=float),
        color="#EE6677",
        s=24.0,
        label="Tuned roots",
    )
    axes[1].set_yscale("log")
    axes[1].set_xlabel(r"$|\chi| T / 2\pi$")
    axes[1].set_ylabel(r"$|K-L|$")
    axes[1].set_title("Aligned-x accidental tuned roots")
    axes[1].grid(True, alpha=0.25)
    axes[1].legend(frameon=False, loc="upper right")

    fig.tight_layout()
    fig.savefig(FIGURES_DIR / f"{TUNED_MAP_FIGURE_BASENAME}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / f"{TUNED_MAP_FIGURE_BASENAME}.pdf", bbox_inches="tight")
    plt.close(fig)


def extension_requests(tuned_chi_t_over_2pi: float) -> list[CaseRequest]:
    return [
        CaseRequest(
            model_variant="chi_only",
            include_chi_prime=False,
            family="equal_angle_aligned_x_tuned",
            n_active=2,
            chi_t_over_2pi=float(tuned_chi_t_over_2pi),
        ),
        CaseRequest(
            model_variant="chi_only",
            include_chi_prime=False,
            family="equal_angle_aligned_x_off_minus",
            n_active=2,
            chi_t_over_2pi=float(tuned_chi_t_over_2pi - OFF_TUNED_OFFSET_CHI_T),
        ),
        CaseRequest(
            model_variant="chi_only",
            include_chi_prime=False,
            family="equal_angle_aligned_x_off_plus",
            n_active=2,
            chi_t_over_2pi=float(tuned_chi_t_over_2pi + OFF_TUNED_OFFSET_CHI_T),
        ),
        CaseRequest(
            model_variant="chi_plus_chiprime",
            include_chi_prime=True,
            family="equal_angle_aligned_x_tuned_chiprime",
            n_active=2,
            chi_t_over_2pi=float(tuned_chi_t_over_2pi),
        ),
    ]


def target_for_request(request: CaseRequest) -> TargetSpec:
    return equal_angle_aligned_x_target(
        theta_rad=EQUAL_ANGLE_X_THETA_RAD,
        n_active=request.n_active,
        family=request.family,
        description="Equal-angle aligned-x target used for tuned-set checkpoints.",
    )


def checkpoint_analytic_terms(
    *,
    model: Any,
    levels: tuple[int, ...],
    duration_s: float,
    tone_specs: tuple[Any, ...],
) -> dict[str, float]:
    freqs = manifold_transition_frequencies_rad_s(model, levels, build_frame(model))
    delta_rad_s = float(abs(freqs[1] - freqs[0]))
    x_value = float(delta_rad_s * duration_s)
    amp_map = {int(spec.manifold): float(spec.amp_rad_s) for spec in tone_specs}
    phase_map = {int(spec.manifold): float(spec.phase_rad) for spec in tone_specs}
    lambda_0 = float(amp_map[int(levels[0])])
    lambda_1 = float(amp_map[int(levels[1])])
    delta_phase = float(phase_map[int(levels[0])] - phase_map[int(levels[1])])
    k_value = float((delta_rad_s * duration_s - np.sin(delta_rad_s * duration_s)) / (delta_rad_s * delta_rad_s * duration_s))
    l_value = float(
        (
            delta_rad_s
            * duration_s
            * (np.cos(delta_phase) + np.cos(delta_rad_s * duration_s + delta_phase))
            + 2.0 * np.sin(delta_phase)
            - 2.0 * np.sin(delta_rad_s * duration_s + delta_phase)
        )
        / (delta_rad_s * delta_rad_s * duration_s)
    )
    zeta_0 = float(-lambda_1 * lambda_1 * k_value + lambda_0 * lambda_1 * l_value)
    zeta_1 = float(lambda_0 * lambda_0 * k_value - lambda_0 * lambda_1 * l_value)
    return {
        "analytic_delta_rad_s": float(delta_rad_s),
        "analytic_x_value": float(x_value),
        "analytic_phase_difference_rad": float(delta_phase),
        "analytic_amplitude_ratio": float(lambda_1 / lambda_0) if abs(lambda_0) > 1.0e-12 else float("nan"),
        "analytic_k_value": float(k_value),
        "analytic_l_value": float(l_value),
        "analytic_zeta_0": float(zeta_0),
        "analytic_zeta_1": float(zeta_1),
        "analytic_equal_amplitude_residual": float(abs(k_value - l_value)),
    }


def save_extension_artifact(
    request: CaseRequest,
    *,
    spec: TargetSpec,
    construction: str,
    target_operator: np.ndarray,
    actual_operator: np.ndarray,
    tone_specs: tuple[Any, ...] | None = None,
    corrections_vector: np.ndarray | None = None,
    waveform_samples: dict[str, Any] | None = None,
    metadata: dict[str, Any] | None = None,
) -> Path:
    path = CASE_DIR / f"extension_{request.case_id}_{construction}.json"
    save_case_artifact(
        path,
        request=request,
        spec=spec,
        construction=construction,
        target_operator=target_operator,
        actual_operator=actual_operator,
        tone_specs=tone_specs,
        corrections_vector=corrections_vector,
        waveform_samples=waveform_samples,
        metadata=metadata,
    )
    return path


def run_checkpoint_case(request: CaseRequest, *, n_starts: int, maxiter: int) -> dict[str, Any]:
    spec = target_for_request(request)
    model = build_model(include_chi_prime=request.include_chi_prime, n_active=request.n_active)
    duration_s = duration_from_chi_t(request.chi_t_over_2pi)
    run_config = make_run_config(model, n_active=request.n_active, duration_s=duration_s)
    levels = logical_levels(request.n_active)
    target_operator = build_target_operator(spec, levels)
    optimization = optimize_square_multitone(
        model,
        spec,
        run_config,
        n_starts=n_starts,
        maxiter=maxiter,
        random_seed=14001 + int(round(1000.0 * request.chi_t_over_2pi)),
        label_prefix=f"extension_{request.case_id}",
    )
    full_validation = optimization.validation
    full_restricted = np.asarray(full_validation.restricted_operator, dtype=np.complex128)
    analytic_terms = checkpoint_analytic_terms(
        model=model,
        levels=levels,
        duration_s=duration_s,
        tone_specs=optimization.tone_specs,
    )
    corrections_vector = corrections_to_vector(optimization.corrections, n_active=request.n_active)
    waveform_samples = channel_waveform_samples(full_validation.compiled)
    full_row = row_from_operator(
        request,
        construction="full_shared_line",
        target_operator=target_operator,
        actual_operator=full_restricted,
        validation=full_validation,
        extra={
            "extension_scope": "exact_checkpoint",
            "optimization_success": bool(optimization.success),
            "optimization_message": str(optimization.message),
            "n_optimizer_history": int(len(optimization.history)),
            "corrections_vector": corrections_vector.tolist(),
            **analytic_terms,
        },
    )
    full_artifact = save_extension_artifact(
        request,
        spec=spec,
        construction="full_shared_line",
        target_operator=target_operator,
        actual_operator=full_restricted,
        tone_specs=optimization.tone_specs,
        corrections_vector=corrections_vector,
        waveform_samples=waveform_samples,
        metadata={
            "extension_scope": "exact_checkpoint",
            **analytic_terms,
        },
    )
    save_waveform_npz(WAVEFORM_DIR / f"extension_{request.case_id}_full_shared_line.npz", waveform_samples)

    reduced_operator = reduced_blockwise_operator(
        model,
        full_validation.compiled,
        optimization.waveform,
        run_config,
        levels=levels,
    )
    reduced_row = row_from_operator(
        request,
        construction="blockwise_exact_reduced",
        target_operator=target_operator,
        actual_operator=reduced_operator,
        extra={
            "extension_scope": "exact_checkpoint",
            "source_case": str(request.case_id),
            "reduced_vs_full_restricted_process_fidelity": float(process_fidelity(full_restricted, reduced_operator)),
            "reduced_vs_full_restricted_average_gate_fidelity": float(average_gate_fidelity(full_restricted, reduced_operator)),
        },
    )
    reduced_artifact = save_extension_artifact(
        request,
        spec=spec,
        construction="blockwise_exact_reduced",
        target_operator=target_operator,
        actual_operator=reduced_operator,
        tone_specs=optimization.tone_specs,
        corrections_vector=corrections_vector,
        metadata={
            "extension_scope": "exact_checkpoint",
            "source_case": str(request.case_id),
        },
    )
    return {
        "request": request,
        "spec": spec,
        "model": model,
        "duration_s": float(duration_s),
        "levels": levels,
        "target_operator": target_operator,
        "optimization": optimization,
        "full_row": full_row,
        "reduced_row": reduced_row,
        "full_artifact": str(full_artifact),
        "reduced_artifact": str(reduced_artifact),
    }


def run_echo_followup(bundle: dict[str, Any]) -> list[dict[str, Any]]:
    request: CaseRequest = bundle["request"]
    spec: TargetSpec = bundle["spec"]
    model = bundle["model"]
    duration_s = float(bundle["duration_s"])
    levels = bundle["levels"]
    target_operator = bundle["target_operator"]

    frame = build_frame(model)
    half_spec = half_target(spec)
    half_run_config = make_run_config(model, n_active=request.n_active, duration_s=0.5 * duration_s)
    half_optimization = optimize_square_multitone(
        model,
        half_spec,
        half_run_config,
        n_starts=ECHO_HALF_N_STARTS,
        maxiter=ECHO_HALF_MAXITER,
        random_seed=22001 + int(round(1000.0 * request.chi_t_over_2pi)),
        label_prefix=f"extension_{request.case_id}_half",
    )
    half_corrections = corrections_to_vector(half_optimization.corrections, n_active=request.n_active)

    echo_rows: list[dict[str, Any]] = []
    half_full_operator = np.asarray(half_optimization.validation.full_operator, dtype=np.complex128)
    ideal_x_full = embed_qubit_operator(IDEAL_X_PI, n_cav=int(model.n_cav))
    ideal_echo_full = ideal_x_full @ half_full_operator @ ideal_x_full @ half_full_operator
    ideal_echo_restricted = restricted_operator_from_full(ideal_echo_full, model, levels)
    ideal_echo_row = row_from_operator(
        request,
        construction="echo_ideal_instantaneous",
        target_operator=target_operator,
        actual_operator=ideal_echo_restricted,
        extra={
            "extension_scope": "echo_followup",
            "echo_half_duration_ns": float(0.5 * duration_s * 1.0e9),
            "echo_pi_duration_ns": 0.0,
            "echo_total_duration_ns": float(duration_s * 1.0e9),
            "refocusing_convention": "ideal_instantaneous_x",
            "half_corrections_vector": half_corrections.tolist(),
        },
    )
    save_extension_artifact(
        request,
        spec=spec,
        construction="echo_ideal_instantaneous",
        target_operator=target_operator,
        actual_operator=ideal_echo_restricted,
        tone_specs=half_optimization.tone_specs,
        corrections_vector=half_corrections,
        metadata={
            "extension_scope": "echo_followup",
            "refocusing_convention": "ideal_instantaneous_x",
            "half_corrections_vector": half_corrections.tolist(),
        },
    )
    echo_rows.append(ideal_echo_row)

    half_pulse = half_optimization.waveform.pulse
    gaussian_pi_pulse = make_gaussian_qubit_rotation_pulse(
        model,
        frame,
        theta=np.pi,
        phase=0.0,
        duration_s=PI_PULSE_DURATION_S,
        manifold_level=0,
        label=f"extension_{request.case_id}_gaussian_pi",
    )
    gaussian_pulses = [
        shift_pulse(half_pulse, t0=0.0, label=f"extension_{request.case_id}_half_1"),
        shift_pulse(gaussian_pi_pulse, t0=0.5 * duration_s, label=f"extension_{request.case_id}_gaussian_pi_1"),
        shift_pulse(half_pulse, t0=0.5 * duration_s + PI_PULSE_DURATION_S, label=f"extension_{request.case_id}_half_2"),
        shift_pulse(gaussian_pi_pulse, t0=duration_s + PI_PULSE_DURATION_S, label=f"extension_{request.case_id}_gaussian_pi_2"),
    ]
    gaussian_total_duration = duration_s + 2.0 * PI_PULSE_DURATION_S
    gaussian_compiled = compile_pulse_sequence(
        gaussian_pulses,
        dt_s=DEFAULT_DT,
        total_duration_s=gaussian_total_duration,
    )
    gaussian_full = simulate_full_operator_on_logical_inputs(
        model,
        gaussian_compiled,
        frame=frame,
        drive_ops={"qubit": "qubit"},
        levels=levels,
    )
    gaussian_restricted = restricted_operator_from_full(gaussian_full, model, levels)
    gaussian_samples = channel_waveform_samples(gaussian_compiled)
    gaussian_row = row_from_operator(
        request,
        construction="echo_finite_gaussian",
        target_operator=target_operator,
        actual_operator=gaussian_restricted,
        extra={
            "extension_scope": "echo_followup",
            "echo_half_duration_ns": float(0.5 * duration_s * 1.0e9),
            "echo_pi_duration_ns": float(PI_PULSE_DURATION_S * 1.0e9),
            "echo_total_duration_ns": float(gaussian_total_duration * 1.0e9),
            "refocusing_convention": "vacuum_gaussian_pi",
            "half_corrections_vector": half_corrections.tolist(),
        },
    )
    save_extension_artifact(
        request,
        spec=spec,
        construction="echo_finite_gaussian",
        target_operator=target_operator,
        actual_operator=gaussian_restricted,
        tone_specs=half_optimization.tone_specs,
        corrections_vector=half_corrections,
        waveform_samples=gaussian_samples,
        metadata={
            "extension_scope": "echo_followup",
            "refocusing_convention": "vacuum_gaussian_pi",
            "half_corrections_vector": half_corrections.tolist(),
            "gaussian_pi_duration_s": float(PI_PULSE_DURATION_S),
        },
    )
    save_waveform_npz(WAVEFORM_DIR / f"extension_{request.case_id}_echo_finite_gaussian.npz", gaussian_samples)
    echo_rows.append(gaussian_row)

    refocus_spec = equal_angle_aligned_x_target(
        theta_rad=np.pi,
        n_active=request.n_active,
        family=f"{request.family}_refocus_multitone_xpi",
        description="Manifold-aware multitone square refocusing target for the echo follow-up.",
    )
    refocus_run_config = make_run_config(model, n_active=request.n_active, duration_s=MANIFOLD_AWARE_PI_DURATION_S)
    refocus_optimization = optimize_square_multitone(
        model,
        refocus_spec,
        refocus_run_config,
        n_starts=REFOCUS_N_STARTS,
        maxiter=REFOCUS_MAXITER,
        random_seed=26001 + int(round(1000.0 * request.chi_t_over_2pi)),
        label_prefix=f"extension_{request.case_id}_refocus",
    )
    refocus_corrections = corrections_to_vector(refocus_optimization.corrections, n_active=request.n_active)
    refocus_pulse = refocus_optimization.waveform.pulse
    manifold_pulses = [
        shift_pulse(half_pulse, t0=0.0, label=f"extension_{request.case_id}_half_1"),
        shift_pulse(refocus_pulse, t0=0.5 * duration_s, label=f"extension_{request.case_id}_refocus_1"),
        shift_pulse(half_pulse, t0=0.5 * duration_s + MANIFOLD_AWARE_PI_DURATION_S, label=f"extension_{request.case_id}_half_2"),
        shift_pulse(refocus_pulse, t0=duration_s + MANIFOLD_AWARE_PI_DURATION_S, label=f"extension_{request.case_id}_refocus_2"),
    ]
    manifold_total_duration = duration_s + 2.0 * MANIFOLD_AWARE_PI_DURATION_S
    manifold_compiled = compile_pulse_sequence(
        manifold_pulses,
        dt_s=DEFAULT_DT,
        total_duration_s=manifold_total_duration,
    )
    manifold_full = simulate_full_operator_on_logical_inputs(
        model,
        manifold_compiled,
        frame=frame,
        drive_ops={"qubit": "qubit"},
        levels=levels,
    )
    manifold_restricted = restricted_operator_from_full(manifold_full, model, levels)
    manifold_samples = channel_waveform_samples(manifold_compiled)
    manifold_row = row_from_operator(
        request,
        construction="echo_finite_manifold_aware_multitone",
        target_operator=target_operator,
        actual_operator=manifold_restricted,
        extra={
            "extension_scope": "echo_followup",
            "echo_half_duration_ns": float(0.5 * duration_s * 1.0e9),
            "echo_pi_duration_ns": float(MANIFOLD_AWARE_PI_DURATION_S * 1.0e9),
            "echo_total_duration_ns": float(manifold_total_duration * 1.0e9),
            "refocusing_convention": "manifold_aware_multitone_square",
            "half_corrections_vector": half_corrections.tolist(),
            "refocus_corrections_vector": refocus_corrections.tolist(),
            "refocus_process_fidelity": float(refocus_optimization.validation.restricted_process_fidelity),
            "refocus_average_gate_fidelity": float(
                average_gate_fidelity(
                    build_target_operator(refocus_spec, levels),
                    np.asarray(refocus_optimization.validation.restricted_operator, dtype=np.complex128),
                )
            ),
        },
    )
    save_extension_artifact(
        request,
        spec=spec,
        construction="echo_finite_manifold_aware_multitone",
        target_operator=target_operator,
        actual_operator=manifold_restricted,
        tone_specs=refocus_optimization.tone_specs,
        corrections_vector=refocus_corrections,
        waveform_samples=manifold_samples,
        metadata={
            "extension_scope": "echo_followup",
            "refocusing_convention": "manifold_aware_multitone_square",
            "half_corrections_vector": half_corrections.tolist(),
            "refocus_corrections_vector": refocus_corrections.tolist(),
            "refocus_duration_s": float(MANIFOLD_AWARE_PI_DURATION_S),
        },
    )
    save_waveform_npz(
        WAVEFORM_DIR / f"extension_{request.case_id}_echo_finite_manifold_aware_multitone.npz",
        manifold_samples,
    )
    echo_rows.append(manifold_row)

    return echo_rows


def figure_checkpoint_comparison(checkpoint_rows: list[dict[str, Any]]) -> None:
    apply_plot_style()
    df = pd.DataFrame(checkpoint_rows)
    full_df = df[df["construction"] == "full_shared_line"].copy().sort_values("chi_t_over_2pi")
    reduced_df = df[df["construction"] == "blockwise_exact_reduced"].copy().sort_values("chi_t_over_2pi")
    if full_df.empty:
        return

    case_labels = []
    for family in full_df["family"]:
        family_name = str(family)
        if "chiprime" in family_name:
            case_labels.append("tuned + chi'")
        elif "off_minus" in family_name:
            case_labels.append("off -")
        elif "off_plus" in family_name:
            case_labels.append("off +")
        else:
            case_labels.append("tuned")
    x = np.arange(len(full_df), dtype=float)
    width = 0.36
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))

    axes[0].bar(
        x - 0.5 * width,
        full_df["restricted_average_gate_fidelity"],
        width=width,
        color="#4477AA",
        label="full shared line",
    )
    axes[0].bar(
        x + 0.5 * width,
        reduced_df["restricted_average_gate_fidelity"],
        width=width,
        color="#CCBB44",
        label="exact reduced replay",
    )
    axes[0].set_ylabel("Restricted average gate fidelity")
    axes[0].set_ylim(0.0, 1.02)
    axes[0].set_xticks(x, case_labels)
    axes[0].grid(True, axis="y", alpha=0.25)
    axes[0].legend(frameon=False, loc="upper right")

    axes[1].bar(
        x,
        full_df["max_residual_z_error_rad"],
        color="#EE6677",
        width=0.55,
    )
    axes[1].set_ylabel("Max residual-Z error [rad]")
    axes[1].set_xticks(x, case_labels)
    axes[1].grid(True, axis="y", alpha=0.25)
    axes[1].set_title("Exact checkpoint diagnostics")

    fig.tight_layout()
    fig.savefig(FIGURES_DIR / f"{CHECKPOINT_FIGURE_BASENAME}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / f"{CHECKPOINT_FIGURE_BASENAME}.pdf", bbox_inches="tight")
    plt.close(fig)


def figure_echo_followup(rows: list[dict[str, Any]]) -> None:
    apply_plot_style()
    df = pd.DataFrame(rows)
    if df.empty:
        return
    order = [
        "full_shared_line",
        "echo_ideal_instantaneous",
        "echo_finite_gaussian",
        "echo_finite_manifold_aware_multitone",
    ]
    labels = ["plain", "echo ideal", "echo gaussian", "echo multitone"]
    palette = ["#4477AA", "#CCBB44", "#EE6677", "#228833"]
    families = list(dict.fromkeys(df["family"].tolist()))
    fig, axes = plt.subplots(2, len(families), figsize=(5.0 * len(families), 7.2), sharey="row")
    if len(families) == 1:
        axes = np.asarray(axes).reshape(2, 1)
    for column, family in enumerate(families):
        family_df = df[df["family"] == family].set_index("construction")
        fidelities = [
            float(family_df.loc[key, "restricted_average_gate_fidelity"]) if key in family_df.index else np.nan
            for key in order
        ]
        residuals = [
            float(family_df.loc[key, "max_residual_z_error_rad"]) if key in family_df.index else np.nan
            for key in order
        ]
        axes[0, column].bar(labels, fidelities, color=palette)
        axes[1, column].bar(labels, residuals, color=palette)
        axes[0, column].set_title(family.replace("equal_angle_aligned_x_", "").replace("_", " "))
        axes[0, column].set_ylim(0.0, 1.02)
        axes[0, column].grid(True, axis="y", alpha=0.25)
        axes[1, column].grid(True, axis="y", alpha=0.25)
    axes[0, 0].set_ylabel("Restricted average gate fidelity")
    axes[1, 0].set_ylabel("Max residual-Z error [rad]")
    for column in range(len(families)):
        axes[1, column].tick_params(axis="x", rotation=15)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / f"{ECHO_FIGURE_BASENAME}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / f"{ECHO_FIGURE_BASENAME}.pdf", bbox_inches="tight")
    plt.close(fig)


def summary_payload(
    *,
    tuned_map: dict[str, Any],
    checkpoint_rows: list[dict[str, Any]],
    echo_rows: list[dict[str, Any]],
    runtime_s: float,
) -> dict[str, Any]:
    checkpoint_df = pd.DataFrame(checkpoint_rows)
    echo_df = pd.DataFrame(echo_rows)
    tuned_full = checkpoint_df[
        (checkpoint_df["construction"] == "full_shared_line")
        & (checkpoint_df["family"] == "equal_angle_aligned_x_tuned")
    ]
    off_full = checkpoint_df[
        (checkpoint_df["construction"] == "full_shared_line")
        & (checkpoint_df["family"].isin(["equal_angle_aligned_x_off_minus", "equal_angle_aligned_x_off_plus"]))
    ]
    tuned_echo = echo_df[echo_df["family"] == "equal_angle_aligned_x_tuned"].copy()
    return {
        "runtime_s": float(runtime_s),
        "selected_tuned_root": tuned_map["equal_amplitude_roots"][0],
        "tuned_full_shared_line": None if tuned_full.empty else tuned_full.iloc[0].to_dict(),
        "off_tuned_mean_fidelity": float(off_full["restricted_average_gate_fidelity"].mean()) if not off_full.empty else float("nan"),
        "off_tuned_mean_max_residual_z_error_rad": float(off_full["max_residual_z_error_rad"].mean()) if not off_full.empty else float("nan"),
        "best_echo_construction_tuned": None
        if tuned_echo.empty
        else tuned_echo.sort_values("restricted_average_gate_fidelity", ascending=False).iloc[0].to_dict(),
        "echo_rows_count": int(len(echo_rows)),
        "checkpoint_rows_count": int(len(checkpoint_rows)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-starts", type=int, default=CHECKPOINT_N_STARTS)
    parser.add_argument("--maxiter", type=int, default=CHECKPOINT_MAXITER)
    args = parser.parse_args()

    start_time = time.perf_counter()
    CASE_DIR.mkdir(parents=True, exist_ok=True)
    WAVEFORM_DIR.mkdir(parents=True, exist_ok=True)

    tuned_map, x_values, delta_values, residual, roots = build_tuned_map()
    figure_tuned_map(x_values, delta_values, residual, roots)
    save_json(
        TUNED_MAP_PATH,
        tuned_map,
        description="Extension-pass tuned-set map for the strict two-block aligned-x cancellation conditions.",
    )

    tuned_root_chi_t = float(tuned_map["equal_amplitude_roots"][0]["chi_t_over_2pi"])
    checkpoint_bundles = [
        run_checkpoint_case(request, n_starts=int(args.n_starts), maxiter=int(args.maxiter))
        for request in extension_requests(tuned_root_chi_t)
    ]
    checkpoint_rows = [bundle["full_row"] for bundle in checkpoint_bundles] + [bundle["reduced_row"] for bundle in checkpoint_bundles]
    figure_checkpoint_comparison(checkpoint_rows)

    echo_source_families = {"equal_angle_aligned_x_tuned", "equal_angle_aligned_x_off_plus"}
    echo_rows: list[dict[str, Any]] = []
    for bundle in checkpoint_bundles:
        echo_rows.append(bundle["full_row"])
        if bundle["request"].family in echo_source_families:
            echo_rows.extend(run_echo_followup(bundle))
    figure_echo_followup(echo_rows)

    save_json(
        CHECKPOINT_PATH,
        {"rows": checkpoint_rows},
        description="Exact shared-line and reduced-replay checkpoint results for the extension pass.",
    )
    save_json(
        ECHO_PATH,
        {"rows": echo_rows},
        description="Echo robustness follow-up results for the aligned-x extension cases.",
    )

    summary = summary_payload(
        tuned_map=tuned_map,
        checkpoint_rows=checkpoint_rows,
        echo_rows=echo_rows,
        runtime_s=float(time.perf_counter() - start_time),
    )
    save_json(SUMMARY_PATH, summary, description="Headline summary for the extension pass results.")
    save_json(
        RESULTS_PATH,
        {
            "tuned_map_path": str(TUNED_MAP_PATH),
            "checkpoint_path": str(CHECKPOINT_PATH),
            "echo_path": str(ECHO_PATH),
            "summary_path": str(SUMMARY_PATH),
            "figures": [
                str(FIGURES_DIR / f"{TUNED_MAP_FIGURE_BASENAME}.pdf"),
                str(FIGURES_DIR / f"{CHECKPOINT_FIGURE_BASENAME}.pdf"),
                str(FIGURES_DIR / f"{ECHO_FIGURE_BASENAME}.pdf"),
            ],
            "runtime_s": float(time.perf_counter() - start_time),
        },
        description="Top-level index for the extension pass outputs.",
    )


if __name__ == "__main__":
    main()