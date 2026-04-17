from __future__ import annotations

from typing import Any, Callable, Sequence

import numpy as np
from scipy.optimize import Bounds, minimize

import runtime_compat  # noqa: F401

from cqed_sim.calibration.conditioned_multitone import ConditionedMultitoneCorrections

from analysis import (
    CaseContext,
    correction_bounds,
    quick_candidate_metrics,
    evaluate_candidate_full,
    evaluate_sequence_fast,
)
from common import (
    TWO_PI,
    build_multitone_waveform_from_corrections,
    corrections_from_vector,
    corrections_to_dict,
    corrections_to_vector,
    gaussian_samples,
    make_gaussian_qubit_rotation_pulse,
    make_run_config,
    normalize_complex_samples,
    orthogonal_basis,
    shift_pulse,
    simulate_full_operator_on_logical_inputs,
    single_pulse_sequence,
)


def direct_sequence_from_corrections(
    context: CaseContext,
    corrections: ConditionedMultitoneCorrections,
    *,
    label: str,
    duration_s: float | None = None,
) -> tuple[list[Any], dict[str, str], dict[str, Any]]:
    use_duration = float(context.duration_s if duration_s is None else duration_s)
    run_config = make_run_config(context.model, n_active=context.request.n_active, duration_s=use_duration, dt_s=float(context.run_config.dt_s))
    waveform, tone_specs = build_multitone_waveform_from_corrections(
        context.model,
        context.spec,
        run_config,
        corrections=corrections,
        label=label,
    )
    return [waveform.pulse], waveform.drive_ops, {
        "construction": "direct",
        "duration_s": float(use_duration),
        "active_duration_s": float(use_duration),
        "total_gate_duration_s": float(use_duration),
        "fairness_mode": "fixed_total_duration",
        "corrections": corrections_to_dict(corrections),
        "tone_rows": [
            {
                "manifold": int(tone.manifold),
                "omega_rad_s": float(tone.omega_rad_s),
                "omega_hz": float(tone.omega_rad_s / TWO_PI),
                "amp_rad_s": float(tone.amp_rad_s),
                "phase_rad": float(tone.phase_rad),
            }
            for tone in tone_specs
        ],
    }


def _single_segment_run_config(context: CaseContext, duration_s: float) -> Any:
    return make_run_config(context.model, n_active=context.request.n_active, duration_s=float(duration_s), dt_s=float(context.run_config.dt_s))


def sampled_family_bounds(family: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if family == "symmetric_two_segment":
        return (
            np.asarray([0.0, 0.0, 0.0, 0.0, 1.0, 0.0], dtype=float),
            np.asarray([-0.9, -0.9, -1.2, -1.2, 0.2, -np.pi], dtype=float),
            np.asarray([0.9, 0.9, 1.2, 1.2, 1.6, np.pi], dtype=float),
        )
    if family == "complex_envelope":
        return (
            np.zeros(6, dtype=float),
            np.asarray([-0.9, -0.9, -0.9, -1.4, -1.4, -1.4], dtype=float),
            np.asarray([0.9, 0.9, 0.9, 1.4, 1.4, 1.4], dtype=float),
        )
    if family == "basis_expanded":
        return (
            np.zeros(8, dtype=float),
            np.asarray([-0.8] * 4 + [-1.0] * 4, dtype=float),
            np.asarray([0.8] * 4 + [1.0] * 4, dtype=float),
        )
    raise ValueError(f"Unsupported sampled family '{family}'.")


def build_symmetric_two_segment_sequence(context: CaseContext, params: np.ndarray) -> tuple[list[Any], dict[str, str], dict[str, Any]]:
    half_duration = 0.5 * float(context.duration_s)
    count = max(16, int(round(half_duration / float(context.run_config.dt_s))))
    base = gaussian_samples(count, sigma_fraction=float(context.run_config.sigma_fraction))
    cos_basis = orthogonal_basis(count, 2, kind="cos")
    amp_profile = 1.0 + float(params[0]) * cos_basis[0] + float(params[1]) * cos_basis[1]
    phase_profile = float(params[2]) * cos_basis[0] + float(params[3]) * cos_basis[1]
    samples_1 = normalize_complex_samples(base * amp_profile * np.exp(1.0j * phase_profile))
    samples_2 = normalize_complex_samples(float(params[4]) * samples_1[::-1].conj() * np.exp(1.0j * float(params[5])))
    run_half = _single_segment_run_config(context, half_duration)
    tone_specs = context.direct_optimization.optimized_result.waveform.tone_specs
    pulses, drive_ops = single_pulse_sequence(tone_specs, run_half, base_samples=samples_1, label="two_segment_env_1")
    pulse_a = shift_pulse(pulses[0], t0=0.0, label="two_segment_env_1")
    pulses, drive_ops = single_pulse_sequence(tone_specs, run_half, base_samples=samples_2, label="two_segment_env_2")
    pulse_b = shift_pulse(pulses[0], t0=half_duration, label="two_segment_env_2")
    return [pulse_a, pulse_b], drive_ops, {
        "construction": "symmetric_two_segment",
        "parameters": [float(x) for x in np.asarray(params, dtype=float)],
        "segment_count": 2,
        "segment_duration_s": float(half_duration),
        "active_duration_s": float(context.duration_s),
        "total_gate_duration_s": float(context.duration_s),
        "fairness_mode": "fixed_total_duration",
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
    tone_specs = context.direct_optimization.optimized_result.waveform.tone_specs
    pulses, drive_ops = single_pulse_sequence(tone_specs, context.run_config, base_samples=samples, label="complex_envelope")
    return pulses, drive_ops, {
        "construction": "complex_envelope",
        "parameters": [float(x) for x in np.asarray(params, dtype=float)],
        "active_duration_s": float(context.duration_s),
        "total_gate_duration_s": float(context.duration_s),
        "fairness_mode": "fixed_total_duration",
    }


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
    tone_specs = context.direct_optimization.optimized_result.waveform.tone_specs
    pulses, drive_ops = single_pulse_sequence(tone_specs, context.run_config, base_samples=samples, label="basis_expanded")
    return pulses, drive_ops, {
        "construction": "basis_expanded",
        "parameters": [float(x) for x in np.asarray(params, dtype=float)],
        "active_duration_s": float(context.duration_s),
        "total_gate_duration_s": float(context.duration_s),
        "fairness_mode": "fixed_total_duration",
    }


def sampled_family_builder(family: str) -> Callable[[CaseContext, np.ndarray], tuple[list[Any], dict[str, str], dict[str, Any]]]:
    if family == "symmetric_two_segment":
        return build_symmetric_two_segment_sequence
    if family == "complex_envelope":
        return build_complex_envelope_sequence
    if family == "basis_expanded":
        return build_basis_expanded_sequence
    raise ValueError(f"Unsupported sampled family '{family}'.")


def build_echo_compiled(
    context: CaseContext,
    *,
    corrections_1: ConditionedMultitoneCorrections,
    corrections_2: ConditionedMultitoneCorrections,
    segment_1_duration_s: float,
    segment_2_duration_s: float,
    construction: str,
) -> tuple[list[Any], Any, dict[str, str], dict[str, Any]]:
    run_1 = _single_segment_run_config(context, segment_1_duration_s)
    run_2 = _single_segment_run_config(context, segment_2_duration_s)
    waveform_1, tones_1 = build_multitone_waveform_from_corrections(context.model, context.half_spec, run_1, corrections=corrections_1, label=f"{construction}_seg1")
    waveform_2, tones_2 = build_multitone_waveform_from_corrections(context.model, context.half_spec, run_2, corrections=corrections_2, label=f"{construction}_seg2")
    channel = str(waveform_1.pulse.channel)
    x_first = make_gaussian_qubit_rotation_pulse(
        context.model,
        context.frame,
        theta=np.pi,
        phase=0.0,
        duration_s=40.0e-9,
        channel=channel,
        manifold_level=0,
        sigma_fraction=0.25,
        t0=segment_1_duration_s,
        label=f"{construction}_xpi_1",
    )
    x_second = make_gaussian_qubit_rotation_pulse(
        context.model,
        context.frame,
        theta=np.pi,
        phase=0.0,
        duration_s=40.0e-9,
        channel=channel,
        manifold_level=0,
        sigma_fraction=0.25,
        t0=segment_1_duration_s + 40.0e-9 + segment_2_duration_s,
        label=f"{construction}_xpi_2",
    )
    pulses = [
        shift_pulse(waveform_1.pulse, t0=0.0, label=f"{construction}_seg1"),
        x_first,
        shift_pulse(waveform_2.pulse, t0=segment_1_duration_s + 40.0e-9, label=f"{construction}_seg2"),
        x_second,
    ]
    from common import compile_pulse_sequence

    total_gate_duration_s = float(segment_1_duration_s + segment_2_duration_s + 80.0e-9)
    compiled = compile_pulse_sequence(pulses, dt_s=float(context.run_config.dt_s), total_duration_s=total_gate_duration_s)
    return pulses, compiled, waveform_1.drive_ops, {
        "construction": str(construction),
        "segment_1_duration_s": float(segment_1_duration_s),
        "segment_2_duration_s": float(segment_2_duration_s),
        "active_duration_s": float(segment_1_duration_s + segment_2_duration_s),
        "total_gate_duration_s": float(total_gate_duration_s),
        "duration_asymmetry_eta": float((segment_1_duration_s - segment_2_duration_s) / (segment_1_duration_s + segment_2_duration_s)),
        "fairness_mode": "fixed_active_duration",
        "half_waveforms_identical": bool(
            np.allclose(corrections_to_vector(corrections_1), corrections_to_vector(corrections_2))
            and np.isclose(segment_1_duration_s, segment_2_duration_s)
        ),
        "corrections_segment_1": corrections_to_dict(corrections_1),
        "corrections_segment_2": corrections_to_dict(corrections_2),
        "tone_rows_segment_1": [
            {
                "manifold": int(tone.manifold),
                "omega_rad_s": float(tone.omega_rad_s),
                "omega_hz": float(tone.omega_rad_s / TWO_PI),
                "amp_rad_s": float(tone.amp_rad_s),
                "phase_rad": float(tone.phase_rad),
            }
            for tone in tones_1
        ],
        "tone_rows_segment_2": [
            {
                "manifold": int(tone.manifold),
                "omega_rad_s": float(tone.omega_rad_s),
                "omega_hz": float(tone.omega_rad_s / TWO_PI),
                "amp_rad_s": float(tone.amp_rad_s),
                "phase_rad": float(tone.phase_rad),
            }
            for tone in tones_2
        ],
    }


def optimize_reduced_direct(context: CaseContext) -> tuple[ConditionedMultitoneCorrections, dict[str, Any]]:
    lower, upper = correction_bounds(len(context.levels))
    x0 = corrections_to_vector(context.direct_optimization.optimized_corrections)
    history: list[dict[str, float]] = []

    def objective(vector: np.ndarray) -> float:
        corr = corrections_from_vector(np.asarray(vector, dtype=float), len(context.levels))
        pulses, drive_ops, _ = direct_sequence_from_corrections(context, corr, label="reduced_unitary_direct")
        metrics = evaluate_sequence_fast(context, pulses, drive_ops)
        value = float((1.0 - metrics["strict_reduced_mean"]) + 0.20 * (1.0 - metrics["strict_joint"]) + 5.0e-4 * np.mean(np.asarray(vector, dtype=float) ** 2))
        history.append({"objective": float(value), "strict_joint": float(metrics["strict_joint"]), "strict_reduced_mean": float(metrics["strict_reduced_mean"])})
        return value

    result = minimize(
        objective,
        x0,
        method="Powell",
        bounds=Bounds(lower, upper),
        options={"maxiter": 18 if context.request.stage == "screen" else 26},
    )
    corr = corrections_from_vector(np.asarray(result.x, dtype=float), len(context.levels))
    return corr, {
        "method": "Powell",
        "success": bool(result.success),
        "message": str(result.message),
        "history": history[-120:],
        "optimized_corrections": corrections_to_dict(corr),
    }


def optimize_sampled_direct_family(context: CaseContext, family: str) -> tuple[np.ndarray, dict[str, Any]]:
    x0, lower, upper = sampled_family_bounds(family)
    builder = sampled_family_builder(family)
    history: list[dict[str, float]] = []

    def objective(params: np.ndarray) -> float:
        pulses, drive_ops, _ = builder(context, np.asarray(params, dtype=float))
        metrics = evaluate_sequence_fast(context, pulses, drive_ops)
        value = float(metrics["strict_objective"] + 3.0e-4 * np.mean(np.asarray(params, dtype=float) ** 2))
        history.append({"objective": float(value), "strict_joint": float(metrics["strict_joint"]), "cpsqr_joint": float(metrics["cpsqr_joint"])})
        return value

    result = minimize(
        objective,
        x0,
        method="Powell",
        bounds=Bounds(lower, upper),
        options={"maxiter": 20 if context.request.stage == "screen" else 28},
    )
    return np.asarray(result.x, dtype=float), {
        "method": "Powell",
        "success": bool(result.success),
        "message": str(result.message),
        "parameters": [float(x) for x in np.asarray(result.x, dtype=float)],
        "history": history[-120:],
    }


def optimize_echo_variant(
    context: CaseContext,
    *,
    objective_mode: str,
    asymmetric: bool,
) -> tuple[ConditionedMultitoneCorrections, ConditionedMultitoneCorrections, float, dict[str, Any]]:
    seed = corrections_to_vector(context.half_optimization.optimized_corrections)
    lower, upper = correction_bounds(len(context.levels))
    x0 = np.concatenate([seed, seed, np.asarray([0.0], dtype=float) if asymmetric else np.asarray([], dtype=float)])
    lb = np.concatenate([lower, lower, np.asarray([-0.18], dtype=float) if asymmetric else np.asarray([], dtype=float)])
    ub = np.concatenate([upper, upper, np.asarray([0.18], dtype=float) if asymmetric else np.asarray([], dtype=float)])
    history: list[dict[str, float]] = []

    def objective(vector: np.ndarray) -> float:
        arr = np.asarray(vector, dtype=float)
        corr_1 = corrections_from_vector(arr[: seed.size], len(context.levels))
        corr_2 = corrections_from_vector(arr[seed.size : 2 * seed.size], len(context.levels))
        eta = float(arr[-1]) if asymmetric else 0.0
        seg_1 = 0.5 * float(context.duration_s) * (1.0 + eta)
        seg_2 = float(context.duration_s) - seg_1
        _, compiled, drive_ops, _ = build_echo_compiled(
            context,
            corrections_1=corr_1,
            corrections_2=corr_2,
            segment_1_duration_s=seg_1,
            segment_2_duration_s=seg_2,
            construction="echoed_asymmetric" if asymmetric else "echoed_variant",
        )
        full_operator = simulate_full_operator_on_logical_inputs(context.model, compiled, frame=context.frame, drive_ops=drive_ops, levels=context.levels)
        qm = quick_candidate_metrics(context, full_operator)
        base = qm["strict_objective"] if objective_mode == "strict" else qm["cpsqr_objective"]
        value = float(base + 3.0e-4 * np.mean(arr**2))
        history.append({"objective": float(value), "strict_joint": float(qm["strict_joint"]), "cpsqr_joint": float(qm["cpsqr_joint"])})
        return value

    result = minimize(
        objective,
        x0,
        method="Powell",
        bounds=Bounds(lb, ub),
        options={"maxiter": 18 if context.request.stage == "screen" else 26},
    )
    arr = np.asarray(result.x, dtype=float)
    corr_1 = corrections_from_vector(arr[: seed.size], len(context.levels))
    corr_2 = corrections_from_vector(arr[seed.size : 2 * seed.size], len(context.levels))
    eta = float(arr[-1]) if asymmetric else 0.0
    return corr_1, corr_2, eta, {
        "method": "Powell",
        "objective_mode": str(objective_mode),
        "success": bool(result.success),
        "message": str(result.message),
        "history": history[-120:],
        "optimized_corrections_segment_1": corrections_to_dict(corr_1),
        "optimized_corrections_segment_2": corrections_to_dict(corr_2),
        "duration_asymmetry_eta": float(eta),
    }


def run_gaussian_seed(context: CaseContext) -> tuple[dict[str, Any], dict[str, Any]]:
    corr = ConditionedMultitoneCorrections.zeros(len(context.levels))
    pulses, drive_ops, metadata = direct_sequence_from_corrections(context, corr, label="gaussian_seed")
    return evaluate_candidate_full(context, "gaussian_seed", pulses=pulses, drive_ops=drive_ops, metadata=metadata, optimizer_payload={"optimized_corrections": corrections_to_dict(corr)}, objective_mode="strict")


def run_native_direct_strict(context: CaseContext) -> tuple[dict[str, Any], dict[str, Any]]:
    pulses, drive_ops, metadata = direct_sequence_from_corrections(context, context.direct_optimization.optimized_corrections, label="native_direct_strict")
    return evaluate_candidate_full(context, "native_direct_strict", pulses=pulses, drive_ops=drive_ops, metadata=metadata, optimizer_payload={"optimized_corrections": corrections_to_dict(context.direct_optimization.optimized_corrections), "source": "targeted_subspace_optimizer"}, objective_mode="strict")


def run_reduced_unitary_direct(context: CaseContext) -> tuple[dict[str, Any], dict[str, Any]]:
    corr, payload = optimize_reduced_direct(context)
    pulses, drive_ops, metadata = direct_sequence_from_corrections(context, corr, label="reduced_unitary_direct")
    return evaluate_candidate_full(context, "reduced_unitary_direct", pulses=pulses, drive_ops=drive_ops, metadata=metadata, optimizer_payload=payload, objective_mode="strict")


def run_sampled_direct(context: CaseContext, family: str) -> tuple[dict[str, Any], dict[str, Any]]:
    params, payload = optimize_sampled_direct_family(context, family)
    pulses, drive_ops, metadata = sampled_family_builder(family)(context, params)
    return evaluate_candidate_full(context, family, pulses=pulses, drive_ops=drive_ops, metadata=metadata, optimizer_payload=payload, objective_mode="strict")


def run_echoed_symmetric(context: CaseContext) -> tuple[dict[str, Any], dict[str, Any]]:
    pulses, _, drive_ops, metadata = build_echo_compiled(
        context,
        corrections_1=context.half_optimization.optimized_corrections,
        corrections_2=context.half_optimization.optimized_corrections,
        segment_1_duration_s=0.5 * context.duration_s,
        segment_2_duration_s=0.5 * context.duration_s,
        construction="echoed_symmetric",
    )
    return evaluate_candidate_full(context, "echoed_symmetric", pulses=pulses, drive_ops=drive_ops, metadata=metadata, optimizer_payload={"optimized_corrections": corrections_to_dict(context.half_optimization.optimized_corrections), "source": "half_target_replay"}, objective_mode="strict")


def run_echo_variant(context: CaseContext, *, objective_mode: str, asymmetric: bool, family_name: str) -> tuple[dict[str, Any], dict[str, Any]]:
    corr_1, corr_2, eta, payload = optimize_echo_variant(context, objective_mode=objective_mode, asymmetric=asymmetric)
    seg_1 = 0.5 * float(context.duration_s) * (1.0 + eta)
    seg_2 = float(context.duration_s) - seg_1
    pulses, _, drive_ops, metadata = build_echo_compiled(
        context,
        corrections_1=corr_1,
        corrections_2=corr_2,
        segment_1_duration_s=seg_1,
        segment_2_duration_s=seg_2,
        construction=family_name,
    )
    return evaluate_candidate_full(context, family_name, pulses=pulses, drive_ops=drive_ops, metadata=metadata, optimizer_payload=payload, objective_mode=objective_mode)


def family_runner(family_name: str) -> Callable[[CaseContext], tuple[dict[str, Any], dict[str, Any]]]:
    if family_name == "gaussian_seed":
        return run_gaussian_seed
    if family_name == "native_direct_strict":
        return run_native_direct_strict
    if family_name == "reduced_unitary_direct":
        return run_reduced_unitary_direct
    if family_name in {"symmetric_two_segment", "complex_envelope", "basis_expanded"}:
        return lambda context: run_sampled_direct(context, family_name)
    if family_name == "echoed_symmetric":
        return run_echoed_symmetric
    if family_name == "echoed_independent":
        return lambda context: run_echo_variant(context, objective_mode="strict", asymmetric=False, family_name="echoed_independent")
    if family_name == "echoed_asymmetric":
        return lambda context: run_echo_variant(context, objective_mode="strict", asymmetric=True, family_name="echoed_asymmetric")
    if family_name == "echoed_cpsqr":
        return lambda context: run_echo_variant(context, objective_mode="cpsqr", asymmetric=False, family_name="echoed_cpsqr")
    raise ValueError(f"Unsupported family '{family_name}'.")
