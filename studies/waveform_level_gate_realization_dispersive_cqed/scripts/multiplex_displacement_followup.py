"""Multiplex-drive follow-up for unconditional displacement in dispersive cQED.

This script benchmarks two physically explicit approximations to the best
existing hardware-aware optimal-control waveform:

1. A full-duration multicarrier family built from dominant Fourier tones.
2. A segmented branch-resonant family that keeps only the two branch carriers
    but allows their complex weights to vary across time segments.
3. A low-parameter shaped two-tone family that keeps the calibrated branch-tone
    ratio and jointly optimizes one complex scale per segment.

The goal is to answer the concrete follow-up left open by the main study: can a
more structured and still interpretable multiplex drive recover the broad-state
advantage of the constrained waveform without falling all the way back to full
sampled optimal control?
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

import common
import unconditional_displacement_study as uds
from common import ARTIFACTS_DIR, FIGURES_DIR, TOL_BRIGHT, apply_plot_style, load_json, save_json


apply_plot_style()

ALPHA_TARGET = 1.0
TONE_COUNTS = (2, 3, 4, 5, 6, 8)
SEGMENT_COUNTS = (1, 2, 4, 8)
JOINT_SHAPED_SEGMENT_COUNT = 4
REFERENCE_TWO_TONE_DURATIONS_NS = (20.0, 40.0, 80.0)
OUTPUT_ARTIFACT = ARTIFACTS_DIR / "unconditional_multiplex_followup.json"
OUTPUT_FIGURE_STEM = FIGURES_DIR / "unconditional_multiplex_followup"


def best_optimal_case() -> dict[str, Any]:
    payload = load_json(ARTIFACTS_DIR / "unconditional_optimal_control_summary.json")
    return max(payload["cases"], key=lambda item: item["full_metrics"]["state_test_mean_fidelity"])


def optimal_complex_samples(case: dict[str, Any]) -> tuple[np.ndarray, float]:
    physical = np.asarray(case["physical_values"], dtype=float)
    samples = physical[0] + 1j * physical[1]
    duration_s = float(case["duration_ns"]) * 1.0e-9
    return np.asarray(samples, dtype=np.complex128), duration_s


def dominant_frequency_subset(samples: np.ndarray, duration_s: float, n_tones: int) -> tuple[np.ndarray, np.ndarray]:
    n_samples = int(samples.size)
    dt_s = duration_s / max(n_samples, 1)
    fft_coeffs = np.fft.fft(samples) / max(n_samples, 1)
    freqs = np.fft.fftfreq(n_samples, d=dt_s) * (2.0 * np.pi)
    order = np.argsort(np.abs(fft_coeffs))[::-1]
    selected = order[: int(n_tones)]
    return np.asarray(freqs[selected], dtype=float), np.asarray(fft_coeffs[selected], dtype=np.complex128)


def midpoint_fit_weights(samples: np.ndarray, duration_s: float, freqs_rad_s: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n_samples = int(samples.size)
    dt_s = duration_s / max(n_samples, 1)
    midpoints_s = (np.arange(n_samples, dtype=float) + 0.5) * dt_s
    design = np.exp(1j * np.outer(midpoints_s, freqs_rad_s))
    weights, *_ = np.linalg.lstsq(design, samples, rcond=None)
    fit = design @ weights
    return np.asarray(weights, dtype=np.complex128), np.asarray(fit, dtype=np.complex128)


def branch_fit_frequencies(model, frame) -> np.ndarray:
    return np.array(
        [
            common.cavity_branch_transition_frequency(model, frame, qubit_level=0),
            common.cavity_branch_transition_frequency(model, frame, qubit_level=1),
        ],
        dtype=float,
    )


def segmented_branch_fit(
    samples: np.ndarray,
    duration_s: float,
    fit_freqs_rad_s: np.ndarray,
    segment_count: int,
) -> tuple[list[Any], np.ndarray, np.ndarray]:
    n_samples = int(samples.size)
    if n_samples % int(segment_count) != 0:
        raise ValueError("Segment count must divide the number of optimal-control samples exactly.")
    dt_s = duration_s / max(n_samples, 1)
    midpoints_s = (np.arange(n_samples, dtype=float) + 0.5) * dt_s
    samples_per_segment = n_samples // int(segment_count)
    pulses: list[Any] = []
    weights = np.zeros((int(segment_count), fit_freqs_rad_s.size), dtype=np.complex128)
    fitted = np.zeros(n_samples, dtype=np.complex128)
    for segment_index in range(int(segment_count)):
        start = segment_index * samples_per_segment
        stop = start + samples_per_segment
        seg_times = midpoints_s[start:stop]
        seg_target = samples[start:stop]
        design = np.exp(1j * np.outer(seg_times, fit_freqs_rad_s))
        seg_weights, *_ = np.linalg.lstsq(design, seg_target, rcond=None)
        weights[segment_index, :] = np.asarray(seg_weights, dtype=np.complex128)
        fitted[start:stop] = design @ seg_weights
        segment_t0_s = float(start) * dt_s
        segment_duration_s = float(samples_per_segment) * dt_s
        for tone_index, (fit_freq_rad_s, weight) in enumerate(zip(fit_freqs_rad_s, seg_weights, strict=True)):
            if abs(weight) <= 1.0e-14:
                continue
            pulses.append(
                common.Pulse(
                    channel="storage",
                    t0=segment_t0_s,
                    duration=segment_duration_s,
                    envelope=common.square_envelope,
                    carrier=float(-fit_freq_rad_s),
                    phase=float(np.angle(weight)),
                    amp=float(abs(weight)),
                    drag=0.0,
                    sample_rate=None,
                    label=f"segment_{segment_count}_{segment_index}_{tone_index}",
                )
            )
    return pulses, weights, fitted


def pulses_from_common_scales(
    duration_s: float,
    fit_freqs_rad_s: np.ndarray,
    base_weights: np.ndarray,
    scales: np.ndarray,
) -> list[Any]:
    segment_count = int(scales.size)
    segment_duration_s = float(duration_s) / max(segment_count, 1)
    pulses: list[Any] = []
    for segment_index, seg_scale in enumerate(scales):
        segment_t0_s = float(segment_index) * segment_duration_s
        segment_weights = np.asarray(seg_scale, dtype=np.complex128) * np.asarray(base_weights, dtype=np.complex128)
        for tone_index, (fit_freq_rad_s, weight) in enumerate(zip(fit_freqs_rad_s, segment_weights, strict=True)):
            if abs(weight) <= 1.0e-14:
                continue
            pulses.append(
                common.Pulse(
                    channel="storage",
                    t0=segment_t0_s,
                    duration=segment_duration_s,
                    envelope=common.square_envelope,
                    carrier=float(-fit_freq_rad_s),
                    phase=float(np.angle(weight)),
                    amp=float(abs(weight)),
                    drag=0.0,
                    sample_rate=None,
                    label=f"shaped_two_tone_{segment_count}_{segment_index}_{tone_index}",
                )
            )
    return pulses


def common_scale_segment_fit(
    samples: np.ndarray,
    duration_s: float,
    fit_freqs_rad_s: np.ndarray,
    base_weights: np.ndarray,
    segment_count: int,
) -> tuple[list[Any], np.ndarray, np.ndarray]:
    n_samples = int(samples.size)
    if n_samples % int(segment_count) != 0:
        raise ValueError("Segment count must divide the number of optimal-control samples exactly.")
    dt_s = duration_s / max(n_samples, 1)
    midpoints_s = (np.arange(n_samples, dtype=float) + 0.5) * dt_s
    samples_per_segment = n_samples // int(segment_count)
    scales = np.zeros(int(segment_count), dtype=np.complex128)
    fitted = np.zeros(n_samples, dtype=np.complex128)
    base_waveform = np.exp(1j * np.outer(midpoints_s, fit_freqs_rad_s)) @ base_weights
    for segment_index in range(int(segment_count)):
        start = segment_index * samples_per_segment
        stop = start + samples_per_segment
        seg_basis = base_waveform[start:stop]
        seg_target = samples[start:stop]
        denom = np.vdot(seg_basis, seg_basis)
        seg_scale = 0.0j if abs(denom) <= 1.0e-15 else np.vdot(seg_basis, seg_target) / denom
        scales[segment_index] = np.complex128(seg_scale)
        fitted[start:stop] = seg_scale * seg_basis
    pulses = pulses_from_common_scales(duration_s, fit_freqs_rad_s, base_weights, scales)
    return pulses, scales, fitted


def multiplex_pulses(duration_s: float, freqs_rad_s: np.ndarray, weights: np.ndarray) -> list[Any]:
    pulses: list[Any] = []
    for tone_index, (freq_rad_s, weight) in enumerate(zip(freqs_rad_s, weights, strict=True)):
        if abs(weight) <= 1.0e-14:
            continue
        pulses.append(
            common.Pulse(
                channel="storage",
                t0=0.0,
                duration=duration_s,
                envelope=common.square_envelope,
                carrier=float(-freq_rad_s),
                phase=float(np.angle(weight)),
                amp=float(abs(weight)),
                drag=0.0,
                sample_rate=None,
                label=f"multiplex_{tone_index}",
            )
        )
    return pulses


def evaluate_multiplex_case(case: dict[str, Any], n_tones: int) -> dict[str, Any]:
    samples, duration_s = optimal_complex_samples(case)
    freqs_rad_s, fft_subset = dominant_frequency_subset(samples, duration_s, n_tones)
    weights, fitted_samples = midpoint_fit_weights(samples, duration_s, freqs_rad_s)
    pulses = multiplex_pulses(duration_s, freqs_rad_s, weights)

    model = uds.variant_model("full")
    frame = common.build_frame(model)
    session = common.compile_and_prepare(model, frame, pulses)
    vacuum_metrics = uds.branch_vacuum_metrics(model, session, alpha_target=ALPHA_TARGET)
    full_metrics = uds.protocol_metrics_from_session(
        model,
        session,
        alpha_target=ALPHA_TARGET,
        state_pairs=uds.FULL_STATE_TEST_SET,
    )
    fit_error = samples - fitted_samples
    return {
        "n_tones": int(n_tones),
        "duration_ns": float(case["duration_ns"]),
        "selected_detunings_rad_s": freqs_rad_s,
        "selected_detunings_mhz": freqs_rad_s / (2.0 * np.pi * 1.0e6),
        "raw_carriers_rad_s": -freqs_rad_s,
        "fft_seed_coeffs": fft_subset,
        "weights": weights,
        "sample_fit": {
            "target": samples,
            "fitted": fitted_samples,
            "relative_l2_error": float(np.linalg.norm(fit_error) / max(np.linalg.norm(samples), 1.0e-15)),
        },
        "vacuum_metrics": vacuum_metrics,
        "full_metrics": full_metrics,
    }


def evaluate_segmented_case(case: dict[str, Any], segment_count: int) -> dict[str, Any]:
    samples, duration_s = optimal_complex_samples(case)
    model = uds.variant_model("full")
    frame = common.build_frame(model)
    fit_freqs = branch_fit_frequencies(model, frame)
    pulses, weights, fitted_samples = segmented_branch_fit(samples, duration_s, fit_freqs, segment_count)
    session = common.compile_and_prepare(model, frame, pulses)
    vacuum_metrics = uds.branch_vacuum_metrics(model, session, alpha_target=ALPHA_TARGET)
    full_metrics = uds.protocol_metrics_from_session(
        model,
        session,
        alpha_target=ALPHA_TARGET,
        state_pairs=uds.FULL_STATE_TEST_SET,
    )
    fit_error = samples - fitted_samples
    return {
        "segment_count": int(segment_count),
        "duration_ns": float(case["duration_ns"]),
        "n_tones_per_segment": int(fit_freqs.size),
        "segment_detunings_rad_s": fit_freqs,
        "segment_detunings_mhz": fit_freqs / (2.0 * np.pi * 1.0e6),
        "raw_carriers_rad_s": -fit_freqs,
        "weights": weights,
        "sample_fit": {
            "target": samples,
            "fitted": fitted_samples,
            "relative_l2_error": float(np.linalg.norm(fit_error) / max(np.linalg.norm(samples), 1.0e-15)),
        },
        "vacuum_metrics": vacuum_metrics,
        "full_metrics": full_metrics,
    }


def evaluate_two_tone_reference(duration_ns: float) -> dict[str, Any]:
    model = uds.variant_model("full")
    frame = common.build_frame(model)
    calibration = uds.calibrate_two_tone(
        model=model,
        frame=frame,
        alpha_target=ALPHA_TARGET,
        duration_s=float(duration_ns) * 1.0e-9,
    )
    session = common.compile_and_prepare(model, frame, calibration["pulses"])
    return {
        "duration_ns": float(duration_ns),
        "vacuum_metrics": uds.branch_vacuum_metrics(model, session, alpha_target=ALPHA_TARGET),
        "full_metrics": uds.protocol_metrics_from_session(
            model,
            session,
            alpha_target=ALPHA_TARGET,
            state_pairs=uds.FULL_STATE_TEST_SET,
        ),
        "pulse_meta": calibration["meta"],
    }


def evaluate_joint_shaped_two_tone_case(case: dict[str, Any], segment_count: int) -> dict[str, Any]:
    samples, duration_s = optimal_complex_samples(case)
    model = uds.variant_model("full")
    frame = common.build_frame(model)
    calibration = uds.calibrate_two_tone(
        model=model,
        frame=frame,
        alpha_target=ALPHA_TARGET,
        duration_s=duration_s,
    )
    base_weights = np.asarray(calibration["meta"]["weights"], dtype=np.complex128)
    fit_freqs = -np.asarray(calibration["meta"]["carrier_detunings"], dtype=float)
    initial_pulses, initial_scales, fitted_samples = common_scale_segment_fit(
        samples,
        duration_s,
        fit_freqs,
        base_weights,
        segment_count,
    )
    initial_session = common.compile_and_prepare(model, frame, initial_pulses)
    initial_vacuum = uds.branch_vacuum_metrics(model, initial_session, alpha_target=ALPHA_TARGET)
    initial_full = uds.protocol_metrics_from_session(
        model,
        initial_session,
        alpha_target=ALPHA_TARGET,
        state_pairs=uds.FULL_STATE_TEST_SET,
    )

    def score_from_metrics(vacuum_metrics: dict[str, Any], full_metrics: dict[str, Any]) -> float:
        return (
            float(full_metrics["state_test_mean_fidelity"])
            + 0.25 * float(full_metrics["state_test_min_fidelity"])
            + 0.05 * float(vacuum_metrics["plus_x_fidelity"])
        )

    def unpack_scales(params: np.ndarray) -> np.ndarray:
        return np.asarray(params[0::2] + 1j * params[1::2], dtype=np.complex128)

    def objective(params: np.ndarray) -> float:
        scales = unpack_scales(params)
        pulses = pulses_from_common_scales(duration_s, fit_freqs, base_weights, scales)
        session = common.compile_and_prepare(model, frame, pulses)
        vacuum_metrics = uds.branch_vacuum_metrics(model, session, alpha_target=ALPHA_TARGET)
        full_metrics = uds.protocol_metrics_from_session(
            model,
            session,
            alpha_target=ALPHA_TARGET,
            state_pairs=uds.FULL_STATE_TEST_SET,
        )
        scale_penalty = 0.01 * max(0.0, float(np.max(np.abs(scales))) - 3.0) ** 2
        return -(score_from_metrics(vacuum_metrics, full_metrics)) + scale_penalty

    x0 = np.empty(2 * int(segment_count), dtype=float)
    x0[0::2] = np.real(initial_scales)
    x0[1::2] = np.imag(initial_scales)
    result = minimize(
        objective,
        x0,
        method="Powell",
        options={"maxfev": 600, "maxiter": 60, "xtol": 1.0e-3, "ftol": 1.0e-4},
    )
    optimized_scales = unpack_scales(np.asarray(result.x, dtype=float))
    optimized_pulses = pulses_from_common_scales(duration_s, fit_freqs, base_weights, optimized_scales)
    optimized_session = common.compile_and_prepare(model, frame, optimized_pulses)
    optimized_vacuum = uds.branch_vacuum_metrics(model, optimized_session, alpha_target=ALPHA_TARGET)
    optimized_full = uds.protocol_metrics_from_session(
        model,
        optimized_session,
        alpha_target=ALPHA_TARGET,
        state_pairs=uds.FULL_STATE_TEST_SET,
    )
    fit_error = samples - fitted_samples
    return {
        "segment_count": int(segment_count),
        "duration_ns": float(case["duration_ns"]),
        "base_weights": base_weights,
        "fit_frequencies_rad_s": fit_freqs,
        "fit_frequencies_mhz": fit_freqs / (2.0 * np.pi * 1.0e6),
        "initial_fit": {
            "scales": initial_scales,
            "sample_fit": {
                "target": samples,
                "fitted": fitted_samples,
                "relative_l2_error": float(np.linalg.norm(fit_error) / max(np.linalg.norm(samples), 1.0e-15)),
            },
            "vacuum_metrics": initial_vacuum,
            "full_metrics": initial_full,
        },
        "optimized": {
            "scales": optimized_scales,
            "vacuum_metrics": optimized_vacuum,
            "full_metrics": optimized_full,
        },
        "optimization": {
            "success": bool(result.success),
            "status": int(result.status),
            "message": str(result.message),
            "nfev": int(result.nfev),
            "score_initial": float(score_from_metrics(initial_vacuum, initial_full)),
            "score_final": float(score_from_metrics(optimized_vacuum, optimized_full)),
        },
    }


def generate_figure(payload: dict[str, Any]) -> None:
    cases = sorted(payload["multiplex_cases"], key=lambda item: item["n_tones"])
    tone_counts = np.array([case["n_tones"] for case in cases], dtype=int)
    mean_fidelity = np.array([case["full_metrics"]["state_test_mean_fidelity"] for case in cases], dtype=float)
    min_fidelity = np.array([case["full_metrics"]["state_test_min_fidelity"] for case in cases], dtype=float)
    plus_x_fidelity = np.array([case["vacuum_metrics"]["plus_x_fidelity"] for case in cases], dtype=float)
    delta_alpha = np.array([case["vacuum_metrics"]["delta_alpha"] for case in cases], dtype=float)

    best_case = max(cases, key=lambda item: item["full_metrics"]["state_test_mean_fidelity"])
    best_k = int(best_case["n_tones"])
    segmented_cases = sorted(payload["segmented_cases"], key=lambda item: item["segment_count"])
    segment_counts = np.array([case["segment_count"] for case in segmented_cases], dtype=int)
    segmented_mean_fidelity = np.array(
        [case["full_metrics"]["state_test_mean_fidelity"] for case in segmented_cases],
        dtype=float,
    )
    segmented_min_fidelity = np.array(
        [case["full_metrics"]["state_test_min_fidelity"] for case in segmented_cases],
        dtype=float,
    )
    best_segmented = max(
        segmented_cases,
        key=lambda item: item["full_metrics"]["state_test_mean_fidelity"],
    )
    joint_case = payload["joint_shaped_two_tone_case"]
    joint_metrics = joint_case["optimized"]["full_metrics"]
    best_segment_count = int(best_segmented["segment_count"])
    best_fit = best_segmented["sample_fit"]
    best_target = np.asarray(best_fit["target"], dtype=np.complex128)
    best_fitted = np.asarray(best_fit["fitted"], dtype=np.complex128)
    duration_s = float(best_case["duration_ns"]) * 1.0e-9
    n_samples = int(best_target.size)
    time_ns = (np.arange(n_samples, dtype=float) + 0.5) * duration_s * 1.0e9 / max(n_samples, 1)

    two_tone_20 = payload["two_tone_references"]["20.0"]
    two_tone_40 = payload["two_tone_references"]["40.0"]
    optimal = payload["optimal_reference"]

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.0))

    axes[0, 0].plot(tone_counts, mean_fidelity, "o-", color=TOL_BRIGHT[0], label="Full-duration mean fidelity")
    axes[0, 0].plot(tone_counts, min_fidelity, "s--", color=TOL_BRIGHT[1], label="Full-duration min fidelity")
    axes[0, 0].axhline(two_tone_40["full_metrics"]["state_test_mean_fidelity"], color=TOL_BRIGHT[2], linestyle=":", label="Two-tone 40 ns mean")
    axes[0, 0].axhline(optimal["full_metrics"]["state_test_mean_fidelity"], color=TOL_BRIGHT[3], linestyle="-.", label="Optimal 40 ns mean")
    axes[0, 0].set_xlabel("Tone count")
    axes[0, 0].set_ylabel("Fidelity")
    axes[0, 0].set_title("Full-duration multiplex benchmark")
    axes[0, 0].legend(fontsize=8)

    axes[0, 1].plot(segment_counts, segmented_mean_fidelity, "o-", color=TOL_BRIGHT[4], label="Segmented mean fidelity")
    axes[0, 1].plot(segment_counts, segmented_min_fidelity, "d--", color=TOL_BRIGHT[5], label="Segmented min fidelity")
    axes[0, 1].axhline(two_tone_20["full_metrics"]["state_test_mean_fidelity"], color=TOL_BRIGHT[2], linestyle=":", label="Two-tone 20 ns mean")
    axes[0, 1].axhline(optimal["full_metrics"]["state_test_mean_fidelity"], color=TOL_BRIGHT[3], linestyle="-.", label="Optimal 40 ns mean")
    axes[0, 1].axhline(joint_metrics["state_test_mean_fidelity"], color=TOL_BRIGHT[0], linestyle="--", label="Joint-shaped 4-seg mean")
    axes[0, 1].set_xlabel("Segment count")
    axes[0, 1].set_ylabel("Fidelity")
    axes[0, 1].set_title("Segmented branch-resonant benchmark")
    axes[0, 1].legend(fontsize=8)

    axes[1, 0].plot(time_ns, np.real(best_target), "o-", color=TOL_BRIGHT[0], label="Optimal I")
    axes[1, 0].plot(time_ns, np.real(best_fitted), "--", color=TOL_BRIGHT[1], label=f"Segmented I ({best_segment_count} segments)")
    axes[1, 0].plot(time_ns, np.imag(best_target), "o-", color=TOL_BRIGHT[2], label="Optimal Q")
    axes[1, 0].plot(time_ns, np.imag(best_fitted), "--", color=TOL_BRIGHT[3], label=f"Segmented Q ({best_segment_count} segments)")
    axes[1, 0].set_xlabel("Time (ns)")
    axes[1, 0].set_ylabel(r"Drive coefficient (rad/s)")
    axes[1, 0].set_title("Best segmented fit to the optimal waveform")
    axes[1, 0].legend(fontsize=8, ncol=2)

    optimized_scales = np.asarray(joint_case["optimized"]["scales"], dtype=np.complex128)
    segment_axis = np.arange(optimized_scales.size, dtype=int) + 1
    axes[1, 1].plot(
        segment_axis,
        np.abs(optimized_scales),
        "o-",
        color=TOL_BRIGHT[0],
        label="|scale_s|",
    )
    axes[1, 1].plot(
        segment_axis,
        np.angle(optimized_scales),
        "s--",
        color=TOL_BRIGHT[1],
        label=r"arg(scale_s)",
    )
    axes[1, 1].set_xlabel("Segment index")
    axes[1, 1].set_ylabel("Optimized segment scale")
    axes[1, 1].set_title("Joint-shaped two-tone segment scales")
    axes[1, 1].legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(OUTPUT_FIGURE_STEM.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(OUTPUT_FIGURE_STEM.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    started = time.time()
    optimal = best_optimal_case()
    multiplex_cases = [evaluate_multiplex_case(optimal, n_tones) for n_tones in TONE_COUNTS]
    segmented_cases = [evaluate_segmented_case(optimal, segment_count) for segment_count in SEGMENT_COUNTS]
    joint_shaped_two_tone_case = evaluate_joint_shaped_two_tone_case(optimal, JOINT_SHAPED_SEGMENT_COUNT)
    two_tone_refs = {
        f"{duration_ns:.1f}": evaluate_two_tone_reference(duration_ns)
        for duration_ns in REFERENCE_TWO_TONE_DURATIONS_NS
    }
    best_multiplex = max(
        multiplex_cases,
        key=lambda item: item["full_metrics"]["state_test_mean_fidelity"],
    )
    best_segmented = max(
        segmented_cases,
        key=lambda item: item["full_metrics"]["state_test_mean_fidelity"],
    )
    payload = {
        "alpha_target": ALPHA_TARGET,
        "tone_counts": list(TONE_COUNTS),
        "segment_counts": list(SEGMENT_COUNTS),
        "optimal_reference": optimal,
        "two_tone_references": two_tone_refs,
        "multiplex_cases": multiplex_cases,
        "segmented_cases": segmented_cases,
        "joint_shaped_two_tone_case": joint_shaped_two_tone_case,
        "best_multiplex_case": best_multiplex,
        "best_segmented_case": best_segmented,
        "summary": {
            "best_multiplex_mean_fidelity": best_multiplex["full_metrics"]["state_test_mean_fidelity"],
            "best_multiplex_min_fidelity": best_multiplex["full_metrics"]["state_test_min_fidelity"],
            "best_multiplex_delta_alpha": best_multiplex["vacuum_metrics"]["delta_alpha"],
            "best_multiplex_plus_x_fidelity": best_multiplex["vacuum_metrics"]["plus_x_fidelity"],
            "best_segmented_mean_fidelity": best_segmented["full_metrics"]["state_test_mean_fidelity"],
            "best_segmented_min_fidelity": best_segmented["full_metrics"]["state_test_min_fidelity"],
            "best_segmented_delta_alpha": best_segmented["vacuum_metrics"]["delta_alpha"],
            "best_segmented_plus_x_fidelity": best_segmented["vacuum_metrics"]["plus_x_fidelity"],
            "best_segmented_segment_count": best_segmented["segment_count"],
            "joint_shaped_two_tone_mean_fidelity": joint_shaped_two_tone_case["optimized"]["full_metrics"]["state_test_mean_fidelity"],
            "joint_shaped_two_tone_min_fidelity": joint_shaped_two_tone_case["optimized"]["full_metrics"]["state_test_min_fidelity"],
            "joint_shaped_two_tone_plus_x_fidelity": joint_shaped_two_tone_case["optimized"]["vacuum_metrics"]["plus_x_fidelity"],
            "joint_shaped_two_tone_delta_alpha": joint_shaped_two_tone_case["optimized"]["vacuum_metrics"]["delta_alpha"],
            "two_tone_20_mean_fidelity": two_tone_refs["20.0"]["full_metrics"]["state_test_mean_fidelity"],
            "two_tone_40_mean_fidelity": two_tone_refs["40.0"]["full_metrics"]["state_test_mean_fidelity"],
            "optimal_mean_fidelity": optimal["full_metrics"]["state_test_mean_fidelity"],
            "optimal_min_fidelity": optimal["full_metrics"]["state_test_min_fidelity"],
        },
        "wall_time_s": time.time() - started,
    }
    generate_figure(payload)
    save_json(
        OUTPUT_ARTIFACT,
        payload,
        description="multiplex displacement follow-up benchmark",
        load_instructions="Use common.load_json(Path('unconditional_multiplex_followup.json')) to load the payload.",
    )
    best_case = payload["best_multiplex_case"]
    print("=" * 72)
    print("Multiplex displacement follow-up complete")
    print("=" * 72)
    print(
        f"Best multiplex case: {best_case['n_tones']} tones, mean fidelity="
        f"{best_case['full_metrics']['state_test_mean_fidelity']:.6f}, min fidelity="
        f"{best_case['full_metrics']['state_test_min_fidelity']:.6f}"
    )
    print(
        f"Best segmented case: {best_segmented['segment_count']} segments, mean fidelity="
        f"{best_segmented['full_metrics']['state_test_mean_fidelity']:.6f}, min fidelity="
        f"{best_segmented['full_metrics']['state_test_min_fidelity']:.6f}"
    )
    print(
        f"Joint-shaped two-tone case: {joint_shaped_two_tone_case['segment_count']} segments, mean fidelity="
        f"{joint_shaped_two_tone_case['optimized']['full_metrics']['state_test_mean_fidelity']:.6f}, min fidelity="
        f"{joint_shaped_two_tone_case['optimized']['full_metrics']['state_test_min_fidelity']:.6f}"
    )
    print(f"Total wall time: {payload['wall_time_s']:.1f} s")


if __name__ == "__main__":
    main()