"""Run the Fock-resolved black-box SQR inference study."""

from __future__ import annotations

import os
import time
from pathlib import Path
from typing import Any, Sequence

import matplotlib.pyplot as plt
import numpy as np

import study_lib as lib

from common import ARTIFACTS_DIR, DATA_DIR, FIGURES_DIR, STUDY_DIR, apply_plot_style, save_figure, save_json


RESULTS_JSON = DATA_DIR / "study_results.json"
SUMMARY_JSON = DATA_DIR / "study_summary.json"
SUMMARY_MD = DATA_DIR / "study_summary.md"
VALIDATION_HINTS_JSON = DATA_DIR / "validation_inputs.json"

PROFILES: dict[str, dict[str, Any]] = {
    "quick": {
        "baseline_repeats": 8,
        "benchmark_repeats": 6,
        "robustness_repeats": 5,
        "full_state_restarts": 4,
        "benchmark_shots": (100, 1_000, 10_000),
        "robustness_shots": 1_000,
        "coherence_shots": 4_000,
        "pulse_shots": 2_000,
    },
    "full": {
        "baseline_repeats": 24,
        "benchmark_repeats": 16,
        "robustness_repeats": 12,
        "full_state_restarts": 8,
        "benchmark_shots": lib.DEFAULT_SHOT_GRID,
        "robustness_shots": 1_500,
        "coherence_shots": 8_000,
        "pulse_shots": 3_000,
    },
}


def profile_config() -> dict[str, Any]:
    name = str(os.environ.get("STUDY_PROFILE", "full")).lower()
    return {"profile_name": name, **PROFILES.get(name, PROFILES["full"])}


def deterministic_diagonal_rows(
    settings: Sequence[lib.MeasurementSetting],
    *,
    weighted_transverse: Sequence[complex],
    z_total: float,
    kernel_dim: int = lib.KERNEL_DIM,
) -> list[dict[str, Any]]:
    exact = lib.diagonal_predict_observables(
        settings,
        weighted_transverse=weighted_transverse,
        z_total=z_total,
        kernel_dim=kernel_dim,
    )
    rows: list[dict[str, Any]] = []
    for index, setting in enumerate(settings):
        for axis in lib.MEASUREMENT_AXES:
            value = float(exact[axis][index])
            rows.append(
                {
                    **setting.as_dict(),
                    "axis": str(axis),
                    "exact_expectation": value,
                    "measured_expectation": value,
                    "counts_plus": None,
                    "shots": 1,
                }
            )
    return rows


def deterministic_joint_rows(
    case: lib.StudyCase,
    settings: Sequence[lib.MeasurementSetting],
) -> list[dict[str, Any]]:
    dataset = lib.exact_joint_dataset(
        case.state,
        settings,
        model=case.model,
        frame=case.frame,
        noise=lib.NoiseModel(shots=1, sample_mode="gaussian"),
        rng=np.random.default_rng(0),
    )
    rows = []
    for row in dataset["rows"]:
        rows.append(
            {
                **row,
                "measured_expectation": float(row["exact_expectation"]),
                "counts_plus": None,
                "shots": 1,
            }
        )
    return rows


def aggregate_metric(rows: Sequence[dict[str, Any]], group_keys: Sequence[str], value_key: str) -> list[dict[str, Any]]:
    grouped: dict[tuple[Any, ...], list[float]] = {}
    for row in rows:
        key = tuple(row[name] for name in group_keys)
        grouped.setdefault(key, []).append(float(row[value_key]))
    out: list[dict[str, Any]] = []
    for key, values in grouped.items():
        arr = np.asarray(values, dtype=float)
        payload = {name: key[index] for index, name in enumerate(group_keys)}
        payload.update(
            {
                f"{value_key}_mean": float(np.mean(arr)),
                f"{value_key}_std": float(np.std(arr)),
                f"{value_key}_median": float(np.median(arr)),
                f"{value_key}_min": float(np.min(arr)),
                f"{value_key}_max": float(np.max(arr)),
                "count": int(arr.size),
            }
        )
        out.append(payload)
    return out


def sample_waveform_payload(waveform: Any, *, duration_s: float) -> dict[str, Any]:
    tlist = np.linspace(0.0, float(duration_s), 2_000)
    samples = np.asarray(waveform.sample(tlist), dtype=np.complex128)
    return {
        "time_s": tlist.tolist(),
        "real": np.real(samples).tolist(),
        "imag": np.imag(samples).tolist(),
    }


def save_artifact(name: str, payload: dict[str, Any]) -> None:
    save_json(ARTIFACTS_DIR / f"{name}.json", payload)


def run_protocol_benchmark(
    benchmark_case: lib.StudyCase,
    *,
    settings_by_protocol: dict[str, list[lib.MeasurementSetting]],
    benchmark_shots: Sequence[int],
    repeats: int,
    seed: int,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    rng = np.random.default_rng(int(seed))
    noisy_rows: list[dict[str, Any]] = []
    exact_rows: list[dict[str, Any]] = []
    for protocol, settings in settings_by_protocol.items():
        exact_dataset_rows = deterministic_diagonal_rows(
            settings,
            weighted_transverse=benchmark_case.weighted_transverse_truth,
            z_total=benchmark_case.z_total_truth,
        )
        fit_exact = lib.infer_weighted_ls(exact_dataset_rows, settings)
        exact_rows.append(
            {
                "protocol": str(protocol),
                "fit_method": "ls",
                "shots": 0,
                "trial": -1,
                "residual_rms": float(fit_exact["residual_rms"]),
                **lib.recoverable_error_summary(
                    fit_exact["weighted_transverse"],
                    benchmark_case.weighted_transverse_truth,
                    inferred_z_total=float(fit_exact["z_total"]),
                    true_z_total=float(benchmark_case.z_total_truth),
                    true_populations=benchmark_case.populations_truth,
                ),
            }
        )
        for shots in benchmark_shots:
            for trial in range(int(repeats)):
                noise = lib.NoiseModel(shots=int(shots), sample_mode="binomial")
                dataset = lib.diagonal_dataset_from_weighted_truth(
                    settings,
                    weighted_transverse=benchmark_case.weighted_transverse_truth,
                    z_total=benchmark_case.z_total_truth,
                    noise=noise,
                    rng=rng,
                )
                for fit_method, fitter in (("ls", lib.infer_weighted_ls), ("mle", lib.infer_weighted_mle)):
                    fit = fitter(dataset["rows"], settings)
                    noisy_rows.append(
                        {
                            "protocol": str(protocol),
                            "fit_method": str(fit_method),
                            "shots": int(shots),
                            "trial": int(trial),
                            "residual_rms": float(fit["residual_rms"]),
                            **lib.recoverable_error_summary(
                                fit["weighted_transverse"],
                                benchmark_case.weighted_transverse_truth,
                                inferred_z_total=float(fit["z_total"]),
                                true_z_total=float(benchmark_case.z_total_truth),
                                true_populations=benchmark_case.populations_truth,
                            ),
                        }
                    )
    return noisy_rows, exact_rows


def run_coherence_witness(
    diagonal_case: lib.StudyCase,
    coherent_case: lib.StudyCase,
    *,
    settings_by_protocol: dict[str, list[lib.MeasurementSetting]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for protocol in ("wait_only", "combined"):
        settings = settings_by_protocol[protocol]
        for case in (diagonal_case, coherent_case):
            deterministic = deterministic_joint_rows(case, settings)
            fit = lib.infer_weighted_ls(deterministic, settings, kernel_dim=int(case.model.n_cav))
            rows.append(
                {
                    "protocol": str(protocol),
                    "case_id": str(case.case_id),
                    "model_class": str(case.model_class),
                    "residual_rms": float(fit["residual_rms"]),
                    "residual_max_abs": float(fit["residual_max_abs"]),
                    **lib.recoverable_error_summary(
                        fit["weighted_transverse"],
                        case.weighted_transverse_truth,
                        inferred_z_total=float(fit["z_total"]),
                        true_z_total=float(case.z_total_truth),
                        true_populations=case.populations_truth,
                    ),
                }
            )
    return rows


def run_full_state_diagnostic(
    benchmark_case: lib.StudyCase,
    *,
    settings: Sequence[lib.MeasurementSetting],
    n_restarts: int,
) -> dict[str, Any]:
    deterministic = deterministic_diagonal_rows(
        settings,
        weighted_transverse=benchmark_case.weighted_transverse_truth,
        z_total=benchmark_case.z_total_truth,
    )
    return lib.run_full_state_restart_diagnostic(
        deterministic,
        settings,
        n_restarts=int(n_restarts),
    )


def run_robustness_sweep(
    benchmark_case: lib.StudyCase,
    *,
    settings_by_protocol: dict[str, list[lib.MeasurementSetting]],
    repeats: int,
    shots: int,
    seed: int,
) -> list[dict[str, Any]]:
    rng = np.random.default_rng(int(seed))
    scenarios = {
        "shot_only": lib.NoiseModel(shots=shots, sample_mode="binomial"),
        "rotation_error": lib.NoiseModel(shots=shots, sample_mode="binomial", rotation_sigma_rad=0.02),
        "displacement_error": lib.NoiseModel(
            shots=shots,
            sample_mode="binomial",
            displacement_rel_sigma=0.03,
            displacement_phase_sigma_rad=0.03,
        ),
        "chi_mismatch": lib.NoiseModel(
            shots=shots,
            sample_mode="binomial",
            chi_rel_sigma=0.01,
            chi_prime_rel_sigma=0.05,
        ),
        "wait_decoherence": lib.NoiseModel(
            shots=shots,
            sample_mode="binomial",
            t1_s=40.0e-6,
            t2_s=25.0e-6,
        ),
    }
    rows: list[dict[str, Any]] = []
    for protocol in ("wait_only", "combined"):
        settings = settings_by_protocol[protocol]
        for scenario_name, noise in scenarios.items():
            for trial in range(int(repeats)):
                dataset = lib.diagonal_dataset_from_weighted_truth(
                    settings,
                    weighted_transverse=benchmark_case.weighted_transverse_truth,
                    z_total=benchmark_case.z_total_truth,
                    noise=noise,
                    rng=rng,
                )
                fit = lib.infer_weighted_mle(dataset["rows"], settings)
                rows.append(
                    {
                        "protocol": str(protocol),
                        "scenario": str(scenario_name),
                        "trial": int(trial),
                        "residual_rms": float(fit["residual_rms"]),
                        **lib.recoverable_error_summary(
                            fit["weighted_transverse"],
                            benchmark_case.weighted_transverse_truth,
                            inferred_z_total=float(fit["z_total"]),
                            true_z_total=float(benchmark_case.z_total_truth),
                            true_populations=benchmark_case.populations_truth,
                        ),
                    }
                )
    return rows


def run_pulse_case_benchmarks(
    pulse_cases: Sequence[lib.StudyCase],
    *,
    settings_by_protocol: dict[str, list[lib.MeasurementSetting]],
    shots: int,
    seed: int,
) -> list[dict[str, Any]]:
    rng = np.random.default_rng(int(seed))
    rows: list[dict[str, Any]] = []
    for case in pulse_cases:
        for protocol in ("wait_only", "combined"):
            settings = settings_by_protocol[protocol]
            sample_mode = "gaussian" if int(case.model.n_tr) > 2 else "binomial"
            noise = lib.NoiseModel(shots=int(shots), sample_mode=sample_mode)
            result = lib.protocol_fit_for_case(
                case,
                settings=settings,
                noise=noise,
                rng=rng,
                fit_method="ls",
            )
            rows.append(
                {
                    "case_id": str(case.case_id),
                    "protocol": str(protocol),
                    "model_class": str(case.model_class),
                    "residual_rms": float(result["fit"]["residual_rms"]),
                    "mean_abs_error": float(result["fit_errors"]["mean_abs_error"]),
                    "max_abs_error": float(result["fit_errors"]["max_abs_error"]),
                    "weighted_rmse": float(result["fit_errors"]["weighted_rmse"]),
                    "z_total_error": float(result["fit_errors"]["z_total_error"]),
                }
            )
            save_artifact(f"pulse_fit_{case.case_id}_{protocol}", result)
    return rows


def convergence_checks(
    benchmark_case: lib.StudyCase,
    *,
    settings_by_protocol: dict[str, list[lib.MeasurementSetting]],
) -> dict[str, Any]:
    wait_rows = deterministic_diagonal_rows(
        settings_by_protocol["wait_only"],
        weighted_transverse=benchmark_case.weighted_transverse_truth,
        z_total=benchmark_case.z_total_truth,
        kernel_dim=lib.KERNEL_DIM,
    )
    wait_fit = lib.infer_weighted_ls(wait_rows, settings_by_protocol["wait_only"], kernel_dim=lib.KERNEL_DIM)
    coarse_rows = deterministic_diagonal_rows(
        settings_by_protocol["combined"],
        weighted_transverse=benchmark_case.weighted_transverse_truth,
        z_total=benchmark_case.z_total_truth,
        kernel_dim=lib.KERNEL_DIM,
    )
    fit_k10 = lib.infer_weighted_ls(coarse_rows, settings_by_protocol["combined"], kernel_dim=lib.KERNEL_DIM)
    fit_k12 = lib.infer_weighted_ls(coarse_rows, settings_by_protocol["combined"], kernel_dim=12)
    err_k10 = lib.recoverable_error_summary(
        fit_k10["weighted_transverse"],
        benchmark_case.weighted_transverse_truth,
        inferred_z_total=float(fit_k10["z_total"]),
        true_z_total=float(benchmark_case.z_total_truth),
    )
    err_k12 = lib.recoverable_error_summary(
        fit_k12["weighted_transverse"],
        benchmark_case.weighted_transverse_truth,
        inferred_z_total=float(fit_k12["z_total"]),
        true_z_total=float(benchmark_case.z_total_truth),
    )
    dense_settings = lib.measurement_settings_combined(
        wait_grid_s=tuple(float(x) for x in np.linspace(0.0, 0.35e-6, 17))
    )
    dense_rows = deterministic_diagonal_rows(
        dense_settings,
        weighted_transverse=benchmark_case.weighted_transverse_truth,
        z_total=benchmark_case.z_total_truth,
    )
    dense_fit = lib.infer_weighted_ls(dense_rows, dense_settings)
    return {
        "wait_only_exact_rmse": float(
            lib.recoverable_error_summary(
                wait_fit["weighted_transverse"],
                benchmark_case.weighted_transverse_truth,
                inferred_z_total=float(wait_fit["z_total"]),
                true_z_total=float(benchmark_case.z_total_truth),
            )["weighted_rmse"]
        ),
        "combined_kernel_dim_10_rmse": float(err_k10["weighted_rmse"]),
        "combined_kernel_dim_12_rmse": float(err_k12["weighted_rmse"]),
        "combined_dense_grid_rmse": float(
            lib.recoverable_error_summary(
                dense_fit["weighted_transverse"],
                benchmark_case.weighted_transverse_truth,
                inferred_z_total=float(dense_fit["z_total"]),
                true_z_total=float(benchmark_case.z_total_truth),
            )["weighted_rmse"]
        ),
        "combined_kernel_dim_delta": float(abs(err_k10["weighted_rmse"] - err_k12["weighted_rmse"])),
    }


def plot_protocol_schematic(recommendation_lines: Sequence[str]) -> None:
    apply_plot_style()
    fig, ax = plt.subplots(figsize=(10.0, 4.4))
    ax.axis("off")
    ax.text(0.02, 0.90, "Protocol Flow", fontsize=16, weight="bold")
    ax.text(0.02, 0.72, "black-box SQR  ->  calibrated displacement D(alpha)  ->  known dispersive wait t", fontsize=12)
    ax.text(0.02, 0.58, "->  qubit pre-rotation for X/Y/Z tomography  ->  z readout  ->  inverse fit on recoverable sector data", fontsize=12)
    ax.text(0.02, 0.38, "Key inference target: weighted transverse sector amplitudes u_n = p_n (X_n + i Y_n)", fontsize=12)
    ax.text(0.02, 0.26, "Key limitation: the allowed protocol does not separate p_n and Z_n.", fontsize=12)
    ax.text(0.02, 0.08, "Recommendation:", fontsize=13, weight="bold")
    for index, line in enumerate(recommendation_lines):
        ax.text(0.18, 0.08 - 0.07 * index, line, fontsize=11)
    save_figure(fig, "protocol_schematic")


def plot_single_qubit_baseline(baseline: dict[str, Any]) -> None:
    apply_plot_style()
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    for method, marker in (("ls", "o"), ("mle", "s")):
        rows = [row for row in baseline["summary_rows"] if row["method"] == method]
        shots = [row["shots"] for row in rows]
        means = [row["mean_fidelity"] for row in rows]
        stds = [row["std_fidelity"] for row in rows]
        ax.errorbar(shots, means, yerr=stds, marker=marker, label=method.upper(), capsize=3)
    ax.set_xscale("log")
    ax.set_xlabel("Shots per branch")
    ax.set_ylabel("State fidelity")
    ax.set_title("Single-Qubit Tomography Baseline")
    ax.legend(frameon=False)
    save_figure(fig, "single_qubit_fidelity_scaling")


def plot_protocol_benchmark(rows: Sequence[dict[str, Any]]) -> None:
    apply_plot_style()
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    for protocol, color in (("wait_only", "#4477AA"), ("displacement_only", "#CC6677"), ("combined", "#228833")):
        subset = [row for row in rows if row["protocol"] == protocol and row["fit_method"] == "mle"]
        for summary in aggregate_metric(subset, ("protocol", "shots"), "weighted_rmse"):
            pass
        summary_rows = sorted(
            aggregate_metric(subset, ("protocol", "shots"), "weighted_rmse"),
            key=lambda item: item["shots"],
        )
        ax.plot(
            [item["shots"] for item in summary_rows],
            [item["weighted_rmse_mean"] for item in summary_rows],
            marker="o",
            color=color,
            label=protocol.replace("_", " "),
        )
    ax.set_xscale("log")
    ax.set_xlabel("Shots per setting")
    ax.set_ylabel("Weighted transverse RMSE")
    ax.set_title("Protocol Comparison on the Recoverable Subspace")
    ax.legend(frameon=False)
    save_figure(fig, "protocol_comparison")


def plot_coherence_witness(rows: Sequence[dict[str, Any]]) -> None:
    apply_plot_style()
    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    labels = []
    values = []
    for protocol in ("wait_only", "combined"):
        for case_id in ("ideal_diag_g", "ideal_coherent_ground"):
            row = next(item for item in rows if item["protocol"] == protocol and item["case_id"] == case_id)
            labels.append(f"{protocol}\n{case_id.replace('ideal_', '')}")
            values.append(float(row["residual_rms"]))
    ax.bar(np.arange(len(labels)), values, color=["#4477AA", "#66AADD", "#228833", "#77CC88"])
    ax.set_xticks(np.arange(len(labels)), labels)
    ax.set_ylabel("Diagonal-model fit residual RMS")
    ax.set_title("Cavity-Coherence Witness")
    save_figure(fig, "coherence_witness_residuals")


def plot_robustness(rows: Sequence[dict[str, Any]]) -> None:
    apply_plot_style()
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    scenarios = ["shot_only", "rotation_error", "displacement_error", "chi_mismatch", "wait_decoherence"]
    x = np.arange(len(scenarios))
    width = 0.35
    for offset, protocol in ((-0.5, "wait_only"), (0.5, "combined")):
        summary = {row["scenario"]: row for row in aggregate_metric([item for item in rows if item["protocol"] == protocol], ("scenario",), "weighted_rmse")}
        ax.bar(
            x + offset * width,
            [summary[name]["weighted_rmse_mean"] for name in scenarios],
            width=width,
            label=protocol.replace("_", " "),
        )
    ax.set_xticks(x, [name.replace("_", "\n") for name in scenarios])
    ax.set_ylabel("Weighted transverse RMSE")
    ax.set_title("Noise Robustness")
    ax.legend(frameon=False)
    save_figure(fig, "robustness_sensitivity")


def plot_pulse_case_recovery(rows: Sequence[dict[str, Any]]) -> None:
    apply_plot_style()
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    case_ids = [row["case_id"] for row in rows if row["protocol"] == "wait_only"]
    x = np.arange(len(case_ids))
    width = 0.34
    for offset, protocol in ((-0.5, "wait_only"), (0.5, "combined")):
        subset = {row["case_id"]: row for row in rows if row["protocol"] == protocol}
        ax.bar(
            x + offset * width,
            [subset[case_id]["weighted_rmse"] for case_id in case_ids],
            width=width,
            label=protocol.replace("_", " "),
        )
    ax.set_xticks(x, case_ids, rotation=20, ha="right")
    ax.set_ylabel("Weighted transverse RMSE")
    ax.set_title("Pulse-Level Black-Box Case Recovery")
    ax.legend(frameon=False)
    save_figure(fig, "pulse_case_recovery")


def plot_full_state_nonuniqueness(diagnostic: dict[str, Any]) -> None:
    apply_plot_style()
    rows = diagnostic["rows"]
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.2))
    objectives = [row["objective"] for row in rows]
    axes[0].plot(np.arange(len(rows)), objectives, marker="o")
    axes[0].set_xlabel("Restart")
    axes[0].set_ylabel("Objective value")
    axes[0].set_title("Restart Objectives")
    prob_matrix = np.asarray([row["probabilities"] for row in rows], dtype=float)
    image = axes[1].imshow(prob_matrix, aspect="auto", cmap="viridis")
    axes[1].set_xlabel("Fock sector n")
    axes[1].set_ylabel("Restart")
    axes[1].set_title("Recovered populations p_n")
    fig.colorbar(image, ax=axes[1], fraction=0.045, pad=0.02)
    save_figure(fig, "full_state_nonuniqueness")


def plot_waveforms(pulse_payloads: dict[str, Any]) -> None:
    apply_plot_style()
    optimized = pulse_payloads["optimized_payload"]
    seed = pulse_payloads["seed_payload"]
    optimized_samples = sample_waveform_payload(optimized["waveform"], duration_s=float(optimized["waveform"].duration_s))
    seed_samples = sample_waveform_payload(seed["waveform"], duration_s=float(seed["waveform"].duration_s))
    save_artifact("optimized_waveform", optimized_samples)
    save_artifact("seed_waveform", seed_samples)
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 6.2), sharex="col")
    for column, (title, payload) in enumerate((("Optimized", optimized_samples), ("Seed / short", seed_samples))):
        time_ns = 1.0e9 * np.asarray(payload["time_s"], dtype=float)
        real = np.asarray(payload["real"], dtype=float)
        imag = np.asarray(payload["imag"], dtype=float)
        axes[0, column].plot(time_ns, real, label="I")
        axes[0, column].plot(time_ns, imag, label="Q")
        axes[0, column].set_title(title)
        axes[0, column].set_ylabel("Baseband amplitude")
        axes[0, column].legend(frameon=False)
        dt = float(np.mean(np.diff(np.asarray(payload["time_s"], dtype=float))))
        spectrum = np.fft.fftshift(np.fft.fft(real + 1.0j * imag))
        freqs = np.fft.fftshift(np.fft.fftfreq(real.size, d=dt)) / 1.0e6
        axes[1, column].plot(freqs, np.abs(spectrum) / max(float(np.max(np.abs(spectrum))), 1.0))
        axes[1, column].set_xlabel("Frequency offset (MHz)")
        axes[1, column].set_ylabel("Normalized |FFT|")
    save_figure(fig, "waveform_examples")


def build_summary(
    *,
    identifiability: dict[str, Any],
    baseline: dict[str, Any],
    benchmark_rows: Sequence[dict[str, Any]],
    coherence_rows: Sequence[dict[str, Any]],
    pulse_rows: Sequence[dict[str, Any]],
    robustness_rows: Sequence[dict[str, Any]],
    full_state_diag: dict[str, Any],
    exact_gauge_family: Sequence[dict[str, Any]],
) -> dict[str, Any]:
    best_baseline = max(row["mean_fidelity"] for row in baseline["summary_rows"] if row["method"] == "mle")
    protocol_summary = aggregate_metric(
        [row for row in benchmark_rows if row["fit_method"] == "mle"],
        ("protocol", "shots"),
        "weighted_rmse",
    )
    coherence_wait = next(row for row in coherence_rows if row["protocol"] == "wait_only" and row["case_id"] == "ideal_coherent_ground")
    coherence_combined = next(row for row in coherence_rows if row["protocol"] == "combined" and row["case_id"] == "ideal_coherent_ground")
    pulse_best = min(pulse_rows, key=lambda row: float(row["weighted_rmse"]))
    pulse_worst = max(pulse_rows, key=lambda row: float(row["weighted_rmse"]))
    robustness_summary = aggregate_metric(robustness_rows, ("protocol", "scenario"), "weighted_rmse")
    return {
        "best_single_qubit_mle_mean_fidelity": float(best_baseline),
        "identifiability": identifiability,
        "protocol_summary": protocol_summary,
        "coherence_wait_residual": float(coherence_wait["residual_rms"]),
        "coherence_combined_residual": float(coherence_combined["residual_rms"]),
        "pulse_best_case": pulse_best,
        "pulse_worst_case": pulse_worst,
        "robustness_summary": robustness_summary,
        "full_state_objective_span": float(full_state_diag["objective_span"]),
        "full_state_probability_std_by_sector": full_state_diag["probability_std_by_sector"],
        "exact_gauge_family_count": int(len(exact_gauge_family)),
    }


def write_summary_markdown(summary: dict[str, Any]) -> None:
    lines = [
        "# Fock-Resolved Black-Box SQR Inference",
        "",
        "## Executive Summary",
        f"- Single-qubit baseline MLE reached mean fidelity {summary['best_single_qubit_mle_mean_fidelity']:.4f} at the top of the shot sweep.",
        f"- Wait-only transverse identifiability rank: {summary['identifiability']['wait_only']['transverse_rank']}.",
        f"- Displacement-only transverse identifiability rank: {summary['identifiability']['displacement_only']['transverse_rank']}.",
        f"- Combined coherence residual: {summary['coherence_combined_residual']:.4e} versus wait-only {summary['coherence_wait_residual']:.4e}.",
        f"- Best pulse-level protocol/case: {summary['pulse_best_case']['protocol']} on {summary['pulse_best_case']['case_id']} with weighted RMSE {summary['pulse_best_case']['weighted_rmse']:.4e}.",
        f"- Worst pulse-level protocol/case: {summary['pulse_worst_case']['protocol']} on {summary['pulse_worst_case']['case_id']} with weighted RMSE {summary['pulse_worst_case']['weighted_rmse']:.4e}.",
        f"- Full-state restart objective span: {summary['full_state_objective_span']:.4e}.",
        f"- Explicit exact gauge-family constructions found: {summary['exact_gauge_family_count']}.",
    ]
    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    start = time.perf_counter()
    config = profile_config()
    rng_seed = 20260331
    print(f"Running study profile: {config['profile_name']}", flush=True)

    settings_by_protocol = lib.measurement_settings_by_protocol()
    identifiability = {
        protocol: lib.protocol_identifiability_summary(settings)
        for protocol, settings in settings_by_protocol.items()
    }

    ideal_cases = {case.case_id: case for case in lib.build_ideal_reference_cases()}
    benchmark_case = ideal_cases["ideal_diag_+"]
    diagonal_partner = ideal_cases["ideal_diag_g"]
    coherent_partner = ideal_cases["ideal_coherent_ground"]

    print("Running single-qubit baseline...", flush=True)
    baseline = lib.run_single_qubit_baseline(
        repeats=int(config["baseline_repeats"]),
        sample_mode="binomial",
    )

    print("Running protocol benchmark...", flush=True)
    benchmark_rows, benchmark_exact_rows = run_protocol_benchmark(
        benchmark_case,
        settings_by_protocol=settings_by_protocol,
        benchmark_shots=config["benchmark_shots"],
        repeats=int(config["benchmark_repeats"]),
        seed=rng_seed,
    )

    print("Running full-state non-uniqueness diagnostic...", flush=True)
    full_state_diag = run_full_state_diagnostic(
        diagonal_partner,
        settings=settings_by_protocol["combined"],
        n_restarts=int(config["full_state_restarts"]),
    )
    exact_gauge_family = lib.construct_exact_gauge_family(
        diagonal_partner.weighted_transverse_truth,
        diagonal_partner.z_total_truth,
    )

    print("Running coherence witness benchmark...", flush=True)
    coherence_rows = run_coherence_witness(
        diagonal_partner,
        coherent_partner,
        settings_by_protocol=settings_by_protocol,
    )

    print("Building pulse-level black-box cases...", flush=True)
    pulse_cases, pulse_payloads = lib.build_pulse_level_cases()

    print("Running pulse-level inference benchmarks...", flush=True)
    pulse_rows = run_pulse_case_benchmarks(
        pulse_cases,
        settings_by_protocol=settings_by_protocol,
        shots=int(config["pulse_shots"]),
        seed=rng_seed + 1,
    )

    print("Running robustness sweep...", flush=True)
    robustness_rows = run_robustness_sweep(
        benchmark_case,
        settings_by_protocol=settings_by_protocol,
        repeats=int(config["robustness_repeats"]),
        shots=int(config["robustness_shots"]),
        seed=rng_seed + 2,
    )

    convergence = convergence_checks(
        benchmark_case,
        settings_by_protocol=settings_by_protocol,
    )

    summary = build_summary(
        identifiability=identifiability,
        baseline=baseline,
        benchmark_rows=benchmark_rows,
        coherence_rows=coherence_rows,
        pulse_rows=pulse_rows,
        robustness_rows=robustness_rows,
        full_state_diag=full_state_diag,
        exact_gauge_family=exact_gauge_family,
    )

    recommendation_lines = [
        "Use wait-only when the goal is the cleanest recovery of the weighted transverse sector data.",
        "Add displacement-plus-wait when you need a diagnostic for cavity coherences or model mismatch.",
        "Do not claim recovery of p_n or sector-resolved Z_n without an extra cavity-sensitive measurement primitive.",
    ]

    print("Generating figures...", flush=True)
    plot_protocol_schematic(recommendation_lines)
    plot_single_qubit_baseline(baseline)
    plot_protocol_benchmark(benchmark_rows)
    plot_coherence_witness(coherence_rows)
    plot_robustness(robustness_rows)
    plot_pulse_case_recovery(pulse_rows)
    plot_full_state_nonuniqueness(full_state_diag)
    plot_waveforms(pulse_payloads)

    results_payload = {
        "study": STUDY_DIR.name,
        "profile": config["profile_name"],
        "generated_at": time.strftime("%Y-%m-%d"),
        "identifiability": identifiability,
        "baseline": baseline,
        "benchmark_rows": benchmark_rows,
        "benchmark_exact_rows": benchmark_exact_rows,
        "coherence_rows": coherence_rows,
        "pulse_rows": pulse_rows,
        "full_state_diagnostic": full_state_diag,
        "exact_gauge_family": exact_gauge_family,
        "robustness_rows": robustness_rows,
        "convergence": convergence,
        "case_catalog": [lib.serialize_case(case) for case in ideal_cases.values()] + [lib.serialize_case(case) for case in pulse_cases],
        "summary": summary,
        "recommendation_lines": recommendation_lines,
    }
    validation_inputs = {
        "identifiability": identifiability,
        "coherence_rows": coherence_rows,
        "convergence": convergence,
        "summary": summary,
        "exact_gauge_family_count": len(exact_gauge_family),
    }
    save_json(RESULTS_JSON, results_payload)
    save_json(SUMMARY_JSON, summary)
    save_json(VALIDATION_HINTS_JSON, validation_inputs)
    write_summary_markdown(summary)

    elapsed = time.perf_counter() - start
    print(f"Completed study in {elapsed:.1f} s", flush=True)
    print(f"Saved results to {RESULTS_JSON}", flush=True)


if __name__ == "__main__":
    main()
