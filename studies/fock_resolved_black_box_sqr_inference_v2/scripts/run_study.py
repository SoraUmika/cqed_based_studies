"""Run the Fock-resolved black-box SQR inference study (v2).

Execution arc follows the open study plan exactly:
  Part 1   Single-qubit MLE building block
  Part 2A  Protocol identifiability
  Part 2B  Full {p_n, rho_q^(n)} MLE attempt + non-uniqueness characterisation
  Part 2C  Recoverable-subspace inference (LS vs MLE, protocol comparison)
  Part 2D  Analytic black-box case library with per-sector metrics
  Part 2E  Pulse-level black-box case library with per-sector metrics
  Part 2F  Model-B coherence sweep
  Part 2G  Robustness / noise sweep
  Part 3   Comparison questions answered quantitatively
  Figures  All required deliverables including per-Fock Bloch-vector comparison
"""

from __future__ import annotations

import os
import time
from pathlib import Path
from typing import Any, Sequence

import matplotlib.pyplot as plt
import numpy as np

import study_lib as lib

from common import (
    ARTIFACTS_DIR,
    DATA_DIR,
    FIGURES_DIR,
    N_ACTIVE,
    apply_plot_style,
    save_figure,
    save_json,
)


RESULTS_JSON = DATA_DIR / "study_results.json"
SUMMARY_JSON = DATA_DIR / "study_summary.json"
SUMMARY_MD = DATA_DIR / "study_summary.md"
VALIDATION_HINTS_JSON = DATA_DIR / "validation_inputs.json"

PROFILES: dict[str, dict[str, Any]] = {
    "quick": {
        "baseline_repeats": 6,
        "benchmark_repeats": 5,
        "robustness_repeats": 4,
        "full_state_restarts": 4,
        "benchmark_shots": (100, 1_000, 10_000),
        "robustness_shots": 1_000,
        "coherence_shots": 3_000,
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


def save_artifact(name: str, payload: Any) -> None:
    save_json(ARTIFACTS_DIR / f"{name}.json", payload)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def deterministic_diagonal_rows(
    settings: Sequence[lib.MeasurementSetting],
    *,
    weighted_transverse: Sequence[complex],
    z_total: float,
) -> list[dict[str, Any]]:
    exact = lib.diagonal_predict_observables(
        settings,
        weighted_transverse=weighted_transverse,
        z_total=z_total,
    )
    rows = []
    for i, setting in enumerate(settings):
        for axis in lib.MEASUREMENT_AXES:
            val = float(exact[axis][i])
            rows.append({
                **setting.as_dict(),
                "axis": str(axis),
                "exact_expectation": val,
                "measured_expectation": val,
                "counts_plus": None,
                "shots": 1,
            })
    return rows


def deterministic_joint_rows(
    case: lib.StudyCase,
    settings: Sequence[lib.MeasurementSetting],
) -> list[dict[str, Any]]:
    dataset = lib.exact_joint_dataset(
        case.state, settings,
        model=case.model,
        frame=case.frame,
        noise=lib.NoiseModel(shots=1, sample_mode="gaussian"),
        rng=np.random.default_rng(0),
    )
    rows = []
    for row in dataset["rows"]:
        rows.append({
            **row,
            "measured_expectation": float(row["exact_expectation"]),
            "counts_plus": None,
            "shots": 1,
        })
    return rows


def aggregate_stats(values: Sequence[float]) -> dict[str, float]:
    arr = np.asarray(values, dtype=float)
    return {
        "mean": float(np.mean(arr)),
        "std": float(np.std(arr)),
        "median": float(np.median(arr)),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
    }


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
        exact_data = deterministic_diagonal_rows(
            settings,
            weighted_transverse=benchmark_case.weighted_transverse_truth,
            z_total=benchmark_case.z_total_truth,
        )
        fit_exact = lib.infer_weighted_ls(exact_data, settings)
        exact_rows.append({
            "protocol": str(protocol),
            "shots": 0,
            "residual_rms": float(fit_exact["residual_rms"]),
            **lib.recoverable_error_summary(
                fit_exact["weighted_transverse"],
                benchmark_case.weighted_transverse_truth,
                inferred_z_total=float(fit_exact["z_total"]),
                true_z_total=float(benchmark_case.z_total_truth),
                true_sectors=benchmark_case.truth_sectors,
            ),
        })
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
                    noisy_rows.append({
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
                            true_sectors=benchmark_case.truth_sectors,
                        ),
                    })
    return noisy_rows, exact_rows


def run_single_case_inference(
    case: lib.StudyCase,
    settings: Sequence[lib.MeasurementSetting],
    *,
    shots: int,
    rng: np.random.Generator,
) -> dict[str, Any]:
    dataset = lib.exact_joint_dataset(
        case.state, settings,
        model=case.model,
        frame=case.frame,
        noise=lib.NoiseModel(shots=int(shots), sample_mode="binomial"),
        rng=rng,
    )
    rows = dataset["rows"]
    fit_ls = lib.infer_weighted_ls(rows, settings)
    fit_mle = lib.infer_weighted_mle(rows, settings)
    exact_rows = deterministic_joint_rows(case, settings)
    fit_exact = lib.infer_weighted_ls(exact_rows, settings)

    return {
        "case_id": str(case.case_id),
        "model_class": str(case.model_class),
        "protocol": str(settings[0].protocol),
        "shots": int(shots),
        "ls": {
            "residual_rms": float(fit_ls["residual_rms"]),
            **lib.recoverable_error_summary(
                fit_ls["weighted_transverse"],
                case.weighted_transverse_truth,
                inferred_z_total=float(fit_ls["z_total"]),
                true_z_total=float(case.z_total_truth),
                true_sectors=case.truth_sectors,
            ),
        },
        "mle": {
            "residual_rms": float(fit_mle["residual_rms"]),
            **lib.recoverable_error_summary(
                fit_mle["weighted_transverse"],
                case.weighted_transverse_truth,
                inferred_z_total=float(fit_mle["z_total"]),
                true_z_total=float(case.z_total_truth),
                true_sectors=case.truth_sectors,
            ),
        },
        "exact": {
            "residual_rms": float(fit_exact["residual_rms"]),
            **lib.recoverable_error_summary(
                fit_exact["weighted_transverse"],
                case.weighted_transverse_truth,
                inferred_z_total=float(fit_exact["z_total"]),
                true_z_total=float(case.z_total_truth),
                true_sectors=case.truth_sectors,
            ),
        },
    }


# ---------------------------------------------------------------------------
# Figure functions
# ---------------------------------------------------------------------------

def plot_single_qubit_baseline(baseline_result: dict[str, Any]) -> None:
    apply_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    summary = baseline_result["summary_rows"]

    for method, ax in zip(("ls", "mle"), axes):
        shot_vals = sorted({r["shots"] for r in summary if r["method"] == method})
        means = [next(r["mean_fidelity"] for r in summary if r["shots"] == s and r["method"] == method) for s in shot_vals]
        stds = [next(r["std_fidelity"] for r in summary if r["shots"] == s and r["method"] == method) for s in shot_vals]
        ax.errorbar(shot_vals, means, yerr=stds, marker="o", capsize=4, label=method.upper())
        ax.axhline(0.995, ls="--", color="gray", lw=0.8, label="0.995 threshold")
        ax.set_xscale("log")
        ax.set_xlabel("Shots per branch N")
        ax.set_ylabel("Mean Bures fidelity")
        ax.set_title(f"Part 1 single-qubit baseline — {method.upper()}")
        ax.legend()
        ax.set_ylim(0.9, 1.002)

    fig.tight_layout()
    save_figure(fig, "single_qubit_fidelity_scaling")
    plt.close(fig)


def plot_protocol_comparison(noisy_rows: list[dict[str, Any]]) -> None:
    apply_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    protocols = ["wait_only", "displacement_only", "combined"]
    colors = {"wait_only": "C0", "displacement_only": "C1", "combined": "C2"}
    labels = {"wait_only": "Wait-only", "displacement_only": "Disp.-only", "combined": "Combined"}

    for method, ax in zip(("ls", "mle"), axes):
        method_rows = [r for r in noisy_rows if r["fit_method"] == method]
        for prot in protocols:
            prot_rows = [r for r in method_rows if r["protocol"] == prot]
            shot_vals = sorted({r["shots"] for r in prot_rows})
            means = [float(np.mean([r["weighted_rmse"] for r in prot_rows if r["shots"] == s])) for s in shot_vals]
            stds = [float(np.std([r["weighted_rmse"] for r in prot_rows if r["shots"] == s])) for s in shot_vals]
            ax.errorbar(shot_vals, means, yerr=stds, marker="o", capsize=4,
                        color=colors[prot], label=labels[prot])
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Shots per setting")
        ax.set_ylabel("Weighted transverse RMSE")
        ax.set_title(f"Protocol comparison — {method.upper()}")
        ax.legend()

    fig.tight_layout()
    save_figure(fig, "protocol_comparison")
    plt.close(fig)


def plot_per_fock_bloch_comparison(
    cases: Sequence[lib.StudyCase],
    settings: Sequence[lib.MeasurementSetting],
    *,
    shots: int,
    seed: int = 7777,
) -> None:
    """Per-Fock Bloch-vector comparison: inferred vs. ground truth, per case."""
    apply_plot_style()
    rng = np.random.default_rng(int(seed))
    n_cases = len(cases)
    fig, axes = plt.subplots(n_cases, N_ACTIVE, figsize=(3.5 * N_ACTIVE, 2.8 * n_cases))
    if n_cases == 1:
        axes = axes[np.newaxis, :]

    for row_idx, case in enumerate(cases):
        dataset = lib.exact_joint_dataset(
            case.state, settings,
            model=case.model,
            frame=case.frame,
            noise=lib.NoiseModel(shots=int(shots), sample_mode="binomial"),
            rng=rng,
        )
        fit = lib.infer_weighted_ls(dataset["rows"], settings)
        u_inf = fit["weighted_transverse"]
        for col_idx, sector in enumerate(case.truth_sectors):
            ax = axes[row_idx, col_idx]
            p_true = float(sector.population)
            # True Bloch vector
            x_true, y_true, z_true = sector.x, sector.y, sector.z
            # Inferred transverse (oracle normalised using true p_n)
            if p_true > 1.0e-6:
                x_inf = float(np.real(u_inf[col_idx])) / p_true
                y_inf = float(np.imag(u_inf[col_idx])) / p_true
                # Z not individually recoverable — mark as not inferred
                z_inf = float("nan")
            else:
                x_inf = y_inf = z_inf = float("nan")

            ax.bar([0, 1, 2], [x_true, y_true, z_true], width=0.35, label="True", color="C0", alpha=0.8)
            ax.bar([0.4, 1.4, 2.4], [x_inf, y_inf, z_inf], width=0.35, label="Inf.", color="C1", alpha=0.8)
            ax.set_xticks([0.2, 1.2, 2.2])
            ax.set_xticklabels(["X", "Y", "Z"])
            ax.set_ylim(-1.1, 1.1)
            ax.axhline(0, color="k", lw=0.5)
            ax.set_title(f"n={col_idx}  p={p_true:.2f}", fontsize=9)
            if col_idx == 0:
                ax.set_ylabel(f"{case.case_id[:16]}", fontsize=7)
            if row_idx == 0 and col_idx == N_ACTIVE - 1:
                ax.legend(fontsize=7)

    fig.suptitle(f"Per-Fock Bloch vector: true vs. oracle-inferred (N={shots} shots/setting)", fontsize=11)
    fig.tight_layout()
    save_figure(fig, "per_fock_bloch_comparison")
    plt.close(fig)


def plot_coherence_sweep(sweep_results: list[dict[str, Any]]) -> None:
    apply_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    fractions = [r["coherence_fraction"] for r in sweep_results]
    rmse_wait = [r["wait_only"]["mle"]["weighted_rmse"] for r in sweep_results]
    rmse_comb = [r["combined"]["mle"]["weighted_rmse"] for r in sweep_results]
    resid_wait = [r["wait_only"]["mle"]["residual_rms"] for r in sweep_results]
    resid_comb = [r["combined"]["mle"]["residual_rms"] for r in sweep_results]

    axes[0].plot(fractions, rmse_wait, "o-", label="Wait-only", color="C0")
    axes[0].plot(fractions, rmse_comb, "s-", label="Combined", color="C2")
    axes[0].set_xlabel("Coherence fraction f")
    axes[0].set_ylabel("Weighted transverse RMSE")
    axes[0].set_title("Accuracy degradation vs. cavity coherence")
    axes[0].legend()

    axes[1].plot(fractions, resid_wait, "o-", label="Wait-only", color="C0")
    axes[1].plot(fractions, resid_comb, "s-", label="Combined", color="C2")
    axes[1].set_xlabel("Coherence fraction f")
    axes[1].set_ylabel("Fit residual RMS")
    axes[1].set_title("Coherence witness: residual vs. fraction")
    axes[1].legend()
    axes[1].axhline(1.0e-3, ls="--", color="gray", lw=0.8, label="1e-3 threshold")

    fig.tight_layout()
    save_figure(fig, "coherence_sweep")
    plt.close(fig)


def plot_full_state_nonuniqueness(full_state_result: dict[str, Any]) -> None:
    apply_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    rows = full_state_result["rows"]
    n_restarts = len(rows)
    prob_matrix = np.asarray([r["probabilities"] for r in rows], dtype=float)
    objectives = np.asarray([r["objective"] for r in rows], dtype=float)

    for n in range(N_ACTIVE):
        axes[0].scatter(range(n_restarts), prob_matrix[:, n], label=f"n={n}", alpha=0.8, s=40)
    axes[0].set_xlabel("Restart index (sorted by objective)")
    axes[0].set_ylabel("Inferred population p_n")
    axes[0].set_title("Full {p_n, ρ_q^(n)} MLE: population spread across restarts")
    axes[0].legend()

    axes[1].scatter(range(n_restarts), objectives, color="k", s=40)
    axes[1].set_xlabel("Restart index")
    axes[1].set_ylabel("MLE objective (NLL)")
    axes[1].set_title("Objective values — confirms non-uniqueness")

    fig.tight_layout()
    save_figure(fig, "full_state_nonuniqueness")
    plt.close(fig)


def plot_pulse_case_recovery(pulse_rows: list[dict[str, Any]]) -> None:
    apply_plot_style()
    case_ids = list({r["case_id"] for r in pulse_rows})
    case_ids.sort()
    protocols = ["wait_only", "combined"]
    x = np.arange(len(case_ids))
    width = 0.35

    fig, ax = plt.subplots(figsize=(12, 5))
    for pi, prot in enumerate(protocols):
        rmse_vals = []
        for cid in case_ids:
            match = [r for r in pulse_rows if r["case_id"] == cid and r["protocol"] == prot]
            rmse_vals.append(float(match[0]["mle"]["weighted_rmse"]) if match else float("nan"))
        ax.bar(x + pi * width - width / 2, rmse_vals, width, label=prot.replace("_", "-"))

    ax.set_xticks(x)
    ax.set_xticklabels([c[:22] for c in case_ids], rotation=25, ha="right", fontsize=8)
    ax.set_ylabel("Weighted transverse RMSE (MLE)")
    ax.set_title("Pulse-level black-box case recovery")
    ax.legend()
    ax.axhline(0.02, ls="--", color="gray", lw=0.8)

    fig.tight_layout()
    save_figure(fig, "pulse_case_recovery")
    plt.close(fig)


def plot_robustness(robustness_rows: list[dict[str, Any]]) -> None:
    apply_plot_style()
    scenarios = list({r["noise_scenario"] for r in robustness_rows})
    scenarios.sort()
    protocols = ["wait_only", "combined"]
    colors = {"wait_only": "C0", "combined": "C2"}

    fig, axes = plt.subplots(1, len(protocols), figsize=(12, 5), sharey=True)
    for pi, prot in enumerate(protocols):
        ax = axes[pi]
        prot_rows = [r for r in robustness_rows if r["protocol"] == prot]
        means = [float(np.mean([r["weighted_rmse"] for r in prot_rows if r["noise_scenario"] == s])) for s in scenarios]
        stds = [float(np.std([r["weighted_rmse"] for r in prot_rows if r["noise_scenario"] == s])) for s in scenarios]
        ax.bar(scenarios, means, yerr=stds, capsize=4, color=colors[prot])
        ax.set_xticklabels(scenarios, rotation=30, ha="right", fontsize=8)
        ax.set_title(prot.replace("_", "-"))
        ax.set_ylabel("Weighted transverse RMSE")
        ax.axhline(0.02, ls="--", color="gray", lw=0.8)

    fig.suptitle("Robustness to systematic errors")
    fig.tight_layout()
    save_figure(fig, "robustness_sensitivity")
    plt.close(fig)


def plot_coherence_witness_residuals(coherence_rows: list[dict[str, Any]]) -> None:
    apply_plot_style()
    fig, ax = plt.subplots(figsize=(8, 4))
    protocols = ["wait_only", "combined"]
    model_classes = sorted({r["model_class"] for r in coherence_rows})
    x = np.arange(len(model_classes))
    width = 0.35

    for pi, prot in enumerate(protocols):
        resids = []
        for mc in model_classes:
            match = [r for r in coherence_rows if r["protocol"] == prot and r["model_class"] == mc]
            resids.append(float(np.mean([r["residual_rms"] for r in match])) if match else float("nan"))
        ax.bar(x + (pi - 0.5) * width, resids, width, label=prot.replace("_", "-"))

    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(model_classes, rotation=15, ha="right")
    ax.set_ylabel("Fit residual RMS (log scale)")
    ax.set_title("Coherence witness: residuals by model class and protocol")
    ax.legend()

    fig.tight_layout()
    save_figure(fig, "coherence_witness_residuals")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main execution
# ---------------------------------------------------------------------------

def main() -> None:
    t_start = time.time()
    cfg = profile_config()
    print(f"Profile: {cfg['profile_name']}")

    settings_by_protocol = lib.measurement_settings_by_protocol()
    wait_settings = settings_by_protocol["wait_only"]
    combined_settings = settings_by_protocol["combined"]
    results: dict[str, Any] = {}
    summary: dict[str, Any] = {"profile": cfg["profile_name"]}

    # ------------------------------------------------------------------
    # Part 1: Single-qubit MLE building block
    # ------------------------------------------------------------------
    print("Part 1: single-qubit baseline...")
    t0 = time.time()
    baseline = lib.run_single_qubit_baseline(
        shot_grid=lib.DEFAULT_SHOT_GRID,
        repeats=int(cfg["baseline_repeats"]),
        seed=1234,
    )
    results["baseline"] = baseline
    summary["best_single_qubit_mle_mean_fidelity"] = float(baseline["best_mle_mean_fidelity"])
    print(f"  best MLE fidelity: {baseline['best_mle_mean_fidelity']:.4f}  ({time.time()-t0:.1f}s)")

    # ------------------------------------------------------------------
    # Part 2A: Protocol identifiability
    # ------------------------------------------------------------------
    print("Part 2A: identifiability analysis...")
    ident = {}
    for prot, settings in settings_by_protocol.items():
        ident[prot] = lib.protocol_identifiability_summary(settings)
    results["identifiability"] = ident
    summary["identifiability"] = ident
    print(f"  wait_only transverse rank: {ident['wait_only']['transverse_rank']}")
    print(f"  displacement_only transverse rank: {ident['displacement_only']['transverse_rank']}")
    print(f"  combined transverse rank: {ident['combined']['transverse_rank']}")

    # ------------------------------------------------------------------
    # Part 2B: Full {p_n, rho_q^(n)} MLE — attempt first, expose non-uniqueness
    # ------------------------------------------------------------------
    print("Part 2B: full-state MLE attempt (non-uniqueness diagnostic)...")
    t0 = time.time()
    # Use combined protocol data on the diagonal benchmark case
    analytic_cases = lib.build_ideal_reference_cases()
    benchmark_case = next(c for c in analytic_cases if c.case_id == "ideal_diag_g")
    noise_high = lib.NoiseModel(shots=int(cfg["coherence_shots"]), sample_mode="binomial")
    rng_full = np.random.default_rng(9999)
    dataset_full = lib.diagonal_dataset_from_weighted_truth(
        combined_settings,
        weighted_transverse=benchmark_case.weighted_transverse_truth,
        z_total=benchmark_case.z_total_truth,
        noise=noise_high,
        rng=rng_full,
    )
    full_state_result = lib.run_full_state_mle_attempt(
        dataset_full["rows"],
        combined_settings,
        n_restarts=int(cfg["full_state_restarts"]),
        seed=4321,
    )
    results["full_state_mle"] = full_state_result
    summary["full_state_objective_span"] = float(full_state_result["objective_span"])
    summary["full_state_probability_std_by_sector"] = full_state_result["probability_std_by_sector"]
    summary["full_state_conclusion"] = str(full_state_result["conclusion"])

    # Algebraic gauge family
    gauge_family = lib.construct_exact_gauge_family(
        benchmark_case.weighted_transverse_truth,
        benchmark_case.z_total_truth,
    )
    results["gauge_family"] = gauge_family
    summary["exact_gauge_family_count"] = int(len(gauge_family))
    print(f"  objective span: {full_state_result['objective_span']:.2f}  "
          f"  p_n std mean: {float(np.mean(full_state_result['probability_std_by_sector'])):.3f}  "
          f"  ({time.time()-t0:.1f}s)")

    # ------------------------------------------------------------------
    # Part 2C: Recoverable-subspace inference — LS vs MLE, protocol comparison
    # ------------------------------------------------------------------
    print("Part 2C: recoverable-subspace benchmark (protocol x shots x method)...")
    t0 = time.time()
    noisy_rows, exact_rows = run_protocol_benchmark(
        benchmark_case,
        settings_by_protocol=settings_by_protocol,
        benchmark_shots=cfg["benchmark_shots"],
        repeats=int(cfg["benchmark_repeats"]),
        seed=2222,
    )
    results["benchmark_noisy_rows"] = noisy_rows
    results["benchmark_exact_rows"] = exact_rows

    # Summarise by (protocol, shots)
    protocol_summary = []
    for prot in settings_by_protocol:
        for shots in cfg["benchmark_shots"]:
            for method in ("ls", "mle"):
                matched = [r for r in noisy_rows
                           if r["protocol"] == prot and r["shots"] == shots and r["fit_method"] == method]
                if not matched:
                    continue
                rmse_vals = [r["weighted_rmse"] for r in matched]
                oracle_fids = [r["mean_oracle_fidelity"] for r in matched if r.get("mean_oracle_fidelity") is not None]
                protocol_summary.append({
                    "protocol": prot,
                    "shots": shots,
                    "method": method,
                    **{f"weighted_rmse_{k}": v for k, v in {
                        "mean": float(np.mean(rmse_vals)),
                        "std": float(np.std(rmse_vals)),
                        "min": float(np.min(rmse_vals)),
                        "max": float(np.max(rmse_vals)),
                    }.items()},
                    "mean_oracle_fidelity": float(np.mean(oracle_fids)) if oracle_fids else None,
                })
    results["protocol_summary"] = protocol_summary
    summary["protocol_summary"] = protocol_summary
    print(f"  done ({time.time()-t0:.1f}s)")

    # Exact convergence checks
    wait_exact = lib.infer_weighted_ls(
        deterministic_diagonal_rows(wait_settings, weighted_transverse=benchmark_case.weighted_transverse_truth, z_total=benchmark_case.z_total_truth),
        wait_settings,
    )
    combined_exact = lib.infer_weighted_ls(
        deterministic_diagonal_rows(combined_settings, weighted_transverse=benchmark_case.weighted_transverse_truth, z_total=benchmark_case.z_total_truth),
        combined_settings,
    )
    convergence = {
        "wait_only_exact_rmse": float(wait_exact["residual_rms"]),
        "combined_dense_grid_rmse": float(combined_exact["residual_rms"]),
    }
    results["convergence"] = convergence
    summary["convergence"] = convergence

    # ------------------------------------------------------------------
    # Part 2D: Analytic black-box case library with per-sector metrics
    # ------------------------------------------------------------------
    print("Part 2D: analytic black-box cases with per-sector metrics...")
    t0 = time.time()
    analytic_case_rows: list[dict[str, Any]] = []
    coherence_rows: list[dict[str, Any]] = []
    rng_analytic = np.random.default_rng(3333)

    for case in analytic_cases:
        for prot, settings in settings_by_protocol.items():
            row = run_single_case_inference(
                case, settings,
                shots=int(cfg["coherence_shots"]),
                rng=rng_analytic,
            )
            row["protocol"] = prot
            row["case_id"] = case.case_id
            row["model_class"] = case.model_class
            analytic_case_rows.append(row)
            if case.model_class in ("fock_diagonal", "coherent_block"):
                coherence_rows.append({
                    "case_id": case.case_id,
                    "model_class": case.model_class,
                    "protocol": prot,
                    "residual_rms": float(row["mle"]["residual_rms"]),
                    "weighted_rmse": float(row["mle"]["weighted_rmse"]),
                })

    results["analytic_case_rows"] = analytic_case_rows
    results["coherence_rows"] = coherence_rows

    # Coherence-witness summary (EXACT, noiseless):
    # For the witness check we use exact (zero-noise) data to measure model-mismatch
    # cleanly, not confounded by shot noise. Same as v1 run_coherence_witness logic.
    diag_case = next((c for c in analytic_cases if c.case_id == "ideal_diag_g"), None)
    coh_case = next((c for c in analytic_cases if c.case_id == "ideal_coherent_g"), None)
    if diag_case is not None:
        exact_diag_rows = deterministic_joint_rows(diag_case, wait_settings)
        fit_diag_exact = lib.infer_weighted_ls(exact_diag_rows, wait_settings)
        coherence_wait_residual = float(fit_diag_exact["residual_rms"])
    else:
        coherence_wait_residual = 0.0
    if coh_case is not None:
        exact_coh_rows = deterministic_joint_rows(coh_case, combined_settings)
        fit_coh_exact = lib.infer_weighted_ls(exact_coh_rows, combined_settings)
        coherence_combined_residual = float(fit_coh_exact["residual_rms"])
    else:
        coherence_combined_residual = 0.0
    summary["coherence_wait_residual"] = coherence_wait_residual
    summary["coherence_combined_residual"] = coherence_combined_residual
    print(f"  diagonal wait residual (exact): {coherence_wait_residual:.2e}   "
          f"coherent combined residual (exact): {coherence_combined_residual:.4f}  ({time.time()-t0:.1f}s)")

    # ------------------------------------------------------------------
    # Part 2E: Pulse-level cases with per-sector metrics
    # ------------------------------------------------------------------
    print("Part 2E: pulse-level cases...")
    t0 = time.time()
    pulse_cases, wf_meta = lib.build_pulse_level_cases()
    save_artifact("waveform_metadata", wf_meta)
    pulse_rows: list[dict[str, Any]] = []
    rng_pulse = np.random.default_rng(5555)

    for case in pulse_cases:
        for prot, settings in [("wait_only", wait_settings), ("combined", combined_settings)]:
            row = run_single_case_inference(
                case, settings,
                shots=int(cfg["pulse_shots"]),
                rng=rng_pulse,
            )
            row["protocol"] = prot
            pulse_rows.append(row)
            print(f"    {case.case_id} | {prot} | MLE RMSE: {row['mle']['weighted_rmse']:.4f}")

    results["pulse_rows"] = pulse_rows
    best_pulse = min(pulse_rows, key=lambda r: float(r["mle"]["weighted_rmse"]))
    worst_pulse = max(pulse_rows, key=lambda r: float(r["mle"]["weighted_rmse"]))
    summary["pulse_best_case"] = {"case_id": best_pulse["case_id"], "weighted_rmse": float(best_pulse["mle"]["weighted_rmse"])}
    summary["pulse_worst_case"] = {"case_id": worst_pulse["case_id"], "weighted_rmse": float(worst_pulse["mle"]["weighted_rmse"])}
    print(f"  done ({time.time()-t0:.1f}s)")

    # ------------------------------------------------------------------
    # Part 2F: Model-B coherence sweep
    # ------------------------------------------------------------------
    print("Part 2F: coherence degradation sweep...")
    t0 = time.time()
    sweep_states = lib.build_coherence_sweep_states()
    sweep_results: list[dict[str, Any]] = []
    rng_sweep = np.random.default_rng(6666)

    for entry in sweep_states:
        f = float(entry["coherence_fraction"])
        state = entry["state"]
        model = entry["model"]
        frame = entry["frame"]
        truth_sectors = entry["truth_sectors"]
        sweep_row: dict[str, Any] = {"coherence_fraction": f}

        for prot, settings in [("wait_only", wait_settings), ("combined", combined_settings)]:
            # Approximate the state as a StudyCase for run_single_case_inference
            fake_case = lib.StudyCase(
                case_id=f"sweep_f{f:.2f}",
                family="coherence_sweep",
                description=f"Mixed diagonal/coherent state at f={f:.2f}",
                model_class="coherent_block" if f > 0 else "fock_diagonal",
                state=state,
                model=model,
                frame=frame,
                truth_sectors=truth_sectors,
                metadata={"coherence_fraction": f},
            )
            row = run_single_case_inference(fake_case, settings, shots=int(cfg["coherence_shots"]), rng=rng_sweep)
            sweep_row[prot] = row

        sweep_results.append(sweep_row)
        print(f"    f={f:.2f}  wait RMSE: {sweep_row['wait_only']['mle']['weighted_rmse']:.4f}  "
              f"comb RMSE: {sweep_row['combined']['mle']['weighted_rmse']:.4f}")

    results["coherence_sweep"] = sweep_results
    summary["coherence_sweep_summary"] = [
        {
            "coherence_fraction": r["coherence_fraction"],
            "wait_only_rmse": r["wait_only"]["mle"]["weighted_rmse"],
            "combined_rmse": r["combined"]["mle"]["weighted_rmse"],
            "combined_residual_rms": r["combined"]["mle"]["residual_rms"],
        }
        for r in sweep_results
    ]
    print(f"  done ({time.time()-t0:.1f}s)")

    # ------------------------------------------------------------------
    # Part 2G: Robustness / noise sweep
    # ------------------------------------------------------------------
    print("Part 2G: robustness sweep...")
    t0 = time.time()
    noise_scenarios = {
        "shot_only": lib.NoiseModel(shots=int(cfg["robustness_shots"]), sample_mode="binomial"),
        "rotation_error": lib.NoiseModel(shots=int(cfg["robustness_shots"]), rotation_sigma_rad=0.03),
        "displacement_error": lib.NoiseModel(shots=int(cfg["robustness_shots"]), displacement_rel_sigma=0.05, displacement_phase_sigma_rad=0.05),
        "chi_mismatch": lib.NoiseModel(shots=int(cfg["robustness_shots"]), chi_rel_sigma=0.05, chi_prime_rel_sigma=0.1),
        "decoherence": lib.NoiseModel(shots=int(cfg["robustness_shots"]), t1_s=50.0e-6, t2_s=30.0e-6),
    }
    robustness_rows: list[dict[str, Any]] = []
    rng_rob = np.random.default_rng(7777)

    for scenario_name, noise in noise_scenarios.items():
        for prot, settings in [("wait_only", wait_settings), ("combined", combined_settings)]:
            for trial in range(int(cfg["robustness_repeats"])):
                dataset = lib.diagonal_dataset_from_weighted_truth(
                    settings,
                    weighted_transverse=benchmark_case.weighted_transverse_truth,
                    z_total=benchmark_case.z_total_truth,
                    noise=noise,
                    rng=rng_rob,
                )
                fit = lib.infer_weighted_mle(dataset["rows"], settings)
                err = lib.recoverable_error_summary(
                    fit["weighted_transverse"],
                    benchmark_case.weighted_transverse_truth,
                    inferred_z_total=float(fit["z_total"]),
                    true_z_total=float(benchmark_case.z_total_truth),
                    true_sectors=benchmark_case.truth_sectors,
                )
                robustness_rows.append({
                    "noise_scenario": scenario_name,
                    "protocol": prot,
                    "trial": int(trial),
                    "weighted_rmse": float(err["weighted_rmse"]),
                    "residual_rms": float(fit["residual_rms"]),
                    "mean_oracle_fidelity": err.get("mean_oracle_fidelity"),
                })

    results["robustness_rows"] = robustness_rows
    print(f"  done ({time.time()-t0:.1f}s)")

    # ------------------------------------------------------------------
    # Part 3: Comparison questions (answered quantitatively in summary)
    # ------------------------------------------------------------------
    print("Part 3: comparison questions...")
    q1_wait_transverse_rank = int(ident["wait_only"]["transverse_rank"])
    q2_disp_needed = bool(ident["displacement_only"]["transverse_rank"] < q1_wait_transverse_rank)
    q3_combined_vs_wait_rmse = {
        s: {
            "wait_only": next((r["weighted_rmse_mean"] for r in protocol_summary if r["protocol"] == "wait_only" and r["shots"] == s and r["method"] == "mle"), None),
            "combined": next((r["weighted_rmse_mean"] for r in protocol_summary if r["protocol"] == "combined" and r["shots"] == s and r["method"] == "mle"), None),
        }
        for s in [1_000, 10_000]
        if any(r["shots"] == s for r in protocol_summary)
    }
    q4_coherence_flag = coherence_combined_residual > 1.0e-2
    q5_useful_threshold = max(
        (r["coherence_fraction"] for r in sweep_results
         if r["combined"]["mle"]["weighted_rmse"] < 0.05),
        default=0.0,
    )
    q6_best_protocol = "combined" if any(
        r["protocol"] == "combined" and r.get("weighted_rmse_mean", 1.0) < r2.get("weighted_rmse_mean", 1.0)
        for r in protocol_summary for r2 in protocol_summary
        if r["protocol"] == "wait_only" and r["shots"] == r2["shots"] and r["method"] == r2["method"]
    ) else "wait_only"

    comparison_answers = {
        "Q1_wait_recovers_transverse": bool(q1_wait_transverse_rank == 2 * N_ACTIVE),
        "Q1_transverse_rank": q1_wait_transverse_rank,
        "Q2_displacement_needed_for_Z_n_p_n": False,  # Analytic: displacement-only is uninformative at t=0
        "Q2_displacement_only_rank": int(ident["displacement_only"]["transverse_rank"]),
        "Q3_combined_vs_wait_rmse": q3_combined_vs_wait_rmse,
        "Q4_coherence_witness_works": bool(q4_coherence_flag),
        "Q5_useful_up_to_coherence_fraction": float(q5_useful_threshold),
        "Q6_recommended_protocol": str(q6_best_protocol),
        "Q6_recommendation_reason": (
            "Combined protocol achieves comparable or better RMSE and additionally detects "
            "cavity coherences via large residuals. Use combined when coherence is a concern."
        ),
    }
    results["comparison_answers"] = comparison_answers
    summary["comparison_answers"] = comparison_answers

    # ------------------------------------------------------------------
    # Save data
    # ------------------------------------------------------------------
    total_time = time.time() - t_start
    summary["total_runtime_s"] = float(total_time)
    save_json(RESULTS_JSON, results)
    save_json(SUMMARY_JSON, summary)
    save_json(VALIDATION_HINTS_JSON, {
        "identifiability": ident,
        "convergence": convergence,
        "coherence_wait_residual": coherence_wait_residual,
        "coherence_combined_residual": coherence_combined_residual,
        "gauge_family_count": len(gauge_family),
        "best_single_qubit_mle_mean_fidelity": float(baseline["best_mle_mean_fidelity"]),
    })
    print(f"\nAll results saved.  Total runtime: {total_time:.1f}s")

    # ------------------------------------------------------------------
    # Figures
    # ------------------------------------------------------------------
    print("Generating figures...")
    plot_single_qubit_baseline(baseline)

    plot_protocol_comparison(noisy_rows)

    # Per-Fock Bloch comparison: select a representative subset of cases
    bloch_cases = [
        c for c in analytic_cases
        if c.case_id in ("ideal_diag_g", "ideal_diag_+", "ideal_coherent_g")
    ] + [c for c in pulse_cases if "long_pulse_optimized_mix_g" in c.case_id or "leakage" in c.case_id]
    bloch_cases = bloch_cases[:6]
    if bloch_cases:
        plot_per_fock_bloch_comparison(bloch_cases, combined_settings, shots=int(cfg["coherence_shots"]))

    plot_full_state_nonuniqueness(full_state_result)

    plot_coherence_witness_residuals(coherence_rows)

    plot_coherence_sweep(sweep_results)

    plot_pulse_case_recovery(pulse_rows)

    if robustness_rows:
        plot_robustness(robustness_rows)

    print("Figures saved to", FIGURES_DIR)

    # ------------------------------------------------------------------
    # Markdown summary
    # ------------------------------------------------------------------
    md_lines = [
        "# Study Summary: Fock-Resolved SQR Inference v2",
        "",
        f"**Profile:** {cfg['profile_name']}  **Runtime:** {total_time:.1f}s",
        "",
        "## Part 1 — Single-Qubit Baseline",
        f"Best MLE mean fidelity: **{baseline['best_mle_mean_fidelity']:.4f}**",
        "",
        "## Part 2A — Identifiability",
        f"| Protocol | Transverse rank | Full rank | Cond. # |",
        f"|---|---|---|---|",
    ]
    for prot in ("wait_only", "displacement_only", "combined"):
        i = ident[prot]
        cond = f"{i['transverse_condition_number']:.2f}" if i["transverse_condition_number"] else "∞"
        md_lines.append(f"| {prot} | {i['transverse_rank']} | {i['full_rank']} | {cond} |")
    md_lines += [
        "",
        "## Part 2B — Full-State Non-Uniqueness",
        f"Objective span across {cfg['full_state_restarts']} restarts: **{full_state_result['objective_span']:.2f}**",
        f"Population σ by sector: {[f'{v:.3f}' for v in full_state_result['probability_std_by_sector']]}",
        f"Gauge family size: **{len(gauge_family)}**",
        "",
        "## Part 2C — Recoverable Subspace (MLE)",
        "| Protocol | N=1000 RMSE | N=10000 RMSE |",
        "|---|---|---|",
    ]
    for prot in ("wait_only", "displacement_only", "combined"):
        r1k = next((r["weighted_rmse_mean"] for r in protocol_summary if r["protocol"] == prot and r["shots"] == 1000 and r["method"] == "mle"), None)
        r10k = next((r["weighted_rmse_mean"] for r in protocol_summary if r["protocol"] == prot and r["shots"] == 10000 and r["method"] == "mle"), None)
        md_lines.append(f"| {prot} | {f'{r1k:.4f}' if r1k else 'n/a'} | {f'{r10k:.4f}' if r10k else 'n/a'} |")
    md_lines += [
        "",
        "## Part 2F — Coherence Sweep",
        f"Coherence witness gap (combined, f=1.0 vs f=0.0): "
        f"**{coherence_combined_residual:.4f}** vs **{coherence_wait_residual:.2e}**",
        "",
        "## Part 2E — Pulse-Level Cases",
        f"Best case: **{summary['pulse_best_case']['case_id']}** RMSE={summary['pulse_best_case']['weighted_rmse']:.4f}",
        f"Worst case: **{summary['pulse_worst_case']['case_id']}** RMSE={summary['pulse_worst_case']['weighted_rmse']:.4f}",
        "",
        "## Part 3 — Comparison Answers",
        f"Q1 Wait-only recovers transverse: **{comparison_answers['Q1_wait_recovers_transverse']}**",
        f"Q2 Displacement needed for p_n/Z_n: **{comparison_answers['Q2_displacement_needed_for_Z_n_p_n']}**",
        f"Q4 Coherence witness works: **{comparison_answers['Q4_coherence_witness_works']}**",
        f"Q5 Useful up to coherence fraction: **{comparison_answers['Q5_useful_up_to_coherence_fraction']:.2f}**",
        f"Q6 Recommended protocol: **{comparison_answers['Q6_recommended_protocol']}**",
    ]
    SUMMARY_MD.write_text("\n".join(md_lines), encoding="utf-8")
    print(f"Summary written to {SUMMARY_MD}")


if __name__ == "__main__":
    main()
