"""Plotting helpers for the procedural readout study."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .config import FIG_DIR

PALETTE = {
    "square": "#1b4965",
    "smooth_square": "#5fa8d3",
    "ring_hold": "#ca6702",
    "procedural_segments": "#2a9d8f",
    "nulling_tail": "#e76f51",
    "fourier_basis": "#7b2cbf",
    "piecewise_reference": "#111111",
}


def _load_summary(summary_path: Path) -> dict:
    return json.loads(summary_path.read_text())


def _duration_array(summary: dict) -> np.ndarray:
    return np.array([row["duration_ns"] for row in summary["hierarchy"]], dtype=float)


def _metric_series(summary: dict, family: str, metric: str) -> np.ndarray:
    return np.array(
        [entry["evaluation"]["metrics"][metric] for entry in summary["full_results"][family]],
        dtype=float,
    )


def _save_figure(fig: plt.Figure, stem: str, figure_dir: Path) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figure_dir / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(figure_dir / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_bound_hierarchy(summary: dict, figure_dir: Path) -> None:
    durations = _duration_array(summary)
    ideal = np.array([row["ideal_bound"] for row in summary["hierarchy"]], dtype=float)
    detector = np.array([row["detector_bound"] for row in summary["hierarchy"]], dtype=float)
    t1_bound = np.array([row["t1_bound"] for row in summary["hierarchy"]], dtype=float)
    realistic = np.array([row["realistic_best"] for row in summary["hierarchy"]], dtype=float)

    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    ax.plot(durations, ideal, label="Ideal linear bound", color="#111111", linewidth=2.4)
    ax.plot(durations, detector, label="Detector-limited", color="#264653", linewidth=2.0)
    ax.plot(durations, t1_bound, label="$T_1$-limited", color="#f4a261", linewidth=2.0)
    ax.plot(durations, realistic, label="Realistic balanced replay", color="#2a9d8f", linewidth=2.2)
    ax.set_xlabel("Readout duration (ns)")
    ax.set_ylabel("Assignment fidelity")
    ax.set_ylim(0.45, 1.01)
    ax.legend(frameon=False)
    ax.set_title("Bound hierarchy vs readout duration")
    _save_figure(fig, "fig1_bound_hierarchy", figure_dir)


def plot_frontiers(summary: dict, figure_dir: Path) -> None:
    durations = _duration_array(summary)
    families = ["square", "procedural_segments", "nulling_tail", "fourier_basis"]
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.0))
    metric_triplet = (
        ("fidelity_eta", "Detector-limited fidelity"),
        ("qnd_preservation", "QND preservation"),
        ("residual_photons", "Residual photons"),
    )
    for axis, (metric, ylabel) in zip(axes, metric_triplet, strict=True):
        for family in families:
            axis.plot(
                durations,
                _metric_series(summary, family, metric),
                label=family.replace("_", " "),
                color=PALETTE[family],
                linewidth=2.0,
                marker="o",
                markersize=4.0,
            )
        axis.set_xlabel("Readout duration (ns)")
        axis.set_ylabel(ylabel)
    axes[0].set_ylim(0.45, 1.01)
    axes[1].set_ylim(0.90, 1.01)
    axes[0].legend(frameon=False, loc="lower right")
    fig.suptitle("Full-model performance frontiers")
    _save_figure(fig, "fig2_frontiers", figure_dir)


def plot_representative_diagnostics(summary: dict, traces_path: Path, figure_dir: Path) -> None:
    traces = np.load(traces_path)
    family = "procedural_segments"
    t_ns = traces[f"{family}_t_ns"]
    alpha_g = traces[f"{family}_alpha_g"]
    alpha_e = traces[f"{family}_alpha_e"]
    signal_g = traces[f"{family}_signal_g"]
    signal_e = traces[f"{family}_signal_e"]
    metrics = summary["representative"][family]["metrics"]
    mu = float(np.sqrt(max(metrics["snr2_ideal"] * summary["config"]["eta_nominal"], 1.0e-12)))

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.0))
    axes[0].plot(t_ns[:-1], np.abs(signal_g[:-1]), color=PALETTE["square"], label="|s_g|")
    axes[0].plot(t_ns[:-1], np.abs(signal_e[:-1]), color=PALETTE["procedural_segments"], label="|s_e|")
    axes[0].set_xlabel("Time (ns)")
    axes[0].set_ylabel("Output amplitude (arb.)")
    axes[0].legend(frameon=False)

    axes[1].plot(alpha_g.real, alpha_g.imag, color=PALETTE["square"], label=r"$\alpha_g$")
    axes[1].plot(alpha_e.real, alpha_e.imag, color=PALETTE["procedural_segments"], label=r"$\alpha_e$")
    axes[1].scatter([alpha_g[-1].real, alpha_e[-1].real], [alpha_g[-1].imag, alpha_e[-1].imag], s=18, color="#111111")
    axes[1].set_xlabel("Re($\\alpha$)")
    axes[1].set_ylabel("Im($\\alpha$)")
    axes[1].legend(frameon=False)

    x = np.linspace(-4.0, 4.0, 400)
    gaussian = lambda center: np.exp(-0.5 * (x - center) ** 2) / np.sqrt(2.0 * np.pi)
    axes[2].plot(x, gaussian(-mu), color=PALETTE["square"], label="Ground trace")
    axes[2].plot(x, gaussian(mu), color=PALETTE["procedural_segments"], label="Excited trace")
    axes[2].axvline(0.0, color="#444444", linestyle="--", linewidth=1.0)
    axes[2].set_xlabel("Matched-filter score (normalized)")
    axes[2].set_ylabel("Probability density")
    axes[2].legend(frameon=False)

    fig.suptitle("Representative procedural-pulse diagnostics")
    _save_figure(fig, "fig3_representative_diagnostics", figure_dir)


def plot_tradeoffs(summary: dict, figure_dir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0))
    for family in ("square", "procedural_segments", "nulling_tail", "fourier_basis"):
        entries = summary["tradeoff_slice"][family]
        fidelities = [entries[name]["evaluation"]["metrics"]["fidelity_eta"] for name in ("info", "balanced", "emptying")]
        residuals = [entries[name]["evaluation"]["metrics"]["residual_photons"] for name in ("info", "balanced", "emptying")]
        qnds = [entries[name]["evaluation"]["metrics"]["qnd_preservation"] for name in ("info", "balanced", "emptying")]
        axes[0].plot(residuals, fidelities, marker="o", linewidth=1.8, color=PALETTE[family], label=family.replace("_", " "))
        axes[1].plot(qnds, fidelities, marker="o", linewidth=1.8, color=PALETTE[family], label=family.replace("_", " "))
    axes[0].set_xlabel("Residual photons")
    axes[0].set_ylabel("Detector-limited fidelity")
    axes[1].set_xlabel("QND preservation")
    axes[1].set_ylabel("Detector-limited fidelity")
    axes[0].legend(frameon=False, loc="lower right")
    fig.suptitle("Objective-tradeoff slice at the representative duration")
    _save_figure(fig, "fig4_tradeoffs", figure_dir)


def plot_robustness(summary: dict, figure_dir: Path) -> None:
    families = ["square", "procedural_segments", "nulling_tail", "fourier_basis"]
    nominal = np.array([summary["robustness"][family]["nominal_fidelity_eta"] for family in families], dtype=float)
    mean = np.array([summary["robustness"][family]["mean_fidelity_eta"] for family in families], dtype=float)
    worst = np.array([summary["robustness"][family]["worst_fidelity_eta"] for family in families], dtype=float)
    x = np.arange(len(families))

    fig, ax = plt.subplots(figsize=(8.0, 4.2))
    ax.bar(x - 0.22, nominal, width=0.22, label="Nominal", color="#1b4965")
    ax.bar(x, mean, width=0.22, label="Mean perturbed", color="#5fa8d3")
    ax.bar(x + 0.22, worst, width=0.22, label="Worst perturbed", color="#ca6702")
    ax.set_xticks(x)
    ax.set_xticklabels([family.replace("_", "\n") for family in families])
    ax.set_ylabel("Detector-limited fidelity")
    ax.set_ylim(0.45, 1.01)
    ax.legend(frameon=False)
    ax.set_title("Robustness to parameter and calibration errors")
    _save_figure(fig, "fig5_robustness", figure_dir)


def generate_all_figures(
    *,
    summary_path: Path,
    traces_path: Path,
    figure_dir: Path = FIG_DIR,
) -> None:
    summary = _load_summary(summary_path)
    plot_bound_hierarchy(summary, figure_dir)
    plot_frontiers(summary, figure_dir)
    plot_representative_diagnostics(summary, traces_path, figure_dir)
    plot_tradeoffs(summary, figure_dir)
    plot_robustness(summary, figure_dir)
