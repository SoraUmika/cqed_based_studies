"""Plotting helpers for the nonlinear-QND readout follow-up study."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from .config import FIG_DIR

STYLE_PATH = Path(__file__).resolve().parents[4] / ".github" / "skills" / "publication-figures" / "assets" / "cqed_style.mplstyle"
plt.style.use(str(STYLE_PATH))

PALETTE = {
    "square": "#4477AA",
    "smooth_square": "#66CCEE",
    "ring_hold": "#CCBB44",
    "procedural_segments": "#228833",
    "nulling_tail": "#EE6677",
    "fourier_basis": "#AA3377",
    "piecewise_reference": "#BBBBBB",
}


def _load_summary(summary_path: Path) -> dict:
    return json.loads(summary_path.read_text())


def _duration_array(summary: dict) -> np.ndarray:
    return np.array([row["duration_ns"] for row in summary["hierarchy"]], dtype=float)


def _metric_series(summary: dict, bucket: str, family: str, metric: str) -> np.ndarray:
    records = summary[bucket][family]
    if bucket in {"full_results", "rich_results"}:
        return np.array([entry["evaluation"]["metrics"][metric] for entry in records], dtype=float)
    if bucket == "nominal_rich_replay":
        return np.array([entry["metrics"][metric] for entry in records], dtype=float)
    raise ValueError(f"Unsupported bucket '{bucket}'.")


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
    legacy = np.array([row["legacy_best"] for row in summary["hierarchy"]], dtype=float)
    nominal_rich = np.array([row["nominal_rich_best"] for row in summary["hierarchy"]], dtype=float)
    rich = np.array([row["rich_best"] for row in summary["hierarchy"]], dtype=float)
    qnd = np.array([row["rich_qnd_constrained_best"] for row in summary["hierarchy"]], dtype=float)

    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    ax.plot(durations, ideal, label="Ideal linear reference", color="#222222", linewidth=2.4)
    ax.plot(durations, detector, label="Detector-limited", color="#4477AA", linewidth=2.0)
    ax.plot(durations, t1_bound, label="$T_1$ reference", color="#CCBB44", linewidth=2.0)
    ax.plot(durations, legacy, label="Legacy full-model best", color="#66CCEE", linewidth=2.0, marker="o", markersize=4.0)
    ax.plot(durations, nominal_rich, label="Nominal pulses replayed in rich model", color="#EE6677", linewidth=2.0, marker="s", markersize=4.0)
    ax.plot(durations, rich, label="Rich-model re-optimized best", color="#228833", linewidth=2.3, marker="D", markersize=4.0)
    ax.plot(durations, qnd, label="QND-constrained practical frontier", color="#AA3377", linewidth=1.8, linestyle="--")
    ax.set_xlabel("Readout duration (ns)")
    ax.set_ylabel("Assignment fidelity")
    ax.set_ylim(0.45, 1.01)
    ax.legend(frameon=False, fontsize=8)
    _save_figure(fig, "fig1_bound_hierarchy", figure_dir)


def plot_rich_frontiers(summary: dict, figure_dir: Path) -> None:
    durations = _duration_array(summary)
    families = ["square", "procedural_segments", "nulling_tail", "fourier_basis"]
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))
    metric_triplet = (
        ("fidelity_eta", "Detector-limited fidelity"),
        ("qnd_preservation", "QND preservation"),
        ("measurement_induced_transition", "Induced transition probability"),
    )
    for axis, (metric, ylabel) in zip(axes, metric_triplet, strict=True):
        for family in families:
            axis.plot(
                durations,
                _metric_series(summary, "rich_results", family, metric),
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
    axes[0].legend(frameon=False, fontsize=8, loc="lower right")
    _save_figure(fig, "fig2_rich_frontiers", figure_dir)


def plot_regime_breakdown(summary: dict, figure_dir: Path) -> None:
    families = ["square", "procedural_segments", "nulling_tail", "fourier_basis"]
    regimes = ["full", "hardware", "nonlinear", "rich"]
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.2))
    for family in families:
        fidelity = [summary["representative_breakdown"][family][regime]["metrics"]["fidelity_eta"] for regime in regimes]
        qnd_defect = [summary["representative_breakdown"][family][regime]["metrics"]["measurement_induced_qnd_defect"] for regime in regimes]
        axes[0].plot(regimes, fidelity, marker="o", linewidth=1.8, color=PALETTE[family], label=family.replace("_", " "))
        axes[1].plot(regimes, qnd_defect, marker="o", linewidth=1.8, color=PALETTE[family], label=family.replace("_", " "))
    axes[0].set_ylabel("Assignment fidelity")
    axes[1].set_ylabel("Excess QND defect")
    for axis in axes:
        axis.set_xlabel("Replay regime")
        axis.tick_params(axis="x", rotation=20)
    axes[0].legend(frameon=False, fontsize=8)
    _save_figure(fig, "fig3_regime_breakdown", figure_dir)


def plot_transport_waveforms(summary: dict, traces_path: Path, figure_dir: Path) -> None:
    traces = np.load(traces_path)
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 5.8), sharex="col")
    for col, family in enumerate(("procedural_segments", "nulling_tail")):
        program = traces[f"{family}_program_waveform"]
        distorted = traces[f"{family}_distorted_waveform"]
        t_ns = np.arange(program.size, dtype=float) * summary["config"]["dt"] * 1.0e9
        axes[0, col].plot(t_ns, np.abs(program), color=PALETTE[family], linewidth=1.8, label="Programmed")
        axes[0, col].plot(t_ns, np.abs(distorted), color="#222222", linewidth=1.4, linestyle="--", label="Distorted")
        axes[0, col].set_ylabel("|$\\epsilon$| (rad/s)")
        axes[0, col].set_title(family.replace("_", " "))
        axes[1, col].plot(t_ns, np.unwrap(np.angle(program)), color=PALETTE[family], linewidth=1.8)
        axes[1, col].plot(t_ns, np.unwrap(np.angle(distorted)), color="#222222", linewidth=1.4, linestyle="--")
        axes[1, col].set_ylabel("Phase (rad)")
        axes[1, col].set_xlabel("Time (ns)")
    axes[0, 0].legend(frameon=False, fontsize=8)
    _save_figure(fig, "fig4_transport_waveforms", figure_dir)


def plot_tradeoffs(summary: dict, figure_dir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.0))
    for family in ("square", "procedural_segments", "nulling_tail", "fourier_basis", "piecewise_reference"):
        entries = summary["tradeoff_slice"][family]
        fidelities = [entries[name]["evaluation"]["metrics"]["fidelity_eta"] for name in ("info", "balanced", "emptying")]
        residuals = [entries[name]["evaluation"]["metrics"]["residual_photons"] for name in ("info", "balanced", "emptying")]
        induced = [entries[name]["evaluation"]["metrics"]["measurement_induced_transition"] for name in ("info", "balanced", "emptying")]
        axes[0].plot(residuals, fidelities, marker="o", linewidth=1.8, color=PALETTE[family], label=family.replace("_", " "))
        axes[1].plot(induced, fidelities, marker="o", linewidth=1.8, color=PALETTE[family], label=family.replace("_", " "))
    axes[0].set_xlabel("Residual photons")
    axes[0].set_ylabel("Assignment fidelity")
    axes[1].set_xlabel("Induced transition probability")
    axes[1].set_ylabel("Assignment fidelity")
    axes[0].legend(frameon=False, fontsize=8, loc="lower right")
    _save_figure(fig, "fig5_tradeoffs", figure_dir)


def plot_qnd_stress(summary: dict, figure_dir: Path) -> None:
    families = ["square", "procedural_segments", "nulling_tail", "fourier_basis"]
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.0))
    for family in families:
        records = summary["qnd_stress"][family]
        scales = np.array([row["amp_scale"] for row in records], dtype=float)
        qnd = np.array([row["qnd_preservation"] for row in records], dtype=float)
        induced = np.array([row["measurement_induced_transition"] for row in records], dtype=float)
        axes[0].plot(scales, qnd, marker="o", linewidth=1.8, color=PALETTE[family], label=family.replace("_", " "))
        axes[1].plot(scales, induced, marker="o", linewidth=1.8, color=PALETTE[family], label=family.replace("_", " "))
    axes[0].set_xlabel("Amplitude scale")
    axes[0].set_ylabel("QND preservation")
    axes[1].set_xlabel("Amplitude scale")
    axes[1].set_ylabel("Induced transition probability")
    axes[0].legend(frameon=False, fontsize=8)
    _save_figure(fig, "fig6_qnd_stress", figure_dir)


def plot_robustness(summary: dict, figure_dir: Path) -> None:
    families = ["square", "procedural_segments", "nulling_tail", "fourier_basis"]
    nominal = np.array([summary["robustness"][family]["nominal_fidelity_eta"] for family in families], dtype=float)
    mean = np.array([summary["robustness"][family]["mean_fidelity_eta"] for family in families], dtype=float)
    worst = np.array([summary["robustness"][family]["worst_fidelity_eta"] for family in families], dtype=float)
    x = np.arange(len(families))

    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    ax.bar(x - 0.22, nominal, width=0.22, label="Nominal", color="#4477AA")
    ax.bar(x, mean, width=0.22, label="Mean perturbed", color="#66CCEE")
    ax.bar(x + 0.22, worst, width=0.22, label="Worst perturbed", color="#EE6677")
    ax.set_xticks(x)
    ax.set_xticklabels([family.replace("_", "\n") for family in families])
    ax.set_ylabel("Assignment fidelity")
    ax.set_ylim(0.45, 1.01)
    ax.legend(frameon=False, fontsize=8)
    _save_figure(fig, "fig7_robustness", figure_dir)


def generate_all_figures(
    *,
    summary_path: Path,
    traces_path: Path,
    figure_dir: Path = FIG_DIR,
) -> None:
    summary = _load_summary(summary_path)
    plot_bound_hierarchy(summary, figure_dir)
    plot_rich_frontiers(summary, figure_dir)
    plot_regime_breakdown(summary, figure_dir)
    plot_transport_waveforms(summary, traces_path, figure_dir)
    plot_tradeoffs(summary, figure_dir)
    plot_qnd_stress(summary, figure_dir)
    plot_robustness(summary, figure_dir)
