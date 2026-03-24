"""Plot results for the simultaneous multitone SQR study."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from common import CASE_LABELS, FIGURES_DIR


DATA_PATH = Path(__file__).resolve().parent.parent / "data" / "simultaneous_multitone_sqr_results.npz"
SUMMARY_PATH = Path(__file__).resolve().parent.parent / "data" / "simultaneous_multitone_sqr_summary.json"

TOL = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB"]


def configure_style() -> None:
    plt.style.use("default")
    plt.rcParams.update(
        {
            "figure.figsize": (10.5, 7.5),
            "font.size": 11,
            "axes.grid": False,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "savefig.bbox": "tight",
            "savefig.dpi": 220,
        }
    )


def load_payload():
    if not DATA_PATH.exists():
        raise FileNotFoundError(f"{DATA_PATH} not found. Run run_multitone_simultaneous_sqr_study.py first.")
    if not SUMMARY_PATH.exists():
        raise FileNotFoundError(f"{SUMMARY_PATH} not found. Run run_multitone_simultaneous_sqr_study.py first.")
    data = np.load(DATA_PATH, allow_pickle=True)
    summary = json.loads(SUMMARY_PATH.read_text(encoding="utf-8"))
    return data, summary


def _annotate_heatmap(ax, values: np.ndarray, *, fmt: str = "{:.3f}", threshold: float = 0.55) -> None:
    for row in range(values.shape[0]):
        for col in range(values.shape[1]):
            val = float(values[row, col])
            color = "white" if val < threshold else "black"
            ax.text(col, row, fmt.format(val), ha="center", va="center", color=color, fontsize=8)


def _save_figure(fig, stem: str) -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURES_DIR / f"{stem}.png")
    fig.savefig(FIGURES_DIR / f"{stem}.pdf")
    plt.close(fig)


def plot_heatmaps(data, *, field: str, title: str, stem: str, vmin: float, vmax: float, cmap: str) -> None:
    case_names = [str(x) for x in data["case_names"]]
    chi_t = np.asarray(data["chi_t_values"], dtype=float)
    theta = np.asarray(data["theta_over_pi"], dtype=float)
    values = np.asarray(data[field], dtype=float)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8.2), constrained_layout=True)
    axes = axes.reshape(-1)
    for ax, case_name, panel in zip(axes, case_names, values, strict=True):
        image = ax.imshow(panel, aspect="auto", origin="lower", cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_title(CASE_LABELS[case_name])
        ax.set_xticks(np.arange(theta.size))
        ax.set_xticklabels([f"{x:.3g}" for x in theta])
        ax.set_yticks(np.arange(chi_t.size))
        ax.set_yticklabels([f"{x:.0f}" if abs(x - round(x)) < 1e-9 else f"{x:.2g}" for x in chi_t])
        ax.set_xlabel(r"Target angle $\theta/\pi$")
        ax.set_ylabel(r"$\chi T / 2\pi$")
        _annotate_heatmap(ax, panel, threshold=0.45 if vmax <= 1.0 else 0.001)
    cbar = fig.colorbar(image, ax=axes.tolist(), shrink=0.9)
    cbar.set_label(title)
    _save_figure(fig, stem)


def plot_correction_checks(data, summary: dict[str, object]) -> None:
    amp_scales = np.asarray(data["amplitude_scan_scales"], dtype=float)
    amp_reduced = np.asarray(data["amplitude_scan_reduced_fidelity"], dtype=float)
    amp_target = np.asarray(data["amplitude_scan_target_mean"], dtype=float)
    amp_spectator = np.asarray(data["amplitude_scan_spectator_mean"], dtype=float)

    multistart = summary["multistart"]
    qutrit_rows = summary["qutrit_spotcheck"]["rows"]

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6), constrained_layout=True)

    ax = axes[0]
    ax.plot(amp_scales, amp_reduced, color=TOL[0], lw=2.2, label="Reduced mean fidelity")
    ax.plot(amp_scales, amp_target, color=TOL[1], lw=2.0, label="Target-branch mean")
    ax.plot(amp_scales, amp_spectator, color=TOL[2], lw=2.0, label="Spectator mean")
    ax.axvline(0.0, color="0.5", ls="--", lw=1.0)
    ax.set_xlabel(r"Common target-tone $\Delta \lambda$")
    ax.set_ylabel("Fidelity")
    ax.set_ylim(0.0, 1.02)
    ax.set_title("Simple amplitude correction is flat")
    ax.legend(frameon=False, fontsize=9)

    ax = axes[1]
    baseline = float(multistart["baseline_reduced_fidelity"])
    best = float(multistart["best_reduced_fidelity"])
    seed_rows = multistart["seed_rows"]
    seed_ids = [int(row["seed"]) for row in seed_rows]
    seed_values = [float(row["optimized_weighted_mean_fidelity"]) for row in seed_rows]
    ax.axhline(baseline, color=TOL[0], lw=2.0, label="Baseline")
    ax.bar(np.arange(len(seed_ids)), seed_values, color=TOL[1], alpha=0.85, label="Best per seed")
    ax.axhline(best, color=TOL[2], lw=2.0, ls="--", label="Best multistart")
    ax.set_xticks(np.arange(len(seed_ids)))
    ax.set_xticklabels([str(seed) for seed in seed_ids])
    ax.set_xlabel("Random seed")
    ax.set_ylabel("Reduced mean fidelity")
    ax.set_ylim(0.45, 0.55)
    ax.set_title("Multistart waveform correction")
    ax.legend(frameon=False, fontsize=9)

    ax = axes[2]
    labels = [
        f"{CASE_LABELS[str(row['case_name'])]}\n" + r"$\theta/\pi=$" + f"{float(row['theta_over_pi']):.3g}"
        for row in qutrit_rows
    ]
    values = np.asarray([float(row["max_f_leakage"]) for row in qutrit_rows], dtype=float)
    ax.bar(np.arange(len(labels)), values, color=TOL[3], alpha=0.9)
    ax.set_yscale("log")
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=18, ha="right")
    ax.set_ylabel(r"Max $|f\rangle$ leakage")
    ax.set_title("Qutrit replay spot-check")
    ax.set_ylim(1.0e-12, max(1.0e-6, 5.0 * float(np.max(values))))

    _save_figure(fig, "simultaneous_multitone_correction_checks")


def plot_segmented_vs_direct(summary: dict[str, object]) -> None:
    rows = summary["segmented_check"]["rows"]
    labels = [CASE_LABELS[str(row["case_name"])] for row in rows]
    direct = np.asarray([float(row["direct_pi_strict_fidelity"]) for row in rows], dtype=float)
    segmented = np.asarray([float(row["segmented_operator_fidelity"]) for row in rows], dtype=float)

    fig, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
    xpos = np.arange(len(labels), dtype=float)
    width = 0.34
    ax.bar(xpos - width / 2.0, direct, width=width, color=TOL[0], label=r"Direct $\pi$")
    ax.bar(xpos + width / 2.0, segmented, width=width, color=TOL[1], label=r"$8\times (\pi/8)$ compiled")
    ax.set_xticks(xpos)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Logical process fidelity")
    ax.set_ylim(0.0, max(0.3, 1.1 * float(np.max(np.concatenate([direct, segmented])))))
    ax.set_title("Naive repeated small-angle compilation is not enough")
    ax.legend(frameon=False)
    _save_figure(fig, "simultaneous_multitone_segmented_vs_direct")


def plot_target_angle_response(data) -> None:
    case_names = [str(x) for x in data["case_names"]]
    theta = np.asarray(data["theta_over_pi"], dtype=float)
    theta_sim = np.asarray(data["reduced_target_theta_mean"], dtype=float) / np.pi
    chi_t = np.asarray(data["chi_t_values"], dtype=float)
    chi_index = int(np.where(np.isclose(chi_t, 3.0))[0][0])

    fig, ax = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
    ax.plot(theta, theta, color="0.35", lw=1.6, ls="--", label="Ideal response")
    for index, case_name in enumerate(case_names):
        ax.plot(
            theta,
            theta_sim[index, chi_index, :],
            marker="o",
            lw=2.0,
            color=TOL[index],
            label=CASE_LABELS[case_name],
        )
    ax.set_xlabel(r"Requested target angle $\theta/\pi$")
    ax.set_ylabel(r"Simulated target angle $\bar{\theta}_{\mathrm{sim}}/\pi$")
    ax.set_title(r"At $\chi T/2\pi = 3$, the common waveform barely rotates the targets")
    ax.legend(frameon=False, fontsize=9)
    ax.set_ylim(-0.01, 0.03)
    _save_figure(fig, "simultaneous_multitone_target_angle_response")


def main() -> None:
    configure_style()
    data, summary = load_payload()

    plot_heatmaps(
        data,
        field="reduced_fidelity",
        title="Reduced conditioned fidelity",
        stem="simultaneous_multitone_conditioned_heatmaps",
        vmin=0.0,
        vmax=1.0,
        cmap="viridis",
    )
    plot_heatmaps(
        data,
        field="strict_fidelity",
        title="Strict logical fidelity",
        stem="simultaneous_multitone_strict_heatmaps",
        vmin=0.0,
        vmax=1.0,
        cmap="magma",
    )
    plot_heatmaps(
        data,
        field="compiled_gain",
        title="Best-fit cavity compilation gain",
        stem="simultaneous_multitone_compiled_gain_heatmaps",
        vmin=0.0,
        vmax=0.01,
        cmap="plasma",
    )
    plot_correction_checks(data, summary)
    plot_target_angle_response(data)
    plot_segmented_vs_direct(summary)
    print(f"Saved figures to {FIGURES_DIR}")


if __name__ == "__main__":
    main()
