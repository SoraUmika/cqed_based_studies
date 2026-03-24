"""Generate report-level summary figures for the consolidated readout study.

The underlying numbers come from the validated consolidated report tables and
summary paragraphs. These figures are meant to make the report easier to read
as a standalone experimental guide; they do not replace the archived study data.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


FIG_DIR = Path(__file__).resolve().parents[1] / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

COLORS = {
    "Square": "#4E79A7",
    "Ring-hold": "#F28E2B",
    "Procedural segments": "#59A14F",
    "Nulled tail": "#E15759",
    "Fourier basis": "#76B7B2",
    "Piecewise reference": "#B07AA1",
}

LEVEL1 = [
    {"family": "Square", "fidelity": 0.871, "residual": 1.485, "qnd": 0.853},
    {"family": "Ring-hold", "fidelity": 0.947, "residual": 0.554, "qnd": 0.928},
    {
        "family": "Procedural segments",
        "fidelity": 0.956,
        "residual": 0.620,
        "qnd": 0.950,
    },
    {"family": "Nulled tail", "fidelity": 0.944, "residual": 0.265, "qnd": 0.936},
    {"family": "Fourier basis", "fidelity": 0.953, "residual": 0.721, "qnd": 0.938},
    {
        "family": "Piecewise reference",
        "fidelity": 0.985,
        "residual": 5.572,
        "qnd": 0.979,
    },
]

LEVEL2 = [
    {"family": "Square", "fidelity": 0.9201, "residual": 2.366, "qnd": 0.9931},
    {"family": "Ring-hold", "fidelity": 0.9599, "residual": 0.801, "qnd": 0.9921},
    {
        "family": "Procedural segments",
        "fidelity": 0.9547,
        "residual": 0.713,
        "qnd": 0.9886,
    },
    {"family": "Nulled tail", "fidelity": 0.9216, "residual": 0.887, "qnd": 0.9931},
    {"family": "Fourier basis", "fidelity": 0.9347, "residual": 0.673, "qnd": 0.9931},
    {
        "family": "Piecewise reference",
        "fidelity": 0.9092,
        "residual": 0.416,
        "qnd": 0.9933,
    },
]

LEAKAGE_THRESHOLD = {
    "duration_ns": np.array([120.0, 240.0, 480.0]),
    "amplitude_mhz": np.array([7.5, 5.5, 4.5]),
}

QND_STRESS = {
    "Procedural segments": 0.9618,
    "Fourier basis": 0.9890,
}

EMPTYING_COMPARISON = {
    "Square": 0.216,
    "Gaussian": 0.053,
}


def _configure_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 160,
            "savefig.dpi": 300,
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "legend.fontsize": 9,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
        }
    )


def _plot_tradeoff_panel(ax: plt.Axes, records: list[dict], title: str) -> None:
    for item in records:
        size = 2200.0 * max(item["qnd"] - 0.80, 0.03)
        ax.scatter(
            item["residual"],
            item["fidelity"],
            s=size,
            color=COLORS[item["family"]],
            alpha=0.88,
            edgecolor="black",
            linewidth=0.6,
        )
        ax.annotate(
            item["family"],
            (item["residual"], item["fidelity"]),
            textcoords="offset points",
            xytext=(6, 6),
            fontsize=8,
        )

    ax.set_xscale("log")
    ax.set_xlabel("Residual cavity photons")
    ax.set_ylabel("Detector-limited fidelity")
    ax.set_title(title)
    ax.grid(alpha=0.25, which="both")
    ax.set_ylim(0.86, 0.99 if "Level 1" in title else 0.965)

    ax.text(
        0.03,
        0.05,
        "Larger markers = better repeated-readout / QND preservation",
        transform=ax.transAxes,
        fontsize=8,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85},
    )


def make_tradeoff_figure() -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8), constrained_layout=True)
    _plot_tradeoff_panel(axes[0], LEVEL1, "Level 1: multilevel procedural replay at 240 ns")
    _plot_tradeoff_panel(axes[1], LEVEL2, "Level 2: hardware-realistic replay at 240 ns")
    fig.suptitle(
        "Structured readout pulses move the experiment toward the top-left frontier",
        fontsize=13,
    )
    for ext in ("pdf", "png"):
        fig.savefig(FIG_DIR / f"readout_family_tradeoffs.{ext}", bbox_inches="tight")
    plt.close(fig)


def make_operating_window_figure() -> None:
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.4), constrained_layout=True)

    ax = axes[0]
    ax.plot(
        LEAKAGE_THRESHOLD["duration_ns"],
        LEAKAGE_THRESHOLD["amplitude_mhz"],
        marker="o",
        color="#E15759",
        linewidth=2.2,
    )
    ax.fill_between(
        LEAKAGE_THRESHOLD["duration_ns"],
        LEAKAGE_THRESHOLD["amplitude_mhz"],
        8.0,
        color="#E15759",
        alpha=0.15,
        label="Higher-risk region",
    )
    ax.set_xlabel("Pulse duration (ns)")
    ax.set_ylabel(r"Approx. $P_\mathrm{leak}=10^{-4}$ threshold amplitude (MHz)")
    ax.set_title("High-power boundary moves to lower drive for longer pulses")
    ax.grid(alpha=0.25)
    ax.legend(loc="upper right")

    ax = axes[1]
    names = list(QND_STRESS)
    values = [QND_STRESS[name] for name in names]
    ax.bar(names, values, color=[COLORS[name] for name in names])
    ax.set_ylim(0.94, 1.0)
    ax.set_ylabel(r"Repeated-readout preservation $Q_\mathrm{QND}$")
    ax.set_title("At 1.3x nominal amplitude, smoother bases are safer")
    for idx, value in enumerate(values):
        ax.text(idx, value + 0.001, f"{value:.4f}", ha="center", va="bottom", fontsize=8)
    ax.tick_params(axis="x", rotation=15)

    ax = axes[2]
    names = list(EMPTYING_COMPARISON)
    values = [EMPTYING_COMPARISON[name] for name in names]
    ax.bar(names, values, color=[COLORS.get(name, "#4E79A7") for name in names])
    ax.set_ylabel("Residual photons at matched information level")
    ax.set_title("Smoother envelopes mainly buy faster cavity emptying")
    for idx, value in enumerate(values):
        ax.text(idx, value + 0.005, f"{value:.3f}", ha="center", va="bottom", fontsize=8)

    fig.suptitle(
        "Experiment-facing readout window: ring-hold for balance, smoother pulses near limits",
        fontsize=13,
    )
    for ext in ("pdf", "png"):
        fig.savefig(FIG_DIR / f"readout_operating_window.{ext}", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    _configure_style()
    make_tradeoff_figure()
    make_operating_window_figure()


if __name__ == "__main__":
    main()
