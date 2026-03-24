"""Generate summary figures for the consolidated hybrid-control report."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


FIG_DIR = Path(__file__).resolve().parents[1] / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

LIBRARY_COLORS = {
    "A": "#4E79A7",
    "B": "#F28E2B",
    "C": "#59A14F",
    "D": "#E15759",
    "E": "#76B7B2",
    "F": "#B07AA1",
}

LOCAL_BENCHMARK = [
    {"library": "A", "label": "D-SNAP-D", "duration": 1260, "fidelity": 0.9887, "leakage": 0.0185},
    {"library": "E", "label": "GRAPE", "duration": 320, "fidelity": 0.9618, "leakage": 0.0524},
    {"library": "F", "label": "Native SWAP", "duration": 440, "fidelity": 0.8784, "leakage": 0.0788},
    {"library": "B", "label": "Selective SQR", "duration": 2440, "fidelity": 0.8615, "leakage": 0.1350},
    {"library": "D", "label": "chi-wait probe", "duration": 336, "fidelity": 0.5542, "leakage": 0.4747},
    {"library": "C", "label": "ECD-like", "duration": 440, "fidelity": 0.4996, "leakage": 0.4535},
]

ENTANGLER_BENCHMARK = [
    {"library": "D", "label": "chi-wait + Rq", "duration": 256, "strict": 1.0000, "block": 1.0000, "leakage": 0.0000},
    {"library": "A", "label": "Baseline A", "duration": 256, "strict": 1.0000, "block": 1.0000, "leakage": 0.0000},
    {"library": "E", "label": "GRAPE", "duration": 400, "strict": 0.9458, "block": 0.9467, "leakage": 0.0854},
    {"library": "F", "label": "Native SWAP", "duration": 320, "strict": 0.7510, "block": 0.7510, "leakage": 0.1834},
    {"library": "B", "label": "Single SQR", "duration": 1100, "strict": 0.7071, "block": 1.0000, "leakage": 0.0000},
    {"library": "C", "label": "ECD-like", "duration": 440, "strict": 0.5000, "block": 0.5077, "leakage": 0.0000},
]

UTARGET_DEPTH = np.array([7, 9, 11], dtype=float)
UTARGET_FIDELITY = {
    "Library 1: D + Rq + SQR": np.array([0.8846, 0.8971, 0.9193]),
    "Library 2: D + Rq + ConditionalPhaseSQR": np.array([0.6726, 0.6858, 0.8741]),
}
UTARGET_LEAKAGE = {
    "Library 1: D + Rq + SQR": np.array([0.1226, 0.1225, 0.1014]),
    "Library 2: D + Rq + ConditionalPhaseSQR": np.array([0.0921, 0.3523, 0.1236]),
}

ROBUSTNESS = {
    "chi error": 0.0000,
    "Amplitude": 0.0058,
    "Duration": 0.0000,
    "Phase offset": 0.0012,
}


def _configure_style() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 160,
            "savefig.dpi": 300,
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
        }
    )


def _scatter_panel(ax: plt.Axes, records: list[dict], y_key: str, title: str) -> None:
    for record in records:
        size = 900.0 * max(0.02, 1.0 - record["leakage"])
        ax.scatter(
            record["duration"],
            record[y_key],
            s=size,
            color=LIBRARY_COLORS[record["library"]],
            edgecolor="black",
            linewidth=0.6,
            alpha=0.88,
        )
        ax.annotate(
            f"{record['library']}: {record['label']}",
            (record["duration"], record[y_key]),
            textcoords="offset points",
            xytext=(6, 6),
            fontsize=8,
        )
    ax.set_xlabel("Gate duration (ns)")
    ax.set_ylabel("Strict logical fidelity")
    ax.set_title(title)
    ax.grid(alpha=0.25)


def make_library_summary() -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.6, 4.8), constrained_layout=True)
    _scatter_panel(
        axes[0],
        LOCAL_BENCHMARK,
        "fidelity",
        "Local cavity control: fidelity-duration-leakage tradeoff",
    )
    axes[0].set_ylim(0.45, 1.02)

    _scatter_panel(
        axes[1],
        ENTANGLER_BENCHMARK,
        "strict",
        "Entangler benchmark: native chi-wait is the practical winner",
    )
    axes[1].set_ylim(0.45, 1.02)

    fig.suptitle(
        "Best practical universal-control route is mixed: selective local control plus native entangling",
        fontsize=13,
    )
    for ext in ("pdf", "png"):
        fig.savefig(FIG_DIR / f"hybrid_gate_library_summary.{ext}", bbox_inches="tight")
    plt.close(fig)


def make_utarget_summary() -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.8), constrained_layout=True)

    ax = axes[0]
    for label, values in UTARGET_FIDELITY.items():
        color = "#4E79A7" if "Library 1" in label else "#F28E2B"
        ax.plot(UTARGET_DEPTH, values, marker="o", linewidth=2.2, label=label, color=color)
    ax.set_xlabel("Ansatz depth")
    ax.set_ylabel(r"Projected logical fidelity $F_\mathrm{proj}$")
    ax.set_ylim(0.65, 0.94)
    ax.set_title(r"$U_\mathrm{target}$ fidelity improves with depth, but leakage stays relevant")
    ax.grid(alpha=0.25)

    ax2 = ax.twinx()
    ax2.plot(
        UTARGET_DEPTH,
        UTARGET_LEAKAGE["Library 1: D + Rq + SQR"],
        marker="s",
        linestyle="--",
        linewidth=1.8,
        color="#59A14F",
        label="Library 1 leakage",
    )
    ax2.plot(
        UTARGET_DEPTH,
        UTARGET_LEAKAGE["Library 2: D + Rq + ConditionalPhaseSQR"],
        marker="^",
        linestyle="--",
        linewidth=1.8,
        color="#E15759",
        label="Library 2 leakage",
    )
    ax2.set_ylabel(r"Average leakage $L_\mathrm{avg}$")
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc="lower right")

    ax = axes[1]
    labels = list(ROBUSTNESS)
    values = [ROBUSTNESS[label] for label in labels]
    ax.bar(labels, values, color=["#76B7B2", "#E15759", "#59A14F", "#B07AA1"])
    ax.set_ylabel("RMS fidelity drop")
    ax.set_ylim(0.0, 0.007)
    ax.set_title("The best physical sequence is least sensitive to chi and timing drift")
    for idx, value in enumerate(values):
        ax.text(idx, value + 0.00015, f"{value:.4f}", ha="center", va="bottom", fontsize=8)
    ax.tick_params(axis="x", rotation=15)

    fig.suptitle(
        "Experiment-facing synthesis message: drift is manageable, displacement leakage is not",
        fontsize=13,
    )
    for ext in ("pdf", "png"):
        fig.savefig(FIG_DIR / f"hybrid_utarget_summary.{ext}", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    _configure_style()
    make_library_summary()
    make_utarget_summary()


if __name__ == "__main__":
    main()
