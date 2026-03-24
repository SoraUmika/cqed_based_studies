"""Generate summary figures for the consolidated SQR report.

These figures condense the validated numerical conclusions already described in
the report into experiment-facing visuals.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


FIG_DIR = Path(__file__).resolve().parents[1] / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

MULTIBRANCH_TARGETS = [r"{0,1}", r"{0,1,2}", r"{0,1,2,3}"]
MULTIBRANCH_ANGLES = [r"$\pi/8$", r"$\pi/4$", r"$\pi/2$", r"$\pi$"]
MULTIBRANCH_FIDELITY = np.array(
    [
        [0.9759, 0.9205, 0.7242, 0.2477],
        [0.9664, 0.8843, 0.6051, 0.0616],
        [0.9570, 0.8490, 0.4970, 3.8e-7],
    ]
)

PHASE_COMPILATION = {
    "Gaussian": (0.927, 0.949),
    "Cosine-squared": (0.931, 0.954),
}

FAMILY_MEAN_NOISY = {
    "Square": 0.898,
    "Cosine-squared": 0.891,
    "Gaussian": 0.877,
}

GRAPE_REPLAY = {
    "chi_t": np.array([1, 2, 3, 5], dtype=float),
    "grape": np.array([0.978, 0.958, 0.941, 0.907]),
    "square": np.array([0.944, 0.924, 0.898, 0.850]),
}

CONCURRENT_READOUT = {
    "amplitude_mhz": np.array([0.5, 1.0, 2.5]),
    "fidelity": np.array([0.822, 0.726, 0.632]),
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


def make_control_summary() -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.6), constrained_layout=True)

    ax = axes[0]
    labels = list(PHASE_COMPILATION)
    before = [PHASE_COMPILATION[label][0] for label in labels]
    after = [PHASE_COMPILATION[label][1] for label in labels]
    x = np.arange(len(labels))
    width = 0.34
    ax.bar(x - width / 2, before, width=width, label="Before cavity-phase layer", color="#4E79A7")
    ax.bar(x + width / 2, after, width=width, label="After cavity-phase layer", color="#F28E2B")
    ax.set_xticks(x, labels)
    ax.set_ylim(0.90, 0.96)
    ax.set_ylabel("Strict logical fidelity at |chi|T/2pi = 3")
    ax.set_title("Single-target phase compilation gives a real but limited gain")
    ax.legend(loc="lower right")

    ax = axes[1]
    image = ax.imshow(MULTIBRANCH_FIDELITY, cmap="viridis", vmin=0.0, vmax=1.0, aspect="auto")
    ax.set_xticks(np.arange(len(MULTIBRANCH_ANGLES)), MULTIBRANCH_ANGLES)
    ax.set_yticks(np.arange(len(MULTIBRANCH_TARGETS)), MULTIBRANCH_TARGETS)
    ax.set_xlabel("Requested rotation angle")
    ax.set_ylabel("Targeted cavity branches")
    ax.set_title("Common multitone SQR fails as targets become experimentally useful")
    for i in range(MULTIBRANCH_FIDELITY.shape[0]):
        for j in range(MULTIBRANCH_FIDELITY.shape[1]):
            value = MULTIBRANCH_FIDELITY[i, j]
            text = f"{value:.3f}" if value >= 1e-3 else "<1e-3"
            ax.text(j, i, text, ha="center", va="center", color="white", fontsize=8)
    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04, label="Strict logical fidelity")

    fig.suptitle(
        "Closed-system lesson: single-branch selective control works, common multitone does not",
        fontsize=13,
    )
    for ext in ("pdf", "png"):
        fig.savefig(FIG_DIR / f"sqr_control_summary.{ext}", bbox_inches="tight")
    plt.close(fig)


def make_noise_summary() -> None:
    fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.6), constrained_layout=True)

    ax = axes[0]
    labels = list(FAMILY_MEAN_NOISY)
    values = [FAMILY_MEAN_NOISY[label] for label in labels]
    ax.bar(labels, values, color=["#4E79A7", "#F28E2B", "#59A14F"])
    ax.set_ylim(0.85, 0.91)
    ax.set_ylabel("Mean noisy target fidelity")
    ax.set_title("Noise shifts the practical winner toward shorter pulses")
    for idx, value in enumerate(values):
        ax.text(idx, value + 0.001, f"{value:.3f}", ha="center", va="bottom", fontsize=8)
    ax.tick_params(axis="x", rotation=12)

    ax = axes[1]
    ax.plot(
        GRAPE_REPLAY["chi_t"],
        GRAPE_REPLAY["grape"],
        marker="o",
        linewidth=2.0,
        color="#E15759",
        label="Archived GRAPE replay",
    )
    ax.plot(
        GRAPE_REPLAY["chi_t"],
        GRAPE_REPLAY["square"],
        marker="s",
        linewidth=2.0,
        color="#4E79A7",
        label="Best parametric square",
    )
    ax.set_xlabel(r"$|\chi|T/2\pi$")
    ax.set_ylabel("Noisy target fidelity")
    ax.set_ylim(0.84, 0.99)
    ax.set_title("GRAPE keeps a 3-6 point fidelity margin under replayed noise")
    ax.grid(alpha=0.25)
    ax.legend(loc="lower left")

    ax = axes[2]
    ax.plot(
        CONCURRENT_READOUT["amplitude_mhz"],
        CONCURRENT_READOUT["fidelity"],
        marker="o",
        linewidth=2.2,
        color="#B07AA1",
    )
    ax.set_xlabel(r"Readout amplitude $\varepsilon_r/2\pi$ (MHz)")
    ax.set_ylabel("Reduced logical fidelity")
    ax.set_ylim(0.60, 0.85)
    ax.set_title("Concurrent readout is a real control penalty even with chi_sr = 0")
    ax.grid(alpha=0.25)

    fig.suptitle(
        "Experiment-facing SQR guidance: shorter noisy primitives and temporal separation from readout",
        fontsize=13,
    )
    for ext in ("pdf", "png"):
        fig.savefig(FIG_DIR / f"sqr_noise_summary.{ext}", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    _configure_style()
    make_control_summary()
    make_noise_summary()


if __name__ == "__main__":
    main()
