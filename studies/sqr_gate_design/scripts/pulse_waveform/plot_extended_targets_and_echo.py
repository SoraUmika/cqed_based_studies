"""Generate figures for the extended generalized-target and echoed-SQR study."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
STUDY_DIR = SCRIPT_DIR.parent
DATA_DIR = STUDY_DIR / "data"
FIG_DIR = STUDY_DIR / "figures"

style_path = STUDY_DIR.parents[1] / ".github" / "skills" / "publication-figures" / "assets" / "cqed_style.mplstyle"
if style_path.exists():
    plt.style.use(str(style_path))

sys.path.insert(0, str(SCRIPT_DIR))
from common import TOL_BRIGHT


FAMILY_LABELS = {
    "single_tone_gaussian": "Single-tone Gaussian",
    "multitone_one_segment": "Multitone baseline",
    "echoed_single_tone_gaussian": "Echoed single-tone",
    "echoed_multitone_one_segment": "Echoed multitone",
}
COLOR_BY_FAMILY = {
    "single_tone_gaussian": TOL_BRIGHT[0],
    "multitone_one_segment": TOL_BRIGHT[3],
    "echoed_single_tone_gaussian": TOL_BRIGHT[1],
    "echoed_multitone_one_segment": TOL_BRIGHT[2],
}
MARKER_BY_FAMILY = {
    "single_tone_gaussian": "o",
    "multitone_one_segment": "D",
    "echoed_single_tone_gaussian": "s",
    "echoed_multitone_one_segment": "^",
}


def save_fig(fig, stem: str):
    fig.savefig(FIG_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {stem}.png/.pdf")


def load_results():
    data = np.load(DATA_DIR / "extended_targets_results.npz", allow_pickle=True)
    result = {key: data[key] for key in data.files}
    result["family_names"] = list(result["family_names"])
    result["axis_labels"] = list(result["axis_labels"])
    result["angle_labels"] = list(result["angle_labels"])
    return result


def plot_representative_family_comparison(data):
    chi_t = data["chi_t_values"]
    families = data["family_names"]
    metrics = [
        ("representative_scan_branch_true_mean", "Branch-average true-SQR fidelity"),
        ("representative_scan_branch_cphase_mean", "Branch-average conditional-phase fidelity"),
        ("representative_scan_joint_strict_fidelity", "Full logical-space fidelity"),
        ("representative_scan_joint_best_block_phase_fidelity", "Logical fidelity after cavity block-phase fit"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    axes = axes.flatten()

    for ax, (metric_name, title) in zip(axes, metrics, strict=True):
        for family_index, family_name in enumerate(families):
            ax.plot(
                chi_t,
                data[metric_name][family_index],
                label=FAMILY_LABELS[family_name],
                color=COLOR_BY_FAMILY[family_name],
                marker=MARKER_BY_FAMILY[family_name],
                linewidth=1.7,
                markersize=5,
            )
        ax.set_title(title)
        ax.set_xlabel(r"$\chi T / 2\pi$")
        ax.set_ylim(0.0, 1.02)
        ax.axhline(0.999, color="gray", linestyle=":", linewidth=0.8, alpha=0.7)

    axes[0].set_ylabel("Fidelity")
    axes[2].set_ylabel("Fidelity")
    axes[0].legend(fontsize=8, loc="lower right")
    fig.suptitle("Representative X_pi selective-rotation comparison", fontsize=12)
    fig.tight_layout()
    save_fig(fig, "extended_representative_family_comparison")


def plot_axis_angle_heatmaps(data):
    chi_t = data["chi_t_values"]
    axis_labels = data["axis_labels"]
    angle_labels = data["angle_labels"]
    chi_index = int(np.argmin(np.abs(chi_t - 3.0)))

    single_index = data["family_names"].index("single_tone_gaussian")
    echo_index = data["family_names"].index("echoed_single_tone_gaussian")

    strict_single = data["axis_scan_joint_strict_fidelity"][single_index, :, :, chi_index]
    strict_echo = data["axis_scan_joint_strict_fidelity"][echo_index, :, :, chi_index]
    gain = strict_echo - strict_single

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.2))
    panels = [
        (strict_single, "Single-tone logical fidelity"),
        (strict_echo, "Echoed single-tone logical fidelity"),
        (gain, "Echo gain in logical fidelity"),
    ]

    for ax, (image_data, title) in zip(axes, panels, strict=True):
        cmap = "RdYlGn" if np.min(image_data) >= 0.0 else "coolwarm"
        im = ax.imshow(image_data, aspect="auto", cmap=cmap)
        ax.set_xticks(range(len(axis_labels)))
        ax.set_xticklabels(axis_labels)
        ax.set_yticks(range(len(angle_labels)))
        ax.set_yticklabels([label.replace("_", " ") for label in angle_labels])
        ax.set_xlabel("Rotation axis")
        ax.set_title(title + r"\n(at $\chi T / 2\pi = 3$)")
        for row in range(image_data.shape[0]):
            for col in range(image_data.shape[1]):
                ax.text(col, row, f"{image_data[row, col]:.3f}", ha="center", va="center", fontsize=8)
        fig.colorbar(im, ax=ax, shrink=0.9)

    axes[0].set_ylabel("Rotation angle")
    fig.tight_layout()
    save_fig(fig, "extended_axis_angle_heatmaps")


def plot_branch_and_truncation_sensitivity(data):
    chi_t = data["chi_t_values"]
    branch_cases = data["target_branch_cases"]
    logical_n_cases = data["logical_n_cases"]
    chi_index = int(np.argmin(np.abs(chi_t - 3.0)))
    families = data["family_names"]

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))

    for family_index, family_name in enumerate(families):
        axes[0].plot(
            branch_cases,
            data["branch_scan_joint_strict_fidelity"][family_index, :, chi_index],
            label=FAMILY_LABELS[family_name],
            color=COLOR_BY_FAMILY[family_name],
            marker=MARKER_BY_FAMILY[family_name],
            linewidth=1.7,
        )
        axes[1].plot(
            logical_n_cases,
            data["trunc_scan_joint_strict_fidelity"][family_index, :, chi_index],
            label=FAMILY_LABELS[family_name],
            color=COLOR_BY_FAMILY[family_name],
            marker=MARKER_BY_FAMILY[family_name],
            linewidth=1.7,
        )

    axes[0].set_title(r"Logical fidelity vs target branch\n($\theta=\pi$, X-axis, $\chi T / 2\pi = 3$)")
    axes[0].set_xlabel("Target branch n0")
    axes[0].set_ylabel("Full logical-space fidelity")
    axes[0].set_xticks(branch_cases)
    axes[0].set_ylim(0.0, 1.02)

    axes[1].set_title(r"Logical fidelity vs truncated logical subspace\n($\theta=\pi$, X-axis, $n_0=1$, $\chi T / 2\pi = 3$)")
    axes[1].set_xlabel("Logical Fock levels kept")
    axes[1].set_ylabel("Full logical-space fidelity")
    axes[1].set_xticks(logical_n_cases)
    axes[1].set_ylim(0.0, 1.02)
    axes[1].legend(fontsize=8, loc="lower right")

    fig.tight_layout()
    save_fig(fig, "extended_branch_truncation_sensitivity")


def plot_cavity_phase_effects(data):
    chi_t = data["chi_t_values"]
    families = data["family_names"]

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))

    for family_index, family_name in enumerate(families):
        strict = data["representative_scan_joint_strict_fidelity"][family_index]
        fitted = data["representative_scan_joint_best_block_phase_fidelity"][family_index]
        phase_gap = fitted - strict
        axes[0].plot(
            chi_t,
            data["representative_scan_block_global_phase_spread"][family_index],
            label=FAMILY_LABELS[family_name],
            color=COLOR_BY_FAMILY[family_name],
            marker=MARKER_BY_FAMILY[family_name],
            linewidth=1.7,
        )
        axes[1].plot(
            chi_t,
            phase_gap,
            label=FAMILY_LABELS[family_name],
            color=COLOR_BY_FAMILY[family_name],
            marker=MARKER_BY_FAMILY[family_name],
            linewidth=1.7,
        )

    axes[0].set_title("Branch-global phase spread")
    axes[0].set_xlabel(r"$\chi T / 2\pi$")
    axes[0].set_ylabel("Phase spread (rad)")

    axes[1].set_title("Gain from fitted cavity block phase")
    axes[1].set_xlabel(r"$\chi T / 2\pi$")
    axes[1].set_ylabel(r"$F_{\mathrm{best-fit}} - F_{\mathrm{strict}}$")
    axes[1].legend(fontsize=8, loc="upper right")

    fig.tight_layout()
    save_fig(fig, "extended_cavity_phase_effects")


def main():
    data = load_results()
    plot_representative_family_comparison(data)
    plot_axis_angle_heatmaps(data)
    plot_branch_and_truncation_sensitivity(data)
    plot_cavity_phase_effects(data)
    print("Extended figures complete.")


if __name__ == "__main__":
    main()
