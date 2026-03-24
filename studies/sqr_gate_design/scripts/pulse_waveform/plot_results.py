"""
Generate publication-quality figures from all simulation phases.

Phase 1+2: Baseline scan (ideal model, closed-system)
Phase 4:   Higher-order corrections (χ', K, closed-system)
Phase 5:   Open-system (χ', K, T1, T2)
GRAPE:     Optimal control upper bound (if available)

Usage:
    python scripts/plot_results.py

Input:
    data/phase1_phase2_results.npz
    data/phase4_results.npz
    data/phase5_results.npz
    data/grape_benchmark_results.npz  (optional)

Output:
    figures/fidelity_vs_chi_t.{png,pdf}
    figures/spectator_phase_vs_chi_t.{png,pdf}
    figures/branch_fidelity_heatmap.{png,pdf}
    figures/true_vs_cphase_comparison.{png,pdf}
    figures/higher_order_comparison.{png,pdf}         (Phase 4)
    figures/decoherence_fidelity.{png,pdf}            (Phase 5)
    figures/net_fidelity_vs_chi_t.{png,pdf}           (Phase 5)
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
STUDY_DIR = SCRIPT_DIR.parent
FIG_DIR = STUDY_DIR / "figures"
DATA_DIR = STUDY_DIR / "data"

# Load publication style
style_path = STUDY_DIR.parents[1] / ".github" / "skills" / "publication-figures" / "assets" / "cqed_style.mplstyle"
if style_path.exists():
    plt.style.use(str(style_path))

sys.path.insert(0, str(SCRIPT_DIR))
from common import TOL_BRIGHT

FAMILY_LABELS = {
    "single_tone_gaussian": "Gaussian",
    "square": "Square",
    "cosine_squared": r"Cosine$^2$ (Hann)",
    "multitone_one_segment": "Multitone (Gaussian)",
}
FAMILY_COLORS = {
    "single_tone_gaussian": TOL_BRIGHT[0],
    "square": TOL_BRIGHT[1],
    "cosine_squared": TOL_BRIGHT[2],
    "multitone_one_segment": TOL_BRIGHT[3],
}
FAMILY_MARKERS = {
    "single_tone_gaussian": "o",
    "square": "s",
    "cosine_squared": "^",
    "multitone_one_segment": "D",
}


def save_fig(fig, name):
    fig.savefig(FIG_DIR / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {name}.png/.pdf")


def load_data():
    d = np.load(DATA_DIR / "phase1_phase2_results.npz", allow_pickle=True)
    data = {k: d[k] for k in d.files}
    data["family_names"] = list(data["family_names"])
    return data


def load_grape():
    p = DATA_DIR / "grape_benchmark_results.npz"
    if p.exists():
        d = np.load(p, allow_pickle=True)
        return {k: d[k] for k in d.files}
    return None


def load_phase4():
    p = DATA_DIR / "phase4_results.npz"
    if p.exists():
        d = np.load(p, allow_pickle=True)
        data = {k: d[k] for k in d.files}
        data["family_names"] = list(data["family_names"])
        return data
    return None


def load_phase5():
    p = DATA_DIR / "phase5_results.npz"
    if p.exists():
        d = np.load(p, allow_pickle=True)
        data = {k: d[k] for k in d.files}
        data["family_names"] = list(data["family_names"])
        return data
    return None


# -----------------------------------------------------------------------
# Figure 1: Fidelity vs χT (both targets, all families)
# -----------------------------------------------------------------------
def plot_fidelity_vs_chi_t(data, grape=None):
    chi_t = data["chi_t_values"]
    families = data["family_names"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

    for ax, fid_key, title in zip(
        axes,
        ["true_sqr_fidelity", "cphase_sqr_fidelity"],
        ["True SQR target", "Conditional-phase SQR target"],
    ):
        for i, fam in enumerate(families):
            ax.plot(
                chi_t, data[fid_key][i],
                marker=FAMILY_MARKERS[fam], color=FAMILY_COLORS[fam],
                label=FAMILY_LABELS[fam], linewidth=1.5, markersize=5,
            )

        if grape is not None:
            if fid_key == "true_sqr_fidelity" and "fidelity_true" in grape:
                ax.plot(
                    grape["chi_t_values"], grape["fidelity_true"],
                    marker="*", color=TOL_BRIGHT[5], linewidth=1.5,
                    markersize=8, label="GRAPE (true SQR)", linestyle="--",
                )
            elif fid_key == "cphase_sqr_fidelity" and "fidelity_cphase" in grape:
                ax.plot(
                    grape["chi_t_values"], grape["fidelity_cphase"],
                    marker="*", color=TOL_BRIGHT[5], linewidth=1.5,
                    markersize=8, label="GRAPE (cphase SQR)", linestyle="--",
                )

        ax.axhline(0.999, color="gray", linestyle=":", linewidth=0.8, alpha=0.7)
        ax.set_xlabel(r"$\chi T / 2\pi$")
        ax.set_title(title)
        ax.set_ylim(-0.05, 1.05)

    axes[0].set_ylabel("Subspace process fidelity")
    axes[0].legend(fontsize=8, loc="lower right")

    fig.suptitle(
        r"SQR fidelity vs $\chi T / 2\pi$"
        f"  (target branch $n_0={int(data['target_n0'])}$, "
        r"$\theta = \pi$)",
        fontsize=12,
    )
    fig.tight_layout()
    save_fig(fig, "fidelity_vs_chi_t")


# -----------------------------------------------------------------------
# Figure 2: Spectator phase spread vs χT
# -----------------------------------------------------------------------
def plot_spectator_phase(data):
    chi_t = data["chi_t_values"]
    families = data["family_names"]

    fig, ax = plt.subplots(figsize=(7, 5))
    for i, fam in enumerate(families):
        ax.plot(
            chi_t, data["spectator_phase_spread"][i],
            marker=FAMILY_MARKERS[fam], color=FAMILY_COLORS[fam],
            label=FAMILY_LABELS[fam], linewidth=1.5, markersize=5,
        )
    ax.set_xlabel(r"$\chi T / 2\pi$")
    ax.set_ylabel(r"Spectator phase spread $\Delta\varphi_{\mathrm{spec}}$ (rad)")
    ax.legend(fontsize=8)
    ax.set_title("Spectator phase spread vs gate duration")
    fig.tight_layout()
    save_fig(fig, "spectator_phase_vs_chi_t")


# -----------------------------------------------------------------------
# Figure 3: Branch-resolved fidelity heatmap
# -----------------------------------------------------------------------
def plot_branch_fidelity_heatmap(data):
    chi_t = data["chi_t_values"]
    families = data["family_names"]
    n_fock = int(data["n_fock"])

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    for i, fam in enumerate(families):
        ax = axes[i]
        fid_data = data["branch_fidelities_true"][i]  # shape (n_chi_t, n_fock)

        im = ax.imshow(
            fid_data.T, aspect="auto", origin="lower",
            extent=[chi_t[0], chi_t[-1], -0.5, n_fock - 0.5],
            cmap="RdYlGn", vmin=0.0, vmax=1.0,
        )
        ax.set_xlabel(r"$\chi T / 2\pi$")
        ax.set_ylabel("Fock level $n$")
        ax.set_yticks(range(n_fock))
        ax.set_title(FAMILY_LABELS[fam], fontsize=10)
        fig.colorbar(im, ax=ax, label="Branch fidelity (true SQR)")

    fig.suptitle("Branch-resolved fidelity (true SQR target)", fontsize=12)
    fig.tight_layout()
    save_fig(fig, "branch_fidelity_heatmap")


# -----------------------------------------------------------------------
# Figure 4: True SQR vs conditional-phase SQR comparison
# -----------------------------------------------------------------------
def plot_true_vs_cphase(data):
    chi_t = data["chi_t_values"]
    families = data["family_names"]

    fig, ax = plt.subplots(figsize=(7, 5))

    for i, fam in enumerate(families):
        fid_gap = data["cphase_sqr_fidelity"][i] - data["true_sqr_fidelity"][i]
        ax.plot(
            chi_t, fid_gap,
            marker=FAMILY_MARKERS[fam], color=FAMILY_COLORS[fam],
            label=FAMILY_LABELS[fam], linewidth=1.5, markersize=5,
        )

    ax.axhline(0, color="gray", linestyle=":", linewidth=0.8)
    ax.set_xlabel(r"$\chi T / 2\pi$")
    ax.set_ylabel(r"$F_{\mathrm{cphase}} - F_{\mathrm{true}}$")
    ax.set_title("Fidelity advantage of conditional-phase SQR over true SQR")
    ax.legend(fontsize=8)
    fig.tight_layout()
    save_fig(fig, "true_vs_cphase_comparison")


# -----------------------------------------------------------------------
# Figure 5: Higher-order corrections comparison (Phase 4 vs Phase 1-2)
# -----------------------------------------------------------------------
def plot_higher_order_comparison(data_p12, data_p4):
    chi_t = data_p12["chi_t_values"]
    families = data_p12["family_names"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

    for ax, fid_key, title in zip(
        axes,
        ["cphase_sqr_fidelity", "cphase_sqr_fidelity"],
        [r"Ideal model ($\chi$ only)", r"With $\chi^\prime$, $K$"],
    ):
        src = data_p12 if "only" in title else data_p4
        for i, fam in enumerate(families):
            ax.plot(
                chi_t, src[fid_key][i],
                marker=FAMILY_MARKERS[fam], color=FAMILY_COLORS[fam],
                label=FAMILY_LABELS[fam], linewidth=1.5, markersize=5,
            )
        ax.axhline(0.999, color="gray", linestyle=":", linewidth=0.8, alpha=0.7)
        ax.set_xlabel(r"$\chi T / 2\pi$")
        ax.set_title(title)
        ax.set_ylim(0.85, 1.005)

    axes[0].set_ylabel("Conditional-phase SQR fidelity")
    axes[0].legend(fontsize=8, loc="lower right")
    fig.suptitle(
        r"Effect of higher-order terms on cphase SQR fidelity"
        "\n"
        r"($\chi^\prime = 2\pi \times (-21\,\mathrm{kHz})$, "
        r"$K = 2\pi \times (-28\,\mathrm{kHz})$)",
        fontsize=11,
    )
    fig.tight_layout()
    save_fig(fig, "higher_order_comparison")


# -----------------------------------------------------------------------
# Figure 6: Decoherence fidelity (Phase 5)
# -----------------------------------------------------------------------
def plot_decoherence_fidelity(data_p5):
    chi_t = data_p5["chi_t_values"]
    families = data_p5["family_names"]

    # Duration from χT/(2π): T = chi_t_2pi / (|χ|/(2π))
    chi_rad = float(data_p5["chi_rad_s"])
    f_chi = abs(chi_rad) / (2 * np.pi)  # Hz

    fig, ax = plt.subplots(figsize=(7, 5))
    for i, fam in enumerate(families):
        ax.plot(
            chi_t, data_p5["deco_fid"][i],
            marker=FAMILY_MARKERS[fam], color=FAMILY_COLORS[fam],
            label=FAMILY_LABELS[fam], linewidth=1.5, markersize=5,
        )

    # Analytic decoherence estimate: F ≈ 1 - T/(2T1)
    t1_s = float(data_p5["t1_s"])
    T_s = chi_t / f_chi
    F_analytic = 1 - T_s / (2 * t1_s)
    ax.plot(chi_t, F_analytic, "k--", linewidth=1, alpha=0.5,
            label=r"$1 - T/(2T_1)$")

    ax.set_xlabel(r"$\chi T / 2\pi$")
    ax.set_ylabel("Decoherence fidelity")
    ax.set_title(
        r"State fidelity vs closed-system ideal"
        f"\n$T_1 = T_2 = {t1_s*1e6:.0f}$ μs"
    )
    ax.legend(fontsize=8)
    ax.set_ylim(0.9, 1.005)

    # Add top x-axis showing pulse duration
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    tick_positions = ax.get_xticks()
    tick_positions = [t for t in tick_positions if ax.get_xlim()[0] <= t <= ax.get_xlim()[1]]
    ax2.set_xticks(tick_positions)
    ax2.set_xticklabels([f"{t / f_chi * 1e6:.1f}" for t in tick_positions])
    ax2.set_xlabel(r"Pulse duration ($\mu$s)")

    fig.tight_layout()
    save_fig(fig, "decoherence_fidelity")


# -----------------------------------------------------------------------
# Figure 7: Net cphase fidelity (Phase 5) — selectivity + decoherence
# -----------------------------------------------------------------------
def plot_net_fidelity(data_p12, data_p5, grape=None):
    chi_t = data_p5["chi_t_values"]
    families = data_p5["family_names"]

    fig, ax = plt.subplots(figsize=(8, 5.5))

    for i, fam in enumerate(families):
        # Phase 1-2 baseline (ideal, closed-system)
        ax.plot(
            chi_t, data_p12["cphase_sqr_fidelity"][i],
            marker=FAMILY_MARKERS[fam], color=FAMILY_COLORS[fam],
            linewidth=1.5, markersize=5, linestyle="--", alpha=0.4,
        )
        # Phase 5 net fidelity (with χ', K, T1, T2)
        ax.plot(
            chi_t, data_p5["cphase_fid_net"][i],
            marker=FAMILY_MARKERS[fam], color=FAMILY_COLORS[fam],
            label=FAMILY_LABELS[fam], linewidth=2, markersize=6,
        )

    # GRAPE upper bound (cphase, if available)
    if grape is not None:
        grape_key = "fidelity_cphase" if "fidelity_cphase" in grape else "fidelity"
        ax.plot(
            grape["chi_t_values"], grape[grape_key],
            marker="*", color=TOL_BRIGHT[5], linewidth=1.5,
            markersize=8, label="GRAPE (cphase SQR)", linestyle="--",
        )

    ax.axhline(0.999, color="gray", linestyle=":", linewidth=0.8, alpha=0.7,
               label=r"$F = 0.999$")
    ax.set_xlabel(r"$\chi T / 2\pi$")
    ax.set_ylabel("Net conditional-phase SQR fidelity")
    ax.set_title(
        r"Net fidelity including $\chi^\prime$, $K$, $T_1$, $T_2$"
        "\n(dashed: ideal closed-system; solid: full model)"
    )
    ax.legend(fontsize=8, loc="lower right")
    ax.set_ylim(0.2, 1.02)
    fig.tight_layout()
    save_fig(fig, "net_fidelity_vs_chi_t")


# -----------------------------------------------------------------------
# Figure 8: Leakage to |f⟩ (if data available)
# -----------------------------------------------------------------------
def plot_leakage(data):
    """Plot max leakage (to |f⟩) vs χT/(2π) for all families."""
    if "leakage" not in data:
        return
    chi_t = data["chi_t_values"]
    families = data["family_names"]
    leakage = data["leakage"]  # (n_families, n_chi_t, n_fock)

    fig, ax = plt.subplots(figsize=(7, 5))
    for i, fam in enumerate(families):
        max_leak = np.max(leakage[i], axis=-1)  # max over Fock levels
        ax.semilogy(
            chi_t, max_leak,
            marker=FAMILY_MARKERS[fam], color=FAMILY_COLORS[fam],
            label=FAMILY_LABELS[fam], linewidth=1.5, markersize=5,
        )
    ax.set_xlabel(r"$\chi T / 2\pi$")
    ax.set_ylabel(r"Max leakage to $|f\rangle$")
    ax.set_title(r"Leakage to third transmon level ($n_{\mathrm{tr}}=3$)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    save_fig(fig, "leakage_vs_chi_t")


# -----------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------
def main():
    print("Loading data...")
    data = load_data()
    grape = load_grape()
    data_p4 = load_phase4()
    data_p5 = load_phase5()

    print("Generating figures...")
    plot_fidelity_vs_chi_t(data, grape)
    plot_spectator_phase(data)
    plot_branch_fidelity_heatmap(data)
    plot_true_vs_cphase(data)
    plot_leakage(data)

    if data_p4 is not None:
        plot_higher_order_comparison(data, data_p4)
    if data_p5 is not None:
        plot_decoherence_fidelity(data_p5)
        plot_net_fidelity(data, data_p5, grape)

    print("Done.")


if __name__ == "__main__":
    main()
