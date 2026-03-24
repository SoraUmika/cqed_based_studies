"""
Plot follow-up multitone study results.

Figures produced:
  1. followup_fidelity_comparison   — F_strict and F_block vs χT/2π for all families
  2. followup_optimized_vs_baseline — side-by-side strict vs block fidelity improvement
  3. followup_leakage_comparison    — leakage for optimized families vs baselines
  4. followup_grape_comparison      — GRAPE benchmark overlay
  5. followup_branch_scan           — baseline vs optimized across target branches
  6. followup_angle_scan            — heatmap of optimized improvement across angles
  7. followup_spectator_metrics     — spectator phase spread and transverse error
  8. followup_summary_table         — printed summary table (console + saved to data/)

Usage:
    python scripts/plot_followup_multitone.py

Input:
    data/followup_multitone_results.npz
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
STUDY_DIR = SCRIPT_DIR.parent
FIG_DIR = STUDY_DIR / "figures"
DATA_DIR = STUDY_DIR / "data"

# Load publication style
style_path = (
    STUDY_DIR.parents[1]
    / ".github" / "skills" / "publication-figures" / "assets" / "cqed_style.mplstyle"
)
if style_path.exists():
    plt.style.use(str(style_path))

sys.path.insert(0, str(SCRIPT_DIR))
from common import TOL_BRIGHT

FIG_DIR.mkdir(parents=True, exist_ok=True)

# ===================================================================
# Colors, labels, markers
# ===================================================================
BLUE, RED, GREEN, YELLOW, CYAN, PURPLE, GREY = TOL_BRIGHT

FAMILY_CFG = {
    "single_tone_gaussian": {"label": "Gaussian (baseline)", "color": BLUE, "marker": "o", "ls": "-"},
    "cosine_squared":       {"label": r"Cosine$^2$ (baseline)", "color": GREY, "marker": "s", "ls": "--"},
    "opt_indep":            {"label": "Optimized indep-tone", "color": RED, "marker": "^", "ls": "-"},
    "opt_detuned":          {"label": "Optimized + detuning", "color": GREEN, "marker": "D", "ls": "-"},
    "opt_smooth":           {"label": "Smooth-basis", "color": YELLOW, "marker": "v", "ls": "-"},
    "opt_2seg":             {"label": "2-segment", "color": PURPLE, "marker": "P", "ls": "-"},
    "grape_cphase":         {"label": "GRAPE (cphase)", "color": CYAN, "marker": "*", "ls": ":"},
    "grape_true":           {"label": "GRAPE (true)", "color": CYAN, "marker": "X", "ls": "-."},
}


def save_fig(fig, name):
    fig.savefig(FIG_DIR / f"{name}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {name}.png/.pdf")


# ===================================================================
# Load data
# ===================================================================
def load_followup():
    p = DATA_DIR / "followup_multitone_results.npz"
    if not p.exists():
        raise FileNotFoundError(f"Data file not found: {p}")
    d = np.load(p, allow_pickle=True)
    return {k: d[k] for k in d.files}


# ===================================================================
# Figure 1: Main fidelity comparison (strict + block)
# ===================================================================
def plot_fidelity_comparison(data):
    """Two-panel plot: F_strict and F_block vs χT/2π for all families."""
    chi_t = data["chi_t_values"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True)

    for metric_key, ax, title in [
        ("strict_logical_fid", ax1, r"Strict logical fidelity $\mathcal{F}_{\mathrm{strict}}$"),
        ("block_phase_relaxed_fid", ax2, r"Block-phase-relaxed fidelity $\mathcal{F}_{\mathrm{block}}$"),
    ]:
        # Baselines
        baseline_names = ["single_tone_gaussian", "cosine_squared"]
        for bi, bname in enumerate(baseline_names):
            cfg = FAMILY_CFG[bname]
            y = data[f"baseline_{metric_key}"][bi]
            ax.plot(chi_t, y, color=cfg["color"], marker=cfg["marker"],
                    ls=cfg["ls"], label=cfg["label"], markersize=5, linewidth=1.5)

        # Optimized families
        for prefix in ["opt_indep", "opt_detuned", "opt_smooth", "opt_2seg"]:
            key = f"{prefix}_{metric_key}"
            if key in data:
                cfg = FAMILY_CFG[prefix]
                ax.plot(chi_t, data[key], color=cfg["color"], marker=cfg["marker"],
                        ls=cfg["ls"], label=cfg["label"], markersize=6, linewidth=1.5)

        # GRAPE (on block panel only, using cphase fidelity as proxy for block)
        if metric_key == "block_phase_relaxed_fid" and "grape_chi_t" in data:
            cfg = FAMILY_CFG["grape_cphase"]
            ax.plot(data["grape_chi_t"], data["grape_cphase_fid"],
                    color=cfg["color"], marker=cfg["marker"], ls=cfg["ls"],
                    label=cfg["label"], markersize=8, linewidth=1.5)

        ax.set_xlabel(r"$\chi T / 2\pi$")
        ax.set_title(title)
        ax.set_ylim(-0.02, 1.05)

    ax1.set_ylabel("Fidelity")
    ax2.legend(fontsize=7, loc="lower right", framealpha=0.9)

    fig.tight_layout()
    save_fig(fig, "followup_fidelity_comparison")


# ===================================================================
# Figure 2: Improvement ratio (optimized / baseline)
# ===================================================================
def plot_improvement(data):
    """Bar chart showing fidelity improvement at each χT value."""
    chi_t = data["chi_t_values"]
    n_chi = len(chi_t)
    baseline = data["baseline_block_phase_relaxed_fid"][0]  # Gaussian baseline

    families = ["opt_indep", "opt_detuned", "opt_smooth", "opt_2seg"]
    present = [f for f in families if f"{f}_block_phase_relaxed_fid" in data]

    fig, ax = plt.subplots(figsize=(8, 4))
    x = np.arange(n_chi)
    width = 0.8 / (len(present) + 1)

    # Baseline bars
    ax.bar(x - 0.4 + width / 2, baseline, width, label="Gaussian",
           color=BLUE, alpha=0.7, edgecolor="black", linewidth=0.5)

    for i, fam in enumerate(present):
        cfg = FAMILY_CFG[fam]
        vals = data[f"{fam}_block_phase_relaxed_fid"]
        ax.bar(x - 0.4 + (i + 1.5) * width, vals, width, label=cfg["label"],
               color=cfg["color"], alpha=0.8, edgecolor="black", linewidth=0.5)

    ax.set_xticks(x)
    ax.set_xticklabels([f"{v:.1f}" for v in chi_t])
    ax.set_xlabel(r"$\chi T / 2\pi$")
    ax.set_ylabel(r"$\mathcal{F}_{\mathrm{block}}$")
    ax.set_ylim(0, 1.05)
    ax.legend(fontsize=7, ncol=2, loc="upper left")
    fig.tight_layout()
    save_fig(fig, "followup_optimized_vs_baseline")


# ===================================================================
# Figure 3: Leakage comparison
# ===================================================================
def plot_leakage(data):
    """Leakage vs χT for baselines and optimized families."""
    chi_t = data["chi_t_values"]

    fig, ax = plt.subplots(figsize=(6, 4))

    # Baseline
    cfg = FAMILY_CFG["single_tone_gaussian"]
    ax.plot(chi_t, data["baseline_leakage_max"][0], color=cfg["color"],
            marker=cfg["marker"], ls=cfg["ls"], label=cfg["label"], markersize=5)

    for prefix in ["opt_indep", "opt_detuned", "opt_smooth", "opt_2seg"]:
        key = f"{prefix}_leakage_max"
        if key in data:
            cfg = FAMILY_CFG[prefix]
            ax.plot(chi_t, data[key], color=cfg["color"], marker=cfg["marker"],
                    ls=cfg["ls"], label=cfg["label"], markersize=5)

    ax.set_xlabel(r"$\chi T / 2\pi$")
    ax.set_ylabel("Max leakage")
    ax.set_yscale("log")
    ax.legend(fontsize=7)
    fig.tight_layout()
    save_fig(fig, "followup_leakage_comparison")


# ===================================================================
# Figure 4: GRAPE comparison
# ===================================================================
def plot_grape_comparison(data):
    """Optimized multitone fidelity vs GRAPE upper bound."""
    if "grape_chi_t" not in data:
        print("  Skipping GRAPE comparison (no data)")
        return

    chi_t_all = data["chi_t_values"]
    grape_ct = data["grape_chi_t"]

    fig, ax = plt.subplots(figsize=(6, 4))

    # Baselines (Gaussian strict and block)
    ax.plot(chi_t_all, data["baseline_strict_logical_fid"][0],
            color=BLUE, marker="o", ls="--", label="Gaussian (strict)", alpha=0.6)
    ax.plot(chi_t_all, data["baseline_block_phase_relaxed_fid"][0],
            color=BLUE, marker="s", ls="-", label="Gaussian (block)", alpha=0.6)

    # Best optimized (indep)
    if "opt_indep_block_phase_relaxed_fid" in data:
        ax.plot(chi_t_all, data["opt_indep_block_phase_relaxed_fid"],
                color=RED, marker="^", ls="-", label="Opt indep (block)", linewidth=1.5)
    if "opt_indep_strict_logical_fid" in data:
        ax.plot(chi_t_all, data["opt_indep_strict_logical_fid"],
                color=RED, marker="v", ls="--", label="Opt indep (strict)", linewidth=1.5)

    # GRAPE
    ax.plot(grape_ct, data["grape_cphase_fid"], color=CYAN, marker="*",
            ls=":", label="GRAPE (cphase)", markersize=10, linewidth=2)
    ax.plot(grape_ct, data["grape_true_fid"], color=CYAN, marker="X",
            ls="-.", label="GRAPE (true)", markersize=8, linewidth=2)

    ax.set_xlabel(r"$\chi T / 2\pi$")
    ax.set_ylabel("Fidelity")
    ax.set_ylim(-0.02, 1.05)
    ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    save_fig(fig, "followup_grape_comparison")


# ===================================================================
# Figure 5: Branch scan
# ===================================================================
def plot_branch_scan(data):
    """Baseline vs optimized block fidelity across target branches."""
    if "scan_branches" not in data:
        print("  Skipping branch scan (no data)")
        return

    branches = data["scan_branches"]
    scan_ct = data["scan_chi_t"]

    fig, axes = plt.subplots(1, len(branches), figsize=(3.5 * len(branches), 4),
                             sharey=True)
    if len(branches) == 1:
        axes = [axes]

    for bi, n0 in enumerate(branches):
        ax = axes[bi]
        base_key = "branch_scan_baseline_gauss_block_phase_relaxed_fid"
        opt_key = "branch_scan_opt_indep_block_phase_relaxed_fid"
        if base_key in data:
            ax.plot(scan_ct, data[base_key][bi], color=BLUE, marker="o",
                    ls="-", label="Gaussian")
        if opt_key in data:
            ax.plot(scan_ct, data[opt_key][bi], color=RED, marker="^",
                    ls="-", label="Opt indep")
        ax.set_xlabel(r"$\chi T / 2\pi$")
        ax.set_title(f"Target branch $n_0 = {int(n0)}$")
        ax.set_ylim(-0.02, 1.05)

    axes[0].set_ylabel(r"$\mathcal{F}_{\mathrm{block}}$")
    axes[-1].legend(fontsize=7)
    fig.tight_layout()
    save_fig(fig, "followup_branch_scan")


# ===================================================================
# Figure 6: Angle scan heatmap
# ===================================================================
def plot_angle_scan(data):
    """Heatmap of improvement factor across (θ, φ) at best χT."""
    if "scan_thetas" not in data:
        print("  Skipping angle scan (no data)")
        return

    thetas = data["scan_thetas"]
    phis = data["scan_phis"]
    scan_ct = data["scan_chi_t"]

    base_key = "angle_scan_baseline_gauss_block_phase_relaxed_fid"
    opt_key = "angle_scan_opt_indep_block_phase_relaxed_fid"
    if base_key not in data or opt_key not in data:
        print("  Skipping angle scan (missing keys)")
        return

    base = data[base_key]  # (n_theta, n_phi, n_chi)
    opt = data[opt_key]    # (n_theta, n_phi, n_chi)

    fig, axes = plt.subplots(1, len(scan_ct), figsize=(3 * len(scan_ct), 3.5),
                             sharey=True)
    if len(scan_ct) == 1:
        axes = [axes]

    theta_labels = [f"{t/np.pi:.2f}π" for t in thetas]
    phi_labels = [f"{p/np.pi:.2f}π" for p in phis]

    for ci, chi_t in enumerate(scan_ct):
        ax = axes[ci]
        improvement = opt[:, :, ci] - base[:, :, ci]
        im = ax.imshow(improvement, aspect="auto", cmap="RdBu_r",
                       vmin=-0.3, vmax=0.3, origin="lower")
        ax.set_xticks(range(len(phis)))
        ax.set_xticklabels(phi_labels, fontsize=7)
        ax.set_xlabel(r"$\varphi$")
        ax.set_title(rf"$\chi T/2\pi = {chi_t:.0f}$", fontsize=9)
        if ci == 0:
            ax.set_yticks(range(len(thetas)))
            ax.set_yticklabels(theta_labels, fontsize=7)
            ax.set_ylabel(r"$\theta$")

    cbar = fig.colorbar(im, ax=axes, shrink=0.9, pad=0.02)
    cbar.set_label(r"$\Delta\mathcal{F}_{\mathrm{block}}$ (opt $-$ baseline)", fontsize=8)
    fig.tight_layout()
    save_fig(fig, "followup_angle_scan")


# ===================================================================
# Figure 7: Spectator metrics
# ===================================================================
def plot_spectator_metrics(data):
    """Spectator phase spread and transverse error for all families."""
    chi_t = data["chi_t_values"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    for metric_key, ax, ylabel in [
        ("spectator_phase_spread", ax1, "Spectator phase spread (rad)"),
        ("spectator_max_transverse", ax2, "Max spectator transverse error"),
    ]:
        # Gaussian baseline
        cfg = FAMILY_CFG["single_tone_gaussian"]
        ax.plot(chi_t, data[f"baseline_{metric_key}"][0], color=cfg["color"],
                marker=cfg["marker"], ls=cfg["ls"], label=cfg["label"], markersize=5)

        for prefix in ["opt_indep", "opt_detuned", "opt_smooth", "opt_2seg"]:
            key = f"{prefix}_{metric_key}"
            if key in data:
                cfg = FAMILY_CFG[prefix]
                ax.plot(chi_t, data[key], color=cfg["color"], marker=cfg["marker"],
                        ls=cfg["ls"], label=cfg["label"], markersize=5)

        ax.set_xlabel(r"$\chi T / 2\pi$")
        ax.set_ylabel(ylabel)

    ax2.legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    save_fig(fig, "followup_spectator_metrics")


# ===================================================================
# Summary table (printed + saved)
# ===================================================================
def print_summary(data):
    """Print a formatted summary table to console and save to file."""
    chi_t = data["chi_t_values"]

    lines = []
    lines.append("=" * 100)
    lines.append("Follow-up Multitone SQR Study: Summary Table")
    lines.append("=" * 100)
    lines.append(f"{'Family':<25s} {'χT/2π':>6s}  {'F_strict':>9s}  {'F_block':>9s}  "
                 f"{'F_cphase':>9s}  {'Leak_max':>9s}  {'Spec_φ':>8s}")
    lines.append("-" * 100)

    # Baselines
    for bi, bname in enumerate(["single_tone_gaussian", "cosine_squared"]):
        label = FAMILY_CFG[bname]["label"]
        for ci, ct in enumerate(chi_t):
            strict = data["baseline_strict_logical_fid"][bi, ci]
            block = data["baseline_block_phase_relaxed_fid"][bi, ci]
            cph = data["baseline_branch_cphase_mean"][bi, ci]
            leak = data["baseline_leakage_max"][bi, ci]
            spec = data["baseline_spectator_phase_spread"][bi, ci]
            lines.append(f"{label:<25s} {ct:6.1f}  {strict:9.4f}  {block:9.4f}  "
                         f"{cph:9.4f}  {leak:9.5f}  {spec:8.4f}")

    lines.append("-" * 100)

    # Optimized
    for prefix in ["opt_indep", "opt_detuned", "opt_smooth", "opt_2seg"]:
        key_check = f"{prefix}_strict_logical_fid"
        if key_check not in data:
            continue
        label = FAMILY_CFG[prefix]["label"]
        for ci, ct in enumerate(chi_t):
            strict = data[f"{prefix}_strict_logical_fid"][ci]
            block = data[f"{prefix}_block_phase_relaxed_fid"][ci]
            cph = data[f"{prefix}_branch_cphase_mean"][ci]
            leak = data[f"{prefix}_leakage_max"][ci]
            spec = data[f"{prefix}_spectator_phase_spread"][ci]
            lines.append(f"{label:<25s} {ct:6.1f}  {strict:9.4f}  {block:9.4f}  "
                         f"{cph:9.4f}  {leak:9.5f}  {spec:8.4f}")

    lines.append("-" * 100)

    # GRAPE
    if "grape_chi_t" in data:
        for ci, ct in enumerate(data["grape_chi_t"]):
            cp_fid = data["grape_cphase_fid"][ci]
            tr_fid = data["grape_true_fid"][ci]
            lines.append(f"{'GRAPE (cphase)':<25s} {ct:6.1f}  {'---':>9s}  "
                         f"{cp_fid:9.6f}  {'---':>9s}  {'---':>9s}  {'---':>8s}")
            lines.append(f"{'GRAPE (true)':<25s} {ct:6.1f}  {tr_fid:9.6f}  "
                         f"{'---':>9s}  {'---':>9s}  {'---':>9s}  {'---':>8s}")

    lines.append("=" * 100)

    text = "\n".join(lines)
    print(text)

    summary_path = DATA_DIR / "followup_summary.txt"
    summary_path.write_text(text, encoding="utf-8")
    print(f"\n  Summary saved to {summary_path}")


# ===================================================================
# Main
# ===================================================================
def main():
    print("Loading follow-up multitone results...")
    data = load_followup()
    print(f"  Keys: {sorted(data.keys())[:10]}... ({len(data)} total)\n")

    print("Generating figures...")
    plot_fidelity_comparison(data)
    plot_improvement(data)
    plot_leakage(data)
    plot_grape_comparison(data)
    plot_branch_scan(data)
    plot_angle_scan(data)
    plot_spectator_metrics(data)
    print("\nGenerating summary table...")
    print_summary(data)
    print("\nAll figures generated.")


if __name__ == "__main__":
    main()
