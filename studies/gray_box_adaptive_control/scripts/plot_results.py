"""
plot_results.py — Figure generation for the gray-box adaptive control study.

Generates all publication-quality figures from the data files produced by
study_phase4.py and study_phase5.py.

Required data files (in data/):
  phase4_results.npz
  phase5_1_chi_mismatch.npz
  phase5_2_noise_sweep.npz
  phase5_3_readout_sweep.npz
  phase5_4_probe_budget.npz
  phase5_5_drift.npz
  phase5_6_omission.npz

Outputs (in figures/):
  phase4_main_comparison.pdf
  phase4_per_fock.pdf
  phase5_1_chi_mismatch.pdf
  phase5_2_noise.pdf
  phase5_3_readout.pdf
  phase5_4_probe_budget.pdf
  phase5_5_drift.pdf
  phase5_6_omission.pdf
  viability_summary.pdf

Usage (from study directory):
    python scripts/plot_results.py
    python scripts/plot_results.py --figure phase4  # only phase4 figures
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")   # non-interactive backend for script use
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    from matplotlib.gridspec import GridSpec
    HAVE_MPL = True
except ImportError:
    HAVE_MPL = False
    warnings.warn("matplotlib not available. Install with: pip install matplotlib --user")

STUDY_DIR = Path(__file__).parent.parent
DATA_DIR = STUDY_DIR / "data"
FIG_DIR = STUDY_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Colorblind-friendly palette (Paul Tol's muted set)
# ---------------------------------------------------------------------------

COLORS = {
    "nominal": "#CC6677",    # rose / reddish
    "gray_box": "#44AA99",   # teal / green
    "perfect": "#332288",    # deep indigo
    "black_box": "#DDCC77",  # sand / gold
    "noisy": "#AA4499",      # purple (for noisy variants)
}

LINESTYLES = {
    "nominal": "-",
    "gray_box": "-",
    "perfect": "--",
    "black_box": ":",
    "noisy": "-.",
}

MARKERS = {
    "nominal": "o",
    "gray_box": "s",
    "perfect": "^",
    "black_box": "D",
}

LABELS = {
    "nominal": "Nominal (prior model)",
    "gray_box": "Gray-box (chi corrected)",
    "perfect": "Perfect knowledge",
    "black_box": "Black-box (model-free)",
}


# ---------------------------------------------------------------------------
# Figure style
# ---------------------------------------------------------------------------


def setup_style():
    if not HAVE_MPL:
        return
    plt.rcParams.update({
        "figure.dpi": 150,
        "font.size": 11,
        "axes.labelsize": 12,
        "axes.titlesize": 12,
        "legend.fontsize": 10,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "lines.linewidth": 2.0,
        "lines.markersize": 6,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


def save_fig(fig, name: str):
    path = FIG_DIR / name
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")


# ---------------------------------------------------------------------------
# Figure 1: Phase 4 main comparison
# ---------------------------------------------------------------------------


def plot_phase4_main():
    if not HAVE_MPL:
        return
    data_path = DATA_DIR / "phase4_results.npz"
    if not data_path.exists():
        print(f"  [SKIP] {data_path} not found")
        return

    d = np.load(data_path, allow_pickle=True)
    mismatches = d["mismatch_fractions"] * 100   # percent
    nominal = d["nominal_fidelities"]
    gray_box = d["gray_box_fidelities"]
    perfect = d["perfect_fidelities"]
    nominal_noisy = d["nominal_fidelities_noisy"]
    gray_box_noisy = d["gray_box_fidelities_noisy"]
    perfect_noisy = d["perfect_fidelities_noisy"]
    bb_fidelity = float(d["bb_fidelity"])
    bb_mismatch = float(d["bb_mismatch"]) * 100

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Left: infidelity vs mismatch (noiseless)
    ax = axes[0]
    ax.semilogy(mismatches, 1 - nominal, color=COLORS["nominal"],
                marker=MARKERS["nominal"], ls=LINESTYLES["nominal"], label=LABELS["nominal"])
    ax.semilogy(mismatches, 1 - gray_box, color=COLORS["gray_box"],
                marker=MARKERS["gray_box"], ls=LINESTYLES["gray_box"], label=LABELS["gray_box"])
    ax.semilogy(mismatches, 1 - perfect, color=COLORS["perfect"],
                marker=MARKERS["perfect"], ls=LINESTYLES["perfect"], label=LABELS["perfect"])
    if not np.isnan(bb_fidelity):
        ax.semilogy([bb_mismatch], [1 - bb_fidelity], color=COLORS["black_box"],
                    marker=MARKERS["black_box"], ms=10, ls="none", label=LABELS["black_box"])
    ax.set_xlabel("Chi mismatch (%)")
    ax.set_ylabel("Gate infidelity (1 - F)")
    ax.set_title("Noiseless evaluation on truth model")
    ax.legend(loc="upper left")
    ax.set_xlim(-2, max(mismatches) + 2)

    # Right: infidelity vs mismatch (noisy)
    ax = axes[1]
    ax.semilogy(mismatches, 1 - nominal_noisy, color=COLORS["nominal"],
                marker=MARKERS["nominal"], ls=LINESTYLES["noisy"],
                label=f"{LABELS['nominal']} + noise")
    ax.semilogy(mismatches, 1 - gray_box_noisy, color=COLORS["gray_box"],
                marker=MARKERS["gray_box"], ls=LINESTYLES["noisy"],
                label=f"{LABELS['gray_box']} + noise")
    ax.semilogy(mismatches, 1 - perfect_noisy, color=COLORS["perfect"],
                marker=MARKERS["perfect"], ls=LINESTYLES["noisy"],
                label=f"{LABELS['perfect']} + noise")
    ax.set_xlabel("Chi mismatch (%)")
    ax.set_ylabel("Gate infidelity (1 - F)")
    ax.set_title(f"Noisy evaluation (T1={50} μs, Tφ={40} μs)")
    ax.legend(loc="upper left")
    ax.set_xlim(-2, max(mismatches) + 2)

    fig.suptitle("Phase 4: Gray-box adaptive control comparison", fontsize=13)
    fig.tight_layout()
    save_fig(fig, "phase4_main_comparison.pdf")


# ---------------------------------------------------------------------------
# Figure 2: Phase 4 per-Fock fidelity (30% mismatch case)
# ---------------------------------------------------------------------------


def plot_phase4_per_fock():
    if not HAVE_MPL:
        return
    data_path = DATA_DIR / "phase4_results.npz"
    if not data_path.exists():
        print(f"  [SKIP] {data_path} not found")
        return

    d = np.load(data_path, allow_pickle=True)
    nom_per = d["per_state_nominal_30"]
    gb_per = d["per_state_gray_box_30"]
    perf_per = d["per_state_perfect_30"]

    if np.all(np.isnan(nom_per)):
        print("  [SKIP] Per-state fidelities not available (30% case may not have run)")
        return

    n_states = len(nom_per)
    x = np.arange(n_states)

    fig, ax = plt.subplots(figsize=(8, 5))
    width = 0.25
    ax.bar(x - width, 1 - nom_per, width, color=COLORS["nominal"], alpha=0.8, label=LABELS["nominal"])
    ax.bar(x, 1 - gb_per, width, color=COLORS["gray_box"], alpha=0.8, label=LABELS["gray_box"])
    ax.bar(x + width, 1 - perf_per, width, color=COLORS["perfect"], alpha=0.8, label=LABELS["perfect"])
    ax.set_xticks(x)
    ax.set_xticklabels([f"probe {i}" for i in range(n_states)], rotation=45, ha="right")
    ax.set_ylabel("Per-state infidelity (1 - F)")
    ax.set_title("Per-probe-state fidelity at 30% chi mismatch")
    ax.legend()
    fig.tight_layout()
    save_fig(fig, "phase4_per_fock.pdf")


# ---------------------------------------------------------------------------
# Figure 3: Phase 5.1 chi mismatch sweep
# ---------------------------------------------------------------------------


def plot_phase5_1():
    if not HAVE_MPL:
        return
    data_path = DATA_DIR / "phase5_1_chi_mismatch.npz"
    if not data_path.exists():
        print(f"  [SKIP] {data_path} not found")
        return

    d = np.load(data_path)
    ratios = d["ratios"]
    mismatches_pct = (1.0 - ratios) * 100

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.semilogy(mismatches_pct, 1 - d["nominal_fidelities"], color=COLORS["nominal"],
                marker=MARKERS["nominal"], label=LABELS["nominal"])
    ax.semilogy(mismatches_pct, 1 - d["gray_box_fidelities"], color=COLORS["gray_box"],
                marker=MARKERS["gray_box"], label=LABELS["gray_box"])
    ax.semilogy(mismatches_pct, 1 - d["perfect_fidelities"], color=COLORS["perfect"],
                marker=MARKERS["perfect"], ls="--", label=LABELS["perfect"])
    ax.set_xlabel("Chi mismatch (1 - chi_prior/chi_true) [%]")
    ax.set_ylabel("Gate infidelity (1 - F)")
    ax.set_title("Phase 5.1: Chi mismatch sweep (noiseless)")
    ax.legend()
    fig.tight_layout()
    save_fig(fig, "phase5_1_chi_mismatch.pdf")


# ---------------------------------------------------------------------------
# Figure 4: Phase 5.2 noise sweep
# ---------------------------------------------------------------------------


def plot_phase5_2():
    if not HAVE_MPL:
        return
    data_path = DATA_DIR / "phase5_2_noise_sweep.npz"
    if not data_path.exists():
        print(f"  [SKIP] {data_path} not found")
        return

    d = np.load(data_path, allow_pickle=True)
    t1_us = np.array([
        float("inf") if not np.isfinite(t) else t * 1e6
        for t in d["t1_values"]
    ])
    x = np.arange(len(t1_us))
    labels_x = [f"∞" if not np.isfinite(t) else f"{t:.0f}" for t in t1_us]

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.semilogy(x, 1 - d["nominal_fidelities"], color=COLORS["nominal"],
                marker=MARKERS["nominal"], label=LABELS["nominal"])
    ax.semilogy(x, 1 - d["gray_box_fidelities"], color=COLORS["gray_box"],
                marker=MARKERS["gray_box"], label=LABELS["gray_box"])
    ax.semilogy(x, 1 - d["perfect_fidelities"], color=COLORS["perfect"],
                marker=MARKERS["perfect"], ls="--", label=LABELS["perfect"])
    ax.set_xticks(x)
    ax.set_xticklabels([f"T1={v} μs" for v in labels_x])
    ax.set_ylabel("Gate infidelity (1 - F)")
    ax.set_title("Phase 5.2: Noise strength sweep (30% chi mismatch)")
    ax.legend()
    fig.tight_layout()
    save_fig(fig, "phase5_2_noise.pdf")


# ---------------------------------------------------------------------------
# Figure 5: Phase 5.3 readout imperfection
# ---------------------------------------------------------------------------


def plot_phase5_3():
    if not HAVE_MPL:
        return
    data_path = DATA_DIR / "phase5_3_readout_sweep.npz"
    if not data_path.exists():
        print(f"  [SKIP] {data_path} not found")
        return

    d = np.load(data_path)
    c01 = d["confusion_01_values"] * 100   # percent

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Left: chi inference error vs confusion rate
    ax = axes[0]
    ax.plot(c01, d["chi_errors"] * 100, color=COLORS["gray_box"], marker="s")
    ax.set_xlabel("Readout error rate (%)")
    ax.set_ylabel("Chi inference error |Δchi/chi_true| (%)")
    ax.set_title("(a) Chi inference error vs readout error")

    # Right: gate fidelity vs confusion rate
    ax = axes[1]
    ax.plot(c01, 1 - d["nominal_fidelities"], color=COLORS["nominal"],
            marker=MARKERS["nominal"], label=LABELS["nominal"])
    ax.plot(c01, 1 - d["gray_box_fidelities"], color=COLORS["gray_box"],
            marker=MARKERS["gray_box"], label=LABELS["gray_box"])
    ax.plot(c01, 1 - d["perfect_fidelities"], color=COLORS["perfect"],
            marker=MARKERS["perfect"], ls="--", label=LABELS["perfect"])
    ax.set_xlabel("Readout error rate (%)")
    ax.set_ylabel("Gate infidelity (1 - F)")
    ax.set_title("(b) Gate fidelity vs readout error")
    ax.legend()

    fig.suptitle("Phase 5.3: Readout imperfection sweep (30% chi mismatch, n_shots=1000)")
    fig.tight_layout()
    save_fig(fig, "phase5_3_readout.pdf")


# ---------------------------------------------------------------------------
# Figure 6: Phase 5.4 probe budget
# ---------------------------------------------------------------------------


def plot_phase5_4():
    if not HAVE_MPL:
        return
    data_path = DATA_DIR / "phase5_4_probe_budget.npz"
    if not data_path.exists():
        print(f"  [SKIP] {data_path} not found")
        return

    d = np.load(data_path)
    n_shots = d["n_shots_values"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Left: chi error vs n_shots
    ax = axes[0]
    ax.loglog(n_shots, d["chi_errors"] * 100, color=COLORS["gray_box"], marker="s", label="Chi error")
    # Theoretical shot-noise scaling: 1/sqrt(n_shots)
    n_ref = np.logspace(np.log10(n_shots[0]), np.log10(n_shots[-1]), 100)
    scale = d["chi_errors"][0] * 100 * np.sqrt(n_shots[0])
    ax.loglog(n_ref, scale / np.sqrt(n_ref), "k--", alpha=0.5, label="∝ 1/√N (shot noise)")
    ax.set_xlabel("Probe shots per time point (N)")
    ax.set_ylabel("Chi inference error (%)")
    ax.set_title("(a) Chi inference error vs probe budget")
    ax.legend()

    # Right: gate fidelity vs n_shots
    ax = axes[1]
    nom_f = float(d["nominal_fidelity"])
    perf_f = float(d["perfect_fidelity"])
    ax.semilogx(n_shots, 1 - d["gray_box_fidelities"], color=COLORS["gray_box"],
                marker="s", label=LABELS["gray_box"])
    ax.axhline(1 - nom_f, color=COLORS["nominal"], ls="--", label=f"{LABELS['nominal']} (F={nom_f:.4f})")
    ax.axhline(1 - perf_f, color=COLORS["perfect"], ls=":", label=f"{LABELS['perfect']} (F={perf_f:.4f})")
    ax.set_xlabel("Probe shots per time point (N)")
    ax.set_ylabel("Gate infidelity (1 - F)")
    ax.set_title("(b) Gate fidelity vs probe budget")
    ax.legend()

    fig.suptitle("Phase 5.4: Probe budget sweep (30% chi mismatch)")
    fig.tight_layout()
    save_fig(fig, "phase5_4_probe_budget.pdf")


# ---------------------------------------------------------------------------
# Figure 7: Phase 5.5 drift study
# ---------------------------------------------------------------------------


def plot_phase5_5():
    if not HAVE_MPL:
        return
    data_path = DATA_DIR / "phase5_5_drift.npz"
    if not data_path.exists():
        print(f"  [SKIP] {data_path} not found")
        return

    d = np.load(data_path, allow_pickle=True)
    traces = d["fidelity_traces"]  # (n_drift, n_recal, n_cycles)
    drift_labels = [str(x) for x in d["drift_labels"]]
    recal_labels = [str(x) for x in d["recal_labels"]]
    n_cycles = int(d["n_cycles"])
    cycles = np.arange(n_cycles)

    n_drift = len(drift_labels)
    n_recal = len(recal_labels)

    recal_colors = ["#4477AA", "#EE6677", "#228833"]  # blue, red, green
    recal_styles = ["-", "--", "-."]

    fig, axes = plt.subplots(1, n_drift, figsize=(4 * n_drift, 5), sharey=True)
    if n_drift == 1:
        axes = [axes]

    for di, d_label in enumerate(drift_labels):
        ax = axes[di]
        for ri, (r_label, color, ls) in enumerate(zip(recal_labels, recal_colors, recal_styles)):
            ax.plot(cycles, traces[di, ri, :], color=color, ls=ls, label=r_label, alpha=0.9)
        ax.set_xlabel("Experiment cycle")
        ax.set_title(f"Drift: {d_label}/cycle")
        if di == 0:
            ax.set_ylabel("Gate fidelity F")
        ax.legend(title="Recalibration", fontsize=9)
        ax.set_ylim(0, 1.05)

    fig.suptitle("Phase 5.5: Chi drift study – recalibration strategies", fontsize=13)
    fig.tight_layout()
    save_fig(fig, "phase5_5_drift.pdf")


# ---------------------------------------------------------------------------
# Figure 8: Phase 5.6 Hamiltonian omission
# ---------------------------------------------------------------------------


def plot_phase5_6():
    if not HAVE_MPL:
        return
    data_path = DATA_DIR / "phase5_6_omission.npz"
    if not data_path.exists():
        print(f"  [SKIP] {data_path} not found")
        return

    d = np.load(data_path)
    mult = d["chi_higher_multipliers"]
    chi_h_khz = d["chi_higher_values"] / (2 * np.pi) / 1e3

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.semilogy(mult, 1 - d["fidelities_chi_correct"], color=COLORS["gray_box"],
                marker="s", label="chi correct, chi_higher varied")
    ax.semilogy(mult, 1 - d["fidelities_chi_wrong"], color=COLORS["nominal"],
                marker="o", label="chi wrong (30%), chi_higher varied")
    ax.axvline(1.0, color="gray", ls=":", alpha=0.7, label="True chi_higher value")
    ax.axvline(0.0, color="red", ls=":", alpha=0.5, label="Learner omits chi_higher")
    ax.set_xlabel("chi_higher multiplier (0=omitted, 1=true, 2=2x true)")
    ax.set_ylabel("Gate infidelity (1 - F)")
    ax.set_title("Phase 5.6: Hamiltonian omission – chi_higher effect")
    ax.legend()

    # Secondary x-axis showing chi_higher in kHz
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(mult)
    ax2.set_xticklabels([f"{x:.1f}" for x in chi_h_khz])
    ax2.set_xlabel("chi_higher / (2π) [kHz]")

    fig.tight_layout()
    save_fig(fig, "phase5_6_omission.pdf")


# ---------------------------------------------------------------------------
# Figure 9: Viability summary heatmap
# ---------------------------------------------------------------------------


def plot_viability_summary():
    """
    Generate a summary heatmap: regime × method → viability judgment.

    Viability is based on fidelity thresholds:
      Green  (viable):     F > 0.99
      Yellow (marginal):   0.95 < F <= 0.99
      Red    (not viable): F <= 0.95
    """
    if not HAVE_MPL:
        return

    # Attempt to load data; generate representative placeholders if missing
    def load_or_nan(path, key, fallback=np.nan):
        if Path(path).exists():
            d = np.load(path, allow_pickle=True)
            return float(d[key]) if key in d else fallback
        return fallback

    def load_array(path, key):
        if Path(path).exists():
            d = np.load(path, allow_pickle=True)
            return d[key] if key in d else None
        return None

    # Methods: nominal, gray_box, perfect
    # Regimes: 0% mismatch, 30% mismatch (noiseless), 30% mismatch (noisy),
    #          30% + 5% readout err, 30% + 50 shots, 70% mismatch

    regime_labels = [
        "0% mismatch\n(noiseless)",
        "30% mismatch\n(noiseless)",
        "30% mismatch\n(noisy)",
        "30% + 5%\nreadout err",
        "30%\n50 shots",
        "70% mismatch\n(noiseless)",
    ]
    method_labels = ["Nominal", "Gray-box", "Perfect"]

    # Fill from data where available, else nan
    fidelity_matrix = np.full((len(method_labels), len(regime_labels)), np.nan)

    p4 = load_array(DATA_DIR / "phase4_results.npz", "nominal_fidelities")
    if p4 is not None:
        mismatches = load_array(DATA_DIR / "phase4_results.npz", "mismatch_fractions")
        if mismatches is not None:
            for col, target_mismatch in zip([0, 1], [0.0, 0.30]):
                idx = np.argmin(np.abs(np.array(mismatches) - target_mismatch))
                fidelity_matrix[0, col] = float(
                    load_array(DATA_DIR / "phase4_results.npz", "nominal_fidelities")[idx]
                )
                fidelity_matrix[1, col] = float(
                    load_array(DATA_DIR / "phase4_results.npz", "gray_box_fidelities")[idx]
                )
                fidelity_matrix[2, col] = float(
                    load_array(DATA_DIR / "phase4_results.npz", "perfect_fidelities")[idx]
                )
            # Noisy 30%
            nom_noisy = load_array(DATA_DIR / "phase4_results.npz", "nominal_fidelities_noisy")
            gb_noisy = load_array(DATA_DIR / "phase4_results.npz", "gray_box_fidelities_noisy")
            perf_noisy = load_array(DATA_DIR / "phase4_results.npz", "perfect_fidelities_noisy")
            if nom_noisy is not None:
                idx_30 = np.argmin(np.abs(np.array(mismatches) - 0.30))
                fidelity_matrix[0, 2] = float(nom_noisy[idx_30])
                fidelity_matrix[1, 2] = float(gb_noisy[idx_30])
                fidelity_matrix[2, 2] = float(perf_noisy[idx_30])

    p5_3 = DATA_DIR / "phase5_3_readout_sweep.npz"
    if p5_3.exists():
        d53 = np.load(p5_3)
        c01 = d53["confusion_01_values"]
        idx_5pct = np.argmin(np.abs(c01 - 0.05))
        fidelity_matrix[0, 3] = float(d53["nominal_fidelities"][idx_5pct])
        fidelity_matrix[1, 3] = float(d53["gray_box_fidelities"][idx_5pct])
        fidelity_matrix[2, 3] = float(d53["perfect_fidelities"][idx_5pct])

    p5_4 = DATA_DIR / "phase5_4_probe_budget.npz"
    if p5_4.exists():
        d54 = np.load(p5_4)
        n_shots = d54["n_shots_values"]
        idx_50 = np.argmin(np.abs(n_shots - 50))
        fidelity_matrix[0, 4] = float(d54["nominal_fidelity"])
        fidelity_matrix[1, 4] = float(d54["gray_box_fidelities"][idx_50])
        fidelity_matrix[2, 4] = float(d54["perfect_fidelity"])

    p5_1 = DATA_DIR / "phase5_1_chi_mismatch.npz"
    if p5_1.exists():
        d51 = np.load(p5_1)
        ratios = d51["ratios"]
        idx_70 = np.argmin(np.abs(ratios - 0.30))  # 70% mismatch = ratio 0.30
        fidelity_matrix[0, 5] = float(d51["nominal_fidelities"][idx_70])
        fidelity_matrix[1, 5] = float(d51["gray_box_fidelities"][idx_70])
        fidelity_matrix[2, 5] = float(d51["perfect_fidelities"][idx_70])

    # Viability colors
    def viability_color(f):
        if np.isnan(f):
            return [0.7, 0.7, 0.7]   # gray = no data
        if f > 0.99:
            return [0.2, 0.7, 0.2]   # green
        if f > 0.95:
            return [0.95, 0.85, 0.1] # yellow
        return [0.85, 0.2, 0.2]      # red

    color_matrix = np.zeros((*fidelity_matrix.shape, 3))
    for i in range(fidelity_matrix.shape[0]):
        for j in range(fidelity_matrix.shape[1]):
            color_matrix[i, j] = viability_color(fidelity_matrix[i, j])

    fig, ax = plt.subplots(figsize=(len(regime_labels) * 1.6 + 1, len(method_labels) * 1.5 + 1.5))
    ax.imshow(color_matrix, aspect="auto")

    for i in range(len(method_labels)):
        for j in range(len(regime_labels)):
            f = fidelity_matrix[i, j]
            text = f"{f:.3f}" if not np.isnan(f) else "N/A"
            ax.text(j, i, text, ha="center", va="center", fontsize=11, fontweight="bold",
                    color="black")

    ax.set_xticks(np.arange(len(regime_labels)))
    ax.set_xticklabels(regime_labels, fontsize=10)
    ax.set_yticks(np.arange(len(method_labels)))
    ax.set_yticklabels(method_labels, fontsize=11)
    ax.set_title("Viability Summary: Gate Fidelity by Regime and Method\n"
                 "(Green F>0.99, Yellow 0.95<F≤0.99, Red F≤0.95, Gray=No data)",
                 fontsize=12)
    ax.set_xlabel("Operating Regime")
    ax.set_ylabel("Control Strategy")
    fig.tight_layout()
    save_fig(fig, "viability_summary.pdf")


# ---------------------------------------------------------------------------
# CLI entrypoint
# ---------------------------------------------------------------------------


def main():
    if not HAVE_MPL:
        print("matplotlib not available; cannot generate figures.")
        sys.exit(1)

    parser = argparse.ArgumentParser(description="Generate figures for gray-box adaptive control study")
    parser.add_argument("--figure", type=str, default="all",
                        help="Which figure to plot: phase4, phase5_1, phase5_2, phase5_3, "
                             "phase5_4, phase5_5, phase5_6, viability, or 'all'")
    args = parser.parse_args()

    setup_style()
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    figure_map = {
        "phase4": [plot_phase4_main, plot_phase4_per_fock],
        "phase5_1": [plot_phase5_1],
        "phase5_2": [plot_phase5_2],
        "phase5_3": [plot_phase5_3],
        "phase5_4": [plot_phase5_4],
        "phase5_5": [plot_phase5_5],
        "phase5_6": [plot_phase5_6],
        "viability": [plot_viability_summary],
    }

    if args.figure == "all":
        for key, fns in figure_map.items():
            for fn in fns:
                try:
                    fn()
                except Exception as exc:
                    print(f"  [ERROR] {fn.__name__}: {exc}")
    elif args.figure in figure_map:
        for fn in figure_map[args.figure]:
            try:
                fn()
            except Exception as exc:
                print(f"  [ERROR] {fn.__name__}: {exc}")
    else:
        print(f"Unknown figure '{args.figure}'. Options: {', '.join(list(figure_map.keys()) + ['all'])}")
        sys.exit(1)

    print("Done.")


if __name__ == "__main__":
    main()
