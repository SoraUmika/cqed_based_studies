"""
Plot all-branch multitone SQR study results.

Generates publication-quality figures comparing baseline, optimized multitone,
and GRAPE fidelities for the all-branch (n=0..3) simultaneous R_X(pi) gate
as a function of chi*T/2pi.

Input:
    data/allbranch_multitone_results.npz

Output:
    figures/allbranch_fidelity_comparison.{png,pdf}
    figures/allbranch_branch_breakdown.{png,pdf}
    figures/allbranch_leakage.{png,pdf}
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
style_path = (
    STUDY_DIR.parents[1]
    / ".github"
    / "skills"
    / "publication-figures"
    / "assets"
    / "cqed_style.mplstyle"
)
if style_path.exists():
    plt.style.use(str(style_path))

sys.path.insert(0, str(SCRIPT_DIR))
from common import TOL_BRIGHT

FIG_DIR.mkdir(exist_ok=True)


def main():
    data_path = DATA_DIR / "allbranch_multitone_results.npz"
    if not data_path.exists():
        print(f"ERROR: {data_path} not found. Run run_allbranch_multitone.py first.")
        sys.exit(1)

    d = np.load(str(data_path), allow_pickle=True)
    chi_t = d["chi_t_values"]

    # ---------------------------------------------------------------
    # Figure 1: Fidelity comparison (F_block and F_strict)
    # ---------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5), sharey=False)

    # Block fidelity
    ax1.plot(chi_t, d["baseline_fid_block"], "o-", color=TOL_BRIGHT[0],
             label="Baseline (common env.)", markersize=5)
    ax1.plot(chi_t, d["opt_ap_fid_block"], "s-", color=TOL_BRIGHT[1],
             label="Opt. amp+phase (8p)", markersize=5)
    ax1.plot(chi_t, d["opt_apd_fid_block"], "^-", color=TOL_BRIGHT[2],
             label="Opt. amp+phase+det (12p)", markersize=5)
    ax1.plot(chi_t, d["opt_de_fid_block"], "D-", color=TOL_BRIGHT[3],
             label="DE + NM (12p)", markersize=5)

    # GRAPE if available
    grape_fid = d["grape_fid"]
    grape_valid = ~np.isnan(grape_fid)
    if np.any(grape_valid):
        ax1.plot(chi_t[grape_valid], grape_fid[grape_valid], "*-",
                 color=TOL_BRIGHT[5], label="GRAPE", markersize=8)

    ax1.set_xlabel(r"$\chi T / 2\pi$")
    ax1.set_ylabel(r"$F_\mathrm{block}$")
    ax1.set_title("Block fidelity (per-branch phase freedom)")
    ax1.legend(fontsize=7, loc="lower right")
    ax1.set_ylim(-0.05, 1.05)
    ax1.axhline(0.999, ls="--", color="gray", alpha=0.4, lw=0.8)

    # Strict fidelity
    ax2.plot(chi_t, d["baseline_fid_strict"], "o-", color=TOL_BRIGHT[0],
             label="Baseline", markersize=5)
    ax2.plot(chi_t, d["opt_ap_fid_strict"], "s-", color=TOL_BRIGHT[1],
             label="Opt. amp+phase", markersize=5)
    ax2.plot(chi_t, d["opt_apd_fid_strict"], "^-", color=TOL_BRIGHT[2],
             label="Opt. amp+phase+det", markersize=5)
    ax2.plot(chi_t, d["opt_de_fid_strict"], "D-", color=TOL_BRIGHT[3],
             label="DE + NM", markersize=5)
    if np.any(grape_valid):
        ax2.plot(chi_t[grape_valid], grape_fid[grape_valid], "*-",
                 color=TOL_BRIGHT[5], label="GRAPE", markersize=8)

    ax2.set_xlabel(r"$\chi T / 2\pi$")
    ax2.set_ylabel(r"$F_\mathrm{strict}$")
    ax2.set_title("Strict fidelity (global phase only)")
    ax2.legend(fontsize=7, loc="lower right")
    ax2.set_ylim(-0.05, 1.05)
    ax2.axhline(0.999, ls="--", color="gray", alpha=0.4, lw=0.8)

    fig.tight_layout()
    fig.savefig(str(FIG_DIR / "allbranch_fidelity_comparison.png"),
                dpi=300, bbox_inches="tight")
    fig.savefig(str(FIG_DIR / "allbranch_fidelity_comparison.pdf"),
                bbox_inches="tight")
    plt.close(fig)
    print("Saved allbranch_fidelity_comparison.{png,pdf}")

    # ---------------------------------------------------------------
    # Figure 2: Per-branch fidelity breakdown (best optimizer: DE+NM)
    # ---------------------------------------------------------------
    branch_labels = [f"$n={n}$" for n in range(4)]
    fig, ax = plt.subplots(figsize=(7, 4.5))

    n_branches = d["opt_de_branch_fids"].shape[1]
    x = np.arange(len(chi_t))
    width = 0.18

    for b in range(n_branches):
        ax.bar(x + b * width, d["opt_de_branch_fids"][:, b],
               width, label=branch_labels[b], color=TOL_BRIGHT[b])

    ax.set_xticks(x + 1.5 * width)
    ax.set_xticklabels([f"{ct}" for ct in chi_t])
    ax.set_xlabel(r"$\chi T / 2\pi$")
    ax.set_ylabel(r"Per-branch $F_\mathrm{block}$")
    ax.set_title("DE+NM per-branch fidelity (all-branch target)")
    ax.legend(fontsize=8)
    ax.set_ylim(0, 1.05)

    fig.tight_layout()
    fig.savefig(str(FIG_DIR / "allbranch_branch_breakdown.png"),
                dpi=300, bbox_inches="tight")
    fig.savefig(str(FIG_DIR / "allbranch_branch_breakdown.pdf"),
                bbox_inches="tight")
    plt.close(fig)
    print("Saved allbranch_branch_breakdown.{png,pdf}")

    # ---------------------------------------------------------------
    # Figure 3: Leakage comparison
    # ---------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(6, 4))

    ax.plot(chi_t, d["baseline_leakage"], "o-", color=TOL_BRIGHT[0],
            label="Baseline", markersize=5)
    ax.plot(chi_t, d["opt_ap_leakage"], "s-", color=TOL_BRIGHT[1],
            label="Opt. amp+phase", markersize=5)
    ax.plot(chi_t, d["opt_apd_leakage"], "^-", color=TOL_BRIGHT[2],
            label="Opt. amp+phase+det", markersize=5)
    ax.plot(chi_t, d["opt_de_leakage"], "D-", color=TOL_BRIGHT[3],
            label="DE + NM", markersize=5)

    ax.set_xlabel(r"$\chi T / 2\pi$")
    ax.set_ylabel("Leakage")
    ax.set_yscale("log")
    ax.legend(fontsize=8)
    ax.set_title("Leakage to non-computational subspace")

    fig.tight_layout()
    fig.savefig(str(FIG_DIR / "allbranch_leakage.png"),
                dpi=300, bbox_inches="tight")
    fig.savefig(str(FIG_DIR / "allbranch_leakage.pdf"),
                bbox_inches="tight")
    plt.close(fig)
    print("Saved allbranch_leakage.{png,pdf}")

    # ---------------------------------------------------------------
    # Summary table to stdout
    # ---------------------------------------------------------------
    print("\n=== All-branch summary ===")
    print(f"{'chiT/2pi':>8} | {'Base_Fb':>8} {'OptAP_Fb':>9} {'OptAPD_Fb':>10} "
          f"{'DE_Fb':>8} {'DE_Fs':>8} {'GRAPE':>8}")
    print("-" * 75)
    for i, ct in enumerate(chi_t):
        grape_str = f"{grape_fid[i]:.4f}" if not np.isnan(grape_fid[i]) else "  N/A  "
        print(f"{ct:>8.1f} | {d['baseline_fid_block'][i]:>8.4f} "
              f"{d['opt_ap_fid_block'][i]:>9.4f} "
              f"{d['opt_apd_fid_block'][i]:>10.4f} "
              f"{d['opt_de_fid_block'][i]:>8.4f} "
              f"{d['opt_de_fid_strict'][i]:>8.4f} "
              f"{grape_str:>8}")

    print("\nDone.")


if __name__ == "__main__":
    main()
