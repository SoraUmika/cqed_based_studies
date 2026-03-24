"""
Plotting for the SQR phase-compilation follow-up study.

Reads ``data/phase_compilation_results.npz`` and generates dedicated figures
for the report:

- single-target strict-vs-compiled recovery,
- extracted cavity-phase profiles and linear fits,
- coherent superposition benchmark improvement,
- all-branch short-gate comparison.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from common import TOL_BRIGHT

SCRIPT_DIR = Path(__file__).resolve().parent
STUDY_DIR = SCRIPT_DIR.parent
DATA_PATH = STUDY_DIR / "data" / "phase_compilation_results.npz"
FIG_DIR = STUDY_DIR / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)


def save_figure(fig, stem: str) -> None:
    png_path = FIG_DIR / f"{stem}.png"
    pdf_path = FIG_DIR / f"{stem}.pdf"
    fig.savefig(png_path, dpi=220, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    print(f"Saved {png_path.name} and {pdf_path.name}")


def family_label(name: str) -> str:
    return {
        "single_tone_gaussian": "Gaussian",
        "cosine_squared": "Cosine-squared",
    }.get(name, name.replace("_", " "))


def plot_single_target_recovery(data) -> None:
    chi_t = np.asarray(data["single_chi_t_values"], dtype=float)
    family_names = [str(name) for name in data["single_family_names"]]
    raw = np.asarray(data["single_raw_global_z_fid"], dtype=float)
    compiled = np.asarray(data["single_linear_cavity_compiled_fid"], dtype=float)
    branchz = np.asarray(data["single_branch_local_z_relaxed_fid"], dtype=float)

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0), sharey=True)
    for family_index, ax in enumerate(axes):
        ax.plot(chi_t, raw[family_index], "o-", lw=2.0, color=TOL_BRIGHT[0], label="Raw strict + global Z")
        ax.plot(
            chi_t,
            compiled[family_index],
            "s-",
            lw=2.0,
            color=TOL_BRIGHT[2],
            label="Linear cavity phase + global Z",
        )
        ax.plot(
            chi_t,
            branchz[family_index],
            "^-",
            lw=2.0,
            color=TOL_BRIGHT[1],
            label="Per-branch local Z upper bound",
        )
        ax.set_title(family_label(family_names[family_index]))
        ax.set_xlabel(r"$\chi T / 2\pi$")
        ax.set_ylim(0.58, 1.01)
        ax.grid(alpha=0.25, linestyle=":")
    axes[0].set_ylabel("Logical fidelity")
    axes[0].legend(loc="lower right", fontsize=9)
    fig.suptitle("Single-target SQR on enlarged logical window (N = 8)", y=1.02, fontsize=13)
    save_figure(fig, "phase_compilation_single_target_recovery")
    plt.close(fig)


def plot_phase_profiles(data) -> None:
    family_names = [str(name) for name in data["single_family_names"]]
    chi_t = np.asarray(data["single_chi_t_values"], dtype=float)
    exact = np.asarray(data["single_exact_cavity_phases"], dtype=float)
    linear = np.asarray(data["single_linear_phase_profile"], dtype=float)
    n = np.arange(exact.shape[-1], dtype=int)
    rep_index = int(np.where(np.isclose(chi_t, 3.0))[0][0])

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0), sharey=True)
    for family_index, ax in enumerate(axes):
        ax.plot(
            n,
            exact[family_index, rep_index],
            "o",
            ms=6,
            color=TOL_BRIGHT[0],
            label="Exact fitted cavity phase",
        )
        ax.plot(
            n,
            linear[family_index, rep_index],
            "-",
            lw=2.0,
            color=TOL_BRIGHT[2],
            label="Linear fit",
        )
        ax.set_title(f"{family_label(family_names[family_index])}, " + r"$\chi T/2\pi = 3$")
        ax.set_xlabel("Fock level n")
        ax.grid(alpha=0.25, linestyle=":")
    axes[0].set_ylabel("Gauge-fixed phase (rad)")
    axes[0].legend(loc="upper right", fontsize=9)
    fig.suptitle("Extracted cavity-only phase profiles are essentially linear in n", y=1.02, fontsize=13)
    save_figure(fig, "phase_compilation_phase_profiles")
    plt.close(fig)


def plot_superposition_benchmarks(data) -> None:
    family_names = [str(name) for name in data["single_family_names"]]
    chi_t = np.asarray(data["single_chi_t_values"], dtype=float)
    rep_index = int(np.where(np.isclose(chi_t, 3.0))[0][0])
    raw_mean = np.asarray(data["single_pair_superposition_raw_mean"], dtype=float)[:, rep_index]
    raw_min = np.asarray(data["single_pair_superposition_raw_min"], dtype=float)[:, rep_index]
    compiled_mean = np.asarray(data["single_pair_superposition_compiled_mean"], dtype=float)[:, rep_index]
    compiled_min = np.asarray(data["single_pair_superposition_compiled_min"], dtype=float)[:, rep_index]

    x = np.arange(len(family_names), dtype=float)
    width = 0.18
    fig, ax = plt.subplots(figsize=(8.0, 4.2))
    ax.bar(x - 1.5 * width, raw_mean, width, color=TOL_BRIGHT[0], label="Raw mean")
    ax.bar(x - 0.5 * width, compiled_mean, width, color=TOL_BRIGHT[2], label="Compiled mean")
    ax.bar(x + 0.5 * width, raw_min, width, color=TOL_BRIGHT[4], label="Raw min")
    ax.bar(x + 1.5 * width, compiled_min, width, color=TOL_BRIGHT[1], label="Compiled min")
    ax.set_xticks(x, [family_label(name) for name in family_names])
    ax.set_ylabel("Pair-superposition fidelity")
    ax.set_ylim(0.6, 1.01)
    ax.grid(alpha=0.25, linestyle=":", axis="y")
    ax.legend(ncol=2, fontsize=9)
    ax.set_title(r"Coherence benchmark over $|g,n\rangle \pm |g,m\rangle$, $|e,n\rangle \pm |e,m\rangle$, and $|g,n\rangle \pm |e,m\rangle$ probes")
    save_figure(fig, "phase_compilation_superposition_benchmarks")
    plt.close(fig)


def plot_allbranch_short_gate(data) -> None:
    chi_t = np.asarray(data["allbranch_chi_t_values"], dtype=float)
    raw = np.asarray(data["allbranch_raw_global_z_fid"], dtype=float)
    cavity = np.asarray(data["allbranch_exact_cavity_compiled_fid"], dtype=float)
    branchz = np.asarray(data["allbranch_branch_local_z_relaxed_fid"], dtype=float)
    grape = np.asarray(data["allbranch_grape_reference_fid"], dtype=float)

    x = np.arange(chi_t.size, dtype=float)
    width = 0.2
    fig, ax = plt.subplots(figsize=(8.0, 4.2))
    ax.bar(x - 1.5 * width, raw, width, color=TOL_BRIGHT[0], label="Raw strict + global Z")
    ax.bar(x - 0.5 * width, cavity, width, color=TOL_BRIGHT[2], label="Cavity-only compiled")
    ax.bar(x + 0.5 * width, branchz, width, color=TOL_BRIGHT[1], label="Per-branch local Z")
    ax.bar(x + 1.5 * width, grape, width, color=TOL_BRIGHT[3], label="Saved GRAPE reference")
    ax.set_xticks(x, [f"{value:.1f}" for value in chi_t])
    ax.set_xlabel(r"$\chi T / 2\pi$")
    ax.set_ylabel("Logical fidelity")
    ax.set_ylim(0.0, 1.05)
    ax.grid(alpha=0.25, linestyle=":", axis="y")
    ax.legend(fontsize=9)
    ax.set_title("All-branch short-gate structured multitone: cavity-only compilation is insufficient")
    save_figure(fig, "phase_compilation_allbranch_short_gate")
    plt.close(fig)


def main() -> None:
    if not DATA_PATH.exists():
        raise FileNotFoundError(f"Missing data file: {DATA_PATH}")

    data = np.load(DATA_PATH, allow_pickle=True)
    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.titlesize": 12,
            "axes.labelsize": 11,
            "legend.fontsize": 10,
            "figure.titlesize": 13,
        }
    )

    plot_single_target_recovery(data)
    plot_phase_profiles(data)
    plot_superposition_benchmarks(data)
    plot_allbranch_short_gate(data)


if __name__ == "__main__":
    main()
