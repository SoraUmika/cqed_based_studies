"""
Plot pulse waveforms for all families and the HW-constrained GRAPE pulse.

Reads data/pulse_waveforms.npz (produced by run_hwgrape_and_waveforms.py)
and generates:
  - figures/pulse_waveforms_families.{png,pdf}  — envelope comparison
  - figures/pulse_waveform_grape.{png,pdf}       — GRAPE I/Q waveform

Usage:
    python scripts/plot_pulse_waveforms.py
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

# Publication style
style_path = SCRIPT_DIR.parents[1] / ".github" / "skills" / "publication-figures" / "assets" / "cqed_style.mplstyle"
if style_path.exists():
    plt.style.use(str(style_path))

TOL_BRIGHT = ['#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB']


def to_mhz(arr_rad_s):
    """Convert rad/s amplitude to MHz."""
    return arr_rad_s / (2 * np.pi * 1e6)


def plot_family_waveforms(d):
    """Plot I/Q envelopes for all four parametric families side by side."""
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 4.5), sharex=True)
    chi_t = float(d["waveform_chi_t"])
    dur_ns = float(d["duration_ns"])

    families = [
        ("Gaussian", d["gauss_t_ns"], d["gauss_I"], d["gauss_Q"]),
        ("Square", d["sq_t_ns"], d["sq_I"], d["sq_Q"]),
        ("Cosine-squared", d["cos2_t_ns"], d["cos2_I"], d["cos2_Q"]),
        ("Multitone", d["mt_t_ns"], d["mt_I"], d["mt_Q"]),
    ]

    for idx, (name, t, I, Q) in enumerate(families):
        ax = axes.flat[idx]
        amp = np.sqrt(I**2 + Q**2)
        ax.plot(t, to_mhz(I), color=TOL_BRIGHT[0], lw=0.8, label="I", alpha=0.85)
        ax.plot(t, to_mhz(Q), color=TOL_BRIGHT[1], lw=0.8, label="Q", alpha=0.85)
        ax.plot(t, to_mhz(amp), color=TOL_BRIGHT[2], lw=1.0, ls="--",
                label=r"$|\varepsilon|$", alpha=0.9)
        ax.set_title(name, fontsize=9)
        if idx >= 2:
            ax.set_xlabel("Time (ns)")
        if idx % 2 == 0:
            ax.set_ylabel(r"Amplitude (MHz)")
        if idx == 0:
            ax.legend(fontsize=7, loc="upper right", framealpha=0.8)
        ax.set_xlim(0, dur_ns)

    fig.suptitle(
        rf"Pulse waveforms at $\chi T/2\pi = {chi_t:.0f}$ "
        rf"($T = {dur_ns:.0f}$ ns)",
        fontsize=10, y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    for ext in ("png", "pdf"):
        fig.savefig(str(FIG_DIR / f"pulse_waveforms_families.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved pulse_waveforms_families.{png,pdf}")


def plot_grape_waveform(d):
    """Plot the HW-constrained GRAPE I/Q pulse waveform."""
    t = d["grape_t_ns"]
    I = d["grape_I"]
    Q = d["grape_Q"]
    if len(t) == 0:
        print("No GRAPE waveform data — skipping grape plot.")
        return

    chi_t = float(d["waveform_chi_t"])
    dur_ns = float(t[-1] + t[1] - t[0]) if len(t) > 1 else float(d["duration_ns"])

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(3.375, 3.5), sharex=True)

    # I/Q channels
    ax1.step(t, to_mhz(I), where="post", color=TOL_BRIGHT[0], lw=0.8, label="I")
    ax1.step(t, to_mhz(Q), where="post", color=TOL_BRIGHT[1], lw=0.8, label="Q")
    ax1.set_ylabel(r"Amplitude (MHz)")
    ax1.legend(fontsize=7, loc="upper right")
    ax1.set_title(
        rf"HW-constrained GRAPE pulse ($\chi T/2\pi = {chi_t:.0f}$)",
        fontsize=9)

    # Amplitude envelope
    amp = np.sqrt(I**2 + Q**2)
    ax2.step(t, to_mhz(amp), where="post", color=TOL_BRIGHT[2], lw=1.0)
    ax2.axhline(100.0, color="gray", ls=":", lw=0.6, label="100 MHz limit")
    ax2.set_ylabel(r"$|\varepsilon|$ (MHz)")
    ax2.set_xlabel("Time (ns)")
    ax2.legend(fontsize=7, loc="upper right")
    ax2.set_xlim(0, dur_ns)

    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(str(FIG_DIR / f"pulse_waveform_grape.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved pulse_waveform_grape.{png,pdf}")


def plot_grape_fidelity_comparison():
    """Plot HW-constrained vs unconstrained GRAPE fidelity comparison."""
    fpath = DATA_DIR / "hwgrape_results.npz"
    if not fpath.exists():
        print("No hwgrape_results.npz — skipping fidelity comparison.")
        return

    d = np.load(str(fpath))
    ct = d["grape_chi_t"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 2.8), sharey=True)

    # Left: cphase SQR
    ax1.plot(ct, d["hw_fid_cphase"], "o-", color=TOL_BRIGHT[0], ms=4,
             label="HW-constrained", lw=1.2)
    ax1.plot(ct, d["unc_fid_cphase"], "s--", color=TOL_BRIGHT[2], ms=4,
             label="Unconstrained", lw=1.2)
    ax1.axhline(0.999, color="gray", ls=":", lw=0.5)
    ax1.set_xlabel(r"$\chi T / 2\pi$")
    ax1.set_ylabel("Fidelity")
    ax1.set_title("Cphase SQR (GRAPE)", fontsize=9)
    ax1.legend(fontsize=7)
    ax1.set_ylim(0.78, 1.005)

    # Right: true SQR
    ax2.plot(ct, d["hw_fid_true"], "o-", color=TOL_BRIGHT[0], ms=4,
             label="HW-constrained", lw=1.2)
    ax2.plot(ct, d["unc_fid_true"], "s--", color=TOL_BRIGHT[2], ms=4,
             label="Unconstrained", lw=1.2)
    ax2.axhline(0.999, color="gray", ls=":", lw=0.5)
    ax2.set_xlabel(r"$\chi T / 2\pi$")
    ax2.set_title("True SQR (GRAPE)", fontsize=9)
    ax2.legend(fontsize=7)

    fig.suptitle(
        "GRAPE fidelity: HW-constrained (1 GHz AWG, 100 MHz) vs unconstrained",
        fontsize=9, y=1.01)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(str(FIG_DIR / f"grape_hw_vs_unconstrained.{ext}"),
                    dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved grape_hw_vs_unconstrained.{png,pdf}")


def main():
    FIG_DIR.mkdir(exist_ok=True)

    # Plot pulse waveforms
    wf_path = DATA_DIR / "pulse_waveforms.npz"
    if wf_path.exists():
        d = np.load(str(wf_path))
        plot_family_waveforms(d)
        plot_grape_waveform(d)
    else:
        print(f"Missing {wf_path} — run run_hwgrape_and_waveforms.py first.")

    # Plot GRAPE comparison
    plot_grape_fidelity_comparison()

    print("All plots done.")


if __name__ == "__main__":
    main()
