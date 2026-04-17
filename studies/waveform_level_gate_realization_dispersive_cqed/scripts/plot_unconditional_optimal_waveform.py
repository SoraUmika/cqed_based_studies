"""Plot the best hardware-aware optimal-control waveform for the report."""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import common
from common import ARTIFACTS_DIR, FIGURES_DIR, apply_plot_style


def main() -> None:
    apply_plot_style()
    payload = common.load_json(ARTIFACTS_DIR / "unconditional_optimal_control_summary.json")
    best = max(
        payload["cases"],
        key=lambda item: item["full_metrics"]["state_test_mean_fidelity"],
    )

    physical = np.asarray(best["physical_values"], dtype=float)
    duration_ns = float(best["duration_ns"])
    n_steps = physical.shape[1]
    dt_ns = duration_ns / n_steps
    times_ns = (np.arange(n_steps, dtype=float) + 0.5) * dt_ns

    complex_wave = physical[0] + 1j * physical[1]
    freqs_mhz = np.fft.fftshift(np.fft.fftfreq(n_steps, d=dt_ns * 1.0e-9)) / 1.0e6
    spectrum = np.fft.fftshift(np.abs(np.fft.fft(complex_wave)))

    fig, axes = plt.subplots(2, 1, figsize=(6.5, 5.5))
    axes[0].step(times_ns, physical[0] / 1.0e6, where="mid", label="I", linewidth=2.0)
    axes[0].step(times_ns, physical[1] / 1.0e6, where="mid", label="Q", linewidth=2.0)
    axes[0].set_xlabel("Time (ns)")
    axes[0].set_ylabel(r"Amplitude (Mrad/s)")
    axes[0].set_title("Physical optimal-control waveform")
    axes[0].legend()

    axes[1].plot(freqs_mhz, spectrum, linewidth=2.0)
    axes[1].set_xlabel("Frequency (MHz)")
    axes[1].set_ylabel("Arbitrary units")
    axes[1].set_title("Magnitude spectrum")

    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "unconditional_optimal_waveform.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "unconditional_optimal_waveform.pdf", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
