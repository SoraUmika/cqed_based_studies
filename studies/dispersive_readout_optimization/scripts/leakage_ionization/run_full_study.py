"""Run the full measurement-induced leakage and ionization study workflow."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import time

import matplotlib.pyplot as plt
import numpy as np

try:
    from .common import (
        DATA_DIR,
        FIGURES_DIR,
        STUDY_DIR,
        TWO_PI,
        StudyConfig,
        build_protocol,
        ensure_directories,
        run_continuous_replay,
        save_json,
        simulate_point,
    )
except ImportError:  # pragma: no cover - direct script execution path
    from common import (
        DATA_DIR,
        FIGURES_DIR,
        STUDY_DIR,
        TWO_PI,
        StudyConfig,
        build_protocol,
        ensure_directories,
        run_continuous_replay,
        save_json,
        simulate_point,
    )

LEAK_THRESHOLD = 1.0e-4


def save_figure(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIGURES_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def classify_records(g_records: np.ndarray, e_records: np.ndarray) -> dict[str, float]:
    mean_g = float(np.mean(g_records))
    mean_e = float(np.mean(e_records))
    threshold = 0.5 * (mean_g + mean_e)
    if mean_g <= mean_e:
        correct_g = float(np.mean(g_records <= threshold))
        correct_e = float(np.mean(e_records >= threshold))
    else:
        correct_g = float(np.mean(g_records >= threshold))
        correct_e = float(np.mean(e_records <= threshold))
    fidelity = 0.5 * (correct_g + correct_e)

    def skewness(values: np.ndarray) -> float:
        centered = values - np.mean(values)
        sigma = np.std(centered)
        if sigma <= 1.0e-15:
            return 0.0
        return float(np.mean((centered / sigma) ** 3))

    return {
        "mean_g": mean_g,
        "mean_e": mean_e,
        "threshold": threshold,
        "correct_g": correct_g,
        "correct_e": correct_e,
        "assignment_fidelity": fidelity,
        "skew_g": skewness(g_records),
        "skew_e": skewness(e_records),
        "var_g": float(np.var(g_records)),
        "var_e": float(np.var(e_records)),
    }


def sample_histogram_proxy(
    protocol,
    record: dict[str, object],
    *,
    shots: int = 4096,
    seed: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    trace_g = protocol.readout_chain.simulate_waveform(
        "g",
        protocol.distorted_waveform,
        dt=protocol.readout_chain.dt,
        include_noise=False,
    )
    trace_e = protocol.readout_chain.simulate_waveform(
        "e",
        protocol.distorted_waveform,
        dt=protocol.readout_chain.dt,
        include_noise=False,
    )
    center_g = np.asarray(trace_g.iq_sample, dtype=float)
    center_e = np.asarray(trace_e.iq_sample, dtype=float)
    axis = center_e - center_g
    norm = np.linalg.norm(axis)
    if norm <= 1.0e-15:
        axis = np.array([1.0, 0.0], dtype=float)
    else:
        axis = axis / norm
    sigma = protocol.readout_chain.integrated_noise_sigma(duration=protocol.duration, dt=protocol.readout_chain.dt)
    rng = np.random.default_rng(seed)

    def draw(state_payload: dict[str, float]) -> np.ndarray:
        p_g = float(state_payload["p_g"])
        p_e_like = float(max(0.0, 1.0 - p_g))
        latent = rng.choice(np.array([0, 1], dtype=int), size=shots, p=[p_g, p_e_like])
        centers = np.where(latent[:, None] == 0, center_g[None, :], center_e[None, :])
        iq = centers + rng.normal(scale=sigma, size=(shots, 2))
        return np.asarray(iq @ axis, dtype=float)

    return draw(record["g"]), draw(record["e"])


def choose_representative_indices(records: list[dict[str, object]]) -> tuple[int, int, int]:
    leakage = np.asarray([float(record["mean_p_leak"]) for record in records], dtype=float)
    above = np.flatnonzero(leakage > LEAK_THRESHOLD)
    if above.size == 0:
        return 0, len(records) // 2, len(records) - 1
    near = int(above[0])
    below = max(0, near - 1)
    higher = min(len(records) - 1, near + 1)
    return below, near, higher


def threshold_from_scan(records: list[dict[str, object]]) -> float:
    for record in records:
        if float(record["mean_p_leak"]) > LEAK_THRESHOLD:
            return float(record["amplitude_mhz"])
    return float(records[-1]["amplitude_mhz"])


def plot_regime_map(cfg: StudyConfig, sweep: list[dict[str, object]]) -> None:
    durations = sorted({float(record["duration_ns"]) for record in sweep})
    amplitudes = sorted({float(record["amplitude_mhz"]) for record in sweep})
    leak = np.zeros((len(durations), len(amplitudes)), dtype=float)
    ion = np.zeros_like(leak)
    qnd = np.zeros_like(leak)
    for record in sweep:
        i = durations.index(float(record["duration_ns"]))
        j = amplitudes.index(float(record["amplitude_mhz"]))
        leak[i, j] = float(record["mean_p_leak"])
        ion[i, j] = float(record["mean_p_ion"])
        qnd[i, j] = float(record["qnd_defect"])

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.8), constrained_layout=True)
    titles = [
        r"$\log_{10} P_{\mathrm{leak}}$",
        r"$\log_{10} P_{\mathrm{ion}}$",
        r"$Q_{\mathrm{QND}}$ defect",
    ]
    panels = [
        np.log10(np.maximum(leak, 1.0e-6)),
        np.log10(np.maximum(ion, 1.0e-7)),
        qnd,
    ]
    cmaps = ["viridis", "magma", "cividis"]
    for ax, panel, title, cmap in zip(axes, panels, titles, cmaps, strict=True):
        image = ax.imshow(panel, origin="lower", aspect="auto", cmap=cmap)
        ax.set_xticks(range(len(amplitudes)), [f"{value:.1f}" for value in amplitudes], rotation=0)
        ax.set_yticks(range(len(durations)), [f"{value:.0f}" for value in durations])
        ax.set_xlabel(r"Amplitude $|\epsilon|/2\pi$ (MHz)")
        ax.set_ylabel("Duration (ns)")
        ax.set_title(title)
        fig.colorbar(image, ax=ax, shrink=0.85)
    save_figure(fig, "regime_maps")


def plot_detuning_scan(detuning_records: list[dict[str, object]]) -> None:
    detuning = np.asarray([record["detuning_mhz"] for record in detuning_records], dtype=float)
    leak = np.asarray([record["mean_p_leak"] for record in detuning_records], dtype=float)
    qnd = np.asarray([record["qnd_defect"] for record in detuning_records], dtype=float)
    assign = np.asarray([record["assignment_proxy"] for record in detuning_records], dtype=float)

    fig, ax = plt.subplots(figsize=(6.5, 4.0), constrained_layout=True)
    ax.plot(detuning, leak, marker="o", label=r"$P_{\mathrm{leak}}$")
    ax.plot(detuning, qnd, marker="s", label=r"$1-Q_{\mathrm{QND}}$")
    ax.plot(detuning, 1.0 - assign, marker="^", label=r"$1-F_{\mathrm{assign}}$")
    ax.set_xlabel("Drive detuning from cavity (MHz)")
    ax.set_ylabel("Metric value")
    ax.set_title("Detuning Dependence Near Threshold")
    ax.legend()
    save_figure(fig, "detuning_dependence")


def plot_bandwidth_scan(bandwidth_records: list[dict[str, object]]) -> None:
    bw = np.asarray([record["lowpass_bw_mhz"] for record in bandwidth_records], dtype=float)
    leak = np.asarray([record["mean_p_leak"] for record in bandwidth_records], dtype=float)
    ion = np.asarray([record["mean_p_ion"] for record in bandwidth_records], dtype=float)
    rms = np.asarray([record["hardware_rms_error"] for record in bandwidth_records], dtype=float)

    fig, ax1 = plt.subplots(figsize=(6.5, 4.0), constrained_layout=True)
    ax1.plot(bw, leak, marker="o", label=r"$P_{\mathrm{leak}}$", color="tab:blue")
    ax1.plot(bw, ion, marker="s", label=r"$P_{\mathrm{ion}}$", color="tab:red")
    ax1.set_xlabel("Low-pass bandwidth (MHz)")
    ax1.set_ylabel("Population metric")
    ax1.set_title("Filter-Chain Sensitivity")
    ax2 = ax1.twinx()
    ax2.plot(bw, rms, marker="^", linestyle="--", color="tab:green", label="Waveform RMS error")
    ax2.set_ylabel("Transport RMS error")
    handles_1, labels_1 = ax1.get_legend_handles_labels()
    handles_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(handles_1 + handles_2, labels_1 + labels_2, loc="upper right")
    save_figure(fig, "bandwidth_sensitivity")


def plot_mitigation_tradeoff(mitigation_records: list[dict[str, object]]) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 4.2), constrained_layout=True)
    for record in mitigation_records:
        ax.scatter(
            float(record["assignment_proxy"]),
            float(record["mean_p_leak"]),
            s=120.0 * max(float(record["mean_residual_n_c"]), 0.02),
            label=f"{record['shape']} @ {record['amplitude_mhz']:.1f} MHz",
        )
    ax.set_xlabel(r"Assignment proxy $F_{\mathrm{assign}}$")
    ax.set_ylabel(r"$P_{\mathrm{leak}}$")
    ax.set_title("Mitigation Tradeoff at Fixed-Fidelity Search")
    ax.legend(fontsize=8)
    save_figure(fig, "mitigation_tradeoff")


def plot_representative_histograms(representative: dict[str, dict[str, object]]) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 3.8), constrained_layout=True, sharey=True)
    for ax, (label, payload) in zip(axes, representative.items(), strict=True):
        g = np.asarray(payload["g_records"], dtype=float)
        e = np.asarray(payload["e_records"], dtype=float)
        bins = np.linspace(min(g.min(), e.min()), max(g.max(), e.max()), 25)
        ax.hist(g, bins=bins, alpha=0.6, label="prep |g>")
        ax.hist(e, bins=bins, alpha=0.6, label="prep |e>")
        ax.axvline(float(payload["classification"]["threshold"]), color="black", linestyle="--", linewidth=1.0)
        ax.set_title(
            f"{label}: {payload['record']['amplitude_mhz']:.1f} MHz\n"
            f"F={payload['classification']['assignment_fidelity']:.3f}"
        )
        ax.set_xlabel("Integrated homodyne record")
    axes[0].set_ylabel("Counts")
    axes[0].legend()
    save_figure(fig, "representative_histograms")


def plot_convergence(convergence_records: list[dict[str, object]]) -> None:
    labels = [record["label"] for record in convergence_records]
    leak = np.asarray([record["metrics"]["mean_p_leak"] for record in convergence_records], dtype=float)
    ion = np.asarray([record["metrics"]["mean_p_ion"] for record in convergence_records], dtype=float)
    qnd = np.asarray([record["metrics"]["qnd_defect"] for record in convergence_records], dtype=float)

    fig, ax = plt.subplots(figsize=(7.2, 4.0), constrained_layout=True)
    ax.plot(labels, leak, marker="o", label=r"$P_{\mathrm{leak}}$")
    ax.plot(labels, ion, marker="s", label=r"$P_{\mathrm{ion}}$")
    ax.plot(labels, qnd, marker="^", label=r"$1-Q_{\mathrm{QND}}$")
    ax.set_ylabel("Metric value")
    ax.set_title("Convergence Checks Near Threshold")
    ax.legend()
    save_figure(fig, "convergence_checks")


def main() -> None:
    start = time.perf_counter()
    ensure_directories()
    plt.style.use("seaborn-v0_8-whitegrid")
    cfg = StudyConfig()
    print("Running study from", STUDY_DIR)

    regime_sweep: list[dict[str, object]] = []
    for duration in cfg.duration_values:
        for amplitude in cfg.amplitude_values:
            protocol = build_protocol(
                cfg,
                amplitude=float(amplitude),
                duration=float(duration),
                detuning=0.0,
                shape=cfg.baseline_shape,
                lowpass_bw_hz=float(cfg.baseline_bandwidth_mhz * 1.0e6),
            )
            regime_sweep.append(simulate_point(protocol, cfg))
            print(
                "regime point",
                f"T={duration * 1.0e9:.0f} ns",
                f"A={amplitude / (TWO_PI * 1.0e6):.1f} MHz",
                f"leak={regime_sweep[-1]['mean_p_leak']:.3e}",
            )

    representative_duration_ns = float(cfg.duration_values_ns[1])
    representative_records = [
        record for record in regime_sweep if abs(float(record["duration_ns"]) - representative_duration_ns) < 1.0e-9
    ]
    representative_records.sort(key=lambda record: float(record["amplitude_mhz"]))
    threshold_amplitude_mhz = threshold_from_scan(representative_records)
    representative_indices = choose_representative_indices(representative_records)
    near_record = representative_records[representative_indices[1]]

    detuning_records: list[dict[str, object]] = []
    for detuning in cfg.detuning_values:
        protocol = build_protocol(
            cfg,
            amplitude=float(near_record["amplitude"]),
            duration=float(near_record["duration"]),
            detuning=float(detuning),
            shape=cfg.baseline_shape,
            lowpass_bw_hz=float(cfg.baseline_bandwidth_mhz * 1.0e6),
        )
        detuning_records.append(simulate_point(protocol, cfg))

    bandwidth_records: list[dict[str, object]] = []
    for bandwidth in cfg.bandwidth_values:
        protocol = build_protocol(
            cfg,
            amplitude=float(near_record["amplitude"]),
            duration=float(near_record["duration"]),
            detuning=0.0,
            shape=cfg.baseline_shape,
            lowpass_bw_hz=float(bandwidth),
        )
        bandwidth_records.append(simulate_point(protocol, cfg))

    mitigation_candidates: list[dict[str, object]] = []
    mitigation_records: list[dict[str, object]] = []
    for shape in cfg.mitigation_shapes:
        shape_records = []
        for amplitude in cfg.amplitude_values:
            protocol = build_protocol(
                cfg,
                amplitude=float(amplitude),
                duration=float(near_record["duration"]),
                detuning=0.0,
                shape=shape,
                lowpass_bw_hz=float(cfg.baseline_bandwidth_mhz * 1.0e6),
            )
            metrics = simulate_point(protocol, cfg)
            shape_records.append(metrics)
            mitigation_candidates.append(metrics)
        best = min(shape_records, key=lambda record: abs(float(record["assignment_proxy"]) - cfg.target_assignment))
        mitigation_records.append(best)

    convergence_cases = [
        ("baseline", cfg, None),
        ("more_transmon_levels", replace(cfg, n_tr=7), None),
        ("more_cavity_levels", replace(cfg, n_cav=24), None),
        ("smaller_solver_step", cfg, cfg.dt / 2.0),
    ]
    convergence_records: list[dict[str, object]] = []
    for label, variant_cfg, max_step in convergence_cases:
        protocol = build_protocol(
            variant_cfg,
            amplitude=float(near_record["amplitude"]),
            duration=float(near_record["duration"]),
            detuning=0.0,
            shape=variant_cfg.baseline_shape,
            lowpass_bw_hz=float(variant_cfg.baseline_bandwidth_mhz * 1.0e6),
        )
        convergence_records.append(
            {
                "label": label,
                "config": variant_cfg.as_dict(),
                "metrics": simulate_point(protocol, variant_cfg, max_step=max_step),
            }
        )

    representative_payload: dict[str, dict[str, object]] = {}
    stochastic_cfg = replace(
        cfg,
        mixing_ge_scale=0.0,
        mixing_ef_scale=0.0,
        mixing_slew_ge_scale=0.0,
        mixing_slew_ef_scale=0.0,
    )
    representative_labels = ["below", "near", "above"]
    for label, index in zip(representative_labels, representative_indices, strict=True):
        record = representative_records[index]
        protocol = build_protocol(
            stochastic_cfg,
            amplitude=float(record["amplitude"]),
            duration=float(record["duration"]),
            detuning=0.0,
            shape=stochastic_cfg.baseline_shape,
            lowpass_bw_hz=float(stochastic_cfg.baseline_bandwidth_mhz * 1.0e6),
        )
        replay_mode = "continuous_replay"
        note = "Representative measurement-record histograms were generated with the cavity-monitored SME replay."
        g_records = run_continuous_replay(protocol, stochastic_cfg, initial_level=0)
        e_records = run_continuous_replay(protocol, stochastic_cfg, initial_level=1)
        if np.isnan(g_records).any() or np.isnan(e_records).any():
            replay_mode = "mixture_proxy"
            note = (
                "Native continuous replay became numerically unstable at this operating point, so the "
                "histogram was generated from a deterministic final-state mixture model using the same "
                "readout-chain waveform and integrated-noise scale."
            )
            g_records, e_records = sample_histogram_proxy(protocol, record, shots=4096, seed=cfg.seed + index)
        representative_payload[label] = {
            "record": record,
            "replay_mode": replay_mode,
            "stochastic_note": note,
            "g_records": g_records.tolist(),
            "e_records": e_records.tolist(),
            "classification": classify_records(g_records, e_records),
        }

    amplitude_curve = np.asarray([record["amplitude_mhz"] for record in representative_records], dtype=float)
    leak_curve = np.asarray([record["mean_p_leak"] for record in representative_records], dtype=float)
    occupancy_curve = np.asarray([record["peak_mean_occupancy"] for record in representative_records], dtype=float)
    monotonic_leak_fraction = float(np.mean(np.diff(leak_curve) >= -5.0e-4))
    monotonic_occupancy_fraction = float(np.mean(np.diff(occupancy_curve) >= -5.0e-4))

    baseline_record = representative_records[0]
    validation = {
        "sanity": {
            "low_power_mean_leak": float(baseline_record["mean_p_leak"]),
            "low_power_qnd_consistency": float(baseline_record["qnd_consistency"]),
            "amplitude_monotonic_leak_fraction": monotonic_leak_fraction,
            "amplitude_monotonic_occupancy_fraction": monotonic_occupancy_fraction,
        },
        "convergence": {
            "max_leak_shift": float(
                max(
                    abs(record["metrics"]["mean_p_leak"] - convergence_records[0]["metrics"]["mean_p_leak"])
                    for record in convergence_records[1:]
                )
            ),
            "max_ion_shift": float(
                max(
                    abs(record["metrics"]["mean_p_ion"] - convergence_records[0]["metrics"]["mean_p_ion"])
                    for record in convergence_records[1:]
                )
            ),
        },
        "literature": {
            "applicable": False,
            "note": "This ANA/DES/OPT study is not a direct literature reproduction, so validation is against internal sanity and convergence criteria rather than a paper benchmark.",
        },
    }

    summary = {
        "config": cfg.as_dict(),
        "thresholds": {
            "representative_duration_ns": representative_duration_ns,
            "representative_threshold_amplitude_mhz": threshold_amplitude_mhz,
        },
        "validation": validation,
        "regime_sweep": regime_sweep,
        "detuning_scan": detuning_records,
        "bandwidth_scan": bandwidth_records,
        "mitigation_scan": mitigation_records,
        "convergence": convergence_records,
        "representative_histograms": representative_payload,
        "wall_clock_s": float(time.perf_counter() - start),
    }

    save_json(DATA_DIR / "summary.json", summary)
    np.savez(
        DATA_DIR / "representative_histograms.npz",
        below_g=np.asarray(representative_payload["below"]["g_records"], dtype=float),
        below_e=np.asarray(representative_payload["below"]["e_records"], dtype=float),
        near_g=np.asarray(representative_payload["near"]["g_records"], dtype=float),
        near_e=np.asarray(representative_payload["near"]["e_records"], dtype=float),
        above_g=np.asarray(representative_payload["above"]["g_records"], dtype=float),
        above_e=np.asarray(representative_payload["above"]["e_records"], dtype=float),
    )

    plot_regime_map(cfg, regime_sweep)
    plot_detuning_scan(detuning_records)
    plot_bandwidth_scan(bandwidth_records)
    plot_mitigation_tradeoff(mitigation_records)
    plot_representative_histograms(representative_payload)
    plot_convergence(convergence_records)

    print("Study complete.")
    print("Representative leakage threshold:", f"{threshold_amplitude_mhz:.2f} MHz at {representative_duration_ns:.0f} ns")
    print("Wall clock:", f"{summary['wall_clock_s']:.1f} s")


if __name__ == "__main__":
    main()
