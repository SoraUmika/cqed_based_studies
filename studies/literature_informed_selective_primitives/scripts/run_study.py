"""Run the literature-informed selective-primitives study.

This study compares literature-backed pulse families for:
1. selective qubit rotation (SQR / cphase-SQR target),
2. geometric SNAP,
3. noisy replay of the same optimized primitives.
"""

from __future__ import annotations

import time
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import (
    CHI,
    COLOR_MAP,
    DATA_DIR,
    DT,
    FIGURES_DIR,
    LOGICAL_N,
    SNAP_COLOR_MAP,
    SNAP_PHASE,
    SQR_PHI,
    SQR_THETA,
    TARGET_BRANCH,
    build_frame,
    build_model,
    build_noise_spec,
    build_session,
    build_snap_pulses,
    build_sqr_pulses,
    cavity_ground_indices,
    duration_from_chi_t,
    extract_restricted_operator,
    json_dump,
    logical_indices,
    qobj_probability_in_indices,
    sample_total_waveform,
    snap_probe_state_vectors,
    snap_target_operator,
    sqr_probe_state_vectors,
    sqr_relaxed_target_operator_from_actual,
    sqr_strict_target_operator,
    average_target_state_fidelity,
)


SQR_FAMILY_GRIDS: dict[str, dict[str, list[float]]] = {
    "gaussian": {
        "chi_t": [0.8, 1.0, 1.3, 1.7, 2.2, 3.0],
        "shape": [0.16, 0.20, 0.24],
        "amp_scale": [0.90, 1.00, 1.10],
    },
    "cosine_squared": {
        "chi_t": [0.8, 1.0, 1.3, 1.7, 2.2, 3.0],
        "shape": [0.18],
        "amp_scale": [0.90, 1.00, 1.10],
    },
    "flat_top_gaussian": {
        "chi_t": [0.8, 1.0, 1.3, 1.7, 2.2, 3.0],
        "shape": [0.12, 0.18, 0.24],
        "amp_scale": [0.90, 1.00, 1.10],
    },
}

SNAP_FAMILY_GRIDS: dict[str, dict[str, list[float]]] = {
    "gaussian": {
        "chi_t": [0.8, 1.0, 1.3, 1.7, 2.2],
        "shape": [0.16, 0.20, 0.24],
        "amp_scale": [0.90, 1.00, 1.10],
    },
    "flat_top_gaussian": {
        "chi_t": [0.8, 1.0, 1.3, 1.7, 2.2],
        "shape": [0.12, 0.18, 0.24],
        "amp_scale": [0.90, 1.00, 1.10],
    },
}


def _json_ready(value):
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    return value


def _strict_process_fidelity(target: np.ndarray, actual: np.ndarray) -> float:
    dim = float(target.shape[0])
    return float(np.clip(abs(np.trace(target.conj().T @ actual)) ** 2 / (dim * dim), 0.0, 1.0))


def evaluate_sqr_candidate(
    model,
    frame,
    noise_spec,
    *,
    family: str,
    chi_t: float,
    shape_parameter: float,
    amplitude_scale: float,
    dt: float = DT,
) -> dict[str, object]:
    duration_s = duration_from_chi_t(float(chi_t))
    pulses, drive_ops, total_duration_s = build_sqr_pulses(
        model,
        frame,
        family=family,
        branch=TARGET_BRANCH,
        theta=SQR_THETA,
        phi=SQR_PHI,
        duration_s=duration_s,
        shape_parameter=shape_parameter,
        amplitude_scale=amplitude_scale,
    )

    logical_idx = logical_indices(model, LOGICAL_N)
    closed_session = build_session(
        model,
        frame,
        pulses,
        drive_ops,
        total_duration_s=total_duration_s,
        noise=None,
        dt=dt,
    )

    basis_states = [
        model.basis_state(qubit_level, storage_level)
        for storage_level in range(LOGICAL_N)
        for qubit_level in (0, 1)
    ]
    closed_basis_results = closed_session.run_many(basis_states)
    closed_basis_states = [item.final_state for item in closed_basis_results]
    restricted_operator = extract_restricted_operator(closed_basis_states, model, logical_n=LOGICAL_N)

    strict_target = sqr_strict_target_operator(LOGICAL_N, TARGET_BRANCH, SQR_THETA, SQR_PHI)
    relaxed_target, relaxed_meta = sqr_relaxed_target_operator_from_actual(
        restricted_operator,
        logical_n=LOGICAL_N,
        target_branch=TARGET_BRANCH,
        theta=SQR_THETA,
        phi=SQR_PHI,
    )

    strict_process = _strict_process_fidelity(strict_target, restricted_operator)
    relaxed_process = _strict_process_fidelity(relaxed_target, restricted_operator)
    relaxed_branch_mean = float(relaxed_meta["branch_cphase_mean"])
    basis_logical_pops = [qobj_probability_in_indices(state, logical_idx) for state in closed_basis_states]
    closed_basis_in_logical = float(np.mean(basis_logical_pops))
    closed_basis_leakage = float(1.0 - closed_basis_in_logical)

    noisy_session = build_session(
        model,
        frame,
        pulses,
        drive_ops,
        total_duration_s=total_duration_s,
        noise=noise_spec,
        dt=dt,
    )
    probes = sqr_probe_state_vectors(LOGICAL_N)
    noisy_strict_avg, noisy_strict_details = average_target_state_fidelity(
        noisy_session,
        probes,
        strict_target,
        model=model,
        indices=logical_idx,
    )
    noisy_relaxed_avg, noisy_relaxed_details = average_target_state_fidelity(
        noisy_session,
        probes,
        relaxed_target,
        model=model,
        indices=logical_idx,
    )

    return {
        "family": family,
        "chi_t": float(chi_t),
        "duration_us": float(total_duration_s * 1.0e6),
        "shape_parameter": float(shape_parameter),
        "amplitude_scale": float(amplitude_scale),
        "closed_strict_process_fidelity": strict_process,
        "closed_relaxed_branch_mean": relaxed_branch_mean,
        "closed_relaxed_process_fidelity": relaxed_process,
        "closed_basis_in_logical": closed_basis_in_logical,
        "closed_basis_leakage": closed_basis_leakage,
        "relaxed_target_meta": relaxed_meta,
        "noisy_strict_avg_state_fidelity": noisy_strict_avg,
        "noisy_relaxed_avg_state_fidelity": noisy_relaxed_avg,
        "noisy_strict_details": noisy_strict_details,
        "noisy_relaxed_details": noisy_relaxed_details,
    }


def evaluate_snap_candidate(
    model,
    frame,
    noise_spec,
    *,
    family: str,
    chi_t: float,
    shape_parameter: float,
    amplitude_scale: float,
    dt: float = DT,
) -> dict[str, object]:
    pi_pulse_duration_s = duration_from_chi_t(float(chi_t))
    pulses, drive_ops, total_duration_s = build_snap_pulses(
        model,
        frame,
        family=family,
        branch=TARGET_BRANCH,
        phase_angle=SNAP_PHASE,
        duration_s=pi_pulse_duration_s,
        shape_parameter=shape_parameter,
        amplitude_scale=amplitude_scale,
    )

    cavity_idx = cavity_ground_indices(model, LOGICAL_N)
    closed_session = build_session(
        model,
        frame,
        pulses,
        drive_ops,
        total_duration_s=total_duration_s,
        noise=None,
        dt=dt,
    )
    basis_states = [model.basis_state(0, storage_level) for storage_level in range(LOGICAL_N)]
    closed_basis_results = closed_session.run_many(basis_states)
    closed_basis_states = [item.final_state for item in closed_basis_results]
    restricted_operator = extract_restricted_operator(closed_basis_states, model, logical_n=LOGICAL_N)
    restricted_operator = restricted_operator[0::2, :]
    snap_target = snap_target_operator(LOGICAL_N, TARGET_BRANCH, SNAP_PHASE)
    strict_process = _strict_process_fidelity(snap_target, restricted_operator)

    qubit_ground_population = [qobj_probability_in_indices(state, cavity_idx) for state in closed_basis_states]
    mean_ground_population = float(np.mean(qubit_ground_population))
    mean_leakage = float(1.0 - mean_ground_population)

    noisy_session = build_session(
        model,
        frame,
        pulses,
        drive_ops,
        total_duration_s=total_duration_s,
        noise=noise_spec,
        dt=dt,
    )
    probes = snap_probe_state_vectors(LOGICAL_N)
    noisy_avg, noisy_details = average_target_state_fidelity(
        noisy_session,
        probes,
        snap_target,
        model=model,
        indices=cavity_idx,
    )

    return {
        "family": family,
        "chi_t": float(chi_t),
        "pi_pulse_duration_us": float(pi_pulse_duration_s * 1.0e6),
        "total_duration_us": float(total_duration_s * 1.0e6),
        "shape_parameter": float(shape_parameter),
        "amplitude_scale": float(amplitude_scale),
        "closed_process_fidelity": strict_process,
        "closed_ground_manifold_population": mean_ground_population,
        "closed_leakage": mean_leakage,
        "noisy_avg_state_fidelity": noisy_avg,
        "noisy_details": noisy_details,
    }


def run_sqr_grid(model, frame, noise_spec) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    total = sum(
        len(config["chi_t"]) * len(config["shape"]) * len(config["amp_scale"])
        for config in SQR_FAMILY_GRIDS.values()
    )
    done = 0
    for family, config in SQR_FAMILY_GRIDS.items():
        for chi_t in config["chi_t"]:
            for shape in config["shape"]:
                for amp_scale in config["amp_scale"]:
                    done += 1
                    print(
                        f"[SQR {done:03d}/{total}] family={family:>17s} "
                        f"chiT={chi_t:>4.1f} shape={shape:.2f} amp={amp_scale:.2f}",
                        flush=True,
                    )
                    rows.append(
                        evaluate_sqr_candidate(
                            model,
                            frame,
                            noise_spec,
                            family=family,
                            chi_t=chi_t,
                            shape_parameter=shape,
                            amplitude_scale=amp_scale,
                        )
                    )
    return rows


def run_snap_grid(model, frame, noise_spec) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    total = sum(
        len(config["chi_t"]) * len(config["shape"]) * len(config["amp_scale"])
        for config in SNAP_FAMILY_GRIDS.values()
    )
    done = 0
    for family, config in SNAP_FAMILY_GRIDS.items():
        for chi_t in config["chi_t"]:
            for shape in config["shape"]:
                for amp_scale in config["amp_scale"]:
                    done += 1
                    print(
                        f"[SNAP {done:03d}/{total}] family={family:>17s} "
                        f"chiT={chi_t:>4.1f} shape={shape:.2f} amp={amp_scale:.2f}",
                        flush=True,
                    )
                    rows.append(
                        evaluate_snap_candidate(
                            model,
                            frame,
                            noise_spec,
                            family=family,
                            chi_t=chi_t,
                            shape_parameter=shape,
                            amplitude_scale=amp_scale,
                        )
                    )
    return rows


def _best_row(rows: list[dict[str, object]], key: str) -> dict[str, object]:
    return max(rows, key=lambda row: float(row[key]))


def _best_by_family(rows: list[dict[str, object]], families: list[str], key: str) -> dict[str, dict[str, object]]:
    out: dict[str, dict[str, object]] = {}
    for family in families:
        subset = [row for row in rows if row["family"] == family]
        out[family] = _best_row(subset, key)
    return out


def make_sqr_figure(rows: list[dict[str, object]]) -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "grid.linestyle": ":",
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), sharex=True)
    for family, color in COLOR_MAP.items():
        subset = [row for row in rows if row["family"] == family]
        subset = sorted(subset, key=lambda row: float(row["duration_us"]))
        durations = np.asarray([row["duration_us"] for row in subset], dtype=float)
        closed = np.asarray([row["closed_relaxed_branch_mean"] for row in subset], dtype=float)
        noisy = np.asarray([row["noisy_relaxed_avg_state_fidelity"] for row in subset], dtype=float)
        axes[0].scatter(durations, closed, color=color, s=36, alpha=0.8, label=family.replace("_", " "))
        axes[1].scatter(durations, noisy, color=color, s=36, alpha=0.8, label=family.replace("_", " "))
    axes[0].set_title("SQR Closed-System cphase Fidelity")
    axes[1].set_title("SQR Noisy Average Fidelity")
    for ax in axes:
        ax.set_xlabel("Total duration (us)")
        ax.set_ylabel("Fidelity")
        ax.set_ylim(0.70, 1.01)
    axes[0].legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "sqr_family_scan.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "sqr_family_scan.pdf", bbox_inches="tight")
    plt.close(fig)


def make_snap_figure(rows: list[dict[str, object]]) -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "grid.linestyle": ":",
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.2), sharex=True)
    for family, color in SNAP_COLOR_MAP.items():
        subset = [row for row in rows if row["family"] == family]
        subset = sorted(subset, key=lambda row: float(row["total_duration_us"]))
        durations = np.asarray([row["total_duration_us"] for row in subset], dtype=float)
        closed = np.asarray([row["closed_process_fidelity"] for row in subset], dtype=float)
        noisy = np.asarray([row["noisy_avg_state_fidelity"] for row in subset], dtype=float)
        axes[0].scatter(durations, closed, color=color, s=36, alpha=0.85, label=family.replace("_", " "))
        axes[1].scatter(durations, noisy, color=color, s=36, alpha=0.85, label=family.replace("_", " "))
    axes[0].set_title("SNAP Closed Process Fidelity")
    axes[1].set_title("SNAP Noisy Average Fidelity")
    for ax in axes:
        ax.set_xlabel("Total duration (us)")
        ax.set_ylabel("Fidelity")
        ax.set_ylim(0.70, 1.01)
    axes[0].legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "snap_family_scan.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "snap_family_scan.pdf", bbox_inches="tight")
    plt.close(fig)


def make_waveform_summary_figure(
    model,
    frame,
    best_sqr_noisy: dict[str, object],
    best_snap_noisy: dict[str, object],
) -> None:
    sqr_pulses, _, _ = build_sqr_pulses(
        model,
        frame,
        family=str(best_sqr_noisy["family"]),
        branch=TARGET_BRANCH,
        theta=SQR_THETA,
        phi=SQR_PHI,
        duration_s=duration_from_chi_t(float(best_sqr_noisy["chi_t"])),
        shape_parameter=float(best_sqr_noisy["shape_parameter"]),
        amplitude_scale=float(best_sqr_noisy["amplitude_scale"]),
    )
    snap_pulses, _, _ = build_snap_pulses(
        model,
        frame,
        family=str(best_snap_noisy["family"]),
        branch=TARGET_BRANCH,
        phase_angle=SNAP_PHASE,
        duration_s=duration_from_chi_t(float(best_snap_noisy["chi_t"])),
        shape_parameter=float(best_snap_noisy["shape_parameter"]),
        amplitude_scale=float(best_snap_noisy["amplitude_scale"]),
    )

    t_sqr, waveform_sqr = sample_total_waveform(sqr_pulses)
    t_snap, waveform_snap = sample_total_waveform(snap_pulses)

    fig, axes = plt.subplots(2, 2, figsize=(10.8, 6.8))
    axes[0, 0].plot(t_sqr * 1.0e9, np.real(waveform_sqr), color="#1f77b4", label="I")
    axes[0, 0].plot(t_sqr * 1.0e9, np.imag(waveform_sqr), color="#d62728", label="Q")
    axes[0, 0].set_title("Best Noisy SQR Waveform")
    axes[0, 0].set_xlabel("Time (ns)")
    axes[0, 0].set_ylabel("Drive amplitude (rad/s)")
    axes[0, 0].legend(frameon=False)

    axes[0, 1].plot(t_snap * 1.0e9, np.real(waveform_snap), color="#4c78a8", label="I")
    axes[0, 1].plot(t_snap * 1.0e9, np.imag(waveform_snap), color="#f58518", label="Q")
    axes[0, 1].set_title("Best Noisy SNAP Waveform")
    axes[0, 1].set_xlabel("Time (ns)")
    axes[0, 1].set_ylabel("Drive amplitude (rad/s)")
    axes[0, 1].legend(frameon=False)

    summary_labels = ["SQR closed", "SQR noisy", "SNAP closed", "SNAP noisy"]
    summary_values = [
        float(best_sqr_noisy["closed_relaxed_branch_mean"]),
        float(best_sqr_noisy["noisy_relaxed_avg_state_fidelity"]),
        float(best_snap_noisy["closed_process_fidelity"]),
        float(best_snap_noisy["noisy_avg_state_fidelity"]),
    ]
    colors = ["#1f77b4", "#1f77b4", "#f58518", "#f58518"]
    axes[1, 0].bar(summary_labels, summary_values, color=colors, alpha=0.85)
    axes[1, 0].set_ylim(0.70, 1.01)
    axes[1, 0].set_ylabel("Fidelity")
    axes[1, 0].set_title("Best Practical Operating Points")
    axes[1, 0].tick_params(axis="x", rotation=15)

    axes[1, 1].axis("off")
    lines = [
        "Representative best-noisy points",
        f"SQR: {best_sqr_noisy['family']} at chiT={best_sqr_noisy['chi_t']:.1f}",
        f"  closed cphase = {best_sqr_noisy['closed_relaxed_branch_mean']:.4f}",
        f"  noisy avg = {best_sqr_noisy['noisy_relaxed_avg_state_fidelity']:.4f}",
        f"  total duration = {best_sqr_noisy['duration_us']:.3f} us",
        f"SNAP: {best_snap_noisy['family']} at chiT={best_snap_noisy['chi_t']:.1f}",
        f"  closed process = {best_snap_noisy['closed_process_fidelity']:.4f}",
        f"  noisy avg = {best_snap_noisy['noisy_avg_state_fidelity']:.4f}",
        f"  total duration = {best_snap_noisy['total_duration_us']:.3f} us",
        f"Model chi/2pi = {CHI / (2.0 * np.pi) / 1.0e6:.2f} MHz",
    ]
    axes[1, 1].text(
        0.0,
        0.98,
        "\n".join(lines),
        va="top",
        ha="left",
        family="monospace",
        fontsize=10,
    )

    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "best_waveforms_and_summary.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "best_waveforms_and_summary.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    t_start = time.time()
    model = build_model()
    frame = build_frame(model)
    noise_spec = build_noise_spec()

    sqr_rows = run_sqr_grid(model, frame, noise_spec)
    snap_rows = run_snap_grid(model, frame, noise_spec)

    sqr_best_closed = _best_by_family(
        sqr_rows,
        list(SQR_FAMILY_GRIDS.keys()),
        "closed_relaxed_branch_mean",
    )
    sqr_best_noisy = _best_by_family(
        sqr_rows,
        list(SQR_FAMILY_GRIDS.keys()),
        "noisy_relaxed_avg_state_fidelity",
    )
    snap_best_closed = _best_by_family(
        snap_rows,
        list(SNAP_FAMILY_GRIDS.keys()),
        "closed_process_fidelity",
    )
    snap_best_noisy = _best_by_family(
        snap_rows,
        list(SNAP_FAMILY_GRIDS.keys()),
        "noisy_avg_state_fidelity",
    )

    best_sqr_noisy = _best_row(sqr_rows, "noisy_relaxed_avg_state_fidelity")
    best_snap_noisy = _best_row(snap_rows, "noisy_avg_state_fidelity")

    make_sqr_figure(sqr_rows)
    make_snap_figure(snap_rows)
    make_waveform_summary_figure(model, frame, best_sqr_noisy, best_snap_noisy)

    payload = {
        "metadata": {
            "logical_window_size": LOGICAL_N,
            "target_branch": TARGET_BRANCH,
            "sqr_theta_rad": SQR_THETA,
            "snap_phase_rad": SNAP_PHASE,
            "dt_s": DT,
            "wall_clock_seconds": time.time() - t_start,
        },
        "sqr": {
            "candidates": sqr_rows,
            "best_closed_by_family": sqr_best_closed,
            "best_noisy_by_family": sqr_best_noisy,
            "best_noisy_overall": best_sqr_noisy,
        },
        "snap": {
            "candidates": snap_rows,
            "best_closed_by_family": snap_best_closed,
            "best_noisy_by_family": snap_best_noisy,
            "best_noisy_overall": best_snap_noisy,
        },
    }
    json_dump(DATA_DIR / "study_results.json", _json_ready(payload))
    print(f"Saved {DATA_DIR / 'study_results.json'}")
    print(f"Total runtime: {time.time() - t_start:.1f} s")


if __name__ == "__main__":
    main()
