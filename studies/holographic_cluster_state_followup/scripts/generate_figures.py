from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import common as c


STYLE_PATH = c.STUDY_ROOT.parents[1] / ".github" / "skills" / "publication-figures" / "assets" / "cqed_style.mplstyle"
if STYLE_PATH.exists():
    plt.style.use(str(STYLE_PATH))


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def save_figure(fig: plt.Figure, stem: str) -> None:
    png_path = c.FIG_DIR / f"{stem}.png"
    pdf_path = c.FIG_DIR / f"{stem}.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)


def min_wigner_fidelity(block: dict[str, Any]) -> float:
    return float(min(float(values["cavity_fidelity"]) for values in block.values()))


def plot_grape_duration_frontier(grape: dict[str, Any], probe: dict[str, Any]) -> None:
    durations = []
    best = []
    median = []
    worst = []
    open_process = []
    for key, payload in sorted(grape["durations"].items(), key=lambda item: item[1]["duration_ns"]):
        durations.append(float(payload["duration_ns"]))
        best.append(float(payload["summary"]["best_replay_fidelity"]))
        median.append(float(payload["summary"]["median_replay_fidelity"]))
        worst.append(float(payload["summary"]["worst_replay_fidelity"]))
        open_process.append(float(payload["open_process"]["process_fidelity"]))

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5))

    ax = axes[0]
    ax.plot(durations, best, "o-", label="Best replay", color="#4477AA")
    ax.plot(durations, median, "s--", label="Median replay", color="#228833")
    ax.plot(durations, worst, "d:", label="Worst replay", color="#CC6677")
    ax.fill_between(durations, worst, best, color="#4477AA", alpha=0.10)
    ax.set_xlabel("Duration (ns)")
    ax.set_ylabel("Restricted replay fidelity")
    ax.set_title("GRAPE seed statistics")
    ax.set_ylim(0.45, 1.01)
    ax.legend(loc="lower right", fontsize=8)

    ax = axes[1]
    ax.plot(durations, open_process, "o-", color="#AA3377", label="Open-system process")
    probe_durations = sorted(float(key.replace("ns", "")) for key in probe)
    probe_values = [float(probe[f"{int(d)}ns"]["open_process"]["process_fidelity"]) for d in probe_durations]
    ax.plot(probe_durations, probe_values, "s--", color="#EE7733", label="Open-system process (N_cav=12 probe)")
    ax.set_xlabel("Duration (ns)")
    ax.set_ylabel("Process fidelity")
    ax.set_title("Closed vs open-system frontier")
    ax.set_ylim(0.60, 0.97)
    ax.legend(loc="lower right", fontsize=8)

    fig.suptitle("GRAPE fidelity vs duration", y=1.02)
    save_figure(fig, "grape_duration_frontier")


def plot_decomposition_validation(decomp: dict[str, Any]) -> None:
    keys = list(decomp["candidates"].keys())
    labels = [decomp["candidates"][key]["label"] for key in keys]
    ideal = [float(decomp["candidates"][key]["embedded_evaluations"]["2"]["fidelity"]) for key in keys]
    embedded = [float(decomp["candidates"][key]["embedded_evaluations"]["12"]["fidelity"]) for key in keys]
    pulse = [
        float(decomp["candidates"][key]["pulse_replay"].get("12", {}).get("fidelity", np.nan))
        for key in keys
    ]
    leak = [float(decomp["candidates"][key]["embedded_evaluations"]["12"]["leakage_average"]) for key in keys]

    x = np.arange(len(keys))
    width = 0.25

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))

    ax = axes[0]
    ax.bar(x - width, ideal, width=width, label="Ideal N_cav=2", color="#4477AA")
    ax.bar(x, embedded, width=width, label="Embedded N_cav=12", color="#228833")
    ax.bar(x + width, pulse, width=width, label="Pulse replay N_cav=12", color="#CC6677")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel("Fidelity")
    ax.set_title("Structured-candidate validation")
    ax.set_ylim(0.0, 1.05)
    ax.legend(fontsize=8)

    ax = axes[1]
    ax.bar(x, leak, color="#AA3377")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right")
    ax.set_ylabel("Average leakage at N_cav=12")
    ax.set_title("Embedded leakage")
    ax.set_ylim(0.0, 0.9)

    fig.suptitle("Decomposition study results", y=1.02)
    save_figure(fig, "decomposition_validation")


def plot_truncation_convergence(grape: dict[str, Any], probe: dict[str, Any]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5))
    ax_fid, ax_leak = axes

    ncavs = [8, 10, 12, 15]
    for key in ("300ns", "400ns"):
        payload = grape["durations"][key]
        ax_fid.plot(
            ncavs,
            [float(payload["truncation_replays"][str(n)]["fidelity"]) for n in ncavs],
            "o--",
            label=f"{key} replay from N_cav=8 opt",
        )
        ax_leak.plot(
            ncavs,
            [float(payload["truncation_replays"][str(n)]["leakage_average"]) for n in ncavs],
            "o--",
            label=f"{key} replay from N_cav=8 opt",
        )

    for key in ("300ns", "400ns"):
        payload = probe[key]
        ax_fid.plot([12], [float(payload["best_replay_fidelity"])], "s", ms=9, label=f"{key} direct N_cav=12 opt")
        best_seed = max(payload["seed_rows"], key=lambda row: row["replay_fidelity"])
        ax_leak.plot([12], [float(best_seed["replay_leakage_average"])], "s", ms=9, label=f"{key} direct N_cav=12 opt")

    ax_fid.set_xlabel("Replay cavity truncation N_cav")
    ax_fid.set_ylabel("Replay fidelity")
    ax_fid.set_title("GRAPE truncation convergence")
    ax_fid.set_ylim(0.0, 1.0)
    ax_fid.legend(fontsize=7, loc="lower left")

    ax_leak.set_xlabel("Replay cavity truncation N_cav")
    ax_leak.set_ylabel("Average leakage")
    ax_leak.set_title("GRAPE leakage vs truncation")
    ax_leak.set_ylim(0.0, 1.0)
    ax_leak.legend(fontsize=7, loc="upper left")

    fig.suptitle("GRAPE truncation replay checks", y=1.02)
    save_figure(fig, "grape_truncation_convergence")


def plot_wigner_grid(wigner: dict[str, Any]) -> None:
    candidates = [
        ("target", "Target"),
        ("best_decomposition", "Best decomposition"),
        ("best_sqr_like", "Best SQR-like"),
        ("best_grape", "Best GRAPE"),
    ]
    inputs = list(wigner["target"].keys())
    fig, axes = plt.subplots(len(candidates), len(inputs), figsize=(12, 11))
    for row_idx, (key, title) in enumerate(candidates):
        for col_idx, input_label in enumerate(inputs):
            ax = axes[row_idx, col_idx]
            entry = wigner[key][input_label]
            data = np.asarray(entry["target_wigner"] if key == "target" else entry["candidate_wigner"], dtype=float)
            xvec = np.asarray(entry["xvec"], dtype=float)
            yvec = np.asarray(entry["yvec"], dtype=float)
            mesh = ax.contourf(xvec, yvec, data, levels=41, cmap="RdBu_r")
            if row_idx == 0:
                ax.set_title(input_label)
            if col_idx == 0:
                ax.set_ylabel(title)
            ax.set_xticks([])
            ax.set_yticks([])
    cbar = fig.colorbar(mesh, ax=axes.ravel().tolist(), shrink=0.85)
    cbar.set_label("Wigner value")
    fig.suptitle("Wigner-function comparisons on logical basis inputs", y=0.93)
    save_figure(fig, "wigner_comparison_grid")


def plot_wigner_differences(wigner: dict[str, Any]) -> None:
    candidates = [
        ("best_decomposition", "Best decomposition - target"),
        ("best_sqr_like", "Best SQR-like - target"),
        ("best_grape", "Best GRAPE - target"),
    ]
    inputs = list(wigner["target"].keys())
    fig, axes = plt.subplots(len(candidates), len(inputs), figsize=(12, 8.5))
    for row_idx, (key, title) in enumerate(candidates):
        for col_idx, input_label in enumerate(inputs):
            ax = axes[row_idx, col_idx]
            entry = wigner[key][input_label]
            data = np.asarray(entry["difference_wigner"], dtype=float)
            xvec = np.asarray(entry["xvec"], dtype=float)
            yvec = np.asarray(entry["yvec"], dtype=float)
            mesh = ax.contourf(xvec, yvec, data, levels=41, cmap="PiYG")
            if row_idx == 0:
                ax.set_title(input_label)
            if col_idx == 0:
                ax.set_ylabel(title)
            ax.set_xticks([])
            ax.set_yticks([])
    cbar = fig.colorbar(mesh, ax=axes.ravel().tolist(), shrink=0.85)
    cbar.set_label("Delta Wigner")
    fig.suptitle("Difference Wigner functions", y=0.95)
    save_figure(fig, "wigner_difference_grid")


def plot_active_tones(decomp: dict[str, Any]) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    for key, payload in decomp["candidates"].items():
        if payload["category"] not in {"decomposition", "sqr_like", "entangler_assisted"}:
            continue
        tone_count = int(payload["gate_summary"]["max_active_tones"])
        if payload["pulse_replay"].get("12"):
            fidelity = float(payload["pulse_replay"]["12"]["fidelity"])
        else:
            fidelity = float(payload["embedded_evaluations"]["12"]["fidelity"])
        category = payload["category"]
        marker = {"decomposition": "o", "sqr_like": "s", "entangler_assisted": "D"}[category]
        ax.scatter(tone_count, fidelity, s=90, marker=marker, label=payload["label"])
        ax.annotate(payload["label"], (tone_count, fidelity), textcoords="offset points", xytext=(6, 6), fontsize=8)
    ax.set_xlabel("Max active tones per selective gate")
    ax.set_ylabel("Physically relevant fidelity")
    ax.set_title("Active-tone count vs fidelity")
    ax.set_xlim(-0.2, 2.4)
    ax.set_ylim(0.15, 0.7)
    save_figure(fig, "active_tones_vs_fidelity")


def plot_waveforms_and_timeline(probe: dict[str, Any], decomp: dict[str, Any]) -> None:
    best_grape_key = max(probe, key=lambda key: float(probe[key]["best_replay_fidelity"]))
    best_grape_seed = int(probe[best_grape_key]["best_seed"])
    grape_payload = load_json(c.ARTIFACT_DIR / f"grape_nc12_{best_grape_key.replace('ns', '')}ns_seed{best_grape_seed}_probe.json")
    result_payload = grape_payload["result_payload"]
    time_ns = 1.0e9 * np.asarray(result_payload["time_grid_s"], dtype=float)
    command_values = np.asarray(result_payload["command_values"], dtype=float)
    labels = result_payload["control_terms"]

    best_sqr_key = decomp["summary"]["best_sqr_key"]
    sqr_payload = decomp["candidates"][best_sqr_key]

    fig, axes = plt.subplots(2, 1, figsize=(11.0, 6.5), gridspec_kw={"height_ratios": [2.2, 1.3]})

    ax = axes[0]
    for idx, label in enumerate(labels):
        ax.plot(time_ns, command_values[idx] / (2.0 * np.pi * 1.0e6), label=label)
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Command amplitude / 2pi (MHz)")
    ax.set_title("Best N_cav=12 GRAPE waveform schedule")
    ax.legend(ncol=2, fontsize=8)

    ax = axes[1]
    seq = sqr_payload["sequence"]
    cursor = 0.0
    for row in seq:
        duration = float(row["duration_ns"])
        gate_type = row["type"]
        color = {
            "Displacement": "#4477AA",
            "QubitRotation": "#228833",
            "SQR": "#CC6677",
            "ConditionalPhaseSQR": "#AA3377",
        }.get(gate_type, "#999999")
        ax.barh([0], [duration], left=[cursor], color=color, edgecolor="black", height=0.5)
        ax.text(cursor + 0.5 * duration, 0, row["name"], ha="center", va="center", fontsize=7, color="white")
        cursor += duration
    ax.set_xlim(0, max(cursor, 1.0))
    ax.set_yticks([])
    ax.set_xlabel("Time (ns)")
    ax.set_title(f"Best SQR-like sequence timeline ({sqr_payload['label']})")

    fig.suptitle("Waveform and sequence views", y=1.01)
    save_figure(fig, "waveforms_and_timeline")


def plot_final_summary(summary: dict[str, Any], decomp: dict[str, Any], wigner: dict[str, Any]) -> None:
    rows = []
    rows.append(
        (
            "Best decomposition",
            float(decomp["candidates"][summary["best_decomposition_key"]]["embedded_evaluations"]["12"]["fidelity"]),
            min_wigner_fidelity(wigner["best_decomposition"]),
            np.nan,
            float(decomp["candidates"][summary["best_decomposition_key"]]["gate_summary"]["total_duration_ns"]),
        )
    )
    rows.append(
        (
            "Best SQR-like",
            float(summary["best_sqr_nc12_metric"]["fidelity"]),
            min_wigner_fidelity(wigner["best_sqr_like"]),
            np.nan,
            float(decomp["candidates"][summary["best_sqr_key"]]["gate_summary"]["total_duration_ns"]),
        )
    )
    rows.append(
        (
            "Best GRAPE",
            float(summary["best_grape_replay_fidelity"]),
            min_wigner_fidelity(wigner["best_grape"]),
            float(summary["best_grape_open_process_fidelity"]),
            float(summary["best_grape_key"].replace("ns", "")),
        )
    )

    labels = [row[0] for row in rows]
    replay = [row[1] for row in rows]
    wigner_min = [row[2] for row in rows]
    open_process = [row[3] for row in rows]
    durations = [row[4] for row in rows]
    x = np.arange(len(rows))

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))

    ax = axes[0]
    ax.bar(x - 0.18, replay, width=0.18, label="Replay / embedded fidelity", color="#4477AA")
    ax.bar(x, wigner_min, width=0.18, label="Min Wigner cavity fidelity", color="#228833")
    ax.bar(x + 0.18, [0.0 if np.isnan(v) else v for v in open_process], width=0.18, label="Open-system process", color="#CC6677")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel("Metric value")
    ax.set_title("Shortlisted candidate comparison")
    ax.legend(fontsize=8)

    ax = axes[1]
    ax.bar(np.arange(len(labels)), durations, color=["#999999", "#CC79A7", "#0072B2"])
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel("Active time (ns)")
    ax.set_title("Duration comparison")

    fig.suptitle("Final recommendation summary", y=1.02)
    save_figure(fig, "final_summary_comparison")


def main() -> None:
    decomp = load_json(c.DATA_DIR / "decomposition_results.json")
    grape = load_json(c.DATA_DIR / "grape_results.json")
    probe = load_json(c.DATA_DIR / "grape_large_truncation_probe.json")
    wigner = load_json(c.DATA_DIR / "wigner_results.json")
    summary = load_json(c.DATA_DIR / "followup_summary.json")

    plot_grape_duration_frontier(grape, probe)
    plot_decomposition_validation(decomp)
    plot_truncation_convergence(grape, probe)
    plot_wigner_grid(wigner)
    plot_wigner_differences(wigner)
    plot_active_tones(decomp)
    plot_waveforms_and_timeline(probe, decomp)
    plot_final_summary(summary, decomp, wigner)


if __name__ == "__main__":
    main()
