"""Run the hybrid 2x2 qubit+cavity gate-set comparison study.

This script executes the first-pass benchmark defined in the study README:

1. decomposition-based synthesis for gate sets A-D and F,
2. waveform-level GRAPE references for gate set E,
3. lightweight robustness replay for decomposition libraries,
4. publication-quality summary figures and machine-readable outputs.
"""

from __future__ import annotations

import csv
import json
import math
import time
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from common import (
    DATA_DIR,
    FIGURES_DIR,
    GrapeCase,
    ensure_dirs,
    grape_cases,
    json_dump,
    replay_sequence_record,
    run_grape_case,
    run_sequence_case,
    sequence_cases,
)


STUDY_ROOT = Path(__file__).resolve().parents[2]

SCENARIOS: dict[str, dict[str, float]] = {
    "amp_minus_2pct": {"amplitude_scale": 0.98},
    "amp_plus_2pct": {"amplitude_scale": 1.02},
    "chi_minus_1pct": {"chi_scale": 0.99},
    "chi_plus_1pct": {"chi_scale": 1.01},
    "dur_plus_2pct": {"duration_scale": 1.02},
}

ROW_COLUMNS = [
    "label",
    "category",
    "target_key",
    "target_name",
    "description",
    "approach_group",
    "strict_fidelity",
    "block_fidelity",
    "leakage_average",
    "leakage_worst",
    "duration_ns",
    "gate_count",
    "score_strict",
    "score_block",
    "decoherence_proxy",
    "objective",
    "success",
]

LIBRARY_COLORS = {
    "A": "#1b9e77",
    "B": "#d95f02",
    "C": "#7570b3",
    "D": "#e7298a",
    "E": "#66a61e",
    "F": "#e6ab02",
}
CATEGORY_MARKERS = {"decomposition": "o", "grape": "s"}


def _json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, complex):
        return {"real": float(value.real), "imag": float(value.imag)}
    if isinstance(value, (float, int, str, bool)) or value is None:
        return value
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return repr(value)


def _row_to_csvable(row: dict[str, Any]) -> dict[str, Any]:
    csv_row: dict[str, Any] = {}
    for key in ROW_COLUMNS:
        value = row.get(key)
        if isinstance(value, (dict, list, tuple)):
            csv_row[key] = json.dumps(_json_ready(value), sort_keys=True)
        else:
            csv_row[key] = value
    return csv_row


def _library_key(label: str) -> str:
    return label.split("_", 1)[0]


def _approach_group(label: str) -> str:
    if label in {"A_local", "B_local", "B_ent"}:
        return "selective"
    if label in {"A_ent", "C_local", "C_ent", "D_local", "D_ent", "F_local", "F_ent"}:
        return "native"
    return "grape"


def _annotate_row(row: dict[str, Any]) -> dict[str, Any]:
    enriched = dict(row)
    enriched["approach_group"] = _approach_group(str(row["label"]))
    return enriched


def evaluate_robustness(record: dict[str, Any]) -> dict[str, Any]:
    if record["category"] != "decomposition":
        return {
            "label": record["label"],
            "nominal_strict_fidelity": float(record["strict_fidelity"]),
            "nominal_block_fidelity": float(record["block_fidelity"]),
            "min_strict_fidelity": float(record["strict_fidelity"]),
            "min_block_fidelity": float(record["block_fidelity"]),
            "max_leakage_average": float(record["leakage_average"]),
            "strict_drop": 0.0,
            "block_drop": 0.0,
            "scenarios": {},
        }

    scenario_rows: dict[str, dict[str, float]] = {}
    for scenario_name, kwargs in SCENARIOS.items():
        scenario_rows[scenario_name] = replay_sequence_record(record, **kwargs)

    strict_values = [float(row["strict_fidelity"]) for row in scenario_rows.values()]
    block_values = [float(row["block_fidelity"]) for row in scenario_rows.values()]
    leak_values = [float(row["leakage_average"]) for row in scenario_rows.values()]
    min_strict = min(strict_values)
    min_block = min(block_values)
    return {
        "label": record["label"],
        "nominal_strict_fidelity": float(record["strict_fidelity"]),
        "nominal_block_fidelity": float(record["block_fidelity"]),
        "min_strict_fidelity": float(min_strict),
        "min_block_fidelity": float(min_block),
        "max_leakage_average": float(max(leak_values)),
        "strict_drop": float(record["strict_fidelity"] - min_strict),
        "block_drop": float(record["block_fidelity"] - min_block),
        "scenarios": scenario_rows,
    }


def _sort_records(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return sorted(
        rows,
        key=lambda row: (
            row["target_key"],
            -float(row["score_strict"]),
            float(row["duration_ns"]),
            row["label"],
        ),
    )


def _winner_map(rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    winners: dict[str, dict[str, Any]] = {}
    for target_key in sorted({row["target_key"] for row in rows}):
        subset = [row for row in rows if row["target_key"] == target_key]
        strict_best = max(subset, key=lambda row: float(row["score_strict"]))
        block_best = max(subset, key=lambda row: float(row["score_block"]))
        shortest_high_fid = min(
            subset,
            key=lambda row: (
                0 if float(row["block_fidelity"]) >= 0.99 else 1,
                float(row["duration_ns"]),
                -float(row["block_fidelity"]),
            ),
        )
        winners[target_key] = {
            "best_score_strict": strict_best["label"],
            "best_score_block": block_best["label"],
            "fastest_block_fid_ge_0p99": shortest_high_fid["label"],
        }
    return winners


def _set_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "grid.linestyle": ":",
            "figure.facecolor": "white",
        }
    )


def _save_figure(fig: plt.Figure, stem: str) -> None:
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def make_candidate_tradeoff_figure(rows: list[dict[str, Any]], *, target_key: str, stem: str, title: str) -> None:
    subset = [row for row in rows if row["target_key"] == target_key]
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    for row in subset:
        library = _library_key(row["label"])
        ax.scatter(
            row["duration_ns"],
            row["strict_fidelity"],
            s=200.0 * max(0.15, 1.0 - float(row["leakage_average"])),
            marker=CATEGORY_MARKERS[row["category"]],
            color=LIBRARY_COLORS[library],
            edgecolor="black",
            linewidth=0.7,
            alpha=0.9,
        )
        if float(row["block_fidelity"]) > float(row["strict_fidelity"]) + 5.0e-4:
            ax.plot(
                [row["duration_ns"], row["duration_ns"]],
                [row["strict_fidelity"], row["block_fidelity"]],
                color=LIBRARY_COLORS[library],
                linewidth=1.1,
                alpha=0.8,
            )
            ax.scatter(
                row["duration_ns"],
                row["block_fidelity"],
                s=55,
                marker="_",
                color="black",
                linewidth=1.2,
                zorder=3,
            )
        ax.annotate(
            row["label"],
            (row["duration_ns"], row["strict_fidelity"]),
            textcoords="offset points",
            xytext=(6, 5),
            fontsize=8,
        )
    ax.set_xlabel("Duration (ns)")
    ax.set_ylabel("Logical fidelity")
    ax.set_ylim(0.40, 1.02)
    ax.set_title(title)
    legend_items = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#444444", markeredgecolor="black", label="Decomposition", markersize=8),
        plt.Line2D([0], [0], marker="s", color="w", markerfacecolor="#444444", markeredgecolor="black", label="GRAPE", markersize=8),
        plt.Line2D([0], [0], marker="_", color="black", label="Block-gauge fidelity", markersize=12, linewidth=0),
    ]
    ax.legend(handles=legend_items, loc="lower right", frameon=False)
    _save_figure(fig, stem)


def make_robustness_figure(rows: list[dict[str, Any]], robustness: dict[str, dict[str, Any]]) -> None:
    decomposition_rows = [row for row in rows if row["category"] == "decomposition"]
    labels = [row["label"] for row in decomposition_rows]
    nominal = [float(robustness[label]["nominal_strict_fidelity"]) for label in labels]
    worst = [float(robustness[label]["min_strict_fidelity"]) for label in labels]
    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    x = np.arange(len(labels), dtype=float)
    width = 0.38
    ax.bar(x - width / 2.0, nominal, width=width, color="#7fb3d5", edgecolor="black", linewidth=0.7, label="Nominal")
    ax.bar(x + width / 2.0, worst, width=width, color="#f5b041", edgecolor="black", linewidth=0.7, label="Worst perturbed")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel("Strict logical fidelity")
    ax.set_title("Robustness under amplitude, drift, and duration perturbations")
    ax.legend(frameon=False, loc="lower right")
    _save_figure(fig, "fig3_robustness")


def make_score_figure(rows: list[dict[str, Any]]) -> None:
    ordered = _sort_records(rows)
    fig, ax = plt.subplots(figsize=(9.0, 4.8))
    x = np.arange(len(ordered), dtype=float)
    colors = [LIBRARY_COLORS[_library_key(row["label"])] for row in ordered]
    ax.bar(x, [float(row["score_strict"]) for row in ordered], color=colors, edgecolor="black", linewidth=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels([row["label"] for row in ordered], rotation=35, ha="right")
    ax.set_ylabel(r"Composite score $J$")
    ax.set_title("Leakage- and duration-penalized score across all candidates")
    _save_figure(fig, "fig4_scores")


def make_grouped_route_figure(rows: list[dict[str, Any]]) -> None:
    groups = ("selective", "native", "grape")
    group_colors = {"selective": "#4c78a8", "native": "#f58518", "grape": "#54a24b"}
    target_keys = ("local_h", "cx_c_to_q")
    target_labels = {"local_h": r"$I_q \otimes H_c$", "cx_c_to_q": r"$\mathrm{CX}(c\rightarrow q)$"}

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.4), sharey=True)
    for axis, target_key in zip(axes, target_keys, strict=True):
        heights: list[float] = []
        annotations: list[str] = []
        for group in groups:
            subset = [row for row in rows if row["target_key"] == target_key and row["approach_group"] == group]
            best = max(subset, key=lambda row: float(row["score_strict"]))
            heights.append(float(best["score_strict"]))
            annotations.append(str(best["label"]))
        x = np.arange(len(groups), dtype=float)
        axis.bar(x, heights, color=[group_colors[group] for group in groups], edgecolor="black", linewidth=0.7)
        axis.set_xticks(x)
        axis.set_xticklabels([group.title() for group in groups])
        axis.set_title(target_labels[target_key])
        for xpos, height, label in zip(x, heights, annotations, strict=True):
            axis.text(xpos, height + 0.015, label, ha="center", va="bottom", fontsize=8)

    axes[0].set_ylabel(r"Best route score $J$")
    axes[0].set_ylim(0.35, 1.05)
    fig.suptitle("Best selective, native-interaction, and GRAPE routes by target")
    _save_figure(fig, "fig5_route_groups")


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=ROW_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow(_row_to_csvable(row))


def main() -> None:
    ensure_dirs()
    _set_plot_style()
    start = time.perf_counter()

    results: list[dict[str, Any]] = []
    for label, case in sequence_cases().items():
        print(f"[sequence] {label}")
        results.append(_annotate_row(run_sequence_case(case)))
    for label, case in grape_cases().items():
        print(f"[grape] {label}")
        results.append(_annotate_row(run_grape_case(case)))

    robustness = {row["label"]: evaluate_robustness(row) for row in results}
    sorted_results = _sort_records(results)
    winners = _winner_map(sorted_results)

    write_csv(DATA_DIR / "candidate_results.csv", sorted_results)
    json_dump(DATA_DIR / "robustness_summary.json", _json_ready(robustness))

    make_candidate_tradeoff_figure(
        sorted_results,
        target_key="local_h",
        stem="fig1_local_tradeoffs",
        title="Hybrid 2x2 local cavity-action benchmark",
    )
    make_candidate_tradeoff_figure(
        sorted_results,
        target_key="cx_c_to_q",
        stem="fig2_entangler_tradeoffs",
        title="Hybrid entangler benchmark on the logical 2x2 block",
    )
    make_robustness_figure(sorted_results, robustness)
    make_score_figure(sorted_results)
    make_grouped_route_figure(sorted_results)

    elapsed_s = time.perf_counter() - start
    summary = {
        "study_name": "hybrid_universal_control_gate_set_comparison",
        "results": sorted_results,
        "robustness": robustness,
        "winners": winners,
        "scenarios": SCENARIOS,
        "runtime_s": elapsed_s,
        "generated_files": {
            "csv": str(DATA_DIR / "candidate_results.csv"),
            "robustness_json": str(DATA_DIR / "robustness_summary.json"),
            "figures": [
                str(FIGURES_DIR / "fig1_local_tradeoffs.pdf"),
                str(FIGURES_DIR / "fig2_entangler_tradeoffs.pdf"),
                str(FIGURES_DIR / "fig3_robustness.pdf"),
                str(FIGURES_DIR / "fig4_scores.pdf"),
                str(FIGURES_DIR / "fig5_route_groups.pdf"),
            ],
        },
    }
    json_dump(DATA_DIR / "study_summary.json", _json_ready(summary))
    print(f"Completed study run in {elapsed_s:.1f} s")


if __name__ == "__main__":
    main()
