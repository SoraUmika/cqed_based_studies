"""Bootstrap a cost-weighted candidate frontier from the earlier hybrid study.

This script does not rerun the entire legacy synthesis campaign. Instead, it
collects the validated candidate sequences and replay metrics from the previous
hybrid-qubit-cavity study, merges them into a single bookkeeping layer, and
re-ranks the candidates using an entangling-aware objective.

Outputs
-------
- data/phase1_candidate_bootstrap.json
- figures/phase1_candidate_bootstrap.png
- figures/phase1_candidate_bootstrap.pdf
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Mapping

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import (
    DATA_DIR,
    FIG_DIR,
    LEGACY_STUDY_ROOT,
    NATIVE_WAIT_DURATION_NS,
    CostWeights,
    apply_publication_style,
    as_plain_dict,
    average_gate_fidelity_from_process,
    candidate_weighted_cost,
    count_gate_types,
    dump_json,
    entangling_time_from_sequence,
    implementation_complexity_score,
    infer_family,
    is_entangling_gate,
    load_json,
    local_gate_count_from_sequence,
    normalize_label,
    best_available_fidelity,
)

sys.stdout.reconfigure(encoding="utf-8")

OUTPUT_JSON = DATA_DIR / "phase1_candidate_bootstrap.json"
OUTPUT_PNG = FIG_DIR / "phase1_candidate_bootstrap.png"
OUTPUT_PDF = FIG_DIR / "phase1_candidate_bootstrap.pdf"

IDEAL_FRONTIER_PATH = LEGACY_STUDY_ROOT / "data" / "speed_limit_feasibility" / "ideal_frontier.json"
REFINED_SUMMARY_PATH = LEGACY_STUDY_ROOT / "data" / "speed_limit_feasibility" / "strategy_summary_refined.json"
FOLLOWUP_PATH = LEGACY_STUDY_ROOT / "data" / "followup_optimization" / "followup_results.json"
REPLAY_PATH = LEGACY_STUDY_ROOT / "data" / "extension_pass" / "sequence_replay_results.json"
CPSQR_REPLAY_PATH = LEGACY_STUDY_ROOT / "data" / "iter10_l2d_replay" / "cpsqr_replay_results.json"

FAMILY_COLORS = {
    "native": "#228833",
    "SNAP": "#4477AA",
    "SQR": "#EE6677",
    "CPSQR": "#CCBB44",
    "mixed": "#AA3377",
    "unknown": "#BBBBBB",
}

LEVEL_MARKERS = {
    "pulse": "o",
    "compiled-estimate": "s",
    "ideal-restart": "^",
    "ideal-frontier": "D",
    "missing": "x",
}


def _base_candidate(label: str) -> dict[str, Any]:
    return {
        "label": label,
        "source_labels": [],
        "source_files": [],
        "family": None,
        "target_key": None,
        "sequence": None,
        "gate_types": [],
        "depth": None,
        "local_gate_count": None,
        "entangling_gate_count": None,
        "entangling_time_ns": None,
        "total_duration_ns": None,
        "active_tones": 0,
        "implementation_complexity": None,
        "ideal_fidelity": None,
        "ideal_leakage": None,
        "process_fidelity": None,
        "average_gate_fidelity": None,
        "leakage_average": None,
        "leakage_worst": None,
        "compiled_estimated_fidelity": None,
        "pulse_fidelity": None,
        "pulse_leakage_average": None,
        "pulse_leakage_worst": None,
        "comparison_fidelity": None,
        "comparison_level": None,
        "weighted_cost": None,
        "notes": [],
    }


def _ensure_candidate(store: dict[str, dict[str, Any]], raw_label: str) -> dict[str, Any]:
    label = normalize_label(raw_label)
    candidate = store.setdefault(label, _base_candidate(label))
    if raw_label not in candidate["source_labels"]:
        candidate["source_labels"].append(raw_label)
    return candidate


def _remember_source(candidate: dict[str, Any], source_name: str) -> None:
    if source_name not in candidate["source_files"]:
        candidate["source_files"].append(source_name)


def _set(candidate: dict[str, Any], key: str, value: Any, *, overwrite: bool = False) -> None:
    if value is None:
        return
    if overwrite or candidate.get(key) in (None, [], ""):
        candidate[key] = value


def _sequence_gate_counts(sequence: list[dict[str, Any]]) -> tuple[int, int, float, list[str]]:
    gate_type_counts = count_gate_types(sequence)
    entangling_gate_count = sum(count for gate_type, count in gate_type_counts.items() if is_entangling_gate(gate_type))
    local_gate_count = local_gate_count_from_sequence(sequence)
    entangling_time_ns = entangling_time_from_sequence(sequence)
    return local_gate_count, entangling_gate_count, entangling_time_ns, sorted(gate_type_counts)


def update_from_ideal_frontier(store: dict[str, dict[str, Any]]) -> None:
    entries = load_json(IDEAL_FRONTIER_PATH)
    for entry in entries:
        candidate = _ensure_candidate(store, str(entry.get("label", "unknown")))
        _remember_source(candidate, IDEAL_FRONTIER_PATH.name)
        sequence = list(entry.get("sequence") or [])
        _set(candidate, "family", entry.get("family"), overwrite=True)
        _set(candidate, "target_key", entry.get("target_key"))
        if sequence:
            local_count, ent_count, ent_time_ns, gate_types = _sequence_gate_counts(sequence)
            _set(candidate, "sequence", sequence, overwrite=True)
            _set(candidate, "depth", len(sequence), overwrite=True)
            _set(candidate, "local_gate_count", local_count, overwrite=True)
            _set(candidate, "entangling_gate_count", ent_count, overwrite=True)
            _set(candidate, "entangling_time_ns", ent_time_ns, overwrite=True)
            _set(candidate, "gate_types", gate_types, overwrite=True)
            _set(candidate, "family", infer_family(sequence, entry.get("family")), overwrite=True)
            total_duration_ns = sum(1.0e9 * float(gate.get("duration", 0.0)) for gate in sequence)
            _set(candidate, "total_duration_ns", total_duration_ns, overwrite=True)
        _set(candidate, "ideal_fidelity", entry.get("ideal_fidelity"), overwrite=True)
        _set(candidate, "ideal_leakage", entry.get("ideal_leakage"), overwrite=True)


def update_from_followup_results(store: dict[str, dict[str, Any]]) -> None:
    payload = load_json(FOLLOWUP_PATH)
    structured_restarts = payload.get("structured_restarts", {})
    for raw_label, body in structured_restarts.items():
        candidate = _ensure_candidate(store, raw_label)
        _remember_source(candidate, FOLLOWUP_PATH.name)
        best_restart = None
        restarts = body.get("restarts") or []
        if restarts:
            best_restart = max(restarts, key=lambda row: float(row.get("process_fidelity", -1.0)))
        summary = body.get("summary") or {}
        sequence = list((best_restart or {}).get("sequence") or [])
        if sequence:
            local_count, ent_count, ent_time_ns, gate_types = _sequence_gate_counts(sequence)
            _set(candidate, "sequence", sequence, overwrite=True)
            _set(candidate, "depth", len(sequence), overwrite=True)
            _set(candidate, "local_gate_count", local_count, overwrite=True)
            _set(candidate, "entangling_gate_count", ent_count, overwrite=True)
            _set(candidate, "entangling_time_ns", ent_time_ns, overwrite=True)
            _set(candidate, "gate_types", gate_types, overwrite=True)
            _set(candidate, "family", infer_family(sequence, candidate.get("family")), overwrite=True)
        _set(candidate, "process_fidelity", summary.get("best_process_fidelity"), overwrite=True)
        if summary.get("best_process_fidelity") is not None:
            _set(
                candidate,
                "average_gate_fidelity",
                average_gate_fidelity_from_process(float(summary["best_process_fidelity"]), 4),
                overwrite=True,
            )
        _set(candidate, "leakage_average", summary.get("best_leakage_average"), overwrite=True)
        _set(candidate, "total_duration_ns", summary.get("best_duration_ns"), overwrite=True)
        if best_restart is not None:
            _set(candidate, "leakage_worst", best_restart.get("leakage_worst"), overwrite=True)


def update_from_refined_summary(store: dict[str, dict[str, Any]]) -> None:
    entries = load_json(REFINED_SUMMARY_PATH)
    for entry in entries:
        candidate = _ensure_candidate(store, str(entry.get("label", "unknown")))
        _remember_source(candidate, REFINED_SUMMARY_PATH.name)
        _set(candidate, "family", entry.get("family"), overwrite=True)
        _set(candidate, "target_key", entry.get("target_key"))
        _set(candidate, "compiled_estimated_fidelity", entry.get("estimated_sequence_fidelity"), overwrite=True)
        _set(candidate, "total_duration_ns", entry.get("total_duration_ns"), overwrite=True)
        _set(candidate, "active_tones", int(entry.get("total_active_tones", 0)), overwrite=True)
        gate_details = list(entry.get("gate_details") or [])
        if gate_details:
            entangling_time_ns = sum(float(detail.get("duration_ns", 0.0)) for detail in gate_details)
            _set(candidate, "entangling_gate_count", len(gate_details), overwrite=True)
            _set(candidate, "entangling_time_ns", entangling_time_ns, overwrite=True)
        elif candidate["label"] == "A_ent":
            _set(candidate, "entangling_gate_count", 1, overwrite=True)
            _set(candidate, "entangling_time_ns", float(entry.get("total_duration_ns", 0.0)), overwrite=True)
            _set(candidate, "depth", 1, overwrite=True)
            _set(candidate, "gate_types", ["FreeEvolveCondPhase"], overwrite=True)
            candidate["notes"].append("Native baseline inferred as a single free-evolution entangler.")


def update_from_replay(store: dict[str, dict[str, Any]], replay_path: Path) -> None:
    payload = load_json(replay_path)
    for raw_label, body in payload.items():
        candidate = _ensure_candidate(store, raw_label)
        _remember_source(candidate, replay_path.name)
        _set(candidate, "pulse_fidelity", body.get("pulse_fidelity"), overwrite=True)
        _set(candidate, "pulse_leakage_average", body.get("leakage_average"), overwrite=True)
        _set(candidate, "pulse_leakage_worst", body.get("leakage_worst"), overwrite=True)
        _set(candidate, "ideal_fidelity", body.get("ideal_fidelity"), overwrite=False)
        _set(candidate, "depth", body.get("gate_count"), overwrite=False)
        gate_types = list(body.get("gate_types") or [])
        if gate_types:
            _set(candidate, "gate_types", sorted(set(gate_types)), overwrite=False)
            if candidate.get("family") in (None, "unknown"):
                pseudo_sequence = [{"type": gate_type} for gate_type in gate_types]
                _set(candidate, "family", infer_family(pseudo_sequence, None), overwrite=False)


def finalize_candidates(store: dict[str, dict[str, Any]], weights: CostWeights) -> list[dict[str, Any]]:
    ordered: list[dict[str, Any]] = []
    for candidate in store.values():
        candidate["family"] = candidate.get("family") or infer_family(candidate.get("sequence"), None)
        candidate["active_tones"] = int(candidate.get("active_tones") or 0)
        candidate["gate_types"] = sorted(set(candidate.get("gate_types") or []))
        if candidate.get("sequence") and candidate.get("depth") is None:
            candidate["depth"] = len(candidate["sequence"])
        if candidate.get("sequence") and candidate.get("local_gate_count") is None:
            candidate["local_gate_count"] = local_gate_count_from_sequence(candidate["sequence"])
        if candidate.get("sequence") and candidate.get("entangling_gate_count") is None:
            counts = count_gate_types(candidate["sequence"])
            candidate["entangling_gate_count"] = sum(
                count for gate_type, count in counts.items() if is_entangling_gate(gate_type)
            )
        if candidate.get("sequence") and candidate.get("entangling_time_ns") is None:
            candidate["entangling_time_ns"] = entangling_time_from_sequence(candidate["sequence"])
        if candidate["label"] == "A_ent" and candidate.get("entangling_gate_count") is None:
            candidate["entangling_gate_count"] = 1
            candidate["entangling_time_ns"] = float(candidate.get("total_duration_ns") or NATIVE_WAIT_DURATION_NS)
            candidate["depth"] = int(candidate.get("depth") or 1)
            candidate["family"] = "native"
            candidate["gate_types"] = ["FreeEvolveCondPhase"]
        candidate["implementation_complexity"] = implementation_complexity_score(
            sequence=candidate.get("sequence"),
            active_tones=candidate.get("active_tones") or 0,
            gate_types=candidate.get("gate_types") or [],
        )
        comparison_fidelity, comparison_level = best_available_fidelity(candidate)
        candidate["comparison_fidelity"] = comparison_fidelity
        candidate["comparison_level"] = comparison_level
        candidate["weighted_cost"] = candidate_weighted_cost(candidate, weights)
        candidate["entangling_efficiency"] = comparison_fidelity / max(
            float(candidate.get("entangling_time_ns") or NATIVE_WAIT_DURATION_NS),
            1.0,
        )
        ordered.append(candidate)
    ordered.sort(key=lambda row: (row["weighted_cost"], -row["comparison_fidelity"], row["label"]))
    return ordered


def plot_frontier(candidates: list[dict[str, Any]]) -> None:
    apply_publication_style()
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8))
    panel_specs = (
        ("depth", "Gate depth"),
        ("entangling_time_ns", "Entangling time (ns)"),
        ("weighted_cost", "Weighted physical cost"),
    )
    for ax, (x_key, x_label) in zip(axes, panel_specs):
        for candidate in candidates:
            family = str(candidate.get("family") or "unknown")
            level = str(candidate.get("comparison_level") or "missing")
            color = FAMILY_COLORS.get(family, FAMILY_COLORS["unknown"])
            marker = LEVEL_MARKERS.get(level, "x")
            x_value = float(candidate.get(x_key) or 0.0)
            y_value = float(candidate.get("comparison_fidelity") or 0.0)
            size = 45.0 + 10.0 * float(candidate.get("entangling_gate_count") or 0.0)
            ax.scatter(x_value, y_value, s=size, color=color, marker=marker, edgecolor="black", linewidth=0.5)
            ax.annotate(
                candidate["label"],
                (x_value, y_value),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=8,
            )
        ax.axhline(0.99, linestyle="--", linewidth=1.0, color="0.5")
        ax.axhline(0.90, linestyle=":", linewidth=1.0, color="0.6")
        ax.set_xlabel(x_label)
        ax.set_ylabel("Best available comparison fidelity")
        ax.set_ylim(-0.02, 1.05)
        ax.grid(True, alpha=0.25)
    handles = []
    labels = []
    for family, color in FAMILY_COLORS.items():
        handles.append(plt.Line2D([0], [0], marker="o", linestyle="", color=color, markeredgecolor="black"))
        labels.append(family)
    axes[2].legend(handles, labels, title="Family", loc="lower left", fontsize=8)
    fig.suptitle("Phase 1 bootstrap frontier: fidelity vs depth, entangling time, and weighted cost")
    fig.tight_layout()
    fig.savefig(OUTPUT_PNG, dpi=300, bbox_inches="tight")
    fig.savefig(OUTPUT_PDF, bbox_inches="tight")
    plt.close(fig)


def print_summary(candidates: list[dict[str, Any]]) -> None:
    print("Top cost-weighted candidates")
    print("=" * 96)
    print(
        f"{'Label':<10} {'Family':<8} {'Depth':>5} {'N_ent':>6} {'T_ent(ns)':>10} "
        f"{'F_cmp':>8} {'Level':<16} {'Cost':>8}"
    )
    print("-" * 96)
    for candidate in candidates[:10]:
        print(
            f"{candidate['label']:<10} {str(candidate['family']):<8} "
            f"{int(candidate.get('depth') or 0):>5d} "
            f"{int(candidate.get('entangling_gate_count') or 0):>6d} "
            f"{float(candidate.get('entangling_time_ns') or 0.0):>10.1f} "
            f"{float(candidate.get('comparison_fidelity') or 0.0):>8.4f} "
            f"{str(candidate.get('comparison_level') or 'missing'):<16} "
            f"{float(candidate.get('weighted_cost') or 0.0):>8.4f}"
        )


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    store: dict[str, dict[str, Any]] = {}
    weights = CostWeights()

    update_from_ideal_frontier(store)
    update_from_followup_results(store)
    update_from_refined_summary(store)
    update_from_replay(store, REPLAY_PATH)
    update_from_replay(store, CPSQR_REPLAY_PATH)

    candidates = finalize_candidates(store, weights)
    payload = {
        "metadata": {
            "source_study": str(LEGACY_STUDY_ROOT),
            "native_wait_duration_ns": NATIVE_WAIT_DURATION_NS,
            "cost_weights": as_plain_dict(weights),
            "notes": [
                "Comparison fidelity prefers pulse replay when available, then compiled estimate, then ideal restart, then ideal frontier.",
                "This bootstrap is a candidate-ranking baseline, not the final study conclusion.",
            ],
        },
        "candidates": candidates,
    }
    dump_json(OUTPUT_JSON, payload)
    plot_frontier(candidates)
    print_summary(candidates)
    print(f"\nWrote {OUTPUT_JSON}")
    print(f"Wrote {OUTPUT_PNG}")
    print(f"Wrote {OUTPUT_PDF}")


if __name__ == "__main__":
    main()