"""Phase 2: native-block composition search for the full U_target.

This phase moves beyond archive-only ranking. It reuses the best validated local
and entangling blocks from the earlier hybrid study and composes full-U_target
candidates that explicitly prioritize native chi-wait entangling evolution.

Search space
------------
- One-native-entangler constructions
- Two-native-entangler constructions derived from the local-equivalence identity
  CNOT_{q->c} = (H_q otimes H_c) CNOT_{c->q} (H_q otimes H_c)

Local cavity blocks considered
------------------------------
- exact_hc: ideal logical cavity Hadamard lower bound
- A_local: optimized D-SNAP-D local block from the earlier study
- B_local: SQR-based local comparison block from the earlier study

Entangler block
---------------
- D_ent: optimized native R-W-R entangler from the earlier study
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import runtime_compat  # noqa: F401
from common import (
    ALPHA,
    ARTIFACT_DIR,
    CHI,
    CHIP,
    DATA_DIR,
    FIG_DIR,
    KERR,
    LEGACY_STUDY_ROOT,
    LOGICAL_BASIS_LABELS,
    NATIVE_WAIT_DURATION_NS,
    CostWeights,
    apply_publication_style,
    average_gate_fidelity_from_process,
    candidate_weighted_cost,
    count_gate_types,
    dump_json,
    ensure_sim_root_on_path,
    entangling_time_from_sequence,
    implementation_complexity_score,
    is_entangling_gate,
    logical_subspace_indices,
    target_unitary_matrix,
)

sys.stdout.reconfigure(encoding="utf-8")

ensure_sim_root_on_path()

from cqed_sim.unitary_synthesis import (  # noqa: E402
    ConditionalPhaseSQR,
    Displacement,
    DriftPhaseModel,
    FreeEvolveCondPhase,
    GateSequence,
    PrimitiveGate,
    QubitRotation,
    SNAP,
    SQR,
    Subspace,
    simulate_sequence,
    subspace_unitary_fidelity,
)

N_CAV = 8
FULL_DIM = 2 * N_CAV
DRIFT = DriftPhaseModel(chi=CHI, chi2=CHIP, kerr=KERR)
TARGET = target_unitary_matrix()
SUBSPACE = Subspace.custom(FULL_DIM, [0, 1, N_CAV, N_CAV + 1], LOGICAL_BASIS_LABELS)

FOLLOWUP_PATH = LEGACY_STUDY_ROOT / "data" / "followup_optimization" / "followup_results.json"
IDEAL_FRONTIER_PATH = LEGACY_STUDY_ROOT / "data" / "speed_limit_feasibility" / "ideal_frontier.json"

OUTPUT_JSON = DATA_DIR / "phase2_native_block_search.json"
OUTPUT_PNG = FIG_DIR / "phase2_native_block_search.png"
OUTPUT_PDF = FIG_DIR / "phase2_native_block_search.pdf"

LOCAL_KIND_COLORS = {
    "exact_hc": "#228833",
    "A_local": "#4477AA",
    "B_local": "#EE6677",
}

WAIT_MARKERS = {
    1: "o",
    2: "s",
}


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def best_restart_sequence(label: str) -> list[dict[str, Any]]:
    payload = load_json(FOLLOWUP_PATH)
    restarts = payload["structured_restarts"][label]["restarts"]
    best = max(restarts, key=lambda row: float(row["process_fidelity"]))
    return list(best["sequence"])


def ideal_frontier_sequence(label: str) -> list[dict[str, Any]]:
    for entry in load_json(IDEAL_FRONTIER_PATH):
        if str(entry.get("label")) == label:
            return list(entry.get("sequence") or [])
    raise KeyError(f"ideal frontier label not found: {label}")


def deserialize_gate(spec: dict[str, Any], suffix: str) -> Any:
    gate_type = str(spec["type"])
    name = f"{spec['name']}_{suffix}"
    duration = float(spec.get("duration", 0.0))
    params = list(spec.get("parameters") or [])

    if gate_type == "Displacement":
        real = float(params[0]) if len(params) >= 1 else 0.0
        imag = float(params[1]) if len(params) >= 2 else 0.0
        return Displacement(name=name, alpha=complex(real, imag), duration=duration, optimize_time=False)
    if gate_type == "QubitRotation":
        theta = float(params[0]) if len(params) >= 1 else 0.0
        phi = float(params[1]) if len(params) >= 2 else 0.0
        return QubitRotation(name=name, theta=theta, phi=phi, duration=duration, optimize_time=False)
    if gate_type == "SNAP":
        return SNAP(name=name, phases=[float(x) for x in params], duration=duration, optimize_time=False)
    if gate_type == "SQR":
        half = len(params) // 2
        theta_n = [float(x) for x in params[:half]]
        phi_n = [float(x) for x in params[half:]]
        return SQR(name=name, theta_n=theta_n, phi_n=phi_n, drift_model=DRIFT, duration=duration, optimize_time=False)
    if gate_type == "ConditionalPhaseSQR":
        return ConditionalPhaseSQR(
            name=name,
            phases_n=[float(x) for x in params],
            drift_model=DRIFT,
            duration=duration,
            optimize_time=False,
        )
    if gate_type == "FreeEvolveCondPhase":
        return FreeEvolveCondPhase(name=name, drift_model=DRIFT, duration=duration, optimize_time=False)
    raise ValueError(f"Unsupported gate type: {gate_type}")


def deserialize_sequence(sequence: list[dict[str, Any]], suffix: str) -> list[Any]:
    return [deserialize_gate(spec, f"{suffix}_{index}") for index, spec in enumerate(sequence)]


def exact_qubit_h(name: str, *, n_cav: int = N_CAV) -> PrimitiveGate:
    hadamard = np.array([[1.0, 1.0], [1.0, -1.0]], dtype=np.complex128) / np.sqrt(2.0)
    matrix = np.kron(hadamard, np.eye(int(n_cav), dtype=np.complex128))
    return PrimitiveGate(
        name=name,
        duration=40e-9,
        matrix=matrix,
        optimize_time=False,
        hilbert_dim=2 * int(n_cav),
    )


def exact_cavity_h(name: str, *, n_cav: int = N_CAV) -> PrimitiveGate:
    n_cav = int(n_cav)
    full_dim = 2 * n_cav
    matrix = np.eye(full_dim, dtype=np.complex128)
    hadamard = np.array([[1.0, 1.0], [1.0, -1.0]], dtype=np.complex128) / np.sqrt(2.0)
    for qubit_index in range(2):
        base = qubit_index * n_cav
        matrix[base : base + 2, base : base + 2] = hadamard
    return PrimitiveGate(name=name, duration=40e-9, matrix=matrix, optimize_time=False, hilbert_dim=full_dim)


def local_hc_block(kind: str, suffix: str, *, n_cav: int = N_CAV) -> list[Any]:
    if kind == "exact_hc":
        return [exact_cavity_h(f"Hc_{suffix}", n_cav=n_cav)]
    if kind == "A_local":
        return deserialize_sequence(best_restart_sequence("A_local"), suffix)
    if kind == "B_local":
        return deserialize_sequence(ideal_frontier_sequence("B_local"), suffix)
    raise ValueError(f"Unknown local block kind: {kind}")


def native_entangler_block(suffix: str) -> list[Any]:
    return deserialize_sequence(best_restart_sequence("D_ent"), f"D_ent_{suffix}")


def build_candidate_sequence(*, waits: int, inner_kind: str, outer_kind: str, n_cav: int = N_CAV) -> GateSequence:
    n_cav = int(n_cav)
    gates: list[Any] = []
    if waits == 1:
        gates.extend([exact_qubit_h("Hq_in_1", n_cav=n_cav)])
        gates.extend(local_hc_block(inner_kind, "inner_1", n_cav=n_cav))
        gates.extend(native_entangler_block("1"))
        gates.extend(local_hc_block(outer_kind, "outer_final", n_cav=n_cav))
    elif waits == 2:
        gates.extend([exact_qubit_h("Hq_in_1", n_cav=n_cav)])
        gates.extend(local_hc_block(inner_kind, "inner_1", n_cav=n_cav))
        gates.extend(native_entangler_block("1"))
        gates.extend([exact_qubit_h("Hq_in_2", n_cav=n_cav)])
        gates.extend(local_hc_block(inner_kind, "inner_2", n_cav=n_cav))
        gates.extend(native_entangler_block("2"))
        gates.extend(local_hc_block(outer_kind, "outer_final", n_cav=n_cav))
    else:
        raise ValueError("waits must be 1 or 2")
    return GateSequence(gates=gates, n_cav=n_cav, full_dim=2 * n_cav)


def sequence_cost_summary(serialized: list[dict[str, Any]]) -> tuple[int, int, float, float, int, list[str]]:
    gate_type_counts = count_gate_types(serialized)
    depth = len(serialized)
    entangling_count = sum(count for gate_type, count in gate_type_counts.items() if is_entangling_gate(gate_type))
    entangling_time_ns = entangling_time_from_sequence(serialized)
    total_duration_ns = sum(1.0e9 * float(spec.get("duration", 0.0)) for spec in serialized)
    active_tones = 0
    for spec in serialized:
        if spec.get("type") == "SQR":
            params = list(spec.get("parameters") or [])
            half = len(params) // 2
            active_tones += sum(1 for value in params[:half] if abs(float(value)) > 1.0e-9)
        if spec.get("type") == "ConditionalPhaseSQR":
            active_tones += sum(1 for value in spec.get("parameters") or [] if abs(float(value)) > 1.0e-9)
    return depth, entangling_count, entangling_time_ns, total_duration_ns, active_tones, sorted(gate_type_counts)


def evaluate_candidate(*, waits: int, inner_kind: str, outer_kind: str, weights: CostWeights) -> dict[str, Any]:
    sequence = build_candidate_sequence(waits=waits, inner_kind=inner_kind, outer_kind=outer_kind, n_cav=N_CAV)
    simulation = simulate_sequence(sequence, subspace=SUBSPACE, backend="ideal")
    fidelity = float(subspace_unitary_fidelity(np.asarray(simulation.subspace_operator), TARGET, gauge="global"))
    serialized = sequence.serialize()
    depth, entangling_count, entangling_time_ns, total_duration_ns, active_tones, gate_types = sequence_cost_summary(serialized)
    candidate = {
        "label": f"N{waits}_{inner_kind}_to_{outer_kind}",
        "family": "native_block",
        "waits": waits,
        "inner_local_kind": inner_kind,
        "outer_local_kind": outer_kind,
        "sequence": serialized,
        "depth": depth,
        "entangling_gate_count": entangling_count,
        "entangling_time_ns": entangling_time_ns,
        "total_duration_ns": total_duration_ns,
        "active_tones": active_tones,
        "gate_types": gate_types,
        "ideal_fidelity": fidelity,
        "compiled_estimated_fidelity": fidelity,
        "average_gate_fidelity": average_gate_fidelity_from_process(fidelity, 4),
        "implementation_complexity": implementation_complexity_score(
            sequence=serialized,
            active_tones=active_tones,
            gate_types=gate_types,
        ),
        "notes": [
            "Built from optimized legacy blocks plus ideal local qubit Hadamard references.",
            "Candidate fidelity is the ideal cqed_sim logical-subspace fidelity of the composed sequence.",
        ],
    }
    candidate["weighted_cost"] = candidate_weighted_cost(candidate, weights)
    return candidate


def make_figure(candidates: list[dict[str, Any]]) -> None:
    apply_publication_style()
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 4.8))

    ax = axes[0]
    for candidate in candidates:
        color = LOCAL_KIND_COLORS[str(candidate["inner_local_kind"])]
        marker = WAIT_MARKERS[int(candidate["waits"])]
        x_value = float(candidate["entangling_time_ns"])
        y_value = float(candidate["ideal_fidelity"])
        ax.scatter(x_value, y_value, color=color, marker=marker, s=75, edgecolor="black", linewidth=0.5)
        ax.annotate(candidate["label"], (x_value, y_value), xytext=(4, 4), textcoords="offset points", fontsize=7)
    ax.axvline(NATIVE_WAIT_DURATION_NS, linestyle=":", color="0.5", linewidth=1.0)
    ax.axhline(0.99, linestyle="--", color="0.6", linewidth=1.0)
    ax.set_xlabel("Total entangling time (ns)")
    ax.set_ylabel("Ideal logical process fidelity")
    ax.set_title("Native-block candidates")
    ax.grid(True, alpha=0.25)

    ax = axes[1]
    ordered = sorted(candidates, key=lambda row: (row["weighted_cost"], -row["ideal_fidelity"]))[:10]
    labels = [row["label"] for row in ordered]
    costs = [float(row["weighted_cost"]) for row in ordered]
    fidelities = [float(row["ideal_fidelity"]) for row in ordered]
    colors = [LOCAL_KIND_COLORS[str(row["inner_local_kind"])] for row in ordered]
    y_positions = np.arange(len(ordered))
    bars = ax.barh(y_positions, costs, color=colors, alpha=0.9, edgecolor="black", linewidth=0.5)
    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Weighted physical cost")
    ax.set_title("Top candidates by cost")
    ax.grid(True, axis="x", alpha=0.25)
    for bar, fidelity in zip(bars, fidelities, strict=True):
        ax.text(bar.get_width() + 0.03, bar.get_y() + bar.get_height() / 2.0, f"F={fidelity:.3f}", va="center", fontsize=8)

    handles = [
        plt.Line2D([0], [0], marker="o", linestyle="", color=LOCAL_KIND_COLORS["exact_hc"], markeredgecolor="black", label="exact_hc"),
        plt.Line2D([0], [0], marker="o", linestyle="", color=LOCAL_KIND_COLORS["A_local"], markeredgecolor="black", label="A_local"),
        plt.Line2D([0], [0], marker="o", linestyle="", color=LOCAL_KIND_COLORS["B_local"], markeredgecolor="black", label="B_local"),
        plt.Line2D([0], [0], marker=WAIT_MARKERS[1], linestyle="", color="black", label="1 native wait"),
        plt.Line2D([0], [0], marker=WAIT_MARKERS[2], linestyle="", color="black", label="2 native waits"),
    ]
    axes[0].legend(handles=handles, loc="lower right", fontsize=8, frameon=False)

    fig.suptitle("Phase 2 native-entangler-biased ideal search")
    fig.tight_layout()
    fig.savefig(OUTPUT_PNG, dpi=300, bbox_inches="tight")
    fig.savefig(OUTPUT_PDF, bbox_inches="tight")
    plt.close(fig)


def print_summary(candidates: list[dict[str, Any]]) -> None:
    print("Phase 2 native-block search summary")
    print("=" * 110)
    print(
        f"{'Label':<24} {'Waits':>5} {'Inner':<10} {'Outer':<10} {'Depth':>5} {'T_ent(ns)':>10} {'F_ideal':>8} {'Cost':>8}"
    )
    print("-" * 110)
    for candidate in sorted(candidates, key=lambda row: (row["weighted_cost"], -row["ideal_fidelity"])):
        print(
            f"{candidate['label']:<24} {int(candidate['waits']):>5d} {candidate['inner_local_kind']:<10} {candidate['outer_local_kind']:<10} "
            f"{int(candidate['depth']):>5d} {float(candidate['entangling_time_ns']):>10.1f} {float(candidate['ideal_fidelity']):>8.4f} {float(candidate['weighted_cost']):>8.4f}"
        )


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    weights = CostWeights()

    candidates: list[dict[str, Any]] = []
    for waits in (1, 2):
        for inner_kind in ("exact_hc", "A_local", "B_local"):
            for outer_kind in ("exact_hc", "A_local", "B_local"):
                candidates.append(evaluate_candidate(waits=waits, inner_kind=inner_kind, outer_kind=outer_kind, weights=weights))

    candidates.sort(key=lambda row: (row["weighted_cost"], -row["ideal_fidelity"], row["label"]))
    payload = {
        "metadata": {
            "description": "Phase 2 native-block composition search for full U_target candidates.",
            "legacy_sources": [str(FOLLOWUP_PATH), str(IDEAL_FRONTIER_PATH)],
            "weights": {
                "infidelity": weights.infidelity,
                "entangling_gate_count": weights.entangling_gate_count,
                "entangling_time": weights.entangling_time,
                "depth": weights.depth,
                "implementation_complexity": weights.implementation_complexity,
                "leakage": weights.leakage,
            },
        },
        "candidates": candidates,
    }
    dump_json(OUTPUT_JSON, payload)
    make_figure(candidates)
    print_summary(candidates)
    print(f"\nWrote {OUTPUT_JSON}")
    print(f"Wrote {OUTPUT_PNG}")
    print(f"Wrote {OUTPUT_PDF}")


if __name__ == "__main__":
    main()