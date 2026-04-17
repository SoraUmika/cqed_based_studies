"""Phase 4: depth-resolved diagnostics for native-heavy hybrid-unitary candidates.

This phase uses checkpointed cqed_sim gate-sequence propagation on the ideal
gate models to quantify how representative Bloch observables and cavity Wigner
structure evolve with gate depth for the leading native-entangling candidates.

It does not attempt pulse-bridge replay. The output therefore answers the
depth-diagnostics part of the study while keeping waveform-backed validation as
an explicitly separate limitation.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import qutip as qt

import runtime_compat  # noqa: F401
from common import (
    DATA_DIR,
    FIG_DIR,
    N_TR,
    apply_publication_style,
    build_model,
    dump_json,
    embed_logical_state,
    ensure_sim_root_on_path,
    load_json,
    target_unitary_matrix,
)

sys.stdout.reconfigure(encoding="utf-8")

ensure_sim_root_on_path()

from cqed_sim.observables import bloch_trajectory_from_states, wigner_negativity  # noqa: E402
from cqed_sim.plotting.bloch_plots import GATE_COLORS, plot_bloch_track  # noqa: E402
from cqed_sim.plotting.wigner_grids import plot_wigner_grid  # noqa: E402
from cqed_sim.sim.extractors import cavity_wigner, reduced_cavity_state  # noqa: E402

from phase2_native_block_search import build_candidate_sequence  # noqa: E402

PHASE2_JSON = DATA_DIR / "phase2_native_block_search.json"
OUTPUT_JSON = DATA_DIR / "phase4_depth_diagnostics.json"
SUMMARY_PNG = FIG_DIR / "phase4_probe_fidelity_summary.png"
SUMMARY_PDF = FIG_DIR / "phase4_probe_fidelity_summary.pdf"

SELECTED_LABELS = (
    "N1_exact_hc_to_exact_hc",
    "N2_exact_hc_to_exact_hc",
    "N2_A_local_to_A_local",
)

ROLE_BY_LABEL = {
    "N1_exact_hc_to_exact_hc": "one_wait_baseline",
    "N2_exact_hc_to_exact_hc": "ideal_two_wait_upper_bound",
    "N2_A_local_to_A_local": "best_local_calibrated_candidate",
}

SUMMARY_PROBES: tuple[tuple[str, np.ndarray], ...] = (
    ("g0", np.array([1.0, 0.0, 0.0, 0.0], dtype=np.complex128)),
    ("g1", np.array([0.0, 1.0, 0.0, 0.0], dtype=np.complex128)),
    ("e0", np.array([0.0, 0.0, 1.0, 0.0], dtype=np.complex128)),
    ("e1", np.array([0.0, 0.0, 0.0, 1.0], dtype=np.complex128)),
    ("qx_plus_0", np.array([1.0, 0.0, 1.0, 0.0], dtype=np.complex128) / np.sqrt(2.0)),
    ("g_c_plus", np.array([1.0, 1.0, 0.0, 0.0], dtype=np.complex128) / np.sqrt(2.0)),
)

BLOCH_PROBE_LABEL = "qx_plus_0"
WIGNER_PROBE_LABEL = "g0"
CONDITIONED_LEVELS = (0, 1)

GATE_COLORS.update(
    {
        "PrimitiveGate": "0.25",
        "QubitRotation": "tab:orange",
        "SNAP": "tab:purple",
        "FreeEvolveCondPhase": "tab:red",
        "ConditionalPhaseSQR": "tab:green",
    }
)


def selected_rows() -> list[dict[str, Any]]:
    payload = load_json(PHASE2_JSON)
    by_label = {str(row["label"]): row for row in payload["candidates"]}
    return [dict(by_label[label]) for label in SELECTED_LABELS]


def short_gate_label(gate: Any) -> str:
    gate_type = str(gate.type)
    if gate_type == "PrimitiveGate":
        return str(gate.name)
    aliases = {
        "QubitRotation": "R",
        "Displacement": "D",
        "SNAP": "SNAP",
        "SQR": "SQR",
        "ConditionalPhaseSQR": "CPSQR",
        "FreeEvolveCondPhase": "W",
    }
    return aliases.get(gate_type, gate_type)


def pure_state_fidelity(target: qt.Qobj, actual: qt.Qobj) -> float:
    if not target.isket or not actual.isket:
        raise ValueError("pure_state_fidelity expects ket inputs.")
    overlap = complex(target.overlap(actual))
    return float(abs(overlap) ** 2)


def propagate_with_checkpoints(sequence: Any, initial_state: qt.Qobj, *, model: Any) -> list[qt.Qobj]:
    checkpoints = list(range(len(sequence.gates) + 1))
    history = sequence.propagate_states_with_checkpoints(
        [initial_state],
        checkpoints,
        backend="ideal",
        backend_settings={"model": model},
    )
    return [history[index][0] for index in checkpoints]


def build_track(sequence: Any, states: list[qt.Qobj], *, case: str, target_state: qt.Qobj) -> dict[str, Any]:
    trajectory = bloch_trajectory_from_states(
        states,
        conditioned_n_levels=list(CONDITIONED_LEVELS),
        probability_threshold=1.0e-8,
    )
    snapshots: list[dict[str, Any]] = []
    wigner_snapshots: list[dict[str, Any]] = []
    target_overlap: list[float] = []
    negativity: list[float] = []

    for index, state in enumerate(states):
        if index == 0:
            gate_type = "INIT"
            gate_name = "INIT"
            top_label = "INIT"
        else:
            gate = sequence.gates[index - 1]
            gate_type = str(gate.type)
            gate_name = str(gate.name)
            top_label = short_gate_label(gate)

        target_overlap.append(pure_state_fidelity(target_state, state))
        snapshots.append(
            {
                "index": index,
                "gate_type": gate_type,
                "gate_name": gate_name,
                "top_label": top_label,
            }
        )
        rho_c = reduced_cavity_state(state)
        xvec, yvec, w = cavity_wigner(rho_c, n_points=41, extent=2.0)
        wigner_snapshot = {
            "index": index,
            "gate_type": gate_type,
            "gate_name": gate_name,
            "top_label": top_label,
            "wigner": {
                "xvec": np.asarray(xvec, dtype=float),
                "yvec": np.asarray(yvec, dtype=float),
                "w": np.asarray(w, dtype=float),
            },
        }
        wigner_snapshots.append(wigner_snapshot)
        negativity.append(float(wigner_negativity(wigner_snapshot)))

    return {
        "case": case,
        "indices": np.arange(len(states), dtype=int),
        "x": np.asarray(trajectory["x"], dtype=float),
        "y": np.asarray(trajectory["y"], dtype=float),
        "z": np.asarray(trajectory["z"], dtype=float),
        "conditioned": trajectory["conditioned"],
        "snapshots": snapshots,
        "wigner_snapshots": wigner_snapshots,
        "target_overlap": np.asarray(target_overlap, dtype=float),
        "wigner_negativity": np.asarray(negativity, dtype=float),
    }


def save_bloch_figure(track: dict[str, Any], *, title: str, stem: str) -> tuple[Path, Path]:
    figure = plot_bloch_track(track, title=title, label_stride=1)
    png_path = FIG_DIR / f"{stem}.png"
    pdf_path = FIG_DIR / f"{stem}.pdf"
    figure.savefig(png_path, dpi=300, bbox_inches="tight")
    figure.savefig(pdf_path, bbox_inches="tight")
    plt.close(figure)
    return png_path, pdf_path


def save_wigner_figure(track: dict[str, Any], *, title: str, stem: str) -> tuple[Path, Path]:
    figure = plot_wigner_grid(track, title=title, stride=1, show_colorbar=False)
    if figure is None:
        raise RuntimeError(f"plot_wigner_grid returned None for {title}")
    png_path = FIG_DIR / f"{stem}.png"
    pdf_path = FIG_DIR / f"{stem}.pdf"
    figure.savefig(png_path, dpi=300, bbox_inches="tight")
    figure.savefig(pdf_path, bbox_inches="tight")
    plt.close(figure)
    return png_path, pdf_path


def serialize_conditioned(conditioned: dict[int, dict[str, np.ndarray]]) -> dict[str, dict[str, list[float] | list[bool]]]:
    payload: dict[str, dict[str, list[float] | list[bool]]] = {}
    for level, row in conditioned.items():
        payload[str(level)] = {
            "x": np.asarray(row["x"], dtype=float).tolist(),
            "y": np.asarray(row["y"], dtype=float).tolist(),
            "z": np.asarray(row["z"], dtype=float).tolist(),
            "probability": np.asarray(row["probability"], dtype=float).tolist(),
            "valid": np.asarray(row["valid"], dtype=bool).tolist(),
        }
    return payload


def make_summary_figure(final_fidelities: dict[str, dict[str, float]]) -> None:
    apply_publication_style()
    candidates = list(SELECTED_LABELS)
    probes = [label for label, _ in SUMMARY_PROBES]
    matrix = np.asarray(
        [[float(final_fidelities[candidate][probe]) for probe in probes] for candidate in candidates],
        dtype=float,
    )

    fig, ax = plt.subplots(figsize=(9.6, 3.8))
    image = ax.imshow(matrix, cmap="viridis", vmin=0.0, vmax=1.0, aspect="auto")
    ax.set_xticks(np.arange(len(probes)))
    ax.set_xticklabels(probes, rotation=30, ha="right")
    ax.set_yticks(np.arange(len(candidates)))
    ax.set_yticklabels([ROLE_BY_LABEL[candidate] for candidate in candidates])
    ax.set_title("Final target-state fidelity across logical probes")
    ax.set_xlabel("Probe state")
    ax.set_ylabel("Candidate")
    for row in range(matrix.shape[0]):
        for col in range(matrix.shape[1]):
            value = matrix[row, col]
            ax.text(
                col,
                row,
                f"{value:.3f}",
                ha="center",
                va="center",
                color="white" if value < 0.72 else "black",
                fontsize=8,
            )
    fig.colorbar(image, ax=ax, shrink=0.88, label="State fidelity")
    fig.tight_layout()
    fig.savefig(SUMMARY_PNG, dpi=300, bbox_inches="tight")
    fig.savefig(SUMMARY_PDF, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    apply_publication_style()

    target = target_unitary_matrix()
    model = build_model(n_cav=8, n_tr=N_TR)
    phase2_rows = selected_rows()

    probes = {label: vector for label, vector in SUMMARY_PROBES}
    final_fidelities: dict[str, dict[str, float]] = {}
    candidate_payloads: list[dict[str, Any]] = []

    for row in phase2_rows:
        label = str(row["label"])
        waits = int(row["waits"])
        inner_kind = str(row["inner_local_kind"])
        outer_kind = str(row["outer_local_kind"])
        sequence = build_candidate_sequence(waits=waits, inner_kind=inner_kind, outer_kind=outer_kind)

        probe_fidelities: dict[str, float] = {}
        for probe_label, logical_vector in probes.items():
            initial_state = embed_logical_state(logical_vector, n_cav=8)
            target_state = embed_logical_state(target @ logical_vector, n_cav=8)
            final_state = sequence.propagate_states(
                [initial_state],
                backend="ideal",
                backend_settings={"model": model},
            )[0]
            probe_fidelities[probe_label] = pure_state_fidelity(target_state, final_state)
        final_fidelities[label] = probe_fidelities

        bloch_initial = embed_logical_state(probes[BLOCH_PROBE_LABEL], n_cav=8)
        bloch_target = embed_logical_state(target @ probes[BLOCH_PROBE_LABEL], n_cav=8)
        bloch_states = propagate_with_checkpoints(sequence, bloch_initial, model=model)
        bloch_track = build_track(
            sequence,
            bloch_states,
            case=f"{label}::{BLOCH_PROBE_LABEL}",
            target_state=bloch_target,
        )
        bloch_png, bloch_pdf = save_bloch_figure(
            bloch_track,
            title=f"{label}: Bloch trajectory for {BLOCH_PROBE_LABEL}",
            stem=f"phase4_bloch_{label}_{BLOCH_PROBE_LABEL}",
        )

        wigner_initial = embed_logical_state(probes[WIGNER_PROBE_LABEL], n_cav=8)
        wigner_target = embed_logical_state(target @ probes[WIGNER_PROBE_LABEL], n_cav=8)
        wigner_states = propagate_with_checkpoints(sequence, wigner_initial, model=model)
        wigner_track = build_track(
            sequence,
            wigner_states,
            case=f"{label}::{WIGNER_PROBE_LABEL}",
            target_state=wigner_target,
        )
        wigner_png, wigner_pdf = save_wigner_figure(
            wigner_track,
            title=f"{label}: cavity Wigner snapshots for {WIGNER_PROBE_LABEL}",
            stem=f"phase4_wigner_{label}_{WIGNER_PROBE_LABEL}",
        )

        candidate_payloads.append(
            {
                "label": label,
                "role": ROLE_BY_LABEL[label],
                "waits": waits,
                "inner_local_kind": inner_kind,
                "outer_local_kind": outer_kind,
                "phase2_ideal_fidelity": float(row["ideal_fidelity"]),
                "phase2_weighted_cost": float(row["weighted_cost"]),
                "final_probe_fidelities": probe_fidelities,
                "bloch_probe": {
                    "label": BLOCH_PROBE_LABEL,
                    "indices": bloch_track["indices"].astype(int).tolist(),
                    "x": np.asarray(bloch_track["x"], dtype=float).tolist(),
                    "y": np.asarray(bloch_track["y"], dtype=float).tolist(),
                    "z": np.asarray(bloch_track["z"], dtype=float).tolist(),
                    "target_overlap": np.asarray(bloch_track["target_overlap"], dtype=float).tolist(),
                    "conditioned": serialize_conditioned(bloch_track["conditioned"]),
                    "gate_types": [str(snapshot["gate_type"]) for snapshot in bloch_track["snapshots"]],
                    "gate_names": [str(snapshot["gate_name"]) for snapshot in bloch_track["snapshots"]],
                    "top_labels": [str(snapshot["top_label"]) for snapshot in bloch_track["snapshots"]],
                    "figures": {
                        "png": str(bloch_png),
                        "pdf": str(bloch_pdf),
                    },
                },
                "wigner_probe": {
                    "label": WIGNER_PROBE_LABEL,
                    "indices": wigner_track["indices"].astype(int).tolist(),
                    "target_overlap": np.asarray(wigner_track["target_overlap"], dtype=float).tolist(),
                    "wigner_negativity": np.asarray(wigner_track["wigner_negativity"], dtype=float).tolist(),
                    "gate_types": [str(snapshot["gate_type"]) for snapshot in wigner_track["snapshots"]],
                    "gate_names": [str(snapshot["gate_name"]) for snapshot in wigner_track["snapshots"]],
                    "top_labels": [str(snapshot["top_label"]) for snapshot in wigner_track["snapshots"]],
                    "figures": {
                        "png": str(wigner_png),
                        "pdf": str(wigner_pdf),
                    },
                },
            }
        )

    make_summary_figure(final_fidelities)

    payload = {
        "metadata": {
            "description": "Phase 4 gate-depth diagnostics using checkpointed cqed_sim ideal gate propagation.",
            "note": "These diagnostics use the symbolic cqed_sim gate models, not the waveform bridge. Pulse-backed replay remains a separate validation task.",
            "phase2_source": str(PHASE2_JSON),
            "selected_labels": list(SELECTED_LABELS),
            "bloch_probe": BLOCH_PROBE_LABEL,
            "wigner_probe": WIGNER_PROBE_LABEL,
            "summary_probe_labels": [label for label, _ in SUMMARY_PROBES],
            "conditioned_levels": list(CONDITIONED_LEVELS),
            "summary_figures": {
                "png": str(SUMMARY_PNG),
                "pdf": str(SUMMARY_PDF),
            },
        },
        "candidates": candidate_payloads,
    }
    dump_json(OUTPUT_JSON, payload)

    print("Phase 4 depth diagnostics summary")
    print("=" * 100)
    print(f"{'Label':<24} {'Role':<31} {'Avg probe F':>11} {'Bloch final':>11} {'Wigner final':>12}")
    print("-" * 100)
    for row in candidate_payloads:
        avg_probe = float(np.mean(list(row["final_probe_fidelities"].values())))
        bloch_final = float(row["bloch_probe"]["target_overlap"][-1])
        wigner_final = float(row["wigner_probe"]["target_overlap"][-1])
        print(
            f"{row['label']:<24} {row['role']:<31} {avg_probe:>11.4f} {bloch_final:>11.4f} {wigner_final:>12.4f}"
        )
    print(f"\nWrote {OUTPUT_JSON}")
    print(f"Wrote {SUMMARY_PNG}")
    print(f"Wrote {SUMMARY_PDF}")


if __name__ == "__main__":
    main()