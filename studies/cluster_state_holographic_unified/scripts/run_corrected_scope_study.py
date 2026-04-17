from __future__ import annotations

import argparse
import copy
import math
import statistics
import sys
import time
from pathlib import Path
from typing import Any, Iterable

import matplotlib.pyplot as plt
import numpy as np
import qutip as qt


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import common as c


STYLE_PATH = SCRIPT_DIR.parents[2] / ".github" / "skills" / "publication-figures" / "assets" / "cqed_style.mplstyle"
if STYLE_PATH.exists():
    plt.style.use(str(STYLE_PATH))


SCREEN_MAXITER = 20
REFINE_MAXITER = 60
FINAL_MAXITER = 100
ACTIVE_SAMPLE_POINTS = 12
PRELIM_TRUNC = c.DECOMP_N_CAV
FINAL_TRUNCATIONS = (10, 12, 14)
DEFAULT_FINAL_TRUNC = 12
ACTIVE_THRESHOLD = 1.0e-3
ACTIVE_CAPTURE_TARGET = 0.999
RETENTION_FIDELITY_MIN = 0.90
RETENTION_LEAKAGE_MAX = 0.05
LOW_CONFIDENCE_LEAKAGE_MAX = 0.10
RETENTION_CONVERGENCE_DELTA = 0.01
LOW_CONFIDENCE_CONVERGENCE_DELTA = 0.05

BUILDERS = {
    "drsqr_3sqr": c.build_drsqr_3sqr_sequence,
    "drsqr_4sqr": c.build_drsqr_4sqr_sequence,
    "drcpsqr": c.build_drcpsqr_sequence,
}


def case_id(case: dict[str, Any]) -> str:
    parts = [str(case["family_key"]), str(case["variant_key"])]
    if case.get("pattern"):
        parts.append(str(case["pattern"]))
    if case.get("levels"):
        parts.append("lv" + "-".join(str(level) for level in case["levels"]))
    if case.get("blocks") is not None:
        parts.append(f"blk{int(case['blocks'])}")
    return "_".join(parts)


def ranking_key(record: dict[str, Any]) -> tuple[float, float, float, float, float, float]:
    summary = record["summary"]
    expensive_blocks = int(summary["selective_block_count"]) + int(summary["entangling_block_count"])
    return (
        float(record["fidelity"]),
        -float(summary["total_duration_ns"]),
        -float(summary["gate_depth"]),
        -float(expensive_blocks),
        -float(summary["total_active_tones"]),
        -float(summary["max_active_tones"]),
    )


def best_record(records: Iterable[dict[str, Any]]) -> dict[str, Any]:
    return max(records, key=ranking_key)


def top_records(records: list[dict[str, Any]], count: int) -> list[dict[str, Any]]:
    return sorted(records, key=ranking_key, reverse=True)[: int(count)]


def record_without_fit(record: dict[str, Any]) -> dict[str, Any]:
    return c.json_ready({key: value for key, value in record.items() if key != "fit"})


def build_sequence_from_case(case: dict[str, Any], *, n_cav: int) -> c.GateSequence:
    builder = BUILDERS[str(case["builder_name"])]
    kwargs = dict(case.get("builder_kwargs", {}))
    kwargs["n_cav"] = int(n_cav)
    if case.get("levels") is not None:
        kwargs["levels"] = tuple(int(level) for level in case["levels"])
    if case.get("pattern") is not None:
        kwargs["pattern"] = str(case["pattern"])
    if case.get("blocks") is not None and str(case["builder_name"]) == "drcpsqr":
        kwargs["blocks"] = int(case["blocks"])
    return builder(**kwargs)


def apply_solution_to_case(record: dict[str, Any], *, n_cav: int) -> c.GateSequence:
    sequence = build_sequence_from_case(record, n_cav=int(n_cav))
    sequence.set_parameter_vector(np.asarray(record["parameter_vector"], dtype=float))
    sequence.set_time_vector(np.asarray(record["time_vector"], dtype=float), active_only=False)
    return sequence


def run_synthesis_trial(
    case: dict[str, Any],
    *,
    search_phase: str,
    seed: int,
    init_guess: str,
    maxiter: int,
    multistart: int = 1,
    warm_start: Any | None = None,
) -> dict[str, Any]:
    started = time.perf_counter()
    sequence = build_sequence_from_case(case, n_cav=PRELIM_TRUNC)
    fit = c.fit_sequence(
        sequence,
        n_cav=PRELIM_TRUNC,
        seed=int(seed),
        init_guess=str(init_guess),
        multistart=int(multistart),
        maxiter=int(maxiter),
        duration_weight=float(case.get("duration_weight", 0.0)),
        gate_count_weight=float(case.get("gate_count_weight", 0.0)),
        warm_start=warm_start,
    )
    summary = dict(fit["summary"])
    record = {
        "case_id": case_id(case),
        "family_key": str(case["family_key"]),
        "family_label": str(case["family_label"]),
        "variant_key": str(case["variant_key"]),
        "variant_label": str(case["variant_label"]),
        "builder_name": str(case["builder_name"]),
        "builder_kwargs": dict(case.get("builder_kwargs", {})),
        "pattern": case.get("pattern"),
        "levels": None if case.get("levels") is None else [int(level) for level in case["levels"]],
        "max_tones": None if case.get("max_tones") is None else int(case["max_tones"]),
        "blocks": None if case.get("blocks") is None else int(case["blocks"]),
        "search_phase": str(search_phase),
        "seed": int(seed),
        "init_guess": str(init_guess),
        "maxiter": int(maxiter),
        "multistart": int(multistart),
        "fidelity": float(fit["fidelity"]),
        "objective": float(fit["objective"]),
        "success": bool(fit["success"]),
        "message": str(fit["message"]),
        "summary": summary,
        "metrics": dict(fit["metrics"]),
        "sequence": fit["sequence_payload"],
        "parameter_vector": fit["result"].sequence.get_parameter_vector().tolist(),
        "time_vector": fit["result"].sequence.get_time_vector(active_only=False).tolist(),
        "fit": fit,
        "elapsed_s": float(time.perf_counter() - started),
    }
    return record


def screen_and_refine_cases(
    *,
    cases: list[dict[str, Any]],
    shortlist_per_group: int,
    group_key: str,
    include_random_refine: bool = True,
) -> dict[str, Any]:
    screened: list[dict[str, Any]] = []
    for case in cases:
        print(f"[screen] {case_id(case)}", flush=True)
        screened.append(
            run_synthesis_trial(
                case,
                search_phase="screen",
                seed=17,
                init_guess="heuristic",
                maxiter=SCREEN_MAXITER,
                multistart=1,
            )
        )

    grouped: dict[str, list[dict[str, Any]]] = {}
    for record in screened:
        grouped.setdefault(str(record[group_key]), []).append(record)

    shortlist: dict[str, dict[str, Any]] = {}
    for records in grouped.values():
        for record in top_records(records, shortlist_per_group):
            shortlist[str(record["case_id"])] = {
                "family_key": str(record["family_key"]),
                "family_label": str(record["family_label"]),
                "variant_key": str(record["variant_key"]),
                "variant_label": str(record["variant_label"]),
                "builder_name": str(record["builder_name"]),
                "builder_kwargs": dict(record["builder_kwargs"]),
                "pattern": record.get("pattern"),
                "levels": None if record.get("levels") is None else tuple(int(level) for level in record["levels"]),
                "max_tones": record.get("max_tones"),
                "blocks": record.get("blocks"),
            }

    refined: list[dict[str, Any]] = []
    for case in shortlist.values():
        print(f"[refine] {case_id(case)}", flush=True)
        refined.append(
            run_synthesis_trial(
                case,
                search_phase="refine",
                seed=17,
                init_guess="heuristic",
                maxiter=REFINE_MAXITER,
                multistart=1,
            )
        )
        if include_random_refine:
            refined.append(
                run_synthesis_trial(
                    case,
                    search_phase="refine",
                    seed=42,
                    init_guess="random",
                    maxiter=REFINE_MAXITER,
                    multistart=1,
                )
            )

    return {"screened": screened, "refined": refined}


def make_sqr_cases() -> list[dict[str, Any]]:
    cases: list[dict[str, Any]] = []
    for max_tones in (1, 2, 3, 4):
        for levels in c.ordered_level_subsets(PRELIM_TRUNC, max_tones):
            cases.append(
                {
                    "family_key": "drsqr",
                    "family_label": "D + R + SQR",
                    "variant_key": f"3sqr_t{max_tones}",
                    "variant_label": "3 selective blocks",
                    "builder_name": "drsqr_3sqr",
                    "pattern": "drs",
                    "levels": tuple(levels),
                    "max_tones": int(max_tones),
                    "blocks": 3,
                }
            )
    for max_tones in (1, 2, 3, 4):
        for levels in c.ordered_level_subsets(PRELIM_TRUNC, max_tones):
            cases.append(
                {
                    "family_key": "drsqr",
                    "family_label": "D + R + SQR",
                    "variant_key": f"4sqr_t{max_tones}",
                    "variant_label": "4 selective blocks",
                    "builder_name": "drsqr_4sqr",
                    "levels": tuple(levels),
                    "max_tones": int(max_tones),
                    "blocks": 4,
                }
            )
    return cases


def make_cpsqr_cases() -> list[dict[str, Any]]:
    cases: list[dict[str, Any]] = []
    for blocks in (2, 3, 4):
        for max_tones in (1, 2, 3, 4):
            for levels in c.ordered_level_subsets(PRELIM_TRUNC, max_tones):
                cases.append(
                    {
                        "family_key": "drcpsqr",
                        "family_label": "D + R + CPSQR",
                        "variant_key": f"cpsqr_blk{blocks}_t{max_tones}",
                        "variant_label": f"{blocks} conditional-phase blocks",
                        "builder_name": "drcpsqr",
                        "levels": tuple(levels),
                        "max_tones": int(max_tones),
                        "blocks": int(blocks),
                    }
                )
    return cases


def cavity_populations_from_state(state: qt.Qobj, *, n_cav: int) -> np.ndarray:
    if state.isket:
        vec = np.asarray(state.full()).reshape(-1)
        populations = np.zeros(int(n_cav), dtype=float)
        for level in range(int(n_cav)):
            populations[level] = float(np.abs(vec[level]) ** 2 + np.abs(vec[int(n_cav) + level]) ** 2)
        return populations
    rho = state if state.isoper else qt.ket2dm(state)
    populations = np.zeros(int(n_cav), dtype=float)
    for level in range(int(n_cav)):
        populations[level] = float(np.real(rho[level, level] + rho[int(n_cav) + level, int(n_cav) + level]))
    return populations


def logical_input_states(n_cav: int) -> list[tuple[str, qt.Qobj]]:
    labels = list(c.LOGICAL_LABELS)
    indices = c.logical_indices(int(n_cav))
    basis = [qt.basis(2 * int(n_cav), idx) for idx in indices]
    return list(zip(labels, basis))


def _flat_operator(op: qt.Qobj, *, n_cav: int) -> qt.Qobj:
    return qt.Qobj(op.full(), dims=[[2 * int(n_cav)], [2 * int(n_cav)]])


def _partial_gate_unitary(gate: Any, *, n_cav: int, fraction: float) -> qt.Qobj:
    if fraction <= 0.0:
        return _flat_operator(qt.qeye(2 * int(n_cav)), n_cav=int(n_cav))
    if fraction >= 1.0:
        return _flat_operator(gate.ideal_unitary(int(n_cav)), n_cav=int(n_cav))

    partial_gate = copy.deepcopy(gate)
    gate_type = type(gate).__name__
    if gate_type == "Displacement":
        params = np.asarray(gate.get_parameters(int(n_cav)), dtype=float)
        partial_gate.set_parameters(params * float(fraction), int(n_cav))
        partial_gate.duration = float(gate.duration) * float(fraction)
        return _flat_operator(partial_gate.ideal_unitary(int(n_cav)), n_cav=int(n_cav))
    if gate_type == "QubitRotation":
        params = np.asarray(gate.get_parameters(int(n_cav)), dtype=float)
        params[0] *= float(fraction)
        partial_gate.set_parameters(params, int(n_cav))
        partial_gate.duration = float(gate.duration) * float(fraction)
        return _flat_operator(partial_gate.ideal_unitary(int(n_cav)), n_cav=int(n_cav))
    if gate_type == "PrimitiveGate":
        metadata = getattr(gate, "metadata", {})
        params = copy.deepcopy(gate.parameters)
        kind = metadata.get("ideal_kind")
        if kind == "MaskedSQR":
            params["theta"] = np.asarray(params["theta"], dtype=float) * float(fraction)
        elif kind == "MaskedCPSQR":
            params["phases"] = np.asarray(params["phases"], dtype=float) * float(fraction)
            params["duration"] = float(params["duration"]) * float(fraction)
        else:
            return _flat_operator(gate.ideal_unitary(int(n_cav)), n_cav=int(n_cav))
        partial_gate.parameters = params
        partial_gate.duration = float(gate.duration) * float(fraction)
        return _flat_operator(partial_gate.ideal_unitary(int(n_cav)), n_cav=int(n_cav))
    return _flat_operator(gate.ideal_unitary(int(n_cav)), n_cav=int(n_cav))


def active_subspace_analysis(sequence: c.GateSequence, *, n_cav: int, sample_points: int = ACTIVE_SAMPLE_POINTS) -> dict[str, Any]:
    snapshots: list[dict[str, Any]] = []
    peak_by_level = np.zeros(int(n_cav), dtype=float)
    per_input_peak: dict[str, np.ndarray] = {}
    per_input_summaries: dict[str, Any] = {}

    for label, state in logical_input_states(int(n_cav)):
        running = state
        local_peak = np.zeros(int(n_cav), dtype=float)
        for gate_index, gate in enumerate(sequence.gates):
            for fraction in np.linspace(1.0 / sample_points, 1.0, int(sample_points)):
                partial = _partial_gate_unitary(gate, n_cav=int(n_cav), fraction=float(fraction))
                evolved = partial * running
                populations = cavity_populations_from_state(evolved, n_cav=int(n_cav))
                local_peak = np.maximum(local_peak, populations)
                peak_by_level = np.maximum(peak_by_level, populations)
                snapshots.append(
                    {
                        "input": label,
                        "gate_index": int(gate_index),
                        "fraction": float(fraction),
                        "populations": populations,
                    }
                )
            running = _flat_operator(gate.ideal_unitary(int(n_cav)), n_cav=int(n_cav)) * running

        initial_active = {idx for idx, pop in enumerate(local_peak) if pop >= ACTIVE_THRESHOLD}
        ordered_levels = list(np.argsort(local_peak)[::-1])
        active = set(initial_active)
        if not active and ordered_levels:
            active.add(int(ordered_levels[0]))
        capture = 0.0
        for level in ordered_levels:
            if capture >= ACTIVE_CAPTURE_TARGET:
                break
            active.add(int(level))
            capture = float(np.sum(local_peak[list(sorted(active))]))
        worst_captured = 1.0
        for snapshot in snapshots:
            if snapshot["input"] != label:
                continue
            captured = float(np.sum(snapshot["populations"][list(sorted(active))]))
            worst_captured = min(worst_captured, captured)
            if worst_captured >= ACTIVE_CAPTURE_TARGET:
                continue
            for level in ordered_levels:
                active.add(int(level))
                captured = float(np.sum(snapshot["populations"][list(sorted(active))]))
                worst_captured = min(worst_captured, captured)
                if captured >= ACTIVE_CAPTURE_TARGET:
                    break

        per_input_peak[label] = local_peak.copy()
        per_input_summaries[label] = {
            "active_levels": [int(level) for level in sorted(active)],
            "peak_population_by_level": local_peak.tolist(),
            "worst_captured_population": float(worst_captured),
            "max_active_level": int(max(active)) if active else -1,
        }

    candidate_active = {idx for idx, pop in enumerate(peak_by_level) if pop >= ACTIVE_THRESHOLD}
    ordered_candidate_levels = list(np.argsort(peak_by_level)[::-1])
    if not candidate_active and ordered_candidate_levels:
        candidate_active.add(int(ordered_candidate_levels[0]))
    worst_candidate_capture = 1.0
    for snapshot in snapshots:
        captured = float(np.sum(snapshot["populations"][list(sorted(candidate_active))]))
        while captured < ACTIVE_CAPTURE_TARGET:
            added = False
            for level in ordered_candidate_levels:
                if int(level) not in candidate_active:
                    candidate_active.add(int(level))
                    captured = float(np.sum(snapshot["populations"][list(sorted(candidate_active))]))
                    added = True
                    break
            if not added:
                break
        worst_candidate_capture = min(worst_candidate_capture, captured)

    boundary_touch = bool(candidate_active and max(candidate_active) >= int(n_cav) - 1)
    return {
        "threshold": ACTIVE_THRESHOLD,
        "capture_target": ACTIVE_CAPTURE_TARGET,
        "sample_points_per_gate": int(sample_points),
        "candidate_active_levels": [int(level) for level in sorted(candidate_active)],
        "candidate_peak_population_by_level": peak_by_level.tolist(),
        "candidate_worst_captured_population": float(worst_candidate_capture),
        "candidate_max_active_level": int(max(candidate_active)) if candidate_active else -1,
        "touches_truncation_boundary": boundary_touch,
        "per_input": per_input_summaries,
    }


def classify_candidate(evaluations: dict[str, Any]) -> tuple[str, list[str]]:
    notes: list[str] = []
    default_eval = evaluations[str(DEFAULT_FINAL_TRUNC)]
    fidelity_12 = float(default_eval["fidelity"])
    leakage_12 = float(default_eval["leakage_worst"])
    active_12 = default_eval["active_subspace"]
    delta_12_14 = abs(float(evaluations[str(14)]["fidelity"]) - fidelity_12)

    if fidelity_12 < RETENTION_FIDELITY_MIN:
        notes.append(f"N_cav=12 fidelity {fidelity_12:.4f} is below the retained threshold {RETENTION_FIDELITY_MIN:.2f}.")
    if leakage_12 > LOW_CONFIDENCE_LEAKAGE_MAX:
        notes.append(f"N_cav=12 worst leakage {leakage_12:.4f} is above {LOW_CONFIDENCE_LEAKAGE_MAX:.2f}.")
    if active_12["touches_truncation_boundary"]:
        notes.append("Active cavity subspace reaches the truncation boundary at N_cav=12.")
    if delta_12_14 > LOW_CONFIDENCE_CONVERGENCE_DELTA:
        notes.append(f"Fidelity drift |F_14-F_12|={delta_12_14:.4f} exceeds {LOW_CONFIDENCE_CONVERGENCE_DELTA:.2f}.")

    if not notes and leakage_12 <= RETENTION_LEAKAGE_MAX and delta_12_14 <= RETENTION_CONVERGENCE_DELTA:
        return "retained", [f"Retained at N_cav=12 with fidelity {fidelity_12:.4f} and leakage {leakage_12:.4f}."]

    if fidelity_12 >= RETENTION_FIDELITY_MIN and leakage_12 <= LOW_CONFIDENCE_LEAKAGE_MAX:
        if not notes:
            notes.append("Candidate passes N_cav=12 fidelity but fails a stricter convergence or support check.")
        return "low_confidence", notes

    if not notes:
        notes.append("Candidate fails the corrected physical-retention rules.")
    return "discarded", notes


def evaluate_physical_candidate(record: dict[str, Any]) -> dict[str, Any]:
    evaluations: dict[str, Any] = {}
    for n_cav in FINAL_TRUNCATIONS:
        sequence = apply_solution_to_case(record, n_cav=int(n_cav))
        evaluation = c.evaluate_sequence(sequence, n_cav=int(n_cav))
        active = active_subspace_analysis(sequence, n_cav=int(n_cav))
        evaluations[str(n_cav)] = {
            "fidelity": float(evaluation["fidelity"]),
            "leakage_average": float(evaluation["leakage_average"]),
            "leakage_worst": float(evaluation["leakage_worst"]),
            "unitarity_error": float(evaluation["unitarity_error"]),
            "active_subspace": active,
        }
    status, notes = classify_candidate(evaluations)
    return {
        "status": status,
        "notes": notes,
        "by_n_cav": evaluations,
    }


def reduced_cavity_density(state: qt.Qobj, *, n_cav: int) -> qt.Qobj:
    if state.isket:
        ket = state
        if ket.dims != [[2, int(n_cav)], [1]] and ket.dims != [[2, int(n_cav)], [1, 1]]:
            ket = qt.Qobj(ket.full(), dims=[[2, int(n_cav)], [1, 1]])
        rho = qt.ket2dm(ket)
    else:
        rho = state
        if rho.dims != [[2, int(n_cav)], [2, int(n_cav)]]:
            rho = qt.Qobj(rho.full(), dims=[[2, int(n_cav)], [2, int(n_cav)]])
    return rho.ptrace(1)


def target_output_state(label: str, *, n_cav: int) -> qt.Qobj:
    index = list(c.LOGICAL_LABELS).index(label)
    target_column = c.TARGET_UNITARY[:, index]
    vec = np.zeros(2 * int(n_cav), dtype=np.complex128)
    indices = c.logical_indices(int(n_cav))
    for offset, amplitude in zip(indices, target_column):
        vec[int(offset)] = amplitude
    return qt.Qobj(vec)


def candidate_output_state(record: dict[str, Any], *, n_cav: int, input_label: str) -> qt.Qobj:
    sequence = apply_solution_to_case(record, n_cav=int(n_cav))
    state = dict(logical_input_states(int(n_cav)))[str(input_label)]
    for gate in sequence.gates:
        state = _flat_operator(gate.ideal_unitary(int(n_cav)), n_cav=int(n_cav)) * state
    return state


def save_record_artifact(record: dict[str, Any], *, stem: str) -> str:
    payload = {
        "study_name": "cluster_state_holographic_unified",
        "date_created": time.strftime("%Y-%m-%d"),
        "description": f"Corrected-scope best candidate: {record['family_label']} / {record['variant_label']}",
        "parameters": {
            "family_key": record["family_key"],
            "variant_key": record["variant_key"],
            "levels": record.get("levels"),
            "blocks": record.get("blocks"),
            "max_tones": record.get("max_tones"),
            "search_phase": record.get("search_phase"),
            "seed": record.get("seed"),
        },
        "load_instructions": "import json; from pathlib import Path; payload = json.loads(Path(filename).read_text(encoding='utf-8'))",
        "record": record_without_fit(record),
    }
    path = c.ARTIFACT_DIR / f"{stem}.json"
    c.save_json(path, payload)
    return str(path)


def load_grape_reference() -> dict[str, Any]:
    summary = c.load_json(c.DATA_DIR / "consolidated_summary.json")
    return {
        "family_label": "GRAPE reference",
        "duration_ns": int(summary["grape_physical_validation"]["duration_ns"]),
        "n_cav_optimization": int(summary["grape_physical_validation"]["n_cav_optimization"]),
        "replay_fidelity": float(summary["grape_physical_validation"]["replay_fidelity"]),
        "open_system_process_fidelity": float(summary["grape_physical_validation"]["open_system_process_fidelity"]),
    }


def physical_eval(record: dict[str, Any], n_cav: int) -> dict[str, Any]:
    return record["physical"]["by_n_cav"][str(int(n_cav))]


def plot_family_tradeoff(summary: dict[str, Any]) -> None:
    prelim = summary["preliminary_best"]
    retained = summary["retained_candidates"]
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0))

    prelim_rows = [prelim["drsqr"], prelim["drcpsqr"]]
    labels = [row["family_label"] for row in prelim_rows]
    x = np.arange(len(labels))
    axes[0].bar(x - 0.18, [row["fidelity"] for row in prelim_rows], width=0.36, label=f"preliminary N={PRELIM_TRUNC}")
    axes[0].bar(x + 0.18, [physical_eval(row, DEFAULT_FINAL_TRUNC)["fidelity"] for row in prelim_rows], width=0.36, label=f"physical N={DEFAULT_FINAL_TRUNC}")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels, rotation=10)
    axes[0].set_ylabel("Fidelity")
    axes[0].set_ylim(0.0, 1.05)
    axes[0].set_title("Preliminary vs retained physical fidelity")
    axes[0].legend(frameon=False)

    retained_rows = [retained.get("drsqr"), retained.get("drcpsqr")]
    retained_rows = [row for row in retained_rows if row is not None]
    if retained_rows:
        labels2 = [row["family_label"] for row in retained_rows]
        x2 = np.arange(len(labels2))
        axes[1].bar(x2 - 0.2, [physical_eval(row, 10)["fidelity"] for row in retained_rows], width=0.2, label="N=10")
        axes[1].bar(x2, [physical_eval(row, 12)["fidelity"] for row in retained_rows], width=0.2, label="N=12")
        axes[1].bar(x2 + 0.2, [physical_eval(row, 14)["fidelity"] for row in retained_rows], width=0.2, label="N=14")
        axes[1].set_xticks(x2)
        axes[1].set_xticklabels(labels2, rotation=10)
        axes[1].set_ylabel("Fidelity")
        axes[1].set_ylim(0.0, 1.05)
        axes[1].set_title("Shortlisted truncation convergence")
        axes[1].legend(frameon=False)
    else:
        axes[1].axis("off")

    fig.tight_layout()
    fig.savefig(c.FIG_DIR / "corrected_family_tradeoff.pdf", bbox_inches="tight")
    fig.savefig(c.FIG_DIR / "corrected_family_tradeoff.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_active_subspace(summary: dict[str, Any]) -> None:
    rows = []
    for family_key in ("drsqr", "drcpsqr"):
        candidate = summary["preliminary_best"][family_key]
        physical = physical_eval(candidate, DEFAULT_FINAL_TRUNC)
        active = physical["active_subspace"]
        rows.append(
            {
                "label": candidate["family_label"],
                "max_level": active["candidate_max_active_level"],
                "captured": active["candidate_worst_captured_population"],
                "leakage": physical["leakage_worst"],
            }
        )

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.0))
    x = np.arange(len(rows))
    axes[0].bar(x, [row["max_level"] for row in rows], color=["#4477AA", "#EE6677"])
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([row["label"] for row in rows], rotation=10)
    axes[0].set_ylabel("Max active Fock level")
    axes[0].set_title(f"Active subspace at N={DEFAULT_FINAL_TRUNC}")

    axes[1].bar(x - 0.18, [row["captured"] for row in rows], width=0.36, label="captured population")
    axes[1].bar(x + 0.18, [1.0 - row["leakage"] for row in rows], width=0.36, label="1 - worst leakage")
    axes[1].axhline(ACTIVE_CAPTURE_TARGET, color="0.3", linestyle="--", linewidth=1.0)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([row["label"] for row in rows], rotation=10)
    axes[1].set_ylim(0.8, 1.01)
    axes[1].set_ylabel("Population fraction")
    axes[1].set_title("Support coverage and leakage")
    axes[1].legend(frameon=False)

    fig.tight_layout()
    fig.savefig(c.FIG_DIR / "corrected_active_subspace_summary.pdf", bbox_inches="tight")
    fig.savefig(c.FIG_DIR / "corrected_active_subspace_summary.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_wigner_panels(summary: dict[str, Any]) -> None:
    selected_rows: list[tuple[str, Any]] = [("Target", None)]
    sqr = summary["retained_candidates"].get("drsqr")
    cpsqr = summary["retained_candidates"].get("drcpsqr")
    if sqr is not None:
        selected_rows.append(("Best D + R + SQR", sqr))
    if cpsqr is not None:
        selected_rows.append(("Best D + R + CPSQR", cpsqr))

    xvec = np.linspace(-3.0, 3.0, 151)
    labels = list(c.LOGICAL_LABELS)
    fig = plt.figure(figsize=(9.6, 9.8), constrained_layout=True)
    grid = fig.add_gridspec(
        len(labels),
        len(selected_rows) + 1,
        width_ratios=[1.0] * len(selected_rows) + [0.08],
        wspace=0.08,
        hspace=0.10,
    )
    axes = np.empty((len(labels), len(selected_rows)), dtype=object)
    vmax = 0.35

    shared_axis = None
    for row_index, input_label in enumerate(labels):
        for col_index, (title, record) in enumerate(selected_rows):
            if shared_axis is None:
                ax = fig.add_subplot(grid[row_index, col_index])
                shared_axis = ax
            else:
                ax = fig.add_subplot(grid[row_index, col_index], sharex=shared_axis, sharey=shared_axis)
            axes[row_index, col_index] = ax
            if record is None:
                state = target_output_state(input_label, n_cav=DEFAULT_FINAL_TRUNC)
            else:
                state = candidate_output_state(record, n_cav=DEFAULT_FINAL_TRUNC, input_label=input_label)
            rho_cav = reduced_cavity_density(state, n_cav=DEFAULT_FINAL_TRUNC)
            wig = qt.wigner(rho_cav, xvec, xvec)
            mesh = ax.pcolormesh(xvec, xvec, wig, cmap="RdBu_r", shading="auto", vmin=-vmax, vmax=vmax, rasterized=True)
            ax.set_aspect("equal")
            if row_index == 0:
                ax.set_title(title)
            if col_index == 0:
                ax.set_ylabel(f"{input_label}\nIm(α)")
            if row_index == len(labels) - 1:
                ax.set_xlabel("Re(α)")

    for ax in axes.ravel():
        ax.label_outer()

    cax = fig.add_subplot(grid[:, -1])
    colorbar = fig.colorbar(mesh, cax=cax)
    colorbar.set_label("W(α)")
    cax.yaxis.set_ticks_position("right")
    cax.yaxis.set_label_position("right")

    fig.savefig(c.FIG_DIR / "appendix_wigner.pdf", bbox_inches="tight")
    fig.savefig(c.FIG_DIR / "appendix_wigner.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def run_searches() -> dict[str, Any]:
    sqr_runs = screen_and_refine_cases(cases=make_sqr_cases(), shortlist_per_group=2, group_key="variant_key")
    cpsqr_runs = screen_and_refine_cases(cases=make_cpsqr_cases(), shortlist_per_group=1, group_key="blocks")

    sqr_pool = sqr_runs["refined"] if sqr_runs["refined"] else sqr_runs["screened"]
    cpsqr_pool = cpsqr_runs["refined"] if cpsqr_runs["refined"] else cpsqr_runs["screened"]
    sqr_best = best_record(sqr_pool)
    cpsqr_best = best_record(cpsqr_pool)

    sqr_best["physical"] = evaluate_physical_candidate(sqr_best)
    cpsqr_best["physical"] = evaluate_physical_candidate(cpsqr_best)

    sqr_best["artifact_path"] = save_record_artifact(sqr_best, stem="corrected_best_sqr")
    cpsqr_best["artifact_path"] = save_record_artifact(cpsqr_best, stem="corrected_best_cpsqr")

    retained_candidates = {
        "drsqr": sqr_best if sqr_best["physical"]["status"] != "discarded" else None,
        "drcpsqr": cpsqr_best if cpsqr_best["physical"]["status"] != "discarded" else None,
    }

    summary = {
        "study_name": "cluster_state_holographic_unified",
        "date_created": time.strftime("%Y-%m-%d"),
        "scope": {
            "exclusive_families": ["D + R + SQR", "D + R + CPSQR"],
            "grape_role": "separate reference",
            "default_final_truncation": DEFAULT_FINAL_TRUNC,
            "truncation_checks": list(FINAL_TRUNCATIONS),
        },
        "active_subspace_definition": {
            "threshold": ACTIVE_THRESHOLD,
            "capture_target": ACTIVE_CAPTURE_TARGET,
            "rule": (
                "For each candidate and each logical input, sample the full sequence gate-by-gate with intermediate fractional steps. "
                "Define the active cavity subspace as the set of Fock levels whose peak population is at least 1e-3, then expand that set until the worst instantaneous captured population is at least 99.9%."
            ),
            "candidate_reporting": "Report per-input active levels and the candidate-wide union at N_cav = 10, 12, 14.",
        },
        "preliminary_search": {
            "drsqr": {
                "screened": [record_without_fit(row) for row in sqr_runs["screened"]],
                "refined": [record_without_fit(row) for row in sqr_runs["refined"]],
            },
            "drcpsqr": {
                "screened": [record_without_fit(row) for row in cpsqr_runs["screened"]],
                "refined": [record_without_fit(row) for row in cpsqr_runs["refined"]],
            },
        },
        "preliminary_best": {
            "drsqr": record_without_fit(sqr_best),
            "drcpsqr": record_without_fit(cpsqr_best),
        },
        "retained_candidates": {
            key: None if value is None else record_without_fit(value)
            for key, value in retained_candidates.items()
        },
        "discarded_candidates": {
            "drsqr": None if retained_candidates["drsqr"] is not None else record_without_fit(sqr_best),
            "drcpsqr": None if retained_candidates["drcpsqr"] is not None else record_without_fit(cpsqr_best),
        },
        "grape_reference": load_grape_reference(),
    }
    c.save_json(c.DATA_DIR / "corrected_scope_summary.json", summary)
    return summary


def write_short_change_summary(summary: dict[str, Any]) -> None:
    sqr = summary["preliminary_best"]["drsqr"]
    cpsqr = summary["preliminary_best"]["drcpsqr"]
    lines = [
        "# What Changed From the Previous Version",
        "",
        "- Restricted the gate-family comparison to exclusive D + R + SQR and D + R + CPSQR families.",
        "- Re-ranked candidates using N_cav = 10, 12, 14 with N_cav = 12 as the default final evidence.",
        "- Added an explicit active-subspace rule based on sampled full-evolution cavity populations.",
        f"- Best preliminary D + R + SQR candidate: fidelity {physical_eval(sqr, DEFAULT_FINAL_TRUNC)['fidelity']:.4f} at N_cav = {DEFAULT_FINAL_TRUNC}.",
        f"- Best preliminary D + R + CPSQR candidate: fidelity {physical_eval(cpsqr, DEFAULT_FINAL_TRUNC)['fidelity']:.4f} at N_cav = {DEFAULT_FINAL_TRUNC}.",
        "- Removed mixed SQR/CPSQR families from the final narrative, tables, and conclusions.",
        "- Replaced the appendix Wigner figure with a local corrected-scope figure generated in this study's figures/ directory.",
    ]
    (c.DATA_DIR / "what_changed_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Rerun the unified holographic study under the corrected scope restrictions.")
    parser.add_argument("--skip-search", action="store_true", help="Reuse corrected_scope_summary.json if it already exists.")
    args = parser.parse_args()

    summary_path = c.DATA_DIR / "corrected_scope_summary.json"
    if args.skip_search and summary_path.exists():
        summary = c.load_json(summary_path)
    else:
        summary = run_searches()
    plot_family_tradeoff(summary)
    plot_active_subspace(summary)
    plot_wigner_panels(summary)
    write_short_change_summary(summary)
    print("[done] corrected_scope_summary.json and corrected figures written", flush=True)


if __name__ == "__main__":
    main()