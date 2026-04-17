from __future__ import annotations

import argparse
import os
import copy
import csv
import itertools
import statistics
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np
import runtime_compat  # noqa: F401
import qutip as qt


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import common as c


STYLE_PATH = SCRIPT_DIR.parents[2] / ".github" / "skills" / "publication-figures" / "assets" / "cqed_style.mplstyle"
if STYLE_PATH.exists():
    plt.style.use(str(STYLE_PATH))


SCREEN_MAXITER = 2
PRELIM_REFINE_MAXITER = 20
PHYSICAL_REFINE_MAXITER = 36
ACTIVE_SAMPLE_POINTS = 12
PRELIM_TRUNC = c.DECOMP_N_CAV
FINAL_TRUNCATIONS = (10, 12, 14)
DEFAULT_FINAL_TRUNC = 12
TARGET_FIDELITY = 0.99
ACTIVE_THRESHOLD = 1.0e-3
ACTIVE_CAPTURE_TARGET = 0.999
RETENTION_FIDELITY_MIN = 0.90
RETENTION_LEAKAGE_MAX = 0.05
LOW_CONFIDENCE_LEAKAGE_MAX = 0.10
RETENTION_CONVERGENCE_DELTA = 0.01
LOW_CONFIDENCE_CONVERGENCE_DELTA = 0.05
STRUCTURAL_FINALISTS_PER_FAMILY = 8
PHYSICAL_FINALISTS_PER_FAMILY = 6
MAX_WORKERS = max(1, min(8, os.cpu_count() or 1))
SQR_BLOCK_OPTIONS = (2, 3, 4, 5)
CPSQR_BLOCK_OPTIONS = (1, 2, 3, 4, 5)
N_ACTIVE_OPTIONS = tuple(range(1, PRELIM_TRUNC + 1))

SQR_ORDERINGS = tuple(itertools.permutations(("D", "R", "SQR")))
CPSQR_ORDERINGS = tuple(itertools.permutations(("D", "R", "CPSQR")))

ORDER_ABBREVIATIONS = {"D": "D", "R": "R", "SQR": "S", "CPSQR": "CP"}

BUILDERS = {
    "ordered_sqr": c.build_ordered_sqr_sequence,
    "ordered_cpsqr": c.build_ordered_cpsqr_sequence,
}


def compact_order_label(order: Sequence[str]) -> str:
    return "".join(ORDER_ABBREVIATIONS[str(token).upper()] for token in order)


def pretty_order_label(order: Sequence[str]) -> str:
    return " > ".join(str(token).upper() for token in order)


def case_id(case: dict[str, Any]) -> str:
    parts = [str(case["family_key"]), str(case["variant_key"])]
    if case.get("order_label"):
        parts.append(str(case["order_label"]))
    if case.get("levels"):
        parts.append("lv" + "-".join(str(level) for level in case["levels"]))
    if case.get("blocks") is not None:
        parts.append(f"blk{int(case['blocks'])}")
    if case.get("max_tones") is not None:
        parts.append(f"t{int(case['max_tones'])}")
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


def replay_sort_key(record: dict[str, Any]) -> tuple[float, float, float, float, float, float]:
    replay = record["replay"][str(DEFAULT_FINAL_TRUNC)]
    summary = record["summary"]
    return (
        float(replay["fidelity"]),
        -float(replay["leakage_worst"]),
        float(record["fidelity"]),
        -float(summary["total_duration_ns"]),
        -float(summary["gate_depth"]),
        -float(summary["max_active_tones"]),
    )


def best_record(records: Iterable[dict[str, Any]], *, key=ranking_key) -> dict[str, Any]:
    return max(records, key=key)


def record_without_fit(record: dict[str, Any]) -> dict[str, Any]:
    return c.json_ready({key: value for key, value in record.items() if key not in {"fit", "warm_start_payload"}})


def case_from_record(record: dict[str, Any]) -> dict[str, Any]:
    return {
        "family_key": str(record["family_key"]),
        "family_label": str(record["family_label"]),
        "variant_key": str(record["variant_key"]),
        "variant_label": str(record["variant_label"]),
        "builder_name": str(record["builder_name"]),
        "builder_kwargs": dict(record.get("builder_kwargs", {})),
        "order_tokens": tuple(record.get("order_tokens", [])),
        "order_label": str(record.get("order_label", "")),
        "levels": None if record.get("levels") is None else tuple(int(level) for level in record["levels"]),
        "max_tones": None if record.get("max_tones") is None else int(record["max_tones"]),
        "blocks": None if record.get("blocks") is None else int(record["blocks"]),
    }


def build_sequence_from_case(case: dict[str, Any], *, n_cav: int) -> c.GateSequence:
    builder = BUILDERS[str(case["builder_name"])]
    kwargs = dict(case.get("builder_kwargs", {}))
    kwargs["n_cav"] = int(n_cav)
    if case.get("levels") is not None:
        kwargs["levels"] = tuple(int(level) for level in case["levels"])
    if case.get("blocks") is not None:
        kwargs["blocks"] = int(case["blocks"])
    if case.get("order_tokens"):
        kwargs["order"] = tuple(str(token) for token in case["order_tokens"])
    return builder(**kwargs)


def apply_solution_to_case(record: dict[str, Any], *, n_cav: int) -> c.GateSequence:
    sequence = build_sequence_from_case(record, n_cav=int(n_cav))
    sequence.set_parameter_vector(np.asarray(record["parameter_vector"], dtype=float))
    sequence.set_time_vector(np.asarray(record["time_vector"], dtype=float), active_only=False)
    return sequence


def run_synthesis_trial(
    case: dict[str, Any],
    *,
    n_cav: int,
    search_phase: str,
    seed: int,
    init_guess: str,
    maxiter: int,
    multistart: int = 1,
    warm_start: Any | None = None,
) -> dict[str, Any]:
    started = time.perf_counter()
    sequence = build_sequence_from_case(case, n_cav=int(n_cav))
    fit = c.fit_sequence(
        sequence,
        n_cav=int(n_cav),
        seed=int(seed),
        init_guess=str(init_guess),
        multistart=int(multistart),
        maxiter=int(maxiter),
        duration_weight=float(case.get("duration_weight", 0.0)),
        gate_count_weight=float(case.get("gate_count_weight", 0.0)),
        warm_start=warm_start,
    )
    record = {
        "case_id": case_id(case),
        "family_key": str(case["family_key"]),
        "family_label": str(case["family_label"]),
        "variant_key": str(case["variant_key"]),
        "variant_label": str(case["variant_label"]),
        "builder_name": str(case["builder_name"]),
        "builder_kwargs": dict(case.get("builder_kwargs", {})),
        "order_tokens": [str(token) for token in case.get("order_tokens", ())],
        "order_label": str(case.get("order_label", "")),
        "levels": None if case.get("levels") is None else [int(level) for level in case["levels"]],
        "max_tones": None if case.get("max_tones") is None else int(case["max_tones"]),
        "blocks": None if case.get("blocks") is None else int(case["blocks"]),
        "search_phase": str(search_phase),
        "optimization_n_cav": int(n_cav),
        "seed": int(seed),
        "init_guess": str(init_guess),
        "maxiter": int(maxiter),
        "multistart": int(multistart),
        "fidelity": float(fit["fidelity"]),
        "objective": float(fit["objective"]),
        "success": bool(fit["success"]),
        "message": str(fit["message"]),
        "summary": dict(fit["summary"]),
        "metrics": dict(fit["metrics"]),
        "sequence": fit["sequence_payload"],
        "parameter_vector": fit["result"].sequence.get_parameter_vector().tolist(),
        "time_vector": fit["result"].sequence.get_time_vector(active_only=False).tolist(),
        "warm_start_payload": fit["result"].to_payload(include_history=False),
        "elapsed_s": float(time.perf_counter() - started),
    }
    return record


def _trial_worker(job: dict[str, Any]) -> dict[str, Any]:
    record_updates = dict(job.get("record_updates", {}))
    record = run_synthesis_trial(
        job["case"],
        n_cav=int(job["n_cav"]),
        search_phase=str(job["search_phase"]),
        seed=int(job["seed"]),
        init_guess=str(job["init_guess"]),
        maxiter=int(job["maxiter"]),
        multistart=int(job.get("multistart", 1)),
        warm_start=job.get("warm_start"),
    )
    record.update(record_updates)
    return record


def multiprocessing_available() -> bool:
    main_module = sys.modules.get("__main__")
    main_file = getattr(main_module, "__file__", "")
    return bool(main_file) and "<" not in str(main_file)


def run_trial_batch(jobs: list[dict[str, Any]], *, progress_label: str) -> list[dict[str, Any]]:
    if not jobs:
        return []

    total = len(jobs)
    if MAX_WORKERS <= 1 or total == 1 or not multiprocessing_available():
        records: list[dict[str, Any]] = []
        for index, job in enumerate(jobs, start=1):
            print(f"[{progress_label} {index}/{total}] {case_id(job['case'])}", flush=True)
            records.append(_trial_worker(job))
        return records

    records = []
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(_trial_worker, job): job for job in jobs}
        completed = 0
        for future in as_completed(futures):
            job = futures[future]
            try:
                record = future.result()
            except Exception as exc:
                raise RuntimeError(f"{progress_label} failed for {case_id(job['case'])}") from exc
            completed += 1
            print(f"[{progress_label} {completed}/{total}] {record['case_id']}", flush=True)
            records.append(record)
    return records


def design_space_key(record: dict[str, Any]) -> tuple[str, int, int, str]:
    return (
        str(record["family_key"]),
        int(record["blocks"]),
        int(record["max_tones"]),
        str(record["order_label"]),
    )


def replay_record(record: dict[str, Any], *, n_cav_values: Sequence[int] = (DEFAULT_FINAL_TRUNC,)) -> dict[str, Any]:
    replay: dict[str, Any] = {}
    for n_cav in n_cav_values:
        evaluation = c.evaluate_sequence(apply_solution_to_case(record, n_cav=int(n_cav)), n_cav=int(n_cav))
        replay[str(int(n_cav))] = {
            "fidelity": float(evaluation["fidelity"]),
            "leakage_average": float(evaluation["leakage_average"]),
            "leakage_worst": float(evaluation["leakage_worst"]),
            "unitarity_error": float(evaluation["unitarity_error"]),
        }
    record["replay"] = replay
    return record


def screen_level_choices(max_tones: int) -> list[tuple[int, ...]]:
    return [tuple(range(int(max_tones)))]


def refinement_level_choices(max_tones: int) -> list[tuple[int, ...]]:
    return [tuple(levels) for levels in c.ordered_level_subsets(PRELIM_TRUNC, int(max_tones))]


def make_sqr_cases() -> list[dict[str, Any]]:
    cases: list[dict[str, Any]] = []
    for blocks in SQR_BLOCK_OPTIONS:
        for order in SQR_ORDERINGS:
            order_label = compact_order_label(order)
            for max_tones in N_ACTIVE_OPTIONS:
                for levels in screen_level_choices(max_tones):
                    cases.append(
                        {
                            "family_key": "drsqr",
                            "family_label": "D + R + SQR",
                            "variant_key": f"blk{blocks}_a{max_tones}_ord{order_label}",
                            "variant_label": f"{blocks} blocks / {order_label} / n_active={max_tones}",
                            "builder_name": "ordered_sqr",
                            "builder_kwargs": {"blocks": int(blocks), "order": tuple(order)},
                            "order_tokens": tuple(order),
                            "order_label": order_label,
                            "levels": tuple(levels),
                            "max_tones": int(max_tones),
                            "blocks": int(blocks),
                        }
                    )
    return cases


def make_cpsqr_cases() -> list[dict[str, Any]]:
    cases: list[dict[str, Any]] = []
    for blocks in CPSQR_BLOCK_OPTIONS:
        for order in CPSQR_ORDERINGS:
            order_label = compact_order_label(order)
            for max_tones in N_ACTIVE_OPTIONS:
                for levels in screen_level_choices(max_tones):
                    cases.append(
                        {
                            "family_key": "drcpsqr",
                            "family_label": "D + R + CPSQR",
                            "variant_key": f"blk{blocks}_a{max_tones}_ord{order_label}",
                            "variant_label": f"{blocks} blocks / {order_label} / n_active={max_tones}",
                            "builder_name": "ordered_cpsqr",
                            "builder_kwargs": {"blocks": int(blocks), "order": tuple(order)},
                            "order_tokens": tuple(order),
                            "order_label": order_label,
                            "levels": tuple(levels),
                            "max_tones": int(max_tones),
                            "blocks": int(blocks),
                        }
                    )
    return cases


def screen_design_space(cases: list[dict[str, Any]]) -> dict[str, Any]:
    best_by_config: dict[tuple[str, int, int, str], dict[str, Any]] = {}
    family_counts: dict[str, int] = {}
    jobs = [
        {
            "case": case,
            "n_cav": PRELIM_TRUNC,
            "search_phase": "screen",
            "seed": 17,
            "init_guess": "heuristic",
            "maxiter": SCREEN_MAXITER,
            "multistart": 1,
        }
        for case in cases
    ]
    for record in run_trial_batch(jobs, progress_label="screen"):
        family_counts[str(record["family_key"])] = family_counts.get(str(record["family_key"]), 0) + 1
        key = design_space_key(record)
        current = best_by_config.get(key)
        if current is None or ranking_key(record) > ranking_key(current):
            best_by_config[key] = record
    return {
        "case_count": len(cases),
        "family_counts": family_counts,
        "best_by_config": list(best_by_config.values()),
    }


def select_finalists(records: list[dict[str, Any]], *, count: int) -> list[dict[str, Any]]:
    ordered = sorted(records, key=replay_sort_key, reverse=True)
    selected: dict[str, dict[str, Any]] = {}
    for field in ("blocks", "max_tones", "order_label"):
        best_by_value: dict[Any, dict[str, Any]] = {}
        for record in ordered:
            best_by_value.setdefault(record[field], record)
        for record in best_by_value.values():
            selected.setdefault(str(record["case_id"]), record)
    for record in ordered:
        if len(selected) >= int(count):
            break
        selected.setdefault(str(record["case_id"]), record)
    finalists = sorted(selected.values(), key=replay_sort_key, reverse=True)
    return finalists[: int(count)]


def refine_level_subsets(records: list[dict[str, Any]], *, maxiter: int) -> list[dict[str, Any]]:
    expanded: dict[str, tuple[dict[str, Any], Any | None, str]] = {}
    for record in records:
        base_case = case_from_record(record)
        warm_start = record.get("warm_start_payload")
        for levels in refinement_level_choices(int(record["max_tones"])):
            case = dict(base_case)
            case["levels"] = tuple(levels)
            expanded.setdefault(case_id(case), (case, warm_start, str(record["case_id"])))

    jobs = [
        {
            "case": case,
            "n_cav": PRELIM_TRUNC,
            "search_phase": "level_refine",
            "seed": 17,
            "init_guess": "heuristic",
            "maxiter": int(maxiter),
            "multistart": 1,
            "warm_start": warm_start,
            "record_updates": {"structure_seed": structure_seed},
        }
        for case, warm_start, structure_seed in expanded.values()
    ]
    return run_trial_batch(jobs, progress_label="level_refine")


def refine_finalists(records: list[dict[str, Any]], *, n_cav: int, phase: str, maxiter: int) -> list[dict[str, Any]]:
    jobs = []
    for record in records:
        record_updates: dict[str, Any] = {}
        if record.get("replay"):
            record_updates["replay_seed"] = c.json_ready(record["replay"])
        jobs.append(
            {
                "case": case_from_record(record),
                "n_cav": int(n_cav),
                "search_phase": str(phase),
                "seed": 17,
                "init_guess": "heuristic",
                "maxiter": int(maxiter),
                "multistart": 1,
                "warm_start": record.get("warm_start_payload"),
                "record_updates": record_updates,
            }
        )
    return run_trial_batch(jobs, progress_label=str(phase))


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
    return {"status": status, "notes": notes, "by_n_cav": evaluations}


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
        "description": f"Best comprehensive candidate: {record['family_label']} / {record['variant_label']}",
        "parameters": {
            "family_key": record["family_key"],
            "variant_key": record["variant_key"],
            "levels": record.get("levels"),
            "blocks": record.get("blocks"),
            "max_tones": record.get("max_tones"),
            "order_tokens": record.get("order_tokens"),
            "order_label": record.get("order_label"),
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


def _median(values: Sequence[float]) -> float:
    return float(statistics.median(values)) if values else float("nan")


def slim_design_space_row(record: dict[str, Any]) -> dict[str, Any]:
    replay = record.get("replay", {}).get(str(DEFAULT_FINAL_TRUNC), {})
    return {
        "family_key": str(record["family_key"]),
        "family_label": str(record["family_label"]),
        "blocks": int(record["blocks"]),
        "n_active": int(record["max_tones"]),
        "order_label": str(record["order_label"]),
        "order_pretty": pretty_order_label(record.get("order_tokens", [])),
        "levels": "-".join(str(level) for level in record.get("levels", [])),
        "screen_fidelity": float(record["fidelity"]),
        "screen_duration_ns": float(record["summary"]["total_duration_ns"]),
        "screen_gate_depth": int(record["summary"]["gate_depth"]),
        "replay_n12_fidelity": float(replay.get("fidelity", np.nan)),
        "replay_n12_leakage": float(replay.get("leakage_worst", np.nan)),
        "replay_n12_unitarity_error": float(replay.get("unitarity_error", np.nan)),
    }


def write_design_space_tables(compressed_records: list[dict[str, Any]], final_records: list[dict[str, Any]]) -> dict[str, str]:
    design_path = c.DATA_DIR / "corrected_design_space_best.csv"
    with design_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = list(slim_design_space_row(compressed_records[0]).keys()) if compressed_records else []
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        if fieldnames:
            writer.writeheader()
            for record in sorted(compressed_records, key=replay_sort_key, reverse=True):
                writer.writerow(slim_design_space_row(record))

    finalist_path = c.DATA_DIR / "corrected_target_99_finalists.csv"
    with finalist_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "family_key",
            "family_label",
            "blocks",
            "n_active",
            "order_label",
            "levels",
            "n12_optimized_fidelity",
            "n12_physical_fidelity",
            "n12_leakage_worst",
            "status",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for record in sorted(final_records, key=lambda row: physical_eval(row, DEFAULT_FINAL_TRUNC)["fidelity"], reverse=True):
            writer.writerow(
                {
                    "family_key": record["family_key"],
                    "family_label": record["family_label"],
                    "blocks": record["blocks"],
                    "n_active": record["max_tones"],
                    "order_label": record["order_label"],
                    "levels": "-".join(str(level) for level in record.get("levels", [])),
                    "n12_optimized_fidelity": float(record["fidelity"]),
                    "n12_physical_fidelity": float(physical_eval(record, DEFAULT_FINAL_TRUNC)["fidelity"]),
                    "n12_leakage_worst": float(physical_eval(record, DEFAULT_FINAL_TRUNC)["leakage_worst"]),
                    "status": record["physical"]["status"],
                }
            )
    return {"design_space": str(design_path), "finalists": str(finalist_path)}


def build_factor_summary(records: list[dict[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for field in ("max_tones", "blocks", "order_label"):
        grouped: dict[str, list[dict[str, Any]]] = {}
        for record in records:
            grouped.setdefault(str(record[field]), []).append(record)
        rows = []
        best_values: list[float] = []
        for value, group in sorted(grouped.items(), key=lambda item: item[0]):
            fidelities = [float(item["replay"][str(DEFAULT_FINAL_TRUNC)]["fidelity"]) for item in group]
            best = best_record(group, key=replay_sort_key)
            best_fidelity = max(fidelities)
            best_values.append(best_fidelity)
            rows.append(
                {
                    "value": value,
                    "best_fidelity": float(best_fidelity),
                    "mean_fidelity": float(statistics.fmean(fidelities)),
                    "median_fidelity": _median(fidelities),
                    "best_case_id": str(best["case_id"]),
                    "best_levels": list(best.get("levels", [])),
                    "best_order_label": str(best.get("order_label", "")),
                    "best_blocks": int(best.get("blocks", 0)),
                    "best_n_active": int(best.get("max_tones", 0)),
                }
            )
        effect_range = max(best_values) - min(best_values) if best_values else 0.0
        result[field] = {"rows": rows, "effect_range": float(effect_range)}
    return result


def build_heatmap_rows(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[int, int], list[dict[str, Any]]] = {}
    for record in records:
        grouped.setdefault((int(record["blocks"]), int(record["max_tones"])), []).append(record)
    rows = []
    for (blocks, n_active), group in sorted(grouped.items()):
        best = best_record(group, key=replay_sort_key)
        rows.append(
            {
                "blocks": int(blocks),
                "n_active": int(n_active),
                "best_replay_fidelity": float(best["replay"][str(DEFAULT_FINAL_TRUNC)]["fidelity"]),
                "best_order_label": str(best["order_label"]),
                "best_case_id": str(best["case_id"]),
            }
        )
    return rows


def build_family_analysis(records: list[dict[str, Any]]) -> dict[str, Any]:
    grouped: dict[str, list[dict[str, Any]]] = {}
    for record in records:
        grouped.setdefault(str(record["family_key"]), []).append(record)
    analysis: dict[str, Any] = {}
    for family_key, family_records in grouped.items():
        factor_summary = build_factor_summary(family_records)
        dominant_factor = max(
            factor_summary.items(),
            key=lambda item: float(item[1]["effect_range"]),
        )[0]
        ordering_rows = factor_summary["order_label"]["rows"]
        consistent_ordering = max(ordering_rows, key=lambda row: (float(row["median_fidelity"]), float(row["best_fidelity"])))
        best_replay = best_record(family_records, key=replay_sort_key)
        analysis[family_key] = {
            "family_label": str(best_replay["family_label"]),
            "dominant_factor": str(dominant_factor),
            "factor_summary": factor_summary,
            "consistent_best_ordering": consistent_ordering,
            "heatmap_rows": build_heatmap_rows(family_records),
            "best_replay_candidate": record_without_fit(best_replay),
        }
    return analysis


def build_target_99_analysis(final_records: list[dict[str, Any]]) -> dict[str, Any]:
    grouped: dict[str, list[dict[str, Any]]] = {}
    for record in final_records:
        grouped.setdefault(str(record["family_key"]), []).append(record)
    analysis: dict[str, Any] = {}
    for family_key, family_records in grouped.items():
        best_any = best_record(family_records, key=lambda row: physical_eval(row, DEFAULT_FINAL_TRUNC)["fidelity"])
        realistic_pool = [row for row in family_records if row["physical"]["status"] != "discarded"]
        best_realistic = None if not realistic_pool else best_record(realistic_pool, key=lambda row: physical_eval(row, DEFAULT_FINAL_TRUNC)["fidelity"])
        hits = [
            row for row in realistic_pool if float(physical_eval(row, DEFAULT_FINAL_TRUNC)["fidelity"]) >= TARGET_FIDELITY
        ]
        best_hit = None if not hits else best_record(hits, key=lambda row: physical_eval(row, DEFAULT_FINAL_TRUNC)["fidelity"])
        baseline = best_realistic if best_realistic is not None else best_any
        best_fidelity = float(physical_eval(baseline, DEFAULT_FINAL_TRUNC)["fidelity"])
        analysis[family_key] = {
            "family_label": str(best_any["family_label"]),
            "achieved": best_hit is not None,
            "gap_to_target": float(max(0.0, TARGET_FIDELITY - best_fidelity)),
            "best_any": record_without_fit(best_any),
            "best_realistic": None if best_realistic is None else record_without_fit(best_realistic),
            "best_hit": None if best_hit is None else record_without_fit(best_hit),
            "closest_configuration": {
                "blocks": int(baseline["blocks"]),
                "n_active": int(baseline["max_tones"]),
                "order_label": str(baseline["order_label"]),
                "levels": list(baseline.get("levels", [])),
                "n12_fidelity": float(physical_eval(baseline, DEFAULT_FINAL_TRUNC)["fidelity"]),
                "n12_leakage_worst": float(physical_eval(baseline, DEFAULT_FINAL_TRUNC)["leakage_worst"]),
                "status": str(baseline["physical"]["status"]),
            },
        }
    return analysis


def plot_family_tradeoff(summary: dict[str, Any]) -> None:
    best_candidates = summary["best_candidates"]
    target_analysis = summary["target_99_analysis"]
    labels = [best_candidates[key]["family_label"] for key in ("drsqr", "drcpsqr")]
    prelim = [summary["preliminary_best"][key] for key in ("drsqr", "drcpsqr")]
    best_replay = [summary["family_analysis"][key]["best_replay_candidate"] for key in ("drsqr", "drcpsqr")]
    best_final = [best_candidates[key] for key in ("drsqr", "drcpsqr")]

    x = np.arange(len(labels))
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.1))
    axes[0].bar(x - 0.24, [row["fidelity"] for row in prelim], width=0.24, label=f"screen N={PRELIM_TRUNC}")
    axes[0].bar(
        x,
        [row["replay"][str(DEFAULT_FINAL_TRUNC)]["fidelity"] for row in best_replay],
        width=0.24,
        label=f"replay N={DEFAULT_FINAL_TRUNC}",
    )
    axes[0].bar(
        x + 0.24,
        [physical_eval(row, DEFAULT_FINAL_TRUNC)["fidelity"] for row in best_final],
        width=0.24,
        label=f"optimized N={DEFAULT_FINAL_TRUNC}",
    )
    axes[0].axhline(TARGET_FIDELITY, color="0.25", linestyle="--", linewidth=1.0)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels, rotation=8)
    axes[0].set_ylabel("Fidelity")
    axes[0].set_ylim(0.0, 1.03)
    axes[0].set_title("Best fidelity by search stage")
    axes[0].legend(frameon=False)

    gaps = [float(target_analysis[key]["gap_to_target"]) for key in ("drsqr", "drcpsqr")]
    axes[1].bar(x, gaps, color=["#4477AA", "#EE6677"])
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=8)
    axes[1].set_ylabel("Gap to 99% target")
    axes[1].set_title("Closest realistic miss to 99%")
    for idx, gap in enumerate(gaps):
        axes[1].text(idx, gap + 0.002, f"{gap:.3f}", ha="center", va="bottom", fontsize=8)

    fig.tight_layout()
    fig.savefig(c.FIG_DIR / "corrected_family_tradeoff.pdf", bbox_inches="tight")
    fig.savefig(c.FIG_DIR / "corrected_family_tradeoff.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_design_space_heatmaps(summary: dict[str, Any]) -> None:
    family_analysis = summary["family_analysis"]
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.3), sharey=True)
    families = ("drsqr", "drcpsqr")
    for axis, family_key in zip(axes, families):
        rows = family_analysis[family_key]["heatmap_rows"]
        block_values = sorted({int(row["blocks"]) for row in rows})
        n_active_values = sorted({int(row["n_active"]) for row in rows})
        matrix = np.full((len(block_values), len(n_active_values)), np.nan, dtype=float)
        annotations: dict[tuple[int, int], tuple[str, float]] = {}
        for row in rows:
            block_idx = block_values.index(int(row["blocks"]))
            tone_idx = n_active_values.index(int(row["n_active"]))
            matrix[block_idx, tone_idx] = float(row["best_replay_fidelity"])
            annotations[(block_idx, tone_idx)] = (str(row["best_order_label"]), float(row["best_replay_fidelity"]))
        image = axis.imshow(matrix, vmin=0.85, vmax=1.0, cmap="viridis", aspect="auto", origin="lower")
        axis.set_xticks(np.arange(len(n_active_values)))
        axis.set_xticklabels(n_active_values)
        axis.set_yticks(np.arange(len(block_values)))
        axis.set_yticklabels(block_values)
        axis.set_xlabel("n_active")
        axis.set_title(f"{family_analysis[family_key]['family_label']} replay landscape")
        for (block_idx, tone_idx), (order_label, fidelity) in annotations.items():
            axis.text(tone_idx, block_idx, f"{order_label}\n{fidelity:.3f}", ha="center", va="center", fontsize=7, color="white")
    axes[0].set_ylabel("Blocks")
    colorbar = fig.colorbar(image, ax=axes.ravel().tolist(), shrink=0.9)
    colorbar.set_label(f"Best replay fidelity at N={DEFAULT_FINAL_TRUNC}")
    fig.tight_layout()
    fig.savefig(c.FIG_DIR / "corrected_design_space_heatmaps.pdf", bbox_inches="tight")
    fig.savefig(c.FIG_DIR / "corrected_design_space_heatmaps.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_ordering_tradeoff(summary: dict[str, Any]) -> None:
    family_analysis = summary["family_analysis"]
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.2), sharey=True)
    for axis, family_key in zip(axes, ("drsqr", "drcpsqr")):
        rows = sorted(
            family_analysis[family_key]["factor_summary"]["order_label"]["rows"],
            key=lambda row: float(row["median_fidelity"]),
            reverse=True,
        )
        x = np.arange(len(rows))
        axis.bar(x - 0.18, [float(row["best_fidelity"]) for row in rows], width=0.36, label="best")
        axis.bar(x + 0.18, [float(row["median_fidelity"]) for row in rows], width=0.36, label="median")
        axis.set_xticks(x)
        axis.set_xticklabels([row["value"] for row in rows], rotation=20)
        axis.set_ylim(0.85, 1.01)
        axis.set_title(f"{family_analysis[family_key]['family_label']} ordering comparison")
        axis.axhline(TARGET_FIDELITY, color="0.25", linestyle="--", linewidth=1.0)
        axis.legend(frameon=False)
    axes[0].set_ylabel(f"Replay fidelity at N={DEFAULT_FINAL_TRUNC}")
    fig.tight_layout()
    fig.savefig(c.FIG_DIR / "corrected_ordering_tradeoff.pdf", bbox_inches="tight")
    fig.savefig(c.FIG_DIR / "corrected_ordering_tradeoff.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_active_subspace(summary: dict[str, Any]) -> None:
    rows = []
    for family_key in ("drsqr", "drcpsqr"):
        candidate = summary["best_candidates"][family_key]
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
    sqr = summary.get("best_candidates", {}).get("drsqr") or summary.get("retained_candidates", {}).get("drsqr")
    cpsqr = summary.get("best_candidates", {}).get("drcpsqr") or summary.get("retained_candidates", {}).get("drcpsqr")
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
                ax.set_ylabel(f"{input_label}\nIm(alpha)")
            if row_index == len(labels) - 1:
                ax.set_xlabel("Re(alpha)")

    for ax in axes.ravel():
        ax.label_outer()

    cax = fig.add_subplot(grid[:, -1])
    colorbar = fig.colorbar(mesh, cax=cax)
    colorbar.set_label("W(alpha)")
    cax.yaxis.set_ticks_position("right")
    cax.yaxis.set_label_position("right")

    fig.savefig(c.FIG_DIR / "appendix_wigner.pdf", bbox_inches="tight")
    fig.savefig(c.FIG_DIR / "appendix_wigner.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_short_change_summary(summary: dict[str, Any]) -> None:
    sqr = summary["target_99_analysis"]["drsqr"]
    cpsqr = summary["target_99_analysis"]["drcpsqr"]
    lines = [
        "# Comprehensive Structured Design-Space Update",
        "",
        f"- Systematically screened {summary['design_space']['screen_case_count']} structural cases at N_cav = {PRELIM_TRUNC} across blocks, n_active, and gate ordering.",
        f"- The structural screen produced {summary['design_space']['compressed_case_count']} unique (family, blocks, n_active, ordering) configurations, then refined {summary['design_space']['level_refine_case_count']} level-subset variants for the strongest structures.",
        f"- Best realistic D + R + SQR configuration: blocks={sqr['closest_configuration']['blocks']}, n_active={sqr['closest_configuration']['n_active']}, ordering={sqr['closest_configuration']['order_label']}, N_cav=12 fidelity={sqr['closest_configuration']['n12_fidelity']:.4f}.",
        f"- Best realistic D + R + CPSQR configuration: blocks={cpsqr['closest_configuration']['blocks']}, n_active={cpsqr['closest_configuration']['n_active']}, ordering={cpsqr['closest_configuration']['order_label']}, N_cav=12 fidelity={cpsqr['closest_configuration']['n12_fidelity']:.4f}.",
        f"- D + R + SQR reaches the 99% target: {'yes' if sqr['achieved'] else 'no'}.",
        f"- D + R + CPSQR reaches the 99% target: {'yes' if cpsqr['achieved'] else 'no'}.",
        f"- Dominant replay-space driver for D + R + SQR: {summary['family_analysis']['drsqr']['dominant_factor']}.",
        f"- Dominant replay-space driver for D + R + CPSQR: {summary['family_analysis']['drcpsqr']['dominant_factor']}.",
    ]
    (c.DATA_DIR / "what_changed_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_searches() -> dict[str, Any]:
    all_cases = make_sqr_cases() + make_cpsqr_cases()
    screened = screen_design_space(all_cases)
    compressed_records = screened["best_by_config"]
    for record in compressed_records:
        replay_record(record)

    family_groups: dict[str, list[dict[str, Any]]] = {"drsqr": [], "drcpsqr": []}
    for record in compressed_records:
        family_groups[str(record["family_key"])].append(record)

    structural_finalists: list[dict[str, Any]] = []
    for family_key in ("drsqr", "drcpsqr"):
        structural_finalists.extend(select_finalists(family_groups[family_key], count=STRUCTURAL_FINALISTS_PER_FAMILY))

    level_refined = refine_level_subsets(structural_finalists, maxiter=PRELIM_REFINE_MAXITER)
    for record in level_refined:
        replay_record(record)

    level_refined_groups: dict[str, list[dict[str, Any]]] = {"drsqr": [], "drcpsqr": []}
    for record in level_refined:
        level_refined_groups[str(record["family_key"])].append(record)

    physical_finalists: list[dict[str, Any]] = []
    for family_key in ("drsqr", "drcpsqr"):
        physical_finalists.extend(select_finalists(level_refined_groups[family_key], count=PHYSICAL_FINALISTS_PER_FAMILY))

    physical_refined = refine_finalists(physical_finalists, n_cav=DEFAULT_FINAL_TRUNC, phase="physical_refine", maxiter=PHYSICAL_REFINE_MAXITER)

    for record in physical_refined:
        record["physical"] = evaluate_physical_candidate(record)

    best_candidates: dict[str, dict[str, Any]] = {}
    retained_candidates: dict[str, dict[str, Any] | None] = {}
    for family_key in ("drsqr", "drcpsqr"):
        family_records = [row for row in physical_refined if str(row["family_key"]) == family_key]
        best_candidate = best_record(family_records, key=lambda row: physical_eval(row, DEFAULT_FINAL_TRUNC)["fidelity"])
        best_candidates[family_key] = best_candidate
        realistic = [row for row in family_records if row["physical"]["status"] != "discarded"]
        retained_candidates[family_key] = None if not realistic else best_record(realistic, key=lambda row: physical_eval(row, DEFAULT_FINAL_TRUNC)["fidelity"])

    for family_key, stem in (("drsqr", "corrected_best_sqr"), ("drcpsqr", "corrected_best_cpsqr")):
        best_candidates[family_key]["artifact_path"] = save_record_artifact(best_candidates[family_key], stem=stem)

    table_paths = write_design_space_tables(compressed_records, physical_refined)
    family_analysis = build_family_analysis(compressed_records)
    target_analysis = build_target_99_analysis(physical_refined)

    preliminary_best = {
        family_key: record_without_fit(best_record(level_refined_groups[family_key], key=replay_sort_key))
        for family_key in ("drsqr", "drcpsqr")
    }

    summary = {
        "study_name": "cluster_state_holographic_unified",
        "date_created": time.strftime("%Y-%m-%d"),
        "scope": {
            "exclusive_families": ["D + R + SQR", "D + R + CPSQR"],
            "grape_role": "separate reference",
            "default_final_truncation": DEFAULT_FINAL_TRUNC,
            "truncation_checks": list(FINAL_TRUNCATIONS),
            "target_fidelity": TARGET_FIDELITY,
        },
        "search_config": {
            "screen_maxiter": SCREEN_MAXITER,
            "prelim_refine_maxiter": PRELIM_REFINE_MAXITER,
            "physical_refine_maxiter": PHYSICAL_REFINE_MAXITER,
            "prelim_truncation": PRELIM_TRUNC,
            "sqr_block_options": list(SQR_BLOCK_OPTIONS),
            "cpsqr_block_options": list(CPSQR_BLOCK_OPTIONS),
            "n_active_options": list(N_ACTIVE_OPTIONS),
            "sqr_orderings": [compact_order_label(order) for order in SQR_ORDERINGS],
            "cpsqr_orderings": [compact_order_label(order) for order in CPSQR_ORDERINGS],
            "screen_level_strategy": "low_contiguous_only",
            "level_refine_strategy": "all ordered level subsets for selected structural finalists",
            "structural_finalists_per_family": STRUCTURAL_FINALISTS_PER_FAMILY,
            "physical_finalists_per_family": PHYSICAL_FINALISTS_PER_FAMILY,
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
        "design_space": {
            "screen_case_count": int(screened["case_count"]),
            "screen_family_counts": screened["family_counts"],
            "compressed_case_count": len(compressed_records),
            "level_refine_case_count": len(level_refined),
            "compressed_records": [record_without_fit(row) for row in sorted(compressed_records, key=replay_sort_key, reverse=True)],
            "table_paths": table_paths,
        },
        "preliminary_search": {
            family_key: {
                "best_by_configuration": [
                    record_without_fit(row)
                    for row in sorted(family_groups[family_key], key=replay_sort_key, reverse=True)
                ],
                "selected_structural_finalists": [
                    record_without_fit(row)
                    for row in sorted(
                        [candidate for candidate in structural_finalists if str(candidate["family_key"]) == family_key],
                        key=replay_sort_key,
                        reverse=True,
                    )
                ],
                "level_refined_candidates": [
                    record_without_fit(row)
                    for row in sorted(level_refined_groups[family_key], key=replay_sort_key, reverse=True)
                ],
                "selected_physical_finalists": [
                    record_without_fit(row)
                    for row in sorted(
                        [candidate for candidate in physical_finalists if str(candidate["family_key"]) == family_key],
                        key=replay_sort_key,
                        reverse=True,
                    )
                ],
            }
            for family_key in ("drsqr", "drcpsqr")
        },
        "preliminary_best": preliminary_best,
        "best_candidates": {key: record_without_fit(value) for key, value in best_candidates.items()},
        "retained_candidates": {key: None if value is None else record_without_fit(value) for key, value in retained_candidates.items()},
        "discarded_candidates": {
            key: None if retained_candidates[key] is not None else record_without_fit(best_candidates[key])
            for key in ("drsqr", "drcpsqr")
        },
        "family_analysis": family_analysis,
        "target_99_analysis": target_analysis,
        "grape_reference": load_grape_reference(),
    }
    c.save_json(c.DATA_DIR / "corrected_scope_summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the comprehensive structured design-space search for the unified holographic study.")
    parser.add_argument("--skip-search", action="store_true", help="Reuse corrected_scope_summary.json if it already exists.")
    args = parser.parse_args()

    summary_path = c.DATA_DIR / "corrected_scope_summary.json"
    if args.skip_search and summary_path.exists():
        summary = c.load_json(summary_path)
    else:
        summary = run_searches()
    plot_family_tradeoff(summary)
    plot_design_space_heatmaps(summary)
    plot_ordering_tradeoff(summary)
    plot_active_subspace(summary)
    plot_wigner_panels(summary)
    write_short_change_summary(summary)
    print("[done] corrected_scope_summary.json and comprehensive design-space figures written", flush=True)


if __name__ == "__main__":
    main()