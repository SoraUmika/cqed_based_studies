from __future__ import annotations

import argparse
import copy
import csv
import math
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Sequence

import matplotlib.pyplot as plt
import numpy as np
import runtime_compat  # noqa: F401
import qutip as qt
from scipy.optimize import minimize


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import common as c
import run_design_space_study as ds


STYLE_PATH = SCRIPT_DIR.parents[2] / ".github" / "skills" / "publication-figures" / "assets" / "cqed_style.mplstyle"
if STYLE_PATH.exists():
    plt.style.use(str(STYLE_PATH))


BEST_ARTIFACTS = {
    "drsqr": ("corrected_best_sqr.json", "D + R + SQR"),
    "drcpsqr": ("corrected_best_cpsqr.json", "D + R + CPSQR"),
}

FAMILY_ORDERINGS = {
    "drsqr": ds.SQR_ORDERINGS,
    "drcpsqr": ds.CPSQR_ORDERINGS,
}

FAMILY_BLOCK_OPTIONS = {
    "drsqr": ds.SQR_BLOCK_OPTIONS,
    "drcpsqr": ds.CPSQR_BLOCK_OPTIONS,
}

ANGLE_TOL = 1.0e-3
DISPLACEMENT_TOL = 1.0e-3
FIDELITY_THRESHOLD = 0.99
NEAR_BEST_DELTA = 1.0e-3
REPAIR_OBJECTIVE_LEAKAGE_WEIGHT = 0.02
WIGNER_XVEC = np.linspace(-3.0, 3.0, 151)
TAIL_TIME_BOUNDS = (20.0e-9, 2.0e-6)
TAIL_REPAIR_MAXITER = 60


def _load_best_record(family_key: str) -> dict[str, Any]:
    artifact_name, _label = BEST_ARTIFACTS[str(family_key)]
    payload = c.load_json(c.ARTIFACT_DIR / artifact_name)
    return dict(payload["record"])


def _cross_check_values(base: Sequence[int], extra: Sequence[int]) -> tuple[int, ...]:
    values = [int(value) for value in base]
    for value in extra:
        value_int = int(value)
        if value_int not in values:
            values.append(value_int)
    return tuple(values)


def _sequence_from_record(record: dict[str, Any], *, n_cav: int) -> c.GateSequence:
    return c.sequence_from_payload(record["sequence"], n_cav=int(n_cav))


def _warm_start_from_sequence(sequence: c.GateSequence) -> dict[str, Any]:
    return {
        "parameter_vector": sequence.get_parameter_vector().tolist(),
        "time_raw_vector": sequence.get_time_raw_vector(active_only=True).tolist(),
    }


def _copy_timing(source: Any, target: Any) -> Any:
    target.optimize_time = bool(getattr(source, "optimize_time", False))
    target.time_bounds = getattr(source, "time_bounds", None)
    target.duration_ref = float(getattr(source, "duration_ref", getattr(source, "duration", 0.0)))
    target.time_group = getattr(source, "time_group", None)
    target.time_policy_locked = bool(getattr(source, "time_policy_locked", False))
    return target


def _enable_tail_optimization(gate: Any) -> Any:
    gate.optimize_time = True
    gate.time_bounds = TAIL_TIME_BOUNDS
    gate.duration_ref = float(getattr(gate, "duration", 0.0))
    gate.time_group = f"tail:{gate.name}"
    gate.time_policy_locked = False
    return gate


def _wrap_centered(value: float, period: float) -> float:
    return ((float(value) + 0.5 * float(period)) % float(period)) - 0.5 * float(period)


def _canonicalize_xy(theta: float, phi: float) -> tuple[float, float]:
    theta_wrapped = _wrap_centered(theta, 4.0 * np.pi)
    phi_wrapped = _wrap_centered(phi, 2.0 * np.pi)
    if theta_wrapped < 0.0:
        theta_wrapped = -theta_wrapped
        phi_wrapped = _wrap_centered(phi_wrapped + np.pi, 2.0 * np.pi)
    return float(theta_wrapped), float(phi_wrapped)


def _canonicalize_cpsqr_phase(phase: float) -> float:
    return float(_wrap_centered(phase, 4.0 * np.pi))


def _make_displacement(name: str, alpha: complex, template: Any) -> Any:
    gate = c.displacement_gate(name)
    gate.alpha = complex(alpha)
    gate.duration = float(getattr(template, "duration", gate.duration))
    return _copy_timing(template, gate)


def _make_rotation(name: str, theta: float, phi: float, template: Any) -> Any:
    gate = c.rotation_gate(name, phi=float(phi))
    gate.theta = float(theta)
    gate.phi = float(phi)
    gate.duration = float(getattr(template, "duration", gate.duration))
    return _copy_timing(template, gate)


def _make_masked_sqr(name: str, *, levels: Sequence[int], theta: Sequence[float], phi: Sequence[float], duration: float, template: Any) -> Any:
    gate = c.make_masked_sqr_gate(
        name=name,
        levels=tuple(int(level) for level in levels),
        n_cav=int(template.hilbert_dim // 2),
        duration_s=float(duration),
        include_conditional_phase=bool(getattr(template, "metadata", {}).get("include_conditional_phase", False)),
    )
    gate.parameters = {
        "theta": np.asarray(theta, dtype=float),
        "phi": np.asarray(phi, dtype=float),
        "duration": float(duration),
    }
    gate.duration = float(duration)
    return _copy_timing(template, gate)


def _make_masked_cpsqr(name: str, *, levels: Sequence[int], phases: Sequence[float], duration: float, template: Any) -> Any:
    gate = c.make_masked_cpsqr_gate(
        name=name,
        levels=tuple(int(level) for level in levels),
        n_cav=int(template.hilbert_dim // 2),
        duration_s=float(duration),
        include_drift=bool(getattr(template, "metadata", {}).get("include_drift", True)),
    )
    gate.parameters = {
        "phases": np.asarray(phases, dtype=float),
        "duration": float(duration),
    }
    gate.duration = float(duration)
    return _copy_timing(template, gate)


def _canonicalize_gate(gate: Any) -> tuple[Any | None, dict[str, int]]:
    stats = {
        "wrapped": 0,
        "dropped_gates": 0,
        "dropped_tones": 0,
    }
    kind = c.gate_kind(gate)
    if kind == "Displacement":
        alpha = complex(getattr(gate, "alpha", 0.0 + 0.0j))
        if abs(alpha) < DISPLACEMENT_TOL:
            stats["dropped_gates"] += 1
            return None, stats
        return _make_displacement(str(gate.name), alpha, gate), stats

    if kind == "QubitRotation":
        theta, phi = _canonicalize_xy(float(getattr(gate, "theta", 0.0)), float(getattr(gate, "phi", 0.0)))
        if abs(theta) < ANGLE_TOL:
            stats["dropped_gates"] += 1
            return None, stats
        if abs(theta - float(getattr(gate, "theta", 0.0))) > 1.0e-12 or abs(phi - float(getattr(gate, "phi", 0.0))) > 1.0e-12:
            stats["wrapped"] += 1
        return _make_rotation(str(gate.name), theta, phi, gate), stats

    if kind == "MaskedSQR":
        levels = list(getattr(gate, "metadata", {}).get("levels", []))
        theta_raw = np.asarray(gate.parameters.get("theta", []), dtype=float)
        phi_raw = np.asarray(gate.parameters.get("phi", []), dtype=float)
        kept_levels: list[int] = []
        kept_theta: list[float] = []
        kept_phi: list[float] = []
        for level, theta_value, phi_value in zip(levels, theta_raw, phi_raw, strict=True):
            theta_can, phi_can = _canonicalize_xy(float(theta_value), float(phi_value))
            if abs(theta_can) < ANGLE_TOL:
                stats["dropped_tones"] += 1
                continue
            if abs(theta_can - float(theta_value)) > 1.0e-12 or abs(phi_can - float(phi_value)) > 1.0e-12:
                stats["wrapped"] += 1
            kept_levels.append(int(level))
            kept_theta.append(theta_can)
            kept_phi.append(phi_can)
        if not kept_levels:
            stats["dropped_gates"] += 1
            return None, stats
        return _make_masked_sqr(
            str(gate.name),
            levels=kept_levels,
            theta=kept_theta,
            phi=kept_phi,
            duration=float(getattr(gate, "duration", gate.parameters.get("duration", c.SQR_S))),
            template=gate,
        ), stats

    if kind == "MaskedCPSQR":
        levels = list(getattr(gate, "metadata", {}).get("levels", []))
        phases_raw = np.asarray(gate.parameters.get("phases", []), dtype=float)
        kept_levels: list[int] = []
        kept_phases: list[float] = []
        for level, phase_value in zip(levels, phases_raw, strict=True):
            phase_can = _canonicalize_cpsqr_phase(float(phase_value))
            if abs(phase_can) < ANGLE_TOL:
                stats["dropped_tones"] += 1
                continue
            if abs(phase_can - float(phase_value)) > 1.0e-12:
                stats["wrapped"] += 1
            kept_levels.append(int(level))
            kept_phases.append(phase_can)
        if not kept_levels:
            stats["dropped_gates"] += 1
            return None, stats
        return _make_masked_cpsqr(
            str(gate.name),
            levels=kept_levels,
            phases=kept_phases,
            duration=float(getattr(gate, "duration", gate.parameters.get("duration", c.CPSQR_S))),
            template=gate,
        ), stats

    return copy.deepcopy(gate), stats


def _same_levels(lhs: Any, rhs: Any) -> bool:
    return tuple(lhs.metadata.get("levels", [])) == tuple(rhs.metadata.get("levels", []))


def _same_phi(lhs: float, rhs: float) -> bool:
    return abs(_wrap_centered(float(lhs) - float(rhs), 2.0 * np.pi)) <= ANGLE_TOL


def _merge_pair(lhs: Any, rhs: Any) -> tuple[Any | None, bool]:
    lhs_kind = c.gate_kind(lhs)
    rhs_kind = c.gate_kind(rhs)
    if lhs_kind != rhs_kind:
        return None, False
    if lhs_kind == "Displacement":
        merged = _make_displacement(str(lhs.name), complex(lhs.alpha) + complex(rhs.alpha), lhs)
        merged.duration = float(lhs.duration) + float(rhs.duration)
        merged.duration_ref = float(getattr(lhs, "duration_ref", lhs.duration)) + float(getattr(rhs, "duration_ref", rhs.duration))
        return merged, True
    if lhs_kind == "QubitRotation" and _same_phi(lhs.phi, rhs.phi):
        merged = _make_rotation(str(lhs.name), float(lhs.theta) + float(rhs.theta), float(lhs.phi), lhs)
        merged.duration = float(lhs.duration) + float(rhs.duration)
        merged.duration_ref = float(getattr(lhs, "duration_ref", lhs.duration)) + float(getattr(rhs, "duration_ref", rhs.duration))
        return merged, True
    if lhs_kind == "MaskedCPSQR" and _same_levels(lhs, rhs):
        merged = _make_masked_cpsqr(
            str(lhs.name),
            levels=lhs.metadata.get("levels", []),
            phases=np.asarray(lhs.parameters["phases"], dtype=float) + np.asarray(rhs.parameters["phases"], dtype=float),
            duration=float(lhs.duration) + float(rhs.duration),
            template=lhs,
        )
        merged.duration_ref = float(getattr(lhs, "duration_ref", lhs.duration)) + float(getattr(rhs, "duration_ref", rhs.duration))
        return merged, True
    if lhs_kind == "MaskedSQR" and _same_levels(lhs, rhs):
        lhs_phi = np.asarray(lhs.parameters["phi"], dtype=float)
        rhs_phi = np.asarray(rhs.parameters["phi"], dtype=float)
        if lhs_phi.shape == rhs_phi.shape and np.all(np.abs([_wrap_centered(a - b, 2.0 * np.pi) for a, b in zip(lhs_phi, rhs_phi, strict=True)]) <= ANGLE_TOL):
            merged = _make_masked_sqr(
                str(lhs.name),
                levels=lhs.metadata.get("levels", []),
                theta=np.asarray(lhs.parameters["theta"], dtype=float) + np.asarray(rhs.parameters["theta"], dtype=float),
                phi=lhs_phi,
                duration=float(lhs.duration) + float(rhs.duration),
                template=lhs,
            )
            merged.duration_ref = float(getattr(lhs, "duration_ref", lhs.duration)) + float(getattr(rhs, "duration_ref", rhs.duration))
            return merged, True
    return None, False


def _renumber_gates(gates: Sequence[Any]) -> list[Any]:
    counters = {"Displacement": 0, "QubitRotation": 0, "MaskedSQR": 0, "MaskedCPSQR": 0}
    renamed: list[Any] = []
    for index, gate in enumerate(gates):
        kind = c.gate_kind(gate)
        if kind == "Displacement":
            gate.name = f"D{counters[kind]}"
        elif kind == "QubitRotation":
            gate.name = f"R{counters[kind]}"
        elif kind == "MaskedSQR":
            gate.name = f"S{counters[kind]}"
        elif kind == "MaskedCPSQR":
            gate.name = f"CP{counters[kind]}"
        counters[kind] = counters.get(kind, 0) + 1
        if getattr(gate, "optimize_time", False):
            gate.time_group = f"instance:{index}:{gate.name}"
        renamed.append(gate)
    return renamed


def canonicalize_sequence(sequence: c.GateSequence) -> tuple[c.GateSequence, dict[str, Any]]:
    gates: list[Any] = []
    stats = {"wrapped": 0, "dropped_gates": 0, "dropped_tones": 0}
    for gate in sequence.gates:
        canonical_gate, gate_stats = _canonicalize_gate(gate)
        for key, value in gate_stats.items():
            stats[key] = int(stats.get(key, 0)) + int(value)
        if canonical_gate is not None:
            gates.append(canonical_gate)
    canonical = c.GateSequence(gates=_renumber_gates(gates), n_cav=int(sequence.n_cav))
    stats["gate_depth"] = len(canonical.gates)
    return canonical, stats


def compress_sequence(sequence: c.GateSequence) -> tuple[c.GateSequence, dict[str, Any]]:
    canonical, canonical_stats = canonicalize_sequence(sequence)
    merged: list[Any] = []
    merge_count = 0
    for gate in canonical.gates:
        if not merged:
            merged.append(copy.deepcopy(gate))
            continue
        merged_gate, changed = _merge_pair(merged[-1], gate)
        if changed and merged_gate is not None:
            merge_count += 1
            canonical_gate, gate_stats = _canonicalize_gate(merged_gate)
            canonical_stats["wrapped"] = int(canonical_stats.get("wrapped", 0)) + int(gate_stats.get("wrapped", 0))
            canonical_stats["dropped_gates"] = int(canonical_stats.get("dropped_gates", 0)) + int(gate_stats.get("dropped_gates", 0))
            canonical_stats["dropped_tones"] = int(canonical_stats.get("dropped_tones", 0)) + int(gate_stats.get("dropped_tones", 0))
            if canonical_gate is None:
                merged.pop()
            else:
                merged[-1] = canonical_gate
        else:
            merged.append(copy.deepcopy(gate))
    compressed = c.GateSequence(gates=_renumber_gates(merged), n_cav=int(sequence.n_cav))
    canonical_stats["merged_pairs"] = int(merge_count)
    canonical_stats["gate_depth"] = len(compressed.gates)
    return compressed, canonical_stats


def gate_parameter_count(gate: Any, *, n_cav: int) -> int:
    count = len(gate.parameter_names(int(n_cav)))
    if bool(getattr(gate, "optimize_time", False)):
        count += 1
    return int(count)


def sequence_complexity(sequence: c.GateSequence) -> dict[str, Any]:
    parameter_count = 0
    total_nonzero_tones = 0
    max_nonzero_tones = 0
    for gate in sequence.gates:
        parameter_count += gate_parameter_count(gate, n_cav=int(sequence.n_cav))
        kind = c.gate_kind(gate)
        if kind == "MaskedSQR":
            active = int(np.count_nonzero(np.abs(np.asarray(gate.parameters["theta"], dtype=float)) >= ANGLE_TOL))
            total_nonzero_tones += active
            max_nonzero_tones = max(max_nonzero_tones, active)
        elif kind == "MaskedCPSQR":
            active = int(np.count_nonzero(np.abs(np.asarray(gate.parameters["phases"], dtype=float)) >= ANGLE_TOL))
            total_nonzero_tones += active
            max_nonzero_tones = max(max_nonzero_tones, active)
    return {
        "parameter_count": int(parameter_count),
        "gate_depth": int(len(sequence.gates)),
        "total_nonzero_tones": int(total_nonzero_tones),
        "max_nonzero_tones": int(max_nonzero_tones),
    }


def _candidate_key(record: dict[str, Any]) -> tuple[str, int, int, str, tuple[int, ...]]:
    return (
        str(record["family_key"]),
        int(record["blocks"]),
        int(record["max_tones"]),
        str(record["order_label"]),
        tuple(int(level) for level in record.get("levels", [])),
    )


def _base_case_match(case: dict[str, Any], base_record: dict[str, Any]) -> bool:
    return _candidate_key(case) == _candidate_key(base_record)


def _nearest_level_subsets(base_levels: Sequence[int], tone_count: int, *, n_cav: int, limit: int) -> list[tuple[int, ...]]:
    if int(tone_count) < 1 or int(tone_count) > int(n_cav):
        return []
    base_set = set(int(level) for level in base_levels)
    base_mean = float(np.mean(tuple(int(level) for level in base_levels))) if base_levels else 0.0
    ordered = c.ordered_level_subsets(int(n_cav), int(tone_count))
    scored = sorted(
        ordered,
        key=lambda subset: (
            len(base_set.symmetric_difference(set(subset))),
            0 if 0 in subset and 1 in subset else 1,
            abs(float(np.mean(subset)) - base_mean),
            tuple(subset),
        ),
    )
    return [tuple(int(level) for level in subset) for subset in scored[: int(limit)]]


def _make_case(base_record: dict[str, Any], *, blocks: int, max_tones: int, levels: Sequence[int], order_tokens: Sequence[str]) -> dict[str, Any]:
    case = ds.case_from_record(base_record)
    case["blocks"] = int(blocks)
    case["max_tones"] = int(max_tones)
    case["levels"] = tuple(int(level) for level in levels)
    case["order_tokens"] = tuple(str(token) for token in order_tokens)
    case["order_label"] = ds.compact_order_label(order_tokens)
    case["variant_key"] = f"blk{int(blocks)}_a{int(max_tones)}_ord{case['order_label']}"
    case["variant_label"] = f"{int(blocks)} blocks / {case['order_label']} / n_active={int(max_tones)}"
    builder_kwargs = dict(case.get("builder_kwargs", {}))
    builder_kwargs["blocks"] = int(blocks)
    builder_kwargs["order"] = tuple(str(token) for token in order_tokens)
    case["builder_kwargs"] = builder_kwargs
    return case


def _generate_neighborhood_cases(base_record: dict[str, Any], *, n_cav: int, level_limit: int) -> list[dict[str, Any]]:
    family_key = str(base_record["family_key"])
    base_order = tuple(str(token) for token in base_record.get("order_tokens", []))
    base_levels = tuple(int(level) for level in base_record.get("levels", []))
    base_blocks = int(base_record["blocks"])
    base_tones = int(base_record["max_tones"])
    block_options = tuple(int(value) for value in FAMILY_BLOCK_OPTIONS[family_key])

    cases: dict[tuple[str, int, int, str, tuple[int, ...]], dict[str, Any]] = {}

    def add(case: dict[str, Any], tag: str) -> None:
        key = _candidate_key(case)
        if key in cases:
            return
        case_copy = dict(case)
        case_copy["neighborhood_tag"] = str(tag)
        cases[key] = case_copy

    add(_make_case(base_record, blocks=base_blocks, max_tones=base_tones, levels=base_levels, order_tokens=base_order), "exact")

    for order_tokens in FAMILY_ORDERINGS[family_key]:
        add(_make_case(base_record, blocks=base_blocks, max_tones=base_tones, levels=base_levels, order_tokens=order_tokens), "ordering")

    for blocks in block_options:
        if blocks == base_blocks:
            continue
        if abs(int(blocks) - base_blocks) <= 1:
            add(_make_case(base_record, blocks=blocks, max_tones=base_tones, levels=base_levels, order_tokens=base_order), "block_neighbor")

    tone_candidates = sorted({base_tones, max(1, base_tones - 1), min(int(n_cav), base_tones + 1)})
    for tone_count in tone_candidates:
        for levels in _nearest_level_subsets(base_levels, tone_count, n_cav=int(n_cav), limit=int(level_limit)):
            add(_make_case(base_record, blocks=base_blocks, max_tones=tone_count, levels=levels, order_tokens=base_order), "level_neighbor")

    return list(cases.values())


def _replay_n12(record: dict[str, Any], *, n_cav: int) -> dict[str, Any]:
    sequence = c.sequence_from_payload(record["sequence"], n_cav=int(n_cav))
    evaluation = c.evaluate_sequence(sequence, n_cav=int(n_cav))
    record["replay"] = {
        str(int(n_cav)): {
            "fidelity": float(evaluation["fidelity"]),
            "leakage_average": float(evaluation["leakage_average"]),
            "leakage_worst": float(evaluation["leakage_worst"]),
            "unitarity_error": float(evaluation["unitarity_error"]),
        }
    }
    return record


def run_synthesis_job(job: dict[str, Any]) -> dict[str, Any]:
    record_updates = dict(job.get("record_updates", {}))
    record = ds.run_synthesis_trial(
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
    return _replay_n12(record, n_cav=int(job["n_cav"]))


def run_synthesis_jobs(jobs: Sequence[dict[str, Any]], *, progress_label: str) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    total = len(jobs)
    with tempfile.TemporaryDirectory(prefix="closed_followup_jobs_") as temp_dir_name:
        temp_dir = Path(temp_dir_name)
        for index, job in enumerate(jobs, start=1):
            print(f"[{progress_label} {index}/{total}] {ds.case_id(job['case'])}", flush=True)
            job_path = temp_dir / f"job_{index:03d}.json"
            out_path = temp_dir / f"result_{index:03d}.json"
            c.save_json(job_path, job)
            completed = subprocess.run(
                [
                    sys.executable,
                    str(Path(__file__).resolve()),
                    "--job-json",
                    str(job_path),
                    "--job-output",
                    str(out_path),
                ],
                capture_output=True,
                text=True,
                encoding="utf-8",
            )
            if completed.returncode != 0:
                raise RuntimeError(
                    f"{progress_label} failed for {ds.case_id(job['case'])}\n"
                    f"stdout:\n{completed.stdout}\n"
                    f"stderr:\n{completed.stderr}"
                )
            records.append(c.load_json(out_path))
    return records


def _run_search_for_family(
    base_record: dict[str, Any],
    *,
    n_cav: int,
    screen_maxiter: int,
    refine_maxiter: int,
    level_limit: int,
    shortlist_count: int,
) -> dict[str, Any]:
    neighborhood = _generate_neighborhood_cases(base_record, n_cav=int(n_cav), level_limit=int(level_limit))
    screen_jobs: list[dict[str, Any]] = []
    for case in neighborhood:
        screen_jobs.append(
            {
                "case": case,
                "n_cav": int(n_cav),
                "search_phase": "followup_screen",
                "seed": 17,
                "init_guess": "heuristic",
                "maxiter": int(screen_maxiter),
                "multistart": 1,
                "record_updates": {"neighborhood_tag": case.get("neighborhood_tag", "screen")},
            }
        )
    screen_records = run_synthesis_jobs(screen_jobs, progress_label=f"{base_record['family_key']}_followup_screen")
    promoted = ds.select_finalists(screen_records, count=int(shortlist_count))

    refine_jobs: list[dict[str, Any]] = []
    for record in promoted:
        for seed in (17, 42, 73):
            refine_jobs.append(
                {
                    "case": ds.case_from_record(record),
                    "n_cav": int(n_cav),
                    "search_phase": "followup_refine",
                    "seed": int(seed),
                    "init_guess": "heuristic",
                    "maxiter": int(refine_maxiter),
                    "multistart": 1,
                    "record_updates": {"neighborhood_tag": record.get("neighborhood_tag", "screen")},
                }
            )
    refine_records = run_synthesis_jobs(refine_jobs, progress_label=f"{base_record['family_key']}_followup_refine")

    best_by_case: dict[tuple[str, int, int, str, tuple[int, ...]], dict[str, Any]] = {}
    for record in refine_records:
        key = _candidate_key(record)
        current = best_by_case.get(key)
        if current is None or ds.replay_sort_key(record) > ds.replay_sort_key(current):
            best_by_case[key] = record
    finalists = ds.select_finalists(list(best_by_case.values()), count=int(shortlist_count))
    return {
        "screen_records": screen_records,
        "refine_records": refine_records,
        "finalists": finalists,
    }


def _sequence_output_state(sequence: c.GateSequence, *, n_cav: int, input_label: str) -> qt.Qobj:
    state = dict(ds.logical_input_states(int(n_cav)))[str(input_label)]
    for gate in sequence.gates:
        state = ds._flat_operator(gate.ideal_unitary(int(n_cav)), n_cav=int(n_cav)) * state
    return state


def _full_state_fidelity(target_state: qt.Qobj, candidate_state: qt.Qobj) -> float:
    target_vec = np.asarray(target_state.full(), dtype=np.complex128).reshape(-1)
    candidate_vec = np.asarray(candidate_state.full(), dtype=np.complex128).reshape(-1)
    return float(abs(np.vdot(target_vec, candidate_vec)) ** 2)


def _quadrature_stats(rho: qt.Qobj, *, n_cav: int) -> dict[str, float]:
    a = qt.destroy(int(n_cav))
    x_op = (a + a.dag()) / math.sqrt(2.0)
    p_op = (a - a.dag()) / (1j * math.sqrt(2.0))
    n_op = a.dag() * a
    parity_diag = np.array([(-1) ** level for level in range(int(n_cav))], dtype=float)
    parity = qt.Qobj(np.diag(parity_diag), dims=[[int(n_cav)], [int(n_cav)]])
    return {
        "x_mean": float(np.real(qt.expect(x_op, rho))),
        "p_mean": float(np.real(qt.expect(p_op, rho))),
        "n_mean": float(np.real(qt.expect(n_op, rho))),
        "parity": float(np.real(qt.expect(parity, rho))),
        "purity": float(np.real((rho * rho).tr())),
        "entropy_vn": float(qt.entropy_vn(rho)),
    }


def sequence_wigner_diagnostics(payload: Sequence[dict[str, Any]], *, n_cav: int) -> dict[str, Any]:
    sequence = c.sequence_from_payload(payload, n_cav=int(n_cav))
    by_input: dict[str, Any] = {}
    reduced_fidelities: list[float] = []
    full_fidelities: list[float] = []
    wigner_rms_values: list[float] = []
    entropy_values: list[float] = []
    for label in c.LOGICAL_LABELS:
        target_state = ds.target_output_state(label, n_cav=int(n_cav))
        candidate_state = _sequence_output_state(sequence, n_cav=int(n_cav), input_label=label)
        target_rho = ds.reduced_cavity_density(target_state, n_cav=int(n_cav))
        candidate_rho = ds.reduced_cavity_density(candidate_state, n_cav=int(n_cav))
        target_wigner = qt.wigner(target_rho, WIGNER_XVEC, WIGNER_XVEC)
        candidate_wigner = qt.wigner(candidate_rho, WIGNER_XVEC, WIGNER_XVEC)
        delta = candidate_wigner - target_wigner
        target_pop = ds.cavity_populations_from_state(target_rho, n_cav=int(n_cav))
        candidate_pop = ds.cavity_populations_from_state(candidate_rho, n_cav=int(n_cav))
        target_stats = _quadrature_stats(target_rho, n_cav=int(n_cav))
        candidate_stats = _quadrature_stats(candidate_rho, n_cav=int(n_cav))
        reduced_fidelity = float(qt.fidelity(target_rho, candidate_rho))
        full_fidelity = _full_state_fidelity(target_state, candidate_state)
        wigner_rms = float(np.sqrt(np.mean(np.square(delta))))
        by_input[str(label)] = {
            "full_state_fidelity": full_fidelity,
            "reduced_state_fidelity": reduced_fidelity,
            "wigner_rms": wigner_rms,
            "wigner_l1": float(np.mean(np.abs(delta))),
            "fock_l1": float(np.mean(np.abs(candidate_pop - target_pop))),
            "target_fock_populations": target_pop.tolist(),
            "candidate_fock_populations": candidate_pop.tolist(),
            "target_stats": target_stats,
            "candidate_stats": candidate_stats,
        }
        reduced_fidelities.append(reduced_fidelity)
        full_fidelities.append(full_fidelity)
        wigner_rms_values.append(wigner_rms)
        entropy_values.append(candidate_stats["entropy_vn"])
    return {
        "by_input": by_input,
        "mean_full_state_fidelity": float(np.mean(full_fidelities)),
        "mean_reduced_state_fidelity": float(np.mean(reduced_fidelities)),
        "min_reduced_state_fidelity": float(np.min(reduced_fidelities)),
        "mean_wigner_rms": float(np.mean(wigner_rms_values)),
        "max_wigner_rms": float(np.max(wigner_rms_values)),
        "mean_cavity_entropy": float(np.mean(entropy_values)),
        "max_cavity_entropy": float(np.max(entropy_values)),
    }


def support_trajectory_analysis(payload: Sequence[dict[str, Any]], *, n_cav: int, sample_points: int) -> dict[str, Any]:
    sequence = c.sequence_from_payload(payload, n_cav=int(n_cav))
    active = ds.active_subspace_analysis(sequence, n_cav=int(n_cav), sample_points=int(sample_points))
    candidate_levels = [int(level) for level in active["candidate_active_levels"]]
    logical_levels = [0, 1]
    peak_by_level = np.zeros(int(n_cav), dtype=float)
    per_input: dict[str, Any] = {}
    for label, state in ds.logical_input_states(int(n_cav)):
        running = state
        trace_rows: list[dict[str, Any]] = []
        max_outside_logical = 0.0
        max_outside_active = 0.0
        max_boundary_population = 0.0
        for gate_index, gate in enumerate(sequence.gates):
            for fraction in np.linspace(1.0 / int(sample_points), 1.0, int(sample_points)):
                partial = ds._partial_gate_unitary(gate, n_cav=int(n_cav), fraction=float(fraction))
                evolved = partial * running
                populations = ds.cavity_populations_from_state(evolved, n_cav=int(n_cav))
                peak_by_level = np.maximum(peak_by_level, populations)
                outside_logical = 1.0 - float(np.sum(populations[logical_levels]))
                outside_active = 1.0 - float(np.sum(populations[candidate_levels])) if candidate_levels else 1.0
                boundary_population = float(populations[-1]) if populations.size else 0.0
                max_outside_logical = max(max_outside_logical, outside_logical)
                max_outside_active = max(max_outside_active, outside_active)
                max_boundary_population = max(max_boundary_population, boundary_population)
                trace_rows.append(
                    {
                        "gate_index": int(gate_index),
                        "gate_name": str(gate.name),
                        "fraction": float(fraction),
                        "outside_logical": float(outside_logical),
                        "outside_active": float(outside_active),
                        "boundary_population": float(boundary_population),
                        "mean_level": float(np.dot(np.arange(int(n_cav), dtype=float), populations)),
                        "peak_level": int(np.argmax(populations)),
                    }
                )
            running = ds._flat_operator(gate.ideal_unitary(int(n_cav)), n_cav=int(n_cav)) * running
        per_input[str(label)] = {
            "max_outside_logical": float(max_outside_logical),
            "max_outside_active": float(max_outside_active),
            "max_boundary_population": float(max_boundary_population),
            "trace": trace_rows,
        }
    return {
        "active_summary": active,
        "candidate_active_width": int(len(candidate_levels)),
        "candidate_peak_population_by_level": peak_by_level.tolist(),
        "per_input": per_input,
    }


def evaluate_variant(payload: Sequence[dict[str, Any]], *, n_cav_values: Sequence[int], support_sample_points: int, include_wigner: bool = False) -> dict[str, Any]:
    by_n: dict[str, Any] = {}
    for n_cav in n_cav_values:
        sequence = c.sequence_from_payload(payload, n_cav=int(n_cav))
        evaluation = c.evaluate_sequence(sequence, n_cav=int(n_cav))
        active = ds.active_subspace_analysis(sequence, n_cav=int(n_cav), sample_points=int(support_sample_points))
        by_n[str(int(n_cav))] = {
            "fidelity": float(evaluation["fidelity"]),
            "leakage_average": float(evaluation["leakage_average"]),
            "leakage_worst": float(evaluation["leakage_worst"]),
            "unitarity_error": float(evaluation["unitarity_error"]),
            "active_subspace": active,
        }
    payload_summary = c.sequence_summary(c.sequence_from_payload(payload, n_cav=int(n_cav_values[0])))
    result = {
        "summary": payload_summary,
        "complexity": sequence_complexity(c.sequence_from_payload(payload, n_cav=int(n_cav_values[0]))),
        "by_n_cav": by_n,
    }
    if include_wigner:
        default_n = int(n_cav_values[0]) if 12 not in n_cav_values else 12
        result["wigner"] = sequence_wigner_diagnostics(payload, n_cav=int(default_n))
        result["support_trajectory"] = support_trajectory_analysis(payload, n_cav=int(default_n), sample_points=int(support_sample_points))
    return result


def _record_rank(record: dict[str, Any], *, n_cav: int) -> tuple[float, float, float, float, float]:
    evaluation = record["evaluated"]["by_n_cav"][str(int(n_cav))]
    complexity = record["evaluated"]["complexity"]
    return (
        float(evaluation["fidelity"]),
        -float(evaluation["leakage_worst"]),
        -float(record["evaluated"]["summary"]["total_duration_ns"]),
        -float(complexity["parameter_count"]),
        -float(complexity["total_nonzero_tones"]),
    )


def _canonical_complexity_for_payload(payload: Sequence[dict[str, Any]], *, n_cav: int) -> dict[str, Any]:
    sequence = c.sequence_from_payload(payload, n_cav=int(n_cav))
    canonical, _stats = canonicalize_sequence(sequence)
    return sequence_complexity(canonical)


def _slim_search_row(record: dict[str, Any], *, n_cav: int) -> dict[str, Any]:
    evaluation = record["evaluated"]["by_n_cav"][str(int(n_cav))]
    complexity = record["evaluated"]["complexity"]
    canonical_complexity = record["canonical_complexity"]
    return {
        "family_key": str(record["family_key"]),
        "family_label": str(record["family_label"]),
        "case_id": str(record["case_id"]),
        "neighborhood_tag": str(record.get("neighborhood_tag", "")),
        "blocks": int(record["blocks"]),
        "max_tones": int(record["max_tones"]),
        "order_label": str(record["order_label"]),
        "levels": "-".join(str(level) for level in record.get("levels", [])),
        "fidelity_n12": float(evaluation["fidelity"]),
        "leakage_worst_n12": float(evaluation["leakage_worst"]),
        "duration_ns": float(record["evaluated"]["summary"]["total_duration_ns"]),
        "support_width_levels": int(len(evaluation["active_subspace"]["candidate_active_levels"])),
        "canonical_parameter_count": int(canonical_complexity["parameter_count"]),
        "canonical_total_nonzero_tones": int(canonical_complexity["total_nonzero_tones"]),
        "gate_depth": int(complexity["gate_depth"]),
    }


def _dominates(lhs: dict[str, Any], rhs: dict[str, Any]) -> bool:
    comparisons = [
        lhs["fidelity_n12"] >= rhs["fidelity_n12"],
        lhs["leakage_worst_n12"] <= rhs["leakage_worst_n12"],
        lhs["duration_ns"] <= rhs["duration_ns"],
        lhs["blocks"] <= rhs["blocks"],
        lhs["max_tones"] <= rhs["max_tones"],
        lhs["support_width_levels"] <= rhs["support_width_levels"],
        lhs["canonical_parameter_count"] <= rhs["canonical_parameter_count"],
    ]
    strict = [
        lhs["fidelity_n12"] > rhs["fidelity_n12"],
        lhs["leakage_worst_n12"] < rhs["leakage_worst_n12"],
        lhs["duration_ns"] < rhs["duration_ns"],
        lhs["blocks"] < rhs["blocks"],
        lhs["max_tones"] < rhs["max_tones"],
        lhs["support_width_levels"] < rhs["support_width_levels"],
        lhs["canonical_parameter_count"] < rhs["canonical_parameter_count"],
    ]
    return all(comparisons) and any(strict)


def _pareto_rows(rows: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    frontier: list[dict[str, Any]] = []
    for row in rows:
        if any(_dominates(other, row) for other in rows if other is not row):
            continue
        frontier.append(dict(row))
    return sorted(frontier, key=lambda row: (-row["fidelity_n12"], row["duration_ns"], row["canonical_parameter_count"]))


def _write_csv(path: Path, rows: Sequence[dict[str, Any]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames))
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key) for key in fieldnames})


def _family_color(family_key: str) -> str:
    return {"drsqr": "#1b6ca8", "drcpsqr": "#c84c09"}.get(str(family_key), "#444444")


def plot_pareto_frontier(rows: Sequence[dict[str, Any]], *, grape_reference: dict[str, Any] | None) -> str:
    frontier_keys = {row["case_id"] for row in _pareto_rows(rows)}
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    for family_key in BEST_ARTIFACTS:
        family_rows = [row for row in rows if row["family_key"] == family_key]
        if not family_rows:
            continue
        ax.scatter(
            [row["duration_ns"] for row in family_rows],
            [row["fidelity_n12"] for row in family_rows],
            s=[40 + 2.0 * row["canonical_parameter_count"] for row in family_rows],
            color=_family_color(family_key),
            alpha=0.25,
            edgecolor="none",
            label=BEST_ARTIFACTS[family_key][1],
        )
        family_frontier = [row for row in family_rows if row["case_id"] in frontier_keys]
        ax.scatter(
            [row["duration_ns"] for row in family_frontier],
            [row["fidelity_n12"] for row in family_frontier],
            s=[65 + 2.5 * row["canonical_parameter_count"] for row in family_frontier],
            color=_family_color(family_key),
            edgecolor="black",
            linewidth=0.7,
        )
    if grape_reference is not None:
        ax.scatter(
            [float(grape_reference["duration_ns"])],
            [float(grape_reference["replay_fidelity"])],
            marker="*",
            s=220,
            color="#222222",
            edgecolor="white",
            linewidth=0.7,
            label="GRAPE reference",
        )
    ax.set_xlabel("Total duration (ns)")
    ax.set_ylabel(r"Validated $F_{12}$")
    ax.set_title("Closed-system structured Pareto frontier")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    stem = "closed_system_pareto_frontier"
    fig.savefig(c.FIG_DIR / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(c.FIG_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return stem


def plot_support_diagnostics(family_payloads: dict[str, dict[str, Any]], *, n_cav: int) -> str:
    fig, axes = plt.subplots(2, 2, figsize=(10.2, 6.8), constrained_layout=True)
    for col, family_key in enumerate(("drsqr", "drcpsqr")):
        payload = family_payloads[family_key]["support_trajectory"]
        peak = np.asarray(payload["candidate_peak_population_by_level"], dtype=float)
        axes[0, col].bar(np.arange(peak.size), peak, color=_family_color(family_key), alpha=0.85)
        axes[0, col].set_title(BEST_ARTIFACTS[family_key][1])
        axes[0, col].set_xlabel("Cavity level")
        axes[0, col].set_ylabel("Peak population")
        for label, per_input in payload["per_input"].items():
            trace = per_input["trace"]
            xs = np.arange(len(trace), dtype=float)
            ys = [float(row["outside_active"]) for row in trace]
            axes[1, col].plot(xs, ys, linewidth=1.1, label=label)
        axes[1, col].set_xlabel("Snapshot index")
        axes[1, col].set_ylabel("Population outside active support")
        axes[1, col].set_title(f"{BEST_ARTIFACTS[family_key][1]} at N_cav={int(n_cav)}")
        axes[1, col].legend(frameon=False, fontsize=8)
    stem = "closed_system_support_diagnostics"
    fig.savefig(c.FIG_DIR / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(c.FIG_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return stem


def plot_repair_tests(repair_summary: dict[str, Any]) -> str:
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2), constrained_layout=True)
    for index, family_key in enumerate(("drsqr", "drcpsqr")):
        tests = repair_summary[family_key]["tests"]
        labels = [row["label"] for row in tests]
        delta_reduced = [row["mean_reduced_state_fidelity"] - repair_summary[family_key]["baseline"]["mean_reduced_state_fidelity"] for row in tests]
        delta_wigner = [repair_summary[family_key]["baseline"]["mean_wigner_rms"] - row["mean_wigner_rms"] for row in tests]
        xs = np.arange(len(labels), dtype=float)
        axes[index].bar(xs - 0.18, delta_reduced, width=0.36, color=_family_color(family_key), alpha=0.9, label=r"$\Delta \bar F_{\rho_c}$")
        axes[index].bar(xs + 0.18, delta_wigner, width=0.36, color="#444444", alpha=0.7, label=r"$-\Delta \overline{\mathrm{RMS}}_W$")
        axes[index].axhline(0.0, color="#333333", linewidth=0.8)
        axes[index].set_xticks(xs)
        axes[index].set_xticklabels(labels, rotation=15)
        axes[index].set_title(BEST_ARTIFACTS[family_key][1])
        axes[index].legend(frameon=False, fontsize=8)
    stem = "closed_system_repair_tests"
    fig.savefig(c.FIG_DIR / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(c.FIG_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    return stem


def _tail_gate_bounds(gate: Any, *, n_cav: int) -> list[tuple[float, float]]:
    kind = c.gate_kind(gate)
    if kind == "Displacement":
        return [(-1.5, 1.5), (-1.5, 1.5)]
    if kind == "QubitRotation":
        return [(-2.0 * np.pi, 2.0 * np.pi), (-np.pi, np.pi)]
    if kind == "MaskedCPSQR":
        return [(-2.0 * np.pi, 2.0 * np.pi)] * len(gate.metadata.get("levels", []))
    raise ValueError(f"Unsupported tail-gate kind '{kind}'.")


def _tail_gate_initial(gate: Any, *, n_cav: int) -> np.ndarray:
    kind = c.gate_kind(gate)
    if kind in {"Displacement", "QubitRotation"}:
        return np.asarray(gate.get_parameters(int(n_cav)), dtype=float)
    if kind == "MaskedCPSQR":
        return np.asarray(gate.parameters.get("phases", []), dtype=float)
    raise ValueError(f"Unsupported tail-gate kind '{kind}'.")


def _set_tail_gate_parameters(gate: Any, params: Sequence[float], *, n_cav: int) -> None:
    kind = c.gate_kind(gate)
    values = np.asarray(params, dtype=float)
    if kind in {"Displacement", "QubitRotation"}:
        gate.set_parameters(values, int(n_cav))
        return
    if kind == "MaskedCPSQR":
        gate.parameters = {
            "phases": np.asarray(values, dtype=float),
            "duration": float(gate.duration),
        }
        return
    raise ValueError(f"Unsupported tail-gate kind '{kind}'.")


def _optimize_tail_gate(base_payload: Sequence[dict[str, Any]], tail_gate: Any, *, n_cav: int) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    base_sequence = c.sequence_from_payload(base_payload, n_cav=int(n_cav))
    initial = _tail_gate_initial(tail_gate, n_cav=int(n_cav))
    bounds = _tail_gate_bounds(tail_gate, n_cav=int(n_cav))

    def objective(params: np.ndarray) -> float:
        gate = copy.deepcopy(tail_gate)
        _set_tail_gate_parameters(gate, params, n_cav=int(n_cav))
        trial = c.GateSequence(gates=copy.deepcopy(base_sequence.gates) + [gate], n_cav=int(n_cav))
        evaluation = c.evaluate_sequence(trial, n_cav=int(n_cav))
        return 1.0 - float(evaluation["fidelity"]) + REPAIR_OBJECTIVE_LEAKAGE_WEIGHT * float(evaluation["leakage_worst"])

    result = minimize(objective, initial, method="L-BFGS-B", bounds=bounds, options={"maxiter": TAIL_REPAIR_MAXITER})
    optimized_gate = copy.deepcopy(tail_gate)
    _set_tail_gate_parameters(optimized_gate, result.x, n_cav=int(n_cav))
    payload = c.GateSequence(gates=copy.deepcopy(base_sequence.gates) + [optimized_gate], n_cav=int(n_cav)).serialize()
    return payload, {"success": bool(result.success), "message": str(result.message), "objective": float(result.fun)}


def run_repair_tests(family_key: str, base_payload: Sequence[dict[str, Any]], *, n_cav: int, levels: Sequence[int]) -> dict[str, Any]:
    baseline = sequence_wigner_diagnostics(base_payload, n_cav=int(n_cav))
    tests: list[dict[str, Any]] = []
    probes = [
        ("disp", _enable_tail_optimization(c.displacement_gate("D_tail"))),
        ("rot", _enable_tail_optimization(c.rotation_gate("R_tail", phi=0.0))),
        (
            "diag",
            _enable_tail_optimization(
                c.make_masked_cpsqr_gate(name="CP_tail", levels=tuple(int(level) for level in levels), n_cav=int(n_cav), include_drift=True)
            ),
        ),
    ]
    for label, gate in probes:
        payload, optimization = _optimize_tail_gate(base_payload, gate, n_cav=int(n_cav))
        diagnostics = sequence_wigner_diagnostics(payload, n_cav=int(n_cav))
        tests.append(
            {
                "label": str(label),
                "optimization": optimization,
                "mean_full_state_fidelity": float(diagnostics["mean_full_state_fidelity"]),
                "mean_reduced_state_fidelity": float(diagnostics["mean_reduced_state_fidelity"]),
                "mean_wigner_rms": float(diagnostics["mean_wigner_rms"]),
                "payload": payload,
            }
        )
    best_test = min(tests, key=lambda row: (row["mean_wigner_rms"], -row["mean_reduced_state_fidelity"])) if tests else None
    mechanism = "no_simple_cleanup"
    if best_test is not None and best_test["mean_wigner_rms"] < baseline["mean_wigner_rms"]:
        mechanism = {
            "disp": "residual cavity displacement",
            "rot": "residual qubit-frame mismatch",
            "diag": "residual conditional phase / spectator-phase error",
        }.get(best_test["label"], "no_simple_cleanup")
    return {
        "baseline": {
            "mean_full_state_fidelity": float(baseline["mean_full_state_fidelity"]),
            "mean_reduced_state_fidelity": float(baseline["mean_reduced_state_fidelity"]),
            "mean_wigner_rms": float(baseline["mean_wigner_rms"]),
        },
        "tests": tests,
        "likely_mechanism": mechanism,
        "best_test": None if best_test is None else best_test["label"],
    }


def load_grape_reference() -> dict[str, Any] | None:
    path = c.DATA_DIR / "grape_frontier_extension.json"
    if not path.exists():
        return None
    payload = c.load_json(path)
    rows = [dict(row) for row in payload.values() if str(row.get("status", "")) == "complete"]
    if not rows:
        return None
    best = max(rows, key=lambda row: (float(row.get("best_replay_fidelity", np.nan)), -float(row.get("duration_ns", np.nan))))
    return {
        "duration_ns": int(best["duration_ns"]),
        "replay_fidelity": float(best["best_replay_fidelity"]),
        "replay_leakage_worst": float(best["best_replay_leakage_worst"]),
        "artifact_path": str(best.get("artifact_path", "")),
    }


def _save_artifact(path: Path, payload: dict[str, Any]) -> str:
    c.save_json(path, payload)
    return str(path)


def _variant_artifact_payload(
    *,
    family_key: str,
    family_label: str,
    variant_name: str,
    sequence_payload: Sequence[dict[str, Any]],
    evaluation: dict[str, Any],
    notes: dict[str, Any],
) -> dict[str, Any]:
    return {
        "study_name": "cluster_state_holographic_unified",
        "date_created": time.strftime("%Y-%m-%d"),
        "description": f"{family_label} {variant_name} closed-system follow-up artifact.",
        "family_key": str(family_key),
        "family_label": str(family_label),
        "variant_name": str(variant_name),
        "load_instructions": "import common as c; payload = c.load_json(path); seq = c.sequence_from_payload(payload['sequence_payload'], n_cav=12)",
        "sequence_payload": sequence_payload,
        "evaluation": evaluation,
        "notes": notes,
    }


def _family_minimality_summary(rows: Sequence[dict[str, Any]], *, family_key: str) -> dict[str, Any]:
    family_rows = [row for row in rows if row["family_key"] == family_key]
    best = max(family_rows, key=lambda row: (row["fidelity_n12"], -row["duration_ns"]))
    above_threshold = [row for row in family_rows if row["fidelity_n12"] >= FIDELITY_THRESHOLD]
    near_best = [row for row in family_rows if row["fidelity_n12"] >= best["fidelity_n12"] - NEAR_BEST_DELTA]
    return {
        "best_case_id": best["case_id"],
        "best_fidelity_n12": float(best["fidelity_n12"]),
        "minimum_blocks_above_threshold": None if not above_threshold else int(min(row["blocks"] for row in above_threshold)),
        "minimum_duration_above_threshold_ns": None if not above_threshold else float(min(row["duration_ns"] for row in above_threshold)),
        "minimum_canonical_parameter_count_near_best": int(min(row["canonical_parameter_count"] for row in near_best)),
        "near_best_case_ids": [str(row["case_id"]) for row in near_best],
    }


def _summary_table_rows(rows: Sequence[dict[str, Any]]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for criterion, selector in (
        ("max_fidelity", lambda family_rows: max(family_rows, key=lambda row: (row["fidelity_n12"], -row["duration_ns"]))),
        (
            "min_duration_F>=0.99",
            lambda family_rows: min([row for row in family_rows if row["fidelity_n12"] >= FIDELITY_THRESHOLD], key=lambda row: (row["duration_ns"], -row["fidelity_n12"])),
        ),
        (
            "min_blocks_F>=0.99",
            lambda family_rows: min([row for row in family_rows if row["fidelity_n12"] >= FIDELITY_THRESHOLD], key=lambda row: (row["blocks"], row["duration_ns"])),
        ),
        (
            "min_canonical_complexity_near_best",
            lambda family_rows: min(
                [row for row in family_rows if row["fidelity_n12"] >= max(item["fidelity_n12"] for item in family_rows) - NEAR_BEST_DELTA],
                key=lambda row: (row["canonical_parameter_count"], row["duration_ns"]),
            ),
        ),
    ):
        for family_key in BEST_ARTIFACTS:
            family_rows = [row for row in rows if row["family_key"] == family_key]
            try:
                winner = selector(family_rows)
            except ValueError:
                continue
            output.append(
                {
                    "criterion": criterion,
                    "family_key": family_key,
                    "family_label": BEST_ARTIFACTS[family_key][1],
                    "case_id": winner["case_id"],
                    "fidelity_n12": float(winner["fidelity_n12"]),
                    "leakage_worst_n12": float(winner["leakage_worst_n12"]),
                    "duration_ns": float(winner["duration_ns"]),
                    "blocks": int(winner["blocks"]),
                    "max_tones": int(winner["max_tones"]),
                    "canonical_parameter_count": int(winner["canonical_parameter_count"]),
                }
            )
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the closed-system structured follow-up study.")
    parser.add_argument("--job-json", type=str, default="")
    parser.add_argument("--job-output", type=str, default="")
    parser.add_argument("--n-cav", type=int, default=12)
    parser.add_argument("--cross-checks", type=int, nargs="*", default=(10, 12, 14))
    parser.add_argument("--extra-truncations", type=int, nargs="*", default=(16,))
    parser.add_argument("--screen-maxiter", type=int, default=32)
    parser.add_argument("--refine-maxiter", type=int, default=64)
    parser.add_argument("--level-limit", type=int, default=4)
    parser.add_argument("--shortlist-count", type=int, default=8)
    parser.add_argument("--support-sample-points", type=int, default=12)
    parser.add_argument("--max-workers", type=int, default=1)
    args = parser.parse_args()

    if args.job_json and args.job_output:
        job = c.load_json(Path(args.job_json))
        result = run_synthesis_job(job)
        c.save_json(Path(args.job_output), result)
        return

    n_cav = int(args.n_cav)
    truncations = _cross_check_values(args.cross_checks, args.extra_truncations)
    grape_reference = load_grape_reference()
    ds.MAX_WORKERS = max(1, int(args.max_workers))
    print(f"[closed-followup] using max_workers={ds.MAX_WORKERS}", flush=True)

    family_search: dict[str, Any] = {}
    slim_rows: list[dict[str, Any]] = []
    final_family_payloads: dict[str, Any] = {}
    repair_summary: dict[str, Any] = {}

    for family_key, (_artifact, family_label) in BEST_ARTIFACTS.items():
        print(f"[closed-followup] search start for {family_key}", flush=True)
        base_record = _load_best_record(family_key)
        search = _run_search_for_family(
            base_record,
            n_cav=n_cav,
            screen_maxiter=int(args.screen_maxiter),
            refine_maxiter=int(args.refine_maxiter),
            level_limit=int(args.level_limit),
            shortlist_count=int(args.shortlist_count),
        )
        evaluated_records: list[dict[str, Any]] = []
        for record in search["finalists"]:
            record["evaluated"] = evaluate_variant(
                record["sequence"],
                n_cav_values=truncations,
                support_sample_points=int(args.support_sample_points),
                include_wigner=False,
            )
            record["canonical_complexity"] = _canonical_complexity_for_payload(record["sequence"], n_cav=int(n_cav))
            evaluated_records.append(record)
            slim_rows.append(_slim_search_row(record, n_cav=int(n_cav)))

        raw_best = max(evaluated_records, key=lambda row: _record_rank(row, n_cav=int(n_cav)))
        raw_payload = raw_best["sequence"]
        raw_variant = evaluate_variant(raw_payload, n_cav_values=truncations, support_sample_points=int(args.support_sample_points), include_wigner=True)

        raw_sequence = c.sequence_from_payload(raw_payload, n_cav=int(n_cav))
        canonical_sequence, canonical_stats = canonicalize_sequence(raw_sequence)
        canonical_payload = canonical_sequence.serialize()
        canonical_variant = evaluate_variant(canonical_payload, n_cav_values=truncations, support_sample_points=int(args.support_sample_points), include_wigner=True)

        compressed_sequence, compression_stats = compress_sequence(raw_sequence)
        compressed_warm_start = _warm_start_from_sequence(compressed_sequence)
        compressed_fit = c.fit_sequence(
            compressed_sequence,
            n_cav=int(n_cav),
            seed=17,
            init_guess="heuristic",
            multistart=1,
            maxiter=max(16, int(args.screen_maxiter)),
            warm_start=compressed_warm_start,
        )
        compressed_payload = compressed_fit["sequence_payload"]
        compressed_variant = evaluate_variant(compressed_payload, n_cav_values=truncations, support_sample_points=int(args.support_sample_points), include_wigner=True)

        raw_artifact = _save_artifact(
            c.ARTIFACT_DIR / f"{family_key}_followup_raw_best.json",
            _variant_artifact_payload(
                family_key=family_key,
                family_label=family_label,
                variant_name="raw_best",
                sequence_payload=raw_payload,
                evaluation=raw_variant,
                notes={"source_case_id": raw_best["case_id"], "search_phase": raw_best["search_phase"]},
            ),
        )
        canonical_artifact = _save_artifact(
            c.ARTIFACT_DIR / f"{family_key}_followup_canonicalized.json",
            _variant_artifact_payload(
                family_key=family_key,
                family_label=family_label,
                variant_name="canonicalized",
                sequence_payload=canonical_payload,
                evaluation=canonical_variant,
                notes={"source_case_id": raw_best["case_id"], "canonicalization": canonical_stats},
            ),
        )
        compressed_artifact = _save_artifact(
            c.ARTIFACT_DIR / f"{family_key}_followup_compressed.json",
            _variant_artifact_payload(
                family_key=family_key,
                family_label=family_label,
                variant_name="compressed",
                sequence_payload=compressed_payload,
                evaluation=compressed_variant,
                notes={
                    "source_case_id": raw_best["case_id"],
                    "canonicalization": canonical_stats,
                    "compression": compression_stats,
                    "compressed_fit_fidelity": float(compressed_fit["fidelity"]),
                },
            ),
        )

        repair_summary[family_key] = run_repair_tests(
            family_key,
            raw_payload,
            n_cav=int(n_cav),
            levels=tuple(int(level) for level in raw_best.get("levels", [])),
        )
        family_search[family_key] = {
            "family_label": family_label,
            "saved_best_case_id": str(base_record["case_id"]),
            "saved_best_fidelity_n12": float(base_record["fidelity"]),
            "screen_count": int(len(search["screen_records"])),
            "refine_count": int(len(search["refine_records"])),
            "raw_best_case_id": str(raw_best["case_id"]),
            "raw_best_search_phase": str(raw_best["search_phase"]),
            "raw_best_artifact": raw_artifact,
            "canonical_artifact": canonical_artifact,
            "compressed_artifact": compressed_artifact,
            "minimality": _family_minimality_summary([row for row in slim_rows if row["family_key"] == family_key], family_key=family_key),
            "raw_best": raw_variant,
            "canonicalized": canonical_variant,
            "compressed": compressed_variant,
            "canonicalization": canonical_stats,
            "compression": compression_stats,
        }
        final_family_payloads[family_key] = raw_variant
        print(
            f"[closed-followup] {family_key}: raw F12={raw_variant['by_n_cav'][str(int(n_cav))]['fidelity']:.6f} "
            f"canonical F12={canonical_variant['by_n_cav'][str(int(n_cav))]['fidelity']:.6f} "
            f"compressed F12={compressed_variant['by_n_cav'][str(int(n_cav))]['fidelity']:.6f}",
            flush=True,
        )

    pareto_rows = _pareto_rows(slim_rows)
    pareto_stem = plot_pareto_frontier(slim_rows, grape_reference=grape_reference)
    support_stem = plot_support_diagnostics({family_key: family_search[family_key]["raw_best"] for family_key in BEST_ARTIFACTS}, n_cav=int(n_cav))
    repair_stem = plot_repair_tests(repair_summary)

    candidates_csv = c.DATA_DIR / "closed_system_search_candidates.csv"
    _write_csv(
        candidates_csv,
        slim_rows,
        fieldnames=[
            "family_key",
            "family_label",
            "case_id",
            "neighborhood_tag",
            "blocks",
            "max_tones",
            "order_label",
            "levels",
            "fidelity_n12",
            "leakage_worst_n12",
            "duration_ns",
            "support_width_levels",
            "canonical_parameter_count",
            "canonical_total_nonzero_tones",
            "gate_depth",
        ],
    )
    pareto_csv = c.DATA_DIR / "closed_system_pareto_frontier.csv"
    _write_csv(
        pareto_csv,
        pareto_rows,
        fieldnames=[
            "family_key",
            "family_label",
            "case_id",
            "blocks",
            "max_tones",
            "order_label",
            "levels",
            "fidelity_n12",
            "leakage_worst_n12",
            "duration_ns",
            "support_width_levels",
            "canonical_parameter_count",
            "canonical_total_nonzero_tones",
            "gate_depth",
        ],
    )
    summary_rows = _summary_table_rows(slim_rows)
    summary_csv = c.DATA_DIR / "closed_system_final_summary_table.csv"
    _write_csv(
        summary_csv,
        summary_rows,
        fieldnames=[
            "criterion",
            "family_key",
            "family_label",
            "case_id",
            "fidelity_n12",
            "leakage_worst_n12",
            "duration_ns",
            "blocks",
            "max_tones",
            "canonical_parameter_count",
        ],
    )

    summary_payload = {
        "study_name": "cluster_state_holographic_unified",
        "date_created": time.strftime("%Y-%m-%d"),
        "scope": {
            "families": {family_key: BEST_ARTIFACTS[family_key][1] for family_key in BEST_ARTIFACTS},
            "n_cav_final": int(n_cav),
            "cross_checks": list(truncations),
            "closed_system_only": True,
        },
        "config": {
            "screen_maxiter": int(args.screen_maxiter),
            "refine_maxiter": int(args.refine_maxiter),
            "level_limit": int(args.level_limit),
            "shortlist_count": int(args.shortlist_count),
            "support_sample_points": int(args.support_sample_points),
        },
        "grape_reference": grape_reference,
        "families": family_search,
        "repair_tests": {
            family_key: {
                "baseline": repair_summary[family_key]["baseline"],
                "likely_mechanism": repair_summary[family_key]["likely_mechanism"],
                "best_test": repair_summary[family_key]["best_test"],
                "tests": [
                    {
                        "label": row["label"],
                        "mean_full_state_fidelity": row["mean_full_state_fidelity"],
                        "mean_reduced_state_fidelity": row["mean_reduced_state_fidelity"],
                        "mean_wigner_rms": row["mean_wigner_rms"],
                        "optimization": row["optimization"],
                    }
                    for row in repair_summary[family_key]["tests"]
                ],
            }
            for family_key in repair_summary
        },
        "search_candidates_csv": str(candidates_csv),
        "pareto_frontier_csv": str(pareto_csv),
        "final_summary_table_csv": str(summary_csv),
        "pareto_frontier": pareto_rows,
        "figures": {
            "pareto": pareto_stem,
            "support": support_stem,
            "repair": repair_stem,
        },
    }
    c.save_json(c.DATA_DIR / "closed_system_followup_summary.json", summary_payload)
    print("[closed-followup] wrote data/closed_system_followup_summary.json", flush=True)


if __name__ == "__main__":
    main()