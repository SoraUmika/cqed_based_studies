"""Phase 5: runtime-validated follow-up for native-entangling hybrid synthesis."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
import csv
import hashlib
import math
import sys
from pathlib import Path
from typing import Any, Iterable, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import qutip as qt
from scipy.linalg import polar

import runtime_compat  # noqa: F401
from common import (
    ARTIFACT_DIR,
    CHI,
    DATA_DIR,
    FIG_DIR,
    LOGICAL_BASIS_LABELS,
    NATIVE_WAIT_DURATION_NS,
    OMEGA_C,
    OMEGA_Q,
    STUDY_ROOT,
    CostWeights,
    apply_publication_style,
    build_frame,
    build_model,
    candidate_weighted_cost,
    dump_json,
    embed_logical_state,
    ensure_sim_root_on_path,
    implementation_complexity_score,
    load_json,
    logical_subspace_indices,
    target_unitary_matrix,
)

sys.stdout.reconfigure(encoding="utf-8")

ensure_sim_root_on_path()

from cqed_sim import (  # noqa: E402
    GrapeConfig,
    GrapeSolver,
    ModelControlChannelSpec,
    NoiseSpec,
    PiecewiseConstantTimeGrid,
    UnitaryObjective,
    build_control_problem_from_model,
)
from cqed_sim.observables import bloch_trajectory_from_states  # noqa: E402
from cqed_sim.pulses.envelopes import square_envelope  # noqa: E402
from cqed_sim.pulses.pulse import Pulse  # noqa: E402
from cqed_sim.sequence import SequenceCompiler  # noqa: E402
from cqed_sim.sim import SimulationConfig, prepare_simulation, simulate_batch  # noqa: E402
from cqed_sim.sim.extractors import cavity_wigner, reduced_cavity_state  # noqa: E402
from cqed_sim.sim.noise import pure_dephasing_time_from_t1_t2  # noqa: E402
from cqed_sim.unitary_synthesis import (  # noqa: E402
    GateSequence,
    QubitRotation,
    Subspace,
    simulate_sequence as synth_simulate_sequence,
    subspace_unitary_fidelity,
)
from cqed_sim.unitary_synthesis.waveform_bridge import waveform_primitive_from_gate  # noqa: E402

import phase2_native_block_search as phase2  # noqa: E402


REFERENCE_N_CAV = 12
TRUNCATION_SWEEP = (10, 12, 14)
RUNTIME_N_TR = 3
RUNTIME_DT_S = 2.0e-9
SURROGATE_STEPS = 40
SURROGATE_DT_S = 32.0e-9
SURROGATE_RESTARTS = 3
SURROGATE_MAXITER = 60
SURROGATE_STORAGE_AMP = 2.0e7
SURROGATE_QUBIT_AMP = 8.0e7
WIGNER_POINTS = 51
WIGNER_EXTENT = 2.5

PHASE5_JSON = DATA_DIR / "phase5_runtime_validation.json"
COMPARISON_CSV = DATA_DIR / "phase5_candidate_comparison.csv"
CONVERGENCE_CSV = DATA_DIR / "phase5_convergence_table.csv"
WEIGHT_SENSITIVITY_CSV = DATA_DIR / "phase5_weight_sensitivity.csv"
SYMBOLIC_CSV = DATA_DIR / "phase5_symbolic_metrics.csv"
RUNTIME_CSV = DATA_DIR / "phase5_runtime_metrics.csv"

COMPARISON_PNG = FIG_DIR / "phase5_candidate_comparison.png"
COMPARISON_PDF = FIG_DIR / "phase5_candidate_comparison.pdf"
CONVERGENCE_PNG = FIG_DIR / "phase5_convergence.png"
CONVERGENCE_PDF = FIG_DIR / "phase5_convergence.pdf"
SENSITIVITY_PNG = FIG_DIR / "phase5_weight_sensitivity.png"
SENSITIVITY_PDF = FIG_DIR / "phase5_weight_sensitivity.pdf"


QUBIT_T1_S = 30.0e-6
QUBIT_T2_S = 20.0e-6
CAVITY_T1_S = 250.0e-6
QUBIT_TPHI_S = pure_dephasing_time_from_t1_t2(t1_s=QUBIT_T1_S, t2_s=QUBIT_T2_S)
NOMINAL_NOISE = NoiseSpec(t1=QUBIT_T1_S, tphi=QUBIT_TPHI_S, kappa=1.0 / CAVITY_T1_S)

PROBE_STATES: tuple[tuple[str, np.ndarray], ...] = (
    ("g0", np.array([1.0, 0.0, 0.0, 0.0], dtype=np.complex128)),
    ("g1", np.array([0.0, 1.0, 0.0, 0.0], dtype=np.complex128)),
    ("e0", np.array([0.0, 0.0, 1.0, 0.0], dtype=np.complex128)),
    ("e1", np.array([0.0, 0.0, 0.0, 1.0], dtype=np.complex128)),
    ("qx_plus_0", np.array([1.0, 0.0, 1.0, 0.0], dtype=np.complex128) / np.sqrt(2.0)),
    ("g_c_plus", np.array([1.0, 1.0, 0.0, 0.0], dtype=np.complex128) / np.sqrt(2.0)),
)

BLOCH_PROBE = "qx_plus_0"
WIGNER_PROBE = "g0"
DIAGNOSTIC_RUNTIME_LABELS = ("R2_exact_runtime_to_exact_runtime", "R2_A_runtime_to_A_runtime")


@dataclass(frozen=True)
class SymbolicCandidateDef:
    label: str
    waits: int
    inner_kind: str
    outer_kind: str
    role: str


@dataclass(frozen=True)
class RuntimeCandidateDef:
    label: str
    waits: int
    symbolic_inner_kind: str
    symbolic_outer_kind: str
    runtime_inner_kind: str
    runtime_outer_kind: str
    role: str
    symbolic_reference: str


@dataclass
class ReplaySegment:
    name: str
    kind: str
    pulses: list[Pulse]
    drive_ops: dict[str, Any]
    duration_s: float
    metadata: dict[str, Any]


@dataclass
class SurrogateBlock:
    label: str
    target_kind: str
    target_matrix: np.ndarray
    pulses: list[Pulse]
    drive_ops: dict[str, Any]
    pulse_metadata: dict[str, Any]
    artifact_path: Path
    result_path: Path
    nominal_fidelity: float
    objective_value: float
    total_duration_s: float
    implementation_complexity: float
    restart_summaries: list[dict[str, Any]] = field(default_factory=list)


@dataclass
class CompiledReplay:
    candidate_label: str
    pulses: list[Pulse]
    drive_ops: dict[str, Any]
    segment_names: list[str]
    segment_kinds: list[str]
    segment_end_times_s: list[float]
    segment_metadata: list[dict[str, Any]]
    total_duration_s: float
    compiled: Any
    prefix_cache: dict[int, Any] = field(default_factory=dict)


SYMBOLIC_CANDIDATES: tuple[SymbolicCandidateDef, ...] = (
    SymbolicCandidateDef(
        label="N1_exact_hc_to_exact_hc",
        waits=1,
        inner_kind="exact_hc",
        outer_kind="exact_hc",
        role="one_wait_exact_baseline",
    ),
    SymbolicCandidateDef(
        label="N2_exact_hc_to_exact_hc",
        waits=2,
        inner_kind="exact_hc",
        outer_kind="exact_hc",
        role="two_wait_symbolic_upper_bound",
    ),
    SymbolicCandidateDef(
        label="N2_A_local_to_A_local",
        waits=2,
        inner_kind="A_local",
        outer_kind="A_local",
        role="two_wait_archive_symbolic_candidate",
    ),
    SymbolicCandidateDef(
        label="N2_B_local_to_B_local",
        waits=2,
        inner_kind="B_local",
        outer_kind="B_local",
        role="two_wait_replayable_archive_reference",
    ),
)

RUNTIME_CANDIDATES: tuple[RuntimeCandidateDef, ...] = (
    RuntimeCandidateDef(
        label="R1_exact_runtime_to_exact_runtime",
        waits=1,
        symbolic_inner_kind="exact_hc",
        symbolic_outer_kind="exact_hc",
        runtime_inner_kind="exact_hc_runtime",
        runtime_outer_kind="exact_hc_runtime",
        role="one_wait_replay_surrogate",
        symbolic_reference="N1_exact_hc_to_exact_hc",
    ),
    RuntimeCandidateDef(
        label="R2_exact_runtime_to_exact_runtime",
        waits=2,
        symbolic_inner_kind="exact_hc",
        symbolic_outer_kind="exact_hc",
        runtime_inner_kind="exact_hc_runtime",
        runtime_outer_kind="exact_hc_runtime",
        role="two_wait_replay_exact_surrogate",
        symbolic_reference="N2_exact_hc_to_exact_hc",
    ),
    RuntimeCandidateDef(
        label="R2_A_runtime_to_A_runtime",
        waits=2,
        symbolic_inner_kind="A_local",
        symbolic_outer_kind="A_local",
        runtime_inner_kind="A_local_runtime",
        runtime_outer_kind="A_local_runtime",
        role="two_wait_replay_archive_surrogate",
        symbolic_reference="N2_A_local_to_A_local",
    ),
    RuntimeCandidateDef(
        label="R2_B_local_replay",
        waits=2,
        symbolic_inner_kind="B_local",
        symbolic_outer_kind="B_local",
        runtime_inner_kind="B_local_direct",
        runtime_outer_kind="B_local_direct",
        role="two_wait_direct_replay_archive_candidate",
        symbolic_reference="N2_B_local_to_B_local",
    ),
)


def logical_subspace(n_cav: int, n_tr: int) -> Subspace:
    full_dim = int(n_cav) * int(n_tr)
    return Subspace.custom(full_dim=full_dim, indices=logical_subspace_indices(n_cav), labels=LOGICAL_BASIS_LABELS)


def embed_logical_state_multilevel(logical_vector: Sequence[complex], *, n_cav: int, n_tr: int) -> qt.Qobj:
    vector = np.asarray(logical_vector, dtype=np.complex128).reshape(-1)
    if vector.shape != (4,):
        raise ValueError(f"logical_vector must have shape (4,), got {vector.shape}")
    full = np.zeros(int(n_cav) * int(n_tr), dtype=np.complex128)
    full[0] = vector[0]
    full[1] = vector[1]
    full[int(n_cav)] = vector[2]
    full[int(n_cav) + 1] = vector[3]
    return qt.Qobj(full, dims=[[int(n_tr), int(n_cav)], [1, 1]])


def pure_state_fidelity(target_state: qt.Qobj, actual_state: qt.Qobj) -> float:
    overlap = complex(target_state.overlap(actual_state))
    return float(abs(overlap) ** 2)


def state_leakage(state: qt.Qobj, subspace_indices: Sequence[int]) -> float:
    vector = np.asarray(state.full(), dtype=np.complex128).reshape(-1)
    logical = vector[list(subspace_indices)]
    return float(max(0.0, 1.0 - np.vdot(logical, logical).real))


def serialize_pulse(pulse: Pulse) -> dict[str, Any]:
    return {
        "channel": str(pulse.channel),
        "t0": float(pulse.t0),
        "duration": float(pulse.duration),
        "carrier": float(pulse.carrier),
        "phase": float(pulse.phase),
        "amp": float(pulse.amp),
        "drag": float(pulse.drag),
        "sample_rate": None if pulse.sample_rate is None else float(pulse.sample_rate),
        "label": None if pulse.label is None else str(pulse.label),
    }


def hash_pulses(pulses: Sequence[Pulse]) -> str:
    payload = json_dumps_stable([serialize_pulse(pulse) for pulse in pulses])
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def json_dumps_stable(payload: Any) -> str:
    import json

    return json.dumps(payload, indent=2, sort_keys=True)


def to_jsonable(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, Pulse):
        return serialize_pulse(value)
    if isinstance(value, np.ndarray):
        return to_jsonable(value.tolist())
    if isinstance(value, (np.floating, np.integer, np.bool_)):
        return value.item()
    if isinstance(value, np.complexfloating):
        return {"real": float(value.real), "imag": float(value.imag)}
    if isinstance(value, complex):
        return {"real": float(value.real), "imag": float(value.imag)}
    if isinstance(value, qt.Qobj):
        array = np.asarray(value.full(), dtype=np.complex128)
        return {
            "real": array.real.tolist(),
            "imag": array.imag.tolist(),
            "dims": value.dims,
        }
    if isinstance(value, dict):
        return {str(key): to_jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_jsonable(item) for item in value]
    return value


def unpack_waveform_payload(result: Any) -> tuple[list[Pulse], dict[str, Any], dict[str, Any]]:
    if isinstance(result, dict):
        return list(result.get("pulses", [])), dict(result.get("drive_ops", {})), dict(result.get("meta", {}))
    if isinstance(result, tuple):
        if len(result) == 2:
            pulses, drive_ops = result
            return list(pulses), dict(drive_ops), {}
        if len(result) == 3:
            pulses, drive_ops, meta = result
            return list(pulses), dict(drive_ops), dict(meta)
        raise ValueError("Unsupported waveform payload tuple length.")
    return list(result), {}, {}


def shift_pulses(pulses: Sequence[Pulse], *, offset_s: float, prefix: str) -> list[Pulse]:
    shifted: list[Pulse] = []
    for pulse in pulses:
        label = None if pulse.label is None else f"{prefix}:{pulse.label}"
        shifted.append(replace(pulse, t0=float(pulse.t0) + float(offset_s), label=label))
    return shifted


def exact_cavity_h_target() -> np.ndarray:
    hadamard = np.array([[1.0, 1.0], [1.0, -1.0]], dtype=np.complex128) / np.sqrt(2.0)
    zeros = np.zeros((2, 2), dtype=np.complex128)
    return np.block([[hadamard, zeros], [zeros, hadamard]])


def symbolic_local_target(kind: str, *, n_cav: int) -> np.ndarray:
    if kind == "exact_hc":
        return exact_cavity_h_target()
    subspace = Subspace.custom(full_dim=2 * int(n_cav), indices=logical_subspace_indices(n_cav), labels=LOGICAL_BASIS_LABELS)
    sequence = GateSequence(gates=phase2.local_hc_block(kind, f"target_{kind}", n_cav=n_cav), n_cav=int(n_cav), full_dim=2 * int(n_cav))
    projected = np.asarray(synth_simulate_sequence(sequence, subspace=subspace, backend="ideal").subspace_operator, dtype=np.complex128)
    return np.asarray(polar(projected)[0], dtype=np.complex128)


def surrogate_channel_specs() -> tuple[ModelControlChannelSpec, ...]:
    return (
        ModelControlChannelSpec(
            name="storage_i",
            target="storage",
            quadratures=("I",),
            amplitude_bounds=(-SURROGATE_STORAGE_AMP, SURROGATE_STORAGE_AMP),
            export_channel="storage_i",
        ),
        ModelControlChannelSpec(
            name="storage_q",
            target="storage",
            quadratures=("Q",),
            amplitude_bounds=(-SURROGATE_STORAGE_AMP, SURROGATE_STORAGE_AMP),
            export_channel="storage_q",
        ),
        ModelControlChannelSpec(
            name="qubit_i",
            target="qubit",
            quadratures=("I",),
            amplitude_bounds=(-SURROGATE_QUBIT_AMP, SURROGATE_QUBIT_AMP),
            export_channel="qubit_i",
        ),
        ModelControlChannelSpec(
            name="qubit_q",
            target="qubit",
            quadratures=("Q",),
            amplitude_bounds=(-SURROGATE_QUBIT_AMP, SURROGATE_QUBIT_AMP),
            export_channel="qubit_q",
        ),
    )


def initial_schedule_for_restart(restart_index: int, *, steps: int, seed: int) -> np.ndarray:
    channels = surrogate_channel_specs()
    if restart_index == 0:
        return np.zeros((len(channels), steps), dtype=float)
    rng = np.random.default_rng(seed + restart_index)
    rows: list[np.ndarray] = []
    for spec in channels:
        low, high = spec.amplitude_bounds
        span = 0.20 * max(abs(float(low)), abs(float(high)))
        rows.append(rng.uniform(-span, span, size=steps))
    return np.vstack(rows)


def surrogate_extra_complexity(steps: int, channels: int) -> float:
    return float(0.5 * channels + 0.02 * steps)


def optimize_local_surrogate(target_kind: str) -> SurrogateBlock:
    label = f"{target_kind}_runtime"
    print(f"[phase5] Optimizing surrogate {target_kind} at n_cav={REFERENCE_N_CAV}, n_tr={RUNTIME_N_TR}", flush=True)
    model = build_model(n_cav=REFERENCE_N_CAV, n_tr=RUNTIME_N_TR)
    frame = build_frame(model)
    subspace = logical_subspace(REFERENCE_N_CAV, RUNTIME_N_TR)
    target = symbolic_local_target(target_kind, n_cav=REFERENCE_N_CAV)

    problem = build_control_problem_from_model(
        model,
        frame=frame,
        time_grid=PiecewiseConstantTimeGrid.uniform(steps=SURROGATE_STEPS, dt_s=SURROGATE_DT_S),
        channel_specs=surrogate_channel_specs(),
        objectives=(
            UnitaryObjective(
                target_operator=target,
                subspace=subspace,
                ignore_global_phase=True,
                name=label,
            ),
        ),
        metadata={
            "study_name": "hybrid_unitary_native_entangling_evolution",
            "description": f"Replayable local surrogate for {target_kind}",
            "reference_n_cav": REFERENCE_N_CAV,
            "reference_n_tr": RUNTIME_N_TR,
        },
    )

    best_result = None
    best_fidelity = -math.inf
    restart_summaries: list[dict[str, Any]] = []
    for restart_index in range(SURROGATE_RESTARTS):
        print(f"[phase5]   restart {restart_index + 1}/{SURROGATE_RESTARTS} for {target_kind}", flush=True)
        initial_schedule = initial_schedule_for_restart(restart_index, steps=SURROGATE_STEPS, seed=17)
        result = GrapeSolver(GrapeConfig(maxiter=SURROGATE_MAXITER, seed=21 + restart_index)).solve(
            problem,
            initial_schedule=initial_schedule,
        )
        fidelity = float(result.metrics.get("nominal_fidelity", 0.0))
        print(
            f"[phase5]   restart {restart_index + 1} complete: fidelity={fidelity:.6f}, success={bool(result.success)}",
            flush=True,
        )
        restart_summaries.append(
            {
                "restart_index": restart_index,
                "success": bool(result.success),
                "objective_value": float(result.objective_value),
                "nominal_fidelity": fidelity,
                "message": str(result.message),
            }
        )
        if fidelity > best_fidelity:
            best_fidelity = fidelity
            best_result = result

    if best_result is None:
        raise RuntimeError(f"No GRAPE result was produced for surrogate {target_kind}.")

    pulses, drive_ops, pulse_meta = unpack_waveform_payload(best_result.to_pulses())
    result_path = ARTIFACT_DIR / f"phase5_local_surrogate_{target_kind}_grape.json"
    artifact_path = ARTIFACT_DIR / f"phase5_local_surrogate_{target_kind}.json"
    dump_json(result_path, to_jsonable(best_result.to_payload()))
    artifact_payload = {
        "study_name": "hybrid_unitary_native_entangling_evolution",
        "date_created": "2026-03-24",
        "description": f"Replayable local surrogate for {target_kind} at n_cav={REFERENCE_N_CAV}, n_tr={RUNTIME_N_TR}.",
        "target_kind": target_kind,
        "target_matrix": serialize_complex_matrix(target),
        "reference_n_cav": REFERENCE_N_CAV,
        "reference_n_tr": RUNTIME_N_TR,
        "command_values_shape": list(best_result.schedule.values.shape),
        "nominal_fidelity": float(best_result.metrics.get("nominal_fidelity", 0.0)),
        "objective_value": float(best_result.objective_value),
        "pulse_hash": hash_pulses(pulses),
        "pulses": [serialize_pulse(pulse) for pulse in pulses],
        "drive_ops": {str(key): str(value) for key, value in drive_ops.items()},
        "pulse_metadata": to_jsonable(pulse_meta),
        "restart_summaries": restart_summaries,
        "load_instructions": f"Load {result_path.name} for the full GRAPE payload and {artifact_path.name} for exported runtime pulses.",
        "result_artifact": str(result_path),
    }
    dump_json(artifact_path, to_jsonable(artifact_payload))
    print(
        f"[phase5] Surrogate {target_kind} saved with nominal fidelity {float(best_result.metrics.get('nominal_fidelity', 0.0)):.6f}",
        flush=True,
    )
    return SurrogateBlock(
        label=label,
        target_kind=target_kind,
        target_matrix=np.asarray(target, dtype=np.complex128),
        pulses=pulses,
        drive_ops={str(key): value for key, value in drive_ops.items()},
        pulse_metadata=pulse_meta,
        artifact_path=artifact_path,
        result_path=result_path,
        nominal_fidelity=float(best_result.metrics.get("nominal_fidelity", 0.0)),
        objective_value=float(best_result.objective_value),
        total_duration_s=float(problem.time_grid.duration_s),
        implementation_complexity=surrogate_extra_complexity(SURROGATE_STEPS, len(surrogate_channel_specs())),
        restart_summaries=restart_summaries,
    )


def hq_runtime_segment(name: str, *, model: Any, frame: Any) -> ReplaySegment:
    gates = [
        QubitRotation(name=f"{name}_ry", theta=np.pi / 2.0, phi=np.pi / 2.0, duration=20.0e-9, optimize_time=False),
        QubitRotation(name=f"{name}_rx", theta=np.pi, phi=0.0, duration=20.0e-9, optimize_time=False),
    ]
    return gates_to_segment(name, gates, model=model, frame=frame, kind="hq_runtime")


def idle_wait_gate_to_pulse(name: str, duration_s: float) -> tuple[list[Pulse], dict[str, Any], dict[str, Any]]:
    pulse = Pulse("qubit", 0.0, float(duration_s), square_envelope, amp=0.0, label=name)
    return [pulse], {"qubit": "qubit"}, {"mapping": "Explicit zero-amplitude idle pulse for FreeEvolveCondPhase replay."}


def gates_to_segment(name: str, gates: Sequence[Any], *, model: Any, frame: Any, kind: str) -> ReplaySegment:
    pulses: list[Pulse] = []
    drive_ops: dict[str, Any] = {}
    source_gates: list[dict[str, Any]] = []
    current_t = 0.0
    n_cav = int(model.subsystem_dims[1])
    for gate_index, gate in enumerate(gates):
        if str(gate.type) == "FreeEvolveCondPhase":
            gate_pulses, gate_ops, gate_meta = idle_wait_gate_to_pulse(str(gate.name), float(gate.duration))
        elif str(gate.type) == "PrimitiveGate":
            raise TypeError(f"Unsupported primitive gate in direct replay segment {name}: {gate.name}")
        else:
            primitive = waveform_primitive_from_gate(
                gate,
                index=gate_index,
                frame=frame,
                hilbert_dim=int(np.prod(model.subsystem_dims)),
            )
            gate_pulses, gate_ops, gate_meta = unpack_waveform_payload(primitive.waveform(primitive.runtime_parameters(), model))
        shifted = shift_pulses(gate_pulses, offset_s=current_t, prefix=f"{name}:{gate.name}")
        current_t = max(current_t, max(float(pulse.t1) for pulse in shifted))
        pulses.extend(shifted)
        drive_ops.update(gate_ops)
        parameters = gate.get_parameters(n_cav) if hasattr(gate, "get_parameters") else np.asarray([], dtype=float)
        source_gates.append(
            {
                "gate_name": str(gate.name),
                "gate_type": str(gate.type),
                "duration_s": float(gate.duration),
                "parameters": np.asarray(parameters, dtype=float).tolist(),
                "pulse_metadata": to_jsonable(gate_meta),
            }
        )
    return ReplaySegment(
        name=name,
        kind=kind,
        pulses=pulses,
        drive_ops={str(key): value for key, value in drive_ops.items()},
        duration_s=float(current_t),
        metadata={
            "source_gates": source_gates,
            "pulse_hash": hash_pulses(pulses),
            "pulse_count": len(pulses),
        },
    )


def surrogate_segment(name: str, surrogate: SurrogateBlock) -> ReplaySegment:
    return ReplaySegment(
        name=name,
        kind=surrogate.label,
        pulses=[replace(pulse) for pulse in surrogate.pulses],
        drive_ops=dict(surrogate.drive_ops),
        duration_s=float(surrogate.total_duration_s),
        metadata={
            "surrogate_label": surrogate.label,
            "target_kind": surrogate.target_kind,
            "artifact_path": str(surrogate.artifact_path),
            "result_path": str(surrogate.result_path),
            "nominal_fidelity": float(surrogate.nominal_fidelity),
            "pulse_hash": hash_pulses(surrogate.pulses),
            "pulse_metadata": to_jsonable(surrogate.pulse_metadata),
        },
    )


def local_segment(kind: str, name: str, *, model: Any, frame: Any, surrogates: Mapping[str, SurrogateBlock]) -> ReplaySegment:
    if kind == "exact_hc_runtime":
        return surrogate_segment(name, surrogates["exact_hc"])
    if kind == "A_local_runtime":
        return surrogate_segment(name, surrogates["A_local"])
    if kind == "B_local_direct":
        return gates_to_segment(name, phase2.local_hc_block("B_local", name, n_cav=REFERENCE_N_CAV), model=model, frame=frame, kind=kind)
    raise ValueError(f"Unknown runtime local kind: {kind}")


def entangler_segment(name: str, *, model: Any, frame: Any) -> ReplaySegment:
    return gates_to_segment(name, phase2.native_entangler_block(name), model=model, frame=frame, kind="native_entangler")


def build_runtime_segments(candidate: RuntimeCandidateDef, *, reference_model: Any, frame: Any, surrogates: Mapping[str, SurrogateBlock]) -> list[ReplaySegment]:
    segments: list[ReplaySegment] = []
    if candidate.waits == 1:
        segments.append(hq_runtime_segment("Hq_in_1", model=reference_model, frame=frame))
        segments.append(local_segment(candidate.runtime_inner_kind, "local_inner_1", model=reference_model, frame=frame, surrogates=surrogates))
        segments.append(entangler_segment("D_ent_1", model=reference_model, frame=frame))
        segments.append(local_segment(candidate.runtime_outer_kind, "local_outer_final", model=reference_model, frame=frame, surrogates=surrogates))
        return segments
    if candidate.waits == 2:
        segments.append(hq_runtime_segment("Hq_in_1", model=reference_model, frame=frame))
        segments.append(local_segment(candidate.runtime_inner_kind, "local_inner_1", model=reference_model, frame=frame, surrogates=surrogates))
        segments.append(entangler_segment("D_ent_1", model=reference_model, frame=frame))
        segments.append(hq_runtime_segment("Hq_in_2", model=reference_model, frame=frame))
        segments.append(local_segment(candidate.runtime_inner_kind, "local_inner_2", model=reference_model, frame=frame, surrogates=surrogates))
        segments.append(entangler_segment("D_ent_2", model=reference_model, frame=frame))
        segments.append(local_segment(candidate.runtime_outer_kind, "local_outer_final", model=reference_model, frame=frame, surrogates=surrogates))
        return segments
    raise ValueError(f"Unsupported wait count: {candidate.waits}")


def compile_replay_segments(candidate_label: str, segments: Sequence[ReplaySegment]) -> CompiledReplay:
    pulses: list[Pulse] = []
    drive_ops: dict[str, Any] = {}
    segment_names: list[str] = []
    segment_kinds: list[str] = []
    segment_metadata: list[dict[str, Any]] = []
    segment_end_times: list[float] = []
    current_t = 0.0
    for segment in segments:
        shifted = shift_pulses(segment.pulses, offset_s=current_t, prefix=segment.name)
        end_time = max(float(pulse.t1) for pulse in shifted) if shifted else float(current_t)
        pulses.extend(shifted)
        drive_ops.update(segment.drive_ops)
        segment_names.append(segment.name)
        segment_kinds.append(segment.kind)
        segment_metadata.append(dict(segment.metadata))
        current_t = end_time
        segment_end_times.append(end_time)
    compiler = SequenceCompiler(dt=RUNTIME_DT_S)
    compiled = compiler.compile(pulses, t_end=float(current_t) + RUNTIME_DT_S)
    return CompiledReplay(
        candidate_label=candidate_label,
        pulses=pulses,
        drive_ops={str(key): value for key, value in drive_ops.items()},
        segment_names=segment_names,
        segment_kinds=segment_kinds,
        segment_end_times_s=segment_end_times,
        segment_metadata=segment_metadata,
        total_duration_s=float(current_t),
        compiled=compiled,
    )


def compiled_prefix(bundle: CompiledReplay, segment_count: int) -> Any:
    if segment_count <= 0:
        raise ValueError("segment_count must be positive.")
    if segment_count not in bundle.prefix_cache:
        end_time = float(bundle.segment_end_times_s[segment_count - 1])
        prefix_pulses = [pulse for pulse in bundle.pulses if float(pulse.t1) <= end_time + 1.0e-15]
        compiler = SequenceCompiler(dt=RUNTIME_DT_S)
        bundle.prefix_cache[segment_count] = compiler.compile(prefix_pulses, t_end=end_time + RUNTIME_DT_S)
    return bundle.prefix_cache[segment_count]


def batch_simulate(compiled: Any, drive_ops: Mapping[str, Any], states: Sequence[qt.Qobj], *, model: Any, frame: Any, noise: NoiseSpec | None) -> list[qt.Qobj]:
    active_drive_ops = {str(channel): target for channel, target in drive_ops.items() if str(channel) in compiled.channels}
    session = prepare_simulation(
        model,
        compiled,
        active_drive_ops,
        config=SimulationConfig(frame=frame),
        noise=noise,
        e_ops={},
    )
    results = simulate_batch(session, list(states), max_workers=1)
    return [row.final_state for row in results]


def projected_subspace_operator(states: Sequence[qt.Qobj], indices: Sequence[int]) -> np.ndarray:
    columns: list[np.ndarray] = []
    for state in states:
        vector = np.asarray(state.full(), dtype=np.complex128).reshape(-1)
        columns.append(vector[list(indices)])
    return np.column_stack(columns)


def photon_expectation(state: qt.Qobj) -> float:
    rho_c = reduced_cavity_state(state)
    dim = int(rho_c.shape[0])
    return float(qt.expect(qt.num(dim), rho_c))


def wigner_snapshot(state: qt.Qobj) -> dict[str, np.ndarray]:
    rho_c = reduced_cavity_state(state)
    xvec, yvec, w = cavity_wigner(rho_c, n_points=WIGNER_POINTS, extent=WIGNER_EXTENT)
    return {
        "xvec": np.asarray(xvec, dtype=float),
        "yvec": np.asarray(yvec, dtype=float),
        "w": np.asarray(w, dtype=float),
    }


def wigner_overlap(a: Mapping[str, np.ndarray], b: Mapping[str, np.ndarray]) -> float:
    left = np.asarray(a["w"], dtype=float).reshape(-1)
    right = np.asarray(b["w"], dtype=float).reshape(-1)
    denom = float(np.linalg.norm(left) * np.linalg.norm(right))
    if denom <= 0.0:
        return 0.0
    return float(np.dot(left, right) / denom)


def average_probe_fidelity(final_states: Mapping[str, qt.Qobj], target_states: Mapping[str, qt.Qobj]) -> float:
    values = [pure_state_fidelity(target_states[label], final_state) for label, final_state in final_states.items()]
    return float(np.mean(values)) if values else 0.0


def average_probe_leakage(final_states: Mapping[str, qt.Qobj], indices: Sequence[int]) -> float:
    values = [state_leakage(state, indices) for state in final_states.values()]
    return float(np.mean(values)) if values else 0.0


def target_probe_states(*, n_cav: int, n_tr: int) -> tuple[dict[str, qt.Qobj], dict[str, qt.Qobj]]:
    target = target_unitary_matrix()
    initial: dict[str, qt.Qobj] = {}
    target_states: dict[str, qt.Qobj] = {}
    for label, logical_vector in PROBE_STATES:
        initial[label] = embed_logical_state_multilevel(logical_vector, n_cav=n_cav, n_tr=n_tr)
        target_states[label] = embed_logical_state_multilevel(target @ logical_vector, n_cav=n_cav, n_tr=n_tr)
    return initial, target_states


def symbolic_probe_states(*, n_cav: int) -> tuple[dict[str, qt.Qobj], dict[str, qt.Qobj]]:
    target = target_unitary_matrix()
    initial: dict[str, qt.Qobj] = {}
    target_states: dict[str, qt.Qobj] = {}
    for label, logical_vector in PROBE_STATES:
        initial[label] = embed_logical_state(logical_vector, n_cav=n_cav)
        target_states[label] = embed_logical_state(target @ logical_vector, n_cav=n_cav)
    return initial, target_states


def serialize_complex_matrix(matrix: np.ndarray) -> dict[str, Any]:
    arr = np.asarray(matrix, dtype=np.complex128)
    return {
        "real": arr.real.tolist(),
        "imag": arr.imag.tolist(),
    }


def candidate_base_metrics(*, waits: int, inner_kind: str, outer_kind: str) -> dict[str, Any]:
    sequence = phase2.build_candidate_sequence(waits=waits, inner_kind=inner_kind, outer_kind=outer_kind, n_cav=REFERENCE_N_CAV)
    serialized = sequence.serialize()
    depth, entangling_count, entangling_time_ns, total_duration_ns, active_tones, gate_types = phase2.sequence_cost_summary(serialized)
    return {
        "serialized_sequence": serialized,
        "depth": int(depth),
        "entangling_gate_count": int(entangling_count),
        "entangling_time_ns": float(entangling_time_ns),
        "symbolic_total_duration_ns": float(total_duration_ns),
        "active_tones": int(active_tones),
        "gate_types": gate_types,
        "base_complexity": float(
            implementation_complexity_score(sequence=serialized, active_tones=active_tones, gate_types=gate_types)
        ),
    }


def evaluate_symbolic_candidate(candidate: SymbolicCandidateDef, *, n_cav: int) -> dict[str, Any]:
    subspace = Subspace.custom(full_dim=2 * int(n_cav), indices=logical_subspace_indices(n_cav), labels=LOGICAL_BASIS_LABELS)
    model = build_model(n_cav=n_cav, n_tr=2)
    sequence = phase2.build_candidate_sequence(
        waits=candidate.waits,
        inner_kind=candidate.inner_kind,
        outer_kind=candidate.outer_kind,
        n_cav=n_cav,
    )
    basis_states = [embed_logical_state(logical_vector, n_cav=n_cav) for _, logical_vector in PROBE_STATES[:4]]
    final_basis = sequence.propagate_states(basis_states, backend="ideal", backend_settings={"model": model})
    sub_operator = projected_subspace_operator(final_basis, logical_subspace_indices(n_cav))

    probe_initial, probe_targets = symbolic_probe_states(n_cav=n_cav)
    final_probes = {
        label: sequence.propagate_states([state], backend="ideal", backend_settings={"model": model})[0]
        for label, state in probe_initial.items()
    }
    probe_fidelities = {label: pure_state_fidelity(probe_targets[label], state) for label, state in final_probes.items()}
    probe_leakages = {label: state_leakage(state, logical_subspace_indices(n_cav)) for label, state in final_probes.items()}

    checkpoint_indices = list(range(len(sequence.gates) + 1))
    checkpoint_histories = {
        label: sequence.propagate_states_with_checkpoints([state], checkpoint_indices, backend="ideal", backend_settings={"model": model})
        for label, state in probe_initial.items()
    }
    max_transient_photon = 0.0
    for history in checkpoint_histories.values():
        for step in checkpoint_indices:
            max_transient_photon = max(max_transient_photon, photon_expectation(history[step][0]))

    wigner_target = wigner_snapshot(probe_targets[WIGNER_PROBE])
    wigner_final = wigner_snapshot(final_probes[WIGNER_PROBE])
    bloch_history = [checkpoint_histories[BLOCH_PROBE][index][0] for index in checkpoint_indices]
    bloch = bloch_trajectory_from_states(bloch_history, conditioned_n_levels=[0, 1], probability_threshold=1.0e-8)

    base = candidate_base_metrics(waits=candidate.waits, inner_kind=candidate.inner_kind, outer_kind=candidate.outer_kind)
    return {
        "candidate_label": candidate.label,
        "role": candidate.role,
        "evaluation_level": "symbolic",
        "n_cav": int(n_cav),
        "n_tr": 2,
        "process_fidelity": float(subspace_unitary_fidelity(sub_operator, target_unitary_matrix(), gauge="global")),
        "average_probe_fidelity": float(np.mean(list(probe_fidelities.values()))),
        "average_probe_leakage": float(np.mean(list(probe_leakages.values()))),
        "max_transient_photon": float(max_transient_photon),
        "wigner_target_overlap": float(wigner_overlap(wigner_final, wigner_target)),
        "probe_fidelities": probe_fidelities,
        "probe_leakages": probe_leakages,
        "bloch_final_target_overlap": float(pure_state_fidelity(probe_targets[BLOCH_PROBE], final_probes[BLOCH_PROBE])),
        "wigner_final_target_overlap": float(pure_state_fidelity(probe_targets[WIGNER_PROBE], final_probes[WIGNER_PROBE])),
        "bloch_xyz_final": {
            "x": float(np.asarray(bloch["x"], dtype=float)[-1]),
            "y": float(np.asarray(bloch["y"], dtype=float)[-1]),
            "z": float(np.asarray(bloch["z"], dtype=float)[-1]),
        },
        **base,
        "weighted_cost": float(
            candidate_weighted_cost(
                {
                    "pulse_fidelity": float(subspace_unitary_fidelity(sub_operator, target_unitary_matrix(), gauge="global")),
                    "entangling_gate_count": base["entangling_gate_count"],
                    "entangling_time_ns": base["entangling_time_ns"],
                    "depth": base["depth"],
                    "implementation_complexity": base["base_complexity"],
                    "pulse_leakage_average": float(np.mean(list(probe_leakages.values()))),
                },
                CostWeights(),
            )
        ),
    }


def evaluate_runtime_candidate(
    candidate: RuntimeCandidateDef,
    bundle: CompiledReplay,
    *,
    n_cav: int,
    surrogate_complexity: float,
    base_metrics: Mapping[str, Any],
    include_noise: bool,
) -> dict[str, Any]:
    model = build_model(n_cav=n_cav, n_tr=RUNTIME_N_TR)
    frame = build_frame(model)
    indices = logical_subspace_indices(n_cav)
    probe_initial, probe_targets = target_probe_states(n_cav=n_cav, n_tr=RUNTIME_N_TR)

    basis_states = [probe_initial[label] for label, _ in PROBE_STATES[:4]]
    closed_basis = batch_simulate(bundle.compiled, bundle.drive_ops, basis_states, model=model, frame=frame, noise=None)
    sub_operator = projected_subspace_operator(closed_basis, indices)
    process_fidelity = float(subspace_unitary_fidelity(sub_operator, target_unitary_matrix(), gauge="global"))

    closed_probe_states = batch_simulate(
        bundle.compiled,
        bundle.drive_ops,
        list(probe_initial.values()),
        model=model,
        frame=frame,
        noise=None,
    )
    closed_probe_map = {label: state for (label, _), state in zip(PROBE_STATES, closed_probe_states, strict=True)}
    closed_probe_fidelities = {
        label: pure_state_fidelity(probe_targets[label], state) for label, state in closed_probe_map.items()
    }
    closed_probe_leakages = {label: state_leakage(state, indices) for label, state in closed_probe_map.items()}

    noisy_probe_fidelities: dict[str, float] = {}
    noisy_probe_leakages: dict[str, float] = {}
    if include_noise:
        noisy_probe_states = batch_simulate(
            bundle.compiled,
            bundle.drive_ops,
            list(probe_initial.values()),
            model=model,
            frame=frame,
            noise=NOMINAL_NOISE,
        )
        noisy_probe_map = {label: state for (label, _), state in zip(PROBE_STATES, noisy_probe_states, strict=True)}
        noisy_probe_fidelities = {
            label: pure_state_fidelity(probe_targets[label], state) for label, state in noisy_probe_map.items()
        }
        noisy_probe_leakages = {label: state_leakage(state, indices) for label, state in noisy_probe_map.items()}

    checkpoint_histories: dict[str, list[qt.Qobj]] = {}
    max_transient_photon = 0.0
    for label, initial_state in probe_initial.items():
        history: list[qt.Qobj] = [initial_state]
        for segment_count in range(1, len(bundle.segment_names) + 1):
            prefix = compiled_prefix(bundle, segment_count)
            final_state = batch_simulate(prefix, bundle.drive_ops, [initial_state], model=model, frame=frame, noise=None)[0]
            history.append(final_state)
        checkpoint_histories[label] = history
        for state in history:
            max_transient_photon = max(max_transient_photon, photon_expectation(state))

    bloch_history = checkpoint_histories[BLOCH_PROBE]
    bloch = bloch_trajectory_from_states(bloch_history, conditioned_n_levels=[0, 1], probability_threshold=1.0e-8)
    wigner_target = wigner_snapshot(probe_targets[WIGNER_PROBE])
    wigner_final = wigner_snapshot(closed_probe_map[WIGNER_PROBE])

    implementation_complexity = float(base_metrics["base_complexity"]) + float(surrogate_complexity)
    weighted_cost = float(
        candidate_weighted_cost(
            {
                "pulse_fidelity": process_fidelity,
                "entangling_gate_count": base_metrics["entangling_gate_count"],
                "entangling_time_ns": base_metrics["entangling_time_ns"],
                "depth": base_metrics["depth"],
                "implementation_complexity": implementation_complexity,
                "pulse_leakage_average": float(np.mean(list(closed_probe_leakages.values()))),
            },
            CostWeights(),
        )
    )
    return {
        "candidate_label": candidate.label,
        "role": candidate.role,
        "symbolic_reference": candidate.symbolic_reference,
        "evaluation_level": "runtime",
        "n_cav": int(n_cav),
        "n_tr": RUNTIME_N_TR,
        "process_fidelity": process_fidelity,
        "average_probe_fidelity": float(np.mean(list(closed_probe_fidelities.values()))),
        "average_probe_leakage": float(np.mean(list(closed_probe_leakages.values()))),
        "noisy_average_probe_fidelity": None if not noisy_probe_fidelities else float(np.mean(list(noisy_probe_fidelities.values()))),
        "noisy_average_probe_leakage": None if not noisy_probe_leakages else float(np.mean(list(noisy_probe_leakages.values()))),
        "max_transient_photon": float(max_transient_photon),
        "wigner_target_overlap": float(wigner_overlap(wigner_final, wigner_target)),
        "probe_fidelities": closed_probe_fidelities,
        "probe_leakages": closed_probe_leakages,
        "noisy_probe_fidelities": noisy_probe_fidelities,
        "noisy_probe_leakages": noisy_probe_leakages,
        "bloch_final_target_overlap": float(pure_state_fidelity(probe_targets[BLOCH_PROBE], closed_probe_map[BLOCH_PROBE])),
        "wigner_final_target_overlap": float(pure_state_fidelity(probe_targets[WIGNER_PROBE], closed_probe_map[WIGNER_PROBE])),
        "bloch_xyz_final": {
            "x": float(np.asarray(bloch["x"], dtype=float)[-1]),
            "y": float(np.asarray(bloch["y"], dtype=float)[-1]),
            "z": float(np.asarray(bloch["z"], dtype=float)[-1]),
        },
        "checkpoint_histories": checkpoint_histories,
        "wigner_final": wigner_final,
        "wigner_target": wigner_target,
        "compiled_artifact": {
            "pulse_hash": hash_pulses(bundle.pulses),
            "pulse_count": len(bundle.pulses),
            "active_channels": sorted(str(key) for key in bundle.drive_ops),
            "total_duration_ns": 1.0e9 * float(bundle.total_duration_s),
            "segment_names": list(bundle.segment_names),
            "segment_kinds": list(bundle.segment_kinds),
        },
        **base_metrics,
        "runtime_total_duration_ns": 1.0e9 * float(bundle.total_duration_s),
        "implementation_complexity": implementation_complexity,
        "weighted_cost": weighted_cost,
    }


def weight_sensitivity_rows(runtime_rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    base = CostWeights()
    rows: list[dict[str, Any]] = []
    factors = (0.75, 1.0, 1.25)
    for infidelity_factor in factors:
        for time_factor in factors:
            for complexity_factor in factors:
                weights = CostWeights(
                    infidelity=base.infidelity * infidelity_factor,
                    entangling_gate_count=base.entangling_gate_count,
                    entangling_time=base.entangling_time * time_factor,
                    depth=base.depth,
                    implementation_complexity=base.implementation_complexity * complexity_factor,
                    leakage=base.leakage,
                )
                costs = []
                for row in runtime_rows:
                    payload = {
                        "pulse_fidelity": float(row["process_fidelity"]),
                        "entangling_gate_count": float(row["entangling_gate_count"]),
                        "entangling_time_ns": float(row["entangling_time_ns"]),
                        "depth": float(row["depth"]),
                        "implementation_complexity": float(row["implementation_complexity"]),
                        "pulse_leakage_average": float(row["average_probe_leakage"]),
                    }
                    costs.append((float(candidate_weighted_cost(payload, weights)), str(row["candidate_label"])))
                costs.sort(key=lambda item: item[0])
                rows.append(
                    {
                        "infidelity_factor": infidelity_factor,
                        "time_factor": time_factor,
                        "complexity_factor": complexity_factor,
                        "winner": costs[0][1],
                        "ordered_labels": [label for _, label in costs],
                    }
                )
    return rows


def save_csv(path: Path, rows: Sequence[Mapping[str, Any]], *, fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames))
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key) for key in fieldnames})


def convergence_summary(rows: Sequence[Mapping[str, Any]]) -> dict[str, dict[str, float]]:
    by_label: dict[str, list[Mapping[str, Any]]] = {}
    for row in rows:
        by_label.setdefault(str(row["candidate_label"]), []).append(row)
    summary: dict[str, dict[str, float]] = {}
    for label, group in by_label.items():
        process_values = [float(row["process_fidelity"]) for row in group]
        probe_values = [float(row["average_probe_fidelity"]) for row in group]
        leak_values = [float(row["average_probe_leakage"]) for row in group]
        wigner_values = [float(row["wigner_target_overlap"]) for row in group]
        summary[label] = {
            "process_fidelity_span": float(max(process_values) - min(process_values)),
            "probe_fidelity_span": float(max(probe_values) - min(probe_values)),
            "leakage_span": float(max(leak_values) - min(leak_values)),
            "wigner_overlap_span": float(max(wigner_values) - min(wigner_values)),
        }
    return summary


def comparison_rows(
    symbolic_nominal: Mapping[str, Mapping[str, Any]],
    runtime_nominal: Mapping[str, Mapping[str, Any]],
    runtime_convergence: Mapping[str, Mapping[str, float]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for candidate in SYMBOLIC_CANDIDATES:
        symbolic = symbolic_nominal.get(candidate.label)
        runtime = None
        for runtime_candidate in RUNTIME_CANDIDATES:
            if runtime_candidate.symbolic_reference == candidate.label:
                runtime = runtime_nominal.get(runtime_candidate.label)
        convergence = None if runtime is None else runtime_convergence.get(str(runtime["candidate_label"]))
        rows.append(
            {
                "candidate_label": candidate.label,
                "role": candidate.role,
                "symbolic_process_fidelity": None if symbolic is None else float(symbolic["process_fidelity"]),
                "symbolic_avg_probe_fidelity": None if symbolic is None else float(symbolic["average_probe_fidelity"]),
                "symbolic_leakage": None if symbolic is None else float(symbolic["average_probe_leakage"]),
                "runtime_candidate": None if runtime is None else str(runtime["candidate_label"]),
                "runtime_process_fidelity": None if runtime is None else float(runtime["process_fidelity"]),
                "runtime_avg_probe_fidelity": None if runtime is None else float(runtime["average_probe_fidelity"]),
                "runtime_leakage": None if runtime is None else float(runtime["average_probe_leakage"]),
                "runtime_noisy_avg_probe_fidelity": None if runtime is None else runtime["noisy_average_probe_fidelity"],
                "runtime_max_transient_photon": None if runtime is None else float(runtime["max_transient_photon"]),
                "runtime_wigner_target_overlap": None if runtime is None else float(runtime["wigner_target_overlap"]),
                "entangling_gate_count": None if symbolic is None else int(symbolic["entangling_gate_count"]),
                "entangling_time_ns": None if symbolic is None else float(symbolic["entangling_time_ns"]),
                "symbolic_total_duration_ns": None if symbolic is None else float(symbolic["symbolic_total_duration_ns"]),
                "runtime_total_duration_ns": None if runtime is None else float(runtime["runtime_total_duration_ns"]),
                "implementation_complexity": None if runtime is None else float(runtime["implementation_complexity"]),
                "truncation_process_span": None if convergence is None else float(convergence["process_fidelity_span"]),
                "truncation_probe_span": None if convergence is None else float(convergence["probe_fidelity_span"]),
                "truncation_leakage_span": None if convergence is None else float(convergence["leakage_span"]),
                "truncation_wigner_span": None if convergence is None else float(convergence["wigner_overlap_span"]),
            }
        )
    return rows


def save_comparison_figure(rows: Sequence[Mapping[str, Any]]) -> None:
    apply_publication_style()
    labels = [str(row["candidate_label"]) for row in rows]
    symbolic = [float(row["symbolic_process_fidelity"] or 0.0) for row in rows]
    runtime = [float(row["runtime_process_fidelity"] or 0.0) for row in rows]
    runtime_probe = [float(row["runtime_avg_probe_fidelity"] or 0.0) for row in rows]
    x = np.arange(len(rows), dtype=float)
    width = 0.26

    fig, ax = plt.subplots(figsize=(11.5, 4.8))
    ax.bar(x - width, symbolic, width=width, label="symbolic process", color="#4477AA")
    ax.bar(x, runtime, width=width, label="runtime process", color="#228833")
    ax.bar(x + width, runtime_probe, width=width, label="runtime probe average", color="#EE6677")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right")
    ax.set_ylabel("Fidelity")
    ax.set_ylim(0.0, 1.05)
    ax.set_title("Nominal candidate comparison at n_cav = 12")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(COMPARISON_PNG, dpi=300, bbox_inches="tight")
    fig.savefig(COMPARISON_PDF, bbox_inches="tight")
    plt.close(fig)


def save_convergence_figure(rows: Sequence[Mapping[str, Any]]) -> None:
    apply_publication_style()
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.8))

    runtime_rows = [row for row in rows if str(row["evaluation_level"]) == "runtime"]
    labels = sorted({str(row["candidate_label"]) for row in runtime_rows})
    palette = ["#4477AA", "#EE6677", "#228833", "#CCBB44"]

    for index, label in enumerate(labels):
        group = sorted((row for row in runtime_rows if str(row["candidate_label"]) == label), key=lambda row: int(row["n_cav"]))
        axes[0].plot(
            [int(row["n_cav"]) for row in group],
            [float(row["process_fidelity"]) for row in group],
            marker="o",
            color=palette[index % len(palette)],
            label=label,
        )
        axes[1].plot(
            [int(row["n_cav"]) for row in group],
            [float(row["average_probe_leakage"]) for row in group],
            marker="o",
            color=palette[index % len(palette)],
            label=label,
        )

    axes[0].set_xlabel("Cavity truncation n_cav")
    axes[0].set_ylabel("Runtime process fidelity")
    axes[0].set_ylim(0.0, 1.02)
    axes[0].grid(True, alpha=0.25)

    axes[1].set_xlabel("Cavity truncation n_cav")
    axes[1].set_ylabel("Average probe leakage")
    axes[1].grid(True, alpha=0.25)
    axes[1].legend(frameon=False)
    fig.suptitle("Runtime truncation convergence")
    fig.tight_layout()
    fig.savefig(CONVERGENCE_PNG, dpi=300, bbox_inches="tight")
    fig.savefig(CONVERGENCE_PDF, bbox_inches="tight")
    plt.close(fig)


def save_weight_sensitivity_figure(rows: Sequence[Mapping[str, Any]]) -> None:
    apply_publication_style()
    counts: dict[str, int] = {}
    for row in rows:
        counts[str(row["winner"])] = counts.get(str(row["winner"]), 0) + 1
    ordered = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    labels = [label for label, _ in ordered]
    values = [count for _, count in ordered]
    fig, ax = plt.subplots(figsize=(8.6, 4.2))
    ax.bar(labels, values, color="#4477AA", edgecolor="black", linewidth=0.5)
    ax.set_ylabel("Winner count")
    ax.set_title("Weighted-score sensitivity over runtime candidates")
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(SENSITIVITY_PNG, dpi=300, bbox_inches="tight")
    fig.savefig(SENSITIVITY_PDF, bbox_inches="tight")
    plt.close(fig)


def save_bloch_figure(candidate_label: str, history: Sequence[qt.Qobj], segment_names: Sequence[str]) -> dict[str, str]:
    apply_publication_style()
    bloch = bloch_trajectory_from_states(list(history), conditioned_n_levels=[0, 1], probability_threshold=1.0e-8)
    x = np.arange(len(history), dtype=int)
    labels = ["init", *list(segment_names)]
    fig, ax = plt.subplots(figsize=(9.8, 4.0))
    ax.plot(x, np.asarray(bloch["x"], dtype=float), marker="o", label="X", color="#4477AA")
    ax.plot(x, np.asarray(bloch["y"], dtype=float), marker="s", label="Y", color="#EE6677")
    ax.plot(x, np.asarray(bloch["z"], dtype=float), marker="^", label="Z", color="#228833")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=25, ha="right")
    ax.set_ylim(-1.05, 1.05)
    ax.set_ylabel("Bloch expectation")
    ax.set_title(f"{candidate_label}: Bloch trajectory for {BLOCH_PROBE}")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    png = FIG_DIR / f"phase5_bloch_{candidate_label}_{BLOCH_PROBE}.png"
    pdf = FIG_DIR / f"phase5_bloch_{candidate_label}_{BLOCH_PROBE}.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return {"png": str(png), "pdf": str(pdf)}


def save_wigner_grid(candidate_label: str, history: Sequence[qt.Qobj], segment_names: Sequence[str]) -> dict[str, str]:
    apply_publication_style()
    snapshots = [wigner_snapshot(state) for state in history]
    vmax = max(float(np.max(np.abs(snapshot["w"]))) for snapshot in snapshots)
    labels = ["init", *list(segment_names)]
    ncols = min(4, len(snapshots))
    nrows = int(math.ceil(len(snapshots) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.1 * ncols, 2.8 * nrows), squeeze=False)
    for index, snapshot in enumerate(snapshots):
        row = index // ncols
        col = index % ncols
        ax = axes[row][col]
        image = ax.imshow(
            snapshot["w"],
            origin="lower",
            extent=[snapshot["xvec"][0], snapshot["xvec"][-1], snapshot["yvec"][0], snapshot["yvec"][-1]],
            cmap="RdBu_r",
            vmin=-vmax,
            vmax=vmax,
            aspect="equal",
        )
        ax.set_title(labels[index], fontsize=8)
        ax.set_xlabel("Re(alpha)")
        ax.set_ylabel("Im(alpha)")
    for index in range(len(snapshots), nrows * ncols):
        row = index // ncols
        col = index % ncols
        axes[row][col].axis("off")
    fig.colorbar(image, ax=axes.ravel().tolist(), shrink=0.88, label="Wigner")
    fig.suptitle(f"{candidate_label}: Wigner snapshots for {WIGNER_PROBE}")
    fig.tight_layout()
    png = FIG_DIR / f"phase5_wigner_{candidate_label}_{WIGNER_PROBE}.png"
    pdf = FIG_DIR / f"phase5_wigner_{candidate_label}_{WIGNER_PROBE}.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return {"png": str(png), "pdf": str(pdf)}


def save_final_wigner_compare(candidate_label: str, target_snapshot: Mapping[str, np.ndarray], realized_snapshot: Mapping[str, np.ndarray]) -> dict[str, str]:
    apply_publication_style()
    difference = np.asarray(realized_snapshot["w"], dtype=float) - np.asarray(target_snapshot["w"], dtype=float)
    vmax = max(
        float(np.max(np.abs(np.asarray(target_snapshot["w"], dtype=float)))),
        float(np.max(np.abs(np.asarray(realized_snapshot["w"], dtype=float)))),
    )
    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.6), squeeze=False)
    panels = [
        (target_snapshot["w"], "ideal target"),
        (realized_snapshot["w"], "realized"),
        (difference, "difference"),
    ]
    for col, (data, title) in enumerate(panels):
        ax = axes[0][col]
        if title == "difference":
            local_vmax = float(np.max(np.abs(data)))
        else:
            local_vmax = vmax
        image = ax.imshow(
            data,
            origin="lower",
            extent=[target_snapshot["xvec"][0], target_snapshot["xvec"][-1], target_snapshot["yvec"][0], target_snapshot["yvec"][-1]],
            cmap="RdBu_r",
            vmin=-local_vmax,
            vmax=local_vmax,
            aspect="equal",
        )
        ax.set_title(title)
        ax.set_xlabel("Re(alpha)")
        ax.set_ylabel("Im(alpha)")
        fig.colorbar(image, ax=ax, shrink=0.80)
    fig.suptitle(f"{candidate_label}: final Wigner comparison for {WIGNER_PROBE}")
    fig.tight_layout()
    png = FIG_DIR / f"phase5_wigner_compare_{candidate_label}_{WIGNER_PROBE}.png"
    pdf = FIG_DIR / f"phase5_wigner_compare_{candidate_label}_{WIGNER_PROBE}.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return {"png": str(png), "pdf": str(pdf)}


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)

    print("[phase5] Starting runtime validation pipeline", flush=True)

    surrogates = {
        "exact_hc": optimize_local_surrogate("exact_hc"),
        "A_local": optimize_local_surrogate("A_local"),
    }
    print("[phase5] Built replayable local surrogates", flush=True)

    reference_model = build_model(n_cav=REFERENCE_N_CAV, n_tr=RUNTIME_N_TR)
    reference_frame = build_frame(reference_model)

    runtime_bundles: dict[str, CompiledReplay] = {}
    runtime_base_metrics: dict[str, dict[str, Any]] = {}
    runtime_surrogate_complexity: dict[str, float] = {}

    for candidate in RUNTIME_CANDIDATES:
        print(f"[phase5] Compiling runtime candidate {candidate.label}", flush=True)
        segments = build_runtime_segments(candidate, reference_model=reference_model, frame=reference_frame, surrogates=surrogates)
        bundle = compile_replay_segments(candidate.label, segments)
        runtime_bundles[candidate.label] = bundle
        base = candidate_base_metrics(
            waits=candidate.waits,
            inner_kind=candidate.symbolic_inner_kind,
            outer_kind=candidate.symbolic_outer_kind,
        )
        runtime_base_metrics[candidate.label] = base
        extra_complexity = 0.0
        if candidate.runtime_inner_kind.endswith("_runtime"):
            extra_complexity += surrogates[candidate.runtime_inner_kind.replace("_runtime", "")].implementation_complexity
        if candidate.runtime_outer_kind.endswith("_runtime"):
            extra_complexity += surrogates[candidate.runtime_outer_kind.replace("_runtime", "")].implementation_complexity
        runtime_surrogate_complexity[candidate.label] = extra_complexity
        candidate_payload = {
            "study_name": "hybrid_unitary_native_entangling_evolution",
            "date_created": "2026-03-24",
            "description": f"Runtime replay decomposition for {candidate.label}.",
            "symbolic_reference": candidate.symbolic_reference,
            "decomposition": [
                {
                    "segment_name": name,
                    "segment_kind": kind,
                    "end_time_s": float(end_time),
                    "metadata": metadata,
                }
                for name, kind, end_time, metadata in zip(
                    bundle.segment_names,
                    bundle.segment_kinds,
                    bundle.segment_end_times_s,
                    bundle.segment_metadata,
                    strict=True,
                )
            ],
            "pulses": [serialize_pulse(pulse) for pulse in bundle.pulses],
            "drive_ops": {str(key): str(value) for key, value in bundle.drive_ops.items()},
            "pulse_hash": hash_pulses(bundle.pulses),
            "parameters": {
                "reference_n_cav": REFERENCE_N_CAV,
                "reference_n_tr": RUNTIME_N_TR,
                "runtime_dt_s": RUNTIME_DT_S,
            },
            "load_instructions": f"Rebuild pulses from {candidate.label} or replay the serialized pulse list through cqed_sim.sim.prepare_simulation.",
        }
        dump_json(ARTIFACT_DIR / f"phase5_runtime_candidate_{candidate.label}.json", to_jsonable(candidate_payload))

    symbolic_rows: list[dict[str, Any]] = []
    for candidate in SYMBOLIC_CANDIDATES:
        for n_cav in TRUNCATION_SWEEP:
            print(f"[phase5] Evaluating symbolic candidate {candidate.label} at n_cav={n_cav}", flush=True)
            symbolic_rows.append(evaluate_symbolic_candidate(candidate, n_cav=n_cav))

    runtime_rows: list[dict[str, Any]] = []
    for candidate in RUNTIME_CANDIDATES:
        for n_cav in TRUNCATION_SWEEP:
            print(f"[phase5] Evaluating runtime candidate {candidate.label} at n_cav={n_cav}", flush=True)
            runtime_rows.append(
                evaluate_runtime_candidate(
                    candidate,
                    runtime_bundles[candidate.label],
                    n_cav=n_cav,
                    surrogate_complexity=runtime_surrogate_complexity[candidate.label],
                    base_metrics=runtime_base_metrics[candidate.label],
                    include_noise=(n_cav == REFERENCE_N_CAV),
                )
            )

    symbolic_by_label = {row["candidate_label"]: row for row in symbolic_rows if int(row["n_cav"]) == REFERENCE_N_CAV}
    runtime_by_label = {row["candidate_label"]: row for row in runtime_rows if int(row["n_cav"]) == REFERENCE_N_CAV}
    runtime_convergence = convergence_summary(runtime_rows)
    comparison = comparison_rows(symbolic_by_label, runtime_by_label, runtime_convergence)
    sensitivity = weight_sensitivity_rows(list(runtime_by_label.values()))

    save_csv(
        SYMBOLIC_CSV,
        symbolic_rows,
        fieldnames=(
            "candidate_label",
            "role",
            "evaluation_level",
            "n_cav",
            "n_tr",
            "process_fidelity",
            "average_probe_fidelity",
            "average_probe_leakage",
            "max_transient_photon",
            "wigner_target_overlap",
            "entangling_gate_count",
            "entangling_time_ns",
            "symbolic_total_duration_ns",
            "base_complexity",
            "weighted_cost",
        ),
    )
    save_csv(
        RUNTIME_CSV,
        runtime_rows,
        fieldnames=(
            "candidate_label",
            "role",
            "evaluation_level",
            "symbolic_reference",
            "n_cav",
            "n_tr",
            "process_fidelity",
            "average_probe_fidelity",
            "average_probe_leakage",
            "noisy_average_probe_fidelity",
            "noisy_average_probe_leakage",
            "max_transient_photon",
            "wigner_target_overlap",
            "entangling_gate_count",
            "entangling_time_ns",
            "runtime_total_duration_ns",
            "implementation_complexity",
            "weighted_cost",
        ),
    )
    save_csv(
        COMPARISON_CSV,
        comparison,
        fieldnames=(
            "candidate_label",
            "role",
            "symbolic_process_fidelity",
            "symbolic_avg_probe_fidelity",
            "symbolic_leakage",
            "runtime_candidate",
            "runtime_process_fidelity",
            "runtime_avg_probe_fidelity",
            "runtime_leakage",
            "runtime_noisy_avg_probe_fidelity",
            "runtime_max_transient_photon",
            "runtime_wigner_target_overlap",
            "entangling_gate_count",
            "entangling_time_ns",
            "symbolic_total_duration_ns",
            "runtime_total_duration_ns",
            "implementation_complexity",
            "truncation_process_span",
            "truncation_probe_span",
            "truncation_leakage_span",
            "truncation_wigner_span",
        ),
    )

    convergence_rows: list[dict[str, Any]] = []
    for row in runtime_rows:
        spans = runtime_convergence[str(row["candidate_label"])]
        convergence_rows.append(
            {
                "candidate_label": row["candidate_label"],
                "n_cav": row["n_cav"],
                "process_fidelity": row["process_fidelity"],
                "average_probe_fidelity": row["average_probe_fidelity"],
                "average_probe_leakage": row["average_probe_leakage"],
                "max_transient_photon": row["max_transient_photon"],
                "wigner_target_overlap": row["wigner_target_overlap"],
                "process_fidelity_span": spans["process_fidelity_span"],
                "probe_fidelity_span": spans["probe_fidelity_span"],
                "leakage_span": spans["leakage_span"],
                "wigner_overlap_span": spans["wigner_overlap_span"],
            }
        )
    save_csv(
        CONVERGENCE_CSV,
        convergence_rows,
        fieldnames=(
            "candidate_label",
            "n_cav",
            "process_fidelity",
            "average_probe_fidelity",
            "average_probe_leakage",
            "max_transient_photon",
            "wigner_target_overlap",
            "process_fidelity_span",
            "probe_fidelity_span",
            "leakage_span",
            "wigner_overlap_span",
        ),
    )
    save_csv(
        WEIGHT_SENSITIVITY_CSV,
        sensitivity,
        fieldnames=("infidelity_factor", "time_factor", "complexity_factor", "winner", "ordered_labels"),
    )

    save_comparison_figure(comparison)
    save_convergence_figure(runtime_rows)
    save_weight_sensitivity_figure(sensitivity)

    diagnostic_figures: dict[str, dict[str, dict[str, str]]] = {}
    for label in DIAGNOSTIC_RUNTIME_LABELS:
        row = runtime_by_label[label]
        history = row["checkpoint_histories"][BLOCH_PROBE]
        wigner_history = row["checkpoint_histories"][WIGNER_PROBE]
        diagnostic_figures[label] = {
            "bloch": save_bloch_figure(label, history, runtime_bundles[label].segment_names),
            "wigner_grid": save_wigner_grid(label, wigner_history, runtime_bundles[label].segment_names),
            "wigner_compare": save_final_wigner_compare(label, row["wigner_target"], row["wigner_final"]),
        }

    payload = {
        "metadata": {
            "description": "Phase 5 runtime validation with replayable surrogates, truncation convergence, and noisy replay.",
            "reference_n_cav": REFERENCE_N_CAV,
            "truncation_sweep": list(TRUNCATION_SWEEP),
            "runtime_n_tr": RUNTIME_N_TR,
            "runtime_dt_s": RUNTIME_DT_S,
            "surrogate_steps": SURROGATE_STEPS,
            "surrogate_dt_s": SURROGATE_DT_S,
            "surrogate_restarts": SURROGATE_RESTARTS,
            "surrogate_maxiter": SURROGATE_MAXITER,
            "nominal_noise": {
                "t1_s": QUBIT_T1_S,
                "tphi_s": QUBIT_TPHI_S,
                "kappa": 1.0 / CAVITY_T1_S,
            },
            "comparison_csv": str(COMPARISON_CSV),
            "convergence_csv": str(CONVERGENCE_CSV),
            "weight_sensitivity_csv": str(WEIGHT_SENSITIVITY_CSV),
            "comparison_figures": {"png": str(COMPARISON_PNG), "pdf": str(COMPARISON_PDF)},
            "convergence_figures": {"png": str(CONVERGENCE_PNG), "pdf": str(CONVERGENCE_PDF)},
            "sensitivity_figures": {"png": str(SENSITIVITY_PNG), "pdf": str(SENSITIVITY_PDF)},
        },
        "surrogates": {
            key: {
                "label": block.label,
                "target_kind": block.target_kind,
                "artifact_path": str(block.artifact_path),
                "result_path": str(block.result_path),
                "nominal_fidelity": float(block.nominal_fidelity),
                "objective_value": float(block.objective_value),
                "total_duration_s": float(block.total_duration_s),
                "implementation_complexity": float(block.implementation_complexity),
                "pulse_hash": hash_pulses(block.pulses),
                "pulse_count": len(block.pulses),
                "pulse_metadata": to_jsonable(block.pulse_metadata),
                "restart_summaries": block.restart_summaries,
            }
            for key, block in surrogates.items()
        },
        "symbolic_rows": symbolic_rows,
        "runtime_rows": [
            {
                **{key: value for key, value in row.items() if key not in {"checkpoint_histories", "wigner_final", "wigner_target"}},
                "checkpoint_histories": {
                    label: {
                        "segment_count": len(history) - 1,
                        "max_photon": max(photon_expectation(state) for state in history),
                    }
                    for label, history in row["checkpoint_histories"].items()
                },
                "wigner_final": to_jsonable(row["wigner_final"]),
                "wigner_target": to_jsonable(row["wigner_target"]),
            }
            for row in runtime_rows
        ],
        "comparison_rows": comparison,
        "runtime_convergence_summary": runtime_convergence,
        "weight_sensitivity": sensitivity,
        "diagnostic_figures": diagnostic_figures,
    }
    dump_json(PHASE5_JSON, to_jsonable(payload))

    print("Phase 5 runtime validation summary")
    print("=" * 110)
    print(f"{'Candidate':<32} {'n_cav':>5} {'F_proc':>10} {'F_probe':>10} {'Leak':>9} {'Noisy':>10}")
    print("-" * 110)
    for row in sorted(runtime_rows, key=lambda item: (int(item["n_cav"]), -float(item["process_fidelity"]))):
        noisy = row["noisy_average_probe_fidelity"]
        noisy_text = "n/a" if noisy is None else f"{float(noisy):.4f}"
        print(
            f"{row['candidate_label']:<32} {int(row['n_cav']):>5d} {float(row['process_fidelity']):>10.4f} "
            f"{float(row['average_probe_fidelity']):>10.4f} {float(row['average_probe_leakage']):>9.4f} {noisy_text:>10}"
        )
    print(f"\nWrote {PHASE5_JSON}")
    print(f"Wrote {COMPARISON_CSV}")
    print(f"Wrote {CONVERGENCE_CSV}")
    print(f"Wrote {WEIGHT_SENSITIVITY_CSV}")


if __name__ == "__main__":
    main()