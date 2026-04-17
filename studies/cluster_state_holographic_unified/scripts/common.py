"""
common.py — Shared infrastructure for the unified cluster-state holographic study.

This module defines all physical parameters, the target unitary, gate builders,
UnitarySynthesizer wrappers, and GRAPE helpers used by both the
reproducibility notebook and any future re-optimization scripts.

Adapted from holographic_cluster_state_ideal_gate_followup/scripts/common.py
with paths updated to point to the unified study root.
"""

from __future__ import annotations

import copy
import importlib
import itertools
import json
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np
import runtime_compat  # noqa: F401
import qutip as qt

from cqed_sim import (
    DispersiveTransmonCavityModel,
    FrameSpec,
    ModelControlChannelSpec,
    PiecewiseConstantTimeGrid,
    SequenceCompiler,
    SimulationConfig,
    simulate_sequence as pulse_simulate_sequence,
)
from cqed_sim.core.conventions import qubit_cavity_block_indices
from cqed_sim.core.ideal_gates import qubit_rotation_xy
from cqed_sim.optimal_control import (
    GrapeConfig,
    GrapeSolver,
    LeakagePenalty as OCLeakagePenalty,
    UnitaryObjective as OCUnitaryObjective,
    build_control_problem_from_model,
)
from cqed_sim.sim.noise import NoiseSpec, pure_dephasing_time_from_t1_t2
from cqed_sim.unitary_synthesis import (
    Displacement,
    ExecutionOptions,
    FreeEvolveCondPhase,
    GateSequence,
    LeakagePenalty,
    MultiObjective,
    PrimitiveGate,
    QubitRotation,
    Subspace,
    TargetUnitary,
    UnitarySynthesizer,
    simulate_sequence as synth_simulate_sequence,
    subspace_unitary_fidelity,
)
from cqed_sim.unitary_synthesis.sequence import DriftPhaseModel, drift_phase_table, drift_phase_unitary
from cqed_sim.unitary_synthesis.targets import make_target

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
STUDY_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = STUDY_ROOT / "data"
FIG_DIR = STUDY_ROOT / "figures"
ARTIFACT_DIR = STUDY_ROOT / "artifacts"
REPORT_DIR = STUDY_ROOT / "report"

# Paths to the source studies for loading pre-computed results
SOURCE_IDEAL_GATE = STUDY_ROOT.parent / "holographic_cluster_state_ideal_gate_followup"
SOURCE_FOLLOWUP = STUDY_ROOT.parent / "holographic_cluster_state_followup"

for _path in (DATA_DIR, FIG_DIR, ARTIFACT_DIR, REPORT_DIR):
    _path.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Physical constants (2π × SI value)
# ---------------------------------------------------------------------------
TWO_PI = 2.0 * np.pi
OMEGA_Q = TWO_PI * 6.150e9      # transmon frequency (rad/s)
OMEGA_C = TWO_PI * 5.241e9      # cavity frequency (rad/s)
ALPHA   = TWO_PI * (-255.0e6)   # transmon anharmonicity (rad/s)
CHI     = TWO_PI * (-2.84e6)    # dispersive coupling (rad/s)
CHI_PRIME = TWO_PI * (-21.0e3)  # higher-order dispersive (rad/s)
KERR    = TWO_PI * (-28.0e3)    # cavity self-Kerr (rad/s)

# ---------------------------------------------------------------------------
# Truncation levels
# ---------------------------------------------------------------------------
GRAPE_N_CAV  = 8    # used for GRAPE benchmark
DECOMP_N_CAV = 4    # used for ideal-gate synthesis
N_TR         = 2    # transmon levels retained

# ---------------------------------------------------------------------------
# Logical basis labels and indices
# ---------------------------------------------------------------------------
LOGICAL_LABELS = ("|g,0>", "|g,1>", "|e,0>", "|e,1>")
LOGICAL_BASIS_ORDERING = tuple(LOGICAL_LABELS)
GROUND_SECTOR_LABELS = ("|g,0>", "|g,1>")

# ---------------------------------------------------------------------------
# Nominal gate durations (seconds)
# ---------------------------------------------------------------------------
DISPLACEMENT_S   = 80.0e-9          # Displacement gate
ROTATION_S       = 40.0e-9          # QubitRotation gate
SQR_S            = 1.10e-6          # Selective qubit rotation (SQR)
CPSQR_S          = 160.0e-9         # Conditional-phase SQR (CPSQR)
FE_DEFAULT_S     = np.pi / abs(CHI) # Default free-evolution (≈ 176 ns)

# ---------------------------------------------------------------------------
# GRAPE optimization settings
# ---------------------------------------------------------------------------
GRAPE_DT_S      = 4.0e-9                # Control resolution (4 ns)
GRAPE_AMP_BOUND = TWO_PI * 50.0e6       # Amplitude bound (rad/s)
GRAPE_DURATIONS_NS = (100, 150, 200, 250, 300, 400)
GRAPE_SEEDS = (17, 42, 73, 91, 103, 127, 211, 307, 401, 509)

# ---------------------------------------------------------------------------
# Noise model used for replay validation
# ---------------------------------------------------------------------------
QUBIT_T1_S = 30.0e-6
QUBIT_T2_S = 20.0e-6
QUBIT_TPHI_S = pure_dephasing_time_from_t1_t2(t1_s=QUBIT_T1_S, t2_s=QUBIT_T2_S)
CAVITY_T1_S = 250.0e-6

# ---------------------------------------------------------------------------
# Target unitary
# ---------------------------------------------------------------------------
# The cluster / transfer unitary equals SWAP · CZ · (H ⊗ I) in the
# logical basis {|g,0>, |g,1>, |e,0>, |e,1>}.
TARGET_UNITARY    = np.asarray(make_target("cluster", n_match=1), dtype=np.complex128)
TARGET_UNITARY_U2 = np.asarray(make_target("cluster", n_match=1, which="u2"), dtype=np.complex128)
TARGET = TargetUnitary(TARGET_UNITARY, ignore_global_phase=True)

# Ground-sector cavity-only target: |g> ⊗ |psi> -> |g> ⊗ H_c |psi> on
# {|g,0>, |g,1>}, where H_c is the logical cavity Hadamard.
GROUND_SECTOR_TARGET_UNITARY = np.asarray(
    [[1.0, 1.0], [1.0, -1.0]],
    dtype=np.complex128,
) / np.sqrt(2.0)

# ---------------------------------------------------------------------------
# Drift models
# ---------------------------------------------------------------------------
IDEAL_DRIFT    = DriftPhaseModel(chi=0.0, chi2=0.0, kerr=0.0)
PHYSICAL_DRIFT = DriftPhaseModel(chi=CHI, chi2=CHI_PRIME, kerr=KERR)
PHYSICAL_FE_DRIFT = DriftPhaseModel(chi=abs(CHI), chi2=CHI_PRIME, kerr=KERR)


# ---------------------------------------------------------------------------
# Utility: JSON serialization
# ---------------------------------------------------------------------------
def json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        if np.iscomplexobj(value):
            return {"real": value.real.tolist(), "imag": value.imag.tolist(), "shape": list(value.shape)}
        return value.tolist()
    if isinstance(value, complex):
        return {"real": float(value.real), "imag": float(value.imag)}
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return repr(value)


def save_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(json_ready(payload), indent=2), encoding="utf-8")


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _restore_timing_metadata(gate: Any, row: dict[str, Any]) -> Any:
    gate.optimize_time = bool(row.get("optimize_time", getattr(gate, "optimize_time", False)))
    time_bounds_raw = row.get("time_bounds")
    gate.time_bounds = None if time_bounds_raw is None else (float(time_bounds_raw[0]), float(time_bounds_raw[1]))
    gate.duration_ref = float(row.get("duration", getattr(gate, "duration", 0.0)))
    gate.time_group = row.get("time_group")
    gate.time_policy_locked = bool(row.get("time_policy_locked", getattr(gate, "time_policy_locked", False)))
    if "time_param_id" in row:
        gate.time_param_id = row.get("time_param_id")
    return gate


def sequence_from_payload(payload: Sequence[dict[str, Any]], *, n_cav: int) -> GateSequence:
    gates: list[Any] = []
    for row in payload:
        gate_type = str(row["type"])
        duration = float(row["duration"])
        name = str(row["name"])
        parameters = [float(value) for value in row.get("parameters", [])]
        metadata = dict(row.get("metadata", {})) if isinstance(row.get("metadata"), dict) else {}
        time_bounds_raw = row.get("time_bounds")
        time_bounds = None if time_bounds_raw is None else (float(time_bounds_raw[0]), float(time_bounds_raw[1]))
        common_kwargs = {
            "name": name,
            "duration": duration,
            "optimize_time": bool(row.get("optimize_time", False)),
            "time_bounds": time_bounds,
            "duration_ref": duration,
            "time_group": row.get("time_group"),
            "time_policy_locked": bool(row.get("time_policy_locked", False)),
        }
        if gate_type == "Displacement":
            gate = Displacement(alpha=complex(parameters[0], parameters[1]), **common_kwargs)
        elif gate_type == "QubitRotation":
            gate = QubitRotation(theta=float(parameters[0]), phi=float(parameters[1]), **common_kwargs)
        elif gate_type == "FreeEvolveCondPhase":
            gate = FreeEvolveCondPhase(drift_model=PHYSICAL_FE_DRIFT, **common_kwargs)
        elif gate_type == "PrimitiveGate":
            levels = tuple(int(level) for level in metadata.get("levels", []))
            ideal_kind = str(metadata.get("ideal_kind", ""))
            if ideal_kind == "MaskedSQR":
                tone_count = len(levels)
                gate = make_masked_sqr_gate(
                    name=name,
                    levels=levels,
                    n_cav=int(n_cav),
                    duration_s=duration,
                    include_conditional_phase=bool(metadata.get("include_conditional_phase", False)),
                )
                gate.parameters = {
                    "theta": np.asarray(parameters[:tone_count], dtype=float),
                    "phi": np.asarray(parameters[tone_count : 2 * tone_count], dtype=float),
                    "duration": duration,
                }
            elif ideal_kind == "MaskedCPSQR":
                gate = make_masked_cpsqr_gate(
                    name=name,
                    levels=levels,
                    n_cav=int(n_cav),
                    duration_s=duration,
                    include_drift=bool(metadata.get("include_drift", True)),
                )
                gate.parameters = {
                    "phases": np.asarray(parameters[: len(levels)], dtype=float),
                    "duration": duration,
                }
            else:
                raise ValueError(f"Unsupported PrimitiveGate payload kind '{ideal_kind}'.")
            gate.duration = duration
            _restore_timing_metadata(gate, row)
        else:
            raise ValueError(f"Unsupported gate payload type '{gate_type}'.")
        gates.append(gate)
    return GateSequence(gates=gates, n_cav=int(n_cav))


# ---------------------------------------------------------------------------
# Hilbert-space helpers
# ---------------------------------------------------------------------------
def logical_indices(n_cav: int) -> tuple[int, int, int, int]:
    """Return (i_g0, i_g1, i_e0, i_e1) within the full 2*n_cav Hilbert space."""
    return (0, 1, int(n_cav), int(n_cav) + 1)


def logical_subspace(n_cav: int) -> Subspace:
    return Subspace.custom(2 * int(n_cav), logical_indices(int(n_cav)), LOGICAL_LABELS)


def ground_sector_indices(n_cav: int) -> tuple[int, int]:
    """Return the {|g,0>, |g,1>} indices within the full 2*n_cav Hilbert space."""
    _ = int(n_cav)
    return (0, 1)


def ground_sector_subspace(n_cav: int) -> Subspace:
    return Subspace.custom(2 * int(n_cav), ground_sector_indices(int(n_cav)), GROUND_SECTOR_LABELS)


def ground_sector_input_states(n_cav: int) -> list[tuple[str, qt.Qobj]]:
    basis = [qt.basis(2 * int(n_cav), idx) for idx in ground_sector_indices(int(n_cav))]
    return list(zip(GROUND_SECTOR_LABELS, basis))


def level_subspace(levels: Sequence[int], *, n_cav: int) -> Subspace:
    ordered_levels = tuple(int(level) for level in levels)
    logical_levels = tuple(level for level in ordered_levels if level in (0, 1))
    spectator_levels = tuple(level for level in ordered_levels if level not in (0, 1))
    indices: list[int] = []
    labels: list[str] = []
    for level in logical_levels:
        indices.append(int(level))
        labels.append(f"|g,{int(level)}>")
    for level in logical_levels:
        indices.append(int(n_cav) + int(level))
        labels.append(f"|e,{int(level)}>")
    for level in spectator_levels:
        indices.append(int(level))
        labels.append(f"|g,{int(level)}>")
        indices.append(int(n_cav) + int(level))
        labels.append(f"|e,{int(level)}>")
    return Subspace.custom(2 * int(n_cav), tuple(indices), tuple(labels))


# ---------------------------------------------------------------------------
# Model and frame builders
# ---------------------------------------------------------------------------
def build_model(*, n_cav: int = GRAPE_N_CAV, n_tr: int = N_TR) -> DispersiveTransmonCavityModel:
    """Construct the dispersive transmon-cavity model with given truncations."""
    return DispersiveTransmonCavityModel(
        omega_q=OMEGA_Q,
        omega_c=OMEGA_C,
        alpha=ALPHA,
        chi=CHI,
        chi_higher=(CHI_PRIME,),
        kerr=KERR,
        n_cav=int(n_cav),
        n_tr=int(n_tr),
    )


def build_frame(model: DispersiveTransmonCavityModel | None = None) -> FrameSpec:
    if model is None:
        return FrameSpec(omega_q_frame=OMEGA_Q, omega_c_frame=OMEGA_C)
    return FrameSpec(omega_q_frame=float(model.omega_q), omega_c_frame=float(model.omega_c))


# ---------------------------------------------------------------------------
# Gate-time reference table
# ---------------------------------------------------------------------------
def gate_time_table_ns() -> dict[str, float]:
    return {
        "Displacement": DISPLACEMENT_S * 1.0e9,
        "QubitRotation": ROTATION_S * 1.0e9,
        "SQR": SQR_S * 1.0e9,
        "ConditionalPhaseSQR": CPSQR_S * 1.0e9,
        "FreeEvolveCondPhase_default": FE_DEFAULT_S * 1.0e9,
    }


# ---------------------------------------------------------------------------
# Target definition payload (for documentation / notebooks)
# ---------------------------------------------------------------------------
def target_definition_payload() -> dict[str, Any]:
    return {
        "logical_basis_ordering": LOGICAL_BASIS_ORDERING,
        "target_matrix_u1": TARGET_UNITARY,
        "target_matrix_u2": TARGET_UNITARY_U2,
        "ground_sector_basis_ordering": GROUND_SECTOR_LABELS,
        "ground_sector_target_matrix": GROUND_SECTOR_TARGET_UNITARY,
        "target_relation_note": (
            "The study optimizes make_target('cluster', n_match=1), "
            "which equals SWAP·CZ·(H⊗I). An equivalent convention is available "
            "via which='u2'."
        ),
        "ground_sector_relation_note": (
            "The cavity-only follow-up uses the restricted ground-sector transfer "
            "target |g>⊗|psi> -> |g>⊗H_c|psi> on {|g,0>, |g,1>}, where H_c is the "
            "logical cavity Hadamard."
        ),
    }


# ---------------------------------------------------------------------------
# Level-subset helpers (for SQR tone-budget sweep)
# ---------------------------------------------------------------------------
def ordered_level_subsets(n_cav: int, max_tones: int) -> list[tuple[int, ...]]:
    levels = range(int(n_cav))
    subsets = list(itertools.combinations(levels, int(max_tones)))
    subsets.sort(key=lambda row: (0 if 0 in row and 1 in row else 1, sum(row), row))
    return [tuple(int(level) for level in row) for row in subsets]


# ---------------------------------------------------------------------------
# Ideal-gate matrix constructors (used by PrimitiveGate wrappers)
# ---------------------------------------------------------------------------
def _identity_block_matrix(n_cav: int) -> np.ndarray:
    return np.eye(2 * int(n_cav), dtype=np.complex128)


def _masked_sqr_matrix(
    *,
    theta: Sequence[float],
    phi: Sequence[float],
    levels: Sequence[int],
    n_cav: int,
    duration: float,
    include_conditional_phase: bool,
    drift_model: DriftPhaseModel,
) -> np.ndarray:
    full = _identity_block_matrix(int(n_cav))
    for offset, level in enumerate(levels):
        block = np.asarray(qubit_rotation_xy(float(theta[offset]), float(phi[offset])).full(), dtype=np.complex128)
        idx = qubit_cavity_block_indices(int(n_cav), int(level))
        full[np.ix_(idx, idx)] = block
    if include_conditional_phase:
        drift = np.asarray(drift_phase_unitary(int(n_cav), float(duration), drift_model).full(), dtype=np.complex128)
        full = drift @ full
    return full


def _masked_cpsqr_matrix(
    *,
    phases: Sequence[float],
    levels: Sequence[int],
    n_cav: int,
    duration: float,
    include_drift: bool,
    drift_model: DriftPhaseModel,
) -> np.ndarray:
    full = _identity_block_matrix(int(n_cav))
    for offset, level in enumerate(levels):
        g_idx, e_idx = qubit_cavity_block_indices(int(n_cav), int(level))
        phase = float(phases[offset])
        full[g_idx, g_idx] = np.exp(-0.5j * phase)
        full[e_idx, e_idx] = np.exp(0.5j * phase)
    if include_drift:
        drift = np.asarray(drift_phase_unitary(int(n_cav), float(duration), drift_model).full(), dtype=np.complex128)
        full = drift @ full
    return full


# ---------------------------------------------------------------------------
# PrimitiveGate factory functions
# ---------------------------------------------------------------------------
def make_masked_sqr_gate(
    *,
    name: str,
    levels: Sequence[int],
    n_cav: int,
    duration_s: float = SQR_S,
    include_conditional_phase: bool = False,
    drift_model: DriftPhaseModel = IDEAL_DRIFT,
) -> PrimitiveGate:
    """Build an ideal Fock-selective qubit rotation on the specified levels."""
    level_tuple = tuple(int(level) for level in levels)
    theta_bounds: Any = [(-2.0 * np.pi, 2.0 * np.pi)] * len(level_tuple)
    phi_bounds: Any   = [(-np.pi, np.pi)] * len(level_tuple)
    if len(level_tuple) == 1:
        theta_bounds = (-2.0 * np.pi, 2.0 * np.pi)
        phi_bounds   = (-np.pi, np.pi)

    def matrix_fn(parameters: dict[str, Any], _model: Any | None = None) -> np.ndarray:
        return _masked_sqr_matrix(
            theta=parameters["theta"],
            phi=parameters["phi"],
            levels=level_tuple,
            n_cav=int(n_cav),
            duration=float(parameters["duration"]),
            include_conditional_phase=bool(include_conditional_phase),
            drift_model=drift_model,
        )

    return PrimitiveGate(
        name=name,
        duration=float(duration_s),
        optimize_time=False,
        matrix=matrix_fn,
        parameters={
            "theta": np.zeros(len(level_tuple), dtype=float),
            "phi":   np.zeros(len(level_tuple), dtype=float),
            "duration": float(duration_s),
        },
        parameter_bounds={"theta": theta_bounds, "phi": phi_bounds},
        hilbert_dim=2 * int(n_cav),
        metadata={
            "ideal_kind": "MaskedSQR",
            "levels": list(level_tuple),
            "max_tones": len(level_tuple),
            "include_conditional_phase": bool(include_conditional_phase),
        },
    )


def make_masked_cpsqr_gate(
    *,
    name: str,
    levels: Sequence[int],
    n_cav: int,
    duration_s: float = CPSQR_S,
    include_drift: bool = True,
    drift_model: DriftPhaseModel = PHYSICAL_DRIFT,
) -> PrimitiveGate:
    """Build an ideal Fock-selective conditional-phase gate on the specified levels."""
    level_tuple = tuple(int(level) for level in levels)
    phase_bounds: Any = [(-2.0 * np.pi, 2.0 * np.pi)] * len(level_tuple)
    if len(level_tuple) == 1:
        phase_bounds = (-2.0 * np.pi, 2.0 * np.pi)

    def matrix_fn(parameters: dict[str, Any], _model: Any | None = None) -> np.ndarray:
        return _masked_cpsqr_matrix(
            phases=parameters["phases"],
            levels=level_tuple,
            n_cav=int(n_cav),
            duration=float(parameters["duration"]),
            include_drift=bool(include_drift),
            drift_model=drift_model,
        )

    return PrimitiveGate(
        name=name,
        duration=float(duration_s),
        optimize_time=False,
        matrix=matrix_fn,
        parameters={"phases": np.zeros(len(level_tuple), dtype=float), "duration": float(duration_s)},
        parameter_bounds={"phases": phase_bounds},
        hilbert_dim=2 * int(n_cav),
        metadata={
            "ideal_kind": "MaskedCPSQR",
            "levels": list(level_tuple),
            "max_tones": len(level_tuple),
            "include_drift": bool(include_drift),
        },
    )


def displacement_gate(name: str) -> Displacement:
    return Displacement(name=name, alpha=0.10 + 0.0j, duration=DISPLACEMENT_S, optimize_time=False)


def rotation_gate(name: str, *, phi: float = 0.0) -> QubitRotation:
    return QubitRotation(name=name, theta=np.pi / 2.0, phi=float(phi), duration=ROTATION_S, optimize_time=False)


# ---------------------------------------------------------------------------
# Gate-sequence builders
# ---------------------------------------------------------------------------
def _normalize_order(order: Sequence[str] | str, *, expected: Sequence[str]) -> tuple[str, ...]:
    if isinstance(order, str):
        normalized = str(order).replace("-", "").replace("_", "").upper()
        tokens: list[str] = []
        idx = 0
        while idx < len(normalized):
            if normalized.startswith("CPSQR", idx):
                tokens.append("CPSQR")
                idx += len("CPSQR")
            elif normalized.startswith("SQR", idx):
                tokens.append("SQR")
                idx += len("SQR")
            else:
                tokens.append(normalized[idx])
                idx += 1
    else:
        tokens = [str(token).upper() for token in order]

    expected_tuple = tuple(str(token).upper() for token in expected)
    token_tuple = tuple(tokens)
    if sorted(token_tuple) != sorted(expected_tuple):
        raise ValueError(f"Order {token_tuple!r} does not match expected tokens {expected_tuple!r}.")
    return token_tuple


def _structured_rotation_phi(rotation_index: int) -> float:
    return 0.0 if int(rotation_index) % 2 == 0 else np.pi / 2.0


def _append_structured_gate(
    gates: list[Any],
    *,
    token: str,
    counters: dict[str, int],
    levels: Sequence[int],
    n_cav: int,
) -> None:
    normalized = str(token).upper()
    if normalized == "D":
        gates.append(displacement_gate(f"D{counters['D']}"))
    elif normalized == "R":
        gates.append(rotation_gate(f"R{counters['R']}", phi=_structured_rotation_phi(counters["R"])))
    elif normalized == "SQR":
        gates.append(make_masked_sqr_gate(name=f"S{counters['SQR']}", levels=levels, n_cav=n_cav))
    elif normalized == "CPSQR":
        gates.append(make_masked_cpsqr_gate(name=f"CP{counters['CPSQR']}", levels=levels, n_cav=n_cav))
    else:
        raise ValueError(f"Unknown structured gate token '{token}'.")
    counters[normalized] += 1


def build_ordered_sqr_sequence(
    *,
    levels: Sequence[int],
    n_cav: int,
    blocks: int,
    order: Sequence[str] | str = ("D", "R", "SQR"),
    close_with_final_displacement: bool = True,
) -> GateSequence:
    """Build a repeated ordered D/R/SQR block sequence (+ optional final D)."""
    token_order = _normalize_order(order, expected=("D", "R", "SQR"))
    counters = {"D": 0, "R": 0, "SQR": 0, "CPSQR": 0}
    gates: list[Any] = []
    for _ in range(int(blocks)):
        for token in token_order:
            _append_structured_gate(gates, token=token, counters=counters, levels=levels, n_cav=int(n_cav))
    if close_with_final_displacement:
        _append_structured_gate(gates, token="D", counters=counters, levels=levels, n_cav=int(n_cav))
    return GateSequence(gates=gates, n_cav=int(n_cav))


def build_ordered_cpsqr_sequence(
    *,
    levels: Sequence[int],
    n_cav: int,
    blocks: int,
    order: Sequence[str] | str = ("D", "R", "CPSQR"),
    close_with_final_displacement: bool = True,
) -> GateSequence:
    """Build a repeated ordered D/R/CPSQR block sequence (+ optional final D)."""
    token_order = _normalize_order(order, expected=("D", "R", "CPSQR"))
    counters = {"D": 0, "R": 0, "SQR": 0, "CPSQR": 0}
    gates: list[Any] = []
    for _ in range(int(blocks)):
        for token in token_order:
            _append_structured_gate(gates, token=token, counters=counters, levels=levels, n_cav=int(n_cav))
    if close_with_final_displacement:
        _append_structured_gate(gates, token="D", counters=counters, levels=levels, n_cav=int(n_cav))
    return GateSequence(gates=gates, n_cav=int(n_cav))


def build_drsqr_3sqr_sequence(*, levels: Sequence[int], n_cav: int, pattern: str = "drs") -> GateSequence:
    """Build a D-R-SQR repeated three times (+ final D), the 3-SQR ansatz."""
    if pattern.lower() == "drs":
        order = ("D", "R", "SQR")
    elif pattern.lower() == "dsr":
        order = ("D", "SQR", "R")
    else:
        raise ValueError(f"Unknown 3-SQR pattern '{pattern}'.")
    return build_ordered_sqr_sequence(levels=levels, n_cav=int(n_cav), blocks=3, order=order)


def build_drsqr_4sqr_sequence(*, levels: Sequence[int], n_cav: int) -> GateSequence:
    """Build a D-R-SQR repeated four times (+ final D), the 4-SQR ansatz."""
    return build_ordered_sqr_sequence(levels=levels, n_cav=int(n_cav), blocks=4, order=("D", "R", "SQR"))


def build_drcpsqr_sequence(*, levels: Sequence[int], n_cav: int, blocks: int = 3) -> GateSequence:
    """Build a D-R-CPSQR repeated `blocks` times (+ final D), the pure CPSQR ansatz."""
    return build_ordered_cpsqr_sequence(levels=levels, n_cav=int(n_cav), blocks=int(blocks), order=("D", "R", "CPSQR"))


def build_dsqrcpsqr_sequence(*, levels: Sequence[int], n_cav: int, blocks: int = 2) -> GateSequence:
    """Build D-(SQR-CPSQR-SQR) repeated `blocks` times (+ final D)."""
    gates: list[Any] = []
    for block in range(int(blocks)):
        gates.extend([
            displacement_gate(f"D{block}"),
            make_masked_sqr_gate(name=f"S{2 * block}", levels=levels, n_cav=n_cav),
            make_masked_cpsqr_gate(name=f"CP{block}", levels=levels, n_cav=n_cav),
            make_masked_sqr_gate(name=f"S{2 * block + 1}", levels=levels, n_cav=n_cav),
        ])
    gates.append(displacement_gate(f"D{blocks}"))
    return GateSequence(gates=gates, n_cav=int(n_cav))


def build_drsqrcpsqr_sequence(*, levels: Sequence[int], n_cav: int, blocks: int = 2) -> GateSequence:
    """Build D-R-SQR-CPSQR-R repeated `blocks` times (+ final D)."""
    gates: list[Any] = []
    for block in range(int(blocks)):
        gates.extend([
            displacement_gate(f"D{block}"),
            rotation_gate(f"R{2 * block}", phi=0.0),
            make_masked_sqr_gate(name=f"S{block}", levels=levels, n_cav=n_cav),
            make_masked_cpsqr_gate(name=f"CP{block}", levels=levels, n_cav=n_cav),
            rotation_gate(f"R{2 * block + 1}", phi=np.pi / 2.0),
        ])
    gates.append(displacement_gate(f"D{blocks}"))
    return GateSequence(gates=gates, n_cav=int(n_cav))


def build_drfe_sequence(*, n_cav: int, blocks: int = 2, fe_duration_s: float = FE_DEFAULT_S) -> GateSequence:
    """Build D-R-FreeEvolve repeated `blocks` times (+ final D+R)."""
    gates: list[Any] = []
    for block in range(int(blocks)):
        gates.extend([
            displacement_gate(f"D{block}"),
            rotation_gate(f"R{block}", phi=0.0 if block == 0 else np.pi / 2.0),
            FreeEvolveCondPhase(
                name=f"FE{block}",
                duration=float(fe_duration_s),
                optimize_time=True,
                time_bounds=(40.0e-9, 400.0e-9),
                duration_ref=float(fe_duration_s),
                drift_model=PHYSICAL_FE_DRIFT,
            ),
        ])
    gates.extend([displacement_gate(f"D{blocks}"), rotation_gate(f"R{blocks}", phi=np.pi / 2.0)])
    return GateSequence(gates=gates, n_cav=int(n_cav))


# ---------------------------------------------------------------------------
# Sequence analysis
# ---------------------------------------------------------------------------
def gate_kind(gate: Any) -> str:
    if isinstance(gate, PrimitiveGate):
        return str(gate.metadata.get("ideal_kind", "PrimitiveGate"))
    return type(gate).__name__


def sequence_summary(sequence: GateSequence) -> dict[str, Any]:
    """Return a structured summary of a GateSequence."""
    rows: list[dict[str, Any]] = []
    counts: dict[str, int] = {}
    total_ns = 0.0
    selective_blocks = entangling_blocks = total_tones = max_tones = 0
    for gate in sequence.gates:
        kind = gate_kind(gate)
        counts[kind] = counts.get(kind, 0) + 1
        tones = int(gate.metadata.get("max_tones", 0)) if isinstance(gate, PrimitiveGate) else 0
        total_tones += tones
        max_tones = max(max_tones, tones)
        if kind in {"MaskedSQR", "MaskedCPSQR"}:
            selective_blocks += 1
        if kind in {"MaskedCPSQR", "FreeEvolveCondPhase"}:
            entangling_blocks += 1
        total_ns += float(gate.duration) * 1.0e9
        rows.append({"name": str(gate.name), "kind": kind, "duration_ns": float(gate.duration) * 1.0e9,
                     "active_tones": tones, "metadata": getattr(gate, "metadata", {})})
    return {
        "rows": rows, "counts": counts, "gate_depth": len(sequence.gates),
        "total_duration_ns": float(total_ns), "selective_block_count": int(selective_blocks),
        "entangling_block_count": int(entangling_blocks), "total_active_tones": int(total_tones),
        "max_active_tones": int(max_tones),
    }


# ---------------------------------------------------------------------------
# Synthesis and evaluation helpers
# ---------------------------------------------------------------------------
def fit_sequence(
    sequence: GateSequence,
    *,
    n_cav: int,
    seed: int,
    init_guess: str,
    multistart: int,
    maxiter: int,
    duration_weight: float = 0.0,
    gate_count_weight: float = 0.0,
    use_fast_path: bool = True,
    warm_start: Any | None = None,
    target_unitary: np.ndarray | None = None,
    subspace: Subspace | None = None,
) -> dict[str, Any]:
    """Optimize a GateSequence against the chosen target using UnitarySynthesizer."""
    target_matrix = TARGET_UNITARY if target_unitary is None else np.asarray(target_unitary, dtype=np.complex128)
    resolved_subspace = logical_subspace(int(n_cav)) if subspace is None else subspace
    target = TargetUnitary(target_matrix, ignore_global_phase=True)
    synthesizer = UnitarySynthesizer(
        primitives=copy.deepcopy(sequence.gates),
        subspace=resolved_subspace,
        objectives=MultiObjective(
            fidelity_weight=1.0,
            leakage_weight=0.05,
            duration_weight=float(duration_weight),
            gate_count_weight=float(gate_count_weight),
        ),
        leakage_penalty=LeakagePenalty(weight=0.05),
        execution=ExecutionOptions(engine="auto", use_fast_path=bool(use_fast_path)),
        warm_start=warm_start,
        seed=int(seed),
    )
    result = synthesizer.fit(target=target, init_guess=str(init_guess),
                             multistart=int(multistart), maxiter=int(maxiter))
    fidelity = float(subspace_unitary_fidelity(
        np.asarray(result.simulation.subspace_operator, dtype=np.complex128),
        target_matrix, gauge="global"))
    summary = sequence_summary(result.sequence)
    return {
        "result": result, "fidelity": fidelity, "summary": summary,
        "objective": float(getattr(result, "objective", np.nan)),
        "success": bool(getattr(result, "success", True)),
        "message": str(getattr(result, "message", "")),
        "metrics": dict(result.report.get("metrics", {})),
        "sequence_payload": result.sequence.serialize(),
        "target_matrix": target_matrix,
    }


def evaluate_sequence(
    sequence: GateSequence,
    *,
    n_cav: int,
    target_unitary: np.ndarray | None = None,
    subspace: Subspace | None = None,
) -> dict[str, Any]:
    """Simulate a GateSequence and return fidelity / leakage metrics."""
    target_matrix = TARGET_UNITARY if target_unitary is None else np.asarray(target_unitary, dtype=np.complex128)
    resolved_subspace = logical_subspace(int(n_cav)) if subspace is None else subspace
    simulation = synth_simulate_sequence(
        sequence, subspace=resolved_subspace,
        backend="ideal", target_subspace=target_matrix,
        leakage_weight=0.05, gauge="global", need_operator=True,
    )
    subspace_operator = np.asarray(simulation.subspace_operator, dtype=np.complex128)
    full_operator     = np.asarray(simulation.full_operator,     dtype=np.complex128)
    fidelity = float(subspace_unitary_fidelity(subspace_operator, target_matrix, gauge="global"))
    metrics  = dict(simulation.metrics)
    return {
        "simulation": simulation, "fidelity": fidelity,
        "subspace_operator": subspace_operator, "full_operator": full_operator,
        "metrics": metrics,
        "leakage_average": float(metrics.get("leakage_average", np.nan)),
        "leakage_worst":   float(metrics.get("leakage_worst",   np.nan)),
        "unitarity_error": float(metrics.get("unitarity_error", np.nan)),
    }


def _restrict_operator_to_indices(operator: np.ndarray, indices: Sequence[int]) -> np.ndarray:
    matrix = np.asarray(operator, dtype=np.complex128)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("Operator must be a square matrix.")
    index_array = np.asarray(tuple(int(index) for index in indices), dtype=int)
    if matrix.shape == (index_array.size, index_array.size):
        return matrix
    if np.any(index_array < 0) or np.any(index_array >= matrix.shape[0]):
        raise ValueError("Requested restriction indices fall outside the operator dimensions.")
    return matrix[np.ix_(index_array, index_array)]


def restricted_ground_sector_fidelity(
    operator: np.ndarray,
    *,
    n_cav: int,
    target_unitary: np.ndarray | None = None,
) -> float:
    """Return the process fidelity on {|g,0>, |g,1>} for the cavity-only target."""
    _ = int(n_cav)
    target_matrix = (
        GROUND_SECTOR_TARGET_UNITARY
        if target_unitary is None
        else np.asarray(target_unitary, dtype=np.complex128)
    )
    restricted = _restrict_operator_to_indices(operator, ground_sector_indices(int(n_cav)))
    if restricted.shape != target_matrix.shape:
        raise ValueError(
            f"Ground-sector operator shape {restricted.shape} does not match target shape {target_matrix.shape}."
        )
    return float(subspace_unitary_fidelity(restricted, target_matrix, gauge="global"))


def ground_sector_leakage(
    operator: np.ndarray,
    *,
    n_cav: int,
    target_levels: Sequence[int] = (0, 1),
) -> dict[str, Any]:
    """Decompose leakage from {|g,0>, |g,1>} into ancilla excitation and cavity-support spillover."""
    full_dim = 2 * int(n_cav)
    matrix = np.asarray(operator, dtype=np.complex128)
    if matrix.shape != (full_dim, full_dim):
        raise ValueError(
            f"Ground-sector leakage expects a full ({full_dim}, {full_dim}) operator; received {matrix.shape}."
        )

    support_levels = tuple(int(level) for level in target_levels)
    if any(level < 0 or level >= int(n_cav) for level in support_levels):
        raise ValueError(f"Target levels {support_levels} must lie in [0, {int(n_cav) - 1}].")
    support_indices = np.asarray(support_levels, dtype=int)
    ground_slice = slice(0, int(n_cav))
    excited_slice = slice(int(n_cav), full_dim)

    per_input: list[dict[str, Any]] = []
    ancilla_values: list[float] = []
    support_values: list[float] = []
    outside_values: list[float] = []

    for label, input_index in zip(GROUND_SECTOR_LABELS, ground_sector_indices(int(n_cav))):
        column = matrix[:, int(input_index)]
        populations = np.abs(column) ** 2
        ancilla_excitation = float(np.clip(np.sum(populations[excited_slice]), 0.0, 1.0))
        ground_population = float(np.clip(np.sum(populations[ground_slice]), 0.0, 1.0))
        target_support_population = float(np.clip(np.sum(populations[support_indices]), 0.0, 1.0))
        support_leakage = float(np.clip(ground_population - target_support_population, 0.0, 1.0))
        outside_target_leakage = float(np.clip(ancilla_excitation + support_leakage, 0.0, 1.0))
        ancilla_values.append(ancilla_excitation)
        support_values.append(support_leakage)
        outside_values.append(outside_target_leakage)
        per_input.append(
            {
                "label": label,
                "input_index": int(input_index),
                "ancilla_excitation": ancilla_excitation,
                "ground_population": ground_population,
                "target_support_population": target_support_population,
                "support_leakage": support_leakage,
                "outside_target_leakage": outside_target_leakage,
                "population_norm": float(np.sum(populations)),
            }
        )

    return {
        "target_levels": [int(level) for level in support_levels],
        "per_input": per_input,
        "ancilla_excitation_average": float(np.mean(ancilla_values)) if ancilla_values else float("nan"),
        "ancilla_excitation_worst": float(np.max(ancilla_values)) if ancilla_values else float("nan"),
        "support_leakage_average": float(np.mean(support_values)) if support_values else float("nan"),
        "support_leakage_worst": float(np.max(support_values)) if support_values else float("nan"),
        "outside_target_leakage_average": float(np.mean(outside_values)) if outside_values else float("nan"),
        "outside_target_leakage_worst": float(np.max(outside_values)) if outside_values else float("nan"),
    }


def evaluate_ground_sector_transfer(
    sequence: GateSequence,
    *,
    n_cav: int,
    target_unitary: np.ndarray | None = None,
) -> dict[str, Any]:
    """Evaluate a sequence on {|g,0>, |g,1>} against the cavity-only target and report leakage diagnostics."""
    target_matrix = (
        GROUND_SECTOR_TARGET_UNITARY
        if target_unitary is None
        else np.asarray(target_unitary, dtype=np.complex128)
    )
    evaluation = evaluate_sequence(
        sequence,
        n_cav=int(n_cav),
        target_unitary=target_matrix,
        subspace=ground_sector_subspace(int(n_cav)),
    )
    leakage = ground_sector_leakage(evaluation["full_operator"], n_cav=int(n_cav))
    return {
        "simulation": evaluation["simulation"],
        "fidelity": float(evaluation["fidelity"]),
        "restricted_fidelity": float(evaluation["fidelity"]),
        "target_matrix": target_matrix,
        "subspace_operator": evaluation["subspace_operator"],
        "full_operator": evaluation["full_operator"],
        "metrics": dict(evaluation["metrics"]),
        "subspace_leakage_average": float(evaluation["leakage_average"]),
        "subspace_leakage_worst": float(evaluation["leakage_worst"]),
        "unitarity_error": float(evaluation["unitarity_error"]),
        "ancilla_excitation_average": float(leakage["ancilla_excitation_average"]),
        "ancilla_excitation_worst": float(leakage["ancilla_excitation_worst"]),
        "support_leakage_average": float(leakage["support_leakage_average"]),
        "support_leakage_worst": float(leakage["support_leakage_worst"]),
        "outside_target_leakage_average": float(leakage["outside_target_leakage_average"]),
        "outside_target_leakage_worst": float(leakage["outside_target_leakage_worst"]),
        "ground_sector_leakage": leakage,
    }


def best_result_from_trials(trials: Sequence[dict[str, Any]]) -> dict[str, Any]:
    return max(trials, key=lambda r: (float(r["fidelity"]),
                                      -float(r["summary"]["total_duration_ns"]),
                                      -float(r["summary"]["gate_depth"])))


# ---------------------------------------------------------------------------
# GRAPE helpers
# ---------------------------------------------------------------------------
def jax_available() -> bool:
    try:
        importlib.import_module("jax")
    except Exception:
        return False
    return True


def resolve_grape_engine(engine: str = "auto") -> str:
    requested = str(engine).strip().lower()
    if requested not in {"auto", "numpy", "jax"}:
        raise ValueError(f"Unsupported GRAPE engine '{engine}'. Expected one of auto, numpy, jax.")
    if requested == "auto":
        return "jax" if jax_available() else "numpy"
    if requested == "jax" and not jax_available():
        raise RuntimeError("GRAPE engine 'jax' was requested, but the jax package is not available.")
    return requested


def build_grape_problem(
    *,
    model: DispersiveTransmonCavityModel,
    duration_ns: float,
    amp_bound_rad_s: float = GRAPE_AMP_BOUND,
) -> Any:
    """Build a GRAPE control problem for the cluster unitary."""
    frame    = build_frame(model)
    duration_s = float(duration_ns) * 1.0e-9
    steps    = max(10, round(duration_s / GRAPE_DT_S))
    return build_control_problem_from_model(
        model,
        frame=frame,
        time_grid=PiecewiseConstantTimeGrid.uniform(steps=int(steps), dt_s=duration_s / int(steps)),
        channel_specs=(
            ModelControlChannelSpec(name="storage", target="storage", quadratures=("I", "Q"),
                                    amplitude_bounds=(-float(amp_bound_rad_s), float(amp_bound_rad_s)),
                                    export_channel="storage"),
            ModelControlChannelSpec(name="qubit", target="qubit", quadratures=("I", "Q"),
                                    amplitude_bounds=(-float(amp_bound_rad_s), float(amp_bound_rad_s)),
                                    export_channel="qubit"),
        ),
        objectives=(
            OCUnitaryObjective(
                target_operator=TARGET_UNITARY,
                subspace=logical_subspace(int(model.n_cav)),
                ignore_global_phase=True,
                name=f"cluster_{int(round(duration_ns))}ns",
            ),
        ),
        penalties=(OCLeakagePenalty(weight=0.02, subspace=logical_subspace(int(model.n_cav))),),
    )


def run_grape_seed(
    problem: Any,
    *,
    seed: int,
    maxiter: int,
    engine: str = "auto",
    jax_device: str | None = None,
) -> Any:
    resolved_engine = resolve_grape_engine(engine)
    config_kwargs: dict[str, Any] = {
        "maxiter": int(maxiter),
        "seed": int(seed),
        "random_scale": 0.30,
        "history_every": 1,
        "engine": resolved_engine,
    }
    if resolved_engine == "jax" and jax_device:
        config_kwargs["jax_device"] = str(jax_device)
    solver = GrapeSolver(GrapeConfig(**config_kwargs))
    return solver.solve(problem)


def summarise_grape_result(result: Any, *, n_cav: int) -> dict[str, Any]:
    metrics = dict(getattr(result, "metrics", {}))
    return {
        "nominal_fidelity": float(metrics.get("nominal_fidelity", metrics.get("fidelity", np.nan))),
        "message": str(getattr(result, "message", "")),
        "success": bool(getattr(result, "success", True)),
        "iterations": int(getattr(result, "nit", 0)),
        "metrics": metrics,
    }


def logical_basis_states(model: DispersiveTransmonCavityModel) -> list[qt.Qobj]:
    return [
        model.basis_state(0, 0),
        model.basis_state(0, 1),
        model.basis_state(1, 0),
        model.basis_state(1, 1),
    ]


def default_noise_spec() -> NoiseSpec:
    return NoiseSpec(t1=QUBIT_T1_S, tphi=QUBIT_TPHI_S, kappa=1.0 / CAVITY_T1_S)


def replay_compiled_sequence(
    *,
    model: DispersiveTransmonCavityModel,
    compiled: Any,
    drive_ops: dict[str, Any],
    basis_states: Sequence[tuple[str, qt.Qobj]],
    noise: NoiseSpec | None = None,
    store_states: bool = False,
) -> list[dict[str, Any]]:
    frame = build_frame(model)
    rows: list[dict[str, Any]] = []
    for label, initial_state in basis_states:
        result = pulse_simulate_sequence(
            model,
            compiled,
            initial_state,
            drive_ops,
            config=SimulationConfig(frame=frame, store_states=bool(store_states)),
            noise=noise,
        )
        rows.append({"label": str(label), "simulation": result})
    return rows


def replay_grape_subspace(
    *,
    result: Any,
    problem: Any,
    model: DispersiveTransmonCavityModel,
    subspace: Subspace,
    target_unitary: np.ndarray,
    basis_states: Sequence[tuple[str, qt.Qobj]],
    noise: NoiseSpec | None = None,
    store_states: bool = False,
) -> dict[str, Any]:
    pulses, drive_ops, pulse_meta = result.to_pulses()
    compiler = SequenceCompiler(dt=1.0e-9)
    compiled = compiler.compile(pulses, t_end=problem.time_grid.duration_s)
    rows = replay_compiled_sequence(
        model=model,
        compiled=compiled,
        drive_ops=drive_ops,
        basis_states=basis_states,
        noise=noise,
        store_states=store_states,
    )
    indices = np.asarray(subspace.indices, dtype=int)
    sub_operator = np.zeros((indices.size, indices.size), dtype=np.complex128)
    leakage_values: list[float] = []
    for column, row in enumerate(rows):
        state = row["simulation"].final_state
        if noise is None:
            vector = np.asarray(state.full(), dtype=np.complex128).reshape(-1)
            logical = vector[indices]
            sub_operator[:, column] = logical
            leakage_values.append(float(max(0.0, 1.0 - np.vdot(logical, logical).real)))
    payload = {
        "pulses": pulses,
        "drive_ops": drive_ops,
        "pulse_meta": pulse_meta,
        "compiled": compiled,
        "rows": rows,
    }
    if noise is None:
        payload.update(
            {
                "subspace_operator": sub_operator,
                "fidelity": float(subspace_unitary_fidelity(sub_operator, np.asarray(target_unitary, dtype=np.complex128), gauge="global")),
                "leakage_average": float(np.mean(leakage_values)) if leakage_values else float("nan"),
                "leakage_worst": float(np.max(leakage_values)) if leakage_values else float("nan"),
            }
        )
    return payload
