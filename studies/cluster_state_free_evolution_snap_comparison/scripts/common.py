from __future__ import annotations

import copy
import sys
from pathlib import Path
from typing import Any

import numpy as np

WORKSPACE_ROOT = Path(__file__).resolve().parents[3]
SIM_ROOT_CANDIDATES = (
    Path("C:/Users/jl82323/Box/Shyam Shankar Quantum Circuits Group/Users/Users_JianJun/cQED_simulation"),
    Path("C:/Users/dazzl/Box/Shyam Shankar Quantum Circuits Group/Users/Users_JianJun/cQED_simulation"),
)
SIM_ROOT = next((path for path in SIM_ROOT_CANDIDATES if (path / "cqed_sim").exists()), SIM_ROOT_CANDIDATES[0])
for _path in (WORKSPACE_ROOT, SIM_ROOT):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from cqed_sim.unitary_synthesis import Displacement, FreeEvolveCondPhase, GateSequence, QubitRotation, SNAP
from cqed_sim.unitary_synthesis.metrics import logical_block_phase_diagnostics

from studies.cluster_state_holographic_unified.scripts import common as base


STUDY_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = STUDY_ROOT / "data"
FIG_DIR = STUDY_ROOT / "figures"
ARTIFACT_DIR = STUDY_ROOT / "artifacts"
REPORT_DIR = STUDY_ROOT / "report"

for _path in (DATA_DIR, FIG_DIR, ARTIFACT_DIR, REPORT_DIR):
    _path.mkdir(parents=True, exist_ok=True)

SNAP_S = 200.0e-9
SCREEN_SEED = 17
REFINE_SEEDS = (17, 42)
SCREEN_MAXITER = 20
REFINE_MAXITER = 35
FRONTIER_WEIGHTS = (0.0, 0.0025, 0.005, 0.01, 0.02)
TRUNCATION_LEVELS = (4, 6, 8)

FAMILY_SPECS: dict[str, dict[str, Any]] = {
    "native_fe": {
        "title": "D + R + FE",
        "placement": "native",
        "color": "#4477AA",
        "blocks": (1, 2, 3, 4, 5, 6),
    },
    "snap_tail": {
        "title": "D + R + FE + tail SNAP",
        "placement": "tail",
        "color": "#EE6677",
        "blocks": (1, 2, 3, 4),
    },
    "snap_interleaved": {
        "title": "D + R + FE + interleaved SNAP",
        "placement": "interleaved",
        "color": "#228833",
        "blocks": (1, 2, 3, 4, 5, 6),
    },
    "snap_hybrid": {
        "title": "D + R + FE + interleaved/tail SNAP",
        "placement": "hybrid",
        "color": "#CCBB44",
        "blocks": (1, 2, 3, 4),
    },
}

CAVITY_BLOCK_SLICES = ([0, 2], [1, 3])


def case_id(family: str, blocks: int) -> str:
    return f"{family}_b{int(blocks)}"


def snap_gate(name: str, *, n_cav: int, duration_s: float = SNAP_S) -> SNAP:
    return SNAP(name=name, phases=[0.0] * int(n_cav), duration=float(duration_s), optimize_time=False)


def build_family_sequence(
    *,
    family: str,
    blocks: int,
    n_cav: int = base.DECOMP_N_CAV,
    fe_duration_s: float = base.FE_DEFAULT_S,
    snap_duration_s: float = SNAP_S,
) -> GateSequence:
    placement = FAMILY_SPECS[str(family)]["placement"]
    if placement == "native":
        return base.build_drfe_sequence(n_cav=int(n_cav), blocks=int(blocks), fe_duration_s=float(fe_duration_s))

    gates: list[Any] = []
    for block in range(int(blocks)):
        gates.append(base.displacement_gate(f"D{block}"))
        gates.append(base.rotation_gate(f"R{block}", phi=0.0 if block == 0 else np.pi / 2.0))
        gates.append(
            FreeEvolveCondPhase(
                name=f"FE{block}",
                duration=float(fe_duration_s),
                optimize_time=True,
                time_bounds=(40.0e-9, 400.0e-9),
                duration_ref=float(fe_duration_s),
                drift_model=base.PHYSICAL_FE_DRIFT,
            )
        )
        if placement in {"interleaved", "hybrid"}:
            gates.append(snap_gate(f"S{block}", n_cav=int(n_cav), duration_s=float(snap_duration_s)))
    gates.append(base.displacement_gate(f"D{blocks}"))
    if placement in {"tail", "hybrid"}:
        gates.append(snap_gate("S_tail", n_cav=int(n_cav), duration_s=float(snap_duration_s)))
    gates.append(base.rotation_gate(f"R{blocks}", phi=np.pi / 2.0))
    return GateSequence(gates=gates, n_cav=int(n_cav))


def sequence_for_n_cav(sequence: GateSequence, *, n_cav: int) -> GateSequence:
    return GateSequence(gates=copy.deepcopy(sequence.gates), n_cav=int(n_cav))


def sequence_from_payload(payload: list[dict[str, Any]], *, n_cav: int = base.DECOMP_N_CAV) -> GateSequence:
    gates: list[Any] = []
    for row in payload:
        gate_type = str(row["type"])
        params = [float(value) for value in row.get("parameters", [])]
        time_bounds_raw = row.get("time_bounds")
        time_bounds = None if time_bounds_raw is None else (float(time_bounds_raw[0]), float(time_bounds_raw[1]))
        common_kwargs = {
            "name": str(row["name"]),
            "duration": float(row["duration"]),
            "optimize_time": bool(row.get("optimize_time", False)),
            "time_bounds": time_bounds,
            "duration_ref": float(row["duration"]),
            "time_group": row.get("time_group"),
            "time_policy_locked": bool(row.get("time_policy_locked", False)),
        }
        if gate_type == "Displacement":
            gates.append(Displacement(alpha=complex(params[0], params[1]), **common_kwargs))
        elif gate_type == "QubitRotation":
            gates.append(QubitRotation(theta=float(params[0]), phi=float(params[1]), **common_kwargs))
        elif gate_type == "FreeEvolveCondPhase":
            gates.append(FreeEvolveCondPhase(drift_model=base.PHYSICAL_FE_DRIFT, **common_kwargs))
        elif gate_type == "SNAP":
            gates.append(
                SNAP(
                    phases=params,
                    fock_levels=tuple(int(level) for level in row.get("fock_levels", [])),
                    **common_kwargs,
                )
            )
        else:
            raise ValueError(f"Unsupported gate type in local payload reconstruction: {gate_type}")
    return GateSequence(gates=gates, n_cav=int(n_cav))


def gate_parameter_count(gate: Any, *, n_cav: int) -> int:
    count = len(gate.parameter_names(int(n_cav)))
    if bool(getattr(gate, "optimize_time", False)):
        count += 1
    return int(count)


def sequence_complexity_metrics(sequence: GateSequence) -> dict[str, Any]:
    parameter_count = 0
    snap_gate_count = 0
    snap_phase_count = 0
    entangling_gate_count = 0
    wait_gate_count = 0
    for gate in sequence.gates:
        parameter_count += gate_parameter_count(gate, n_cav=int(sequence.n_cav))
        if isinstance(gate, SNAP):
            snap_gate_count += 1
            snap_phase_count += len(gate.parameter_names(int(sequence.n_cav)))
        if isinstance(gate, FreeEvolveCondPhase):
            entangling_gate_count += 1
            wait_gate_count += 1
    return {
        "parameter_count": int(parameter_count),
        "snap_gate_count": int(snap_gate_count),
        "snap_phase_count": int(snap_phase_count),
        "entangling_gate_count": int(entangling_gate_count),
        "wait_gate_count": int(wait_gate_count),
    }


def sequence_phase_budget(sequence: GateSequence, *, n_match: int = 1) -> dict[str, Any]:
    rows = sequence.phase_decomposition(n_match=int(n_match))
    fe_rows = [row for row in rows if row.get("gate_type") == "FreeEvolveCondPhase"]
    snap_rows = [row for row in rows if row.get("gate_type") == "SNAP"]

    total_wait_time_ns = float(sum(float(row.get("wait_time", row.get("duration", 0.0))) for row in fe_rows) * 1.0e9)
    total_fe_logical_delta_phi_rad = float(
        sum(abs(float(row.get("delta_phi", [0.0, 0.0])[1])) for row in fe_rows if len(row.get("delta_phi", [])) >= 2)
    )
    total_snap_logical_relative_phase_rad = float(
        sum(
            abs(float(row.get("relative_phases_rad", [0.0, 0.0])[1]))
            for row in snap_rows
            if len(row.get("relative_phases_rad", [])) >= 2
        )
    )
    return {
        "phase_rows": rows,
        "fe_rows": fe_rows,
        "snap_rows": snap_rows,
        "total_wait_time_ns": total_wait_time_ns,
        "total_fe_logical_delta_phi_rad": total_fe_logical_delta_phi_rad,
        "total_snap_logical_relative_phase_rad": total_snap_logical_relative_phase_rad,
    }


def evaluate_sequence_with_diagnostics(sequence: GateSequence, *, n_cav: int) -> dict[str, Any]:
    evaluation = base.evaluate_sequence(sequence, n_cav=int(n_cav))
    diag = logical_block_phase_diagnostics(
        evaluation["subspace_operator"],
        base.TARGET_UNITARY,
        block_slices=CAVITY_BLOCK_SLICES,
    )
    return {
        **evaluation,
        "block_gauge_fidelity": float(diag.block_gauge_fidelity),
        "best_fit_block_gauge_fidelity": float(diag.best_fit_block_gauge_fidelity),
        "block_phases_rad": [float(value) for value in diag.block_phases_rad],
        "best_fit_correction_phases_rad": [float(value) for value in diag.best_fit_correction_phases_rad],
        "rms_block_phase_error_rad": float(diag.rms_block_phase_error_rad),
    }


def remove_gate_types(sequence: GateSequence, gate_types: tuple[type[Any], ...]) -> GateSequence:
    gates = [copy.deepcopy(gate) for gate in sequence.gates if not isinstance(gate, gate_types)]
    return GateSequence(gates=gates, n_cav=int(sequence.n_cav))


def ablation_sequences(sequence: GateSequence) -> dict[str, GateSequence]:
    return {
        "full": sequence_for_n_cav(sequence, n_cav=int(sequence.n_cav)),
        "without_snap": remove_gate_types(sequence, (SNAP,)),
        "without_free_evolution": remove_gate_types(sequence, (FreeEvolveCondPhase,)),
    }


def summarize_sequence(sequence: GateSequence) -> dict[str, Any]:
    return {
        **base.sequence_summary(sequence),
        **sequence_complexity_metrics(sequence),
        **{k: v for k, v in sequence_phase_budget(sequence).items() if k != "phase_rows" and k != "fe_rows" and k != "snap_rows"},
    }
