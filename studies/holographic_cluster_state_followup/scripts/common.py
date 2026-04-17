from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Any

import numpy as np
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
from cqed_sim.optimal_control import (
    GrapeConfig,
    GrapeSolver,
    LeakagePenalty as OCLeakagePenalty,
    UnitaryObjective as OCUnitaryObjective,
    build_control_problem_from_model,
)
from cqed_sim.quantum_algorithms import HolographicChannel
from cqed_sim.sim.extractors import cavity_wigner
from cqed_sim.sim.noise import NoiseSpec, pure_dephasing_time_from_t1_t2
from cqed_sim.unitary_synthesis import (
    CQEDSystemAdapter,
    ConditionalPhaseSQR,
    Displacement,
    DriftPhaseModel,
    ExecutionOptions,
    FreeEvolveCondPhase,
    GateSequence,
    LeakagePenalty,
    MultiObjective,
    QubitRotation,
    SNAP,
    SQR,
    Subspace,
    TargetUnitary,
    UnitarySynthesizer,
    SynthesisConstraints,
    leakage_metrics,
    simulate_sequence as synth_simulate_sequence,
    subspace_unitary_fidelity,
)
from cqed_sim.unitary_synthesis.targets import make_target
from cqed_sim.unitary_synthesis.waveform_bridge import waveform_sequence_from_gates


STUDY_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = STUDY_ROOT / "data"
FIG_DIR = STUDY_ROOT / "figures"
ARTIFACT_DIR = STUDY_ROOT / "artifacts"
REPORT_DIR = STUDY_ROOT / "report"
PREVIOUS_STUDY_ROOT = STUDY_ROOT.parent / "cluster_state_holographic_sim"

for _path in (DATA_DIR, FIG_DIR, ARTIFACT_DIR, REPORT_DIR):
    _path.mkdir(parents=True, exist_ok=True)

TWO_PI = 2.0 * np.pi
OMEGA_Q = TWO_PI * 6.150e9
OMEGA_C = TWO_PI * 5.241e9
ALPHA = TWO_PI * (-255.0e6)
CHI = TWO_PI * (-2.84e6)
CHI_ABS = abs(CHI)
CHI_PRIME = TWO_PI * (-21.0e3)
KERR = TWO_PI * (-28.0e3)

N_CAV_OPT = 2
N_CAV_DEFAULT = 8
N_TR_DEFAULT = 2

ROTATION_S = 40.0e-9
DISPLACEMENT_S = 80.0e-9
SELECTIVE_S = 1.10e-6
CP_S = 160.0e-9
FREE_EVOLVE_S = np.pi / CHI_ABS

GRAPE_AMP_BOUND = TWO_PI * 50.0e6
GRAPE_DT_S = 4.0e-9

QUBIT_T1_S = 30.0e-6
QUBIT_T2_S = 20.0e-6
QUBIT_TPHI_S = pure_dephasing_time_from_t1_t2(t1_s=QUBIT_T1_S, t2_s=QUBIT_T2_S)
CAVITY_T1_S = 250.0e-6

LOGICAL_LABELS = ("|g,0>", "|g,1>", "|e,0>", "|e,1>")
PAULI_NAMES = ("I", "X", "Y", "Z")

TARGET_UNITARY = np.asarray(make_target("cluster", n_match=1), dtype=np.complex128)
TARGET = TargetUnitary(TARGET_UNITARY, ignore_global_phase=True)
IDEAL_DRIFT = DriftPhaseModel(chi=0.0, chi2=0.0, kerr=0.0)
PHYSICAL_FE_DRIFT = DriftPhaseModel(chi=CHI_ABS, chi2=0.0, kerr=0.0)


def json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        if np.iscomplexobj(value):
            return {
                "real": value.real.tolist(),
                "imag": value.imag.tolist(),
                "shape": list(value.shape),
            }
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


def build_model(*, n_cav: int = N_CAV_DEFAULT, n_tr: int = N_TR_DEFAULT) -> DispersiveTransmonCavityModel:
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


def logical_indices(n_cav: int) -> tuple[int, int, int, int]:
    return (0, 1, int(n_cav), int(n_cav) + 1)


def logical_subspace(n_cav: int) -> Subspace:
    return Subspace.custom(2 * int(n_cav), logical_indices(int(n_cav)), LOGICAL_LABELS)


def logical_basis_states(model: DispersiveTransmonCavityModel) -> list[qt.Qobj]:
    return [
        model.basis_state(0, 0),
        model.basis_state(0, 1),
        model.basis_state(1, 0),
        model.basis_state(1, 1),
    ]


def embed_target_unitary(n_cav: int) -> np.ndarray:
    full_dim = 2 * int(n_cav)
    embedded = np.eye(full_dim, dtype=np.complex128)
    idx = list(logical_indices(int(n_cav)))
    embedded[np.ix_(idx, idx)] = TARGET_UNITARY
    return embedded


def channel_transfer_summary() -> dict[str, Any]:
    channel = HolographicChannel.from_unitary(TARGET_UNITARY, physical_dim=2, bond_dim=2)
    kraus_ops = [np.asarray(op, dtype=np.complex128) for op in channel.kraus_ops]

    vec_super = np.zeros((4, 4), dtype=np.complex128)
    for op in kraus_ops:
        vec_super += np.kron(op.conj(), op)

    paulis = (
        np.eye(2, dtype=np.complex128),
        np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128),
        np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=np.complex128),
        np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128),
    )
    pauli_transfer = np.zeros((4, 4), dtype=np.complex128)
    action_on_pauli: dict[str, Any] = {}
    for col, (name, operator) in enumerate(zip(PAULI_NAMES, paulis, strict=True)):
        mapped = sum(k @ operator @ k.conj().T for k in kraus_ops)
        coeffs = np.array([0.5 * np.trace(base.conj().T @ mapped) for base in paulis], dtype=np.complex128)
        pauli_transfer[:, col] = coeffs
        action_on_pauli[name] = {
            "matrix": mapped,
            "coefficients": {PAULI_NAMES[row]: coeffs[row] for row in range(4)},
        }

    eigvals = np.linalg.eigvals(pauli_transfer)
    order = np.argsort(-np.abs(eigvals))
    eigvals = eigvals[order]
    if abs(eigvals[0]) <= 1.0e-15:
        xi = float("nan")
        ratio = float("nan")
    else:
        ratio = float(abs(eigvals[1] / eigvals[0]))
        if ratio == 0.0:
            xi = 0.0
        elif abs(ratio - 1.0) <= 1.0e-12:
            xi = float("inf")
        else:
            xi = float(-1.0 / np.log(ratio))

    traceless_block = pauli_transfer[1:, 1:]
    nilpotent_depth = None
    for depth in range(1, 5):
        if np.linalg.norm(np.linalg.matrix_power(traceless_block, depth), ord="fro") <= 1.0e-12:
            nilpotent_depth = depth
            break

    return {
        "kraus_ops": kraus_ops,
        "vec_superoperator": vec_super,
        "pauli_transfer": pauli_transfer,
        "eigenvalues": eigvals,
        "correlation_length_formula": "-1 / ln(|lambda_1 / lambda_0|)",
        "correlation_length": xi,
        "subleading_ratio": ratio,
        "nilpotent_depth_traceless": nilpotent_depth,
        "right_canonical_error": float(channel.right_canonical_error()),
        "kraus_completeness_error": float(channel.kraus_completeness_error()),
        "action_on_pauli": action_on_pauli,
        "fixed_point_density_matrix": 0.5 * np.eye(2, dtype=np.complex128),
    }


def selective_gate_active_tones(gate: Any, *, tol: float = 1.0e-6) -> int:
    if isinstance(gate, SNAP):
        return int(sum(abs(float(value)) > tol for value in gate.phases))
    if isinstance(gate, ConditionalPhaseSQR):
        return int(sum(abs(float(value)) > tol for value in gate.phases_n))
    if isinstance(gate, SQR):
        return int(
            sum(
                max(abs(float(theta)), abs(float(phi))) > tol
                for theta, phi in zip(gate.theta_n, gate.phi_n, strict=True)
            )
        )
    return 0


def sequence_gate_summary(sequence: GateSequence, *, tol: float = 1.0e-6) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    total_duration_ns = 0.0
    total_active_tones = 0
    max_active_tones = 0
    selective_gates = 0
    gate_counts: dict[str, int] = {}
    for gate in sequence.gates:
        gate_type = type(gate).__name__
        gate_counts[gate_type] = gate_counts.get(gate_type, 0) + 1
        duration_ns = float(getattr(gate, "duration", 0.0) * 1.0e9)
        active_tones = selective_gate_active_tones(gate, tol=tol)
        if active_tones > 0:
            selective_gates += 1
        max_active_tones = max(max_active_tones, active_tones)
        total_active_tones += active_tones
        total_duration_ns += duration_ns
        rows.append(
            {
                "name": str(gate.name),
                "type": gate_type,
                "duration_ns": duration_ns,
                "active_tones": active_tones,
            }
        )
    return {
        "rows": rows,
        "gate_count": len(sequence.gates),
        "gate_counts": gate_counts,
        "total_duration_ns": total_duration_ns,
        "selective_gate_count": selective_gates,
        "total_active_tones": total_active_tones,
        "max_active_tones": max_active_tones,
    }


def build_dr_snap_sequence(*, n_cav: int, blocks: int = 2) -> GateSequence:
    gates: list[Any] = []
    for block in range(blocks):
        gates.append(
            Displacement(
                name=f"D{block + 1}",
                alpha=0.20 + 0.0j,
                duration=DISPLACEMENT_S,
                optimize_time=False,
            )
        )
        gates.append(
            QubitRotation(
                name=f"R{block + 1}",
                theta=np.pi / 2.0,
                phi=0.0,
                duration=ROTATION_S,
                optimize_time=False,
            )
        )
        gates.append(
            SNAP(
                name=f"SNAP{block + 1}",
                phases=[0.0] * int(n_cav),
                duration=SELECTIVE_S,
                optimize_time=False,
            )
        )
    gates.append(
        Displacement(
            name=f"D{blocks + 1}",
            alpha=0.20 + 0.0j,
            duration=DISPLACEMENT_S,
            optimize_time=False,
        )
    )
    return GateSequence(gates=gates, n_cav=int(n_cav))


def build_dsqr_cp_sequence(*, n_cav: int, blocks: int = 2) -> GateSequence:
    theta_seed = [np.pi / 2.0] * int(n_cav)
    phi_seed = [0.0] * int(n_cav)
    gates: list[Any] = []
    for block in range(blocks):
        gates.append(
            Displacement(
                name=f"D{block + 1}",
                alpha=0.30 + 0.0j,
                duration=DISPLACEMENT_S,
                optimize_time=False,
            )
        )
        gates.append(
            SQR(
                name=f"S{2 * block + 1}",
                theta_n=theta_seed,
                phi_n=phi_seed,
                drift_model=IDEAL_DRIFT,
                duration=SELECTIVE_S,
                optimize_time=False,
            )
        )
        gates.append(
            ConditionalPhaseSQR(
                name=f"CP{block + 1}",
                phases_n=[0.0] * int(n_cav),
                drift_model=IDEAL_DRIFT,
                duration=CP_S,
                optimize_time=False,
            )
        )
        gates.append(
            SQR(
                name=f"S{2 * block + 2}",
                theta_n=theta_seed,
                phi_n=phi_seed,
                drift_model=IDEAL_DRIFT,
                duration=SELECTIVE_S,
                optimize_time=False,
            )
        )
    gates.append(
        Displacement(
            name=f"D{blocks + 1}",
            alpha=0.30 + 0.0j,
            duration=DISPLACEMENT_S,
            optimize_time=False,
        )
    )
    return GateSequence(gates=gates, n_cav=int(n_cav))


def build_dr_sqr_cp_sequence(*, n_cav: int, blocks: int = 2) -> GateSequence:
    gates: list[Any] = []
    for block in range(blocks):
        gates.append(
            Displacement(
                name=f"D{block + 1}",
                alpha=0.25 + 0.0j,
                duration=DISPLACEMENT_S,
                optimize_time=False,
            )
        )
        gates.append(
            QubitRotation(
                name=f"R{2 * block + 1}",
                theta=np.pi / 2.0,
                phi=0.0,
                duration=ROTATION_S,
                optimize_time=False,
            )
        )
        gates.append(
            SQR(
                name=f"S{block + 1}",
                theta_n=[np.pi / 2.0] * int(n_cav),
                phi_n=[0.0] * int(n_cav),
                drift_model=IDEAL_DRIFT,
                duration=SELECTIVE_S,
                optimize_time=False,
            )
        )
        gates.append(
            ConditionalPhaseSQR(
                name=f"CP{block + 1}",
                phases_n=[0.0] * int(n_cav),
                drift_model=IDEAL_DRIFT,
                duration=CP_S,
                optimize_time=False,
            )
        )
        gates.append(
            QubitRotation(
                name=f"R{2 * block + 2}",
                theta=np.pi / 2.0,
                phi=np.pi / 2.0,
                duration=ROTATION_S,
                optimize_time=False,
            )
        )
    gates.append(
        Displacement(
            name=f"D{blocks + 1}",
            alpha=0.15 + 0.0j,
            duration=DISPLACEMENT_S,
            optimize_time=False,
        )
    )
    return GateSequence(gates=gates, n_cav=int(n_cav))


def build_dr_fe_sequence(*, n_cav: int, blocks: int = 2) -> GateSequence:
    gates: list[Any] = []
    for block in range(blocks):
        gates.append(
            Displacement(
                name=f"D{block + 1}",
                alpha=0.20 + 0.0j,
                duration=DISPLACEMENT_S,
                optimize_time=False,
            )
        )
        gates.append(
            QubitRotation(
                name=f"R{block + 1}",
                theta=np.pi / 2.0,
                phi=0.0,
                duration=ROTATION_S,
                optimize_time=False,
            )
        )
        gates.append(
            FreeEvolveCondPhase(
                name=f"FE{block + 1}",
                duration=FREE_EVOLVE_S,
                drift_model=PHYSICAL_FE_DRIFT,
                optimize_time=True,
            )
        )
    gates.append(
        Displacement(
            name=f"D{blocks + 1}",
            alpha=0.20 + 0.0j,
            duration=DISPLACEMENT_S,
            optimize_time=False,
        )
    )
    gates.append(
        QubitRotation(
            name=f"R{blocks + 1}",
            theta=np.pi / 2.0,
            phi=np.pi / 2.0,
            duration=ROTATION_S,
            optimize_time=False,
        )
    )
    return GateSequence(gates=gates, n_cav=int(n_cav))


def extend_sequence_to_n_cav(sequence: GateSequence, target_n_cav: int) -> GateSequence:
    gates: list[Any] = []
    target_n_cav = int(target_n_cav)
    for gate in sequence.gates:
        if isinstance(gate, Displacement):
            gates.append(
                Displacement(
                    name=gate.name,
                    alpha=gate.alpha,
                    duration=float(gate.duration),
                    optimize_time=bool(getattr(gate, "optimize_time", False)),
                )
            )
        elif isinstance(gate, QubitRotation):
            gates.append(
                QubitRotation(
                    name=gate.name,
                    theta=float(gate.theta),
                    phi=float(gate.phi),
                    duration=float(gate.duration),
                    optimize_time=bool(getattr(gate, "optimize_time", False)),
                )
            )
        elif isinstance(gate, SNAP):
            phases = list(gate.phases) + [0.0] * max(0, target_n_cav - len(gate.phases))
            gates.append(
                SNAP(
                    name=gate.name,
                    phases=phases[:target_n_cav],
                    duration=float(gate.duration),
                    optimize_time=bool(getattr(gate, "optimize_time", False)),
                )
            )
        elif isinstance(gate, SQR):
            theta_n = list(gate.theta_n) + [0.0] * max(0, target_n_cav - len(gate.theta_n))
            phi_n = list(gate.phi_n) + [0.0] * max(0, target_n_cav - len(gate.phi_n))
            gates.append(
                SQR(
                    name=gate.name,
                    theta_n=theta_n[:target_n_cav],
                    phi_n=phi_n[:target_n_cav],
                    drift_model=gate.drift_model,
                    duration=float(gate.duration),
                    optimize_time=bool(getattr(gate, "optimize_time", False)),
                )
            )
        elif isinstance(gate, ConditionalPhaseSQR):
            phases_n = list(gate.phases_n) + [0.0] * max(0, target_n_cav - len(gate.phases_n))
            gates.append(
                ConditionalPhaseSQR(
                    name=gate.name,
                    phases_n=phases_n[:target_n_cav],
                    drift_model=gate.drift_model,
                    duration=float(gate.duration),
                    optimize_time=bool(getattr(gate, "optimize_time", False)),
                )
            )
        elif isinstance(gate, FreeEvolveCondPhase):
            gates.append(
                FreeEvolveCondPhase(
                    name=gate.name,
                    duration=float(gate.duration),
                    drift_model=gate.drift_model,
                    optimize_time=bool(getattr(gate, "optimize_time", False)),
                )
            )
        else:
            gates.append(copy.deepcopy(gate))
    return GateSequence(gates=gates, n_cav=target_n_cav)


def fit_sequence(
    sequence: GateSequence,
    *,
    max_amplitude: float | None = None,
    multistart: int = 4,
    maxiter: int = 250,
    duration_weight: float = 0.0,
) -> dict[str, Any]:
    subspace = logical_subspace(sequence.n_cav)
    constraints = None if max_amplitude is None else SynthesisConstraints(max_amplitude=float(max_amplitude))
    objectives = MultiObjective(fidelity_weight=1.0, leakage_weight=0.05, duration_weight=float(duration_weight))
    synthesizer = UnitarySynthesizer(
        primitives=copy.deepcopy(sequence.gates),
        subspace=subspace,
        objectives=objectives,
        leakage_penalty=LeakagePenalty(weight=0.05),
        synthesis_constraints=constraints,
        execution=ExecutionOptions(engine="auto", use_fast_path=True),
    )
    result = synthesizer.fit(target=TARGET, init_guess="heuristic", multistart=int(multistart), maxiter=int(maxiter))
    sub_op = np.asarray(result.simulation.subspace_operator, dtype=np.complex128)
    ideal_fidelity = float(subspace_unitary_fidelity(sub_op, TARGET_UNITARY, gauge="global"))
    gate_info = sequence_gate_summary(result.sequence)
    return {
        "result": result,
        "ideal_fidelity": ideal_fidelity,
        "objective": float(getattr(result, "objective", getattr(result, "objective_value", np.nan))),
        "success": bool(getattr(result, "success", True)),
        "message": str(getattr(result, "message", "")),
        "gate_summary": gate_info,
    }


def evaluate_sequence_ideal(sequence: GateSequence, *, n_cav: int) -> dict[str, Any]:
    seq = extend_sequence_to_n_cav(sequence, int(n_cav))
    subspace = logical_subspace(int(n_cav))
    simulation = synth_simulate_sequence(seq, subspace)
    full_operator = np.asarray(simulation.full_operator, dtype=np.complex128)
    sub_operator = np.asarray(simulation.subspace_operator, dtype=np.complex128)
    leakage = leakage_metrics(full_operator, subspace)
    return {
        "sequence": seq,
        "full_operator": full_operator,
        "subspace_operator": sub_operator,
        "fidelity": float(subspace_unitary_fidelity(sub_operator, TARGET_UNITARY, gauge="global")),
        "block_fidelity": float(
            subspace_unitary_fidelity(
                sub_operator,
                TARGET_UNITARY,
                gauge="block",
                block_slices=((0, 2), (1, 3)),
            )
        ),
        "leakage_average": float(leakage.average),
        "leakage_worst": float(leakage.worst),
    }


def waveform_bridge_supported(sequence: GateSequence) -> bool:
    supported = {QubitRotation, Displacement, SQR, ConditionalPhaseSQR}
    return all(type(gate) in supported for gate in sequence.gates)


def evaluate_sequence_pulse(
    sequence: GateSequence,
    *,
    n_cav: int,
    n_tr: int = N_TR_DEFAULT,
    dt_s: float = GRAPE_DT_S,
) -> dict[str, Any]:
    seq = extend_sequence_to_n_cav(sequence, int(n_cav))
    if not waveform_bridge_supported(seq):
        return {"supported": False, "reason": "waveform_bridge_gap"}
    model = build_model(n_cav=int(n_cav), n_tr=int(n_tr))
    frame = build_frame(model)
    system = CQEDSystemAdapter(model=model)
    waveform_seq = waveform_sequence_from_gates(seq, frame=frame)
    subspace = logical_subspace(int(n_cav))
    simulation = synth_simulate_sequence(
        waveform_seq,
        subspace,
        backend="pulse",
        target_subspace=TARGET_UNITARY,
        system=system,
        dt=dt_s,
        frame=frame,
    )
    full_operator = np.asarray(simulation.full_operator, dtype=np.complex128)
    leakage = leakage_metrics(full_operator, subspace)
    return {
        "supported": True,
        "sequence": seq,
        "full_operator": full_operator,
        "subspace_operator": np.asarray(simulation.subspace_operator, dtype=np.complex128),
        "fidelity": float(simulation.metrics.get("fidelity", np.nan)),
        "block_fidelity": float(simulation.metrics.get("block_fidelity", np.nan)),
        "leakage_average": float(leakage.average),
        "leakage_worst": float(leakage.worst),
        "metrics": dict(simulation.metrics),
    }


def reduced_cavity_density(state_vector: np.ndarray, *, n_cav: int, n_tr: int) -> qt.Qobj:
    psi = np.asarray(state_vector, dtype=np.complex128).reshape(int(n_tr), int(n_cav))
    rho_c = np.einsum("qn,qm->nm", psi, psi.conj())
    return qt.Qobj(rho_c, dims=[[int(n_cav)], [int(n_cav)]])


def candidate_wigner_summary(
    candidate_full_operator: np.ndarray,
    *,
    target_full_operator: np.ndarray,
    n_cav: int,
    n_tr: int = N_TR_DEFAULT,
    n_points: int = 81,
    extent: float = 4.5,
) -> dict[str, Any]:
    outputs: dict[str, Any] = {}
    for label, index in zip(LOGICAL_LABELS, logical_indices(int(n_cav)), strict=True):
        basis = np.zeros(2 * int(n_cav), dtype=np.complex128)
        basis[index] = 1.0
        target_state = target_full_operator @ basis
        candidate_state = candidate_full_operator @ basis
        target_rho = reduced_cavity_density(target_state, n_cav=int(n_cav), n_tr=int(n_tr))
        candidate_rho = reduced_cavity_density(candidate_state, n_cav=int(n_cav), n_tr=int(n_tr))
        xvec, yvec, target_w = cavity_wigner(target_rho, n_points=int(n_points), extent=float(extent))
        _, _, candidate_w = cavity_wigner(candidate_rho, xvec=xvec, yvec=yvec)
        fidelity = float(qt.fidelity(candidate_rho, target_rho) ** 2)
        outputs[label] = {
            "cavity_fidelity": fidelity,
            "wigner_l2": float(np.sqrt(np.mean((candidate_w - target_w) ** 2))),
            "xvec": xvec,
            "yvec": yvec,
            "target_wigner": target_w,
            "candidate_wigner": candidate_w,
            "difference_wigner": candidate_w - target_w,
            "target_density": np.asarray(target_rho.full(), dtype=np.complex128),
            "candidate_density": np.asarray(candidate_rho.full(), dtype=np.complex128),
        }
    return outputs


def build_grape_problem(
    *,
    model: DispersiveTransmonCavityModel,
    duration_ns: float,
    amp_bound_rad_s: float = GRAPE_AMP_BOUND,
    steps: int | None = None,
) -> Any:
    frame = build_frame(model)
    steps = int(steps or max(10, round((float(duration_ns) * 1.0e-9) / GRAPE_DT_S)))
    duration_s = float(duration_ns) * 1.0e-9
    return build_control_problem_from_model(
        model,
        frame=frame,
        time_grid=PiecewiseConstantTimeGrid.uniform(steps=steps, dt_s=duration_s / steps),
        channel_specs=(
            ModelControlChannelSpec(
                name="storage",
                target="storage",
                quadratures=("I", "Q"),
                amplitude_bounds=(-float(amp_bound_rad_s), float(amp_bound_rad_s)),
                export_channel="storage",
            ),
            ModelControlChannelSpec(
                name="qubit",
                target="qubit",
                quadratures=("I", "Q"),
                amplitude_bounds=(-float(amp_bound_rad_s), float(amp_bound_rad_s)),
                export_channel="qubit",
            ),
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


def run_grape_seed(problem: Any, *, seed: int, maxiter: int) -> Any:
    solver = GrapeSolver(
        GrapeConfig(
            maxiter=int(maxiter),
            seed=int(seed),
            random_scale=0.30,
            history_every=1,
        )
    )
    return solver.solve(problem)


def replay_compiled_sequence(
    *,
    model: DispersiveTransmonCavityModel,
    compiled: Any,
    drive_ops: dict[str, Any],
    noise: NoiseSpec | None = None,
    store_states: bool = False,
) -> list[dict[str, Any]]:
    frame = build_frame(model)
    basis_states = logical_basis_states(model)
    rows: list[dict[str, Any]] = []
    for label, initial_state in zip(LOGICAL_LABELS, basis_states, strict=True):
        result = pulse_simulate_sequence(
            model,
            compiled,
            initial_state,
            drive_ops,
            config=SimulationConfig(frame=frame, store_states=bool(store_states)),
            noise=noise,
        )
        rows.append({"label": label, "simulation": result})
    return rows


def replay_grape_operator(
    *,
    result: Any,
    problem: Any,
    model: DispersiveTransmonCavityModel,
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
        noise=noise,
        store_states=store_states,
    )
    indices = np.asarray(logical_indices(int(model.n_cav)), dtype=int)
    if noise is None:
        sub_operator = np.zeros((4, 4), dtype=np.complex128)
        leakage_values: list[float] = []
        photon_peaks: list[float] = []
        for column, row in enumerate(rows):
            state = row["simulation"].final_state
            vector = np.asarray(state.full(), dtype=np.complex128).reshape(-1)
            logical = vector[indices]
            sub_operator[:, column] = logical
            leakage_values.append(float(max(0.0, 1.0 - np.vdot(logical, logical).real)))
            expectations = row["simulation"].expectations
            photon_trace = expectations.get("n_c")
            if photon_trace is None:
                photon_trace = expectations.get("n_s")
            if photon_trace is not None and len(photon_trace) > 0:
                photon_peaks.append(float(np.max(np.real(photon_trace))))
        return {
            "pulses": pulses,
            "drive_ops": drive_ops,
            "pulse_meta": pulse_meta,
            "compiled": compiled,
            "subspace_operator": sub_operator,
            "fidelity": float(subspace_unitary_fidelity(sub_operator, TARGET_UNITARY, gauge="global")),
            "block_fidelity": float(
                subspace_unitary_fidelity(
                    sub_operator,
                    TARGET_UNITARY,
                    gauge="block",
                    block_slices=((0, 2), (1, 3)),
                )
            ),
            "leakage_average": float(np.mean(leakage_values)),
            "leakage_worst": float(np.max(leakage_values)),
            "max_transient_photon_number": float(max(photon_peaks) if photon_peaks else 0.0),
            "rows": rows,
        }

    return {
        "pulses": pulses,
        "drive_ops": drive_ops,
        "pulse_meta": pulse_meta,
        "compiled": compiled,
        "rows": rows,
    }


def default_noise_spec() -> NoiseSpec:
    return NoiseSpec(t1=QUBIT_T1_S, tphi=QUBIT_TPHI_S, kappa=1.0 / CAVITY_T1_S)


def superposition_state(n_cav: int, i: int, j: int, phase: complex = 1.0j) -> np.ndarray:
    vec = np.zeros(2 * int(n_cav), dtype=np.complex128)
    vec[i] = 1.0 / np.sqrt(2.0)
    vec[j] = phase / np.sqrt(2.0)
    return vec


def density_in_logical_block(state: qt.Qobj, *, n_cav: int) -> np.ndarray:
    rho = np.asarray(state.full(), dtype=np.complex128)
    idx = list(logical_indices(int(n_cav)))
    return rho[np.ix_(idx, idx)]
