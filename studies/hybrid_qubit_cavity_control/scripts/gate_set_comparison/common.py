"""Shared utilities for the 2x2 hybrid universal-control study.

This file keeps the study-specific glue on top of `cqed_sim` small and explicit:

- one device model and logical subspace definition,
- one set of target operators,
- fixed-duration primitive templates for each candidate gate library,
- helpers for synthesis, GRAPE replay, robustness, and plotting data.
"""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


STUDY_ROOT = Path(__file__).resolve().parents[2]
COMPONENT_NAME = Path(__file__).resolve().parent.name
DATA_DIR = STUDY_ROOT / "data" / COMPONENT_NAME
FIGURES_DIR = STUDY_ROOT / "figures" / COMPONENT_NAME
NOTES_DIR = STUDY_ROOT / "notes" / COMPONENT_NAME
REPORT_DIR = STUDY_ROOT / "report"

SIM_ROOT = Path(
    "C:/Users/dazzl/Box/Shyam Shankar Quantum Circuits Group/Users/Users_JianJun/cQED_simulation"
)
if str(SIM_ROOT) not in sys.path:
    sys.path.insert(0, str(SIM_ROOT))

from cqed_sim import (  # noqa: E402
    DispersiveTransmonCavityModel,
    FrameSpec,
    ModelControlChannelSpec,
    PiecewiseConstantTimeGrid,
    SequenceCompiler,
    SimulationConfig,
    simulate_sequence,
)
from cqed_sim.optimal_control import (  # noqa: E402
    ControlEvaluationCase,
    GrapeConfig,
    GrapeSolver,
    LeakagePenalty as OCLeakagePenalty,
    UnitaryObjective as OCUnitaryObjective,
    build_control_problem_from_model,
)
from cqed_sim.sim.noise import NoiseSpec, pure_dephasing_time_from_t1_t2  # noqa: E402
from cqed_sim.unitary_synthesis import (  # noqa: E402
    BlueSidebandExchange,
    ConditionalDisplacement,
    ConditionalPhaseSQR,
    Displacement,
    DriftPhaseModel,
    ExecutionOptions,
    FreeEvolveCondPhase,
    GateSequence,
    JaynesCummingsExchange,
    LeakagePenalty,
    MultiObjective,
    QubitRotation,
    SNAP,
    SQR,
    Subspace,
    TargetUnitary,
    UnitarySynthesizer,
    leakage_metrics,
    subspace_unitary_fidelity,
)


# Device parameters used throughout the study.
OMEGA_Q = 2.0 * np.pi * 6.150e9
OMEGA_C = 2.0 * np.pi * 5.241e9
ALPHA = 2.0 * np.pi * (-255.0e6)
CHI = 2.0 * np.pi * (-2.84e6)
CHI_ABS = abs(CHI)
CHI_PRIME = 2.0 * np.pi * (-21.0e3)
KERR = 2.0 * np.pi * (-28.0e3)

N_CAV_DEFAULT = 8
N_TR_DEFAULT = 2

# Duration priors used for library comparison.
QUBIT_ROTATION_S = 40.0e-9
DISPLACEMENT_S = 80.0e-9
SELECTIVE_S = 1.10e-6
CONDITIONAL_DISPLACEMENT_S = 160.0e-9
EXCHANGE_S = 120.0e-9
BLUE_SIDEBAND_S = 120.0e-9
NATIVE_WAIT_S = np.pi / CHI_ABS

# Reporting weights for the composite study-local score.
SCORE_LAMBDA_LEAK = 0.20
SCORE_LAMBDA_T = 0.02
SCORE_T_REF_NS = 1000.0

# Nominal decoherence values used for proxy comparisons and noisy GRAPE replay.
QUBIT_T1_S = 30.0e-6
QUBIT_T2_S = 20.0e-6
CAVITY_T1_S = 250.0e-6
QUBIT_TPHI_S = pure_dephasing_time_from_t1_t2(t1_s=QUBIT_T1_S, t2_s=QUBIT_T2_S)


@dataclass(frozen=True)
class SequenceCase:
    """A named synthesis ansatz."""

    label: str
    sequence: list[Any]
    description: str
    target_key: str
    multistart: int = 3
    maxiter: int = 180


@dataclass(frozen=True)
class GrapeCase:
    """A named GRAPE configuration."""

    label: str
    target_key: str
    duration_s: float
    steps: int
    maxiter: int = 120
    amp_bound_rad_s: float = 2.0 * np.pi * 40.0e6
    seed: int = 17


def ensure_dirs() -> None:
    for path in (DATA_DIR, FIGURES_DIR, NOTES_DIR, REPORT_DIR):
        path.mkdir(parents=True, exist_ok=True)


def build_model(*, n_cav: int = N_CAV_DEFAULT, n_tr: int = N_TR_DEFAULT) -> DispersiveTransmonCavityModel:
    return DispersiveTransmonCavityModel(
        omega_c=OMEGA_C,
        omega_q=OMEGA_Q,
        alpha=ALPHA,
        chi=CHI,
        chi_higher=(CHI_PRIME,),
        kerr=KERR,
        n_cav=int(n_cav),
        n_tr=int(n_tr),
    )


def build_frame(model: DispersiveTransmonCavityModel) -> FrameSpec:
    return FrameSpec(omega_c_frame=model.omega_c, omega_q_frame=model.omega_q)


def logical_subspace(*, n_cav: int = N_CAV_DEFAULT) -> Subspace:
    return Subspace.custom(
        full_dim=2 * int(n_cav),
        indices=(0, 1, int(n_cav), int(n_cav) + 1),
        labels=("|g,0>", "|g,1>", "|e,0>", "|e,1>"),
    )


def logical_indices(*, n_cav: int = N_CAV_DEFAULT) -> tuple[int, int, int, int]:
    return (0, 1, int(n_cav), int(n_cav) + 1)


def per_fock_block_slices() -> tuple[tuple[int, int], tuple[int, int]]:
    # Basis order is |g,0>, |g,1>, |e,0>, |e,1>.
    return ((0, 2), (1, 3))


def hadamard_cavity_target() -> np.ndarray:
    had = (1.0 / np.sqrt(2.0)) * np.array([[1.0, 1.0], [1.0, -1.0]], dtype=np.complex128)
    return np.kron(np.eye(2, dtype=np.complex128), had)


def cnot_cavity_to_qubit_target() -> np.ndarray:
    return np.array(
        [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
        ],
        dtype=np.complex128,
    )


def cz_target() -> np.ndarray:
    return np.diag([1.0, 1.0, 1.0, -1.0]).astype(np.complex128)


def target_matrix(target_key: str) -> np.ndarray:
    if target_key == "local_h":
        return hadamard_cavity_target()
    if target_key == "cx_c_to_q":
        return cnot_cavity_to_qubit_target()
    if target_key == "cz":
        return cz_target()
    raise ValueError(f"Unsupported target_key '{target_key}'.")


def target_name(target_key: str) -> str:
    return {
        "local_h": "I_q \\otimes H_c",
        "cx_c_to_q": "CX(c->q)",
        "cz": "CZ(q,c)",
    }[target_key]


def synthesis_score(
    *,
    fidelity: float,
    leakage_average: float,
    duration_ns: float,
    t_ref_ns: float = SCORE_T_REF_NS,
) -> float:
    return float(
        fidelity - SCORE_LAMBDA_LEAK * leakage_average - SCORE_LAMBDA_T * (duration_ns / t_ref_ns)
    )


def decoherence_proxy(duration_s: float) -> float:
    qubit_factor = np.exp(-float(duration_s) / QUBIT_T2_S)
    cavity_factor = np.exp(-float(duration_s) / CAVITY_T1_S)
    return float(min(qubit_factor, cavity_factor))


def native_wait_gate(name: str = "W") -> FreeEvolveCondPhase:
    return FreeEvolveCondPhase(
        name=name,
        duration=NATIVE_WAIT_S,
        drift_model=DriftPhaseModel(chi=CHI_ABS, chi2=0.0, kerr=0.0),
        optimize_time=False,
    )


def swap_like_exchange_gate(
    name: str,
    *,
    coupling_rad_s: float = 2.0 * np.pi * 2.0e6,
    phase: float = 0.0,
    duration_s: float = EXCHANGE_S,
) -> JaynesCummingsExchange:
    return JaynesCummingsExchange(
        name=name,
        coupling=float(coupling_rad_s),
        phase=float(phase),
        duration=float(duration_s),
        optimize_time=False,
    )


def blue_sideband_gate(
    name: str,
    *,
    coupling_rad_s: float = 2.0 * np.pi * 2.0e6,
    phase: float = 0.0,
    duration_s: float = BLUE_SIDEBAND_S,
) -> BlueSidebandExchange:
    return BlueSidebandExchange(
        name=name,
        coupling=float(coupling_rad_s),
        phase=float(phase),
        duration=float(duration_s),
        optimize_time=False,
    )


def library_a_local_sequence(*, n_cav: int = N_CAV_DEFAULT) -> list[Any]:
    return [
        Displacement(name="D1", alpha=0.20 + 0.0j, duration=DISPLACEMENT_S, optimize_time=False),
        SNAP(name="S1", phases=[0.0] * int(n_cav), duration=SELECTIVE_S, optimize_time=False),
        Displacement(name="D2", alpha=0.20 + 0.0j, duration=DISPLACEMENT_S, optimize_time=False),
    ]


def library_b_local_sequence(*, n_cav: int = N_CAV_DEFAULT) -> list[Any]:
    zeros = [0.0] * int(n_cav)
    return [
        QubitRotation(name="R1", theta=np.pi / 2.0, phi=0.0, duration=QUBIT_ROTATION_S, optimize_time=False),
        SQR(name="S1", theta_n=[0.10] * int(n_cav), phi_n=zeros, duration=SELECTIVE_S, optimize_time=False),
        QubitRotation(name="R2", theta=np.pi / 2.0, phi=np.pi / 2.0, duration=QUBIT_ROTATION_S, optimize_time=False),
        Displacement(name="D1", alpha=0.20 + 0.0j, duration=DISPLACEMENT_S, optimize_time=False),
        QubitRotation(name="R3", theta=np.pi / 2.0, phi=0.0, duration=QUBIT_ROTATION_S, optimize_time=False),
        SQR(name="S2", theta_n=[0.10] * int(n_cav), phi_n=zeros, duration=SELECTIVE_S, optimize_time=False),
        QubitRotation(name="R4", theta=np.pi / 2.0, phi=np.pi / 2.0, duration=QUBIT_ROTATION_S, optimize_time=False),
    ]


def library_c_local_sequence(*, n_cav: int = N_CAV_DEFAULT) -> list[Any]:
    return [
        QubitRotation(name="R1", theta=np.pi / 2.0, phi=0.0, duration=QUBIT_ROTATION_S, optimize_time=False),
        ConditionalDisplacement(name="CD1", alpha=0.10 + 0.0j, duration=CONDITIONAL_DISPLACEMENT_S, optimize_time=False),
        QubitRotation(name="R2", theta=np.pi / 2.0, phi=np.pi / 2.0, duration=QUBIT_ROTATION_S, optimize_time=False),
        ConditionalDisplacement(name="CD2", alpha=0.10 + 0.0j, duration=CONDITIONAL_DISPLACEMENT_S, optimize_time=False),
        QubitRotation(name="R3", theta=np.pi / 2.0, phi=0.0, duration=QUBIT_ROTATION_S, optimize_time=False),
    ]


def library_d_local_sequence(*, n_cav: int = N_CAV_DEFAULT) -> list[Any]:
    return [
        Displacement(name="D1", alpha=0.20 + 0.0j, duration=DISPLACEMENT_S, optimize_time=False),
        native_wait_gate("W1"),
        Displacement(name="D2", alpha=-0.20 + 0.0j, duration=DISPLACEMENT_S, optimize_time=False),
    ]


def library_a_entangler_sequence() -> list[Any]:
    return [
        QubitRotation(name="R1", theta=np.pi / 2.0, phi=0.0, duration=QUBIT_ROTATION_S, optimize_time=False),
        native_wait_gate("W"),
        QubitRotation(name="R2", theta=np.pi / 2.0, phi=0.0, duration=QUBIT_ROTATION_S, optimize_time=False),
    ]


def library_b_entangler_sequence(*, n_cav: int = N_CAV_DEFAULT) -> list[Any]:
    return [
        SQR(
            name="S",
            theta_n=[0.0, np.pi] + [0.0] * (int(n_cav) - 2),
            phi_n=[0.0] * int(n_cav),
            duration=SELECTIVE_S,
            optimize_time=False,
        )
    ]


def library_c_entangler_sequence(*, n_cav: int = N_CAV_DEFAULT) -> list[Any]:
    return [
        QubitRotation(name="R1", theta=np.pi / 2.0, phi=0.0, duration=QUBIT_ROTATION_S, optimize_time=False),
        ConditionalDisplacement(name="CD1", alpha=0.10 + 0.0j, duration=CONDITIONAL_DISPLACEMENT_S, optimize_time=False),
        QubitRotation(name="R2", theta=np.pi / 2.0, phi=np.pi / 2.0, duration=QUBIT_ROTATION_S, optimize_time=False),
        ConditionalDisplacement(name="CD2", alpha=0.10 + 0.0j, duration=CONDITIONAL_DISPLACEMENT_S, optimize_time=False),
        QubitRotation(name="R3", theta=np.pi / 2.0, phi=0.0, duration=QUBIT_ROTATION_S, optimize_time=False),
    ]


def library_d_entangler_sequence() -> list[Any]:
    return library_a_entangler_sequence()


def library_f_local_sequence() -> list[Any]:
    return [
        blue_sideband_gate("BS1"),
        QubitRotation(name="R1", theta=np.pi / 2.0, phi=0.0, duration=QUBIT_ROTATION_S, optimize_time=False),
        swap_like_exchange_gate("JC1"),
        QubitRotation(name="R2", theta=np.pi / 2.0, phi=np.pi / 2.0, duration=QUBIT_ROTATION_S, optimize_time=False),
        blue_sideband_gate("BS2"),
    ]


def library_f_entangler_sequence() -> list[Any]:
    return [
        QubitRotation(name="R1", theta=np.pi / 2.0, phi=0.0, duration=QUBIT_ROTATION_S, optimize_time=False),
        swap_like_exchange_gate("JC1"),
        blue_sideband_gate("BS1"),
        QubitRotation(name="R2", theta=np.pi / 2.0, phi=np.pi / 2.0, duration=QUBIT_ROTATION_S, optimize_time=False),
    ]


def sequence_cases(*, n_cav: int = N_CAV_DEFAULT) -> dict[str, SequenceCase]:
    return {
        "A_local": SequenceCase(
            label="A_local",
            sequence=library_a_local_sequence(n_cav=n_cav),
            description="Baseline dispersive local cavity control via D-SNAP-D.",
            target_key="local_h",
            multistart=3,
            maxiter=180,
        ),
        "B_local": SequenceCase(
            label="B_local",
            sequence=library_b_local_sequence(n_cav=n_cav),
            description="Selective hybrid control using SQR-mediated cavity action.",
            target_key="local_h",
            multistart=3,
            maxiter=180,
        ),
        "C_local": SequenceCase(
            label="C_local",
            sequence=library_c_local_sequence(n_cav=n_cav),
            description="ECD-like native conditional-displacement route for local cavity action.",
            target_key="local_h",
            multistart=3,
            maxiter=180,
        ),
        "D_local": SequenceCase(
            label="D_local",
            sequence=library_d_local_sequence(n_cav=n_cav),
            description="Minimal native library local-control probe.",
            target_key="local_h",
            multistart=3,
            maxiter=180,
        ),
        "F_local": SequenceCase(
            label="F_local",
            sequence=library_f_local_sequence(),
            description="Native SWAP-/sideband-style route using blue-sideband and Jaynes-Cummings exchange.",
            target_key="local_h",
            multistart=3,
            maxiter=220,
        ),
        "A_ent": SequenceCase(
            label="A_ent",
            sequence=library_a_entangler_sequence(),
            description="Dispersive wait plus fast qubit rotations for CX(c->q).",
            target_key="cx_c_to_q",
            multistart=2,
            maxiter=120,
        ),
        "B_ent": SequenceCase(
            label="B_ent",
            sequence=library_b_entangler_sequence(n_cav=n_cav),
            description="Single selective-hybrid entangler (SQR).",
            target_key="cx_c_to_q",
            multistart=1,
            maxiter=80,
        ),
        "C_ent": SequenceCase(
            label="C_ent",
            sequence=library_c_entangler_sequence(n_cav=n_cav),
            description="ECD-like conditional-displacement entangler ansatz.",
            target_key="cx_c_to_q",
            multistart=3,
            maxiter=180,
        ),
        "D_ent": SequenceCase(
            label="D_ent",
            sequence=library_d_entangler_sequence(),
            description="Minimal native entangler library.",
            target_key="cx_c_to_q",
            multistart=2,
            maxiter=120,
        ),
        "F_ent": SequenceCase(
            label="F_ent",
            sequence=library_f_entangler_sequence(),
            description="Native SWAP-/sideband-style entangler ansatz using exchange primitives.",
            target_key="cx_c_to_q",
            multistart=3,
            maxiter=220,
        ),
    }


def grape_cases() -> dict[str, GrapeCase]:
    return {
        "E_local_320": GrapeCase(label="E_local_320", target_key="local_h", duration_s=320.0e-9, steps=16),
        "E_local_480": GrapeCase(label="E_local_480", target_key="local_h", duration_s=480.0e-9, steps=24),
        "E_local_640": GrapeCase(label="E_local_640", target_key="local_h", duration_s=640.0e-9, steps=32),
        "E_ent_240": GrapeCase(label="E_ent_240", target_key="cx_c_to_q", duration_s=240.0e-9, steps=12),
        "E_ent_400": GrapeCase(label="E_ent_400", target_key="cx_c_to_q", duration_s=400.0e-9, steps=20),
        "E_ent_560": GrapeCase(label="E_ent_560", target_key="cx_c_to_q", duration_s=560.0e-9, steps=28),
    }


def _synthesis_target(target_key: str) -> TargetUnitary:
    return TargetUnitary(target_matrix(target_key), ignore_global_phase=True)


def _base_synth_kwargs() -> dict[str, Any]:
    return {
        "optimize_times": False,
        "optimizer": "powell",
        "objectives": MultiObjective(fidelity_weight=1.0, leakage_weight=0.05),
        "leakage_penalty": LeakagePenalty(weight=0.05),
        "execution": ExecutionOptions(engine="auto", use_fast_path=True),
    }


def evaluate_sequence_result(
    *,
    label: str,
    category: str,
    description: str,
    target_key: str,
    result: Any,
) -> dict[str, Any]:
    sim = result.simulation
    target = target_matrix(target_key)
    strict = float(subspace_unitary_fidelity(sim.subspace_operator, target, gauge="global"))
    block = float(
        subspace_unitary_fidelity(
            sim.subspace_operator,
            target,
            gauge="block",
            block_slices=per_fock_block_slices(),
        )
    )
    metrics = sim.metrics
    duration_ns = float(result.sequence.total_duration() * 1.0e9)
    leakage_avg = float(metrics.get("leakage_average", np.nan))
    leakage_worst = float(metrics.get("leakage_worst", np.nan))
    return {
        "label": label,
        "category": category,
        "description": description,
        "target_key": target_key,
        "target_name": target_name(target_key),
        "strict_fidelity": strict,
        "block_fidelity": block,
        "leakage_average": leakage_avg,
        "leakage_worst": leakage_worst,
        "duration_ns": duration_ns,
        "gate_count": len(result.sequence.gates),
        "score_strict": synthesis_score(
            fidelity=strict,
            leakage_average=leakage_avg,
            duration_ns=duration_ns,
        ),
        "score_block": synthesis_score(
            fidelity=block,
            leakage_average=leakage_avg,
            duration_ns=duration_ns,
        ),
        "decoherence_proxy": decoherence_proxy(duration_ns * 1.0e-9),
        "sequence": result.sequence.serialize(),
        "time_parameters": result.sequence.serialize_time_parameters(),
        "simulation_metrics": {key: float(value) if isinstance(value, (np.floating, float, int)) else value for key, value in metrics.items()},
        "objective": float(result.objective),
        "success": bool(result.success),
        "report": result.report,
    }


def run_sequence_case(case: SequenceCase, *, n_cav: int = N_CAV_DEFAULT) -> dict[str, Any]:
    model = build_model(n_cav=n_cav)
    subspace = logical_subspace(n_cav=n_cav)
    synthesizer = UnitarySynthesizer(
        subspace=subspace,
        primitives=case.sequence,
        target=_synthesis_target(case.target_key),
        model=model,
        seed=11,
        **_base_synth_kwargs(),
    )
    result = synthesizer.fit(multistart=case.multistart, maxiter=case.maxiter)
    return evaluate_sequence_result(
        label=case.label,
        category="decomposition",
        description=case.description,
        target_key=case.target_key,
        result=result,
    )


def _logical_basis_states(model: DispersiveTransmonCavityModel) -> list[Any]:
    return [
        model.basis_state(0, 0),
        model.basis_state(0, 1),
        model.basis_state(1, 0),
        model.basis_state(1, 1),
    ]


def replay_grape_operator(
    *,
    result: Any,
    problem: Any,
    model: DispersiveTransmonCavityModel,
    frame: FrameSpec,
    n_cav: int = N_CAV_DEFAULT,
) -> tuple[np.ndarray, list[float]]:
    pulses, drive_ops, _meta = result.to_pulses()
    compiler = SequenceCompiler(dt=1.0e-9)
    compiled = compiler.compile(pulses, t_end=problem.time_grid.duration_s)
    full_basis = _logical_basis_states(model)
    indices = np.asarray(logical_indices(n_cav=n_cav), dtype=int)
    sub_op = np.zeros((4, 4), dtype=np.complex128)
    leakage: list[float] = []
    for column, initial_state in enumerate(full_basis):
        sim = simulate_sequence(
            model,
            compiled,
            initial_state,
            drive_ops,
            config=SimulationConfig(frame=frame),
        )
        final = np.asarray(sim.final_state.full(), dtype=np.complex128).reshape(-1)
        logical = final[indices]
        sub_op[:, column] = logical
        leakage.append(float(max(0.0, 1.0 - np.vdot(logical, logical).real)))
    return sub_op, leakage


def replay_grape_noise_metrics(
    *,
    result: Any,
    problem: Any,
    model: DispersiveTransmonCavityModel,
    frame: FrameSpec,
) -> dict[str, Any]:
    noisy = result.evaluate_with_simulator(
        problem,
        cases=(
            ControlEvaluationCase(
                model=model,
                frame=frame,
                noise=NoiseSpec(t1=QUBIT_T1_S, tphi=QUBIT_TPHI_S, kappa=1.0 / CAVITY_T1_S),
                label="nominal_open_system",
            ),
        ),
    )
    return {
        "aggregate_fidelity": float(noisy.metrics["aggregate_fidelity"]),
        "aggregate_leakage": float(noisy.metrics.get("aggregate_leakage", np.nan)),
    }


def run_grape_case(case: GrapeCase, *, n_cav: int = N_CAV_DEFAULT) -> dict[str, Any]:
    model = build_model(n_cav=n_cav)
    frame = build_frame(model)
    subspace = logical_subspace(n_cav=n_cav)
    problem = build_control_problem_from_model(
        model,
        frame=frame,
        time_grid=PiecewiseConstantTimeGrid.uniform(steps=case.steps, dt_s=case.duration_s / case.steps),
        channel_specs=(
            ModelControlChannelSpec(
                name="storage",
                target="storage",
                quadratures=("I", "Q"),
                amplitude_bounds=(-case.amp_bound_rad_s, case.amp_bound_rad_s),
                export_channel="storage",
            ),
            ModelControlChannelSpec(
                name="qubit",
                target="qubit",
                quadratures=("I", "Q"),
                amplitude_bounds=(-case.amp_bound_rad_s, case.amp_bound_rad_s),
                export_channel="qubit",
            ),
        ),
        objectives=(
            OCUnitaryObjective(
                target_operator=target_matrix(case.target_key),
                subspace=subspace,
                ignore_global_phase=True,
                name=case.label,
            ),
        ),
        penalties=(
            OCLeakagePenalty(weight=0.02, subspace=subspace),
        ),
    )
    solver = GrapeSolver(GrapeConfig(maxiter=case.maxiter, seed=case.seed, random_scale=0.40))
    result = solver.solve(problem)
    sub_op, leakage = replay_grape_operator(result=result, problem=problem, model=model, frame=frame, n_cav=n_cav)
    strict = float(subspace_unitary_fidelity(sub_op, target_matrix(case.target_key), gauge="global"))
    block = float(
        subspace_unitary_fidelity(
            sub_op,
            target_matrix(case.target_key),
            gauge="block",
            block_slices=per_fock_block_slices(),
        )
    )
    noise_metrics = replay_grape_noise_metrics(result=result, problem=problem, model=model, frame=frame)
    duration_ns = float(case.duration_s * 1.0e9)
    leakage_avg = float(np.mean(leakage))
    leakage_worst = float(np.max(leakage))
    return {
        "label": case.label,
        "category": "grape",
        "description": f"Waveform-level GRAPE reference for {target_name(case.target_key)}.",
        "target_key": case.target_key,
        "target_name": target_name(case.target_key),
        "strict_fidelity": strict,
        "block_fidelity": block,
        "leakage_average": leakage_avg,
        "leakage_worst": leakage_worst,
        "duration_ns": duration_ns,
        "gate_count": 1,
        "score_strict": synthesis_score(
            fidelity=strict,
            leakage_average=leakage_avg,
            duration_ns=duration_ns,
        ),
        "score_block": synthesis_score(
            fidelity=block,
            leakage_average=leakage_avg,
            duration_ns=duration_ns,
        ),
        "decoherence_proxy": decoherence_proxy(case.duration_s),
        "objective": float(result.objective_value),
        "success": bool(result.success),
        "message": str(result.message),
        "nominal_fidelity": float(result.metrics["nominal_fidelity"]),
        "control_metrics": {key: float(value) if isinstance(value, (np.floating, float, int)) else value for key, value in result.metrics.items()},
        "noise_metrics": noise_metrics,
        "steps": int(case.steps),
        "maxiter": int(case.maxiter),
        "amp_bound_rad_s": float(case.amp_bound_rad_s),
        "max_abs_command_rad_s": float(result.schedule.max_abs_amplitude()),
    }


def perturb_sequence_record(
    record: dict[str, Any],
    *,
    amplitude_scale: float = 1.0,
    duration_scale: float = 1.0,
    chi_scale: float = 1.0,
    n_cav: int = N_CAV_DEFAULT,
) -> GateSequence:
    gates: list[Any] = []
    for gate in record["sequence"]:
        gate_type = gate["type"]
        duration = float(gate["duration"]) * float(duration_scale)
        params = gate["parameters"]
        if gate_type == "Displacement":
            gates.append(
                Displacement(
                    name=gate["name"],
                    alpha=complex(params[0] * amplitude_scale, params[1] * amplitude_scale),
                    duration=duration,
                    optimize_time=False,
                )
            )
        elif gate_type == "SNAP":
            gates.append(
                SNAP(
                    name=gate["name"],
                    phases=[float(value) for value in params],
                    duration=duration,
                    optimize_time=False,
                )
            )
        elif gate_type == "QubitRotation":
            gates.append(
                QubitRotation(
                    name=gate["name"],
                    theta=float(params[0]) * amplitude_scale,
                    phi=float(params[1]),
                    duration=duration,
                    optimize_time=False,
                )
            )
        elif gate_type == "SQR":
            half = len(params) // 2
            gates.append(
                SQR(
                    name=gate["name"],
                    theta_n=[float(value) * amplitude_scale for value in params[:half]],
                    phi_n=[float(value) for value in params[half:]],
                    duration=duration,
                    optimize_time=False,
                )
            )
        elif gate_type == "FreeEvolveCondPhase":
            gates.append(
                FreeEvolveCondPhase(
                    name=gate["name"],
                    duration=duration,
                    drift_model=DriftPhaseModel(chi=CHI_ABS * chi_scale, chi2=0.0, kerr=0.0),
                    optimize_time=False,
                )
            )
        elif gate_type == "ConditionalDisplacement":
            gates.append(
                ConditionalDisplacement(
                    name=gate["name"],
                    alpha=complex(float(params[0]) * amplitude_scale, float(params[1]) * amplitude_scale),
                    duration=duration,
                    optimize_time=False,
                )
            )
        elif gate_type == "JaynesCummingsExchange":
            gates.append(
                JaynesCummingsExchange(
                    name=gate["name"],
                    coupling=float(params[0]) * amplitude_scale,
                    phase=float(params[1]),
                    duration=duration,
                    optimize_time=False,
                )
            )
        elif gate_type == "BlueSidebandExchange":
            gates.append(
                BlueSidebandExchange(
                    name=gate["name"],
                    coupling=float(params[0]) * amplitude_scale,
                    phase=float(params[1]),
                    duration=duration,
                    optimize_time=False,
                )
            )
        else:
            raise ValueError(f"Unsupported gate_type '{gate_type}' in perturbation helper.")
    return GateSequence(gates=gates, n_cav=int(n_cav))


def replay_sequence_record(
    record: dict[str, Any],
    *,
    amplitude_scale: float = 1.0,
    duration_scale: float = 1.0,
    chi_scale: float = 1.0,
    n_cav: int = N_CAV_DEFAULT,
) -> dict[str, float]:
    sequence = perturb_sequence_record(
        record,
        amplitude_scale=amplitude_scale,
        duration_scale=duration_scale,
        chi_scale=chi_scale,
        n_cav=n_cav,
    )
    subspace = logical_subspace(n_cav=n_cav)
    target = target_matrix(record["target_key"])
    full_operator = np.asarray(sequence.unitary(backend="ideal"), dtype=np.complex128)
    sub_operator = subspace.restrict_operator(full_operator)
    strict = float(subspace_unitary_fidelity(sub_operator, target, gauge="global"))
    block = float(
        subspace_unitary_fidelity(
            sub_operator,
            target,
            gauge="block",
            block_slices=per_fock_block_slices(),
        )
    )
    leak = leakage_metrics(full_operator, subspace)
    return {
        "strict_fidelity": strict,
        "block_fidelity": block,
        "leakage_average": float(leak.average),
        "leakage_worst": float(leak.worst),
        "duration_ns": float(sequence.total_duration() * 1.0e9),
    }


def json_dump(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
