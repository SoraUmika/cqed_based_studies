"""Common utilities for the multitone SQR arbitrary block-rotation study."""

from __future__ import annotations

import json
import os
import time
from functools import partial
from pathlib import Path
from typing import Any, Sequence

import matplotlib
import numpy as np

matplotlib.use("Agg")

from scipy.linalg import norm, polar, sqrtm
from scipy.spatial.transform import Rotation

import runtime_compat  # noqa: F401
import qutip as qt

from cqed_sim.calibration.conditioned_multitone import (
    ConditionedMultitoneRunConfig,
    ConditionedMultitoneWaveform,
    ConditionedQubitTargets,
    build_conditioned_multitone_waveform,
    compile_conditioned_multitone_waveform,
)
from cqed_sim.core import DispersiveTransmonCavityModel, FrameSpec
from cqed_sim.core.conventions import qubit_cavity_block_indices
from cqed_sim.core.frequencies import carrier_for_transition_frequency, manifold_transition_frequency
from cqed_sim.pulses.envelopes import gaussian_area_fraction, gaussian_envelope
from cqed_sim.pulses.pulse import Pulse
from cqed_sim.sequence import SequenceCompiler
from cqed_sim.sim import SimulationConfig, prepare_simulation


STUDY_DIR = Path(__file__).resolve().parent.parent
DATA_DIR = STUDY_DIR / "data"
FIGURES_DIR = STUDY_DIR / "figures"
ARTIFACTS_DIR = STUDY_DIR / "artifacts"
REPORT_DIR = STUDY_DIR / "report"

for _path in (DATA_DIR, FIGURES_DIR, ARTIFACTS_DIR, REPORT_DIR):
    _path.mkdir(parents=True, exist_ok=True)


STYLE_PATH = (
    STUDY_DIR.parent.parent / ".github" / "skills" / "publication-figures" / "assets" / "cqed_style.mplstyle"
)

TWO_PI = 2.0 * np.pi

OMEGA_Q = TWO_PI * 6.150e9
OMEGA_C = TWO_PI * 5.241e9
ALPHA = TWO_PI * (-255.0e6)
CHI = TWO_PI * (-2.84e6)
CHI_PRIME = TWO_PI * (-21.0e3)
KERR = 0.0

N_TR = 2
N_CAV_PADDING = 2
DEFAULT_DT = 4.0e-9
DEFAULT_SIGMA_FRACTION = 1.0 / 6.0

SUCCESS_TIERS = {
    "strong": 0.999,
    "moderate": 0.99,
    "weak": 0.95,
}

PAULI_X = np.asarray([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
PAULI_Y = np.asarray([[0.0, -1.0j], [1.0j, 0.0]], dtype=np.complex128)
PAULI_Z = np.asarray([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)
PAULIS = (PAULI_X, PAULI_Y, PAULI_Z)
IDEAL_X_PI = np.asarray([[0.0, -1.0j], [-1.0j, 0.0]], dtype=np.complex128)


def apply_plot_style() -> None:
    import matplotlib.pyplot as plt

    if STYLE_PATH.exists():
        plt.style.use(str(STYLE_PATH))


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
                "real": np.real(value).tolist(),
                "imag": np.imag(value).tolist(),
                "shape": list(value.shape),
            }
        return value.tolist()
    if isinstance(value, complex):
        return {"real": float(value.real), "imag": float(value.imag)}
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def save_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(json_ready(payload), indent=2, sort_keys=True)
    temp_path = path.with_suffix(path.suffix + ".tmp")
    last_error: OSError | None = None
    for _ in range(5):
        try:
            temp_path.write_text(text, encoding="utf-8")
            os.replace(str(temp_path), str(path))
            return
        except OSError as exc:
            last_error = exc
            time.sleep(0.25)
    if last_error is not None:
        raise last_error


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def save_waveform_npz(path: Path, waveform_samples: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        path,
        time_s=np.asarray(waveform_samples["time_s"], dtype=float),
        baseband_real=np.asarray(waveform_samples["baseband_real"], dtype=float),
        baseband_imag=np.asarray(waveform_samples["baseband_imag"], dtype=float),
        distorted_real=np.asarray(waveform_samples["distorted_real"], dtype=float),
        distorted_imag=np.asarray(waveform_samples["distorted_imag"], dtype=float),
    )


def duration_from_chi_t(chi_t_2pi: float, *, chi: float = CHI) -> float:
    return float(chi_t_2pi) / (abs(float(chi)) / TWO_PI)


def build_model(*, include_chi_prime: bool, n_active: int, n_tr: int = N_TR) -> DispersiveTransmonCavityModel:
    n_cav = int(n_active) + int(N_CAV_PADDING)
    chi_higher = (CHI_PRIME,) if include_chi_prime else ()
    return DispersiveTransmonCavityModel(
        omega_q=OMEGA_Q,
        omega_c=OMEGA_C,
        alpha=ALPHA,
        chi=CHI,
        chi_higher=chi_higher,
        kerr=KERR,
        n_cav=n_cav,
        n_tr=int(n_tr),
    )


def build_frame(model: DispersiveTransmonCavityModel) -> FrameSpec:
    return FrameSpec(omega_q_frame=float(model.omega_q), omega_c_frame=float(model.omega_c))


def logical_levels(n_active: int) -> tuple[int, ...]:
    return tuple(range(int(n_active)))


def logical_indices(model: DispersiveTransmonCavityModel, levels: Sequence[int]) -> tuple[int, ...]:
    indices: list[int] = []
    for level in levels:
        indices.extend(int(index) for index in qubit_cavity_block_indices(int(model.n_cav), int(level)))
    return tuple(indices)


def manifold_transition_frequencies_hz(
    model: DispersiveTransmonCavityModel,
    levels: Sequence[int],
    frame: FrameSpec,
) -> np.ndarray:
    return np.asarray(
        [
            float(manifold_transition_frequency(model, int(level), frame=frame) / TWO_PI)
            for level in levels
        ],
        dtype=float,
    )


def min_transition_spacing_hz(model: DispersiveTransmonCavityModel, levels: Sequence[int], frame: FrameSpec) -> float:
    freqs = np.sort(manifold_transition_frequencies_hz(model, levels, frame))
    if freqs.size < 2:
        return float("inf")
    return float(np.min(np.diff(freqs)))


def make_run_config(
    model: DispersiveTransmonCavityModel,
    *,
    n_active: int,
    duration_s: float,
    dt_s: float = DEFAULT_DT,
    sigma_fraction: float = DEFAULT_SIGMA_FRACTION,
) -> ConditionedMultitoneRunConfig:
    frame = build_frame(model)
    levels = logical_levels(n_active)
    return ConditionedMultitoneRunConfig(
        frame=frame,
        duration_s=float(duration_s),
        dt_s=float(dt_s),
        sigma_fraction=float(sigma_fraction),
        tone_cutoff=1.0e-12,
        include_all_levels=False,
        max_step_s=float(dt_s),
        fock_fqs_hz=tuple(float(freq) for freq in manifold_transition_frequencies_hz(model, levels, frame)),
    )


def rz(angle: float) -> np.ndarray:
    half = 0.5 * float(angle)
    return np.asarray(
        [[np.exp(-1.0j * half), 0.0], [0.0, np.exp(1.0j * half)]],
        dtype=np.complex128,
    )


def ry(angle: float) -> np.ndarray:
    half = 0.5 * float(angle)
    return np.asarray(
        [[np.cos(half), -np.sin(half)], [np.sin(half), np.cos(half)]],
        dtype=np.complex128,
    )


def su2_from_zyz(phi: float, theta: float, lam: float) -> np.ndarray:
    return rz(phi) @ ry(theta) @ rz(lam)


def haar_random_su2(rng: np.random.Generator) -> np.ndarray:
    u1, u2, u3 = rng.random(3)
    a = np.sqrt(1.0 - u1) * np.exp(1.0j * TWO_PI * u2)
    b = np.sqrt(u1) * np.exp(1.0j * TWO_PI * u3)
    return np.asarray([[a, b], [-np.conj(b), np.conj(a)]], dtype=np.complex128)


def block_to_seed_angles(block: np.ndarray) -> tuple[float, float]:
    vector = np.asarray(block, dtype=np.complex128)[:, 0]
    if abs(vector[0]) > 1.0e-14:
        vector = vector * np.exp(-1.0j * np.angle(vector[0]))
    norm_val = float(np.linalg.norm(vector))
    if norm_val > 0.0:
        vector = vector / norm_val
    theta = 2.0 * np.arccos(np.clip(abs(vector[0]), 0.0, 1.0))
    phi = 0.0 if abs(vector[1]) < 1.0e-14 else float(np.angle(vector[1]))
    return float(theta), float(phi)


def conditioned_targets_from_blocks(blocks: Sequence[np.ndarray]) -> ConditionedQubitTargets:
    spec = {
        int(index): block_to_seed_angles(np.asarray(block, dtype=np.complex128))
        for index, block in enumerate(blocks)
    }
    return ConditionedQubitTargets.from_spec(spec, n_levels=len(spec))


def block_diag_target(blocks: Sequence[np.ndarray]) -> np.ndarray:
    n_blocks = len(blocks)
    operator = np.zeros((2 * n_blocks, 2 * n_blocks), dtype=np.complex128)
    for index, block in enumerate(blocks):
        operator[2 * index : 2 * index + 2, 2 * index : 2 * index + 2] = np.asarray(block, dtype=np.complex128)
    return operator


def make_family_blocks(family: str, n_active: int, *, rng: np.random.Generator | None = None) -> tuple[tuple[np.ndarray, ...], dict[str, Any]]:
    blocks: list[np.ndarray] = []
    metadata: dict[str, Any] = {"family": str(family), "n_active": int(n_active)}
    if family == "A":
        target_level = min(max(int(n_active) // 2, 0), int(n_active) - 1)
        metadata["selected_level"] = int(target_level)
        for level in range(int(n_active)):
            if level == target_level:
                blocks.append(su2_from_zyz(0.35 * np.pi, 0.70 * np.pi, -0.28 * np.pi))
            else:
                blocks.append(np.eye(2, dtype=np.complex128))
    elif family == "B":
        unitary = su2_from_zyz(-0.22 * np.pi, 0.62 * np.pi, 0.31 * np.pi)
        blocks = [unitary.copy() for _ in range(int(n_active))]
    elif family == "C":
        angles: list[tuple[float, float, float]] = []
        for level in range(int(n_active)):
            phi = -0.30 * np.pi + 0.16 * np.pi * float(level)
            theta = min(0.18 * np.pi + 0.14 * np.pi * float(level), 0.92 * np.pi)
            lam = 0.27 * np.pi - 0.09 * np.pi * float(level)
            angles.append((phi, theta, lam))
            blocks.append(su2_from_zyz(phi, theta, lam))
        metadata["smooth_euler_angles"] = angles
    elif family == "D":
        if rng is None:
            raise ValueError("Random family D requires an RNG.")
        blocks = [haar_random_su2(rng) for _ in range(int(n_active))]
    else:
        raise ValueError(f"Unsupported target family '{family}'.")

    euler_rows: list[dict[str, float]] = []
    for index, block in enumerate(blocks):
        seed_theta, seed_phi = block_to_seed_angles(block)
        euler_rows.append(
            {
                "level": int(index),
                "seed_theta": float(seed_theta),
                "seed_phi": float(seed_phi),
            }
        )
    metadata["seed_rows"] = euler_rows
    return tuple(np.asarray(block, dtype=np.complex128) for block in blocks), metadata


def nearest_su2(matrix: np.ndarray) -> np.ndarray:
    unitary, _ = polar(np.asarray(matrix, dtype=np.complex128))
    det = np.linalg.det(unitary)
    phase = np.angle(det) / 2.0
    fixed = unitary * np.exp(-1.0j * phase)
    return np.asarray(fixed, dtype=np.complex128)


def su2_matrix_sqrt(matrix: np.ndarray) -> np.ndarray:
    return nearest_su2(np.asarray(sqrtm(np.asarray(matrix, dtype=np.complex128)), dtype=np.complex128))


def echoed_half_target_blocks(blocks: Sequence[np.ndarray]) -> tuple[np.ndarray, ...]:
    half_blocks: list[np.ndarray] = []
    x_dagger = IDEAL_X_PI.conj().T
    for block in blocks:
        root = su2_matrix_sqrt(np.asarray(block, dtype=np.complex128))
        half_blocks.append(nearest_su2(x_dagger @ root))
    return tuple(half_blocks)


def su2_to_bloch_rotation(block: np.ndarray) -> Rotation:
    unitary = nearest_su2(block)
    rot = np.zeros((3, 3), dtype=float)
    dagger = unitary.conj().T
    for i, sigma_i in enumerate(PAULIS):
        for j, sigma_j in enumerate(PAULIS):
            rot[i, j] = 0.5 * float(np.real(np.trace(sigma_i @ unitary @ sigma_j @ dagger)))
    return Rotation.from_matrix(rot)


def process_fidelity(target_operator: np.ndarray, actual_operator: np.ndarray) -> float:
    target = np.asarray(target_operator, dtype=np.complex128)
    actual = np.asarray(actual_operator, dtype=np.complex128)
    dim = float(target.shape[0])
    overlap = np.trace(target.conj().T @ actual)
    return float(np.clip(abs(overlap) ** 2 / (dim * dim), 0.0, 1.0))


def average_gate_fidelity(target_operator: np.ndarray, actual_operator: np.ndarray) -> float:
    dim = float(np.asarray(target_operator).shape[0])
    proc = process_fidelity(target_operator, actual_operator)
    return float((dim * proc + 1.0) / (dim + 1.0))


def frobenius_error(target_operator: np.ndarray, actual_operator: np.ndarray) -> float:
    target = np.asarray(target_operator, dtype=np.complex128)
    actual = np.asarray(actual_operator, dtype=np.complex128)
    denom = float(norm(target, ord="fro"))
    raw = float(norm(actual - target, ord="fro"))
    if denom <= 0.0:
        return raw
    return raw / denom


def operator_norm_error(target_operator: np.ndarray, actual_operator: np.ndarray) -> float:
    target = np.asarray(target_operator, dtype=np.complex128)
    actual = np.asarray(actual_operator, dtype=np.complex128)
    return float(norm(actual - target, ord=2))


def unitarity_error(actual_operator: np.ndarray) -> float:
    actual = np.asarray(actual_operator, dtype=np.complex128)
    ident = np.eye(actual.shape[0], dtype=np.complex128)
    return float(norm(actual.conj().T @ actual - ident, ord=2))


def restricted_blocks(restricted_operator: np.ndarray) -> tuple[np.ndarray, ...]:
    op = np.asarray(restricted_operator, dtype=np.complex128)
    n_blocks = op.shape[0] // 2
    return tuple(op[2 * index : 2 * index + 2, 2 * index : 2 * index + 2] for index in range(n_blocks))


def block_rotation_metrics(target_block: np.ndarray, actual_block: np.ndarray) -> dict[str, float]:
    target = np.asarray(target_block, dtype=np.complex128)
    actual = np.asarray(actual_block, dtype=np.complex128)
    proc = process_fidelity(target, actual)
    avg = average_gate_fidelity(target, actual)
    fro = frobenius_error(target, actual)
    op_err = operator_norm_error(target, actual)
    unit_err = unitarity_error(actual)

    target_rot = su2_to_bloch_rotation(target)
    actual_rot = su2_to_bloch_rotation(actual)
    target_vec = target_rot.as_rotvec()
    actual_vec = actual_rot.as_rotvec()
    target_angle = float(np.linalg.norm(target_vec))
    actual_angle = float(np.linalg.norm(actual_vec))
    angle_error = abs(actual_angle - target_angle)
    if target_angle < 1.0e-9 or actual_angle < 1.0e-9:
        axis_error = float("nan")
    else:
        target_axis = target_vec / target_angle
        actual_axis = actual_vec / actual_angle
        axis_error = float(np.arccos(np.clip(np.dot(target_axis, actual_axis), -1.0, 1.0)))
    error_unitary = nearest_su2(target.conj().T @ actual)
    error_rotvec = su2_to_bloch_rotation(error_unitary).as_rotvec()
    residual_z = float(abs(error_rotvec[2]))
    transverse = float(np.linalg.norm(error_rotvec[:2]))
    total = float(np.linalg.norm(error_rotvec))
    return {
        "process_fidelity": float(proc),
        "average_gate_fidelity": float(avg),
        "frobenius_error": float(fro),
        "operator_norm_error": float(op_err),
        "unitarity_error": float(unit_err),
        "target_rotation_angle_rad": float(target_angle),
        "actual_rotation_angle_rad": float(actual_angle),
        "rotation_angle_error_rad": float(angle_error),
        "rotation_axis_error_rad": float(axis_error),
        "residual_z_error_rad": float(residual_z),
        "transverse_error_rad": float(transverse),
        "error_rotation_angle_rad": float(total),
        "error_rotvec_x_rad": float(error_rotvec[0]),
        "error_rotvec_y_rad": float(error_rotvec[1]),
        "error_rotvec_z_rad": float(error_rotvec[2]),
    }


def active_subspace_metrics(target_operator: np.ndarray, actual_operator: np.ndarray) -> dict[str, float]:
    return {
        "process_fidelity": process_fidelity(target_operator, actual_operator),
        "average_gate_fidelity": average_gate_fidelity(target_operator, actual_operator),
        "frobenius_error": frobenius_error(target_operator, actual_operator),
        "operator_norm_error": operator_norm_error(target_operator, actual_operator),
        "restricted_unitarity_error": unitarity_error(actual_operator),
    }


def basis_state(model: DispersiveTransmonCavityModel, qubit_level: int, cavity_level: int) -> qt.Qobj:
    return model.basis_state(int(qubit_level), int(cavity_level))


def conditioned_qubit_density_matrix(state: qt.Qobj, n_cav: int, level: int) -> tuple[qt.Qobj, float]:
    rho = state if state.isoper else state.proj()
    matrix = np.asarray(rho.full(), dtype=np.complex128)
    idx_g = int(level)
    idx_e = int(n_cav) + int(level)
    block = matrix[np.ix_([idx_g, idx_e], [idx_g, idx_e])]
    population = float(np.real(np.trace(block)))
    if population <= 0.0:
        return qt.Qobj(np.zeros((2, 2), dtype=np.complex128), dims=[[2], [2]]), 0.0
    return qt.Qobj(block / population, dims=[[2], [2]]), population


def bloch_xyz_from_density_matrix(rho_q: qt.Qobj) -> tuple[float, float, float]:
    matrix = np.asarray(rho_q.full(), dtype=np.complex128)
    x = float(np.real(np.trace(matrix @ PAULI_X)))
    y = float(np.real(np.trace(matrix @ PAULI_Y)))
    z = float(np.real(np.trace(matrix @ PAULI_Z)))
    return x, y, z


def simulate_logical_basis_operator(
    model: DispersiveTransmonCavityModel,
    waveform: ConditionedMultitoneWaveform,
    run_config: ConditionedMultitoneRunConfig,
    levels: Sequence[int],
) -> np.ndarray:
    compiled = compile_conditioned_multitone_waveform(waveform, run_config)
    config = SimulationConfig(frame=run_config.frame, store_states=False)
    session = prepare_simulation(model, compiled, waveform.drive_ops, config=config, e_ops={})
    full_dim = int(model.n_tr) * int(model.n_cav)
    basis_states: list[qt.Qobj] = []
    for level in levels:
        basis_states.append(basis_state(model, 0, int(level)))
        basis_states.append(basis_state(model, 1, int(level)))
    logical_dim = len(basis_states)
    columns = np.zeros((full_dim, logical_dim), dtype=np.complex128)
    for column, psi0 in enumerate(basis_states):
        result = session.run(psi0)
        columns[:, column] = np.asarray(result.final_state.full(), dtype=np.complex128).reshape(-1)
    indices = np.asarray(logical_indices(model, levels), dtype=int)
    return np.asarray(columns[indices, :], dtype=np.complex128)


def shift_pulse(pulse: Pulse, *, t0: float, label: str | None = None) -> Pulse:
    return Pulse(
        channel=str(pulse.channel),
        t0=float(t0),
        duration=float(pulse.duration),
        envelope=pulse.envelope,
        carrier=float(pulse.carrier),
        phase=float(pulse.phase),
        amp=float(pulse.amp),
        drag=float(pulse.drag),
        sample_rate=None if pulse.sample_rate is None else float(pulse.sample_rate),
        label=label if label is not None else pulse.label,
    )


def make_gaussian_qubit_rotation_pulse(
    model: DispersiveTransmonCavityModel,
    frame: FrameSpec,
    *,
    theta: float,
    phase: float,
    duration_s: float,
    channel: str = "qubit",
    manifold_level: int = 0,
    sigma_fraction: float = 0.25,
    drag: float = 0.0,
    t0: float = 0.0,
    label: str | None = None,
) -> Pulse:
    transition = manifold_transition_frequency(model, int(manifold_level), frame=frame)
    carrier = carrier_for_transition_frequency(transition)
    area_fraction = gaussian_area_fraction(float(sigma_fraction))
    amplitude = float(theta) / (2.0 * float(duration_s) * float(area_fraction))
    envelope = partial(gaussian_envelope, sigma=float(sigma_fraction), center=0.5)
    return Pulse(
        channel=str(channel),
        t0=float(t0),
        duration=float(duration_s),
        envelope=envelope,
        carrier=float(carrier),
        phase=float(phase),
        amp=float(amplitude),
        drag=float(drag),
        label=label,
    )


def compile_pulse_sequence(pulses: Sequence[Pulse], *, dt_s: float, total_duration_s: float) -> Any:
    compiler = SequenceCompiler(dt=float(dt_s))
    return compiler.compile(list(pulses), t_end=float(total_duration_s + dt_s))


def _simulation_session_for_compiled(
    model: DispersiveTransmonCavityModel,
    compiled: Any,
    *,
    frame: FrameSpec,
    drive_ops: dict[str, str],
) -> Any:
    config = SimulationConfig(frame=frame, store_states=False)
    return prepare_simulation(model, compiled, drive_ops, config=config, e_ops={})


def simulate_full_operator_on_logical_inputs(
    model: DispersiveTransmonCavityModel,
    compiled: Any,
    *,
    frame: FrameSpec,
    drive_ops: dict[str, str],
    levels: Sequence[int],
) -> np.ndarray:
    session = _simulation_session_for_compiled(model, compiled, frame=frame, drive_ops=drive_ops)
    full_dim = int(model.n_tr) * int(model.n_cav)
    operator = np.zeros((full_dim, full_dim), dtype=np.complex128)
    for level in levels:
        indices = qubit_cavity_block_indices(int(model.n_cav), int(level))
        for qubit_level, input_index in enumerate(indices):
            psi0 = basis_state(model, qubit_level, int(level))
            result = session.run(psi0)
            operator[:, int(input_index)] = np.asarray(result.final_state.full(), dtype=np.complex128).reshape(-1)
    return operator


def restricted_operator_from_full(
    full_operator: np.ndarray,
    model: DispersiveTransmonCavityModel,
    levels: Sequence[int],
) -> np.ndarray:
    indices = np.asarray(logical_indices(model, levels), dtype=int)
    return np.asarray(full_operator[np.ix_(indices, indices)], dtype=np.complex128)


def simulate_state_for_compiled(
    model: DispersiveTransmonCavityModel,
    compiled: Any,
    *,
    frame: FrameSpec,
    drive_ops: dict[str, str],
    initial_state: qt.Qobj,
) -> qt.Qobj:
    session = _simulation_session_for_compiled(model, compiled, frame=frame, drive_ops=drive_ops)
    return session.run(initial_state).final_state


def state_validation_summary_for_compiled(
    model: DispersiveTransmonCavityModel,
    compiled: Any,
    *,
    frame: FrameSpec,
    drive_ops: dict[str, str],
    levels: Sequence[int],
    target_operator: np.ndarray,
) -> dict[str, Any]:
    coeffs = np.ones(len(levels), dtype=np.complex128) / np.sqrt(len(levels))
    g_state = logical_state_qobj(model, levels, coeffs, [1.0, 0.0])
    plus_state = logical_state_qobj(model, levels, coeffs, [1.0 / np.sqrt(2.0), 1.0 / np.sqrt(2.0)])
    rows = []
    for label, psi0 in (("cavity_superposition_ground", g_state), ("cavity_superposition_plus", plus_state)):
        actual = simulate_state_for_compiled(model, compiled, frame=frame, drive_ops=drive_ops, initial_state=psi0)
        ideal = apply_target_operator_to_state(model, levels, target_operator, psi0)
        fidelity = float(qt.fidelity(actual, ideal) ** 2)
        populations = []
        bloch_rows = []
        for level in levels:
            rho_q, pop = conditioned_qubit_density_matrix(actual, int(model.n_cav), int(level))
            x, y, z = bloch_xyz_from_density_matrix(rho_q)
            populations.append(float(pop))
            bloch_rows.append({"level": int(level), "population": float(pop), "x": float(x), "y": float(y), "z": float(z)})
        rows.append(
            {
                "label": label,
                "state_fidelity": float(fidelity),
                "cavity_level_populations": populations,
                "conditioned_bloch_vectors": bloch_rows,
            }
        )
    return {"states": rows}


def channel_waveform_samples(compiled: Any, channel: str = "qubit") -> dict[str, Any]:
    compiled_channel = compiled.channels[str(channel)]
    return {
        "time_s": np.asarray(compiled.tlist, dtype=float).tolist(),
        "baseband_real": np.real(np.asarray(compiled_channel.baseband, dtype=np.complex128)).tolist(),
        "baseband_imag": np.imag(np.asarray(compiled_channel.baseband, dtype=np.complex128)).tolist(),
        "distorted_real": np.real(np.asarray(compiled_channel.distorted, dtype=np.complex128)).tolist(),
        "distorted_imag": np.imag(np.asarray(compiled_channel.distorted, dtype=np.complex128)).tolist(),
    }


def crosstalk_matrix(
    model: DispersiveTransmonCavityModel,
    run_config: ConditionedMultitoneRunConfig,
    tone_specs: Sequence[Any],
    levels: Sequence[int],
) -> np.ndarray:
    level_tuple = tuple(int(level) for level in levels)
    matrix = np.zeros((len(level_tuple), len(level_tuple)), dtype=float)
    by_manifold = {int(spec.manifold): spec for spec in tone_specs}
    for column, target_level in enumerate(level_tuple):
        tone = by_manifold.get(int(target_level))
        if tone is None:
            continue
        waveform = build_conditioned_multitone_waveform([tone], run_config, label=f"tone_{target_level}")
        restricted = simulate_logical_basis_operator(model, waveform, run_config, level_tuple)
        blocks = restricted_blocks(restricted)
        for row, block in enumerate(blocks):
            matrix[row, column] = float(abs(block[1, 0]) ** 2)
    return matrix


def crosstalk_summary(matrix: np.ndarray) -> dict[str, float]:
    arr = np.asarray(matrix, dtype=float)
    if arr.size == 0:
        return {
            "diagonal_mean": float("nan"),
            "diagonal_min": float("nan"),
            "offdiag_mean": float("nan"),
            "offdiag_max": float("nan"),
        }
    diag = np.diag(arr)
    mask = ~np.eye(arr.shape[0], dtype=bool)
    offdiag = arr[mask]
    return {
        "diagonal_mean": float(np.mean(diag)) if diag.size else float("nan"),
        "diagonal_min": float(np.min(diag)) if diag.size else float("nan"),
        "offdiag_mean": float(np.mean(offdiag)) if offdiag.size else 0.0,
        "offdiag_max": float(np.max(offdiag)) if offdiag.size else 0.0,
    }


def logical_state_qobj(
    model: DispersiveTransmonCavityModel,
    levels: Sequence[int],
    cavity_coeffs: Sequence[complex],
    qubit_state: Sequence[complex],
) -> qt.Qobj:
    coeffs = np.asarray(cavity_coeffs, dtype=np.complex128)
    qvec = np.asarray(qubit_state, dtype=np.complex128)
    if coeffs.size != len(levels):
        raise ValueError("cavity_coeffs must match the number of logical levels.")
    full = np.zeros(int(model.n_tr) * int(model.n_cav), dtype=np.complex128)
    for amp, level in zip(coeffs, levels, strict=True):
        idx_g, idx_e = qubit_cavity_block_indices(int(model.n_cav), int(level))
        full[int(idx_g)] = amp * qvec[0]
        full[int(idx_e)] = amp * qvec[1]
    norm_val = float(np.linalg.norm(full))
    if norm_val > 0.0:
        full = full / norm_val
    return qt.Qobj(full.reshape((-1, 1)), dims=[[int(model.n_tr), int(model.n_cav)], [1, 1]])


def apply_target_operator_to_state(
    model: DispersiveTransmonCavityModel,
    levels: Sequence[int],
    target_operator: np.ndarray,
    state: qt.Qobj,
) -> qt.Qobj:
    restricted = np.asarray(target_operator, dtype=np.complex128)
    indices = np.asarray(logical_indices(model, levels), dtype=int)
    vector = np.asarray(state.full(), dtype=np.complex128).reshape(-1)
    logical = vector[indices]
    evolved = restricted @ logical
    full = np.zeros_like(vector)
    full[indices] = evolved
    return qt.Qobj(full.reshape((-1, 1)), dims=state.dims)


def simulate_state(
    model: DispersiveTransmonCavityModel,
    waveform: ConditionedMultitoneWaveform,
    run_config: ConditionedMultitoneRunConfig,
    initial_state: qt.Qobj,
) -> qt.Qobj:
    compiled = compile_conditioned_multitone_waveform(waveform, run_config)
    config = SimulationConfig(frame=run_config.frame, store_states=False)
    session = prepare_simulation(model, compiled, waveform.drive_ops, config=config, e_ops={})
    return session.run(initial_state).final_state


def state_validation_summary(
    model: DispersiveTransmonCavityModel,
    waveform: ConditionedMultitoneWaveform,
    run_config: ConditionedMultitoneRunConfig,
    levels: Sequence[int],
    target_operator: np.ndarray,
) -> dict[str, Any]:
    coeffs = np.ones(len(levels), dtype=np.complex128) / np.sqrt(len(levels))
    g_state = logical_state_qobj(model, levels, coeffs, [1.0, 0.0])
    plus_state = logical_state_qobj(model, levels, coeffs, [1.0 / np.sqrt(2.0), 1.0 / np.sqrt(2.0)])
    rows = []
    for label, psi0 in (("cavity_superposition_ground", g_state), ("cavity_superposition_plus", plus_state)):
        actual = simulate_state(model, waveform, run_config, psi0)
        ideal = apply_target_operator_to_state(model, levels, target_operator, psi0)
        fidelity = float(qt.fidelity(actual, ideal) ** 2)
        populations = []
        bloch_rows = []
        for level in levels:
            rho_q, pop = conditioned_qubit_density_matrix(actual, int(model.n_cav), int(level))
            x, y, z = bloch_xyz_from_density_matrix(rho_q)
            populations.append(float(pop))
            bloch_rows.append({"level": int(level), "population": float(pop), "x": float(x), "y": float(y), "z": float(z)})
        rows.append(
            {
                "label": label,
                "state_fidelity": float(fidelity),
                "cavity_level_populations": populations,
                "conditioned_bloch_vectors": bloch_rows,
            }
        )
    return {"states": rows}


def classify_success(avg_gate_fidelity_value: float) -> str:
    value = float(avg_gate_fidelity_value)
    if value >= SUCCESS_TIERS["strong"]:
        return "strong"
    if value >= SUCCESS_TIERS["moderate"]:
        return "moderate"
    if value >= SUCCESS_TIERS["weak"]:
        return "weak"
    return "failure"


def classify_failure_mode(summary_row: dict[str, Any]) -> str:
    if float(summary_row.get("leakage_outside_target_mean", 0.0)) > 0.05:
        return "leakage outside active subspace"
    if float(summary_row.get("crosstalk_offdiag_max", 0.0)) > 0.15:
        return "off-resonant crosstalk"
    if float(summary_row.get("min_transition_spacing_hz", float("inf"))) < 1.5e6:
        return "spectral crowding"
    if float(summary_row.get("restricted_unitarity_error", 0.0)) > 0.1:
        return "nonunitary restricted dynamics"
    return "insufficient waveform degrees of freedom"