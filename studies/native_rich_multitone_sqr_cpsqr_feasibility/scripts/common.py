"""Common helpers for the final native/rich multitone SQR/CPSQR study."""

from __future__ import annotations

import json
import os
import time
from dataclasses import dataclass
from functools import partial
from pathlib import Path
from typing import Any, Sequence

import matplotlib
import numpy as np

matplotlib.use("Agg")

from scipy.linalg import norm, polar
from scipy.spatial.transform import Rotation

import runtime_compat  # noqa: F401
import qutip as qt

from cqed_sim.calibration.conditioned_multitone import (
    ConditionedMultitoneCorrections,
    ConditionedMultitoneRunConfig,
    ConditionedQubitTargets,
    build_conditioned_multitone_tones,
    build_conditioned_multitone_waveform,
    compile_conditioned_multitone_waveform,
)
from cqed_sim.calibration.targeted_subspace_multitone import build_block_rotation_target_operator
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
PI_PULSE_DURATION_S = 40.0e-9
PI_PULSE_SIGMA_FRACTION = 0.25

PAULI_X = np.asarray([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
PAULI_Y = np.asarray([[0.0, -1.0j], [1.0j, 0.0]], dtype=np.complex128)
PAULI_Z = np.asarray([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)
PAULIS = (PAULI_X, PAULI_Y, PAULI_Z)
IDEAL_X_PI = np.asarray([[0.0, -1.0j], [-1.0j, 0.0]], dtype=np.complex128)


@dataclass(frozen=True)
class TargetSpec:
    family: str
    theta_values: tuple[float, ...]
    phi_values: tuple[float, ...]
    metadata: dict[str, Any]


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


def build_model(
    *,
    include_chi_prime: bool,
    n_active: int,
    n_tr: int = N_TR,
    n_cav_padding: int = N_CAV_PADDING,
) -> DispersiveTransmonCavityModel:
    n_cav = int(n_active) + int(n_cav_padding)
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
        [float(manifold_transition_frequency(model, int(level), frame=frame) / TWO_PI) for level in levels],
        dtype=float,
    )


def detuning_matrix_hz(
    model: DispersiveTransmonCavityModel,
    levels: Sequence[int],
    frame: FrameSpec,
) -> np.ndarray:
    freqs = manifold_transition_frequencies_hz(model, levels, frame)
    return freqs[:, None] - freqs[None, :]


def make_run_config(
    model: DispersiveTransmonCavityModel,
    *,
    n_active: int,
    duration_s: float,
    dt_s: float = DEFAULT_DT,
    sigma_fraction: float = DEFAULT_SIGMA_FRACTION,
) -> ConditionedMultitoneRunConfig:
    frame = build_frame(model)
    return ConditionedMultitoneRunConfig(
        frame=frame,
        duration_s=float(duration_s),
        dt_s=float(dt_s),
        sigma_fraction=float(sigma_fraction),
        tone_cutoff=1.0e-12,
        include_all_levels=False,
        max_step_s=float(dt_s),
        # Patched cqed_sim interprets fock_fqs_hz as an absolute transition-frequency
        # override and subtracts the frame internally. Older study-local helpers passed
        # frame-relative frequencies here; that pattern must not be reused.
        fock_fqs_hz=None,
    )


def target_spec(
    family: str,
    n_active: int,
    *,
    rng: np.random.Generator | None = None,
    seed: int | None = None,
) -> TargetSpec:
    n = int(n_active)
    if family == "smooth_x":
        theta = tuple(float(min(0.30 * np.pi + 0.18 * np.pi * level, 0.92 * np.pi)) for level in range(n))
        metadata = {
            "description": "Smooth monotone x-axis SQR profile",
            "formula": "theta_n = min(0.30 pi + 0.18 pi n, 0.92 pi)",
        }
    elif family == "staggered_x":
        theta_rows = [0.22 * np.pi + 0.58 * np.pi * (level % 2) + 0.06 * np.pi * (level // 2) for level in range(n)]
        theta = tuple(float(min(value, 0.94 * np.pi)) for value in theta_rows)
        metadata = {
            "description": "Harder staggered x-axis SQR profile",
            "formula": "theta_n = min(0.22 pi + 0.58 pi (n mod 2) + 0.06 pi floor(n/2), 0.94 pi)",
        }
    elif family == "random_x":
        generator = np.random.default_rng(int(seed)) if rng is None else rng
        theta = tuple(float(generator.uniform(0.18 * np.pi, 0.92 * np.pi)) for _ in range(n))
        metadata = {
            "description": "Random x-axis SQR profile",
            "formula": "theta_n ~ Uniform(0.18 pi, 0.92 pi)",
            "seed": None if seed is None else int(seed),
        }
    else:
        raise ValueError(f"Unsupported target family '{family}'.")
    phi = tuple(0.0 for _ in range(n))
    metadata["family"] = str(family)
    metadata["n_active"] = int(n)
    metadata["theta_values_rad"] = list(theta)
    metadata["phi_values_rad"] = list(phi)
    return TargetSpec(family=str(family), theta_values=theta, phi_values=phi, metadata=metadata)


def conditioned_targets_from_target_spec(spec: TargetSpec) -> ConditionedQubitTargets:
    pairs = list(zip(spec.theta_values, spec.phi_values, strict=True))
    return ConditionedQubitTargets.from_spec(pairs, n_levels=len(pairs))


def gaussian_samples(count: int, sigma_fraction: float = DEFAULT_SIGMA_FRACTION) -> np.ndarray:
    grid = np.linspace(0.0, 1.0, int(count), endpoint=False)
    sigma = float(sigma_fraction)
    samples = np.exp(-0.5 * ((grid - 0.5) / sigma) ** 2).astype(np.complex128)
    peak = float(np.max(np.abs(samples))) if samples.size else 0.0
    if peak > 0.0:
        samples = samples / peak
    return samples


def orthogonal_basis(count: int, order: int, *, kind: str) -> list[np.ndarray]:
    grid = np.linspace(0.0, 1.0, int(count), endpoint=False)
    rows: list[np.ndarray] = []
    for idx in range(1, int(order) + 1):
        if kind == "cos":
            row = np.cos(np.pi * idx * (grid - 0.5))
        elif kind == "sin":
            row = np.sin(np.pi * idx * grid)
        else:
            raise ValueError(f"Unsupported basis kind '{kind}'.")
        row = row - float(np.mean(row))
        norm_val = float(np.linalg.norm(row))
        if norm_val > 0.0:
            row = row / norm_val
        rows.append(np.asarray(row, dtype=float))
    return rows


def normalize_complex_samples(samples: np.ndarray, *, max_norm: float = 2.0) -> np.ndarray:
    arr = np.asarray(samples, dtype=np.complex128).reshape(-1)
    peak = float(np.max(np.abs(arr))) if arr.size else 1.0
    if peak > max_norm and peak > 0.0:
        arr = arr * (float(max_norm) / peak)
    return arr


def build_target_operator(spec: TargetSpec, levels: Sequence[int]) -> np.ndarray:
    targets = conditioned_targets_from_target_spec(spec)
    return build_block_rotation_target_operator(targets, logical_levels=levels)


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
    return raw if denom <= 0.0 else raw / denom


def unitarity_error(actual_operator: np.ndarray) -> float:
    actual = np.asarray(actual_operator, dtype=np.complex128)
    ident = np.eye(actual.shape[0], dtype=np.complex128)
    return float(norm(actual.conj().T @ actual - ident, ord=2))


def nearest_su2(matrix: np.ndarray) -> np.ndarray:
    unitary, _ = polar(np.asarray(matrix, dtype=np.complex128))
    det = np.linalg.det(unitary)
    phase = np.angle(det) / 2.0
    return np.asarray(unitary * np.exp(-1.0j * phase), dtype=np.complex128)


def su2_to_bloch_rotation(block: np.ndarray) -> Rotation:
    unitary = nearest_su2(block)
    rot = np.zeros((3, 3), dtype=float)
    dagger = unitary.conj().T
    for i, sigma_i in enumerate(PAULIS):
        for j, sigma_j in enumerate(PAULIS):
            rot[i, j] = 0.5 * float(np.real(np.trace(sigma_i @ unitary @ sigma_j @ dagger)))
    return Rotation.from_matrix(rot)


def restricted_blocks(restricted_operator: np.ndarray) -> tuple[np.ndarray, ...]:
    op = np.asarray(restricted_operator, dtype=np.complex128)
    return tuple(op[2 * idx : 2 * idx + 2, 2 * idx : 2 * idx + 2] for idx in range(op.shape[0] // 2))


def block_rotation_metrics(target_block: np.ndarray, actual_block: np.ndarray) -> dict[str, float]:
    target = np.asarray(target_block, dtype=np.complex128)
    actual = np.asarray(actual_block, dtype=np.complex128)
    proc = process_fidelity(target, actual)
    avg = average_gate_fidelity(target, actual)
    target_rot = su2_to_bloch_rotation(target)
    actual_rot = su2_to_bloch_rotation(actual)
    target_vec = target_rot.as_rotvec()
    actual_vec = actual_rot.as_rotvec()
    target_angle = float(np.linalg.norm(target_vec))
    actual_angle = float(np.linalg.norm(actual_vec))
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
    return {
        "process_fidelity": float(proc),
        "average_gate_fidelity": float(avg),
        "frobenius_error": float(frobenius_error(target, actual)),
        "unitarity_error": float(unitarity_error(actual)),
        "target_rotation_angle_rad": float(target_angle),
        "actual_rotation_angle_rad": float(actual_angle),
        "rotation_angle_error_rad": float(abs(actual_angle - target_angle)),
        "rotation_axis_error_rad": float(axis_error),
        "residual_z_error_rad": float(residual_z),
        "transverse_error_rad": float(transverse),
        "error_rotation_angle_rad": float(np.linalg.norm(error_rotvec)),
        "error_rotvec_x_rad": float(error_rotvec[0]),
        "error_rotvec_y_rad": float(error_rotvec[1]),
        "error_rotvec_z_rad": float(error_rotvec[2]),
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
    return (
        float(np.real(np.trace(matrix @ PAULI_X))),
        float(np.real(np.trace(matrix @ PAULI_Y))),
        float(np.real(np.trace(matrix @ PAULI_Z))),
    )


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
    sigma_fraction: float = PI_PULSE_SIGMA_FRACTION,
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


def logical_state_qobj(
    model: DispersiveTransmonCavityModel,
    levels: Sequence[int],
    cavity_coeffs: Sequence[complex],
    qubit_state: Sequence[complex],
) -> qt.Qobj:
    coeffs = np.asarray(cavity_coeffs, dtype=np.complex128)
    qvec = np.asarray(qubit_state, dtype=np.complex128)
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
    indices = np.asarray(logical_indices(model, levels), dtype=int)
    vector = np.asarray(state.full(), dtype=np.complex128).reshape(-1)
    logical = vector[indices]
    evolved = np.asarray(target_operator, dtype=np.complex128) @ logical
    full = np.zeros_like(vector)
    full[indices] = evolved
    return qt.Qobj(full.reshape((-1, 1)), dims=state.dims)


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


def corrections_to_dict(corrections: ConditionedMultitoneCorrections) -> dict[str, list[float]]:
    return {
        "d_lambda": [float(x) for x in corrections.d_lambda],
        "d_alpha": [float(x) for x in corrections.d_alpha],
        "d_omega_rad_s": [float(x) for x in corrections.d_omega_rad_s],
        "d_omega_hz": [float(x / TWO_PI) for x in corrections.d_omega_rad_s],
    }


def corrections_from_vector(vector: np.ndarray, n_active: int) -> ConditionedMultitoneCorrections:
    arr = np.asarray(vector, dtype=float).reshape(-1)
    n = int(n_active)
    return ConditionedMultitoneCorrections(
        d_lambda=tuple(float(x) for x in arr[0:n]),
        d_alpha=tuple(float(x) for x in arr[n : 2 * n]),
        d_omega_rad_s=tuple(float(x) for x in arr[2 * n : 3 * n]),
    )


def corrections_to_vector(corrections: ConditionedMultitoneCorrections) -> np.ndarray:
    return np.concatenate(
        [
            np.asarray(corrections.d_lambda, dtype=float),
            np.asarray(corrections.d_alpha, dtype=float),
            np.asarray(corrections.d_omega_rad_s, dtype=float),
        ]
    )


def build_multitone_waveform_from_corrections(
    model: DispersiveTransmonCavityModel,
    spec: TargetSpec,
    run_config: ConditionedMultitoneRunConfig,
    *,
    corrections: ConditionedMultitoneCorrections,
    channel: str = "qubit",
    label: str | None = None,
):
    targets = conditioned_targets_from_target_spec(spec)
    tone_specs = build_conditioned_multitone_tones(model, targets, run_config, corrections=corrections)
    waveform = build_conditioned_multitone_waveform(tone_specs, run_config, channel=channel, drive_target="qubit", label=label)
    return waveform, tone_specs


def single_pulse_sequence(
    tone_specs: Sequence[Any],
    run_config: ConditionedMultitoneRunConfig,
    *,
    base_samples: np.ndarray,
    label: str,
) -> tuple[list[Pulse], dict[str, str]]:
    waveform = build_conditioned_multitone_waveform(
        tone_specs,
        run_config,
        base_samples=np.asarray(base_samples, dtype=np.complex128),
        sample_rate=1.0 / float(run_config.dt_s),
        label=label,
    )
    return [waveform.pulse], waveform.drive_ops
