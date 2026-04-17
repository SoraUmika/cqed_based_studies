"""Shared helpers for the strong arbitrary Fock-conditional rotation study."""

from __future__ import annotations

import json
import os
import time
from dataclasses import dataclass
from functools import partial
from pathlib import Path
from typing import Any, Mapping, Sequence

import matplotlib
import numpy as np

matplotlib.use("Agg")

from scipy.linalg import norm, polar, sqrtm
from scipy.optimize import minimize_scalar

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
from cqed_sim.calibration.targeted_subspace_multitone import (
    TargetedSubspaceObjectiveWeights,
    analyze_targeted_subspace_operator,
    build_spanning_state_transfer_set,
)
from cqed_sim.core import DispersiveTransmonCavityModel, FrameSpec
from cqed_sim.core.conventions import qubit_cavity_block_indices
from cqed_sim.core.frequencies import carrier_for_transition_frequency, manifold_transition_frequency
from cqed_sim.core.ideal_gates import qubit_rotation_xy
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
DEFAULT_KERR = 0.0

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

PROBE_QUBIT_STATES: dict[str, np.ndarray] = {
    "g": np.asarray([1.0, 0.0], dtype=np.complex128),
    "e": np.asarray([0.0, 1.0], dtype=np.complex128),
    "plus_x": np.asarray([1.0, 1.0], dtype=np.complex128) / np.sqrt(2.0),
    "minus_x": np.asarray([1.0, -1.0], dtype=np.complex128) / np.sqrt(2.0),
    "plus_y": np.asarray([1.0, 1.0j], dtype=np.complex128) / np.sqrt(2.0),
    "minus_y": np.asarray([1.0, -1.0j], dtype=np.complex128) / np.sqrt(2.0),
}


@dataclass(frozen=True)
class TargetSpec:
    family: str
    target_class: str
    blocks: tuple[np.ndarray, ...]
    seed_theta_values: tuple[float, ...]
    seed_phi_values: tuple[float, ...]
    block_rows: tuple[dict[str, Any], ...]
    metadata: dict[str, Any]


@dataclass(frozen=True)
class ZGaugeFit:
    delta_rad: float
    block_phase_rad: float
    process_fidelity: float
    target_block: np.ndarray


def apply_plot_style() -> None:
    import matplotlib.pyplot as plt

    if STYLE_PATH.exists():
        plt.style.use(str(STYLE_PATH))


def wrap_pi(value: float | np.ndarray) -> float | np.ndarray:
    wrapped = (np.asarray(value, dtype=float) + np.pi) % (2.0 * np.pi) - np.pi
    if np.ndim(wrapped) == 0:
        return float(wrapped)
    return wrapped


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


def duration_from_chi_t(chi_t_over_2pi: float, *, chi: float = CHI) -> float:
    return float(chi_t_over_2pi) / (abs(float(chi)) / TWO_PI)


def build_model(
    *,
    include_chi_prime: bool,
    n_active: int,
    n_tr: int = N_TR,
    n_cav_padding: int = N_CAV_PADDING,
    kerr: float = DEFAULT_KERR,
) -> DispersiveTransmonCavityModel:
    n_cav = int(n_active) + int(n_cav_padding)
    chi_higher = (CHI_PRIME,) if include_chi_prime else ()
    return DispersiveTransmonCavityModel(
        omega_q=OMEGA_Q,
        omega_c=OMEGA_C,
        alpha=ALPHA,
        chi=CHI,
        chi_higher=chi_higher,
        kerr=float(kerr),
        n_cav=int(n_cav),
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


def make_run_config(
    model: DispersiveTransmonCavityModel,
    *,
    n_active: int,
    duration_s: float,
    dt_s: float = DEFAULT_DT,
    sigma_fraction: float = DEFAULT_SIGMA_FRACTION,
) -> ConditionedMultitoneRunConfig:
    return ConditionedMultitoneRunConfig(
        frame=build_frame(model),
        duration_s=float(duration_s),
        dt_s=float(dt_s),
        sigma_fraction=float(sigma_fraction),
        tone_cutoff=1.0e-12,
        include_all_levels=False,
        max_step_s=float(dt_s),
        # Patched cqed_sim expects absolute frequencies when this is overridden.
        # The study uses the package default frequency model directly.
        fock_fqs_hz=None,
    )


def basis_state(model: DispersiveTransmonCavityModel, qubit_level: int, cavity_level: int) -> qt.Qobj:
    return model.basis_state(int(qubit_level), int(cavity_level))


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


def rx(angle: float) -> np.ndarray:
    half = 0.5 * float(angle)
    return np.asarray(
        [[np.cos(half), -1.0j * np.sin(half)], [-1.0j * np.sin(half), np.cos(half)]],
        dtype=np.complex128,
    )


def xy_rotation(theta: float, phi: float) -> np.ndarray:
    return np.asarray(qubit_rotation_xy(float(theta), float(phi)).full(), dtype=np.complex128)


def su2_from_zyz(phi: float, theta: float, lam: float) -> np.ndarray:
    return rz(phi) @ ry(theta) @ rz(lam)


def haar_random_su2(rng: np.random.Generator) -> np.ndarray:
    u1, u2, u3 = rng.random(3)
    a = np.sqrt(1.0 - u1) * np.exp(1.0j * TWO_PI * u2)
    b = np.sqrt(u1) * np.exp(1.0j * TWO_PI * u3)
    return np.asarray([[a, b], [-np.conj(b), np.conj(a)]], dtype=np.complex128)


def nearest_unitary(matrix: np.ndarray) -> np.ndarray:
    unitary, _ = polar(np.asarray(matrix, dtype=np.complex128))
    det = np.linalg.det(unitary)
    if abs(det) > 1.0e-15:
        unitary = unitary * np.exp(-0.5j * np.angle(det))
    return np.asarray(unitary, dtype=np.complex128)


def su2_matrix_sqrt(matrix: np.ndarray) -> np.ndarray:
    return nearest_unitary(np.asarray(sqrtm(np.asarray(matrix, dtype=np.complex128)), dtype=np.complex128))


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


def block_diag_target(blocks: Sequence[np.ndarray]) -> np.ndarray:
    n_blocks = len(blocks)
    operator = np.zeros((2 * n_blocks, 2 * n_blocks), dtype=np.complex128)
    for index, block in enumerate(blocks):
        operator[2 * index : 2 * index + 2, 2 * index : 2 * index + 2] = np.asarray(block, dtype=np.complex128)
    return operator


def conditioned_seed_targets_from_blocks(blocks: Sequence[np.ndarray]) -> ConditionedQubitTargets:
    spec = {
        int(index): block_to_seed_angles(np.asarray(block, dtype=np.complex128))
        for index, block in enumerate(blocks)
    }
    return ConditionedQubitTargets.from_spec(spec, n_levels=len(spec))


def echoed_half_target_blocks(blocks: Sequence[np.ndarray]) -> tuple[np.ndarray, ...]:
    half_blocks: list[np.ndarray] = []
    x_dagger = IDEAL_X_PI.conj().T
    for block in blocks:
        root = su2_matrix_sqrt(np.asarray(block, dtype=np.complex128))
        half_blocks.append(nearest_unitary(x_dagger @ root))
    return tuple(half_blocks)


def _xy_structured_rows(n_active: int) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for level in range(int(n_active)):
        theta = min(0.24 * np.pi + 0.14 * np.pi * level, 0.90 * np.pi)
        phi = 0.0 if level % 2 == 0 else 0.5 * np.pi
        rows.append(
            {
                "level": int(level),
                "rotation_kind": "xy_plane",
                "theta_rad": float(theta),
                "phi_rad": float(phi),
                "lambda_rad": 0.0,
                "block_matrix": xy_rotation(theta, phi),
            }
        )
    return rows


def _inplane_axis_rows(n_active: int) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for level in range(int(n_active)):
        theta = min(0.18 * np.pi + 0.19 * np.pi * level, 0.93 * np.pi)
        phi = float(-0.32 * np.pi + 0.21 * np.pi * level)
        rows.append(
            {
                "level": int(level),
                "rotation_kind": "inplane_axis",
                "theta_rad": float(theta),
                "phi_rad": float(phi),
                "lambda_rad": 0.0,
                "block_matrix": xy_rotation(theta, phi),
            }
        )
    return rows


def _structured_zyz_rows(n_active: int) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for level in range(int(n_active)):
        phi = float(-0.30 * np.pi + 0.16 * np.pi * level)
        theta = float(min(0.18 * np.pi + 0.14 * np.pi * level, 0.92 * np.pi))
        lam = float(0.27 * np.pi - 0.09 * np.pi * level)
        rows.append(
            {
                "level": int(level),
                "rotation_kind": "zyz",
                "theta_rad": theta,
                "phi_rad": phi,
                "lambda_rad": lam,
                "block_matrix": su2_from_zyz(phi, theta, lam),
            }
        )
    return rows


def _stress_zyz_rows(n_active: int) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for level in range(int(n_active)):
        phi = float((0.65 if level % 2 == 0 else -0.55) * np.pi + 0.03 * np.pi * level)
        theta = float(0.82 * np.pi - 0.08 * np.pi * (level % 3))
        lam = float(((-1) ** level) * 0.72 * np.pi + 0.06 * np.pi * level)
        rows.append(
            {
                "level": int(level),
                "rotation_kind": "stress_zyz",
                "theta_rad": theta,
                "phi_rad": phi,
                "lambda_rad": lam,
                "block_matrix": su2_from_zyz(phi, theta, lam),
            }
        )
    return rows


def _random_su2_rows(n_active: int, *, seed: int) -> list[dict[str, Any]]:
    rng = np.random.default_rng(int(seed))
    rows: list[dict[str, Any]] = []
    for level in range(int(n_active)):
        block = haar_random_su2(rng)
        seed_theta, seed_phi = block_to_seed_angles(block)
        rows.append(
            {
                "level": int(level),
                "rotation_kind": "random_su2",
                "theta_rad": float(seed_theta),
                "phi_rad": float(seed_phi),
                "lambda_rad": float("nan"),
                "block_matrix": block,
            }
        )
    return rows


def target_spec(family: str, n_active: int, *, seed: int | None = None) -> TargetSpec:
    if family == "xy_structured":
        rows = _xy_structured_rows(n_active)
        target_class = "pure_xy"
        description = "Pure conditional X/Y-plane rotations with axis restricted to x or y."
    elif family == "inplane_axes":
        rows = _inplane_axis_rows(n_active)
        target_class = "general_inplane"
        description = "General axis-in-plane rotations with varying azimuth."
    elif family == "structured_zyz":
        rows = _structured_zyz_rows(n_active)
        target_class = "structured_su2"
        description = "Structured arbitrary SU(2) blocks in ZYZ form."
    elif family == "stress_zyz":
        rows = _stress_zyz_rows(n_active)
        target_class = "stress_test"
        description = "Alternating, large-angle SU(2) blocks chosen to expose failure modes."
    elif family == "random_su2":
        if seed is None:
            raise ValueError("random_su2 requires a seed.")
        rows = _random_su2_rows(n_active, seed=int(seed))
        target_class = "random_su2"
        description = "Haar-random blockwise SU(2) targets."
    else:
        raise ValueError(f"Unsupported target family '{family}'.")

    blocks = tuple(np.asarray(row["block_matrix"], dtype=np.complex128) for row in rows)
    seeds = [block_to_seed_angles(block) for block in blocks]
    metadata = {
        "family": str(family),
        "target_class": str(target_class),
        "n_active": int(n_active),
        "seed": None if seed is None else int(seed),
        "description": description,
    }
    return TargetSpec(
        family=str(family),
        target_class=str(target_class),
        blocks=blocks,
        seed_theta_values=tuple(float(item[0]) for item in seeds),
        seed_phi_values=tuple(float(item[1]) for item in seeds),
        block_rows=tuple(
            {
                "level": int(row["level"]),
                "rotation_kind": str(row["rotation_kind"]),
                "theta_rad": float(row["theta_rad"]),
                "phi_rad": float(row["phi_rad"]),
                "lambda_rad": None if not np.isfinite(float(row["lambda_rad"])) else float(row["lambda_rad"]),
            }
            for row in rows
        ),
        metadata=metadata,
    )


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


def corrections_to_dict(corrections: ConditionedMultitoneCorrections) -> dict[str, list[float]]:
    return {
        "d_lambda": [float(x) for x in corrections.d_lambda],
        "d_alpha": [float(x) for x in corrections.d_alpha],
        "d_omega_rad_s": [float(x) for x in corrections.d_omega_rad_s],
        "d_omega_hz": [float(x / TWO_PI) for x in corrections.d_omega_rad_s],
    }


def corrections_to_vector(corrections: ConditionedMultitoneCorrections) -> np.ndarray:
    return np.concatenate(
        [
            np.asarray(corrections.d_lambda, dtype=float),
            np.asarray(corrections.d_alpha, dtype=float),
            np.asarray(corrections.d_omega_rad_s, dtype=float),
        ]
    )


def corrections_from_vector(vector: np.ndarray, n_active: int) -> ConditionedMultitoneCorrections:
    arr = np.asarray(vector, dtype=float).reshape(-1)
    n = int(n_active)
    return ConditionedMultitoneCorrections(
        d_lambda=tuple(float(x) for x in arr[0:n]),
        d_alpha=tuple(float(x) for x in arr[n : 2 * n]),
        d_omega_rad_s=tuple(float(x) for x in arr[2 * n : 3 * n]),
    )


def build_multitone_waveform_from_corrections(
    model: DispersiveTransmonCavityModel,
    seed_targets: ConditionedQubitTargets,
    run_config: ConditionedMultitoneRunConfig,
    *,
    corrections: ConditionedMultitoneCorrections,
    channel: str = "qubit",
    label: str | None = None,
) -> tuple[Any, Sequence[Any]]:
    tone_specs = build_conditioned_multitone_tones(model, seed_targets, run_config, corrections=corrections)
    waveform = build_conditioned_multitone_waveform(
        tone_specs,
        run_config,
        channel=channel,
        drive_target="qubit",
        label=label,
    )
    return waveform, tone_specs


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
    drive_ops: Mapping[str, str],
) -> Any:
    config = SimulationConfig(frame=frame, store_states=False)
    return prepare_simulation(model, compiled, dict(drive_ops), config=config, e_ops={})


def simulate_full_operator_on_logical_inputs(
    model: DispersiveTransmonCavityModel,
    compiled: Any,
    *,
    frame: FrameSpec,
    drive_ops: Mapping[str, str],
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
    drive_ops: Mapping[str, str],
    initial_state: qt.Qobj,
) -> qt.Qobj:
    session = _simulation_session_for_compiled(model, compiled, frame=frame, drive_ops=drive_ops)
    return session.run(initial_state).final_state


def channel_waveform_samples(compiled: Any, channel: str = "qubit") -> dict[str, Any]:
    compiled_channel = compiled.channels[str(channel)]
    return {
        "time_s": np.asarray(compiled.tlist, dtype=float).tolist(),
        "baseband_real": np.real(np.asarray(compiled_channel.baseband, dtype=np.complex128)).tolist(),
        "baseband_imag": np.imag(np.asarray(compiled_channel.baseband, dtype=np.complex128)).tolist(),
        "distorted_real": np.real(np.asarray(compiled_channel.distorted, dtype=np.complex128)).tolist(),
        "distorted_imag": np.imag(np.asarray(compiled_channel.distorted, dtype=np.complex128)).tolist(),
    }


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


def operator_2norm_error(target_operator: np.ndarray, actual_operator: np.ndarray) -> float:
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


def unitary_rotation_parameters(unitary: np.ndarray) -> tuple[float, float, float]:
    u = nearest_unitary(unitary)
    trace = np.trace(u)
    cos_half = float(np.clip(np.real(trace / 2.0), -1.0, 1.0))
    theta = float(2.0 * np.arccos(cos_half))
    if theta < 1.0e-12:
        return 0.0, 0.0, 0.0
    sin_half = float(np.sin(theta / 2.0))
    nx = float(np.real(1.0j * np.trace(PAULI_X @ u) / (2.0 * sin_half)))
    ny = float(np.real(1.0j * np.trace(PAULI_Y @ u) / (2.0 * sin_half)))
    nz = float(np.real(1.0j * np.trace(PAULI_Z @ u) / (2.0 * sin_half)))
    axis_norm = float(np.sqrt(nx * nx + ny * ny + nz * nz))
    if axis_norm > 1.0e-12:
        nx /= axis_norm
        ny /= axis_norm
        nz /= axis_norm
    phi = float(np.mod(np.arctan2(ny, nx), TWO_PI))
    return theta, phi, nz


def coherent_error_decomposition(actual_block: np.ndarray, target_block: np.ndarray) -> dict[str, float]:
    actual = nearest_unitary(actual_block)
    target = nearest_unitary(target_block)
    error = nearest_unitary(target.conj().T @ actual)
    trace = np.trace(error)
    cos_half = float(np.clip(np.real(trace / 2.0), -1.0, 1.0))
    theta = float(2.0 * np.arccos(cos_half))
    if theta < 1.0e-12:
        rotvec = np.zeros(3, dtype=float)
    else:
        sin_half = float(np.sin(theta / 2.0))
        nx = float(np.real(1.0j * np.trace(PAULI_X @ error) / (2.0 * sin_half)))
        ny = float(np.real(1.0j * np.trace(PAULI_Y @ error) / (2.0 * sin_half)))
        nz = float(np.real(1.0j * np.trace(PAULI_Z @ error) / (2.0 * sin_half)))
        axis = np.asarray([nx, ny, nz], dtype=float)
        axis_norm = float(np.linalg.norm(axis))
        if axis_norm > 1.0e-12:
            axis = axis / axis_norm
        rotvec = axis * theta
    return {
        "error_rotation_angle_rad": float(np.linalg.norm(rotvec)),
        "residual_z_error_rad": float(abs(rotvec[2])),
        "transverse_error_rad": float(np.linalg.norm(rotvec[:2])),
        "error_rotvec_x_rad": float(rotvec[0]),
        "error_rotvec_y_rad": float(rotvec[1]),
        "error_rotvec_z_rad": float(rotvec[2]),
    }


def block_rotation_metrics(target_block: np.ndarray, actual_block: np.ndarray) -> dict[str, float]:
    target = np.asarray(target_block, dtype=np.complex128)
    actual = np.asarray(actual_block, dtype=np.complex128)
    target_theta, target_phi, target_axis_z = unitary_rotation_parameters(target)
    actual_theta, actual_phi, actual_axis_z = unitary_rotation_parameters(actual)
    angle_error = abs(actual_theta - target_theta)
    if target_theta < 1.0e-12 or actual_theta < 1.0e-12:
        axis_error = float("nan")
    else:
        target_axis = np.asarray([np.cos(target_phi), np.sin(target_phi), target_axis_z], dtype=float)
        actual_axis = np.asarray([np.cos(actual_phi), np.sin(actual_phi), actual_axis_z], dtype=float)
        target_axis = target_axis / max(float(np.linalg.norm(target_axis)), 1.0e-12)
        actual_axis = actual_axis / max(float(np.linalg.norm(actual_axis)), 1.0e-12)
        axis_error = float(np.arccos(np.clip(np.dot(target_axis, actual_axis), -1.0, 1.0)))
    return {
        "process_fidelity": float(process_fidelity(target, actual)),
        "average_gate_fidelity": float(average_gate_fidelity(target, actual)),
        "frobenius_error": float(frobenius_error(target, actual)),
        "operator_2norm_error": float(operator_2norm_error(target, actual)),
        "unitarity_error": float(unitarity_error(actual)),
        "target_rotation_angle_rad": float(target_theta),
        "actual_rotation_angle_rad": float(actual_theta),
        "rotation_angle_error_rad": float(angle_error),
        "rotation_axis_error_rad": float(axis_error),
        "target_axis_phi_rad": float(target_phi),
        "actual_axis_phi_rad": float(actual_phi),
        "target_axis_z": float(target_axis_z),
        "actual_axis_z": float(actual_axis_z),
        **coherent_error_decomposition(actual, target),
    }


def state_density_matrix(qubit_state: Sequence[complex]) -> np.ndarray:
    vec = np.asarray(qubit_state, dtype=np.complex128).reshape((2, 1))
    return np.asarray(vec @ vec.conj().T, dtype=np.complex128)


def state_fidelity_from_dm(actual_rho: np.ndarray, target_rho: np.ndarray) -> float:
    actual = np.asarray(actual_rho, dtype=np.complex128)
    target = np.asarray(target_rho, dtype=np.complex128)
    return float(np.clip(np.real(np.trace(actual @ target)), 0.0, 1.0))


def full_state_fidelity(actual_state: np.ndarray, target_state: np.ndarray) -> float:
    actual = np.asarray(actual_state, dtype=np.complex128).reshape(-1)
    target = np.asarray(target_state, dtype=np.complex128).reshape(-1)
    return float(np.clip(abs(np.vdot(target, actual)) ** 2, 0.0, 1.0))


def bloch_vector_from_density_matrix(rho: np.ndarray) -> tuple[float, float, float]:
    matrix = np.asarray(rho, dtype=np.complex128)
    return (
        float(np.real(np.trace(matrix @ PAULI_X))),
        float(np.real(np.trace(matrix @ PAULI_Y))),
        float(np.real(np.trace(matrix @ PAULI_Z))),
    )


def qubit_channel_kraus_from_full(full_operator: np.ndarray, n_cav: int, input_level: int) -> tuple[np.ndarray, ...]:
    full = np.asarray(full_operator, dtype=np.complex128)
    kraus: list[np.ndarray] = []
    for output_level in range(int(n_cav)):
        block = np.zeros((2, 2), dtype=np.complex128)
        for q_out in range(2):
            for q_in in range(2):
                row = q_out * int(n_cav) + int(output_level)
                col = q_in * int(n_cav) + int(input_level)
                block[q_out, q_in] = full[row, col]
        kraus.append(block)
    return tuple(kraus)


def same_manifold_block(full_operator: np.ndarray, n_cav: int, level: int) -> np.ndarray:
    return np.asarray(qubit_channel_kraus_from_full(full_operator, n_cav, level)[int(level)], dtype=np.complex128)


def apply_channel(kraus_ops: Sequence[np.ndarray], rho: np.ndarray) -> np.ndarray:
    state = np.asarray(rho, dtype=np.complex128)
    out = np.zeros((2, 2), dtype=np.complex128)
    for op in kraus_ops:
        mat = np.asarray(op, dtype=np.complex128)
        out = out + mat @ state @ mat.conj().T
    return out


def channel_process_fidelity_to_unitary(kraus_ops: Sequence[np.ndarray], target_unitary: np.ndarray) -> float:
    target = np.asarray(target_unitary, dtype=np.complex128)
    accum = 0.0
    for op in kraus_ops:
        accum += abs(np.trace(target.conj().T @ np.asarray(op, dtype=np.complex128))) ** 2
    return float(np.clip(accum / 4.0, 0.0, 1.0))


def fit_zgauge_block(actual_block: np.ndarray, strict_target_block: np.ndarray) -> ZGaugeFit:
    actual = np.asarray(actual_block, dtype=np.complex128)
    strict_target = np.asarray(strict_target_block, dtype=np.complex128)

    def overlap_abs(delta_rad: float) -> float:
        target = rz(delta_rad) @ strict_target
        return float(abs(np.trace(target.conj().T @ actual)))

    result = minimize_scalar(
        lambda delta: -overlap_abs(delta),
        bounds=(-np.pi, np.pi),
        method="bounded",
        options={"xatol": 1.0e-4},
    )
    delta = float(result.x)
    base_target = rz(delta) @ strict_target
    term = np.trace(base_target.conj().T @ actual)
    block_phase = float(np.angle(term))
    target = np.exp(1.0j * block_phase) * base_target
    return ZGaugeFit(
        delta_rad=delta,
        block_phase_rad=block_phase,
        process_fidelity=process_fidelity(target, actual),
        target_block=target,
    )


def fit_zgauge_channel(kraus_ops: Sequence[np.ndarray], strict_target_block: np.ndarray) -> ZGaugeFit:
    strict_target = np.asarray(strict_target_block, dtype=np.complex128)

    def objective(delta_rad: float) -> float:
        target = rz(delta_rad) @ strict_target
        return -channel_process_fidelity_to_unitary(kraus_ops, target)

    result = minimize_scalar(
        objective,
        bounds=(-np.pi, np.pi),
        method="bounded",
        options={"xatol": 1.0e-4},
    )
    delta = float(result.x)
    base_target = rz(delta) @ strict_target
    return ZGaugeFit(
        delta_rad=delta,
        block_phase_rad=0.0,
        process_fidelity=channel_process_fidelity_to_unitary(kraus_ops, base_target),
        target_block=base_target,
    )


def build_zgauge_joint_target(
    restricted_operator: np.ndarray,
    strict_blocks: Sequence[np.ndarray],
) -> tuple[np.ndarray, list[dict[str, float]]]:
    restricted = np.asarray(restricted_operator, dtype=np.complex128)
    target = np.zeros_like(restricted)
    rows: list[dict[str, float]] = []
    for idx, strict_block in enumerate(strict_blocks):
        block = restricted[2 * idx : 2 * idx + 2, 2 * idx : 2 * idx + 2]
        fit = fit_zgauge_block(block, strict_block)
        target[2 * idx : 2 * idx + 2, 2 * idx : 2 * idx + 2] = fit.target_block
        rows.append(
            {
                "level": int(idx),
                "delta_rad": float(fit.delta_rad),
                "block_phase_rad": float(fit.block_phase_rad),
                "process_fidelity": float(fit.process_fidelity),
            }
        )
    return target, rows


def offblock_norm_rows(restricted_operator: np.ndarray) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    blocks = restricted_blocks(restricted_operator)
    for out_idx in range(len(blocks)):
        for in_idx in range(len(blocks)):
            block = np.asarray(restricted_operator[2 * out_idx : 2 * out_idx + 2, 2 * in_idx : 2 * in_idx + 2], dtype=np.complex128)
            rows.append(
                {
                    "out_level": int(out_idx),
                    "in_level": int(in_idx),
                    "fro_norm": float(norm(block, ord="fro")),
                    "operator_2norm": float(norm(block, ord=2)),
                    "is_diagonal": bool(out_idx == in_idx),
                }
            )
    return rows


def addressed_indices(n_cav: int, levels: Sequence[int]) -> np.ndarray:
    rows: list[int] = []
    for level in levels:
        rows.extend([int(level), int(n_cav) + int(level)])
    return np.asarray(rows, dtype=int)


def leakage_outside_indices(state_vector: np.ndarray, keep_indices: Sequence[int]) -> float:
    vec = np.asarray(state_vector, dtype=np.complex128).reshape(-1)
    keep = np.asarray(tuple(int(index) for index in keep_indices), dtype=int)
    keep_prob = float(np.sum(np.abs(vec[keep]) ** 2))
    return float(max(0.0, 1.0 - keep_prob))


def reduced_probe_rows_from_full(
    full_operator: np.ndarray,
    strict_blocks: Sequence[np.ndarray],
    n_cav: int,
    levels: Sequence[int],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    reduced_rows: list[dict[str, Any]] = []
    level_rows: list[dict[str, Any]] = []
    for idx, (level, strict_block) in enumerate(zip(levels, strict_blocks, strict=True)):
        kraus_ops = qubit_channel_kraus_from_full(full_operator, int(n_cav), int(level))
        same_block = same_manifold_block(full_operator, int(n_cav), int(level))
        dominant_unitary = nearest_unitary(same_block)
        relaxed_fit_channel = fit_zgauge_channel(kraus_ops, strict_block)
        relaxed_fit_block = fit_zgauge_block(dominant_unitary, strict_block)
        theta_actual, phi_actual, axis_z_actual = unitary_rotation_parameters(dominant_unitary)
        theta_target, phi_target, axis_z_target = unitary_rotation_parameters(strict_block)
        level_rows.append(
            {
                "level": int(level),
                "strict_process_fidelity": float(process_fidelity(strict_block, dominant_unitary)),
                "relaxed_process_fidelity": float(relaxed_fit_channel.process_fidelity),
                "strict_target_theta_rad": float(theta_target),
                "strict_target_phi_rad": float(phi_target),
                "strict_target_axis_z": float(axis_z_target),
                "actual_theta_rad": float(theta_actual),
                "actual_phi_rad": float(phi_actual),
                "actual_axis_z": float(axis_z_actual),
                "relaxed_delta_rad": float(relaxed_fit_channel.delta_rad),
                **{f"strict_{key}": float(value) for key, value in coherent_error_decomposition(dominant_unitary, strict_block).items()},
                **{f"relaxed_{key}": float(value) for key, value in coherent_error_decomposition(dominant_unitary, relaxed_fit_block.target_block).items()},
            }
        )
        for probe_name, qubit_state in PROBE_QUBIT_STATES.items():
            rho_in = state_density_matrix(qubit_state)
            actual = apply_channel(kraus_ops, rho_in)
            strict_target = strict_block @ np.asarray(qubit_state, dtype=np.complex128).reshape((2, 1))
            relaxed_target = relaxed_fit_channel.target_block @ np.asarray(qubit_state, dtype=np.complex128).reshape((2, 1))
            strict_rho = strict_target @ strict_target.conj().T
            relaxed_rho = relaxed_target @ relaxed_target.conj().T
            reduced_rows.append(
                {
                    "level": int(level),
                    "probe_label": str(probe_name),
                    "strict_fidelity": state_fidelity_from_dm(actual, strict_rho),
                    "relaxed_fidelity": state_fidelity_from_dm(actual, relaxed_rho),
                    "bloch_x": bloch_vector_from_density_matrix(actual)[0],
                    "bloch_y": bloch_vector_from_density_matrix(actual)[1],
                    "bloch_z": bloch_vector_from_density_matrix(actual)[2],
                    "ideal_strict_bloch_x": bloch_vector_from_density_matrix(strict_rho)[0],
                    "ideal_strict_bloch_y": bloch_vector_from_density_matrix(strict_rho)[1],
                    "ideal_strict_bloch_z": bloch_vector_from_density_matrix(strict_rho)[2],
                    "ideal_relaxed_bloch_x": bloch_vector_from_density_matrix(relaxed_rho)[0],
                    "ideal_relaxed_bloch_y": bloch_vector_from_density_matrix(relaxed_rho)[1],
                    "ideal_relaxed_bloch_z": bloch_vector_from_density_matrix(relaxed_rho)[2],
                }
            )
    return reduced_rows, level_rows


def probe_tier_summary(probe_rows: Sequence[dict[str, Any]], metric_key: str) -> dict[str, dict[str, float]]:
    tiers = {
        "basis_pair": ("g", "e"),
        "cartesian_superpositions": ("plus_x", "minus_x", "plus_y", "minus_y"),
        "full_six_state": tuple(PROBE_QUBIT_STATES.keys()),
    }
    summary: dict[str, dict[str, float]] = {}
    for tier_name, labels in tiers.items():
        values = [float(row[metric_key]) for row in probe_rows if str(row["probe_label"]) in labels]
        summary[tier_name] = {
            "mean_fidelity": float(np.mean(values)) if values else float("nan"),
            "min_fidelity": float(np.min(values)) if values else float("nan"),
        }
    return summary


def full_probe_rows_from_full(
    full_operator: np.ndarray,
    strict_blocks: Sequence[np.ndarray],
    relaxed_blocks: Sequence[np.ndarray],
    n_cav: int,
    levels: Sequence[int],
) -> list[dict[str, Any]]:
    full = np.asarray(full_operator, dtype=np.complex128)
    keep_indices = addressed_indices(n_cav, levels)
    rows: list[dict[str, Any]] = []
    for idx, level in enumerate(levels):
        strict_block = np.asarray(strict_blocks[idx], dtype=np.complex128)
        relaxed_block = np.asarray(relaxed_blocks[idx], dtype=np.complex128)
        for probe_name, qubit_state in PROBE_QUBIT_STATES.items():
            psi0 = np.zeros(2 * int(n_cav), dtype=np.complex128)
            psi0[int(level)] = complex(qubit_state[0])
            psi0[int(n_cav) + int(level)] = complex(qubit_state[1])
            actual = full @ psi0
            strict_target = np.zeros_like(actual)
            relaxed_target = np.zeros_like(actual)
            strict_target[int(level)] = strict_block[0, 0] * qubit_state[0] + strict_block[0, 1] * qubit_state[1]
            strict_target[int(n_cav) + int(level)] = strict_block[1, 0] * qubit_state[0] + strict_block[1, 1] * qubit_state[1]
            relaxed_target[int(level)] = relaxed_block[0, 0] * qubit_state[0] + relaxed_block[0, 1] * qubit_state[1]
            relaxed_target[int(n_cav) + int(level)] = relaxed_block[1, 0] * qubit_state[0] + relaxed_block[1, 1] * qubit_state[1]
            same_manifold = float(abs(actual[int(level)]) ** 2 + abs(actual[int(n_cav) + int(level)]) ** 2)
            rows.append(
                {
                    "level": int(level),
                    "probe_label": str(probe_name),
                    "strict_fidelity": full_state_fidelity(actual, strict_target),
                    "relaxed_fidelity": full_state_fidelity(actual, relaxed_target),
                    "same_manifold_population": float(same_manifold),
                    "addressed_window_population": float(np.sum(np.abs(actual[keep_indices]) ** 2)),
                    "leakage_outside_addressed": leakage_outside_indices(actual, keep_indices),
                }
            )
    return rows


def cross_block_superposition_rows(
    full_operator: np.ndarray,
    strict_blocks: Sequence[np.ndarray],
    relaxed_blocks: Sequence[np.ndarray],
    n_cav: int,
    levels: Sequence[int],
) -> list[dict[str, Any]]:
    if len(levels) < 2:
        return []
    full = np.asarray(full_operator, dtype=np.complex128)
    keep_indices = addressed_indices(n_cav, levels)
    pair_set: list[tuple[int, int]] = [(int(levels[0]), int(levels[1]))]
    if len(levels) >= 3:
        pair_set.append((int(levels[0]), int(levels[-1])))
    rows: list[dict[str, Any]] = []
    for left, right in pair_set:
        psi0 = np.zeros(2 * int(n_cav), dtype=np.complex128)
        psi0[int(left)] = 1.0 / np.sqrt(2.0)
        psi0[int(right)] = 1.0 / np.sqrt(2.0)
        actual = full @ psi0
        strict_target = np.zeros_like(actual)
        relaxed_target = np.zeros_like(actual)
        idx_left = tuple(levels).index(int(left))
        idx_right = tuple(levels).index(int(right))
        strict_target[int(left)] = strict_blocks[idx_left][0, 0] / np.sqrt(2.0)
        strict_target[int(n_cav) + int(left)] = strict_blocks[idx_left][1, 0] / np.sqrt(2.0)
        strict_target[int(right)] = strict_blocks[idx_right][0, 0] / np.sqrt(2.0)
        strict_target[int(n_cav) + int(right)] = strict_blocks[idx_right][1, 0] / np.sqrt(2.0)
        relaxed_target[int(left)] = relaxed_blocks[idx_left][0, 0] / np.sqrt(2.0)
        relaxed_target[int(n_cav) + int(left)] = relaxed_blocks[idx_left][1, 0] / np.sqrt(2.0)
        relaxed_target[int(right)] = relaxed_blocks[idx_right][0, 0] / np.sqrt(2.0)
        relaxed_target[int(n_cav) + int(right)] = relaxed_blocks[idx_right][1, 0] / np.sqrt(2.0)
        rows.append(
            {
                "left_level": int(left),
                "right_level": int(right),
                "probe_label": "g_cross_superposition",
                "strict_fidelity": full_state_fidelity(actual, strict_target),
                "relaxed_fidelity": full_state_fidelity(actual, relaxed_target),
                "addressed_window_population": float(np.sum(np.abs(actual[keep_indices]) ** 2)),
                "leakage_outside_addressed": leakage_outside_indices(actual, keep_indices),
            }
        )
    return rows


def summed_target_rows(strict_blocks: Sequence[np.ndarray], relaxed_blocks: Sequence[np.ndarray]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for idx, (strict_block, relaxed_block) in enumerate(zip(strict_blocks, relaxed_blocks, strict=True)):
        rows.append(
            {
                "level": int(idx),
                "strict_block": np.asarray(strict_block, dtype=np.complex128),
                "relaxed_block": np.asarray(relaxed_block, dtype=np.complex128),
            }
        )
    return rows


def state_validation_summary_for_compiled(
    model: DispersiveTransmonCavityModel,
    compiled: Any,
    *,
    frame: FrameSpec,
    drive_ops: Mapping[str, str],
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
        rows.append({"label": label, "state_fidelity": float(fidelity)})
    return {"states": rows}


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


def build_analysis_bundle(
    model: DispersiveTransmonCavityModel,
    seed_targets: ConditionedQubitTargets,
    strict_blocks: Sequence[np.ndarray],
    levels: Sequence[int],
    full_operator: np.ndarray,
    *,
    objective_weights: TargetedSubspaceObjectiveWeights,
    metadata: dict[str, Any] | None = None,
) -> dict[str, Any]:
    target_operator = block_diag_target(strict_blocks)
    strict_validation = analyze_targeted_subspace_operator(
        full_operator,
        model,
        seed_targets,
        logical_levels=levels,
        target_operator=target_operator,
        transfer_set=build_spanning_state_transfer_set(target_operator),
        objective_weights=objective_weights,
        metadata={} if metadata is None else dict(metadata),
    )
    restricted = restricted_operator_from_full(full_operator, model, levels)
    relaxed_target, relaxed_fit_rows = build_zgauge_joint_target(restricted, strict_blocks)
    reduced_probe_rows, reduced_level_rows = reduced_probe_rows_from_full(full_operator, strict_blocks, int(model.n_cav), levels)
    full_probe_rows = full_probe_rows_from_full(
        full_operator,
        strict_blocks,
        restricted_blocks(relaxed_target),
        int(model.n_cav),
        levels,
    )
    cross_rows = cross_block_superposition_rows(
        full_operator,
        strict_blocks,
        restricted_blocks(relaxed_target),
        int(model.n_cav),
        levels,
    )
    return {
        "target_operator": target_operator,
        "strict_validation": strict_validation,
        "restricted_operator": restricted,
        "relaxed_target_operator": relaxed_target,
        "relaxed_fit_rows": relaxed_fit_rows,
        "reduced_probe_rows": reduced_probe_rows,
        "reduced_level_rows": reduced_level_rows,
        "full_probe_rows": full_probe_rows,
        "cross_block_rows": cross_rows,
        "offblock_rows": offblock_norm_rows(restricted),
    }
