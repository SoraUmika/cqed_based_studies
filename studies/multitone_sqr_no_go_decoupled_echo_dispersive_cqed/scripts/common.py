"""Common helpers for the strict no-detuning multitone-SQR study."""

from __future__ import annotations

import json
import math
import os
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

import matplotlib
import numpy as np

matplotlib.use("Agg")

from scipy.linalg import expm, norm, polar
from scipy.optimize import Bounds, minimize
from scipy.spatial.transform import Rotation

import runtime_compat  # noqa: F401
import qutip as qt

from cqed_sim.calibration.conditioned_multitone import (
    ConditionedMultitoneCorrections,
    ConditionedMultitoneRunConfig,
    ConditionedQubitTargets,
    build_conditioned_multitone_tones,
    build_conditioned_multitone_waveform,
)
from cqed_sim.calibration.targeted_subspace_multitone import (
    TargetedSubspaceObjectiveWeights,
    build_block_rotation_target_operator,
    build_spanning_state_transfer_set,
    evaluate_targeted_subspace_multitone,
)
from cqed_sim.core import DispersiveTransmonCavityModel, FrameSpec
from cqed_sim.core.conventions import qubit_cavity_block_indices
from cqed_sim.core.frequencies import carrier_for_transition_frequency, manifold_transition_frequency
from cqed_sim.core.ideal_gates import qubit_rotation_xy
from cqed_sim.pulses.envelopes import gaussian_area_fraction, gaussian_envelope, square_envelope
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
DEFAULT_DT = 2.0e-9
PI_PULSE_DURATION_S = 40.0e-9
PI_PULSE_SIGMA_FRACTION = 0.25

PAULI_X = np.asarray([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
PAULI_Y = np.asarray([[0.0, -1.0j], [1.0j, 0.0]], dtype=np.complex128)
PAULI_Z = np.asarray([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)
PAULIS = (PAULI_X, PAULI_Y, PAULI_Z)
IDENTITY_2 = np.eye(2, dtype=np.complex128)
IDEAL_X_PI = np.asarray(qubit_rotation_xy(np.pi, 0.0).full(), dtype=np.complex128)
_SIGMA_PLUS = qt.create(2)
_SIGMA_MINUS = qt.destroy(2)
_N_Q = qt.num(2)

OBJECTIVE_WEIGHTS = TargetedSubspaceObjectiveWeights(
    qubit_weight=0.2,
    subspace_weight=1.0,
    preservation_weight=0.35,
    leakage_weight=0.35,
)


@dataclass(frozen=True)
class TargetSpec:
    family: str
    theta_values: tuple[float, ...]
    phi_values: tuple[float, ...]
    metadata: dict[str, Any]


@dataclass(frozen=True)
class CaseRequest:
    model_variant: str
    include_chi_prime: bool
    family: str
    n_active: int
    chi_t_over_2pi: float
    seed: int | None = None

    @property
    def case_id(self) -> str:
        seed_part = "" if self.seed is None else f"_seed{int(self.seed)}"
        duration_label = str(self.chi_t_over_2pi).replace(".", "p")
        return f"{self.model_variant}_{self.family}_na{int(self.n_active)}_chiT{duration_label}{seed_part}"


@dataclass
class OptimizationResult:
    corrections: ConditionedMultitoneCorrections
    validation: Any
    tone_specs: tuple[Any, ...]
    waveform: Any
    history: list[dict[str, float]]
    success: bool
    message: str


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


def save_json(path: Path, payload: Any, *, description: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    wrapped = {
        "study_name": STUDY_DIR.name,
        "date_created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "description": str(description),
        "load_instructions": f"Load `{path.name}` as JSON and inspect the documented result rows and metadata.",
        **(payload if isinstance(payload, dict) else {"payload": payload}),
    }
    text = json.dumps(json_ready(wrapped), indent=2, sort_keys=True)
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


def manifold_transition_frequencies_rad_s(
    model: DispersiveTransmonCavityModel,
    levels: Sequence[int],
    frame: FrameSpec,
) -> np.ndarray:
    return np.asarray(
        [float(manifold_transition_frequency(model, int(level), frame=frame)) for level in levels],
        dtype=float,
    )


def make_run_config(
    model: DispersiveTransmonCavityModel,
    *,
    n_active: int,
    duration_s: float,
    dt_s: float = DEFAULT_DT,
) -> ConditionedMultitoneRunConfig:
    frame = build_frame(model)
    levels = logical_levels(n_active)
    return ConditionedMultitoneRunConfig(
        frame=frame,
        duration_s=float(duration_s),
        dt_s=float(dt_s),
        sigma_fraction=1.0 / 6.0,
        tone_cutoff=1.0e-12,
        include_all_levels=False,
        max_step_s=float(dt_s),
        fock_fqs_hz=tuple(
            float(manifold_transition_frequency(model, int(level), frame=frame) / TWO_PI)
            for level in levels
        ),
    )


def target_spec(family: str, n_active: int, *, seed: int | None = None) -> TargetSpec:
    n = int(n_active)
    if family == "aligned_x":
        theta = tuple(float(min(0.28 * np.pi + 0.18 * np.pi * level, 0.92 * np.pi)) for level in range(n))
        phi = tuple(0.0 for _ in range(n))
        metadata = {
            "description": "Aligned x-axis SQR family that gives the echo sequence its best-case symmetry setting.",
            "formula": "theta_n = min(0.28 pi + 0.18 pi n, 0.92 pi), phi_n = 0",
        }
    elif family == "structured_xy":
        theta = tuple(float(min(0.26 * np.pi + 0.13 * np.pi * level, 0.88 * np.pi)) for level in range(n))
        phi = tuple(float((0.20 * np.pi + 0.37 * np.pi * level) % TWO_PI) for level in range(n))
        metadata = {
            "description": "Structured XY family with nontrivial block-dependent azimuths.",
            "formula": "theta_n = min(0.26 pi + 0.13 pi n, 0.88 pi), phi_n = (0.20 + 0.37 n) pi mod 2 pi",
        }
    elif family == "random_xy":
        if seed is None:
            raise ValueError("random_xy family requires a seed.")
        rng = np.random.default_rng(int(seed))
        theta = tuple(float(rng.uniform(0.18 * np.pi, 0.82 * np.pi)) for _ in range(n))
        phi = tuple(float(rng.uniform(0.0, TWO_PI)) for _ in range(n))
        metadata = {
            "description": "Random XY family used for falsification pressure.",
            "seed": int(seed),
        }
    else:
        raise ValueError(f"Unsupported target family '{family}'.")
    metadata["family"] = str(family)
    metadata["n_active"] = int(n)
    metadata["theta_values_rad"] = list(theta)
    metadata["phi_values_rad"] = list(phi)
    return TargetSpec(family=str(family), theta_values=theta, phi_values=phi, metadata=metadata)


def conditioned_targets_from_target_spec(spec: TargetSpec) -> ConditionedQubitTargets:
    pairs = list(zip(spec.theta_values, spec.phi_values, strict=True))
    return ConditionedQubitTargets.from_spec(pairs, n_levels=len(pairs))


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


def axis_angle_unitary(rotvec: np.ndarray) -> np.ndarray:
    vec = np.asarray(rotvec, dtype=float).reshape(3)
    angle = float(np.linalg.norm(vec))
    if angle <= 1.0e-15:
        return np.array(IDENTITY_2, copy=True)
    nx, ny, nz = vec / angle
    generator = nx * PAULI_X + ny * PAULI_Y + nz * PAULI_Z
    return np.asarray(expm(-0.5j * angle * generator), dtype=np.complex128)


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
    return {
        "process_fidelity": float(proc),
        "average_gate_fidelity": float(avg),
        "frobenius_error": float(frobenius_error(target, actual)),
        "unitarity_error": float(unitarity_error(actual)),
        "target_rotation_angle_rad": float(target_angle),
        "actual_rotation_angle_rad": float(actual_angle),
        "rotation_angle_error_rad": float(abs(actual_angle - target_angle)),
        "rotation_axis_error_rad": float(axis_error),
        "actual_rotvec_x_rad": float(actual_vec[0]),
        "actual_rotvec_y_rad": float(actual_vec[1]),
        "actual_rotvec_z_rad": float(actual_vec[2]),
        "error_rotvec_x_rad": float(error_rotvec[0]),
        "error_rotvec_y_rad": float(error_rotvec[1]),
        "error_rotvec_z_rad": float(error_rotvec[2]),
        "residual_z_error_rad": float(abs(error_rotvec[2])),
        "transverse_error_rad": float(np.linalg.norm(error_rotvec[:2])),
        "error_rotation_angle_rad": float(np.linalg.norm(error_rotvec)),
    }


def corrections_to_vector(corrections: ConditionedMultitoneCorrections, *, n_active: int) -> np.ndarray:
    corr = corrections.padded(int(n_active))
    return np.concatenate(
        [
            np.asarray(corr.d_lambda[:n_active], dtype=float),
            np.asarray(corr.d_alpha[:n_active], dtype=float),
        ]
    )


def corrections_from_vector(vector: np.ndarray, *, n_active: int) -> ConditionedMultitoneCorrections:
    arr = np.asarray(vector, dtype=float).reshape(-1)
    n = int(n_active)
    if arr.size != 2 * n:
        raise ValueError(f"Expected vector of length {2 * n}, received {arr.size}.")
    zeros = tuple(0.0 for _ in range(n))
    return ConditionedMultitoneCorrections(
        d_lambda=tuple(float(x) for x in arr[:n]),
        d_alpha=tuple(float(x) for x in arr[n:]),
        d_omega_rad_s=zeros,
    )


def build_square_multitone_waveform(
    model: DispersiveTransmonCavityModel,
    spec: TargetSpec,
    run_config: ConditionedMultitoneRunConfig,
    *,
    corrections: ConditionedMultitoneCorrections | None = None,
    label: str | None = None,
) -> tuple[Any, tuple[Any, ...]]:
    targets = conditioned_targets_from_target_spec(spec)
    corr = ConditionedMultitoneCorrections.zeros(targets.n_levels) if corrections is None else corrections.padded(targets.n_levels)
    tone_specs = tuple(build_conditioned_multitone_tones(model, targets, run_config, corrections=corr))
    waveform = build_conditioned_multitone_waveform(
        tone_specs,
        run_config,
        base_envelope=square_envelope,
        channel="qubit",
        drive_target="qubit",
        label=label,
    )
    return waveform, tone_specs


def evaluate_square_multitone(
    model: DispersiveTransmonCavityModel,
    spec: TargetSpec,
    run_config: ConditionedMultitoneRunConfig,
    *,
    corrections: ConditionedMultitoneCorrections | None = None,
    objective_weights: TargetedSubspaceObjectiveWeights = OBJECTIVE_WEIGHTS,
    label: str | None = None,
) -> tuple[Any, tuple[Any, ...], Any]:
    levels = logical_levels(len(spec.theta_values))
    targets = conditioned_targets_from_target_spec(spec)
    target_operator = build_target_operator(spec, levels)
    transfer_set = build_spanning_state_transfer_set(target_operator)
    waveform, tone_specs = build_square_multitone_waveform(
        model,
        spec,
        run_config,
        corrections=corrections,
        label=label,
    )
    validation = evaluate_targeted_subspace_multitone(
        model,
        targets,
        waveform,
        run_config,
        corrections=ConditionedMultitoneCorrections.zeros(targets.n_levels) if corrections is None else corrections,
        logical_levels=levels,
        objective_weights=objective_weights,
        target_operator=target_operator,
        transfer_set=transfer_set,
    )
    return validation, tone_specs, waveform


def optimize_square_multitone(
    model: DispersiveTransmonCavityModel,
    spec: TargetSpec,
    run_config: ConditionedMultitoneRunConfig,
    *,
    n_starts: int = 4,
    maxiter: int = 90,
    objective_weights: TargetedSubspaceObjectiveWeights = OBJECTIVE_WEIGHTS,
    random_seed: int = 1234,
    label_prefix: str = "square",
) -> OptimizationResult:
    n_active = len(spec.theta_values)
    rng = np.random.default_rng(int(random_seed))
    starts = [np.zeros(2 * n_active, dtype=float)]
    for _ in range(max(0, int(n_starts) - 1)):
        starts.append(
            np.concatenate(
                [
                    rng.uniform(-0.22, 0.22, size=n_active),
                    rng.uniform(-0.40, 0.40, size=n_active),
                ]
            )
        )
    bounds = Bounds(
        np.concatenate([np.full(n_active, -0.75, dtype=float), np.full(n_active, -np.pi, dtype=float)]),
        np.concatenate([np.full(n_active, 0.75, dtype=float), np.full(n_active, np.pi, dtype=float)]),
    )
    history: list[dict[str, float]] = []
    best_result: OptimizationResult | None = None
    best_loss = float("inf")

    for start_index, x0 in enumerate(starts):

        def objective(vector: np.ndarray) -> float:
            corrections = corrections_from_vector(vector, n_active=n_active)
            validation, _tone_specs, _waveform = evaluate_square_multitone(
                model,
                spec,
                run_config,
                corrections=corrections,
                objective_weights=objective_weights,
                label=f"{label_prefix}_start{start_index}",
            )
            loss = float(validation.weighted_loss)
            history.append(
                {
                    "start": float(start_index),
                    "loss": float(loss),
                    "restricted_process_fidelity": float(validation.restricted_process_fidelity),
                    "best_fit_restricted_process_fidelity": float(validation.best_fit_restricted_process_fidelity),
                }
            )
            return loss

        result = minimize(
            objective,
            x0=np.asarray(x0, dtype=float),
            method="Powell",
            bounds=bounds,
            options={"maxiter": int(maxiter), "xtol": 1.0e-3, "ftol": 1.0e-4},
        )
        corrections = corrections_from_vector(result.x, n_active=n_active)
        validation, tone_specs, waveform = evaluate_square_multitone(
            model,
            spec,
            run_config,
            corrections=corrections,
            objective_weights=objective_weights,
            label=f"{label_prefix}_best{start_index}",
        )
        loss = float(validation.weighted_loss)
        if loss < best_loss:
            best_loss = loss
            best_result = OptimizationResult(
                corrections=corrections,
                validation=validation,
                tone_specs=tuple(tone_specs),
                waveform=waveform,
                history=list(history),
                success=bool(result.success),
                message=str(result.message),
            )
    if best_result is None:
        raise RuntimeError("Square multitone optimization did not produce a result.")
    return best_result


def magnus_z_kernel(delta_rad_s: float, duration_s: float) -> float:
    delta = float(delta_rad_s)
    duration = float(duration_s)
    if abs(delta) <= 1.0e-12:
        return 0.0
    return float((1.0 / delta) - np.sin(delta * duration) / (delta * delta * duration))


def magnus_effective_blocks(
    tone_specs: Sequence[Any],
    *,
    model: DispersiveTransmonCavityModel,
    run_config: ConditionedMultitoneRunConfig,
    levels: Sequence[int],
) -> tuple[np.ndarray, ...]:
    frame = run_config.frame
    freqs = manifold_transition_frequencies_rad_s(model, levels, frame)
    amplitude_map = {int(spec.manifold): float(spec.amp_rad_s) for spec in tone_specs}
    phase_map = {int(spec.manifold): float(spec.phase_rad) for spec in tone_specs}
    blocks: list[np.ndarray] = []
    for row, level in enumerate(levels):
        lam = float(amplitude_map.get(int(level), 0.0))
        phi = float(phase_map.get(int(level), 0.0))
        z_coeff = 0.0
        for col, other in enumerate(levels):
            if int(other) == int(level):
                continue
            lam_other = float(amplitude_map.get(int(other), 0.0))
            delta = float(freqs[col] - freqs[row])
            z_coeff -= lam_other * lam_other * magnus_z_kernel(delta, run_config.duration_s)
        rotvec = np.asarray(
            [
                2.0 * run_config.duration_s * lam * np.cos(phi),
                2.0 * run_config.duration_s * lam * np.sin(phi),
                2.0 * run_config.duration_s * z_coeff,
            ],
            dtype=float,
        )
        blocks.append(axis_angle_unitary(rotvec))
    return tuple(blocks)


def decoupled_block_operator(
    tone_specs: Sequence[Any],
    *,
    levels: Sequence[int],
    duration_s: float,
) -> np.ndarray:
    dim = 2 * len(levels)
    operator = np.zeros((dim, dim), dtype=np.complex128)
    amplitude_map = {int(spec.manifold): float(spec.amp_rad_s) for spec in tone_specs}
    phase_map = {int(spec.manifold): float(spec.phase_rad) for spec in tone_specs}
    for block_index, level in enumerate(levels):
        lam = float(amplitude_map.get(int(level), 0.0))
        phi = float(phase_map.get(int(level), 0.0))
        operator[2 * block_index : 2 * block_index + 2, 2 * block_index : 2 * block_index + 2] = np.asarray(
            qubit_rotation_xy(2.0 * lam * float(duration_s), phi).full(),
            dtype=np.complex128,
        )
    return operator


def reduced_sector_unitary(
    model: DispersiveTransmonCavityModel,
    compiled: Any,
    waveform: Any,
    run_config: ConditionedMultitoneRunConfig,
    *,
    level: int,
) -> np.ndarray:
    coeff = np.asarray(compiled.channels[waveform.drive_channel].distorted, dtype=np.complex128)
    detuning = float(manifold_transition_frequency(model, int(level), frame=run_config.frame))
    hamiltonian = [
        detuning * _N_Q,
        [_SIGMA_PLUS, coeff],
        [_SIGMA_MINUS, np.conj(coeff)],
    ]
    options: dict[str, Any] = {
        "atol": 1.0e-8,
        "rtol": 1.0e-7,
        "nsteps": 100000,
    }
    if run_config.max_step_s is not None:
        options["max_step"] = float(run_config.max_step_s)
    propagators = qt.propagator(
        hamiltonian,
        compiled.tlist,
        options=options,
        tlist=compiled.tlist,
    )
    final = propagators[-1] if isinstance(propagators, list) else propagators
    return np.asarray(final.full(), dtype=np.complex128)


def reduced_blockwise_operator(
    model: DispersiveTransmonCavityModel,
    compiled: Any,
    waveform: Any,
    run_config: ConditionedMultitoneRunConfig,
    *,
    levels: Sequence[int],
) -> np.ndarray:
    dim = 2 * len(levels)
    operator = np.zeros((dim, dim), dtype=np.complex128)
    for block_index, level in enumerate(levels):
        operator[2 * block_index : 2 * block_index + 2, 2 * block_index : 2 * block_index + 2] = reduced_sector_unitary(
            model,
            compiled,
            waveform,
            run_config,
            level=int(level),
        )
    return operator


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

    def envelope(t_rel: np.ndarray) -> np.ndarray:
        return gaussian_envelope(np.asarray(t_rel, dtype=float), sigma=float(sigma_fraction), center=0.5) / area_fraction

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
            psi0 = model.basis_state(qubit_level, int(level))
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


def embed_qubit_operator(qubit_operator: np.ndarray, *, n_cav: int) -> np.ndarray:
    return np.asarray(qt.tensor(qt.Qobj(qubit_operator), qt.qeye(int(n_cav))).full(), dtype=np.complex128)


def channel_waveform_samples(compiled: Any, channel: str = "qubit") -> dict[str, Any]:
    compiled_channel = compiled.channels[str(channel)]
    return {
        "time_s": np.asarray(compiled.tlist, dtype=float).tolist(),
        "baseband_real": np.real(np.asarray(compiled_channel.baseband, dtype=np.complex128)).tolist(),
        "baseband_imag": np.imag(np.asarray(compiled_channel.baseband, dtype=np.complex128)).tolist(),
        "distorted_real": np.real(np.asarray(compiled_channel.distorted, dtype=np.complex128)).tolist(),
        "distorted_imag": np.imag(np.asarray(compiled_channel.distorted, dtype=np.complex128)).tolist(),
    }
