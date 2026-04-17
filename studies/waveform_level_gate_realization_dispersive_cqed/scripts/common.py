"""Shared helpers for the waveform-level gate realization study.

Physical setup: dispersive transmon + cavity system with chi, chi', K.
Units: rad/s and seconds throughout.
Tensor ordering: transmon first, cavity second (cqed_sim convention).
"""

from __future__ import annotations

import json
import os
import sys
import time
from functools import partial
from pathlib import Path
from typing import Any, Sequence

import numpy as np

import runtime_compat  # noqa: F401

from scipy.linalg import norm, polar
from scipy.spatial.transform import Rotation

import qutip as qt

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
STUDY_DIR = SCRIPT_DIR.parent
DATA_DIR = STUDY_DIR / "data"
FIGURES_DIR = STUDY_DIR / "figures"
ARTIFACTS_DIR = STUDY_DIR / "artifacts"
REPORT_DIR = STUDY_DIR / "report"
REPO_ROOT = STUDY_DIR.parent.parent
USER_ROOT = REPO_ROOT.parent
CQED_SIM_ROOT = USER_ROOT / "cQED_simulation"

for _d in (DATA_DIR, FIGURES_DIR, ARTIFACTS_DIR, REPORT_DIR):
    _d.mkdir(parents=True, exist_ok=True)

STYLE_PATH = (
    REPO_ROOT / ".github" / "skills" / "publication-figures" / "assets" / "cqed_style.mplstyle"
)


def ensure_cqed_sim_on_path() -> None:
    candidate = str(CQED_SIM_ROOT)
    if candidate not in sys.path:
        sys.path.insert(0, candidate)


ensure_cqed_sim_on_path()

# ---------------------------------------------------------------------------
# cqed_sim imports
# ---------------------------------------------------------------------------
from cqed_sim.core import DispersiveTransmonCavityModel, FrameSpec  # noqa: E402
from cqed_sim.core.conventions import qubit_cavity_block_indices  # noqa: E402
from cqed_sim.core.frequencies import (  # noqa: E402
    carrier_for_transition_frequency,
    manifold_transition_frequency,
)
from cqed_sim.core.ideal_gates import displacement_op, qubit_rotation_xy  # noqa: E402
from cqed_sim.pulses.envelopes import (  # noqa: E402
    cosine_rise_envelope,
    gaussian_area_fraction,
    gaussian_envelope,
    square_envelope,
)
from cqed_sim.pulses.pulse import Pulse  # noqa: E402
from cqed_sim.sequence import SequenceCompiler  # noqa: E402
from cqed_sim.sim import SimulationConfig, prepare_simulation  # noqa: E402

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
TWO_PI = 2.0 * np.pi

OMEGA_Q = TWO_PI * 6.150e9   # Qubit frequency (rad/s)
OMEGA_C = TWO_PI * 5.241e9   # Cavity frequency (rad/s)
ALPHA = TWO_PI * (-255.0e6)  # Transmon anharmonicity (rad/s)
CHI = TWO_PI * (-2.84e6)     # Dispersive shift (rad/s)
CHI_PRIME = TWO_PI * (-21.0e3)  # Second-order dispersive (rad/s)
KERR = TWO_PI * (-28.0e3)    # Cavity self-Kerr (rad/s)

# ---------------------------------------------------------------------------
# Simulation defaults
# ---------------------------------------------------------------------------
N_TR = 3          # Transmon levels (g, e, f) to capture leakage
N_CAV = 15        # Default cavity Fock dimension
DEFAULT_DT = 0.5e-9  # Time step (s) = 0.5 ns
PI_PULSE_DURATION = 40.0e-9   # Default qubit pi-pulse duration (s)
SIGMA_FRACTION = 0.25         # Gaussian sigma/duration

# ---------------------------------------------------------------------------
# Pauli matrices (2x2 numpy)
# ---------------------------------------------------------------------------
PAULI_X = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
PAULI_Y = np.array([[0.0, -1j], [1j, 0.0]], dtype=np.complex128)
PAULI_Z = np.array([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)
PAULIS = (PAULI_X, PAULI_Y, PAULI_Z)

# ---------------------------------------------------------------------------
# Plot style
# ---------------------------------------------------------------------------
TOL_BRIGHT = ['#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB']


def apply_plot_style() -> None:
    import matplotlib.pyplot as plt
    if STYLE_PATH.exists():
        plt.style.use(str(STYLE_PATH))


# ---------------------------------------------------------------------------
# Model and frame builders
# ---------------------------------------------------------------------------
def build_model(
    *,
    n_cav: int = N_CAV,
    n_tr: int = N_TR,
    chi: float = CHI,
    chi_prime: float | None = CHI_PRIME,
    kerr: float = KERR,
) -> DispersiveTransmonCavityModel:
    chi_higher = (chi_prime,) if chi_prime is not None else ()
    return DispersiveTransmonCavityModel(
        omega_q=OMEGA_Q,
        omega_c=OMEGA_C,
        alpha=ALPHA,
        chi=chi,
        chi_higher=chi_higher,
        kerr=kerr,
        n_cav=n_cav,
        n_tr=n_tr,
    )


def build_frame(model: DispersiveTransmonCavityModel) -> FrameSpec:
    return FrameSpec(omega_q_frame=float(model.omega_q), omega_c_frame=float(model.omega_c))


# ---------------------------------------------------------------------------
# Pulse construction
# ---------------------------------------------------------------------------
def make_square_displacement_pulse(
    model: DispersiveTransmonCavityModel,
    frame: FrameSpec,
    *,
    alpha: complex,
    duration_s: float,
    channel: str = "storage",
    t0: float = 0.0,
    label: str | None = None,
) -> Pulse:
    """Square-envelope cavity displacement pulse.

    Amplitude calibration: epsilon = i * alpha / T  (rotating-frame convention).
    The carrier is set to zero detuning from the cavity in the rotating frame.
    The complex amplitude is decomposed: amp = |epsilon|, phase = arg(epsilon).
    """
    epsilon = 1j * alpha / duration_s
    return Pulse(
        channel=channel,
        t0=t0,
        duration=duration_s,
        envelope=square_envelope,
        carrier=0.0,  # On resonance with cavity in bare-freq frame
        phase=float(np.angle(epsilon)),
        amp=float(abs(epsilon)),
        drag=0.0,
        label=label or "displacement",
    )


def envelope_area(envelope_func, *, n_pts: int = 4097) -> float:
    """Return the real-valued area fraction of a normalized envelope function."""
    grid = np.linspace(0.0, 1.0, int(n_pts))
    env = np.asarray(envelope_func(grid), dtype=np.complex128)
    trapezoid = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    return float(trapezoid(np.real(env), grid))


def displacement_envelope(
    family: str,
    *,
    sigma_fraction: float = SIGMA_FRACTION,
    rise_fraction: float = 0.15,
):
    """Return a normalized cavity-displacement envelope by family name."""
    family_key = str(family).strip().lower()
    if family_key == "square":
        return square_envelope
    if family_key == "gaussian":
        return partial(gaussian_envelope, sigma=float(sigma_fraction), center=0.5)
    if family_key in {"cosine", "cosine_rise", "cosine-rise"}:
        return partial(cosine_rise_envelope, rise_fraction=float(rise_fraction))
    raise ValueError(f"Unsupported displacement envelope family '{family}'.")


def make_shaped_displacement_pulse(
    model: DispersiveTransmonCavityModel,
    frame: FrameSpec,
    *,
    alpha: complex,
    duration_s: float,
    family: str = "square",
    sigma_fraction: float = SIGMA_FRACTION,
    rise_fraction: float = 0.15,
    carrier: float = 0.0,
    channel: str = "storage",
    t0: float = 0.0,
    label: str | None = None,
) -> Pulse:
    """General single-tone cavity displacement pulse with calibrated envelope area."""
    envelope = displacement_envelope(
        family,
        sigma_fraction=sigma_fraction,
        rise_fraction=rise_fraction,
    )
    area = envelope_area(envelope)
    epsilon = 1j * alpha / (duration_s * max(area, 1.0e-15))
    return Pulse(
        channel=channel,
        t0=t0,
        duration=duration_s,
        envelope=envelope,
        carrier=float(carrier),
        phase=float(np.angle(epsilon)),
        amp=float(abs(epsilon)),
        drag=0.0,
        label=label or f"disp_{family}",
    )


def cavity_branch_transition_frequency(
    model: DispersiveTransmonCavityModel,
    frame: FrameSpec,
    *,
    qubit_level: int,
    cavity_level: int = 0,
) -> float:
    """Cavity transition frequency in the rotating frame for a fixed qubit branch."""
    n0 = int(cavity_level)
    return float(
        model.basis_energy(int(qubit_level), n0 + 1, frame)
        - model.basis_energy(int(qubit_level), n0, frame)
    )


def make_gaussian_qubit_pulse(
    model: DispersiveTransmonCavityModel,
    frame: FrameSpec,
    *,
    theta: float,
    phase: float,
    duration_s: float,
    manifold_level: int = 0,
    sigma_fraction: float = SIGMA_FRACTION,
    drag: float = 0.0,
    channel: str = "qubit",
    t0: float = 0.0,
    label: str | None = None,
) -> Pulse:
    """Gaussian (optionally DRAG-corrected) qubit rotation pulse.

    Targets the |g,n> <-> |e,n> transition for manifold_level n.
    Amplitude calibrated so the pulse area gives rotation angle theta.
    """
    transition = manifold_transition_frequency(model, manifold_level, frame=frame)
    carrier = carrier_for_transition_frequency(transition)
    area_frac = gaussian_area_fraction(sigma_fraction)
    amplitude = theta / (2.0 * duration_s * area_frac)
    envelope = partial(gaussian_envelope, sigma=sigma_fraction, center=0.5)
    return Pulse(
        channel=channel,
        t0=t0,
        duration=duration_s,
        envelope=envelope,
        carrier=carrier,
        phase=phase,
        amp=amplitude,
        drag=drag,
        label=label or f"qubit_rot_{theta:.2f}",
    )


# ---------------------------------------------------------------------------
# Compilation and simulation
# ---------------------------------------------------------------------------
def compile_and_prepare(
    model: DispersiveTransmonCavityModel,
    frame: FrameSpec,
    pulses: list[Pulse],
    *,
    dt: float = DEFAULT_DT,
    drive_ops: dict[str, str] | None = None,
    store_states: bool = False,
) -> Any:
    """Compile pulses and prepare a SimulationSession for reuse."""
    if drive_ops is None:
        drive_ops = {}
        channels = {p.channel for p in pulses}
        if "storage" in channels:
            drive_ops["storage"] = "cavity"
        if "qubit" in channels:
            drive_ops["qubit"] = "qubit"
    total_duration = max(p.t0 + p.duration for p in pulses)
    compiler = SequenceCompiler(dt=dt)
    compiled = compiler.compile(pulses, t_end=total_duration + dt)
    config = SimulationConfig(frame=frame, store_states=store_states)
    return prepare_simulation(model, compiled, drive_ops, config=config, e_ops={})


def simulate_state(session: Any, initial_state: qt.Qobj) -> qt.Qobj:
    """Run a single simulation and return the final state."""
    return session.run(initial_state).final_state


def simulate_with_trajectory(session: Any, initial_state: qt.Qobj) -> Any:
    """Run a simulation with state storage and return the full result."""
    return session.run(initial_state)


# ---------------------------------------------------------------------------
# Operator extraction
# ---------------------------------------------------------------------------
def propagate_basis_states(
    model: DispersiveTransmonCavityModel,
    session: Any,
    *,
    cavity_levels: Sequence[int],
    qubit_levels: Sequence[int] = (0, 1),
) -> np.ndarray:
    """Propagate computational basis states to build the full propagator matrix.

    Returns the propagator U such that U[:, idx] is the output column for
    basis state |q, n>. Indices follow the cqed_sim flat convention:
    idx = q * n_cav + n.
    """
    full_dim = int(model.n_tr) * int(model.n_cav)
    U = np.zeros((full_dim, full_dim), dtype=np.complex128)
    for n in cavity_levels:
        for q in qubit_levels:
            idx = int(q) * int(model.n_cav) + int(n)
            psi0 = model.basis_state(q, n)
            psi_f = simulate_state(session, psi0)
            U[:, idx] = np.asarray(psi_f.full(), dtype=np.complex128).ravel()
    return U


def extract_qubit_block(
    U: np.ndarray,
    n_cav: int,
    cavity_level: int,
) -> np.ndarray:
    """Extract the 2x2 qubit block for a given cavity Fock level from the full propagator."""
    idx_g = int(cavity_level)
    idx_e = int(n_cav) + int(cavity_level)
    return U[np.ix_([idx_g, idx_e], [idx_g, idx_e])]


def extract_cavity_block(
    U: np.ndarray,
    n_cav: int,
    qubit_level: int,
    cavity_levels: Sequence[int],
) -> np.ndarray:
    """Extract the cavity sub-block for a given qubit level from the full propagator."""
    indices = [int(qubit_level) * int(n_cav) + int(n) for n in cavity_levels]
    return U[np.ix_(indices, indices)]


# ---------------------------------------------------------------------------
# SU(2) analysis
# ---------------------------------------------------------------------------
def nearest_su2(matrix: np.ndarray) -> np.ndarray:
    """Project a 2x2 matrix to the nearest SU(2) element."""
    unitary, _ = polar(np.asarray(matrix, dtype=np.complex128))
    det = np.linalg.det(unitary)
    phase = np.angle(det) / 2.0
    return np.asarray(unitary * np.exp(-1j * phase), dtype=np.complex128)


def nearest_unitary(matrix: np.ndarray) -> np.ndarray:
    """Project a square matrix to the nearest unitary via the polar decomposition."""
    unitary, _ = polar(np.asarray(matrix, dtype=np.complex128))
    return np.asarray(unitary, dtype=np.complex128)


def su2_to_rotation(block: np.ndarray) -> Rotation:
    """Convert a 2x2 SU(2) matrix to a scipy Rotation via the Bloch-sphere SO(3) map."""
    u = nearest_su2(block)
    rot = np.zeros((3, 3), dtype=float)
    ud = u.conj().T
    for i, si in enumerate(PAULIS):
        for j, sj in enumerate(PAULIS):
            rot[i, j] = 0.5 * float(np.real(np.trace(si @ u @ sj @ ud)))
    return Rotation.from_matrix(rot)


def rotation_metrics(block: np.ndarray, target: np.ndarray | None = None) -> dict[str, float]:
    """Extract rotation angle, axis, and error metrics from a 2x2 block.

    If target is given, also compute fidelity and error decomposition.
    """
    u = nearest_su2(block)
    rot = su2_to_rotation(u)
    rotvec = rot.as_rotvec()
    angle = float(np.linalg.norm(rotvec))
    axis = rotvec / angle if angle > 1e-9 else np.array([0.0, 0.0, 1.0])

    result = {
        "rotation_angle_rad": angle,
        "axis_x": float(axis[0]),
        "axis_y": float(axis[1]),
        "axis_z": float(axis[2]),
        "unitarity_error": float(norm(u.conj().T @ u - np.eye(2), ord=2)),
    }

    if target is not None:
        t = nearest_su2(target)
        dim = 2.0
        overlap = np.trace(t.conj().T @ u)
        proc_fid = float(np.clip(abs(overlap) ** 2 / (dim * dim), 0.0, 1.0))
        avg_fid = float((dim * proc_fid + 1.0) / (dim + 1.0))
        result["process_fidelity"] = proc_fid
        result["average_gate_fidelity"] = avg_fid
        result["frobenius_error"] = float(norm(u - t, ord="fro") / norm(t, ord="fro"))

        error_u = nearest_su2(t.conj().T @ u)
        err_rot = su2_to_rotation(error_u)
        err_vec = err_rot.as_rotvec()
        result["residual_z_rad"] = float(abs(err_vec[2]))
        result["transverse_error_rad"] = float(np.linalg.norm(err_vec[:2]))
        result["error_rotation_angle_rad"] = float(np.linalg.norm(err_vec))

        target_rot = su2_to_rotation(t)
        target_vec = target_rot.as_rotvec()
        target_angle = float(np.linalg.norm(target_vec))
        result["target_rotation_angle_rad"] = target_angle
        result["angle_error_rad"] = abs(angle - target_angle)

    return result


# ---------------------------------------------------------------------------
# Displacement analysis
# ---------------------------------------------------------------------------
def displacement_fidelity(
    final_state: qt.Qobj,
    ideal_state: qt.Qobj,
) -> float:
    """State fidelity between the final state and ideal displaced state."""
    if final_state.isket and ideal_state.isket:
        return float(abs(ideal_state.dag() * final_state) ** 2)
    rho = final_state if final_state.isoper else final_state.proj()
    sigma = ideal_state if ideal_state.isoper else ideal_state.proj()
    from qutip import fidelity
    return float(fidelity(rho, sigma) ** 2)


def cavity_state_from_joint(
    joint_state: qt.Qobj,
    n_tr: int,
    n_cav: int,
    qubit_level: int,
) -> tuple[qt.Qobj, float]:
    """Extract the conditioned cavity state for a given qubit level.

    Returns (normalized cavity state, population in that qubit sector).
    """
    rho = joint_state if joint_state.isoper else joint_state.proj()
    full = np.asarray(rho.full(), dtype=np.complex128)
    # Indices for qubit_level sector: q*n_cav + n for n in 0..n_cav-1
    indices = list(range(qubit_level * n_cav, (qubit_level + 1) * n_cav))
    block = full[np.ix_(indices, indices)]
    pop = float(np.real(np.trace(block)))
    if pop <= 0.0:
        return qt.Qobj(np.zeros((n_cav, n_cav), dtype=np.complex128), dims=[[n_cav], [n_cav]]), 0.0
    return qt.Qobj(block / pop, dims=[[n_cav], [n_cav]]), pop


def qubit_state_from_joint(
    joint_state: qt.Qobj,
    n_tr: int,
    n_cav: int,
    cavity_level: int,
) -> tuple[qt.Qobj, float]:
    """Extract the conditioned qubit state for a given cavity Fock level.

    Returns (normalized qubit density matrix, population in that cavity sector).
    """
    rho = joint_state if joint_state.isoper else joint_state.proj()
    full = np.asarray(rho.full(), dtype=np.complex128)
    indices = [q * n_cav + cavity_level for q in range(n_tr)]
    block = full[np.ix_(indices, indices)]
    pop = float(np.real(np.trace(block)))
    if pop <= 0.0:
        return qt.Qobj(np.zeros((n_tr, n_tr), dtype=np.complex128), dims=[[n_tr], [n_tr]]), 0.0
    return qt.Qobj(block / pop, dims=[[n_tr], [n_tr]]), pop


def entanglement_entropy(
    joint_state: qt.Qobj,
    n_tr: int,
    n_cav: int,
) -> float:
    """Von Neumann entropy of the reduced qubit state (entanglement with cavity)."""
    rho = joint_state if joint_state.isoper else joint_state.proj()
    rho_q = rho.ptrace(0)  # Trace out cavity (transmon first)
    evals = np.real(rho_q.eigenenergies())
    evals = evals[evals > 0]
    return float(-np.sum(evals * np.log2(evals + 1e-30)))


def qubit_purity(
    joint_state: qt.Qobj,
) -> float:
    """Purity of the reduced qubit state."""
    rho = joint_state if joint_state.isoper else joint_state.proj()
    rho_q = rho.ptrace(0)
    return float(np.real((rho_q * rho_q).tr()))


def fock_populations(
    joint_state: qt.Qobj,
    n_cav: int,
) -> np.ndarray:
    """Cavity Fock-state populations (summed over qubit levels)."""
    rho = joint_state if joint_state.isoper else joint_state.proj()
    rho_c = rho.ptrace(1)  # Trace out qubit
    return np.real(np.diag(np.asarray(rho_c.full())))


def annihilation_expectation(joint_state: qt.Qobj) -> complex:
    """Expectation value of the cavity annihilation operator."""
    rho = joint_state if joint_state.isoper else joint_state.proj()
    rho_c = rho.ptrace(1)
    n_cav = int(rho_c.shape[0])
    a = qt.destroy(n_cav)
    return complex(qt.expect(a, rho_c))


def state_trace_distance(state_a: qt.Qobj, state_b: qt.Qobj) -> float:
    """Trace distance between two pure or mixed states."""
    rho = state_a if state_a.isoper else state_a.proj()
    sigma = state_b if state_b.isoper else state_b.proj()
    return float(qt.metrics.tracedist(rho, sigma))


# ---------------------------------------------------------------------------
# Bloch vector
# ---------------------------------------------------------------------------
def bloch_vector(rho_q: qt.Qobj) -> tuple[float, float, float]:
    """Bloch vector (x, y, z) from a 2x2 qubit density matrix."""
    m = np.asarray(rho_q.full(), dtype=np.complex128)[:2, :2]
    return (
        float(np.real(np.trace(m @ PAULI_X))),
        float(np.real(np.trace(m @ PAULI_Y))),
        float(np.real(np.trace(m @ PAULI_Z))),
    )


# ---------------------------------------------------------------------------
# Serialization
# ---------------------------------------------------------------------------
def json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): json_ready(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(v) for v in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        if np.iscomplexobj(value):
            return {"real": np.real(value).tolist(), "imag": np.imag(value).tolist(), "shape": list(value.shape)}
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


def _artifact_document(
    path: Path,
    payload: Any,
    *,
    description: str | None = None,
    load_instructions: str | None = None,
) -> dict[str, Any]:
    if isinstance(payload, dict) and "artifact_payload" in payload:
        return payload
    return {
        "study_name": STUDY_DIR.name,
        "date_created": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "description": description or path.stem.replace("_", " "),
        "load_instructions": load_instructions or f"Use common.load_json(Path('{path.name}')) to load the payload.",
        "artifact_payload": payload,
    }


def save_json(
    path: Path,
    payload: Any,
    *,
    description: str | None = None,
    load_instructions: str | None = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    document = _artifact_document(
        path,
        payload,
        description=description,
        load_instructions=load_instructions,
    )
    text = json.dumps(json_ready(document), indent=2, sort_keys=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    last_error: OSError | None = None
    for _ in range(5):
        try:
            temp.write_text(text, encoding="utf-8")
            os.replace(str(temp), str(path))
            return
        except OSError as exc:
            last_error = exc
            time.sleep(0.25)
    if last_error is not None:
        raise last_error


def load_json(path: Path) -> Any:
    document = json.loads(path.read_text(encoding="utf-8-sig"))
    if isinstance(document, dict) and "artifact_payload" in document:
        return document["artifact_payload"]
    return document
