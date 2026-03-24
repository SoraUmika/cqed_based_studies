"""Shared helpers for the literature-informed selective-primitives study.

All Hamiltonian frequencies are in rad/s and all times are in seconds.
The study uses `cqed_sim` for pulse compilation and simulation; only the
missing pulse-family definitions are implemented locally.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Callable

import numpy as np

from runtime_compat import patch_windows_qutip_import

patch_windows_qutip_import()

import qutip as qt


SCRIPT_DIR = Path(__file__).resolve().parent
STUDY_DIR = SCRIPT_DIR.parent
DATA_DIR = STUDY_DIR / "data"
FIGURES_DIR = STUDY_DIR / "figures"
REPORT_DIR = STUDY_DIR / "report"

from cqed_sim.core import DispersiveTransmonCavityModel, FrameSpec
from cqed_sim.core.frequencies import carrier_for_transition_frequency, manifold_transition_frequency
from cqed_sim.pulses import Pulse, normalized_gaussian
from cqed_sim.sequence import SequenceCompiler
from cqed_sim.sim import SimulationConfig, prepare_simulation
from cqed_sim.sim.noise import NoiseSpec


DATA_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
REPORT_DIR.mkdir(parents=True, exist_ok=True)

TWO_PI = 2.0 * np.pi

# Typical dispersive cQED parameters from AGENTS.md.
OMEGA_Q = TWO_PI * 6.150e9
OMEGA_C = TWO_PI * 5.241e9
ALPHA = TWO_PI * (-255.0e6)
CHI = TWO_PI * (-2.84e6)
CHI_PRIME = TWO_PI * (-21.0e3)
KERR = TWO_PI * (-28.0e3)

# Representative simulation settings.
LOGICAL_N = 4
N_CAV_SIM = LOGICAL_N + 3
N_TR_SIM = 3
DT = 2.0e-9
TARGET_BRANCH = 1
SQR_THETA = np.pi
SQR_PHI = 0.0
SNAP_PHASE = 0.5 * np.pi

# Representative noise point.
TRANSMON_T1 = (50.0e-6, 25.0e-6)
QUBIT_TPHI = 80.0e-6
KAPPA_STORAGE = TWO_PI * 1.0e3
NTH_STORAGE = 0.02

COLOR_MAP = {
    "gaussian": "#1f77b4",
    "cosine_squared": "#d95f02",
    "flat_top_gaussian": "#1b9e77",
}

SNAP_COLOR_MAP = {
    "gaussian": "#4c78a8",
    "flat_top_gaussian": "#f58518",
}


EnvelopeFunc = Callable[[np.ndarray], np.ndarray]


def duration_from_chi_t(chi_t_2pi: float, *, chi: float = CHI) -> float:
    """Return T such that |chi| T / (2pi) = chi_t_2pi."""
    return float(chi_t_2pi / (abs(float(chi)) / TWO_PI))


def build_model(*, n_cav: int = N_CAV_SIM, n_tr: int = N_TR_SIM) -> DispersiveTransmonCavityModel:
    """Construct the representative dispersive qubit-storage model."""
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
    """Use the qubit/cavity rotating frame matching the bare frequencies."""
    return FrameSpec(
        omega_c_frame=float(model.omega_c),
        omega_q_frame=float(model.omega_q),
    )


def build_noise_spec() -> NoiseSpec:
    """Representative noise point for noisy replay."""
    return NoiseSpec(
        transmon_t1=tuple(float(value) for value in TRANSMON_T1),
        tphi=float(QUBIT_TPHI),
        kappa=float(KAPPA_STORAGE),
        nth=float(NTH_STORAGE),
    )


def json_dump(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def _normalized_real_envelope(values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    area = float(np.trapezoid(np.real(values), grid))
    if abs(area) < 1.0e-14:
        return np.asarray(values, dtype=np.complex128)
    return np.asarray(values / area, dtype=np.complex128)


def gaussian_envelope_factory(sigma_fraction: float) -> EnvelopeFunc:
    sigma = float(sigma_fraction)

    def envelope(t_rel: np.ndarray) -> np.ndarray:
        return np.asarray(normalized_gaussian(np.asarray(t_rel, dtype=float), sigma_fraction=sigma), dtype=np.complex128)

    return envelope


def cosine_squared_envelope() -> EnvelopeFunc:
    grid = np.linspace(0.0, 1.0, 4001)
    base = 2.0 * np.cos(np.pi * (grid - 0.5)) ** 2
    normalization = float(np.trapezoid(base, grid))

    def envelope(t_rel: np.ndarray) -> np.ndarray:
        t = np.asarray(t_rel, dtype=float)
        values = 2.0 * np.cos(np.pi * (t - 0.5)) ** 2
        return np.asarray(values / normalization, dtype=np.complex128)

    return envelope


def flat_top_gaussian_envelope_factory(ramp_fraction: float) -> EnvelopeFunc:
    ramp = float(ramp_fraction)
    ramp = max(1.0e-3, min(ramp, 0.45))
    sigma = ramp / 3.0
    grid = np.linspace(0.0, 1.0, 4001)

    def _raw(t: np.ndarray) -> np.ndarray:
        x = np.asarray(t, dtype=float)
        out = np.ones_like(x, dtype=float)
        left = x < ramp
        right = x > (1.0 - ramp)
        out[left] = np.exp(-0.5 * ((x[left] - ramp) / sigma) ** 2)
        out[right] = np.exp(-0.5 * ((x[right] - (1.0 - ramp)) / sigma) ** 2)
        return out

    normalization = float(np.trapezoid(_raw(grid), grid))

    def envelope(t_rel: np.ndarray) -> np.ndarray:
        return np.asarray(_raw(np.asarray(t_rel, dtype=float)) / normalization, dtype=np.complex128)

    return envelope


def family_envelope(family: str, shape_parameter: float) -> EnvelopeFunc:
    if family == "gaussian":
        return gaussian_envelope_factory(shape_parameter)
    if family == "cosine_squared":
        return cosine_squared_envelope()
    if family == "flat_top_gaussian":
        return flat_top_gaussian_envelope_factory(shape_parameter)
    raise ValueError(f"Unsupported family '{family}'.")


def build_selective_qubit_pulse(
    model: DispersiveTransmonCavityModel,
    frame: FrameSpec,
    *,
    family: str,
    branch: int,
    theta: float,
    phi: float,
    duration_s: float,
    shape_parameter: float,
    amplitude_scale: float = 1.0,
    t0: float = 0.0,
    label: str | None = None,
) -> Pulse:
    """Create a branch-selective qubit pulse centered on the chosen manifold."""
    omega_branch = manifold_transition_frequency(model, int(branch), frame)
    carrier = carrier_for_transition_frequency(omega_branch)
    amp = float(amplitude_scale) * float(theta) / (2.0 * float(duration_s))
    return Pulse(
        channel="qubit",
        t0=float(t0),
        duration=float(duration_s),
        envelope=family_envelope(family, shape_parameter),
        carrier=float(carrier),
        phase=float(phi),
        amp=float(amp),
        label=label,
    )


def build_sqr_pulses(
    model: DispersiveTransmonCavityModel,
    frame: FrameSpec,
    *,
    family: str,
    branch: int,
    theta: float,
    phi: float,
    duration_s: float,
    shape_parameter: float,
    amplitude_scale: float,
) -> tuple[list[Pulse], dict[str, str], float]:
    pulse = build_selective_qubit_pulse(
        model,
        frame,
        family=family,
        branch=branch,
        theta=theta,
        phi=phi,
        duration_s=duration_s,
        shape_parameter=shape_parameter,
        amplitude_scale=amplitude_scale,
        label=f"sqr_{family}",
    )
    return [pulse], {"qubit": "qubit"}, float(duration_s)


def build_snap_pulses(
    model: DispersiveTransmonCavityModel,
    frame: FrameSpec,
    *,
    family: str,
    branch: int,
    phase_angle: float,
    duration_s: float,
    shape_parameter: float,
    amplitude_scale: float,
) -> tuple[list[Pulse], dict[str, str], float]:
    """Create the geometric SNAP sequence from two number-selective pi pulses.

    Following the literature convention, the second pulse is phase shifted by
    `pi - theta`, where `theta` is the desired cavity phase on the selected
    Fock manifold.
    """
    first = build_selective_qubit_pulse(
        model,
        frame,
        family=family,
        branch=branch,
        theta=np.pi,
        phi=0.0,
        duration_s=duration_s,
        shape_parameter=shape_parameter,
        amplitude_scale=amplitude_scale,
        t0=0.0,
        label=f"snap_{family}_a",
    )
    second = build_selective_qubit_pulse(
        model,
        frame,
        family=family,
        branch=branch,
        theta=np.pi,
        phi=float(np.pi - phase_angle),
        duration_s=duration_s,
        shape_parameter=shape_parameter,
        amplitude_scale=amplitude_scale,
        t0=float(duration_s),
        label=f"snap_{family}_b",
    )
    total_duration = float(2.0 * duration_s)
    return [first, second], {"qubit": "qubit"}, total_duration


def build_session(
    model: DispersiveTransmonCavityModel,
    frame: FrameSpec,
    pulses: list[Pulse],
    drive_ops: dict[str, str],
    *,
    total_duration_s: float,
    noise: NoiseSpec | None = None,
    dt: float = DT,
) -> object:
    compiler = SequenceCompiler(dt=float(dt))
    compiled = compiler.compile(pulses, t_end=float(total_duration_s + 4.0 * dt))
    return prepare_simulation(
        model,
        compiled,
        drive_ops,
        config=SimulationConfig(frame=frame, store_states=False),
        noise=noise,
        e_ops={},
    )


def logical_indices(model: DispersiveTransmonCavityModel, logical_n: int = LOGICAL_N) -> list[int]:
    indices: list[int] = []
    for n in range(int(logical_n)):
        indices.extend([n, int(model.n_cav) + n])
    return indices


def cavity_ground_indices(model: DispersiveTransmonCavityModel, logical_n: int = LOGICAL_N) -> list[int]:
    return [n for n in range(int(logical_n))]


def qobj_probability_in_indices(state: qt.Qobj, indices: list[int]) -> float:
    idx = np.asarray(indices, dtype=int)
    if state.isket:
        vec = np.asarray(state.full(), dtype=np.complex128).reshape(-1)
        return float(np.sum(np.abs(vec[idx]) ** 2))
    arr = np.asarray(state.full(), dtype=np.complex128)
    return float(np.real(np.trace(arr[np.ix_(idx, idx)])))


def state_fidelity_to_pure_target(target_state: qt.Qobj, actual_state: qt.Qobj) -> float:
    target_vec = np.asarray(target_state.full(), dtype=np.complex128).reshape(-1)
    if actual_state.isket:
        actual_vec = np.asarray(actual_state.full(), dtype=np.complex128).reshape(-1)
        return float(np.clip(abs(np.vdot(target_vec, actual_vec)) ** 2, 0.0, 1.0))
    rho = np.asarray(actual_state.full(), dtype=np.complex128)
    return float(np.clip(np.real(target_vec.conj() @ rho @ target_vec), 0.0, 1.0))


def target_qubit_unitary(theta: float, phi: float) -> np.ndarray:
    c = np.cos(theta / 2.0)
    s = np.sin(theta / 2.0)
    return np.array(
        [
            [c, -1j * s * np.exp(-1j * phi)],
            [-1j * s * np.exp(1j * phi), c],
        ],
        dtype=np.complex128,
    )


def sqr_strict_target_operator(logical_n: int, target_branch: int, theta: float, phi: float) -> np.ndarray:
    operator = np.zeros((2 * int(logical_n), 2 * int(logical_n)), dtype=np.complex128)
    target_block = target_qubit_unitary(theta, phi)
    identity = np.eye(2, dtype=np.complex128)
    for n in range(int(logical_n)):
        block = target_block if n == int(target_branch) else identity
        operator[2 * n : 2 * n + 2, 2 * n : 2 * n + 2] = block
    return operator


def z_corrected_target_fidelity(u_actual: np.ndarray, u_target: np.ndarray) -> tuple[float, float]:
    u_conj = u_target.conj()
    a_term = u_conj[0, 0] * u_actual[0, 0] + u_conj[0, 1] * u_actual[0, 1]
    b_term = u_conj[1, 0] * u_actual[1, 0] + u_conj[1, 1] * u_actual[1, 1]
    alpha_opt = float(np.angle(a_term) - np.angle(b_term)) if abs(b_term) > 1.0e-15 else 0.0
    fid = float(np.clip((abs(a_term) + abs(b_term)) ** 2 / 4.0, 0.0, 1.0))
    return fid, alpha_opt


def spectator_z_fidelity(u_actual: np.ndarray) -> tuple[float, float]:
    u00 = complex(u_actual[0, 0])
    u11 = complex(u_actual[1, 1])
    phase_opt = float(np.angle(u11) - np.angle(u00))
    fid = float(np.clip((abs(u00) + abs(u11)) ** 2 / 4.0, 0.0, 1.0))
    return fid, phase_opt


def sqr_relaxed_target_operator_from_actual(
    restricted_operator: np.ndarray,
    *,
    logical_n: int,
    target_branch: int,
    theta: float,
    phi: float,
) -> tuple[np.ndarray, dict[str, object]]:
    target = target_qubit_unitary(theta, phi)
    relaxed = np.zeros_like(restricted_operator, dtype=np.complex128)
    spectator_phases: list[float] = []
    target_alpha = 0.0
    branch_fidelities: list[float] = []
    for n in range(int(logical_n)):
        block = restricted_operator[2 * n : 2 * n + 2, 2 * n : 2 * n + 2]
        if n == int(target_branch):
            fid, target_alpha = z_corrected_target_fidelity(block, target)
            correction = np.diag([1.0, np.exp(-1j * target_alpha)])
            relaxed_block = correction @ target
            spectator_phases.append(None)
        else:
            fid, phase_opt = spectator_z_fidelity(block)
            relaxed_block = np.diag([np.exp(-0.5j * phase_opt), np.exp(0.5j * phase_opt)])
            spectator_phases.append(phase_opt)
        branch_fidelities.append(fid)
        relaxed[2 * n : 2 * n + 2, 2 * n : 2 * n + 2] = relaxed_block
    return relaxed, {
        "target_z_correction": float(target_alpha),
        "spectator_phases": spectator_phases,
        "branch_cphase_fidelities": branch_fidelities,
        "branch_cphase_mean": float(np.mean(np.asarray(branch_fidelities, dtype=float))),
    }


def extract_restricted_operator(
    final_states: list[qt.Qobj],
    model: DispersiveTransmonCavityModel,
    *,
    logical_n: int = LOGICAL_N,
) -> np.ndarray:
    indices = logical_indices(model, logical_n)
    operator = np.zeros((len(indices), len(final_states)), dtype=np.complex128)
    for column, final_state in enumerate(final_states):
        if not final_state.isket:
            raise TypeError("Restricted operator extraction requires pure-state outputs.")
        vec = np.asarray(final_state.full(), dtype=np.complex128).reshape(-1)
        operator[:, column] = vec[np.asarray(indices, dtype=int)]
    return operator


def sqr_probe_state_vectors(logical_n: int = LOGICAL_N) -> list[dict[str, object]]:
    dim = 2 * int(logical_n)
    eye = np.eye(dim, dtype=np.complex128)
    probes: list[dict[str, object]] = []
    for n in range(int(logical_n)):
        g_idx = 2 * n
        e_idx = 2 * n + 1
        probes.append({"label": f"basis_g_{n}", "vector": eye[:, g_idx]})
        probes.append({"label": f"basis_e_{n}", "vector": eye[:, e_idx]})
        probes.append({"label": f"x_super_{n}", "vector": (eye[:, g_idx] + eye[:, e_idx]) / np.sqrt(2.0)})
        probes.append({"label": f"y_super_{n}", "vector": (eye[:, g_idx] + 1.0j * eye[:, e_idx]) / np.sqrt(2.0)})
    return probes


def snap_target_operator(logical_n: int, target_branch: int, phase_angle: float) -> np.ndarray:
    diag = np.ones(int(logical_n), dtype=np.complex128)
    diag[int(target_branch)] = np.exp(1j * float(phase_angle))
    return np.diag(diag)


def snap_probe_state_vectors(logical_n: int = LOGICAL_N) -> list[dict[str, object]]:
    eye = np.eye(int(logical_n), dtype=np.complex128)
    probes: list[dict[str, object]] = []
    for n in range(int(logical_n)):
        probes.append({"label": f"basis_{n}", "vector": eye[:, n]})
    for left in range(int(logical_n)):
        for right in range(left + 1, int(logical_n)):
            probes.append({"label": f"x_pair_{left}_{right}", "vector": (eye[:, left] + eye[:, right]) / np.sqrt(2.0)})
            probes.append({"label": f"y_pair_{left}_{right}", "vector": (eye[:, left] + 1.0j * eye[:, right]) / np.sqrt(2.0)})
    return probes


def embed_model_logical_state(vector: np.ndarray, model: DispersiveTransmonCavityModel, indices: list[int]) -> qt.Qobj:
    full = np.zeros(int(model.n_tr) * int(model.n_cav), dtype=np.complex128)
    full[np.asarray(indices, dtype=int)] = np.asarray(vector, dtype=np.complex128)
    return qt.Qobj(full, dims=[[int(model.n_tr), int(model.n_cav)], [1, 1]])


def average_target_state_fidelity(
    session: object,
    probe_vectors: list[dict[str, object]],
    target_operator: np.ndarray,
    *,
    model: DispersiveTransmonCavityModel,
    indices: list[int],
) -> tuple[float, list[dict[str, float]]]:
    actual_fidelities: list[dict[str, float]] = []
    states = [embed_model_logical_state(np.asarray(item["vector"], dtype=np.complex128), model, indices) for item in probe_vectors]
    results = session.run_many(states)
    for item, result in zip(probe_vectors, results, strict=True):
        final_state = result.final_state
        target_vector = np.asarray(target_operator @ np.asarray(item["vector"], dtype=np.complex128), dtype=np.complex128)
        target_state = embed_model_logical_state(target_vector, model, indices)
        actual_fidelities.append(
            {
                "label": str(item["label"]),
                "fidelity": state_fidelity_to_pure_target(target_state, final_state),
                "in_logical_population": qobj_probability_in_indices(final_state, indices),
            }
        )
    mean_fidelity = float(np.mean([item["fidelity"] for item in actual_fidelities]))
    return mean_fidelity, actual_fidelities


def sample_total_waveform(pulses: list[Pulse], *, dt: float = 1.0e-9) -> tuple[np.ndarray, np.ndarray]:
    t_end = max(float(pulse.t1) for pulse in pulses)
    count = int(np.ceil(t_end / float(dt))) + 1
    tlist = np.arange(count, dtype=float) * float(dt)
    total = np.zeros_like(tlist, dtype=np.complex128)
    for pulse in pulses:
        total += np.asarray(pulse.sample(tlist), dtype=np.complex128)
    return tlist, total
