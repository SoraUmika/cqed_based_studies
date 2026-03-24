"""Shared helpers for the open-system SQR deep-dive study.

This module keeps the study runnable from its own folder while reusing the
physics conventions established in the completed SQR study. All frequencies are
in rad/s and all durations are in seconds.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

from runtime_compat import patch_windows_qutip_import

patch_windows_qutip_import()

SCRIPT_DIR = Path(__file__).resolve().parent
STUDY_DIR = SCRIPT_DIR.parent
DATA_DIR = STUDY_DIR / "data"
REFERENCE_DATA_DIR = STUDY_DIR.parent / "sqr_pulse_waveform_design" / "data"
CQED_SIM_PATH = Path(
    r"C:\Users\dazzl\Box\Shyam Shankar Quantum Circuits Group\Users\Users_JianJun\cQED_simulation"
)

if str(CQED_SIM_PATH) not in sys.path:
    sys.path.insert(0, str(CQED_SIM_PATH))

import qutip as qt

from cqed_sim.core import (
    DispersiveReadoutTransmonStorageModel,
    DispersiveTransmonCavityModel,
    FrameSpec,
)
from cqed_sim.core.frequencies import carrier_for_transition_frequency, manifold_transition_frequency
from cqed_sim.measurement import AmplifierChain, PurcellFilter, ReadoutChain, ReadoutResonator
from cqed_sim.pulses import Pulse
from cqed_sim.pulses.envelopes import MultitoneTone, multitone_gaussian_envelope, normalized_gaussian
from cqed_sim.sequence import SequenceCompiler
from cqed_sim.sim import SimulationConfig, prepare_simulation
from cqed_sim.sim.noise import NoiseSpec

TWO_PI = 2.0 * np.pi

# Device parameters from AGENTS.md / study READMEs.
OMEGA_Q = TWO_PI * 6.150e9
OMEGA_S = TWO_PI * 5.241e9
OMEGA_R = TWO_PI * 8.597e9
ALPHA = TWO_PI * (-255.0e6)
CHI_S = TWO_PI * (-2.84e6)
CHI_PRIME = TWO_PI * (-21.0e3)
KERR_S = TWO_PI * (-28.0e3)

# Three-mode defaults. chi_r defaults to the same dispersive scale as chi_s; chi_sr
# is left at zero unless the caller explicitly requests a cross-Kerr scenario.
CHI_R_DEFAULT = CHI_S
CHI_SR_DEFAULT = 0.0
KERR_R_DEFAULT = 0.0

# Readout-chain defaults used for Purcell and backaction estimates.
KAPPA_R = TWO_PI * 2.4e6
G_READOUT = TWO_PI * 100.0e6

# Legacy SQR noise point from the completed study.
T1_LEGACY = 20.0e-6
T2_LEGACY = 20.0e-6

# Representative realistic-noise point for replay studies.
T1_GE_DEFAULT = 30.0e-6
T1_FE_DEFAULT = 10.0e-6
TPHI_STORAGE_DEFAULT = 80.0e-6
TPHI_READOUT_DEFAULT = 20.0e-6
KAPPA_STORAGE_DEFAULT = TWO_PI * 10.0e3

N_STORAGE_LOGICAL = 4
N_STORAGE_SIM = 8
N_READOUT_SIM = 6
N_CAV_TWO_MODE = N_STORAGE_LOGICAL + 2
N_TR = 3
DT = 2.0e-9

CHI_T_VALUES = np.array([0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0], dtype=float)
FAMILY_NAMES = (
    "single_tone_gaussian",
    "square",
    "cosine_squared",
    "multitone_one_segment",
)
TARGET_STORAGE_LEVEL = 1
THETA_TARGET = np.pi
PHI_TARGET = 0.0

DATA_DIR.mkdir(parents=True, exist_ok=True)


def pure_dephasing_time(t1: float | None, t2: float | None) -> float | None:
    """Return T_phi from T1 and T2, or None if not well-defined."""
    if t1 is None or t2 is None:
        return None
    denom = (1.0 / float(t2)) - (0.5 / float(t1))
    if denom <= 0.0:
        return None
    return float(1.0 / denom)


def combine_t1_limits(*t1_values: float | None) -> float:
    """Combine independent T1 channels into one effective T1."""
    inv_total = 0.0
    for value in t1_values:
        if value is None or not np.isfinite(value) or value <= 0.0:
            continue
        inv_total += 1.0 / float(value)
    return float(np.inf if inv_total <= 0.0 else 1.0 / inv_total)


def duration_from_chi_t(chi_t_2pi: float, chi_s: float = CHI_S) -> float:
    """Return T such that |chi/(2pi)| * T = chi_t_2pi."""
    return float(chi_t_2pi / (abs(float(chi_s)) / TWO_PI))


def build_two_mode_model(
    *,
    chi_s: float = CHI_S,
    chi_prime: float = CHI_PRIME,
    kerr_s: float = KERR_S,
    n_cav: int = N_CAV_TWO_MODE,
    n_tr: int = N_TR,
) -> DispersiveTransmonCavityModel:
    """Build the two-mode dispersive model used in the completed SQR study."""
    chi_higher = (float(chi_prime),) if abs(float(chi_prime)) > 0.0 else ()
    return DispersiveTransmonCavityModel(
        omega_c=OMEGA_S,
        omega_q=OMEGA_Q,
        alpha=ALPHA,
        chi=float(chi_s),
        chi_higher=chi_higher,
        kerr=float(kerr_s),
        n_cav=int(n_cav),
        n_tr=int(n_tr),
    )


def build_three_mode_model(
    *,
    chi_s: float = CHI_S,
    chi_r: float = CHI_R_DEFAULT,
    chi_sr: float = CHI_SR_DEFAULT,
    kerr_s: float = KERR_S,
    kerr_r: float = KERR_R_DEFAULT,
    n_storage: int = N_STORAGE_SIM,
    n_readout: int = N_READOUT_SIM,
    n_tr: int = N_TR,
) -> DispersiveReadoutTransmonStorageModel:
    """Build the three-mode qubit + storage + readout dispersive model."""
    return DispersiveReadoutTransmonStorageModel(
        omega_s=OMEGA_S,
        omega_r=OMEGA_R,
        omega_q=OMEGA_Q,
        alpha=ALPHA,
        chi_s=float(chi_s),
        chi_r=float(chi_r),
        chi_sr=float(chi_sr),
        kerr_s=float(kerr_s),
        kerr_r=float(kerr_r),
        n_storage=int(n_storage),
        n_readout=int(n_readout),
        n_tr=int(n_tr),
    )


def build_frame(model) -> FrameSpec:
    """Construct the matching rotating frame for a two- or three-mode model."""
    if hasattr(model, "omega_s") and hasattr(model, "omega_r"):
        return FrameSpec(
            omega_c_frame=float(model.omega_s),
            omega_q_frame=float(model.omega_q),
            omega_r_frame=float(model.omega_r),
        )
    return FrameSpec(
        omega_c_frame=float(model.omega_c),
        omega_q_frame=float(model.omega_q),
    )


def build_legacy_noise_spec(
    *,
    t1: float = T1_LEGACY,
    t2: float = T2_LEGACY,
    kappa_storage: float | None = None,
    nth_storage: float = 0.0,
) -> NoiseSpec:
    """Legacy single-T1 / single-T2 noise model from the completed SQR study."""
    tphi = pure_dephasing_time(t1, t2)
    kwargs: dict[str, float] = {
        "t1": float(t1),
    }
    if tphi is not None:
        kwargs["tphi"] = float(tphi)
    if kappa_storage is not None:
        kwargs["kappa"] = float(kappa_storage)
        kwargs["nth"] = float(nth_storage)
    return NoiseSpec(**kwargs)


def build_multilevel_noise_spec(
    *,
    transmon_t1: tuple[float | None, ...] = (T1_GE_DEFAULT, T1_FE_DEFAULT),
    t2_qubit: float = T2_LEGACY,
    kappa_storage: float | None = KAPPA_STORAGE_DEFAULT,
    nth_storage: float = 0.0,
    tphi_storage: float | None = TPHI_STORAGE_DEFAULT,
    kappa_readout: float | None = None,
    nth_readout: float = 0.0,
    tphi_readout: float | None = TPHI_READOUT_DEFAULT,
    for_three_mode: bool = False,
) -> NoiseSpec:
    """Explicit multilevel transmon-T1 noise model for A2 sweeps and replay."""
    t1_ge = transmon_t1[0] if transmon_t1 else None
    tphi_qubit = pure_dephasing_time(t1_ge, t2_qubit)
    kwargs: dict[str, object] = {
        "transmon_t1": tuple(transmon_t1),
    }
    if tphi_qubit is not None:
        kwargs["tphi"] = float(tphi_qubit)
    if for_three_mode:
        if kappa_storage is not None:
            kwargs["kappa_storage"] = float(kappa_storage)
            kwargs["nth_storage"] = float(nth_storage)
        if tphi_storage is not None:
            kwargs["tphi_storage"] = float(tphi_storage)
        if kappa_readout is not None:
            kwargs["kappa_readout"] = float(kappa_readout)
            kwargs["nth_readout"] = float(nth_readout)
        if tphi_readout is not None:
            kwargs["tphi_readout"] = float(tphi_readout)
    else:
        if kappa_storage is not None:
            kwargs["kappa"] = float(kappa_storage)
            kwargs["nth"] = float(nth_storage)
    return NoiseSpec(**kwargs)


def build_nominal_readout_chain(
    *,
    epsilon: complex | float,
    chi_r: float = CHI_R_DEFAULT,
    omega_r: float = OMEGA_R,
    omega_q: float = OMEGA_Q,
    kappa_r: float = KAPPA_R,
    g_coupling: float = G_READOUT,
    include_filter: bool = True,
    filter_bandwidth: float | None = None,
    integration_time: float = 400.0e-9,
    dt: float = DT,
    noise_temperature: float = 0.0,
    gain: float = 1.0,
) -> ReadoutChain:
    """Build the nominal readout chain used for Purcell and backaction analysis."""
    resonator = ReadoutResonator(
        omega_r=float(omega_r),
        kappa=float(kappa_r),
        g=float(g_coupling),
        epsilon=epsilon,
        chi=float(chi_r),
        drive_frequency=float(omega_r),
    )
    purcell_filter = None
    if include_filter:
        purcell_filter = PurcellFilter(
            bandwidth=float(kappa_r if filter_bandwidth is None else filter_bandwidth)
        )
    return ReadoutChain(
        resonator=resonator,
        amplifier=AmplifierChain(noise_temperature=float(noise_temperature), gain=float(gain)),
        purcell_filter=purcell_filter,
        integration_time=float(integration_time),
        dt=float(dt),
    )


def basis_state(model, qubit_level: int, storage_level: int, readout_level: int = 0) -> qt.Qobj:
    """Return a basis state for either the two-mode or three-mode model."""
    try:
        return model.basis_state(int(qubit_level), int(storage_level), int(readout_level))
    except TypeError:
        return model.basis_state(int(qubit_level), int(storage_level))


def qubit_branch_transition_frequency(
    model,
    frame: FrameSpec,
    storage_level: int,
    *,
    readout_level: int = 0,
) -> float:
    """Return the qubit transition frequency for the selected storage branch."""
    if hasattr(model, "transmon_transition_frequency") and hasattr(model, "omega_s"):
        return float(
            model.transmon_transition_frequency(
                storage_level=int(storage_level),
                readout_level=int(readout_level),
                lower_level=0,
                upper_level=1,
                frame=frame,
            )
        )
    return float(manifold_transition_frequency(model, int(storage_level), frame))


def readout_midpoint_frequency(
    model,
    frame: FrameSpec,
    *,
    storage_level: int = 0,
    readout_level: int = 0,
) -> float:
    """Return the midpoint readout drive frequency for the three-mode model."""
    if not hasattr(model, "readout_transition_frequency"):
        raise TypeError("Readout midpoint frequency is only defined for the three-mode model.")
    freq_g = float(model.readout_transition_frequency(0, int(storage_level), int(readout_level), frame))
    freq_e = float(model.readout_transition_frequency(1, int(storage_level), int(readout_level), frame))
    return 0.5 * (freq_g + freq_e)


def target_qubit_unitary(theta: float = THETA_TARGET, phi: float = PHI_TARGET) -> np.ndarray:
    """Return the target qubit rotation matrix R(theta, phi)."""
    cos_term = np.cos(theta / 2.0)
    sin_term = np.sin(theta / 2.0)
    return np.array(
        [
            [cos_term, -1j * sin_term * np.exp(-1j * phi)],
            [-1j * sin_term * np.exp(1j * phi), cos_term],
        ],
        dtype=np.complex128,
    )


def state_fidelity_pure_target(target_state: qt.Qobj, state: qt.Qobj) -> float:
    """Return fidelity against a pure target ket for ket or density-matrix inputs."""
    if state.isoper:
        target_vec = target_state.full().reshape(-1)
        rho = state.full()
        return float(np.real(target_vec.conj() @ rho @ target_vec))
    return float(abs(target_state.overlap(state)) ** 2)


def reduce_qubit_storage(state: qt.Qobj) -> qt.Qobj:
    """Trace out the readout mode when present."""
    rho = state if state.isoper else state.proj()
    if len(rho.dims[0]) < 3:
        return rho
    return rho.ptrace([0, 1])


def reduce_storage(state: qt.Qobj) -> qt.Qobj:
    """Return the reduced storage-mode density matrix."""
    rho = state if state.isoper else state.proj()
    if len(rho.dims[0]) < 2:
        return rho
    return rho.ptrace([1])


def storage_coherence(state: qt.Qobj, lower: int = 0, upper: int = 1) -> complex:
    """Return the selected storage-mode off-diagonal coherence element."""
    rho_storage = reduce_storage(state)
    arr = rho_storage.full()
    return complex(arr[int(lower), int(upper)])


def build_storage_superposition_state(
    model,
    *,
    lower: int = 0,
    upper: int = 1,
    qubit_level: int = 0,
    readout_level: int = 0,
) -> qt.Qobj:
    """Return (|lower> + |upper>) / sqrt(2) in the storage mode."""
    psi_lower = basis_state(model, qubit_level, lower, readout_level)
    psi_upper = basis_state(model, qubit_level, upper, readout_level)
    return (psi_lower + psi_upper).unit()


def build_target_state(
    model,
    *,
    storage_level: int,
    qubit_level: int,
    target_storage_level: int = TARGET_STORAGE_LEVEL,
    readout_level: int = 0,
) -> qt.Qobj:
    """Return the ideal basis-state target for an SQR pulse on one storage branch."""
    target_qubit_level = 1 - int(qubit_level) if int(storage_level) == int(target_storage_level) else int(qubit_level)
    return basis_state(model, target_qubit_level, int(storage_level), int(readout_level))


def build_basis_initial_states(
    model,
    *,
    n_storage_levels: int = N_STORAGE_LOGICAL,
    readout_level: int = 0,
) -> tuple[list[tuple[int, int]], list[qt.Qobj]]:
    """Return labels and computational-basis initial states for branch replay."""
    labels: list[tuple[int, int]] = []
    states: list[qt.Qobj] = []
    for storage_level in range(int(n_storage_levels)):
        for qubit_level in (0, 1):
            labels.append((qubit_level, storage_level))
            states.append(basis_state(model, qubit_level, storage_level, readout_level))
    return labels, states


def _gaussian_envelope(sigma_fraction: float = 1.0 / 6.0):
    def env(t_rel):
        return normalized_gaussian(t_rel, sigma_fraction=sigma_fraction)
    return env


def _square_envelope(t_rel):
    return np.ones_like(np.asarray(t_rel, dtype=float))


def _cosine_squared_envelope(t_rel):
    t = np.asarray(t_rel, dtype=float)
    return 2.0 * np.cos(np.pi * (t - 0.5)) ** 2


def build_single_tone_gaussian(
    model,
    frame: FrameSpec,
    *,
    storage_level: int = TARGET_STORAGE_LEVEL,
    theta: float = THETA_TARGET,
    phi: float = PHI_TARGET,
    duration: float,
    readout_level: int = 0,
    channel: str = "q",
) -> tuple[list[Pulse], dict[str, str]]:
    """Build a single-tone Gaussian selective qubit pulse."""
    omega_branch = qubit_branch_transition_frequency(
        model,
        frame,
        int(storage_level),
        readout_level=int(readout_level),
    )
    amp = float(theta) / (2.0 * float(duration))
    pulse = Pulse(
        channel=channel,
        t0=0.0,
        duration=float(duration),
        envelope=_gaussian_envelope(),
        carrier=carrier_for_transition_frequency(omega_branch),
        phase=float(phi),
        amp=amp,
    )
    return [pulse], {channel: "qubit"}


def build_square_pulse(
    model,
    frame: FrameSpec,
    *,
    storage_level: int = TARGET_STORAGE_LEVEL,
    theta: float = THETA_TARGET,
    phi: float = PHI_TARGET,
    duration: float,
    readout_level: int = 0,
    channel: str = "q",
) -> tuple[list[Pulse], dict[str, str]]:
    """Build a square selective qubit pulse."""
    omega_branch = qubit_branch_transition_frequency(
        model,
        frame,
        int(storage_level),
        readout_level=int(readout_level),
    )
    amp = float(theta) / (2.0 * float(duration))
    pulse = Pulse(
        channel=channel,
        t0=0.0,
        duration=float(duration),
        envelope=_square_envelope,
        carrier=carrier_for_transition_frequency(omega_branch),
        phase=float(phi),
        amp=amp,
    )
    return [pulse], {channel: "qubit"}


def build_cosine_squared_pulse(
    model,
    frame: FrameSpec,
    *,
    storage_level: int = TARGET_STORAGE_LEVEL,
    theta: float = THETA_TARGET,
    phi: float = PHI_TARGET,
    duration: float,
    readout_level: int = 0,
    channel: str = "q",
) -> tuple[list[Pulse], dict[str, str]]:
    """Build a cosine-squared (Hann-like) selective qubit pulse."""
    omega_branch = qubit_branch_transition_frequency(
        model,
        frame,
        int(storage_level),
        readout_level=int(readout_level),
    )
    amp = float(theta) / (2.0 * float(duration))
    pulse = Pulse(
        channel=channel,
        t0=0.0,
        duration=float(duration),
        envelope=_cosine_squared_envelope,
        carrier=carrier_for_transition_frequency(omega_branch),
        phase=float(phi),
        amp=amp,
    )
    return [pulse], {channel: "qubit"}


def build_multitone_pulse(
    model,
    frame: FrameSpec,
    *,
    n_storage_levels: int = N_STORAGE_LOGICAL,
    target_storage_level: int = TARGET_STORAGE_LEVEL,
    theta: float = THETA_TARGET,
    phi: float = PHI_TARGET,
    duration: float,
    readout_level: int = 0,
    channel: str = "q",
) -> tuple[list[Pulse], dict[str, str]]:
    """Build a one-segment multitone Gaussian SQR pulse."""
    tone_specs: list[MultitoneTone] = []
    for storage_level in range(int(n_storage_levels)):
        omega_branch = qubit_branch_transition_frequency(
            model,
            frame,
            storage_level,
            readout_level=int(readout_level),
        )
        amp = float(theta) / (2.0 * float(duration)) if storage_level == int(target_storage_level) else 0.0
        tone_specs.append(
            MultitoneTone(
                manifold=storage_level,
                omega_rad_s=carrier_for_transition_frequency(omega_branch),
                amp_rad_s=amp,
                phase_rad=float(phi) if storage_level == int(target_storage_level) else 0.0,
            )
        )

    def env(t_rel):
        return multitone_gaussian_envelope(
            t_rel,
            duration_s=float(duration),
            sigma_fraction=1.0 / 6.0,
            tone_specs=tone_specs,
        )

    pulse = Pulse(
        channel=channel,
        t0=0.0,
        duration=float(duration),
        envelope=env,
        carrier=0.0,
        phase=0.0,
        amp=1.0,
    )
    return [pulse], {channel: "qubit"}


def build_family_pulse(
    model,
    frame: FrameSpec,
    family_name: str,
    *,
    duration: float,
    target_storage_level: int = TARGET_STORAGE_LEVEL,
    theta: float = THETA_TARGET,
    phi: float = PHI_TARGET,
    n_storage_levels: int = N_STORAGE_LOGICAL,
    readout_level: int = 0,
) -> tuple[list[Pulse], dict[str, str]]:
    """Dispatch helper for the supported parametric SQR pulse families."""
    if family_name == "single_tone_gaussian":
        return build_single_tone_gaussian(
            model,
            frame,
            storage_level=target_storage_level,
            theta=theta,
            phi=phi,
            duration=duration,
            readout_level=readout_level,
        )
    if family_name == "square":
        return build_square_pulse(
            model,
            frame,
            storage_level=target_storage_level,
            theta=theta,
            phi=phi,
            duration=duration,
            readout_level=readout_level,
        )
    if family_name == "cosine_squared":
        return build_cosine_squared_pulse(
            model,
            frame,
            storage_level=target_storage_level,
            theta=theta,
            phi=phi,
            duration=duration,
            readout_level=readout_level,
        )
    if family_name == "multitone_one_segment":
        return build_multitone_pulse(
            model,
            frame,
            n_storage_levels=n_storage_levels,
            target_storage_level=target_storage_level,
            theta=theta,
            phi=phi,
            duration=duration,
            readout_level=readout_level,
        )
    raise ValueError(f"Unsupported family '{family_name}'.")


def build_readout_square_pulse(
    model,
    frame: FrameSpec,
    *,
    duration: float,
    amplitude: float,
    storage_level: int = 0,
    readout_level: int = 0,
    phase: float = 0.0,
    channel: str = "readout",
) -> tuple[list[Pulse], dict[str, str]]:
    """Build a square readout pulse on the readout mode midpoint frequency."""
    omega_mid = readout_midpoint_frequency(
        model,
        frame,
        storage_level=int(storage_level),
        readout_level=int(readout_level),
    )
    pulse = Pulse(
        channel=channel,
        t0=0.0,
        duration=float(duration),
        envelope=_square_envelope,
        carrier=carrier_for_transition_frequency(omega_mid),
        phase=float(phase),
        amp=float(amplitude),
    )
    return [pulse], {channel: "readout"}


def build_session(
    model,
    frame: FrameSpec,
    pulses: list[Pulse],
    drive_ops: dict[str, str],
    *,
    duration: float,
    noise: NoiseSpec | None = None,
    store_states: bool = False,
    dt: float = DT,
):
    """Compile pulses and build a reusable simulation session."""
    compiler = SequenceCompiler(dt=float(dt))
    compiled = compiler.compile(pulses, t_end=float(duration) + 4.0 * float(dt))
    config = SimulationConfig(frame=frame, store_states=store_states)
    return prepare_simulation(
        model,
        compiled,
        drive_ops,
        config=config,
        noise=noise,
        e_ops={},
    )


def run_session_over_states(session, initial_states: list[qt.Qobj]) -> list[qt.Qobj]:
    """Run a prepared simulation session over a list of initial states."""
    final_states: list[qt.Qobj] = []
    for psi0 in initial_states:
        final_states.append(session.run(psi0).final_state)
    return final_states


def branch_average(values: np.ndarray, n_storage_levels: int = N_STORAGE_LOGICAL) -> np.ndarray:
    """Average qubit-g and qubit-e values within each storage branch."""
    averages = np.zeros(int(n_storage_levels), dtype=float)
    for storage_level in range(int(n_storage_levels)):
        start = 2 * storage_level
        averages[storage_level] = float(np.mean(values[start:start + 2]))
    return averages


def load_reference_phase5_results() -> dict[str, np.ndarray | list[str]]:
    """Load the completed SQR phase-5 archive for baseline comparisons."""
    payload = np.load(REFERENCE_DATA_DIR / "phase5_results.npz", allow_pickle=True)
    data = {key: payload[key] for key in payload.files}
    data["family_names"] = [str(name) for name in payload["family_names"]]
    return data


def load_reference_grape_benchmark() -> dict[str, np.ndarray]:
    """Load the completed SQR GRAPE benchmark archive."""
    payload = np.load(REFERENCE_DATA_DIR / "grape_benchmark_results.npz", allow_pickle=True)
    return {key: payload[key] for key in payload.files}