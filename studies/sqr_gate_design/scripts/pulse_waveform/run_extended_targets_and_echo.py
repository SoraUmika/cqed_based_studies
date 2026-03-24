"""
Extended closed-system SQR study with generalized targets and echoed sequences.

This extension keeps the same effective dispersive model as the base study, but
widens the target family and evaluates the resulting control as a full logical
qubit-cavity operator on a truncated Hilbert space. The numerical families are:

- single-tone Gaussian selective pulse
- one-segment common-envelope multitone baseline
- echoed single-tone Gaussian sequence
- echoed one-segment multitone sequence

The echoed families implement
    half-SQR -> pi_x -> half-SQR' -> pi_x
with the second half phase-conjugated (phi -> -phi), matching the toggling-frame
condition for an x-axis refocusing pulse.

Usage:
    python scripts/run_extended_targets_and_echo.py

Output:
    data/extended_targets_results.npz
"""

from __future__ import annotations

from dataclasses import dataclass
import sys
import time
from pathlib import Path
from typing import Callable

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import runtime_compat  # noqa: F401

from common import (
    CHI,
    CHI_PRIME,
    DT,
    KERR,
    N_TR,
    build_frame,
    build_model,
    duration_from_chi_t,
    extract_branch_unitaries,
    extract_leakage,
    identity_fidelity_with_z,
    spectator_transverse_error,
    spectator_z_fidelity,
    target_qubit_unitary,
    z_corrected_target_fidelity,
)
from cqed_sim.core.frequencies import carrier_for_transition_frequency, manifold_transition_frequency
from cqed_sim.sequence import SequenceCompiler
from cqed_sim.sim import SimulationConfig, prepare_simulation


EnvelopeFunc = Callable[[np.ndarray], np.ndarray]


@dataclass(frozen=True)
class Pulse:
    channel: str
    t0: float
    duration: float
    envelope: EnvelopeFunc | np.ndarray
    carrier: float = 0.0
    phase: float = 0.0
    amp: float = 1.0
    drag: float = 0.0
    sample_rate: float | None = None
    label: str | None = None

    @property
    def t1(self) -> float:
        return self.t0 + self.duration

    def _sample_analytic(self, t: np.ndarray) -> np.ndarray:
        t_rel = (t - self.t0) / self.duration
        in_support = (t_rel >= 0.0) & (t_rel < 1.0)
        out = np.zeros_like(t, dtype=np.complex128)
        if np.any(in_support):
            env = np.asarray(self.envelope(t_rel[in_support]), dtype=np.complex128)
            phase = np.exp(1j * (self.carrier * t[in_support] + self.phase))
            out[in_support] = self.amp * env * phase
        return out

    def _sample_discrete(self, t: np.ndarray) -> np.ndarray:
        arr = np.asarray(self.envelope, dtype=np.complex128)
        if self.sample_rate is None:
            raise ValueError("sample_rate is required when envelope is sampled.")
        in_support = (t >= self.t0) & (t < self.t1)
        out = np.zeros_like(t, dtype=np.complex128)
        idx = np.floor((t[in_support] - self.t0) * self.sample_rate).astype(int)
        idx = np.clip(idx, 0, arr.size - 1)
        phase = np.exp(1j * (self.carrier * t[in_support] + self.phase))
        out[in_support] = self.amp * arr[idx] * phase
        return out

    def sample(self, t: np.ndarray) -> np.ndarray:
        if callable(self.envelope):
            return self._sample_analytic(t)
        return self._sample_discrete(t)


@dataclass(frozen=True)
class MultitoneTone:
    manifold: int
    omega_rad_s: float
    amp_rad_s: float
    phase_rad: float


def gaussian_envelope(t_rel: np.ndarray, sigma: float, center: float | None = None) -> np.ndarray:
    center = 0.5 if center is None else center
    return np.exp(-0.5 * ((t_rel - center) / sigma) ** 2).astype(np.complex128)


def gaussian_area_fraction(sigma_fraction: float, n_pts: int = 4097) -> float:
    grid = np.linspace(0.0, 1.0, n_pts)
    env = np.asarray(gaussian_envelope(grid, sigma=sigma_fraction), dtype=np.complex128)
    trapezoid = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    return float(trapezoid(np.real(env), grid))


def normalized_gaussian(t_rel: np.ndarray, sigma_fraction: float) -> np.ndarray:
    base = np.asarray(gaussian_envelope(t_rel, sigma=sigma_fraction), dtype=np.complex128)
    area = gaussian_area_fraction(sigma_fraction)
    return base if abs(area) < 1.0e-12 else base / area


def multitone_gaussian_envelope(t_rel: np.ndarray, duration_s: float, sigma_fraction: float, tone_specs: list[MultitoneTone]) -> np.ndarray:
    env = normalized_gaussian(t_rel, sigma_fraction=sigma_fraction)
    t = t_rel * duration_s
    coeff = np.zeros_like(t, dtype=np.complex128)
    for spec in tone_specs:
        coeff += spec.amp_rad_s * np.exp(1j * spec.phase_rad) * np.exp(1j * spec.omega_rad_s * t)
    return env * coeff

CHI_T_EXT_VALUES = np.array([1.0, 2.0, 3.0, 5.0], dtype=float)
AXIS_CASES = (
    ("X", 0.0),
    ("Diag", 0.25 * np.pi),
    ("Y", 0.5 * np.pi),
)
ANGLE_CASES = (
    ("pi_over_2", 0.5 * np.pi),
    ("pi", np.pi),
)
TARGET_BRANCH_CASES = (0, 1, 2)
LOGICAL_N_CASES = (3, 4)
FAMILY_NAMES = (
    "single_tone_gaussian",
    "multitone_one_segment",
    "echoed_single_tone_gaussian",
    "echoed_multitone_one_segment",
)

GENERAL_TARGET_BRANCH = 1
GENERAL_LOGICAL_N = 4
REPRESENTATIVE_TARGET_BRANCH = 1
REPRESENTATIVE_LOGICAL_N = 4
REPRESENTATIVE_THETA = np.pi
REPRESENTATIVE_PHI = 0.0
SIGMA_FRACTION = 1.0 / 6.0
ECHO_PI_DURATION_S = 8.0e-9


def gaussian_envelope_callable(sigma_fraction: float = SIGMA_FRACTION):
    def envelope(t_rel):
        return normalized_gaussian(t_rel, sigma_fraction=sigma_fraction)

    return envelope


def build_single_tone_gaussian_local(model, frame, target_branch: int, theta: float, phi: float, duration_s: float):
    omega_target = manifold_transition_frequency(model, int(target_branch), frame)
    carrier = carrier_for_transition_frequency(omega_target)
    amplitude = float(theta) / (2.0 * float(duration_s))
    pulse = Pulse(
        channel="q",
        t0=0.0,
        duration=float(duration_s),
        envelope=gaussian_envelope_callable(),
        carrier=float(carrier),
        phase=float(phi),
        amp=amplitude,
        label="single_tone_gaussian",
    )
    return [pulse], {"q": "qubit"}


def build_multitone_baseline_local(model, frame, logical_n: int, target_branch: int, theta: float, phi: float, duration_s: float):
    tone_specs = []
    for n in range(int(logical_n)):
        omega_n = manifold_transition_frequency(model, n, frame)
        tone_specs.append(
            MultitoneTone(
                manifold=int(n),
                omega_rad_s=float(carrier_for_transition_frequency(omega_n)),
                amp_rad_s=float(theta / (2.0 * duration_s)) if n == int(target_branch) else 0.0,
                phase_rad=float(phi) if n == int(target_branch) else 0.0,
            )
        )

    def envelope(t_rel):
        return multitone_gaussian_envelope(
            np.asarray(t_rel, dtype=float),
            duration_s=float(duration_s),
            sigma_fraction=SIGMA_FRACTION,
            tone_specs=tone_specs,
        )

    pulse = Pulse(
        channel="q",
        t0=0.0,
        duration=float(duration_s),
        envelope=envelope,
        carrier=0.0,
        phase=0.0,
        amp=1.0,
        label="multitone_one_segment",
    )
    return [pulse], {"q": "qubit"}


def shifted_pulse(base_pulse: Pulse, *, t0: float, phase: float | None = None, label: str | None = None) -> Pulse:
    kwargs = {
        "channel": base_pulse.channel,
        "t0": float(t0),
        "duration": float(base_pulse.duration),
        "envelope": base_pulse.envelope,
        "carrier": float(base_pulse.carrier),
        "phase": float(base_pulse.phase if phase is None else phase),
        "amp": float(base_pulse.amp),
        "label": label if label is not None else getattr(base_pulse, "label", None),
    }
    sample_rate = getattr(base_pulse, "sample_rate", None)
    if sample_rate is not None:
        kwargs["sample_rate"] = float(sample_rate)
    return Pulse(**kwargs)


def build_hard_global_pi_pulse(t0: float, *, phase: float = 0.0) -> Pulse:
    def square_envelope(t_rel):
        return np.ones_like(np.asarray(t_rel, dtype=float))

    amplitude = np.pi / (2.0 * ECHO_PI_DURATION_S)
    return Pulse(
        channel="q",
        t0=float(t0),
        duration=float(ECHO_PI_DURATION_S),
        envelope=square_envelope,
        carrier=0.0,
        phase=float(phase),
        amp=float(amplitude),
        label="hard_global_pi",
    )


def build_family_sequence(family_name: str, model, frame, logical_n: int, target_branch: int, theta: float, phi: float, selective_duration_s: float):
    if family_name == "single_tone_gaussian":
        pulses, drive_ops = build_single_tone_gaussian_local(model, frame, target_branch, theta, phi, selective_duration_s)
        return pulses, drive_ops, float(selective_duration_s)

    if family_name == "multitone_one_segment":
        pulses, drive_ops = build_multitone_baseline_local(model, frame, logical_n, target_branch, theta, phi, selective_duration_s)
        return pulses, drive_ops, float(selective_duration_s)

    if family_name == "echoed_single_tone_gaussian":
        half_duration_s = 0.5 * float(selective_duration_s)
        first_half, drive_ops = build_single_tone_gaussian_local(model, frame, target_branch, 0.5 * theta, phi, half_duration_s)
        second_half, _ = build_single_tone_gaussian_local(model, frame, target_branch, 0.5 * theta, -phi, half_duration_s)
        pulses = [
            shifted_pulse(first_half[0], t0=0.0, label="echo_half_a"),
            build_hard_global_pi_pulse(half_duration_s, phase=0.0),
            shifted_pulse(second_half[0], t0=half_duration_s + ECHO_PI_DURATION_S, label="echo_half_b"),
            build_hard_global_pi_pulse(selective_duration_s + ECHO_PI_DURATION_S, phase=0.0),
        ]
        return pulses, drive_ops, float(selective_duration_s + 2.0 * ECHO_PI_DURATION_S)

    if family_name == "echoed_multitone_one_segment":
        half_duration_s = 0.5 * float(selective_duration_s)
        first_half, drive_ops = build_multitone_baseline_local(model, frame, logical_n, target_branch, 0.5 * theta, phi, half_duration_s)
        second_half, _ = build_multitone_baseline_local(model, frame, logical_n, target_branch, 0.5 * theta, -phi, half_duration_s)
        pulses = [
            shifted_pulse(first_half[0], t0=0.0, label="echo_multi_half_a"),
            build_hard_global_pi_pulse(half_duration_s, phase=0.0),
            shifted_pulse(second_half[0], t0=half_duration_s + ECHO_PI_DURATION_S, label="echo_multi_half_b"),
            build_hard_global_pi_pulse(selective_duration_s + ECHO_PI_DURATION_S, phase=0.0),
        ]
        return pulses, drive_ops, float(selective_duration_s + 2.0 * ECHO_PI_DURATION_S)

    raise ValueError(f"Unknown family {family_name!r}")


def logical_indices(model, logical_n: int) -> list[int]:
    indices: list[int] = []
    for n in range(int(logical_n)):
        indices.extend([n, int(model.n_cav) + n])
    return indices


def restricted_operator_from_full(full_operator: np.ndarray, model, logical_n: int) -> np.ndarray:
    indices = logical_indices(model, logical_n)
    return np.asarray(full_operator[np.ix_(indices, indices)], dtype=np.complex128)


def target_restricted_operator(logical_n: int, target_branch: int, theta: float, phi: float) -> np.ndarray:
    dim = 2 * int(logical_n)
    operator = np.zeros((dim, dim), dtype=np.complex128)
    target_block = target_qubit_unitary(theta, phi)
    identity_block = np.eye(2, dtype=np.complex128)
    for block_index in range(int(logical_n)):
        operator[2 * block_index : 2 * block_index + 2, 2 * block_index : 2 * block_index + 2] = (
            target_block if block_index == int(target_branch) else identity_block
        )
    return operator


def branch_metrics_from_blocks(blocks: list[np.ndarray], target_branch: int, theta: float, phi: float) -> dict[str, np.ndarray | float]:
    logical_n = len(blocks)
    target_gate = target_qubit_unitary(theta, phi)
    target_fid, alpha_opt = z_corrected_target_fidelity(blocks[int(target_branch)], target_gate)

    branch_true = np.zeros(logical_n, dtype=float)
    branch_cphase = np.zeros(logical_n, dtype=float)
    branch_z_phase = np.zeros(logical_n, dtype=float)
    branch_transverse = np.zeros(logical_n, dtype=float)
    branch_global_phase = np.zeros(logical_n, dtype=float)

    for block_index, block in enumerate(blocks):
        determinant = np.linalg.det(block)
        branch_global_phase[block_index] = 0.0 if abs(determinant) < 1.0e-15 else 0.5 * float(np.angle(determinant))
        if block_index == int(target_branch):
            branch_true[block_index] = float(target_fid)
            branch_cphase[block_index] = float(target_fid)
            branch_z_phase[block_index] = float(alpha_opt)
        else:
            branch_true[block_index] = float(identity_fidelity_with_z(block, alpha_opt))
            spectator_fid, spectator_phase = spectator_z_fidelity(block)
            branch_cphase[block_index] = float(spectator_fid)
            branch_z_phase[block_index] = float(spectator_phase)
            branch_transverse[block_index] = float(spectator_transverse_error(block))

    spectator_mask = np.arange(logical_n) != int(target_branch)
    spectator_phase_spread = 0.0 if not np.any(spectator_mask) else float(np.ptp(branch_z_phase[spectator_mask]))
    spectator_max_transverse = 0.0 if not np.any(spectator_mask) else float(np.max(branch_transverse[spectator_mask]))

    return {
        "branch_true": branch_true,
        "branch_cphase": branch_cphase,
        "branch_z_phase": branch_z_phase,
        "branch_transverse": branch_transverse,
        "branch_global_phase": branch_global_phase,
        "branch_true_mean": float(np.mean(branch_true)),
        "branch_cphase_mean": float(np.mean(branch_cphase)),
        "spectator_phase_spread": spectator_phase_spread,
        "spectator_max_transverse": spectator_max_transverse,
        "block_global_phase_spread": float(np.ptp(branch_global_phase)),
    }


def logical_process_metrics(restricted_operator: np.ndarray, logical_n: int, target_branch: int, theta: float, phi: float) -> dict[str, np.ndarray | float]:
    target_operator = target_restricted_operator(logical_n, target_branch, theta, phi)
    dim = float(target_operator.shape[0])

    overlaps = []
    identity_block = np.eye(2, dtype=np.complex128)
    target_block = target_qubit_unitary(theta, phi)
    best_fit_block_phases = np.zeros(int(logical_n), dtype=float)
    for block_index in range(int(logical_n)):
        actual_block = restricted_operator[2 * block_index : 2 * block_index + 2, 2 * block_index : 2 * block_index + 2]
        ideal_block = target_block if block_index == int(target_branch) else identity_block
        overlap = np.trace(ideal_block.conj().T @ actual_block)
        overlaps.append(overlap)
        best_fit_block_phases[block_index] = 0.0 if abs(overlap) < 1.0e-15 else -float(np.angle(overlap))

    best_fit_block_phases = best_fit_block_phases - best_fit_block_phases[0]
    strict_fidelity = abs(np.trace(target_operator.conj().T @ restricted_operator)) ** 2 / (dim * dim)
    best_block_phase_fidelity = (np.sum(np.abs(np.asarray(overlaps, dtype=np.complex128))) ** 2) / (dim * dim)

    return {
        "joint_strict_fidelity": float(np.clip(strict_fidelity, 0.0, 1.0)),
        "joint_best_block_phase_fidelity": float(np.clip(best_block_phase_fidelity, 0.0, 1.0)),
        "best_fit_block_phases": best_fit_block_phases,
    }


def same_block_population_metrics(full_operator: np.ndarray, model, logical_n: int) -> dict[str, float]:
    same_block_values = []
    for n in range(int(logical_n)):
        for qubit_level in (0, 1):
            logical_index = qubit_level * int(model.n_cav) + n
            column = np.asarray(full_operator[:, logical_index], dtype=np.complex128).reshape(-1)
            same_block = float(abs(column[n]) ** 2 + abs(column[int(model.n_cav) + n]) ** 2)
            same_block_values.append(same_block)
    values = np.asarray(same_block_values, dtype=float)
    return {
        "same_block_population_mean": float(np.mean(values)),
        "same_block_population_min": float(np.min(values)),
    }


def embed_logical_state(logical_state: np.ndarray, model, logical_n: int) -> np.ndarray:
    embedded = np.zeros(int(model.n_tr) * int(model.n_cav), dtype=np.complex128)
    embedded[np.asarray(logical_indices(model, logical_n), dtype=int)] = np.asarray(logical_state, dtype=np.complex128)
    return embedded


def build_transfer_state_set(dim: int) -> list[np.ndarray]:
    eye = np.eye(dim, dtype=np.complex128)
    states = [eye[:, index] for index in range(dim)]
    for left in range(dim):
        for right in range(left + 1, dim):
            states.append((eye[:, left] + eye[:, right]) / np.sqrt(2.0))
            states.append((eye[:, left] + 1.0j * eye[:, right]) / np.sqrt(2.0))
    return states


def state_transfer_metrics(full_operator: np.ndarray, model, logical_n: int, target_branch: int, theta: float, phi: float) -> dict[str, float]:
    target_operator = target_restricted_operator(logical_n, target_branch, theta, phi)
    transfer_inputs = build_transfer_state_set(target_operator.shape[0])
    fidelities = []
    for logical_input in transfer_inputs:
        embedded_input = embed_logical_state(logical_input, model, logical_n)
        actual_output = np.asarray(full_operator @ embedded_input, dtype=np.complex128).reshape(-1)
        target_output = embed_logical_state(target_operator @ logical_input, model, logical_n)
        fidelities.append(float(np.clip(abs(np.vdot(target_output, actual_output)) ** 2, 0.0, 1.0)))
    fidelity_values = np.asarray(fidelities, dtype=float)
    return {
        "state_transfer_mean": float(np.mean(fidelity_values)),
        "state_transfer_min": float(np.min(fidelity_values)),
    }


def simulate_logical_inputs(model, frame, pulses, drive_ops, logical_n: int, total_duration_s: float):
    compiler = SequenceCompiler(dt=DT)
    compiled = compiler.compile(pulses, t_end=float(total_duration_s + 4.0 * DT))
    config = SimulationConfig(frame=frame, store_states=False)
    session = prepare_simulation(model, compiled, drive_ops, config=config, e_ops={})

    full_dim = int(model.n_tr) * int(model.n_cav)
    full_operator = np.eye(full_dim, dtype=np.complex128)
    final_states = []
    for n in range(int(logical_n)):
        for qubit_level in (0, 1):
            initial_state = model.basis_state(qubit_level, n)
            result = session.run(initial_state)
            final_state = result.final_state
            final_states.append(final_state)
            full_operator[:, qubit_level * int(model.n_cav) + n] = np.asarray(final_state.full(), dtype=np.complex128).reshape(-1)
    return full_operator, final_states


def initialize_results():
    results: dict[str, np.ndarray | list[str] | np.ndarray] = {
        "chi_t_values": CHI_T_EXT_VALUES,
        "family_names": np.array(FAMILY_NAMES, dtype=object),
        "axis_labels": np.array([label for label, _ in AXIS_CASES], dtype=object),
        "phi_values": np.array([value for _, value in AXIS_CASES], dtype=float),
        "angle_labels": np.array([label for label, _ in ANGLE_CASES], dtype=object),
        "theta_values": np.array([value for _, value in ANGLE_CASES], dtype=float),
        "target_branch_cases": np.array(TARGET_BRANCH_CASES, dtype=int),
        "logical_n_cases": np.array(LOGICAL_N_CASES, dtype=int),
        "general_target_branch": np.array(GENERAL_TARGET_BRANCH, dtype=int),
        "general_logical_n": np.array(GENERAL_LOGICAL_N, dtype=int),
        "representative_target_branch": np.array(REPRESENTATIVE_TARGET_BRANCH, dtype=int),
        "representative_logical_n": np.array(REPRESENTATIVE_LOGICAL_N, dtype=int),
        "representative_theta": np.array(REPRESENTATIVE_THETA, dtype=float),
        "representative_phi": np.array(REPRESENTATIVE_PHI, dtype=float),
        "echo_pi_duration_s": np.array(ECHO_PI_DURATION_S, dtype=float),
        "chi_rad_s": np.array(CHI, dtype=float),
        "chi_prime_rad_s": np.array(CHI_PRIME, dtype=float),
        "kerr_rad_s": np.array(KERR, dtype=float),
    }

    axis_shape = (len(FAMILY_NAMES), len(ANGLE_CASES), len(AXIS_CASES), len(CHI_T_EXT_VALUES))
    branch_shape = (len(FAMILY_NAMES), len(TARGET_BRANCH_CASES), len(CHI_T_EXT_VALUES))
    trunc_shape = (len(FAMILY_NAMES), len(LOGICAL_N_CASES), len(CHI_T_EXT_VALUES))
    representative_shape = (len(FAMILY_NAMES), len(CHI_T_EXT_VALUES))

    metric_names = (
        "branch_true_mean",
        "branch_cphase_mean",
        "joint_strict_fidelity",
        "joint_best_block_phase_fidelity",
        "state_transfer_mean",
        "state_transfer_min",
        "same_block_population_mean",
        "same_block_population_min",
        "leakage_mean",
        "leakage_max",
        "spectator_phase_spread",
        "spectator_max_transverse",
        "block_global_phase_spread",
    )
    for prefix, shape in (("axis_scan", axis_shape), ("branch_scan", branch_shape), ("trunc_scan", trunc_shape), ("representative_scan", representative_shape)):
        for metric_name in metric_names:
            results[f"{prefix}_{metric_name}"] = np.zeros(shape, dtype=float)

    results["representative_best_fit_block_phases"] = np.zeros((len(FAMILY_NAMES), len(CHI_T_EXT_VALUES), REPRESENTATIVE_LOGICAL_N), dtype=float)
    results["representative_branch_true_blocks"] = np.zeros((len(FAMILY_NAMES), len(CHI_T_EXT_VALUES), REPRESENTATIVE_LOGICAL_N), dtype=float)
    results["representative_branch_cphase_blocks"] = np.zeros((len(FAMILY_NAMES), len(CHI_T_EXT_VALUES), REPRESENTATIVE_LOGICAL_N), dtype=float)
    results["representative_branch_z_phases"] = np.zeros((len(FAMILY_NAMES), len(CHI_T_EXT_VALUES), REPRESENTATIVE_LOGICAL_N), dtype=float)

    return results


def collect_metrics(full_operator: np.ndarray, final_states, model, logical_n: int, target_branch: int, theta: float, phi: float):
    blocks = extract_branch_unitaries(final_states, model, int(logical_n))
    restricted = restricted_operator_from_full(full_operator, model, int(logical_n))
    branch_metrics = branch_metrics_from_blocks(blocks, int(target_branch), float(theta), float(phi))
    process_metrics = logical_process_metrics(restricted, int(logical_n), int(target_branch), float(theta), float(phi))
    transfer_metrics = state_transfer_metrics(full_operator, model, int(logical_n), int(target_branch), float(theta), float(phi))
    same_block_metrics = same_block_population_metrics(full_operator, model, int(logical_n))
    leakage = extract_leakage(final_states, model, int(logical_n))

    return {
        **branch_metrics,
        **process_metrics,
        **transfer_metrics,
        **same_block_metrics,
        "leakage_mean": float(np.mean(leakage)),
        "leakage_max": float(np.max(leakage)),
    }


def store_metric_slice(results: dict[str, np.ndarray | list[str] | np.ndarray], prefix: str, family_index: int, slice_index: tuple[int, ...], metrics: dict[str, np.ndarray | float]):
    for metric_name in (
        "branch_true_mean",
        "branch_cphase_mean",
        "joint_strict_fidelity",
        "joint_best_block_phase_fidelity",
        "state_transfer_mean",
        "state_transfer_min",
        "same_block_population_mean",
        "same_block_population_min",
        "leakage_mean",
        "leakage_max",
        "spectator_phase_spread",
        "spectator_max_transverse",
        "block_global_phase_spread",
    ):
        results[f"{prefix}_{metric_name}"][(family_index,) + slice_index] = float(metrics[metric_name])


def run_axis_scan(results: dict[str, np.ndarray | list[str] | np.ndarray]):
    print("Axis/angle scan")
    for family_index, family_name in enumerate(FAMILY_NAMES):
        print(f"  Family: {family_name}")
        for angle_index, (_angle_label, theta) in enumerate(ANGLE_CASES):
            for axis_index, (_axis_label, phi) in enumerate(AXIS_CASES):
                for chi_index, chi_t_value in enumerate(CHI_T_EXT_VALUES):
                    model = build_model(chi_prime=CHI_PRIME, kerr=KERR, n_cav=GENERAL_LOGICAL_N + 2, n_tr=N_TR)
                    frame = build_frame(model)
                    selective_duration_s = duration_from_chi_t(float(chi_t_value))
                    pulses, drive_ops, total_duration_s = build_family_sequence(
                        family_name,
                        model,
                        frame,
                        GENERAL_LOGICAL_N,
                        GENERAL_TARGET_BRANCH,
                        theta,
                        phi,
                        selective_duration_s,
                    )
                    full_operator, final_states = simulate_logical_inputs(model, frame, pulses, drive_ops, GENERAL_LOGICAL_N, total_duration_s)
                    metrics = collect_metrics(full_operator, final_states, model, GENERAL_LOGICAL_N, GENERAL_TARGET_BRANCH, theta, phi)
                    store_metric_slice(results, "axis_scan", family_index, (angle_index, axis_index, chi_index), metrics)
                    if theta == REPRESENTATIVE_THETA and abs(phi - REPRESENTATIVE_PHI) < 1.0e-12:
                        results["representative_best_fit_block_phases"][family_index, chi_index] = metrics["best_fit_block_phases"]
                        results["representative_branch_true_blocks"][family_index, chi_index] = metrics["branch_true"]
                        results["representative_branch_cphase_blocks"][family_index, chi_index] = metrics["branch_cphase"]
                        results["representative_branch_z_phases"][family_index, chi_index] = metrics["branch_z_phase"]
                        store_metric_slice(results, "representative_scan", family_index, (chi_index,), metrics)
                    print(
                        f"    theta={theta/np.pi:.2f}pi phi={phi/np.pi:.2f}pi chiT={chi_t_value:.1f} "
                        f"F_true={metrics['branch_true_mean']:.4f} "
                        f"F_cphase={metrics['branch_cphase_mean']:.4f} "
                        f"F_joint={metrics['joint_strict_fidelity']:.4f}"
                    )


def run_branch_scan(results: dict[str, np.ndarray | list[str] | np.ndarray]):
    print("Target-branch scan")
    for family_index, family_name in enumerate(FAMILY_NAMES):
        print(f"  Family: {family_name}")
        for branch_index, target_branch in enumerate(TARGET_BRANCH_CASES):
            for chi_index, chi_t_value in enumerate(CHI_T_EXT_VALUES):
                model = build_model(chi_prime=CHI_PRIME, kerr=KERR, n_cav=GENERAL_LOGICAL_N + 2, n_tr=N_TR)
                frame = build_frame(model)
                selective_duration_s = duration_from_chi_t(float(chi_t_value))
                pulses, drive_ops, total_duration_s = build_family_sequence(
                    family_name,
                    model,
                    frame,
                    GENERAL_LOGICAL_N,
                    int(target_branch),
                    np.pi,
                    0.0,
                    selective_duration_s,
                )
                full_operator, final_states = simulate_logical_inputs(model, frame, pulses, drive_ops, GENERAL_LOGICAL_N, total_duration_s)
                metrics = collect_metrics(full_operator, final_states, model, GENERAL_LOGICAL_N, int(target_branch), np.pi, 0.0)
                store_metric_slice(results, "branch_scan", family_index, (branch_index, chi_index), metrics)
                print(
                    f"    n0={target_branch} chiT={chi_t_value:.1f} "
                    f"F_true={metrics['branch_true_mean']:.4f} "
                    f"F_joint={metrics['joint_strict_fidelity']:.4f}"
                )


def run_truncation_scan(results: dict[str, np.ndarray | list[str] | np.ndarray]):
    print("Logical-subspace scan")
    for family_index, family_name in enumerate(FAMILY_NAMES):
        print(f"  Family: {family_name}")
        for logical_index, logical_n in enumerate(LOGICAL_N_CASES):
            target_branch = min(1, int(logical_n) - 1)
            for chi_index, chi_t_value in enumerate(CHI_T_EXT_VALUES):
                model = build_model(chi_prime=CHI_PRIME, kerr=KERR, n_cav=int(logical_n) + 2, n_tr=N_TR)
                frame = build_frame(model)
                selective_duration_s = duration_from_chi_t(float(chi_t_value))
                pulses, drive_ops, total_duration_s = build_family_sequence(
                    family_name,
                    model,
                    frame,
                    int(logical_n),
                    int(target_branch),
                    np.pi,
                    0.0,
                    selective_duration_s,
                )
                full_operator, final_states = simulate_logical_inputs(model, frame, pulses, drive_ops, int(logical_n), total_duration_s)
                metrics = collect_metrics(full_operator, final_states, model, int(logical_n), int(target_branch), np.pi, 0.0)
                store_metric_slice(results, "trunc_scan", family_index, (logical_index, chi_index), metrics)
                print(
                    f"    logical_n={logical_n} chiT={chi_t_value:.1f} "
                    f"F_true={metrics['branch_true_mean']:.4f} "
                    f"F_joint={metrics['joint_strict_fidelity']:.4f}"
                )


def main():
    print("Extended SQR study: generalized targets and echoed constructions")
    print(f"Model: chi'=2pi*{CHI_PRIME/(2*np.pi)/1e3:.1f} kHz, K=2pi*{KERR/(2*np.pi)/1e3:.1f} kHz, n_tr={N_TR}")
    print(f"Echo hard-pi duration: {ECHO_PI_DURATION_S*1e9:.1f} ns")
    print()

    results = initialize_results()
    start_time = time.time()

    run_axis_scan(results)
    print()
    run_branch_scan(results)
    print()
    run_truncation_scan(results)

    elapsed_s = time.time() - start_time
    print(f"\nExtended study complete in {elapsed_s:.1f} s")

    data_path = SCRIPT_DIR.parent / "data" / "extended_targets_results.npz"
    np.savez(data_path, **results)
    print(f"Results saved to {data_path}")


if __name__ == "__main__":
    main()
