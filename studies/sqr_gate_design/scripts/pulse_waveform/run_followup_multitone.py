"""
Follow-up study: Optimized multitone SQR gate design.

Investigates whether genuinely optimized multitone waveforms can serve as
practical gate primitives that outperform the single-tone / common-envelope
baselines in speed and fidelity.

Phases covered:
  B — Reproduce baselines (single-tone Gaussian, one-segment multitone)
  C — Independent-tone multitone with optimized amplitudes/phases/detunings
  D — Segment-wise multitone and echo-like constructions
  E — GRAPE comparison at matching χT values
  F — Metric aggregation for synthesis

Usage:
    python scripts/run_followup_multitone.py

Output:
    data/followup_multitone_results.npz
"""
from __future__ import annotations

import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.optimize import minimize
from scipy.linalg import block_diag

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import runtime_compat  # noqa: F401

from common import (
    CHI, N_CAV, N_TR, DT,
    build_frame, build_model, duration_from_chi_t,
    target_qubit_unitary, z_corrected_target_fidelity,
    spectator_z_fidelity, spectator_transverse_error,
    identity_fidelity_with_z, extract_branch_unitaries, extract_leakage,
)
from cqed_sim.core.frequencies import (
    carrier_for_transition_frequency,
    manifold_transition_frequency,
)
from cqed_sim.sequence import SequenceCompiler
from cqed_sim.sim import SimulationConfig, prepare_simulation

# ===================================================================
# Scan parameters
# ===================================================================
CHI_T_VALUES = np.array([0.5, 1.0, 1.5, 2.0, 3.0, 5.0], dtype=float)
TARGET_BRANCHES = [0, 1, 2]
THETA_VALUES = [np.pi / 2, np.pi]
PHI_VALUES = [0.0, np.pi / 4, np.pi / 2]
LOGICAL_N_VALUES = [3, 4]
SIGMA_FRACTION = 1.0 / 6.0
TONE_COUNTS = [2, 3, 4]       # number of active tones for independent-tone study
SEGMENT_COUNTS = [1, 2, 3]    # for segment-wise study

# Representative case for detailed optimization
REP_N0 = 1
REP_THETA = np.pi
REP_PHI = 0.0
REP_LOGICAL_N = 4

# Optimization settings
OPT_MAXITER = 200
OPT_METHOD = "Nelder-Mead"


# ===================================================================
# Local pulse dataclass (reuses existing pattern from extended study)
# ===================================================================
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


@dataclass(frozen=True)
class MultitoneTone:
    manifold: int
    omega_rad_s: float
    amp_rad_s: float
    phase_rad: float


# ===================================================================
# Envelope helpers
# ===================================================================
def gaussian_envelope(t_rel, sigma=SIGMA_FRACTION):
    return np.exp(-0.5 * ((t_rel - 0.5) / sigma) ** 2).astype(np.complex128)


def gaussian_area(sigma=SIGMA_FRACTION, n_pts=4097):
    grid = np.linspace(0.0, 1.0, n_pts)
    env = np.real(gaussian_envelope(grid, sigma))
    trapezoid = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    return float(trapezoid(env, grid))


def normalized_gaussian(t_rel, sigma=SIGMA_FRACTION):
    base = gaussian_envelope(t_rel, sigma)
    area = gaussian_area(sigma)
    return base / area if abs(area) > 1e-12 else base


def cosine_squared_envelope(t_rel):
    return (2.0 * np.cos(np.pi * (t_rel - 0.5)) ** 2).astype(np.complex128)


def multitone_modulated_envelope(t_rel, duration_s, tone_specs):
    """Envelope with multiple independent tones and shared Gaussian window."""
    env = normalized_gaussian(t_rel)
    t_abs = t_rel * duration_s
    coeff = np.zeros_like(t_abs, dtype=np.complex128)
    for spec in tone_specs:
        coeff += spec.amp_rad_s * np.exp(
            1j * spec.phase_rad
        ) * np.exp(1j * spec.omega_rad_s * t_abs)
    return env * coeff


def smooth_basis_modulated_envelope(t_rel, duration_s, tone_specs, basis_coeffs):
    """Each tone gets a smooth time-dependent amplitude via Fourier basis.

    basis_coeffs: list of arrays, one per tone.
    Each array has shape (n_basis,) giving coefficients c_{m,k}
    for basis functions b_k(t) = cos(k * pi * t_rel) (k=0,1,...).
    The effective amplitude for tone m at time t is:
        A_m(t) = sum_k c_{m,k} * cos(k * pi * t_rel)
    """
    window = normalized_gaussian(t_rel)
    t_abs = t_rel * duration_s
    result = np.zeros_like(t_abs, dtype=np.complex128)
    for i, spec in enumerate(tone_specs):
        # Build time-dependent amplitude
        coeffs = basis_coeffs[i]
        amp_t = np.zeros_like(t_rel, dtype=float)
        for k, ck in enumerate(coeffs):
            amp_t += ck * np.cos(k * np.pi * t_rel)
        result += amp_t * np.exp(
            1j * spec.phase_rad
        ) * np.exp(1j * spec.omega_rad_s * t_abs)
    return window * result


# ===================================================================
# Pulse builders
# ===================================================================
def build_single_tone_gaussian(model, frame, n0, theta, phi, T):
    """Single-tone Gaussian baseline."""
    omega = manifold_transition_frequency(model, int(n0), frame)
    carrier = carrier_for_transition_frequency(omega)
    amp = float(theta) / (2.0 * T)
    pulse = Pulse(
        channel="q", t0=0.0, duration=T,
        envelope=lambda tr: normalized_gaussian(tr),
        carrier=carrier, phase=phi, amp=amp,
        label="single_tone_gaussian",
    )
    return [pulse], {"q": "qubit"}, T


def build_cosine_squared_pulse(model, frame, n0, theta, phi, T):
    """Cosine-squared (Hann) baseline."""
    omega = manifold_transition_frequency(model, int(n0), frame)
    carrier = carrier_for_transition_frequency(omega)
    amp = float(theta) / (2.0 * T)
    pulse = Pulse(
        channel="q", t0=0.0, duration=T,
        envelope=cosine_squared_envelope,
        carrier=carrier, phase=phi, amp=amp,
        label="cosine_squared",
    )
    return [pulse], {"q": "qubit"}, T


def build_independent_tone_multitone(model, frame, logical_n, n0, theta, phi, T,
                                      amp_ratios=None, phases=None, detunings=None):
    """Multitone with independent amplitudes, phases, and optional detunings.

    Parameters
    ----------
    amp_ratios : array of shape (logical_n,) or None
        Relative amplitude for each tone. Target branch gets base amplitude * amp_ratios[n0].
        If None, only target branch is active (baseline behavior).
    phases : array of shape (logical_n,) or None
        Phase for each tone (rad). If None, target branch gets phi, others 0.
    detunings : array of shape (logical_n,) or None
        Detuning from branch resonance (rad/s). If None, all zero.
    """
    base_amp = float(theta) / (2.0 * T)
    tone_specs = []
    for n in range(int(logical_n)):
        omega_n = manifold_transition_frequency(model, n, frame)
        carrier_n = carrier_for_transition_frequency(omega_n)
        if detunings is not None:
            carrier_n += detunings[n]
        if amp_ratios is not None:
            tone_amp = base_amp * amp_ratios[n]
        else:
            tone_amp = base_amp if n == int(n0) else 0.0
        if phases is not None:
            tone_phase = phases[n]
        else:
            tone_phase = phi if n == int(n0) else 0.0
        tone_specs.append(MultitoneTone(
            manifold=n,
            omega_rad_s=carrier_n,
            amp_rad_s=tone_amp,
            phase_rad=tone_phase,
        ))

    dur = float(T)
    ts = list(tone_specs)

    def envelope(t_rel):
        return multitone_modulated_envelope(
            np.asarray(t_rel, dtype=float), dur, ts)

    pulse = Pulse(
        channel="q", t0=0.0, duration=T,
        envelope=envelope, carrier=0.0, phase=0.0, amp=1.0,
        label="independent_tone_multitone",
    )
    return [pulse], {"q": "qubit"}, T


def build_segmented_multitone(model, frame, logical_n, n0, theta, phi, T_total,
                               n_segments, segment_params_list):
    """Multi-segment multitone: each segment has its own tone parameters.

    segment_params_list: list of dicts, one per segment.
    Each dict has keys 'amp_ratios', 'phases', 'detunings' (arrays of length logical_n).
    Segments equally divide T_total.
    """
    seg_duration = T_total / n_segments
    base_amp = float(theta) / (2.0 * T_total)  # total rotation spread over all segments
    pulses = []
    for seg_idx in range(n_segments):
        params = segment_params_list[seg_idx]
        amp_ratios = params.get("amp_ratios")
        seg_phases = params.get("phases")
        seg_detunings = params.get("detunings")

        tone_specs = []
        for n in range(int(logical_n)):
            omega_n = manifold_transition_frequency(model, n, frame)
            carrier_n = carrier_for_transition_frequency(omega_n)
            if seg_detunings is not None:
                carrier_n += seg_detunings[n]
            if amp_ratios is not None:
                tone_amp = base_amp * n_segments * amp_ratios[n]
            else:
                tone_amp = base_amp * n_segments if n == int(n0) else 0.0
            if seg_phases is not None:
                tone_phase = seg_phases[n]
            else:
                tone_phase = phi if n == int(n0) else 0.0
            tone_specs.append(MultitoneTone(
                manifold=n, omega_rad_s=carrier_n,
                amp_rad_s=tone_amp, phase_rad=tone_phase,
            ))

        dur = seg_duration
        ts = list(tone_specs)
        t0 = seg_idx * seg_duration

        def make_env(dur_s, specs):
            def env(t_rel):
                return multitone_modulated_envelope(
                    np.asarray(t_rel, dtype=float), dur_s, specs)
            return env

        pulse = Pulse(
            channel="q", t0=t0, duration=seg_duration,
            envelope=make_env(dur, ts),
            carrier=0.0, phase=0.0, amp=1.0,
            label=f"segment_{seg_idx}",
        )
        pulses.append(pulse)

    return pulses, {"q": "qubit"}, T_total


def build_smooth_basis_multitone(model, frame, logical_n, n0, theta, phi, T,
                                  n_basis, all_coeffs, all_phases, detunings=None):
    """Smooth-basis multitone: each tone has time-varying amplitude via cosine basis.

    all_coeffs: array (logical_n, n_basis) — basis coefficients per tone
    all_phases: array (logical_n,) — phases per tone
    detunings:  array (logical_n,) or None — detunings per tone
    """
    tone_specs = []
    basis_coeffs = []
    for n in range(int(logical_n)):
        omega_n = manifold_transition_frequency(model, n, frame)
        carrier_n = carrier_for_transition_frequency(omega_n)
        if detunings is not None:
            carrier_n += detunings[n]
        tone_specs.append(MultitoneTone(
            manifold=n, omega_rad_s=carrier_n,
            amp_rad_s=1.0,  # amplitude handled by basis expansion
            phase_rad=all_phases[n],
        ))
        basis_coeffs.append(all_coeffs[n])

    dur = float(T)
    ts = list(tone_specs)
    bc = [np.array(c, dtype=float) for c in basis_coeffs]

    def envelope(t_rel):
        return smooth_basis_modulated_envelope(
            np.asarray(t_rel, dtype=float), dur, ts, bc)

    pulse = Pulse(
        channel="q", t0=0.0, duration=T,
        envelope=envelope, carrier=0.0, phase=0.0, amp=1.0,
        label="smooth_basis_multitone",
    )
    return [pulse], {"q": "qubit"}, T


# ===================================================================
# Simulation and metric helpers
# ===================================================================
def simulate_and_extract(model, frame, pulses, drive_ops, logical_n, T_total):
    """Simulate all logical basis inputs and return full-space propagator."""
    compiler = SequenceCompiler(dt=DT)
    compiled = compiler.compile(pulses, t_end=float(T_total + 4.0 * DT))
    config = SimulationConfig(frame=frame, store_states=False)
    session = prepare_simulation(model, compiled, drive_ops, config=config, e_ops={})

    full_dim = int(model.n_tr) * int(model.n_cav)
    full_op = np.eye(full_dim, dtype=np.complex128)
    final_states = []
    for n in range(int(logical_n)):
        for q in (0, 1):
            initial = model.basis_state(q, n)
            result = session.run(initial)
            final_states.append(result.final_state)
            full_op[:, q * int(model.n_cav) + n] = (
                result.final_state.full().flatten()
            )
    return full_op, final_states


def compute_all_metrics(full_op, final_states, model, logical_n, n0, theta, phi):
    """Compute all branch + logical metrics for a simulation result."""
    blocks = extract_branch_unitaries(final_states, model, int(logical_n))
    target_gate = target_qubit_unitary(theta, phi)

    # Branch metrics
    target_fid, alpha_opt = z_corrected_target_fidelity(
        blocks[int(n0)], target_gate
    )
    branch_true = np.zeros(int(logical_n))
    branch_cphase = np.zeros(int(logical_n))
    branch_z_phase = np.zeros(int(logical_n))
    branch_transverse = np.zeros(int(logical_n))

    for n, block in enumerate(blocks):
        if n == int(n0):
            branch_true[n] = target_fid
            branch_cphase[n] = target_fid
            branch_z_phase[n] = alpha_opt
        else:
            branch_true[n] = identity_fidelity_with_z(block, alpha_opt)
            sf, sp = spectator_z_fidelity(block)
            branch_cphase[n] = sf
            branch_z_phase[n] = sp
            branch_transverse[n] = spectator_transverse_error(block)

    # Logical indices
    indices = []
    for n in range(int(logical_n)):
        indices.extend([n, int(model.n_cav) + n])
    restricted = full_op[np.ix_(indices, indices)]

    # Target operator
    I2 = np.eye(2, dtype=np.complex128)
    R = target_qubit_unitary(theta, phi)
    tgt_blocks = [R if n == int(n0) else I2 for n in range(int(logical_n))]
    target_op = block_diag(*tgt_blocks)
    dim = float(target_op.shape[0])

    # Strict logical fidelity
    strict_fid = abs(np.trace(target_op.conj().T @ restricted)) ** 2 / (dim * dim)

    # Block-phase-relaxed fidelity
    overlaps = []
    for bn in range(int(logical_n)):
        actual_block = restricted[2*bn:2*bn+2, 2*bn:2*bn+2]
        ideal_block = R if bn == int(n0) else I2
        overlaps.append(np.trace(ideal_block.conj().T @ actual_block))
    block_fid = (sum(abs(o) for o in overlaps)) ** 2 / (dim * dim)

    # Same-block population
    sb_vals = []
    for n in range(int(logical_n)):
        for q in (0, 1):
            col_idx = q * int(model.n_cav) + n
            col = full_op[:, col_idx]
            sb = abs(col[n]) ** 2 + abs(col[int(model.n_cav) + n]) ** 2
            sb_vals.append(float(sb))
    sb_arr = np.array(sb_vals)

    # Leakage
    leakage = extract_leakage(final_states, model, int(logical_n))

    # Spectator aggregates
    spec_mask = np.arange(int(logical_n)) != int(n0)
    spec_phase_spread = float(np.ptp(branch_z_phase[spec_mask])) if np.any(spec_mask) else 0.0
    spec_max_trans = float(np.max(branch_transverse[spec_mask])) if np.any(spec_mask) else 0.0

    return {
        "target_branch_fid": float(target_fid),
        "branch_true_mean": float(np.mean(branch_true)),
        "branch_cphase_mean": float(np.mean(branch_cphase)),
        "strict_logical_fid": float(np.clip(strict_fid, 0, 1)),
        "block_phase_relaxed_fid": float(np.clip(block_fid, 0, 1)),
        "same_block_pop_mean": float(np.mean(sb_arr)),
        "same_block_pop_min": float(np.min(sb_arr)),
        "leakage_mean": float(np.mean(leakage)),
        "leakage_max": float(np.max(leakage)),
        "spectator_phase_spread": spec_phase_spread,
        "spectator_max_transverse": spec_max_trans,
        "branch_z_phases": branch_z_phase.copy(),
        "branch_cphase": branch_cphase.copy(),
        "branch_true": branch_true.copy(),
    }


# ===================================================================
# Optimization wrappers
# ===================================================================
def optimize_independent_tone(model, frame, logical_n, n0, theta, phi, T,
                               objective="block_phase_relaxed_fid"):
    """Optimize independent per-tone amplitudes and phases.

    Free parameters: amp_ratios[n] for n != n0, phases[n] for n != n0.
    Target branch gets amp=1.0, phase=phi.
    """
    n_free = int(logical_n) - 1  # non-target branches
    # x = [amp_0, ..., amp_{n_free-1}, phase_0, ..., phase_{n_free-1}]
    x0 = np.zeros(2 * n_free)  # start with zero spectator amplitudes

    non_target_indices = [n for n in range(int(logical_n)) if n != int(n0)]

    def cost(x):
        amp_ratios = np.zeros(int(logical_n))
        phases_arr = np.zeros(int(logical_n))
        amp_ratios[int(n0)] = 1.0
        phases_arr[int(n0)] = float(phi)
        for i, n in enumerate(non_target_indices):
            amp_ratios[n] = x[i]
            phases_arr[n] = x[n_free + i]
        try:
            pulses, drive_ops, T_tot = build_independent_tone_multitone(
                model, frame, logical_n, n0, theta, phi, T,
                amp_ratios=amp_ratios, phases=phases_arr
            )
            full_op, final_states = simulate_and_extract(
                model, frame, pulses, drive_ops, logical_n, T_tot
            )
            metrics = compute_all_metrics(
                full_op, final_states, model, logical_n, n0, theta, phi
            )
            return -metrics[objective]
        except Exception:
            return 0.0  # worst case fidelity

    result = minimize(cost, x0, method=OPT_METHOD,
                      options={"maxiter": OPT_MAXITER, "xatol": 1e-5, "fatol": 1e-6})

    # Reconstruct best parameters
    amp_ratios = np.zeros(int(logical_n))
    phases_arr = np.zeros(int(logical_n))
    amp_ratios[int(n0)] = 1.0
    phases_arr[int(n0)] = float(phi)
    for i, n in enumerate(non_target_indices):
        amp_ratios[n] = result.x[i]
        phases_arr[n] = result.x[n_free + i]

    return amp_ratios, phases_arr, result


def optimize_detuned_multitone(model, frame, logical_n, n0, theta, phi, T,
                                objective="block_phase_relaxed_fid"):
    """Optimize per-tone amplitudes, phases, AND detunings."""
    n_tones = int(logical_n)
    non_target = [n for n in range(n_tones) if n != int(n0)]
    n_free = len(non_target)
    # x = [amps..., phases..., detunings_all...]
    # detunings for all tones including target
    x0 = np.zeros(2 * n_free + n_tones)

    chi_abs = abs(CHI)

    def cost(x):
        amp_ratios = np.zeros(n_tones)
        phases_arr = np.zeros(n_tones)
        detunings = np.zeros(n_tones)
        amp_ratios[int(n0)] = 1.0
        phases_arr[int(n0)] = float(phi)
        for i, n in enumerate(non_target):
            amp_ratios[n] = x[i]
            phases_arr[n] = x[n_free + i]
        for n in range(n_tones):
            detunings[n] = x[2 * n_free + n] * chi_abs * 0.1  # scale: 10% of chi
        try:
            pulses, drive_ops, T_tot = build_independent_tone_multitone(
                model, frame, logical_n, n0, theta, phi, T,
                amp_ratios=amp_ratios, phases=phases_arr, detunings=detunings
            )
            full_op, final_states = simulate_and_extract(
                model, frame, pulses, drive_ops, logical_n, T_tot
            )
            metrics = compute_all_metrics(
                full_op, final_states, model, logical_n, n0, theta, phi
            )
            return -metrics[objective]
        except Exception:
            return 0.0

    result = minimize(cost, x0, method=OPT_METHOD,
                      options={"maxiter": OPT_MAXITER, "xatol": 1e-5, "fatol": 1e-6})

    amp_ratios = np.zeros(n_tones)
    phases_arr = np.zeros(n_tones)
    detunings = np.zeros(n_tones)
    amp_ratios[int(n0)] = 1.0
    phases_arr[int(n0)] = float(phi)
    for i, n in enumerate(non_target):
        amp_ratios[n] = result.x[i]
        phases_arr[n] = result.x[n_free + i]
    for n in range(n_tones):
        detunings[n] = result.x[2 * n_free + n] * abs(CHI) * 0.1

    return amp_ratios, phases_arr, detunings, result


def optimize_2segment_multitone(model, frame, logical_n, n0, theta, phi, T,
                                 objective="block_phase_relaxed_fid"):
    """Optimize a 2-segment multitone with independent parameters per segment."""
    n_tones = int(logical_n)
    non_target = [n for n in range(n_tones) if n != int(n0)]
    n_free = len(non_target)
    # Each segment: amps (n_free) + phases (n_tones)
    # Total: 2 * (n_free + n_tones)
    n_per_seg = n_free + n_tones
    x0 = np.zeros(2 * n_per_seg)
    # Initialize: segment 1 gets target amp=1 phase=phi; segment 2 same but phase-flipped
    x0[n_free + int(n0)] = float(phi)            # seg 0 target phase
    x0[n_per_seg + n_free + int(n0)] = -float(phi)  # seg 1 target phase (echo-like)

    def cost(x):
        seg_params = []
        for seg in range(2):
            offset = seg * n_per_seg
            amp_ratios = np.zeros(n_tones)
            phases_arr = np.zeros(n_tones)
            amp_ratios[int(n0)] = 1.0
            for i, n in enumerate(non_target):
                amp_ratios[n] = x[offset + i]
            for n in range(n_tones):
                phases_arr[n] = x[offset + n_free + n]
            seg_params.append({
                "amp_ratios": amp_ratios,
                "phases": phases_arr,
                "detunings": None,
            })
        try:
            pulses, drive_ops, T_tot = build_segmented_multitone(
                model, frame, logical_n, n0, theta, phi, T,
                n_segments=2, segment_params_list=seg_params
            )
            full_op, final_states = simulate_and_extract(
                model, frame, pulses, drive_ops, logical_n, T_tot
            )
            metrics = compute_all_metrics(
                full_op, final_states, model, logical_n, n0, theta, phi
            )
            return -metrics[objective]
        except Exception:
            return 0.0

    result = minimize(cost, x0, method=OPT_METHOD,
                      options={"maxiter": OPT_MAXITER, "xatol": 1e-5, "fatol": 1e-6})
    return result


def optimize_smooth_basis(model, frame, logical_n, n0, theta, phi, T,
                           n_basis=3, objective="block_phase_relaxed_fid"):
    """Optimize smooth-basis multitone: each tone has cos-basis amplitude modulation."""
    n_tones = int(logical_n)
    non_target = [n for n in range(n_tones) if n != int(n0)]
    n_free_t = len(non_target)
    # Parameters: per tone: n_basis coefficients + 1 phase
    # Target tone: n_basis coeffs (phase fixed)
    # Free tones: n_basis coeffs + 1 phase each
    # x = [target_c0, ..., target_c{nb-1},
    #       free0_c0, ..., free0_c{nb-1}, free0_phase,
    #       free1_c0, ..., ...]
    dim_x = n_basis + n_free_t * (n_basis + 1)
    x0 = np.zeros(dim_x)
    # Initialize target tone to have constant amplitude
    base_amp = float(theta) / (2.0 * T)
    x0[0] = base_amp  # c_0 for target = constant amplitude

    def cost(x):
        all_coeffs = np.zeros((n_tones, n_basis))
        all_phases = np.zeros(n_tones)
        # Target tone
        all_coeffs[int(n0)] = x[:n_basis]
        all_phases[int(n0)] = float(phi)
        # Free tones
        offset = n_basis
        for i, n in enumerate(non_target):
            all_coeffs[n] = x[offset:offset + n_basis]
            all_phases[n] = x[offset + n_basis]
            offset += n_basis + 1
        try:
            pulses, drive_ops, T_tot = build_smooth_basis_multitone(
                model, frame, logical_n, n0, theta, phi, T,
                n_basis=n_basis, all_coeffs=all_coeffs, all_phases=all_phases
            )
            full_op, final_states = simulate_and_extract(
                model, frame, pulses, drive_ops, logical_n, T_tot
            )
            metrics = compute_all_metrics(
                full_op, final_states, model, logical_n, n0, theta, phi
            )
            return -metrics[objective]
        except Exception:
            return 0.0

    result = minimize(cost, x0, method=OPT_METHOD,
                      options={"maxiter": OPT_MAXITER, "xatol": 1e-5, "fatol": 1e-6})

    # Reconstruct
    all_coeffs = np.zeros((n_tones, n_basis))
    all_phases = np.zeros(n_tones)
    all_coeffs[int(n0)] = result.x[:n_basis]
    all_phases[int(n0)] = float(phi)
    offset = n_basis
    for i, n in enumerate(non_target):
        all_coeffs[n] = result.x[offset:offset + n_basis]
        all_phases[n] = result.x[offset + n_basis]
        offset += n_basis + 1
    return all_coeffs, all_phases, result


# ===================================================================
# GRAPE benchmark (reuse existing pattern)
# ===================================================================
def run_grape_benchmark(model, frame, logical_n, n0, theta, phi, T, cphase=True):
    """Run GRAPE for comparison. Returns fidelity dict."""
    from cqed_sim import (
        GrapeConfig, GrapeSolver, ModelControlChannelSpec,
        PiecewiseConstantTimeGrid, UnitaryObjective,
        build_control_problem_from_model,
    )
    from cqed_sim.unitary_synthesis import Subspace

    n_slices = 48
    dt_grape = T / n_slices
    amp_bound = 2 * np.pi * 50e6

    # Build subspace indices
    indices = []
    labels = []
    for n in range(int(logical_n)):
        indices.append(0 * int(model.n_cav) + n)
        indices.append(1 * int(model.n_cav) + n)
        labels.append(f"|g,{n}>")
        labels.append(f"|e,{n}>")
    sub = Subspace(
        full_dim=int(model.n_tr) * int(model.n_cav),
        indices=tuple(indices), labels=tuple(labels)
    )

    # Build target
    I2 = np.eye(2, dtype=np.complex128)
    R = target_qubit_unitary(theta, phi)
    tgt_blocks = [R if n == int(n0) else I2 for n in range(int(logical_n))]
    target_sub = block_diag(*tgt_blocks)

    obj_kwargs = dict(
        target_operator=target_sub, subspace=sub,
        ignore_global_phase=True,
    )
    if cphase:
        phase_blocks = tuple((2 * n, 2 * n + 1) for n in range(int(logical_n)))
        obj_kwargs["phase_blocks"] = phase_blocks

    problem = build_control_problem_from_model(
        model, frame=frame,
        time_grid=PiecewiseConstantTimeGrid.uniform(steps=n_slices, dt_s=dt_grape),
        channel_specs=(
            ModelControlChannelSpec(
                name="qubit_I", target="qubit", quadratures=("I",),
                amplitude_bounds=(-amp_bound, amp_bound),
            ),
            ModelControlChannelSpec(
                name="qubit_Q", target="qubit", quadratures=("Q",),
                amplitude_bounds=(-amp_bound, amp_bound),
            ),
        ),
        objectives=(UnitaryObjective(**obj_kwargs),),
    )

    config = GrapeConfig(maxiter=300, seed=42)
    result = GrapeSolver(config).solve(problem)

    if "nominal_fidelity" in result.metrics:
        fid = result.metrics["nominal_fidelity"]
    elif "fidelity" in result.metrics:
        fid = result.metrics["fidelity"]
    else:
        fid = 1.0 - result.objective_value if result.objective_value <= 1.0 else 0.0
    return {"fidelity": float(fid), "converged": result.success}


# ===================================================================
# Main study runner
# ===================================================================
def run_study():
    """Execute the full follow-up study."""
    data_dir = SCRIPT_DIR.parent / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    model = build_model()
    frame = build_frame(model)

    results = {}
    results["chi_t_values"] = CHI_T_VALUES
    results["target_branches"] = np.array(TARGET_BRANCHES)
    results["theta_values"] = np.array(THETA_VALUES)
    results["phi_values"] = np.array(PHI_VALUES)
    results["logical_n_values"] = np.array(LOGICAL_N_VALUES)

    n_chi = len(CHI_T_VALUES)
    t_start = time.time()

    # ----------------------------------------------------------------
    # Phase B: Reproduce baselines
    # ----------------------------------------------------------------
    print("=" * 60)
    print("PHASE B: Baseline reproduction")
    print("=" * 60)

    families_B = ["single_tone_gaussian", "cosine_squared", "multitone_baseline"]
    shape_B = (len(families_B), n_chi)
    for metric in ["target_branch_fid", "branch_true_mean", "branch_cphase_mean",
                    "strict_logical_fid", "block_phase_relaxed_fid",
                    "same_block_pop_mean", "leakage_max",
                    "spectator_phase_spread", "spectator_max_transverse"]:
        results[f"baseline_{metric}"] = np.zeros(shape_B)

    for fi, family in enumerate(families_B):
        for ci, chi_t in enumerate(CHI_T_VALUES):
            T = duration_from_chi_t(chi_t)
            if family == "single_tone_gaussian":
                pulses, dops, T_tot = build_single_tone_gaussian(
                    model, frame, REP_N0, REP_THETA, REP_PHI, T)
            elif family == "cosine_squared":
                pulses, dops, T_tot = build_cosine_squared_pulse(
                    model, frame, REP_N0, REP_THETA, REP_PHI, T)
            else:  # multitone_baseline — same as single-tone but via multitone API
                pulses, dops, T_tot = build_independent_tone_multitone(
                    model, frame, REP_LOGICAL_N, REP_N0, REP_THETA, REP_PHI, T)

            full_op, states = simulate_and_extract(
                model, frame, pulses, dops, REP_LOGICAL_N, T_tot)
            m = compute_all_metrics(
                full_op, states, model, REP_LOGICAL_N, REP_N0, REP_THETA, REP_PHI)

            for metric in ["target_branch_fid", "branch_true_mean", "branch_cphase_mean",
                           "strict_logical_fid", "block_phase_relaxed_fid",
                           "same_block_pop_mean", "leakage_max",
                           "spectator_phase_spread", "spectator_max_transverse"]:
                results[f"baseline_{metric}"][fi, ci] = m[metric]

            print(f"  {family:30s} χT/2π={chi_t:4.1f}  "
                  f"F_strict={m['strict_logical_fid']:.4f}  "
                  f"F_block={m['block_phase_relaxed_fid']:.4f}  "
                  f"F_cphase={m['branch_cphase_mean']:.4f}")

    results["baseline_families"] = np.array(families_B, dtype=object)
    print(f"Phase B complete ({time.time() - t_start:.1f}s)\n")

    # ----------------------------------------------------------------
    # Phase C: Optimized independent-tone multitone
    # ----------------------------------------------------------------
    print("=" * 60)
    print("PHASE C: Optimized independent-tone multitone")
    print("=" * 60)

    # C1: Independent amplitudes/phases
    shape_C = (n_chi,)
    for prefix in ["opt_indep", "opt_detuned", "opt_smooth"]:
        for metric in ["target_branch_fid", "branch_true_mean", "branch_cphase_mean",
                        "strict_logical_fid", "block_phase_relaxed_fid",
                        "same_block_pop_mean", "leakage_max",
                        "spectator_phase_spread", "spectator_max_transverse"]:
            results[f"{prefix}_{metric}"] = np.zeros(shape_C)
        results[f"{prefix}_opt_cost"] = np.zeros(shape_C)

    for ci, chi_t in enumerate(CHI_T_VALUES):
        T = duration_from_chi_t(chi_t)
        print(f"\n  χT/2π = {chi_t:.1f}")

        # C1: Independent-tone optimization
        print(f"    Optimizing independent-tone...", end=" ", flush=True)
        t0 = time.time()
        amp_r, ph_r, res = optimize_independent_tone(
            model, frame, REP_LOGICAL_N, REP_N0, REP_THETA, REP_PHI, T,
            objective="block_phase_relaxed_fid"
        )
        # Evaluate at optimum
        pulses, dops, T_tot = build_independent_tone_multitone(
            model, frame, REP_LOGICAL_N, REP_N0, REP_THETA, REP_PHI, T,
            amp_ratios=amp_r, phases=ph_r
        )
        full_op, states = simulate_and_extract(
            model, frame, pulses, dops, REP_LOGICAL_N, T_tot)
        m = compute_all_metrics(
            full_op, states, model, REP_LOGICAL_N, REP_N0, REP_THETA, REP_PHI)
        for metric in ["target_branch_fid", "branch_true_mean", "branch_cphase_mean",
                        "strict_logical_fid", "block_phase_relaxed_fid",
                        "same_block_pop_mean", "leakage_max",
                        "spectator_phase_spread", "spectator_max_transverse"]:
            results[f"opt_indep_{metric}"][ci] = m[metric]
        results["opt_indep_opt_cost"][ci] = res.fun
        print(f"F_block={m['block_phase_relaxed_fid']:.4f} "
              f"F_strict={m['strict_logical_fid']:.4f} ({time.time()-t0:.1f}s)")

        # C2: Detuned multitone optimization
        print(f"    Optimizing detuned multitone...", end=" ", flush=True)
        t0 = time.time()
        amp_r2, ph_r2, det_r, res2 = optimize_detuned_multitone(
            model, frame, REP_LOGICAL_N, REP_N0, REP_THETA, REP_PHI, T,
            objective="block_phase_relaxed_fid"
        )
        pulses, dops, T_tot = build_independent_tone_multitone(
            model, frame, REP_LOGICAL_N, REP_N0, REP_THETA, REP_PHI, T,
            amp_ratios=amp_r2, phases=ph_r2, detunings=det_r
        )
        full_op, states = simulate_and_extract(
            model, frame, pulses, dops, REP_LOGICAL_N, T_tot)
        m = compute_all_metrics(
            full_op, states, model, REP_LOGICAL_N, REP_N0, REP_THETA, REP_PHI)
        for metric in ["target_branch_fid", "branch_true_mean", "branch_cphase_mean",
                        "strict_logical_fid", "block_phase_relaxed_fid",
                        "same_block_pop_mean", "leakage_max",
                        "spectator_phase_spread", "spectator_max_transverse"]:
            results[f"opt_detuned_{metric}"][ci] = m[metric]
        results["opt_detuned_opt_cost"][ci] = res2.fun
        print(f"F_block={m['block_phase_relaxed_fid']:.4f} "
              f"F_strict={m['strict_logical_fid']:.4f} ({time.time()-t0:.1f}s)")

        # C3: Smooth-basis multitone
        print(f"    Optimizing smooth-basis...", end=" ", flush=True)
        t0 = time.time()
        sb_coeffs, sb_phases, res3 = optimize_smooth_basis(
            model, frame, REP_LOGICAL_N, REP_N0, REP_THETA, REP_PHI, T,
            n_basis=3, objective="block_phase_relaxed_fid"
        )
        pulses, dops, T_tot = build_smooth_basis_multitone(
            model, frame, REP_LOGICAL_N, REP_N0, REP_THETA, REP_PHI, T,
            n_basis=3, all_coeffs=sb_coeffs, all_phases=sb_phases
        )
        full_op, states = simulate_and_extract(
            model, frame, pulses, dops, REP_LOGICAL_N, T_tot)
        m = compute_all_metrics(
            full_op, states, model, REP_LOGICAL_N, REP_N0, REP_THETA, REP_PHI)
        for metric in ["target_branch_fid", "branch_true_mean", "branch_cphase_mean",
                        "strict_logical_fid", "block_phase_relaxed_fid",
                        "same_block_pop_mean", "leakage_max",
                        "spectator_phase_spread", "spectator_max_transverse"]:
            results[f"opt_smooth_{metric}"][ci] = m[metric]
        results["opt_smooth_opt_cost"][ci] = res3.fun
        print(f"F_block={m['block_phase_relaxed_fid']:.4f} "
              f"F_strict={m['strict_logical_fid']:.4f} ({time.time()-t0:.1f}s)")

    print(f"\nPhase C complete ({time.time() - t_start:.1f}s)\n")

    # ----------------------------------------------------------------
    # Phase D: Segmented multitone
    # ----------------------------------------------------------------
    print("=" * 60)
    print("PHASE D: Segmented multitone optimization")
    print("=" * 60)

    for metric in ["target_branch_fid", "branch_true_mean", "branch_cphase_mean",
                    "strict_logical_fid", "block_phase_relaxed_fid",
                    "same_block_pop_mean", "leakage_max",
                    "spectator_phase_spread", "spectator_max_transverse"]:
        results[f"opt_2seg_{metric}"] = np.zeros(shape_C)
    results["opt_2seg_opt_cost"] = np.zeros(shape_C)

    for ci, chi_t in enumerate(CHI_T_VALUES):
        T = duration_from_chi_t(chi_t)
        print(f"  χT/2π = {chi_t:.1f}  Optimizing 2-segment...", end=" ", flush=True)
        t0 = time.time()
        res4 = optimize_2segment_multitone(
            model, frame, REP_LOGICAL_N, REP_N0, REP_THETA, REP_PHI, T,
            objective="block_phase_relaxed_fid"
        )
        # Reconstruct and evaluate
        n_tones = REP_LOGICAL_N
        non_target = [n for n in range(n_tones) if n != REP_N0]
        n_free = len(non_target)
        n_per_seg = n_free + n_tones
        seg_params = []
        for seg in range(2):
            offset = seg * n_per_seg
            amp_ratios = np.zeros(n_tones)
            phases_arr = np.zeros(n_tones)
            amp_ratios[REP_N0] = 1.0
            for i, n in enumerate(non_target):
                amp_ratios[n] = res4.x[offset + i]
            for n in range(n_tones):
                phases_arr[n] = res4.x[offset + n_free + n]
            seg_params.append({
                "amp_ratios": amp_ratios,
                "phases": phases_arr,
                "detunings": None,
            })
        pulses, dops, T_tot = build_segmented_multitone(
            model, frame, REP_LOGICAL_N, REP_N0, REP_THETA, REP_PHI, T,
            n_segments=2, segment_params_list=seg_params
        )
        full_op, states = simulate_and_extract(
            model, frame, pulses, dops, REP_LOGICAL_N, T_tot)
        m = compute_all_metrics(
            full_op, states, model, REP_LOGICAL_N, REP_N0, REP_THETA, REP_PHI)
        for metric in ["target_branch_fid", "branch_true_mean", "branch_cphase_mean",
                        "strict_logical_fid", "block_phase_relaxed_fid",
                        "same_block_pop_mean", "leakage_max",
                        "spectator_phase_spread", "spectator_max_transverse"]:
            results[f"opt_2seg_{metric}"][ci] = m[metric]
        results["opt_2seg_opt_cost"][ci] = res4.fun
        print(f"F_block={m['block_phase_relaxed_fid']:.4f} "
              f"F_strict={m['strict_logical_fid']:.4f} ({time.time()-t0:.1f}s)")

    print(f"\nPhase D complete ({time.time() - t_start:.1f}s)\n")

    # ----------------------------------------------------------------
    # Phase E: GRAPE comparison
    # ----------------------------------------------------------------
    print("=" * 60)
    print("PHASE E: GRAPE comparison")
    print("=" * 60)

    grape_chi_t = np.array([1.0, 2.0, 3.0, 5.0])  # subset for speed
    results["grape_chi_t"] = grape_chi_t
    results["grape_cphase_fid"] = np.zeros(len(grape_chi_t))
    results["grape_true_fid"] = np.zeros(len(grape_chi_t))

    for ci, chi_t in enumerate(grape_chi_t):
        T = duration_from_chi_t(chi_t)
        print(f"  χT/2π = {chi_t:.1f}", end=" ", flush=True)

        t0 = time.time()
        r_cp = run_grape_benchmark(model, frame, REP_LOGICAL_N, REP_N0,
                                    REP_THETA, REP_PHI, T, cphase=True)
        results["grape_cphase_fid"][ci] = r_cp["fidelity"]
        print(f"GRAPE cphase F={r_cp['fidelity']:.6f}", end=" ", flush=True)

        r_tr = run_grape_benchmark(model, frame, REP_LOGICAL_N, REP_N0,
                                    REP_THETA, REP_PHI, T, cphase=False)
        results["grape_true_fid"][ci] = r_tr["fidelity"]
        print(f"true F={r_tr['fidelity']:.6f} ({time.time()-t0:.1f}s)")

    print(f"\nPhase E complete ({time.time() - t_start:.1f}s)\n")

    # ----------------------------------------------------------------
    # Phase C extended: Parameter scan (branch, angle, truncation)
    # ----------------------------------------------------------------
    print("=" * 60)
    print("PHASE C-EXT: Parameter scan over branches/angles/truncations")
    print("=" * 60)

    # Scan over target branch at representative angle
    scan_branches = TARGET_BRANCHES
    scan_chi_t = np.array([1.0, 2.0, 3.0, 5.0])
    results["scan_branches"] = np.array(scan_branches)
    results["scan_chi_t"] = scan_chi_t

    for prefix in ["baseline_gauss", "opt_indep"]:
        for metric in ["strict_logical_fid", "block_phase_relaxed_fid",
                        "branch_cphase_mean", "leakage_max"]:
            results[f"branch_scan_{prefix}_{metric}"] = np.zeros(
                (len(scan_branches), len(scan_chi_t)))

    for bi, n0 in enumerate(scan_branches):
        for ci, chi_t in enumerate(scan_chi_t):
            T = duration_from_chi_t(chi_t)
            ln = min(REP_LOGICAL_N, max(n0 + 2, 3))  # ensure enough levels

            # Baseline
            pulses, dops, T_tot = build_single_tone_gaussian(
                model, frame, n0, REP_THETA, REP_PHI, T)
            full_op, states = simulate_and_extract(
                model, frame, pulses, dops, ln, T_tot)
            m = compute_all_metrics(full_op, states, model, ln, n0, REP_THETA, REP_PHI)
            for metric in ["strict_logical_fid", "block_phase_relaxed_fid",
                            "branch_cphase_mean", "leakage_max"]:
                results[f"branch_scan_baseline_gauss_{metric}"][bi, ci] = m[metric]

            # Optimized independent-tone
            amp_r, ph_r, _ = optimize_independent_tone(
                model, frame, ln, n0, REP_THETA, REP_PHI, T,
                objective="block_phase_relaxed_fid"
            )
            pulses, dops, T_tot = build_independent_tone_multitone(
                model, frame, ln, n0, REP_THETA, REP_PHI, T,
                amp_ratios=amp_r, phases=ph_r
            )
            full_op, states = simulate_and_extract(
                model, frame, pulses, dops, ln, T_tot)
            m = compute_all_metrics(full_op, states, model, ln, n0, REP_THETA, REP_PHI)
            for metric in ["strict_logical_fid", "block_phase_relaxed_fid",
                            "branch_cphase_mean", "leakage_max"]:
                results[f"branch_scan_opt_indep_{metric}"][bi, ci] = m[metric]

            print(f"  n0={n0} χT/2π={chi_t:.1f}  "
                  f"gauss_strict={results[f'branch_scan_baseline_gauss_strict_logical_fid'][bi,ci]:.3f}  "
                  f"opt_strict={results[f'branch_scan_opt_indep_strict_logical_fid'][bi,ci]:.3f}")

    # Scan over angles at representative branch
    scan_thetas = THETA_VALUES
    scan_phis = PHI_VALUES
    results["scan_thetas"] = np.array(scan_thetas)
    results["scan_phis"] = np.array(scan_phis)

    for prefix in ["baseline_gauss", "opt_indep"]:
        for metric in ["strict_logical_fid", "block_phase_relaxed_fid"]:
            results[f"angle_scan_{prefix}_{metric}"] = np.zeros(
                (len(scan_thetas), len(scan_phis), len(scan_chi_t)))

    for ti, theta in enumerate(scan_thetas):
        for pi_idx, phi in enumerate(scan_phis):
            for ci, chi_t in enumerate(scan_chi_t):
                T = duration_from_chi_t(chi_t)
                # Baseline
                pulses, dops, T_tot = build_single_tone_gaussian(
                    model, frame, REP_N0, theta, phi, T)
                full_op, states = simulate_and_extract(
                    model, frame, pulses, dops, REP_LOGICAL_N, T_tot)
                m = compute_all_metrics(
                    full_op, states, model, REP_LOGICAL_N, REP_N0, theta, phi)
                for metric in ["strict_logical_fid", "block_phase_relaxed_fid"]:
                    results[f"angle_scan_baseline_gauss_{metric}"][ti, pi_idx, ci] = m[metric]

                # Optimized
                amp_r, ph_r, _ = optimize_independent_tone(
                    model, frame, REP_LOGICAL_N, REP_N0, theta, phi, T,
                    objective="block_phase_relaxed_fid"
                )
                pulses, dops, T_tot = build_independent_tone_multitone(
                    model, frame, REP_LOGICAL_N, REP_N0, theta, phi, T,
                    amp_ratios=amp_r, phases=ph_r
                )
                full_op, states = simulate_and_extract(
                    model, frame, pulses, dops, REP_LOGICAL_N, T_tot)
                m = compute_all_metrics(
                    full_op, states, model, REP_LOGICAL_N, REP_N0, theta, phi)
                for metric in ["strict_logical_fid", "block_phase_relaxed_fid"]:
                    results[f"angle_scan_opt_indep_{metric}"][ti, pi_idx, ci] = m[metric]

                print(f"  θ={theta/np.pi:.2f}π φ={phi/np.pi:.2f}π χT/2π={chi_t:.1f}  "
                      f"gauss={results[f'angle_scan_baseline_gauss_strict_logical_fid'][ti,pi_idx,ci]:.3f}  "
                      f"opt={results[f'angle_scan_opt_indep_strict_logical_fid'][ti,pi_idx,ci]:.3f}")

    print(f"\nPhase C-EXT complete ({time.time() - t_start:.1f}s)\n")

    # ----------------------------------------------------------------
    # Save all results
    # ----------------------------------------------------------------
    outpath = data_dir / "followup_multitone_results.npz"
    np.savez(str(outpath), **results)
    print(f"\nAll results saved to {outpath}")
    print(f"Total time: {time.time() - t_start:.1f}s")


if __name__ == "__main__":
    run_study()
