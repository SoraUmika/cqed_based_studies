"""
Simultaneous multitone SQR gate design: corrected (theta, d_lambda, d_alpha,
delta_omega) waveform parameterization.

Each tone in the multitone pulse is constructed as:

    w_n[k] = s_n * w0[k] * exp(i * phi_eff_n) * exp(i * omega_n * t_k)

where the per-tone parameters are

    s_n         = theta_n / pi  +  d_lambda_n / lambda_0(T)
    lambda_0(T) = pi / (2 * T)       [natural amplitude scale]
    phi_eff_n   = phi_n + d_alpha_n
    omega_n     = omega_det_n + delta_omega_n

A tone is included in the waveform only when |s_n| > S_THRESHOLD, i.e. when
the total amplitude scale is non-negligible — regardless of whether theta_n is
zero or not.

Three free corrections per tone: (d_lambda_n, d_alpha_n, delta_omega_n).

Gates optimised:
    SQR             — strict selective qubit rotation (strict_logical_fid)
    controlPhaseSQR — conditional-phase variant      (block_phase_relaxed_fid)

Cases:
    single-branch — target n0 = 1, spectators start idle (theta_n = 0)
    all-branch    — all Fock branches simultaneously  (theta_n = theta for all n)

Scan: chi_T_values = {0.5, 1.0, 1.5, 2.0, 3.0, 5.0}

Output:
    data/simultaneous_multitone_results.npz
"""
from __future__ import annotations

import gc
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.linalg import block_diag
from scipy.optimize import differential_evolution, minimize

SCRIPT_DIR = Path(__file__).resolve().parent
STUDY_DIR = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

import runtime_compat  # noqa: F401

from common import (
    CHI,
    DT,
    N_FOCK,
    build_frame,
    build_model,
    duration_from_chi_t,
    extract_branch_unitaries,
    extract_leakage,
    identity_fidelity_with_z,
    spectator_z_fidelity,
    target_qubit_unitary,
    z_corrected_target_fidelity,
)
from cqed_sim.core.frequencies import (
    carrier_for_transition_frequency,
    manifold_transition_frequency,
)
from cqed_sim.sequence import SequenceCompiler
from cqed_sim.sim import SimulationConfig, prepare_simulation

# ===================================================================
# Study parameters
# ===================================================================
CHI_T_VALUES = np.array([0.5, 1.0, 1.5, 2.0, 3.0, 5.0], dtype=float)
LOGICAL_N = N_FOCK        # 4 Fock levels (n = 0 .. 3)
TARGET_N0 = 1             # single-branch target Fock level
REP_THETA = np.pi         # pi rotation
REP_PHI = 0.0             # X axis
SIGMA_FRACTION = 1.0 / 6.0
S_THRESHOLD = 1e-4        # |s_n| < S_THRESHOLD -> tone omitted

OPT_MAXITER_NM = 500
OPT_MAXITER_DE = 150

SAVE_PATH = STUDY_DIR / "data" / "simultaneous_multitone_results.npz"


# ===================================================================
# Pulse dataclass (matches existing study pattern)
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
class _ToneDef:
    """Internal: resolved (carrier_freq, amplitude, phase) for one tone."""
    omega_carrier: float  # rad/s, frequency of modulation in the envelope
    amp_rad_s: float      # rad/s, Rabi amplitude = s_n * lambda_0(T)
    phase_rad: float      # rad, phi_eff_n = phi_n + d_alpha_n


# ===================================================================
# Envelope helpers
# ===================================================================
def _gaussian_envelope(t_rel: np.ndarray, sigma: float = SIGMA_FRACTION) -> np.ndarray:
    return np.exp(-0.5 * ((t_rel - 0.5) / sigma) ** 2).astype(np.complex128)


def _gaussian_area(sigma: float = SIGMA_FRACTION, n_pts: int = 4097) -> float:
    grid = np.linspace(0.0, 1.0, n_pts)
    env = np.real(_gaussian_envelope(grid, sigma))
    trapezoid = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    return float(trapezoid(env, grid))


def _normalized_gaussian(t_rel: np.ndarray, sigma: float = SIGMA_FRACTION) -> np.ndarray:
    base = _gaussian_envelope(t_rel, sigma)
    area = _gaussian_area(sigma)
    return base / area if abs(area) > 1e-12 else base


def _multitone_envelope(t_rel: np.ndarray, duration_s: float,
                        tone_defs: list[_ToneDef]) -> np.ndarray:
    """Gaussian-windowed sum of modulated tones (complex IQ)."""
    window = _normalized_gaussian(t_rel)
    t_abs = t_rel * duration_s
    coeff = np.zeros_like(t_abs, dtype=np.complex128)
    for td in tone_defs:
        coeff += (td.amp_rad_s
                  * np.exp(1j * td.phase_rad)
                  * np.exp(1j * td.omega_carrier * t_abs))
    return window * coeff


# ===================================================================
# New parameterization: waveform builder
# ===================================================================
def lambda_0(T: float) -> float:
    """Natural amplitude scale lambda_0 = pi / (2*T).

    When s_n = 1 (i.e. theta_n = pi, d_lambda_n = 0) and the envelope is
    normalized_gaussian, the resulting pulse implements exactly a pi rotation.
    """
    return np.pi / (2.0 * float(T))


def build_corrected_multitone(
    model,
    frame,
    T: float,
    tone_thetas: np.ndarray,    # (n_tones,) nominal rotation angles (rad)
    tone_phis: np.ndarray,      # (n_tones,) nominal phases (rad)
    tone_det_freqs: np.ndarray, # (n_tones,) nominal carrier frequencies (rad/s)
    d_lambdas: np.ndarray,      # (n_tones,) amplitude corrections (rad/s)
    d_alphas: np.ndarray,       # (n_tones,) phase corrections (rad)
    delta_omegas: np.ndarray,   # (n_tones,) frequency corrections (rad/s)
):
    """Build a corrected multitone Gaussian pulse.

    For tone n:
        s_n         = tone_thetas[n] / pi  +  d_lambdas[n] / lambda_0(T)
        amp_n       = s_n * lambda_0(T)    = tone_thetas[n]/(2T) + d_lambdas[n]
        phi_eff_n   = tone_phis[n] + d_alphas[n]
        omega_n     = tone_det_freqs[n] + delta_omegas[n]

    Tone n is included in the waveform only when |s_n| > S_THRESHOLD.
    """
    lam0 = lambda_0(T)
    n_tones = len(tone_thetas)

    active_tones: list[_ToneDef] = []
    for n in range(n_tones):
        s_n = tone_thetas[n] / np.pi + d_lambdas[n] / lam0
        if abs(s_n) < S_THRESHOLD:
            continue  # omit this tone
        amp_n = s_n * lam0           # = tone_thetas[n]/(2T) + d_lambdas[n]
        phi_eff_n = tone_phis[n] + d_alphas[n]
        omega_n = tone_det_freqs[n] + delta_omegas[n]
        active_tones.append(_ToneDef(omega_carrier=omega_n,
                                     amp_rad_s=amp_n,
                                     phase_rad=phi_eff_n))

    # Null pulse if all tones are omitted (shouldn't happen in practice)
    if not active_tones:
        active_tones.append(_ToneDef(0.0, 0.0, 0.0))

    dur = float(T)
    tds = list(active_tones)

    def envelope(t_rel: np.ndarray) -> np.ndarray:
        return _multitone_envelope(np.asarray(t_rel, dtype=float), dur, tds)

    pulse = Pulse(channel="q", t0=0.0, duration=T,
                  envelope=envelope, carrier=0.0, phase=0.0, amp=1.0,
                  label="corrected_multitone")
    return [pulse], {"q": "qubit"}, T


def get_branch_carrier_freqs(model, frame, n_tones: int) -> np.ndarray:
    """Carrier frequency for each Fock branch in the rotating frame."""
    freqs = []
    for n in range(n_tones):
        omega_n = manifold_transition_frequency(model, n, frame)
        freqs.append(carrier_for_transition_frequency(omega_n))
    return np.array(freqs)


# ===================================================================
# Simulation helper
# ===================================================================
def simulate_and_extract(model, frame, pulses, drive_ops, logical_n: int,
                         T_total: float):
    """Simulate all logical basis inputs and return (full_op, final_states)."""
    compiler = SequenceCompiler(dt=DT)
    compiled = compiler.compile(pulses, t_end=float(T_total + 4.0 * DT))
    config = SimulationConfig(frame=frame, store_states=False)
    session = prepare_simulation(model, compiled, drive_ops, config=config, e_ops={})

    full_dim = int(model.n_tr) * int(model.n_cav)
    full_op = np.eye(full_dim, dtype=np.complex128)
    final_states = []
    for n in range(logical_n):
        for q in (0, 1):
            initial = model.basis_state(q, n)
            result = session.run(initial)
            final_states.append(result.final_state)
            full_op[:, q * int(model.n_cav) + n] = (
                result.final_state.full().flatten()
            )
    return full_op, final_states


# ===================================================================
# Metric computation
# ===================================================================
def _logical_restricted(full_op, model, logical_n: int) -> np.ndarray:
    indices = []
    for n in range(logical_n):
        indices.extend([n, int(model.n_cav) + n])
    return full_op[np.ix_(indices, indices)]


def compute_singlebranch_metrics(full_op, final_states, model,
                                 logical_n: int, n0: int,
                                 theta: float, phi: float) -> dict:
    """Metrics for single-branch SQR target.

    Branch n0 should implement R(theta, phi); spectators should be identity
    (SQR) or any Z rotation (controlPhaseSQR).
    """
    blocks = extract_branch_unitaries(final_states, model, logical_n)
    R = target_qubit_unitary(theta, phi)
    I2 = np.eye(2, dtype=np.complex128)

    target_fid, alpha_opt = z_corrected_target_fidelity(blocks[n0], R)

    branch_true = np.zeros(logical_n)
    branch_cphase = np.zeros(logical_n)
    for n, block in enumerate(blocks):
        if n == n0:
            branch_true[n] = target_fid
            branch_cphase[n] = target_fid
        else:
            branch_true[n] = identity_fidelity_with_z(block, alpha_opt)
            sf, _ = spectator_z_fidelity(block)
            branch_cphase[n] = sf

    restricted = _logical_restricted(full_op, model, logical_n)
    tgt_blocks = [R if n == n0 else I2 for n in range(logical_n)]
    target_op = block_diag(*tgt_blocks)
    dim = float(target_op.shape[0])

    strict_fid = (abs(np.trace(target_op.conj().T @ restricted)) ** 2
                  / (dim * dim))

    overlaps = []
    for bn in range(logical_n):
        blk = restricted[2 * bn:2 * bn + 2, 2 * bn:2 * bn + 2]
        ideal = R if bn == n0 else I2
        overlaps.append(np.trace(ideal.conj().T @ blk))
    block_fid = (sum(abs(o) for o in overlaps)) ** 2 / (dim * dim)

    leakage = extract_leakage(final_states, model, logical_n)

    return {
        "target_branch_fid":      float(target_fid),
        "branch_true_mean":       float(np.mean(branch_true)),
        "branch_cphase_mean":     float(np.mean(branch_cphase)),
        "strict_logical_fid":     float(np.clip(strict_fid, 0.0, 1.0)),
        "block_phase_relaxed_fid": float(np.clip(block_fid, 0.0, 1.0)),
        "leakage_mean":           float(np.mean(leakage)),
        "leakage_max":            float(np.max(leakage)),
    }


def compute_allbranch_metrics(full_op, final_states, model,
                              logical_n: int,
                              theta: float, phi: float) -> dict:
    """Metrics for all-branch simultaneous target.

    Every branch should implement R(theta, phi).
    block_phase_relaxed_fid: per-branch Z correction allowed.
    strict_logical_fid:      single global Z only.
    """
    blocks = extract_branch_unitaries(final_states, model, logical_n)
    R = target_qubit_unitary(theta, phi)

    branch_fids = np.array([
        z_corrected_target_fidelity(b, R)[0] for b in blocks
    ])

    # Block-phase-relaxed: per-branch Z freedom (geometric mean of sqrt)
    block_fid = (np.sum(np.sqrt(branch_fids)) / logical_n) ** 2

    # Strict logical fidelity
    restricted = _logical_restricted(full_op, model, logical_n)
    target_op = block_diag(*[R] * logical_n)
    dim = float(target_op.shape[0])
    strict_fid = (abs(np.trace(target_op.conj().T @ restricted)) ** 2
                  / (dim * dim))

    leakage = extract_leakage(final_states, model, logical_n)

    return {
        "branch_fids":            branch_fids.copy(),
        "branch_fid_mean":        float(np.mean(branch_fids)),
        "branch_fid_min":         float(np.min(branch_fids)),
        "block_phase_relaxed_fid": float(np.clip(block_fid, 0.0, 1.0)),
        "strict_logical_fid":     float(np.clip(strict_fid, 0.0, 1.0)),
        "leakage_mean":           float(np.mean(leakage)),
        "leakage_max":            float(np.max(leakage)),
    }


# ===================================================================
# Optimiser: free corrections (d_lambda, d_alpha, delta_omega) per tone
# ===================================================================
def optimize_corrections(
    model,
    frame,
    T: float,
    logical_n: int,
    n0: int,
    theta: float,
    phi: float,
    all_branch: bool = False,
    objective: str = "block_phase_relaxed_fid",
    use_de: bool = True,
    maxiter_nm: int = OPT_MAXITER_NM,
    maxiter_de: int = OPT_MAXITER_DE,
) -> dict:
    """Optimise three corrections per tone: (d_lambda_n, d_alpha_n, delta_omega_n).

    Baseline (zero corrections):
      - single-branch: tone n0 active at theta, spectators silent (theta_n=0)
      - all-branch:    all tones active at theta (standard equal-amplitude multitone)

    The optimizer finds the corrections that maximise the chosen objective.

    Parameters
    ----------
    n0         : target branch for single-branch case (ignored if all_branch)
    all_branch : if True, all tones have theta_n=theta as nominal
    objective  : key in the metric dict to maximise
    use_de     : warm-start Nelder-Mead with differential evolution
    """
    lam0 = lambda_0(T)
    n_tones = logical_n

    # Nominal parameters (no corrections)
    tone_thetas = np.zeros(n_tones)
    tone_phis = np.zeros(n_tones)
    tone_det_freqs = get_branch_carrier_freqs(model, frame, n_tones)

    if all_branch:
        tone_thetas[:] = theta
        tone_phis[:] = phi
    else:
        tone_thetas[n0] = theta
        tone_phis[n0] = phi

    # Choose metric computation
    if all_branch:
        def compute_metrics(fo, fs):
            return compute_allbranch_metrics(fo, fs, model, n_tones, theta, phi)
    else:
        def compute_metrics(fo, fs):
            return compute_singlebranch_metrics(fo, fs, model, n_tones, n0, theta, phi)

    # Parameter layout: x = [d_lambda_scaled (n), d_alpha (n), d_omega_scaled (n)]
    # Scaling: d_lambda = x_lambda * lam0 * 0.2
    #          delta_omega = x_omega * |chi| * 0.2
    lam_scale = lam0 * 0.2
    chi_scale = abs(CHI) * 0.2
    n_params = 3 * n_tones
    x0 = np.zeros(n_params)

    def decode(x: np.ndarray):
        d_lambdas = x[:n_tones] * lam_scale
        d_alphas = x[n_tones:2 * n_tones]
        delta_omegas = x[2 * n_tones:] * chi_scale
        return d_lambdas, d_alphas, delta_omegas

    def cost(x: np.ndarray) -> float:
        d_lam, d_alp, d_omg = decode(x)
        try:
            pulses, dops, T_tot = build_corrected_multitone(
                model, frame, T,
                tone_thetas, tone_phis, tone_det_freqs,
                d_lam, d_alp, d_omg,
            )
            fo, fs = simulate_and_extract(model, frame, pulses, dops, n_tones, T_tot)
            m = compute_metrics(fo, fs)
            return -float(m[objective])
        except Exception:
            return 0.0  # worst-case fidelity

    if use_de:
        bounds_de = [(-6.0, 6.0)] * n_params
        de_res = differential_evolution(
            cost, bounds_de, maxiter=maxiter_de, seed=42,
            tol=1e-6, polish=False, x0=x0,
            mutation=(0.5, 1.0), recombination=0.7,
        )
        x0 = de_res.x

    nm_res = minimize(
        cost, x0, method="Nelder-Mead",
        options={"maxiter": maxiter_nm, "xatol": 1e-5,
                 "fatol": 1e-6, "adaptive": True},
    )

    d_lam, d_alp, d_omg = decode(nm_res.x)
    return {
        "d_lambdas":     d_lam,
        "d_alphas":      d_alp,
        "delta_omegas":  d_omg,
        "s_values":      tone_thetas / np.pi + d_lam / lam0,
        "opt_fidelity":  -nm_res.fun,
        "opt_nfev":      nm_res.nfev,
        "opt_success":   nm_res.success,
    }


# ===================================================================
# GRAPE comparison
# ===================================================================
def run_grape(model, frame, logical_n: int, n0: int,
              theta: float, phi: float, T: float,
              cphase: bool, all_branch: bool = False) -> dict:
    """GRAPE upper-bound for SQR or cPhaseSQR targets."""
    from cqed_sim import (
        GrapeConfig, GrapeSolver,
        ModelControlChannelSpec, PiecewiseConstantTimeGrid,
        UnitaryObjective, build_control_problem_from_model,
    )
    from cqed_sim.unitary_synthesis import Subspace

    n_slices = 48
    amp_bound = 2 * np.pi * 50e6

    indices, labels = [], []
    for n in range(logical_n):
        indices += [0 * int(model.n_cav) + n, 1 * int(model.n_cav) + n]
        labels += [f"|g,{n}>", f"|e,{n}>"]
    sub = Subspace(full_dim=int(model.n_tr) * int(model.n_cav),
                   indices=tuple(indices), labels=tuple(labels))

    R = target_qubit_unitary(theta, phi)
    I2 = np.eye(2, dtype=np.complex128)
    if all_branch:
        tgt_blocks = [R] * logical_n
    else:
        tgt_blocks = [R if n == n0 else I2 for n in range(logical_n)]
    target_sub = block_diag(*tgt_blocks)

    obj_kwargs: dict = dict(
        target_operator=target_sub, subspace=sub, ignore_global_phase=True,
    )
    if cphase:
        obj_kwargs["phase_blocks"] = tuple(
            (2 * n, 2 * n + 1) for n in range(logical_n)
        )

    problem = build_control_problem_from_model(
        model, frame=frame,
        time_grid=PiecewiseConstantTimeGrid.uniform(
            steps=n_slices, dt_s=T / n_slices,
        ),
        channel_specs=(
            ModelControlChannelSpec("qubit_I", target="qubit",
                                    quadratures=("I",),
                                    amplitude_bounds=(-amp_bound, amp_bound)),
            ModelControlChannelSpec("qubit_Q", target="qubit",
                                    quadratures=("Q",),
                                    amplitude_bounds=(-amp_bound, amp_bound)),
        ),
        objectives=(UnitaryObjective(**obj_kwargs),),
    )

    result = GrapeSolver(GrapeConfig(maxiter=300, seed=42)).solve(problem)
    if "nominal_fidelity" in result.metrics:
        fid = float(result.metrics["nominal_fidelity"])
    elif "fidelity" in result.metrics:
        fid = float(result.metrics["fidelity"])
    else:
        fid = 1.0 - float(result.objective_value)
    return {"fidelity": fid, "converged": bool(result.success)}


# ===================================================================
# Main study runner
# ===================================================================
def run_study() -> None:
    STUDY_DIR.joinpath("data").mkdir(exist_ok=True)

    model = build_model()
    frame = build_frame(model)
    nC = len(CHI_T_VALUES)

    # ------------------------------------------------------------------
    # Storage
    # ------------------------------------------------------------------
    results: dict = {"chi_t_values": CHI_T_VALUES, "logical_n": LOGICAL_N}

    def _init_arrays(prefix: str, extra_keys: list[str]) -> None:
        base_keys = ["strict_logical_fid", "block_phase_relaxed_fid", "leakage_max"]
        for k in base_keys + extra_keys:
            results[f"{prefix}_{k}"] = np.full(nC, np.nan)

    # Single-branch arrays
    for obj in ["sqr", "cphase"]:
        for var in ["baseline", "opt"]:
            _init_arrays(f"sb_{obj}_{var}",
                         ["target_branch_fid", "branch_cphase_mean",
                          "branch_true_mean"])

    # All-branch arrays
    for obj in ["sqr", "cphase"]:
        for var in ["baseline", "opt"]:
            _init_arrays(f"ab_{obj}_{var}",
                         ["branch_fid_mean", "branch_fid_min"])

    # GRAPE arrays (subset of chi_t for speed)
    grape_chi_t = np.array([1.0, 2.0, 3.0, 5.0])
    results["grape_chi_t"] = grape_chi_t
    nG = len(grape_chi_t)
    for tag in ["sb_sqr", "sb_cphase", "ab_sqr", "ab_cphase"]:
        results[f"grape_{tag}_fid"] = np.full(nG, np.nan)

    # Correction arrays (for the optimised runs)
    for case in ["sb", "ab"]:
        for obj in ["sqr", "cphase"]:
            results[f"{case}_{obj}_opt_d_lambdas"] = np.full((nC, LOGICAL_N), np.nan)
            results[f"{case}_{obj}_opt_d_alphas"] = np.full((nC, LOGICAL_N), np.nan)
            results[f"{case}_{obj}_opt_delta_omegas"] = np.full((nC, LOGICAL_N), np.nan)
            results[f"{case}_{obj}_opt_s_values"] = np.full((nC, LOGICAL_N), np.nan)

    t_study = time.time()

    # ------------------------------------------------------------------
    # SINGLE-BRANCH cases
    # ------------------------------------------------------------------
    for obj_name, obj_key in [("sqr", "strict_logical_fid"),
                               ("cphase", "block_phase_relaxed_fid")]:
        print(f"\n{'='*64}")
        print(f"SINGLE-BRANCH  objective = {obj_name.upper()} ({obj_key})")
        print(f"{'='*64}")

        tone_thetas = np.zeros(LOGICAL_N)
        tone_thetas[TARGET_N0] = REP_THETA
        tone_phis = np.zeros(LOGICAL_N)
        tone_phis[TARGET_N0] = REP_PHI
        d0 = np.zeros(LOGICAL_N)

        for ci, chi_t in enumerate(CHI_T_VALUES):
            T = duration_from_chi_t(chi_t)
            tone_det_freqs = get_branch_carrier_freqs(model, frame, LOGICAL_N)

            # --- Baseline (zero corrections) ---
            t0 = time.time()
            pulses, dops, T_tot = build_corrected_multitone(
                model, frame, T,
                tone_thetas, tone_phis, tone_det_freqs, d0, d0, d0,
            )
            fo, fs = simulate_and_extract(model, frame, pulses, dops, LOGICAL_N, T_tot)
            mb = compute_singlebranch_metrics(
                fo, fs, model, LOGICAL_N, TARGET_N0, REP_THETA, REP_PHI)
            pref = f"sb_{obj_name}_baseline"
            for k in ["strict_logical_fid", "block_phase_relaxed_fid", "leakage_max",
                      "target_branch_fid", "branch_cphase_mean", "branch_true_mean"]:
                results[f"{pref}_{k}"][ci] = mb[k]
            print(f"  chiT/2pi={chi_t:.1f}  baseline ({time.time()-t0:.1f}s):  "
                  f"F_strict={mb['strict_logical_fid']:.4f}  "
                  f"F_block={mb['block_phase_relaxed_fid']:.4f}")
            gc.collect()

            # --- Optimised corrections ---
            t0 = time.time()
            opt = optimize_corrections(
                model, frame, T, LOGICAL_N, TARGET_N0, REP_THETA, REP_PHI,
                all_branch=False, objective=obj_key, use_de=True,
            )
            pulses, dops, T_tot = build_corrected_multitone(
                model, frame, T,
                tone_thetas, tone_phis, tone_det_freqs,
                opt["d_lambdas"], opt["d_alphas"], opt["delta_omegas"],
            )
            fo, fs = simulate_and_extract(model, frame, pulses, dops, LOGICAL_N, T_tot)
            mo = compute_singlebranch_metrics(
                fo, fs, model, LOGICAL_N, TARGET_N0, REP_THETA, REP_PHI)
            pref = f"sb_{obj_name}_opt"
            for k in ["strict_logical_fid", "block_phase_relaxed_fid", "leakage_max",
                      "target_branch_fid", "branch_cphase_mean", "branch_true_mean"]:
                results[f"{pref}_{k}"][ci] = mo[k]
            results[f"sb_{obj_name}_opt_d_lambdas"][ci] = opt["d_lambdas"]
            results[f"sb_{obj_name}_opt_d_alphas"][ci] = opt["d_alphas"]
            results[f"sb_{obj_name}_opt_delta_omegas"][ci] = opt["delta_omegas"]
            results[f"sb_{obj_name}_opt_s_values"][ci] = opt["s_values"]
            print(f"             optimised ({time.time()-t0:.1f}s):  "
                  f"F_strict={mo['strict_logical_fid']:.4f}  "
                  f"F_block={mo['block_phase_relaxed_fid']:.4f}  "
                  f"(nfev={opt['opt_nfev']})")
            gc.collect()

        # checkpoint
        np.savez_compressed(str(SAVE_PATH), **results)
        print(f"[checkpoint saved, elapsed {time.time()-t_study:.0f}s]")

    # ------------------------------------------------------------------
    # ALL-BRANCH cases
    # ------------------------------------------------------------------
    tone_thetas_ab = np.full(LOGICAL_N, REP_THETA)
    tone_phis_ab = np.full(LOGICAL_N, REP_PHI)
    d0 = np.zeros(LOGICAL_N)

    for obj_name, obj_key in [("sqr", "strict_logical_fid"),
                               ("cphase", "block_phase_relaxed_fid")]:
        print(f"\n{'='*64}")
        print(f"ALL-BRANCH     objective = {obj_name.upper()} ({obj_key})")
        print(f"{'='*64}")

        for ci, chi_t in enumerate(CHI_T_VALUES):
            T = duration_from_chi_t(chi_t)
            tone_det_freqs = get_branch_carrier_freqs(model, frame, LOGICAL_N)

            # --- Baseline ---
            t0 = time.time()
            pulses, dops, T_tot = build_corrected_multitone(
                model, frame, T,
                tone_thetas_ab, tone_phis_ab, tone_det_freqs, d0, d0, d0,
            )
            fo, fs = simulate_and_extract(model, frame, pulses, dops, LOGICAL_N, T_tot)
            mb = compute_allbranch_metrics(fo, fs, model, LOGICAL_N, REP_THETA, REP_PHI)
            pref = f"ab_{obj_name}_baseline"
            for k in ["strict_logical_fid", "block_phase_relaxed_fid", "leakage_max",
                      "branch_fid_mean", "branch_fid_min"]:
                results[f"{pref}_{k}"][ci] = mb[k]
            print(f"  chiT/2pi={chi_t:.1f}  baseline ({time.time()-t0:.1f}s):  "
                  f"F_strict={mb['strict_logical_fid']:.4f}  "
                  f"F_block={mb['block_phase_relaxed_fid']:.4f}")
            gc.collect()

            # --- Optimised corrections ---
            t0 = time.time()
            opt = optimize_corrections(
                model, frame, T, LOGICAL_N, TARGET_N0, REP_THETA, REP_PHI,
                all_branch=True, objective=obj_key, use_de=True,
            )
            pulses, dops, T_tot = build_corrected_multitone(
                model, frame, T,
                tone_thetas_ab, tone_phis_ab, tone_det_freqs,
                opt["d_lambdas"], opt["d_alphas"], opt["delta_omegas"],
            )
            fo, fs = simulate_and_extract(model, frame, pulses, dops, LOGICAL_N, T_tot)
            mo = compute_allbranch_metrics(fo, fs, model, LOGICAL_N, REP_THETA, REP_PHI)
            pref = f"ab_{obj_name}_opt"
            for k in ["strict_logical_fid", "block_phase_relaxed_fid", "leakage_max",
                      "branch_fid_mean", "branch_fid_min"]:
                results[f"{pref}_{k}"][ci] = mo[k]
            results[f"ab_{obj_name}_opt_d_lambdas"][ci] = opt["d_lambdas"]
            results[f"ab_{obj_name}_opt_d_alphas"][ci] = opt["d_alphas"]
            results[f"ab_{obj_name}_opt_delta_omegas"][ci] = opt["delta_omegas"]
            results[f"ab_{obj_name}_opt_s_values"][ci] = opt["s_values"]
            print(f"             optimised ({time.time()-t0:.1f}s):  "
                  f"F_strict={mo['strict_logical_fid']:.4f}  "
                  f"F_block={mo['block_phase_relaxed_fid']:.4f}  "
                  f"(nfev={opt['opt_nfev']})")
            gc.collect()

        np.savez_compressed(str(SAVE_PATH), **results)
        print(f"[checkpoint saved, elapsed {time.time()-t_study:.0f}s]")

    # ------------------------------------------------------------------
    # GRAPE comparison (subset of chi_t values)
    # ------------------------------------------------------------------
    print(f"\n{'='*64}")
    print("GRAPE comparison")
    print(f"{'='*64}")

    for gi, chi_t in enumerate(grape_chi_t):
        T = duration_from_chi_t(chi_t)
        print(f"  chiT/2pi = {chi_t:.1f}  T = {T*1e9:.1f} ns", flush=True)

        for tag, kwargs in [
            ("sb_sqr",    dict(all_branch=False, cphase=False)),
            ("sb_cphase", dict(all_branch=False, cphase=True)),
            ("ab_sqr",    dict(all_branch=True,  cphase=False)),
            ("ab_cphase", dict(all_branch=True,  cphase=True)),
        ]:
            t0 = time.time()
            try:
                gr = run_grape(model, frame, LOGICAL_N, TARGET_N0,
                               REP_THETA, REP_PHI, T, **kwargs)
                results[f"grape_{tag}_fid"][gi] = gr["fidelity"]
                print(f"    {tag:12s}  F={gr['fidelity']:.6f}  "
                      f"conv={gr['converged']}  ({time.time()-t0:.1f}s)")
            except Exception as exc:
                print(f"    {tag:12s}  GRAPE failed: {exc}")
            gc.collect()

    # ------------------------------------------------------------------
    # Save final results
    # ------------------------------------------------------------------
    np.savez_compressed(str(SAVE_PATH), **results)
    elapsed = time.time() - t_study
    print(f"\nAll results saved to {SAVE_PATH}")
    print(f"Total elapsed: {elapsed:.0f} s  ({elapsed/60:.1f} min)")

    # ------------------------------------------------------------------
    # Console summary
    # ------------------------------------------------------------------
    print("\n" + "=" * 64)
    print("SUMMARY — block_phase_relaxed_fid")
    print("=" * 64)
    hdr = f"{'chiT/2pi':>6}  {'SB-base':>8}  {'SB-sqr':>7}  {'SB-cph':>7}  " \
          f"{'AB-base':>8}  {'AB-sqr':>7}  {'AB-cph':>7}"
    print(hdr)
    print("-" * len(hdr))
    for ci, ct in enumerate(CHI_T_VALUES):
        vals = [
            results["sb_sqr_baseline_block_phase_relaxed_fid"][ci],
            results["sb_sqr_opt_block_phase_relaxed_fid"][ci],
            results["sb_cphase_opt_block_phase_relaxed_fid"][ci],
            results["ab_sqr_baseline_block_phase_relaxed_fid"][ci],
            results["ab_sqr_opt_block_phase_relaxed_fid"][ci],
            results["ab_cphase_opt_block_phase_relaxed_fid"][ci],
        ]
        print(f"{ct:6.1f}  " + "  ".join(f"{v:8.4f}" for v in vals))


if __name__ == "__main__":
    run_study()
