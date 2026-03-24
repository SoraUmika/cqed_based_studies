"""
All-branch multitone SQR study: best-optimized Gaussian multitone vs GRAPE.

Targets simultaneous R_X(π) rotations on ALL Fock branches n = 0, 1, 2, 3.
The multitone Gaussian waveform has per-tone (amplitude, phase, detuning)
parameters optimised via Nelder-Mead.  GRAPE provides the upper-bound comparison.

Scan variable: χT/2π ∈ {0.5, 1.0, 1.5, 2.0, 3.0, 5.0}

Usage:
    python scripts/run_allbranch_multitone.py

Output:
    data/allbranch_multitone_results.npz
"""
from __future__ import annotations

import gc
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.optimize import minimize, differential_evolution
from scipy.linalg import block_diag

SCRIPT_DIR = Path(__file__).resolve().parent
STUDY_DIR = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

import runtime_compat  # noqa: F401  (WMI patch for Windows)

from common import (
    CHI, N_CAV, N_TR, N_FOCK, DT,
    build_frame, build_model, duration_from_chi_t,
    target_qubit_unitary, z_corrected_target_fidelity,
    extract_branch_unitaries, extract_leakage,
)
from cqed_sim.core.frequencies import (
    carrier_for_transition_frequency,
    manifold_transition_frequency,
)
from cqed_sim.sequence import SequenceCompiler
from cqed_sim.sim import SimulationConfig, prepare_simulation

# ===================================================================
# Constants
# ===================================================================
CHI_T_VALUES = np.array([0.5, 1.0, 1.5, 2.0, 3.0, 5.0], dtype=float)

# All-branch target: R_X(π) on EVERY Fock branch n=0..3
THETA = np.pi
PHI = 0.0
LOGICAL_N = N_FOCK  # = 4

SIGMA_FRACTION = 1.0 / 6.0

# Optimization
OPT_MAXITER_NM = 300      # Nelder-Mead
OPT_MAXITER_DE = 100      # differential evolution (for global search)

# GRAPE
N_SLICES = 48
GRAPE_MAXITER = 300
AMP_BOUND = 2 * np.pi * 50e6  # rad/s

DATA_DIR = STUDY_DIR / "data"
DATA_DIR.mkdir(exist_ok=True)
SAVE_PATH = DATA_DIR / "allbranch_multitone_results.npz"


# ===================================================================
# Local pulse dataclass (matches existing study pattern)
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


# ===================================================================
# Simulation helper (follows existing simulate_and_extract pattern)
# ===================================================================
def simulate_and_extract(model, frame, pulses, drive_ops, T_total):
    """Simulate all logical basis inputs and return final states.

    Returns list of qutip Qobj final states, ordered as:
    [|g,0⟩_final, |e,0⟩_final, |g,1⟩_final, |e,1⟩_final, ...]
    """
    compiler = SequenceCompiler(dt=DT)
    compiled = compiler.compile(pulses, t_end=float(T_total + 4.0 * DT))
    config = SimulationConfig(frame=frame, store_states=False)
    session = prepare_simulation(model, compiled, drive_ops, config=config, e_ops={})

    final_states = []
    for n in range(LOGICAL_N):
        for q in (0, 1):
            initial = model.basis_state(q, n)
            result = session.run(initial)
            final_states.append(result.final_state)
    return final_states


# ===================================================================
# Build all-branch target unitary
# ===================================================================
def build_allbranch_target():
    """8×8 block-diagonal target: R_X(π) on every Fock branch."""
    R = target_qubit_unitary(THETA, PHI)
    blocks = [R for _ in range(LOGICAL_N)]
    return block_diag(*blocks)


# ===================================================================
# Fidelity computation — block-phase-relaxed
# ===================================================================
def block_fidelity(U_branches, target_branches):
    """Block-phase-relaxed fidelity: per-branch Z-correction allowed.

    U_branches: list of 2×2 arrays (actual branch unitaries)
    target_branches: list of 2×2 target matrices (one per branch)
    """
    total = 0.0
    for Ub, Ut in zip(U_branches, target_branches):
        fid, _ = z_corrected_target_fidelity(Ub, Ut)
        total += np.sqrt(fid)
    return (total / len(U_branches)) ** 2


def strict_fidelity(U_branches, target_branches):
    """Strict logical fidelity with single global Z-correction."""
    d = 2 * len(U_branches)
    U_full = block_diag(*U_branches)
    T_full = block_diag(*target_branches)
    best_fid = 0.0
    for alpha_trial in np.linspace(0, 2 * np.pi, 361):
        Z_global = np.diag(np.tile([1.0, np.exp(1j * alpha_trial)], len(U_branches)))
        fid = abs(np.trace(T_full.conj().T @ Z_global @ U_full)) ** 2 / d ** 2
        if fid > best_fid:
            best_fid = fid
    return best_fid


# ===================================================================
# Multitone pulse builder for all-branch target
# ===================================================================
def build_allbranch_multitone_pulse(model, frame, T,
                                     amp_ratios, phases, detunings=None):
    """Build a multitone Gaussian pulse targeting all branches simultaneously.

    amp_ratios: array (LOGICAL_N,) — amplitude scaling per tone
    phases: array (LOGICAL_N,) — phase per tone (rad)
    detunings: array (LOGICAL_N,) or None — frequency offset per tone (rad/s)

    Each tone is centered at branch-n transition frequency + detuning.
    All tones share a Gaussian window (σ=T/6).
    """
    base_amp = THETA / (2.0 * T)

    tone_specs = []
    for n in range(LOGICAL_N):
        omega_n = manifold_transition_frequency(model, n, frame)
        carrier_n = carrier_for_transition_frequency(omega_n)
        if detunings is not None:
            carrier_n += detunings[n]
        tone_specs.append(MultitoneTone(
            manifold=n,
            omega_rad_s=carrier_n,
            amp_rad_s=base_amp * amp_ratios[n],
            phase_rad=phases[n],
        ))

    dur = float(T)
    ts = list(tone_specs)

    def envelope(t_rel):
        return multitone_modulated_envelope(
            np.asarray(t_rel, dtype=float), dur, ts)

    pulse = Pulse(
        channel="q", t0=0.0, duration=T,
        envelope=envelope, carrier=0.0, phase=0.0, amp=1.0,
        label="allbranch_multitone",
    )
    return [pulse], {"q": "qubit"}, T


# ===================================================================
# Baseline: common-envelope multitone (all branches, equal amp)
# ===================================================================
def build_baseline_multitone(model, frame, T):
    """Baseline multitone: all tones have amp_ratio=1, phase=0, no detuning."""
    amp_ratios = np.ones(LOGICAL_N)
    phases = np.zeros(LOGICAL_N)
    return build_allbranch_multitone_pulse(model, frame, T, amp_ratios, phases)


# ===================================================================
# Optimization: Nelder-Mead + optional differential evolution
# ===================================================================
def optimize_allbranch_multitone(model, frame, T, use_detuning=True,
                                  use_de=False, verbose=True):
    """Optimize per-tone (amp, phase, detuning) to maximize block fidelity
    for the all-branch R_X(π) target."""
    R_target = target_qubit_unitary(THETA, PHI)
    target_branches = [R_target] * LOGICAL_N

    n_params = LOGICAL_N * 2  # amps + phases
    if use_detuning:
        n_params += LOGICAL_N

    eval_count = [0]

    def cost(params):
        eval_count[0] += 1
        amps = np.clip(params[:LOGICAL_N], 0.0, 3.0)
        phs = params[LOGICAL_N:2*LOGICAL_N]
        det = params[2*LOGICAL_N:] if use_detuning else None

        pulses, cmap, Ttot = build_allbranch_multitone_pulse(
            model, frame, T, amps, phs, det)
        try:
            final_states = simulate_and_extract(model, frame, pulses, cmap, Ttot)
            branches = extract_branch_unitaries(final_states, model, LOGICAL_N)
            fid = block_fidelity(branches, target_branches)
        except Exception:
            return 1.0
        return 1.0 - fid

    x0 = np.zeros(n_params)
    x0[:LOGICAL_N] = 1.0

    if use_de:
        bounds_de = (
            [(0.0, 3.0)] * LOGICAL_N +
            [(-np.pi, np.pi)] * LOGICAL_N
        )
        if use_detuning:
            max_det = abs(CHI) * 0.5
            bounds_de += [(-max_det, max_det)] * LOGICAL_N

        if verbose:
            print(f"    DE global search ({n_params} params)...", flush=True)

        de_result = differential_evolution(
            cost, bounds_de, maxiter=OPT_MAXITER_DE,
            seed=42, tol=1e-6, polish=False,
            x0=x0, mutation=(0.5, 1.0), recombination=0.7,
        )
        x0 = de_result.x
        if verbose:
            print(f"    DE done: F_block={1-de_result.fun:.6f} "
                  f"({eval_count[0]} evals)", flush=True)
        eval_count[0] = 0

    if verbose:
        print(f"    Nelder-Mead ({n_params} params)...", flush=True)

    result = minimize(cost, x0, method="Nelder-Mead",
                      options={"maxiter": OPT_MAXITER_NM, "adaptive": True})

    opt_amps = np.clip(result.x[:LOGICAL_N], 0.0, 3.0)
    opt_phases = result.x[LOGICAL_N:2*LOGICAL_N]
    opt_det = result.x[2*LOGICAL_N:] if use_detuning else None
    fid_block = 1.0 - result.fun

    if verbose:
        print(f"    NM done: F_block={fid_block:.6f} "
              f"({eval_count[0]} evals)", flush=True)

    return {
        "fid_block": fid_block,
        "amp_ratios": opt_amps,
        "phases": opt_phases,
        "detunings": opt_det,
        "nfev": result.nfev,
        "success": result.success,
    }


# ===================================================================
# Detailed metrics extraction
# ===================================================================
def compute_detailed_metrics(model, frame, pulses, drive_ops, T):
    """Extract full set of metrics for a pulse configuration."""
    R_target = target_qubit_unitary(THETA, PHI)
    target_branches = [R_target] * LOGICAL_N

    final_states = simulate_and_extract(model, frame, pulses, drive_ops, T)
    branches = extract_branch_unitaries(final_states, model, LOGICAL_N)
    leakages = extract_leakage(final_states, model, LOGICAL_N)

    fid_block = block_fidelity(branches, target_branches)
    fid_strict = strict_fidelity(branches, target_branches)

    branch_fids = []
    for Ub, Ut in zip(branches, target_branches):
        f, _ = z_corrected_target_fidelity(Ub, Ut)
        branch_fids.append(f)

    return {
        "fid_block": fid_block,
        "fid_strict": fid_strict,
        "branch_fids": np.array(branch_fids),
        "leakage_max": float(np.max(leakages)),
        "leakages": leakages,
    }


# ===================================================================
# GRAPE comparison
# ===================================================================
def run_grape_allbranch(model, frame, T, verbose=True):
    """Run GRAPE for the all-branch R_X(π) target with per-branch phase freedom."""
    from cqed_sim import (
        GrapeConfig,
        GrapeSolver,
        ModelControlChannelSpec,
        PiecewiseConstantTimeGrid,
        UnitaryObjective,
        build_control_problem_from_model,
    )
    from cqed_sim.unitary_synthesis import Subspace

    indices = []
    labels = []
    for n in range(LOGICAL_N):
        indices.append(0 * N_CAV + n)
        indices.append(1 * N_CAV + n)
        labels.append(f"|g,{n}>")
        labels.append(f"|e,{n}>")
    sub = Subspace(
        full_dim=N_TR * N_CAV,
        indices=tuple(indices),
        labels=tuple(labels),
    )

    target_sub = build_allbranch_target()
    phase_blocks = tuple((2 * n, 2 * n + 1) for n in range(LOGICAL_N))
    dt_grape = T / N_SLICES

    problem = build_control_problem_from_model(
        model,
        frame=frame,
        time_grid=PiecewiseConstantTimeGrid.uniform(
            steps=N_SLICES, dt_s=dt_grape),
        channel_specs=(
            ModelControlChannelSpec(
                name="qubit_I", target="qubit",
                quadratures=("I",),
                amplitude_bounds=(-AMP_BOUND, AMP_BOUND),
            ),
            ModelControlChannelSpec(
                name="qubit_Q", target="qubit",
                quadratures=("Q",),
                amplitude_bounds=(-AMP_BOUND, AMP_BOUND),
            ),
        ),
        objectives=(
            UnitaryObjective(
                target_operator=target_sub,
                subspace=sub,
                ignore_global_phase=True,
                phase_blocks=phase_blocks,
            ),
        ),
    )

    config = GrapeConfig(maxiter=GRAPE_MAXITER, seed=42)
    if verbose:
        print(f"    GRAPE solving ({N_SLICES} slices, "
              f"{GRAPE_MAXITER} iter)...", flush=True)
    t0 = time.time()
    result = GrapeSolver(config).solve(problem)
    elapsed = time.time() - t0

    if "nominal_fidelity" in result.metrics:
        fid = result.metrics["nominal_fidelity"]
    elif "fidelity" in result.metrics:
        fid = result.metrics["fidelity"]
    else:
        fid = (1.0 - result.objective_value
               if result.objective_value <= 1.0 else 0.0)

    if verbose:
        print(f"    GRAPE done: F={fid:.8f}, converged={result.success}, "
              f"time={elapsed:.1f}s", flush=True)

    return {
        "fidelity": fid,
        "objective_value": float(result.objective_value),
        "converged": bool(result.success),
        "elapsed_s": elapsed,
    }


# ===================================================================
# MAIN
# ===================================================================
def main():
    print("=" * 70)
    print("All-Branch Multitone SQR Study")
    print("=" * 70)

    model = build_model()
    frame = build_frame(model)

    nC = len(CHI_T_VALUES)

    # Storage arrays
    baseline_fid_block = np.full(nC, np.nan)
    baseline_fid_strict = np.full(nC, np.nan)
    baseline_branch_fids = np.full((nC, LOGICAL_N), np.nan)
    baseline_leakage = np.full(nC, np.nan)

    opt_ap_fid_block = np.full(nC, np.nan)
    opt_ap_fid_strict = np.full(nC, np.nan)
    opt_ap_branch_fids = np.full((nC, LOGICAL_N), np.nan)
    opt_ap_leakage = np.full(nC, np.nan)

    opt_apd_fid_block = np.full(nC, np.nan)
    opt_apd_fid_strict = np.full(nC, np.nan)
    opt_apd_branch_fids = np.full((nC, LOGICAL_N), np.nan)
    opt_apd_leakage = np.full(nC, np.nan)

    opt_de_fid_block = np.full(nC, np.nan)
    opt_de_fid_strict = np.full(nC, np.nan)
    opt_de_branch_fids = np.full((nC, LOGICAL_N), np.nan)
    opt_de_leakage = np.full(nC, np.nan)

    grape_fid = np.full(nC, np.nan)
    grape_converged = np.full(nC, False)

    def save_checkpoint():
        np.savez_compressed(
            str(SAVE_PATH),
            chi_t_values=CHI_T_VALUES,
            baseline_fid_block=baseline_fid_block,
            baseline_fid_strict=baseline_fid_strict,
            baseline_branch_fids=baseline_branch_fids,
            baseline_leakage=baseline_leakage,
            opt_ap_fid_block=opt_ap_fid_block,
            opt_ap_fid_strict=opt_ap_fid_strict,
            opt_ap_branch_fids=opt_ap_branch_fids,
            opt_ap_leakage=opt_ap_leakage,
            opt_apd_fid_block=opt_apd_fid_block,
            opt_apd_fid_strict=opt_apd_fid_strict,
            opt_apd_branch_fids=opt_apd_branch_fids,
            opt_apd_leakage=opt_apd_leakage,
            opt_de_fid_block=opt_de_fid_block,
            opt_de_fid_strict=opt_de_fid_strict,
            opt_de_branch_fids=opt_de_branch_fids,
            opt_de_leakage=opt_de_leakage,
            grape_fid=grape_fid,
            grape_converged=grape_converged,
        )
        print(f"  [checkpoint saved to {SAVE_PATH.name}]", flush=True)

    # ---------------------------------------------------------------
    # Phase A: Baseline (common-envelope multitone, all branches)
    # ---------------------------------------------------------------
    print("\n--- Phase A: Baseline common-envelope multitone ---")
    for i, ct in enumerate(CHI_T_VALUES):
        T = duration_from_chi_t(ct)
        print(f"  chiT/2pi = {ct}: T = {T*1e9:.1f} ns ...", end=" ", flush=True)
        t0 = time.time()

        pulses, cmap, Ttot = build_baseline_multitone(model, frame, T)
        metrics = compute_detailed_metrics(model, frame, pulses, cmap, T)

        baseline_fid_block[i] = metrics["fid_block"]
        baseline_fid_strict[i] = metrics["fid_strict"]
        baseline_branch_fids[i] = metrics["branch_fids"]
        baseline_leakage[i] = metrics["leakage_max"]

        elapsed = time.time() - t0
        print(f"F_block={metrics['fid_block']:.4f}, "
              f"F_strict={metrics['fid_strict']:.4f} ({elapsed:.1f}s)", flush=True)
        gc.collect()

    save_checkpoint()

    # ---------------------------------------------------------------
    # Phase B: Optimized amp+phase (Nelder-Mead)
    # ---------------------------------------------------------------
    print("\n--- Phase B: Optimized amp+phase (8 params, NM) ---")
    for i, ct in enumerate(CHI_T_VALUES):
        T = duration_from_chi_t(ct)
        print(f"  chiT/2pi = {ct}: T = {T*1e9:.1f} ns", flush=True)
        t0 = time.time()

        opt = optimize_allbranch_multitone(model, frame, T,
                                            use_detuning=False, use_de=False)

        pulses, cmap, Ttot = build_allbranch_multitone_pulse(
            model, frame, T, opt["amp_ratios"], opt["phases"])
        metrics = compute_detailed_metrics(model, frame, pulses, cmap, T)

        opt_ap_fid_block[i] = metrics["fid_block"]
        opt_ap_fid_strict[i] = metrics["fid_strict"]
        opt_ap_branch_fids[i] = metrics["branch_fids"]
        opt_ap_leakage[i] = metrics["leakage_max"]

        elapsed = time.time() - t0
        print(f"    -> F_block={metrics['fid_block']:.4f}, "
              f"F_strict={metrics['fid_strict']:.4f} ({elapsed:.1f}s)", flush=True)
        gc.collect()

    save_checkpoint()

    # ---------------------------------------------------------------
    # Phase C: Optimized amp+phase+detuning (Nelder-Mead)
    # ---------------------------------------------------------------
    print("\n--- Phase C: Optimized amp+phase+detuning (12 params, NM) ---")
    for i, ct in enumerate(CHI_T_VALUES):
        T = duration_from_chi_t(ct)
        print(f"  chiT/2pi = {ct}: T = {T*1e9:.1f} ns", flush=True)
        t0 = time.time()

        opt = optimize_allbranch_multitone(model, frame, T,
                                            use_detuning=True, use_de=False)

        pulses, cmap, Ttot = build_allbranch_multitone_pulse(
            model, frame, T, opt["amp_ratios"], opt["phases"],
            opt["detunings"])
        metrics = compute_detailed_metrics(model, frame, pulses, cmap, T)

        opt_apd_fid_block[i] = metrics["fid_block"]
        opt_apd_fid_strict[i] = metrics["fid_strict"]
        opt_apd_branch_fids[i] = metrics["branch_fids"]
        opt_apd_leakage[i] = metrics["leakage_max"]

        elapsed = time.time() - t0
        print(f"    -> F_block={metrics['fid_block']:.4f}, "
              f"F_strict={metrics['fid_strict']:.4f} ({elapsed:.1f}s)", flush=True)
        gc.collect()

    save_checkpoint()

    # ---------------------------------------------------------------
    # Phase D: Optimized amp+phase+detuning with DE warm-start
    # ---------------------------------------------------------------
    print("\n--- Phase D: DE warm-start + NM (12 params) ---")
    for i, ct in enumerate(CHI_T_VALUES):
        T = duration_from_chi_t(ct)
        print(f"  chiT/2pi = {ct}: T = {T*1e9:.1f} ns", flush=True)
        t0 = time.time()

        opt = optimize_allbranch_multitone(model, frame, T,
                                            use_detuning=True, use_de=True)

        pulses, cmap, Ttot = build_allbranch_multitone_pulse(
            model, frame, T, opt["amp_ratios"], opt["phases"],
            opt["detunings"])
        metrics = compute_detailed_metrics(model, frame, pulses, cmap, T)

        opt_de_fid_block[i] = metrics["fid_block"]
        opt_de_fid_strict[i] = metrics["fid_strict"]
        opt_de_branch_fids[i] = metrics["branch_fids"]
        opt_de_leakage[i] = metrics["leakage_max"]

        elapsed = time.time() - t0
        print(f"    -> F_block={metrics['fid_block']:.4f}, "
              f"F_strict={metrics['fid_strict']:.4f} ({elapsed:.1f}s)", flush=True)
        gc.collect()

    save_checkpoint()

    # ---------------------------------------------------------------
    # Phase E: GRAPE comparison
    # ---------------------------------------------------------------
    print("\n--- Phase E: GRAPE all-branch comparison ---")
    for i, ct in enumerate(CHI_T_VALUES):
        T = duration_from_chi_t(ct)
        print(f"  chiT/2pi = {ct}: T = {T*1e9:.1f} ns", flush=True)
        t0 = time.time()

        try:
            grape_res = run_grape_allbranch(model, frame, T)
            grape_fid[i] = grape_res["fidelity"]
            grape_converged[i] = grape_res["converged"]
        except Exception as e:
            print(f"    GRAPE failed: {e}", flush=True)

        elapsed = time.time() - t0
        print(f"    -> F_GRAPE={grape_fid[i]:.8f} ({elapsed:.1f}s)", flush=True)
        gc.collect()

    save_checkpoint()

    # ---------------------------------------------------------------
    # Summary
    # ---------------------------------------------------------------
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"{'chiT/2pi':>8}  {'Baseline':>10}  {'Opt(a,p)':>10}  "
          f"{'Opt(a,p,d)':>11}  {'DE+NM':>10}  {'GRAPE':>10}")
    print("-" * 70)
    for i, ct in enumerate(CHI_T_VALUES):
        print(f"{ct:8.1f}  {baseline_fid_block[i]:10.4f}  "
              f"{opt_ap_fid_block[i]:10.4f}  {opt_apd_fid_block[i]:11.4f}  "
              f"{opt_de_fid_block[i]:10.4f}  {grape_fid[i]:10.6f}")
    print("=" * 70)
    print(f"\nAll results saved to {SAVE_PATH}")


if __name__ == "__main__":
    main()
