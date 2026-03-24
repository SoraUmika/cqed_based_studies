"""
All-branch multitone SQR study — accelerated with multiprocessing.

Targets simultaneous R_X(pi) rotations on ALL Fock branches n = 0..3.
Uses multiprocessing.Pool to parallelize:
  - DE population evaluation across CPU cores (Phase D)
  - Independent chi*T points across CPU cores (Phases A/B/C/E)

Usage:
    python scripts/run_allbranch_fast.py

Output:
    data/allbranch_multitone_results.npz
"""
from __future__ import annotations

import gc
import multiprocessing as mp
import os
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

# Ensure scripts dir is on path (needed for spawned workers on Windows)
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import runtime_compat  # noqa: F401  (WMI patch)

from common import (
    CHI,
    N_CAV,
    N_FOCK,
    N_TR,
    DT,
    build_frame,
    build_model,
    duration_from_chi_t,
    extract_branch_unitaries,
    extract_leakage,
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
# Constants
# ===================================================================
CHI_T_VALUES = np.array([0.5, 1.0, 1.5, 2.0, 3.0, 5.0], dtype=float)

THETA = np.pi
PHI = 0.0
LOGICAL_N = N_FOCK  # = 4
SIGMA_FRACTION = 1.0 / 6.0

OPT_MAXITER_NM = 300
OPT_MAXITER_DE = 100

N_SLICES = 48
GRAPE_MAXITER = 300
AMP_BOUND = 2 * np.pi * 50e6  # rad/s

DATA_DIR = STUDY_DIR / "data"
SAVE_PATH = DATA_DIR / "allbranch_multitone_results.npz"

# Leave 2 cores free for system responsiveness
N_WORKERS = max(1, min(os.cpu_count() - 2, 16))

# ===================================================================
# Precomputed target
# ===================================================================
_R_TARGET = target_qubit_unitary(THETA, PHI)
_TARGET_BRANCHES = [_R_TARGET] * LOGICAL_N

# ===================================================================
# Per-process lazy model initialization (for multiprocessing on Windows)
# ===================================================================
_proc_model = None
_proc_frame = None


def _ensure_model():
    """Lazily build model + frame once per worker process."""
    global _proc_model, _proc_frame
    if _proc_model is None:
        # Re-apply WMI patch in spawned process
        if str(SCRIPT_DIR) not in sys.path:
            sys.path.insert(0, str(SCRIPT_DIR))
        import runtime_compat as _rc  # noqa: F401
        _proc_model = build_model()
        _proc_frame = build_frame(_proc_model)
    return _proc_model, _proc_frame


def _worker_init():
    """Pool initializer: build model eagerly in each worker."""
    _ensure_model()


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
    env = normalized_gaussian(t_rel)
    t_abs = t_rel * duration_s
    coeff = np.zeros_like(t_abs, dtype=np.complex128)
    for spec in tone_specs:
        coeff += spec.amp_rad_s * np.exp(
            1j * spec.phase_rad
        ) * np.exp(1j * spec.omega_rad_s * t_abs)
    return env * coeff


# ===================================================================
# Simulation helper
# ===================================================================
def simulate_and_extract(model, frame, pulses, drive_ops, T_total):
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
# Target
# ===================================================================
def build_allbranch_target():
    return block_diag(*[_R_TARGET for _ in range(LOGICAL_N)])


# ===================================================================
# Fidelity metrics
# ===================================================================
def block_fidelity(U_branches, target_branches):
    total = 0.0
    for Ub, Ut in zip(U_branches, target_branches):
        fid, _ = z_corrected_target_fidelity(Ub, Ut)
        total += np.sqrt(fid)
    return (total / len(U_branches)) ** 2


def strict_fidelity(U_branches, target_branches):
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
# Pulse builder
# ===================================================================
def build_allbranch_multitone_pulse(model, frame, T, amp_ratios, phases,
                                    detunings=None):
    base_amp = THETA / (2.0 * T)
    tone_specs = []
    for n in range(LOGICAL_N):
        omega_n = manifold_transition_frequency(model, n, frame)
        carrier_n = carrier_for_transition_frequency(omega_n)
        if detunings is not None:
            carrier_n += detunings[n]
        tone_specs.append(MultitoneTone(
            manifold=n, omega_rad_s=carrier_n,
            amp_rad_s=base_amp * amp_ratios[n], phase_rad=phases[n],
        ))

    dur = float(T)
    ts = list(tone_specs)

    def envelope(t_rel):
        return multitone_modulated_envelope(np.asarray(t_rel, dtype=float), dur, ts)

    pulse = Pulse(channel="q", t0=0.0, duration=T, envelope=envelope,
                  carrier=0.0, phase=0.0, amp=1.0, label="allbranch_multitone")
    return [pulse], {"q": "qubit"}, T


# ===================================================================
# Detailed metrics
# ===================================================================
def compute_detailed_metrics(model, frame, pulses, drive_ops, T):
    final_states = simulate_and_extract(model, frame, pulses, drive_ops, T)
    branches = extract_branch_unitaries(final_states, model, LOGICAL_N)
    leakages = extract_leakage(final_states, model, LOGICAL_N)
    fid_block = block_fidelity(branches, _TARGET_BRANCHES)
    fid_strict = strict_fidelity(branches, _TARGET_BRANCHES)
    branch_fids = []
    for Ub, Ut in zip(branches, _TARGET_BRANCHES):
        f, _ = z_corrected_target_fidelity(Ub, Ut)
        branch_fids.append(f)
    return {
        "fid_block": fid_block,
        "fid_strict": fid_strict,
        "branch_fids": np.array(branch_fids),
        "leakage_max": float(np.max(leakages)),
    }


# ===================================================================
# Module-level cost function for DE parallel evaluation
# ===================================================================
def _de_cost(params, T_s, use_detuning_flag):
    """Cost for differential_evolution — uses per-process model."""
    model, frame = _ensure_model()
    amps = np.clip(params[:LOGICAL_N], 0.0, 3.0)
    phs = params[LOGICAL_N:2 * LOGICAL_N]
    det = params[2 * LOGICAL_N:] if use_detuning_flag else None
    pulses, cmap, Ttot = build_allbranch_multitone_pulse(
        model, frame, T_s, amps, phs, det)
    try:
        final_states = simulate_and_extract(model, frame, pulses, cmap, Ttot)
        branches = extract_branch_unitaries(final_states, model, LOGICAL_N)
        fid = block_fidelity(branches, _TARGET_BRANCHES)
    except Exception:
        return 1.0
    return 1.0 - fid


# ===================================================================
# Module-level worker functions for parallel chi*T evaluation
# ===================================================================
def _worker_baseline(ct):
    """Evaluate baseline for a single chi*T value."""
    model, frame = _ensure_model()
    T = duration_from_chi_t(ct)
    pulses, cmap, Ttot = build_baseline_multitone(model, frame, T)
    return compute_detailed_metrics(model, frame, pulses, cmap, T)


def build_baseline_multitone(model, frame, T):
    return build_allbranch_multitone_pulse(
        model, frame, T, np.ones(LOGICAL_N), np.zeros(LOGICAL_N))


def _worker_opt_nm(args):
    """Optimize with NM for a single chi*T value. args = (ct, use_detuning)."""
    ct, use_detuning = args
    model, frame = _ensure_model()
    T = duration_from_chi_t(ct)

    n_params = LOGICAL_N * 2
    if use_detuning:
        n_params += LOGICAL_N

    def cost(params):
        amps = np.clip(params[:LOGICAL_N], 0.0, 3.0)
        phs = params[LOGICAL_N:2 * LOGICAL_N]
        det = params[2 * LOGICAL_N:] if use_detuning else None
        pulses, cmap, Ttot = build_allbranch_multitone_pulse(
            model, frame, T, amps, phs, det)
        try:
            final_states = simulate_and_extract(model, frame, pulses, cmap, Ttot)
            branches = extract_branch_unitaries(final_states, model, LOGICAL_N)
            fid = block_fidelity(branches, _TARGET_BRANCHES)
        except Exception:
            return 1.0
        return 1.0 - fid

    x0 = np.zeros(n_params)
    x0[:LOGICAL_N] = 1.0
    result = minimize(cost, x0, method="Nelder-Mead",
                      options={"maxiter": OPT_MAXITER_NM, "adaptive": True})

    # Compute detailed metrics at optimum
    opt_amps = np.clip(result.x[:LOGICAL_N], 0.0, 3.0)
    opt_phases = result.x[LOGICAL_N:2 * LOGICAL_N]
    opt_det = result.x[2 * LOGICAL_N:] if use_detuning else None
    pulses, cmap, Ttot = build_allbranch_multitone_pulse(
        model, frame, T, opt_amps, opt_phases, opt_det)
    metrics = compute_detailed_metrics(model, frame, pulses, cmap, T)
    return {"nfev": result.nfev, **metrics}


def _worker_de_nm(ct):
    """DE + NM optimization for a single chi*T value.

    DE uses parallel workers internally via Pool — but since this function
    is already running in the main process sequentially, we create a
    sub-pool for DE.
    """
    model, frame = _ensure_model()
    T = duration_from_chi_t(ct)
    n_params = LOGICAL_N * 3  # amp + phase + detuning

    bounds_de = (
        [(0.0, 3.0)] * LOGICAL_N
        + [(-np.pi, np.pi)] * LOGICAL_N
        + [(-abs(CHI) * 0.5, abs(CHI) * 0.5)] * LOGICAL_N
    )

    x0 = np.zeros(n_params)
    x0[:LOGICAL_N] = 1.0

    # DE with parallel workers
    with mp.Pool(N_WORKERS, initializer=_worker_init) as pool:
        de_result = differential_evolution(
            _de_cost, bounds_de, args=(T, True),
            maxiter=OPT_MAXITER_DE, seed=42, tol=1e-6,
            polish=False, x0=x0, mutation=(0.5, 1.0), recombination=0.7,
            workers=pool.map, updating="deferred",
        )

    # Local NM refinement (sequential)
    def cost_local(params):
        return _de_cost(params, T, True)

    nm_result = minimize(cost_local, de_result.x, method="Nelder-Mead",
                         options={"maxiter": OPT_MAXITER_NM, "adaptive": True})

    # Detailed metrics at optimum
    opt_amps = np.clip(nm_result.x[:LOGICAL_N], 0.0, 3.0)
    opt_phases = nm_result.x[LOGICAL_N:2 * LOGICAL_N]
    opt_det = nm_result.x[2 * LOGICAL_N:]
    pulses, cmap, Ttot = build_allbranch_multitone_pulse(
        model, frame, T, opt_amps, opt_phases, opt_det)
    metrics = compute_detailed_metrics(model, frame, pulses, cmap, T)
    return {
        "de_fid": 1.0 - de_result.fun,
        "de_nfev": de_result.nfev,
        "nm_nfev": nm_result.nfev,
        **metrics,
    }


def _worker_grape(ct):
    """GRAPE for a single chi*T value."""
    model, frame = _ensure_model()
    T = duration_from_chi_t(ct)

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
        full_dim=N_TR * N_CAV, indices=tuple(indices), labels=tuple(labels))

    target_sub = build_allbranch_target()
    phase_blocks = tuple((2 * n, 2 * n + 1) for n in range(LOGICAL_N))
    dt_grape = T / N_SLICES

    problem = build_control_problem_from_model(
        model, frame=frame,
        time_grid=PiecewiseConstantTimeGrid.uniform(steps=N_SLICES, dt_s=dt_grape),
        channel_specs=(
            ModelControlChannelSpec(
                name="qubit_I", target="qubit", quadratures=("I",),
                amplitude_bounds=(-AMP_BOUND, AMP_BOUND)),
            ModelControlChannelSpec(
                name="qubit_Q", target="qubit", quadratures=("Q",),
                amplitude_bounds=(-AMP_BOUND, AMP_BOUND)),
        ),
        objectives=(
            UnitaryObjective(
                target_operator=target_sub, subspace=sub,
                ignore_global_phase=True, phase_blocks=phase_blocks),
        ),
    )

    config = GrapeConfig(maxiter=GRAPE_MAXITER, seed=42)
    result = GrapeSolver(config).solve(problem)

    if "nominal_fidelity" in result.metrics:
        fid = result.metrics["nominal_fidelity"]
    elif "fidelity" in result.metrics:
        fid = result.metrics["fidelity"]
    else:
        fid = (1.0 - result.objective_value
               if result.objective_value <= 1.0 else 0.0)

    return {"fidelity": fid, "converged": bool(result.success)}


# ===================================================================
# MAIN
# ===================================================================
def main():
    print("=" * 70)
    print("All-Branch Multitone SQR Study (Multiprocessing Accelerated)")
    print(f"Workers: {N_WORKERS} / {os.cpu_count()} CPUs")
    print("=" * 70)
    t_total_start = time.time()

    DATA_DIR.mkdir(exist_ok=True)

    nC = len(CHI_T_VALUES)

    # Storage
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
        print(f"  [checkpoint saved]", flush=True)

    # ---------------------------------------------------------------
    # Phase A: Baseline — parallel across chi*T values
    # ---------------------------------------------------------------
    print("\n--- Phase A: Baseline (parallel across chiT) ---")
    t0 = time.time()
    with mp.Pool(nC, initializer=_worker_init) as pool:
        results_a = pool.map(_worker_baseline, CHI_T_VALUES.tolist())
    for i, (ct, m) in enumerate(zip(CHI_T_VALUES, results_a)):
        baseline_fid_block[i] = m["fid_block"]
        baseline_fid_strict[i] = m["fid_strict"]
        baseline_branch_fids[i] = m["branch_fids"]
        baseline_leakage[i] = m["leakage_max"]
        print(f"  chiT/2pi={ct}: F_block={m['fid_block']:.4f}, "
              f"F_strict={m['fid_strict']:.4f}", flush=True)
    print(f"  Phase A total: {time.time()-t0:.1f}s", flush=True)
    save_checkpoint()

    # ---------------------------------------------------------------
    # Phase B: Opt amp+phase (8p NM) — parallel across chi*T
    # ---------------------------------------------------------------
    print("\n--- Phase B: Opt amp+phase 8p NM (parallel across chiT) ---")
    t0 = time.time()
    args_b = [(ct, False) for ct in CHI_T_VALUES.tolist()]
    with mp.Pool(nC, initializer=_worker_init) as pool:
        results_b = pool.map(_worker_opt_nm, args_b)
    for i, (ct, m) in enumerate(zip(CHI_T_VALUES, results_b)):
        opt_ap_fid_block[i] = m["fid_block"]
        opt_ap_fid_strict[i] = m["fid_strict"]
        opt_ap_branch_fids[i] = m["branch_fids"]
        opt_ap_leakage[i] = m["leakage_max"]
        print(f"  chiT/2pi={ct}: F_block={m['fid_block']:.4f}, "
              f"F_strict={m['fid_strict']:.4f} ({m['nfev']} evals)", flush=True)
    print(f"  Phase B total: {time.time()-t0:.1f}s", flush=True)
    save_checkpoint()

    # ---------------------------------------------------------------
    # Phase C: Opt amp+phase+det (12p NM) — parallel across chi*T
    # ---------------------------------------------------------------
    print("\n--- Phase C: Opt amp+phase+det 12p NM (parallel across chiT) ---")
    t0 = time.time()
    args_c = [(ct, True) for ct in CHI_T_VALUES.tolist()]
    with mp.Pool(nC, initializer=_worker_init) as pool:
        results_c = pool.map(_worker_opt_nm, args_c)
    for i, (ct, m) in enumerate(zip(CHI_T_VALUES, results_c)):
        opt_apd_fid_block[i] = m["fid_block"]
        opt_apd_fid_strict[i] = m["fid_strict"]
        opt_apd_branch_fids[i] = m["branch_fids"]
        opt_apd_leakage[i] = m["leakage_max"]
        print(f"  chiT/2pi={ct}: F_block={m['fid_block']:.4f}, "
              f"F_strict={m['fid_strict']:.4f} ({m['nfev']} evals)", flush=True)
    print(f"  Phase C total: {time.time()-t0:.1f}s", flush=True)
    save_checkpoint()

    # ---------------------------------------------------------------
    # Phase D: DE + NM (12p) — sequential over chi*T, parallel DE within
    # ---------------------------------------------------------------
    print(f"\n--- Phase D: DE+NM 12p (DE parallel w/ {N_WORKERS} workers) ---")
    for i, ct in enumerate(CHI_T_VALUES):
        T_s = duration_from_chi_t(ct)
        print(f"  chiT/2pi={ct}: T={T_s*1e9:.1f} ns", flush=True)
        t0 = time.time()
        try:
            m = _worker_de_nm(ct)
            opt_de_fid_block[i] = m["fid_block"]
            opt_de_fid_strict[i] = m["fid_strict"]
            opt_de_branch_fids[i] = m["branch_fids"]
            opt_de_leakage[i] = m["leakage_max"]
            elapsed = time.time() - t0
            print(f"    DE: F_block={m['de_fid']:.6f} ({m['de_nfev']} evals)", flush=True)
            print(f"    NM: F_block={m['fid_block']:.6f} ({m['nm_nfev']} evals)", flush=True)
            print(f"    -> F_block={m['fid_block']:.4f}, "
                  f"F_strict={m['fid_strict']:.4f} ({elapsed:.1f}s)", flush=True)
        except Exception as e:
            print(f"    FAILED: {e}", flush=True)
        gc.collect()
        save_checkpoint()

    # ---------------------------------------------------------------
    # Phase E: GRAPE — parallel across chi*T values
    # ---------------------------------------------------------------
    print("\n--- Phase E: GRAPE (parallel across chiT) ---")
    t0 = time.time()
    with mp.Pool(nC, initializer=_worker_init) as pool:
        results_e = pool.map(_worker_grape, CHI_T_VALUES.tolist())
    for i, (ct, m) in enumerate(zip(CHI_T_VALUES, results_e)):
        grape_fid[i] = m["fidelity"]
        grape_converged[i] = m["converged"]
        print(f"  chiT/2pi={ct}: F_GRAPE={m['fidelity']:.8f}, "
              f"converged={m['converged']}", flush=True)
    print(f"  Phase E total: {time.time()-t0:.1f}s", flush=True)
    save_checkpoint()

    # ---------------------------------------------------------------
    # Summary
    # ---------------------------------------------------------------
    total_elapsed = time.time() - t_total_start
    print(f"\n{'='*70}")
    print(f"SUMMARY  (total wall time: {total_elapsed:.0f}s = "
          f"{total_elapsed/60:.1f} min)")
    print(f"{'='*70}")
    hdr = (f"{'chiT/2pi':>8}  {'Base_Fb':>8}  {'NM8_Fb':>8}  "
           f"{'NM12_Fb':>8}  {'DE_Fb':>8}  {'DE_Fs':>8}  {'GRAPE':>10}")
    print(hdr)
    print("-" * len(hdr))
    for i, ct in enumerate(CHI_T_VALUES):
        g = f"{grape_fid[i]:.6f}" if not np.isnan(grape_fid[i]) else "   N/A   "
        print(f"{ct:8.1f}  {baseline_fid_block[i]:8.4f}  "
              f"{opt_ap_fid_block[i]:8.4f}  {opt_apd_fid_block[i]:8.4f}  "
              f"{opt_de_fid_block[i]:8.4f}  {opt_de_fid_strict[i]:8.4f}  "
              f"{g:>10}")

    print("\nDone.")


if __name__ == "__main__":
    mp.freeze_support()  # Required on Windows for multiprocessing
    main()
