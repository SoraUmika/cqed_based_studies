"""
Re-run GRAPE benchmark with AWG hardware constraints and collect
pulse waveforms for all families.

Changes from the original run_grape_benchmark.py:
  - Qubit frequency: 6.150 GHz (updated in common.py)
  - AWG sampling rate: 1 GHz (1 ns time resolution)
  - Maximum qubit-drive Rabi rate: 2pi x 100 MHz
  - GRAPE amplitude bound set to 100 MHz (was 50 MHz)
  - Saves I/Q waveform arrays for all pulse families and GRAPE

The 1 GHz AWG constraint is automatically satisfied: each GRAPE slice
(~7-37 ns depending on chi*T) is well above the 1 ns minimum.  The
primary hardware limitation is the 100 MHz peak Rabi rate, enforced
via L-BFGS-B box constraints on each I/Q channel.

Usage:
    python scripts/run_hwgrape_and_waveforms.py

Output:
    data/hwgrape_results.npz
    data/pulse_waveforms.npz
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import block_diag

SCRIPT_DIR = Path(__file__).resolve().parent
STUDY_DIR = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

import runtime_compat  # noqa: F401

from common import (
    CHI, DT, N_CAV, N_FOCK, N_TR, TARGET_N0, THETA_TARGET, PHI_TARGET,
    build_frame, build_model, duration_from_chi_t, target_qubit_unitary,
)
from cqed_sim.core.frequencies import (
    carrier_for_transition_frequency,
    manifold_transition_frequency,
)
from cqed_sim.pulses.envelopes import (
    normalized_gaussian,
    multitone_gaussian_envelope,
    MultitoneTone,
)
from cqed_sim.sequence import SequenceCompiler
from cqed_sim.sim import SimulationConfig, prepare_simulation

# ── cqed_sim optimal control ──
from cqed_sim import (
    GrapeConfig,
    GrapeSolver,
    ModelControlChannelSpec,
    PiecewiseConstantTimeGrid,
    UnitaryObjective,
    build_control_problem_from_model,
)
from cqed_sim.unitary_synthesis import Subspace

# ===================================================================
# AWG hardware constants
# ===================================================================
AWG_SAMPLE_RATE_HZ = 1.0e9           # 1 GHz → 1 ns per sample
AWG_DT = 1.0 / AWG_SAMPLE_RATE_HZ   # 1 ns
MAX_RABI_HZ = 100e6                  # 100 MHz peak Rabi rate
AMP_BOUND = 2 * np.pi * MAX_RABI_HZ # rad/s

# ===================================================================
# GRAPE settings
# ===================================================================
GRAPE_CHI_T = np.array([1, 2, 3, 5, 7, 10], dtype=float)
N_SLICES = 48       # piecewise-constant segments
GRAPE_MAXITER = 400
GRAPE_SEED = 42

# ===================================================================
# Waveform-export chi*T (pick one in 97-100% fidelity range)
# ===================================================================
WAVEFORM_CHI_T = 3.0   # typically gives excellent GRAPE fidelity

# ===================================================================
# Helpers
# ===================================================================
SIGMA_FRAC = 1.0 / 6.0


def _build_subspace():
    """Build 8-dim qubit-cavity subspace for n=0..3, g/e."""
    indices, labels = [], []
    for n in range(N_FOCK):
        indices.append(0 * N_CAV + n)
        indices.append(1 * N_CAV + n)
        labels.append(f"|g,{n}>")
        labels.append(f"|e,{n}>")
    return Subspace(full_dim=N_TR * N_CAV, indices=tuple(indices),
                    labels=tuple(labels))


def _build_sqr_target():
    """Block-diagonal target: R_X(pi) on branch n0, identity on spectators."""
    R = target_qubit_unitary(THETA_TARGET, PHI_TARGET)
    I2 = np.eye(2, dtype=np.complex128)
    blocks = [R if n == TARGET_N0 else I2 for n in range(N_FOCK)]
    return block_diag(*blocks)


def _phase_blocks():
    return tuple((2 * n, 2 * n + 1) for n in range(N_FOCK))


def _extract_fidelity(result):
    if "nominal_fidelity" in result.metrics:
        return result.metrics["nominal_fidelity"]
    if "fidelity" in result.metrics:
        return result.metrics["fidelity"]
    return 1.0 - result.objective_value if result.objective_value <= 1.0 else 0.0


# ===================================================================
# GRAPE — hardware-constrained (100 MHz amplitude bound)
# ===================================================================
def run_hw_grape(model, frame, chi_t_2pi, cphase=True):
    """Run GRAPE with AWG hardware constraints for one chi*T value.

    Hardware constraints:
      - Hard amplitude bound: ±2π × 100 MHz on each I/Q channel
        (enforced by L-BFGS-B box constraints)
      - 1 GHz AWG sampling (1 ns) is inherently satisfied since each
        GRAPE slice is ≥ 7 ns.
    """
    T = duration_from_chi_t(chi_t_2pi)
    dt_grape = T / N_SLICES

    sub = _build_subspace()
    target = _build_sqr_target()

    obj_kwargs = dict(
        target_operator=target, subspace=sub,
        ignore_global_phase=True,
    )
    if cphase:
        obj_kwargs["phase_blocks"] = _phase_blocks()

    problem = build_control_problem_from_model(
        model, frame=frame,
        time_grid=PiecewiseConstantTimeGrid.uniform(
            steps=N_SLICES, dt_s=dt_grape),
        channel_specs=(
            ModelControlChannelSpec(
                name="qubit_I", target="qubit", quadratures=("I",),
                amplitude_bounds=(-AMP_BOUND, AMP_BOUND)),
            ModelControlChannelSpec(
                name="qubit_Q", target="qubit", quadratures=("Q",),
                amplitude_bounds=(-AMP_BOUND, AMP_BOUND)),
        ),
        objectives=(UnitaryObjective(**obj_kwargs),),
    )

    config = GrapeConfig(
        maxiter=GRAPE_MAXITER, seed=GRAPE_SEED,
    )
    result = GrapeSolver(config).solve(problem)
    return result


# ===================================================================
# GRAPE — unconstrained (no amplitude bound, for comparison)
# ===================================================================
AMP_BOUND_UNC = 2 * np.pi * 200e6  # 200 MHz (relaxed, for comparison)


def run_unconstrained_grape(model, frame, chi_t_2pi, cphase=True):
    """Run GRAPE with no meaningful amplitude bound for comparison."""
    T = duration_from_chi_t(chi_t_2pi)
    dt_grape = T / N_SLICES
    sub = _build_subspace()
    target = _build_sqr_target()

    obj_kwargs = dict(target_operator=target, subspace=sub,
                      ignore_global_phase=True)
    if cphase:
        obj_kwargs["phase_blocks"] = _phase_blocks()

    problem = build_control_problem_from_model(
        model, frame=frame,
        time_grid=PiecewiseConstantTimeGrid.uniform(
            steps=N_SLICES, dt_s=dt_grape),
        channel_specs=(
            ModelControlChannelSpec(
                name="qubit_I", target="qubit", quadratures=("I",),
                amplitude_bounds=(-AMP_BOUND_UNC, AMP_BOUND_UNC)),
            ModelControlChannelSpec(
                name="qubit_Q", target="qubit", quadratures=("Q",),
                amplitude_bounds=(-AMP_BOUND_UNC, AMP_BOUND_UNC)),
        ),
        objectives=(UnitaryObjective(**obj_kwargs),),
    )
    config = GrapeConfig(maxiter=GRAPE_MAXITER, seed=GRAPE_SEED)
    result = GrapeSolver(config).solve(problem)
    return result


# ===================================================================
# Parametric pulse waveform samplers
# ===================================================================
def sample_gaussian_waveform(model, frame, duration):
    """Sample a Gaussian SQR pulse waveform at 1 ns resolution."""
    omega_n0 = manifold_transition_frequency(model, TARGET_N0, frame)
    carrier = carrier_for_transition_frequency(omega_n0)
    amp = THETA_TARGET / (2 * duration)

    n_pts = int(round(duration / AWG_DT))
    t_rel = np.linspace(0, 1, n_pts, endpoint=False)
    t_ns = np.arange(n_pts) * (AWG_DT * 1e9)

    env = normalized_gaussian(t_rel, sigma_fraction=SIGMA_FRAC)
    # Baseband I/Q: c(t) = amp * env(t) * exp(i*(carrier*t + phase))
    t_s = t_rel * duration
    iq = amp * env * np.exp(1j * (carrier * t_s + PHI_TARGET))
    return t_ns, np.real(iq), np.imag(iq)


def sample_square_waveform(model, frame, duration):
    """Sample a square SQR pulse waveform at 1 ns resolution."""
    omega_n0 = manifold_transition_frequency(model, TARGET_N0, frame)
    carrier = carrier_for_transition_frequency(omega_n0)
    amp = THETA_TARGET / (2 * duration)

    n_pts = int(round(duration / AWG_DT))
    t_ns = np.arange(n_pts) * (AWG_DT * 1e9)
    t_s = np.arange(n_pts) * AWG_DT

    iq = amp * np.exp(1j * (carrier * t_s + PHI_TARGET))
    return t_ns, np.real(iq), np.imag(iq)


def sample_cosine_squared_waveform(model, frame, duration):
    """Sample a cosine-squared (Hann) SQR pulse waveform at 1 ns resolution."""
    omega_n0 = manifold_transition_frequency(model, TARGET_N0, frame)
    carrier = carrier_for_transition_frequency(omega_n0)
    amp = THETA_TARGET / (2 * duration)

    n_pts = int(round(duration / AWG_DT))
    t_rel = np.linspace(0, 1, n_pts, endpoint=False)
    t_ns = np.arange(n_pts) * (AWG_DT * 1e9)
    t_s = t_rel * duration

    env = 2.0 * np.cos(np.pi * (t_rel - 0.5))**2
    iq = amp * env * np.exp(1j * (carrier * t_s + PHI_TARGET))
    return t_ns, np.real(iq), np.imag(iq)


def sample_multitone_waveform(model, frame, duration):
    """Sample a multitone Gaussian SQR pulse with 4 tones at 1 ns."""
    n_pts = int(round(duration / AWG_DT))
    t_rel = np.linspace(0, 1, n_pts, endpoint=False)
    t_ns = np.arange(n_pts) * (AWG_DT * 1e9)

    tone_specs = []
    for n in range(N_FOCK):
        omega_n = manifold_transition_frequency(model, n, frame)
        tone_carrier = carrier_for_transition_frequency(omega_n)
        amp_n = THETA_TARGET / (2 * duration) if n == TARGET_N0 else 0.0
        tone_specs.append(MultitoneTone(
            manifold=n, omega_rad_s=tone_carrier,
            amp_rad_s=amp_n,
            phase_rad=PHI_TARGET if n == TARGET_N0 else 0.0,
        ))

    env = multitone_gaussian_envelope(
        t_rel, duration_s=duration, sigma_fraction=SIGMA_FRAC,
        tone_specs=tone_specs,
    )
    return t_ns, np.real(env), np.imag(env)


# ===================================================================
# MAIN
# ===================================================================
def main():
    print("=" * 70)
    print("Hardware-Constrained GRAPE + Pulse Waveform Export")
    print(f"AWG: {AWG_SAMPLE_RATE_HZ/1e9:.0f} GHz sampling, "
          f"max Rabi = {MAX_RABI_HZ/1e6:.0f} MHz")
    print("=" * 70)

    model = build_model()
    frame = build_frame(model)
    data_dir = STUDY_DIR / "data"
    data_dir.mkdir(exist_ok=True)

    print(f"\nSystem: omega_q/2pi = {model.omega_q / (2*np.pi) / 1e9:.3f} GHz")
    print(f"        omega_c/2pi = {model.omega_c / (2*np.pi) / 1e9:.3f} GHz")
    print(f"        chi/2pi     = {CHI / (2*np.pi) / 1e6:.2f} MHz")

    # ==============================================================
    # Part 1: GRAPE scan — HW-constrained cphase + true SQR
    # ==============================================================
    n_ct = len(GRAPE_CHI_T)
    hw_fid_cphase = np.full(n_ct, np.nan)
    hw_fid_true = np.full(n_ct, np.nan)
    hw_converged_cphase = np.full(n_ct, False)
    hw_converged_true = np.full(n_ct, False)

    unc_fid_cphase = np.full(n_ct, np.nan)
    unc_fid_true = np.full(n_ct, np.nan)
    unc_converged_cphase = np.full(n_ct, False)
    unc_converged_true = np.full(n_ct, False)

    # Store the GRAPE result for the waveform-export chi*T
    grape_result_for_waveform = None
    grape_T_for_waveform = None

    print("\n--- GRAPE sweep (HW-constrained + unconstrained) ---")
    for i, ct in enumerate(GRAPE_CHI_T):
        T = duration_from_chi_t(ct)
        n_steps = max(1, int(round(T / AWG_DT)))
        print(f"\n  chiT/2pi={ct:.0f}: T={T*1e9:.1f} ns, "
              f"AWG samples={n_steps}")

        # ── HW-constrained cphase ──
        print(f"    HW cphase ...", end=" ", flush=True)
        t0 = time.time()
        try:
            res = run_hw_grape(model, frame, ct, cphase=True)
            hw_fid_cphase[i] = _extract_fidelity(res)
            hw_converged_cphase[i] = res.success
            print(f"F={hw_fid_cphase[i]:.8f} ({'ok' if res.success else 'FAIL'}) "
                  f"[{time.time()-t0:.1f}s]")
            if ct == WAVEFORM_CHI_T:
                grape_result_for_waveform = res
                grape_T_for_waveform = T
        except Exception as e:
            print(f"ERROR: {e}")

        # ── HW-constrained true SQR ──
        print(f"    HW true  ...", end=" ", flush=True)
        t0 = time.time()
        try:
            res = run_hw_grape(model, frame, ct, cphase=False)
            hw_fid_true[i] = _extract_fidelity(res)
            hw_converged_true[i] = res.success
            print(f"F={hw_fid_true[i]:.8f} ({'ok' if res.success else 'FAIL'}) "
                  f"[{time.time()-t0:.1f}s]")
        except Exception as e:
            print(f"ERROR: {e}")

        # ── Unconstrained cphase ──
        print(f"    Unc cphase ...", end=" ", flush=True)
        t0 = time.time()
        try:
            res = run_unconstrained_grape(model, frame, ct, cphase=True)
            unc_fid_cphase[i] = _extract_fidelity(res)
            unc_converged_cphase[i] = res.success
            print(f"F={unc_fid_cphase[i]:.8f} ({'ok' if res.success else 'FAIL'}) "
                  f"[{time.time()-t0:.1f}s]")
        except Exception as e:
            print(f"ERROR: {e}")

        # ── Unconstrained true SQR ──
        print(f"    Unc true  ...", end=" ", flush=True)
        t0 = time.time()
        try:
            res = run_unconstrained_grape(model, frame, ct, cphase=False)
            unc_fid_true[i] = _extract_fidelity(res)
            unc_converged_true[i] = res.success
            print(f"F={unc_fid_true[i]:.8f} ({'ok' if res.success else 'FAIL'}) "
                  f"[{time.time()-t0:.1f}s]")
        except Exception as e:
            print(f"ERROR: {e}")

    # ==============================================================
    # Part 2: Extract GRAPE waveform at WAVEFORM_CHI_T
    # ==============================================================
    grape_t_ns = np.array([])
    grape_I = np.array([])
    grape_Q = np.array([])

    if grape_result_for_waveform is not None:
        res = grape_result_for_waveform
        T = grape_T_for_waveform
        # physical_values: shape (n_channels, n_steps)
        # channels: qubit_I (idx 0), qubit_Q (idx 1)
        phys = res.physical_values
        if phys is not None and phys.ndim == 2:
            n_steps = phys.shape[1]
            grape_t_ns = np.arange(n_steps) * (T / n_steps) * 1e9
            grape_I = phys[0, :]
            grape_Q = phys[1, :]
            print(f"\nGRAPE waveform at chiT/2pi={WAVEFORM_CHI_T}: "
                  f"{n_steps} samples, peak I={np.max(np.abs(grape_I))/(2*np.pi*1e6):.1f} MHz, "
                  f"peak Q={np.max(np.abs(grape_Q))/(2*np.pi*1e6):.1f} MHz")
        else:
            # Fall back to command_values
            cmd = res.command_values
            if cmd is not None and cmd.ndim == 2:
                n_steps = cmd.shape[1]
                grape_t_ns = np.arange(n_steps) * (T / n_steps) * 1e9
                grape_I = cmd[0, :]
                grape_Q = cmd[1, :]
                print(f"\nGRAPE waveform (command) at chiT/2pi={WAVEFORM_CHI_T}: "
                      f"{n_steps} samples")

    # ==============================================================
    # Part 3: Sample parametric pulse waveforms
    # ==============================================================
    T_wf = duration_from_chi_t(WAVEFORM_CHI_T)
    print(f"\nSampling parametric waveforms at chiT/2pi={WAVEFORM_CHI_T} "
          f"(T={T_wf*1e9:.1f} ns)")

    gauss_t, gauss_I, gauss_Q = sample_gaussian_waveform(model, frame, T_wf)
    sq_t, sq_I, sq_Q = sample_square_waveform(model, frame, T_wf)
    cos2_t, cos2_I, cos2_Q = sample_cosine_squared_waveform(model, frame, T_wf)
    mt_t, mt_I, mt_Q = sample_multitone_waveform(model, frame, T_wf)

    # ==============================================================
    # Save results
    # ==============================================================
    np.savez_compressed(
        str(data_dir / "hwgrape_results.npz"),
        grape_chi_t=GRAPE_CHI_T,
        hw_fid_cphase=hw_fid_cphase,
        hw_fid_true=hw_fid_true,
        hw_converged_cphase=hw_converged_cphase,
        hw_converged_true=hw_converged_true,
        unc_fid_cphase=unc_fid_cphase,
        unc_fid_true=unc_fid_true,
        unc_converged_cphase=unc_converged_cphase,
        unc_converged_true=unc_converged_true,
        awg_sample_rate_hz=AWG_SAMPLE_RATE_HZ,
        max_rabi_hz=MAX_RABI_HZ,
        amp_bound_rad_s=AMP_BOUND,
        waveform_chi_t=WAVEFORM_CHI_T,
    )
    print(f"\nGRAPE results saved to data/hwgrape_results.npz")

    np.savez_compressed(
        str(data_dir / "pulse_waveforms.npz"),
        waveform_chi_t=WAVEFORM_CHI_T,
        duration_ns=T_wf * 1e9,
        # Gaussian
        gauss_t_ns=gauss_t, gauss_I=gauss_I, gauss_Q=gauss_Q,
        # Square
        sq_t_ns=sq_t, sq_I=sq_I, sq_Q=sq_Q,
        # Cosine-squared
        cos2_t_ns=cos2_t, cos2_I=cos2_I, cos2_Q=cos2_Q,
        # Multitone
        mt_t_ns=mt_t, mt_I=mt_I, mt_Q=mt_Q,
        # GRAPE (HW-constrained)
        grape_t_ns=grape_t_ns, grape_I=grape_I, grape_Q=grape_Q,
    )
    print(f"Waveform data saved to data/pulse_waveforms.npz")

    # ==============================================================
    # Summary table
    # ==============================================================
    print(f"\n{'='*78}")
    print("GRAPE FIDELITY SUMMARY")
    print(f"{'='*78}")
    hdr = (f"{'chiT':>6}  {'HW_cph':>10}  {'HW_true':>10}  "
           f"{'Unc_cph':>10}  {'Unc_true':>10}")
    print(hdr)
    print("-" * len(hdr))
    for i, ct in enumerate(GRAPE_CHI_T):
        def fmt(v):
            return f"{v:.6f}" if not np.isnan(v) else "   N/A   "
        print(f"{ct:6.0f}  {fmt(hw_fid_cphase[i]):>10}  "
              f"{fmt(hw_fid_true[i]):>10}  "
              f"{fmt(unc_fid_cphase[i]):>10}  "
              f"{fmt(unc_fid_true[i]):>10}")

    print("\nDone.")


if __name__ == "__main__":
    main()
