"""
Phase 2 benchmark: GRAPE piecewise-optimized pulse for SQR.

Uses the cqed_sim optimal control module to find the best achievable fidelity
for SQR targets on the truncated g/e subspace, as a function of χT/(2π).

Two GRAPE targets are optimised:
  1. **True SQR** — identity on all spectator branches (global phase only).
  2. **Cphase SQR** — per-branch Z-freedom via ``phase_blocks``, matching the
     metric used for the parameterised pulse families.

The conditional-phase GRAPE establishes the theoretical upper bound against
which analytic waveform families are compared.

Usage:
    python scripts/run_grape_benchmark.py

Output:
    data/grape_benchmark_results.npz
"""

import sys
import time
from pathlib import Path

import numpy as np
from scipy.linalg import block_diag

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import runtime_compat  # noqa: F401

from common import (
    CHI, N_CAV, N_FOCK, N_TR, TARGET_N0, THETA_TARGET, PHI_TARGET,
    build_frame, build_model, duration_from_chi_t, target_qubit_unitary,
)

from cqed_sim import (
    GrapeConfig,
    GrapeSolver,
    ModelControlChannelSpec,
    PiecewiseConstantTimeGrid,
    UnitaryObjective,
    build_control_problem_from_model,
)
from cqed_sim.unitary_synthesis import Subspace

# χT/(2π) values for GRAPE (fewer points since GRAPE is expensive)
GRAPE_CHI_T = np.array([1, 2, 3, 5, 7, 10], dtype=float)
N_SLICES = 48       # piecewise-constant segments
GRAPE_MAXITER = 300
AMP_BOUND = 2 * np.pi * 50e6  # max drive amplitude (rad/s)


def build_sqr_target_subspace():
    """Build the SQR target in the g/e qubit-cavity subspace.

    Returns (target_8x8, subspace, phase_blocks) where:
      - target_8x8 is block_diag(I, I, ..., R, ..., I) with R on branch n0.
      - subspace  selects |g,0⟩,|e,0⟩,...,|g,N_FOCK-1⟩,|e,N_FOCK-1⟩.
      - phase_blocks is a tuple-of-tuples for per-Fock phase freedom.
    """
    # Build subspace manually to support n_tr >= 3
    # Indices: |g,n⟩ → n, |e,n⟩ → N_CAV + n  (qubit-first tensor ordering)
    indices = []
    labels = []
    for n in range(N_FOCK):
        indices.append(0 * N_CAV + n)       # |g,n⟩
        indices.append(1 * N_CAV + n)       # |e,n⟩
        labels.append(f"|g,{n}>")
        labels.append(f"|e,{n}>")
    sub = Subspace(
        full_dim=N_TR * N_CAV,
        indices=tuple(indices),
        labels=tuple(labels),
    )

    # Build block-diagonal target: 2×2 per Fock level
    R_target = target_qubit_unitary(THETA_TARGET, PHI_TARGET)
    I2 = np.eye(2, dtype=np.complex128)
    blocks = [R_target if n == TARGET_N0 else I2 for n in range(N_FOCK)]
    target_sub = block_diag(*blocks)

    # Per-Fock phase blocks for cphase GRAPE
    phase_blocks = tuple((2 * n, 2 * n + 1) for n in range(N_FOCK))

    return target_sub, sub, phase_blocks


def _build_grape_problem(model, frame, dt_grape, target_sub, sub, cphase):
    """Assemble GRAPE control problem for true SQR or cphase SQR."""
    obj_kwargs = dict(
        target_operator=target_sub,
        subspace=sub,
        ignore_global_phase=True,
    )
    if cphase:
        _, _, phase_blocks = build_sqr_target_subspace()
        obj_kwargs["phase_blocks"] = phase_blocks

    return build_control_problem_from_model(
        model,
        frame=frame,
        time_grid=PiecewiseConstantTimeGrid.uniform(steps=N_SLICES, dt_s=dt_grape),
        channel_specs=(
            ModelControlChannelSpec(
                name="qubit_I",
                target="qubit",
                quadratures=("I",),
                amplitude_bounds=(-AMP_BOUND, AMP_BOUND),
            ),
            ModelControlChannelSpec(
                name="qubit_Q",
                target="qubit",
                quadratures=("Q",),
                amplitude_bounds=(-AMP_BOUND, AMP_BOUND),
            ),
        ),
        objectives=(UnitaryObjective(**obj_kwargs),),
    )


def _extract_fidelity(result):
    """Pull best fidelity from GRAPE result."""
    if "nominal_fidelity" in result.metrics:
        return result.metrics["nominal_fidelity"]
    if "fidelity" in result.metrics:
        return result.metrics["fidelity"]
    obj_val = result.objective_value
    return 1.0 - obj_val if obj_val <= 1.0 else 0.0


def run_grape_for_chi_t(chi_t_2pi, cphase=False):
    """Run GRAPE for a single χT/(2π) value.

    Parameters
    ----------
    chi_t_2pi : float
        Dimensionless χT/(2π).
    cphase : bool
        If True, use per-branch phase_blocks (conditional-phase metric).

    Returns dict with metrics.
    """
    T = duration_from_chi_t(chi_t_2pi)
    dt_grape = T / N_SLICES

    model = build_model()
    frame = build_frame(model)
    target_sub, sub, _ = build_sqr_target_subspace()

    problem = _build_grape_problem(model, frame, dt_grape, target_sub, sub, cphase)
    config = GrapeConfig(maxiter=GRAPE_MAXITER, seed=42)
    result = GrapeSolver(config).solve(problem)

    return {
        "chi_t_2pi": chi_t_2pi,
        "duration_s": T,
        "fidelity": _extract_fidelity(result),
        "objective_value": result.objective_value,
        "converged": result.success,
        "message": result.message,
    }


def run_scan():
    """Run GRAPE benchmark (true SQR and cphase SQR) over χT/(2π) values."""
    print(f"GRAPE benchmark: {len(GRAPE_CHI_T)} χT/(2π) values, "
          f"{N_SLICES} slices, maxiter={GRAPE_MAXITER}")
    print(f"Target: n0={TARGET_N0}, θ=π SQR")
    print()

    true_results = []
    cphase_results = []
    t_start = time.time()

    for chi_t in GRAPE_CHI_T:
        # --- True SQR ---
        print(f"  χT/(2π) = {chi_t:.1f}  [true SQR] ...", end=" ", flush=True)
        t0 = time.time()
        res_true = run_grape_for_chi_t(chi_t, cphase=False)
        print(f"F = {res_true['fidelity']:.8f}  ({time.time() - t0:.1f} s)")
        true_results.append(res_true)

        # --- Cphase SQR ---
        print(f"  χT/(2π) = {chi_t:.1f}  [cphase SQR] ...", end=" ", flush=True)
        t0 = time.time()
        res_cp = run_grape_for_chi_t(chi_t, cphase=True)
        print(f"F = {res_cp['fidelity']:.8f}  ({time.time() - t0:.1f} s)")
        cphase_results.append(res_cp)
        print()

    elapsed = time.time() - t_start
    print(f"GRAPE scan complete in {elapsed:.1f} s")

    # Save
    data_dir = SCRIPT_DIR.parent / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    np.savez(
        data_dir / "grape_benchmark_results.npz",
        chi_t_values=GRAPE_CHI_T,
        fidelity_true=np.array([r["fidelity"] for r in true_results]),
        fidelity_cphase=np.array([r["fidelity"] for r in cphase_results]),
        objective_true=np.array([r["objective_value"] for r in true_results]),
        objective_cphase=np.array([r["objective_value"] for r in cphase_results]),
        converged_true=np.array([r["converged"] for r in true_results]),
        converged_cphase=np.array([r["converged"] for r in cphase_results]),
        n_slices=N_SLICES,
        maxiter=GRAPE_MAXITER,
        amp_bound=AMP_BOUND,
        n_fock=N_FOCK,
        target_n0=TARGET_N0,
        theta_target=THETA_TARGET,
    )
    print(f"Results saved to {data_dir / 'grape_benchmark_results.npz'}")


if __name__ == "__main__":
    run_scan()
