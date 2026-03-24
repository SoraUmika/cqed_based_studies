"""
black_box.py — Black-box (model-free) optimizer baseline.

Implements a direct optimization of control amplitudes on the truth model without
any model learning or GRAPE gradients. Uses scipy differential_evolution as the
optimizer, which is gradient-free and suitable for the black-box scenario.

Physical setup:
  The control is a piecewise-constant pulse with n_steps slices of dt_s duration.
  Two quadrature channels (I, Q) are optimized independently.
  Amplitudes are bounded to ±amp_bound.

Evaluation approach:
  For each candidate parameter vector, we:
  1. Build piecewise-constant pulse amplitudes as arrays.
  2. Use evaluate_control_with_simulator directly with the truth model.
  3. Return the negative average fidelity as the cost (minimization).

Note: The black-box optimizer does NOT build a ControlProblem or use gradients.
It only uses cqed_sim for evaluation via evaluate_control_with_simulator.
This is the intended use: the optimizer treats the device as a black box.

Implementation details:
  - n_steps=8 (fewer than GRAPE for tractability: 8*2=16 parameters)
  - dt_s=10ns (same slice duration as GRAPE)
  - Total gate time: 80 ns (half of GRAPE due to fewer slices)
  - differential_evolution: popsize=15, maxiter=20 (~300 evaluations)
  - Each evaluation calls evaluate_control_with_simulator (expensive for large dim)

To keep runtime tractable, we use the direct evaluate_control_with_simulator API
which internally calls the QuTiP solver. The problem is constructed once from the
truth model and reused across the optimizer loop.
"""

from __future__ import annotations

import time

import numpy as np
from scipy.optimize import differential_evolution

from cqed_sim import (
    ControlEvaluationCase,
    DispersiveTransmonCavityModel,
    FrameSpec,
    ModelControlChannelSpec,
    PiecewiseConstantTimeGrid,
    UnitaryObjective,
    build_control_problem_from_model,
    evaluate_control_with_simulator,
)
from cqed_sim.unitary_synthesis import Subspace


def run_black_box(
    truth_model: DispersiveTransmonCavityModel,
    frame: FrameSpec,
    grape_subspace: Subspace,
    target_matrix: np.ndarray,
    n_evaluations: int = 200,
    seed: int = 42,
    n_steps: int = 8,
    dt_s: float = 10e-9,
    amp_bound: float = 2 * np.pi * 50e6,
    popsize: int = 15,
    maxiter_de: int = 20,
) -> dict:
    """
    Run a black-box optimizer on the truth model without model learning.

    Uses scipy.optimize.differential_evolution with the cqed_sim simulator
    as the oracle. No GRAPE gradients are used. The optimizer only sees
    fidelity values returned by the evaluator.

    The parameter space is n_steps * 2 (I and Q amplitudes for each slice),
    bounded to [-amp_bound, +amp_bound].

    Parameters
    ----------
    truth_model : DispersiveTransmonCavityModel
        The true physical model (unknown to the learner, directly probed here).
    frame : FrameSpec
        Rotating frame for simulation.
    grape_subspace : Subspace
        Control subspace (g+e manifold over Fock 0-3).
    target_matrix : np.ndarray, shape (8, 8)
        Target unitary (X gate in subspace).
    n_evaluations : int
        Approximate total oracle evaluations. Used to set popsize/maxiter.
    seed : int
        Random seed for reproducibility.
    n_steps : int
        Number of piecewise-constant time slices. Default 8 (fewer than GRAPE).
    dt_s : float
        Slice duration in seconds. Default 10 ns.
    amp_bound : float
        Maximum drive amplitude (rad/s). Default 2pi*50 MHz.
    popsize : int
        Population size multiplier for differential_evolution. Actual pop size
        = popsize * n_parameters. Default 15.
    maxiter_de : int
        Maximum generations for differential_evolution. Default 20.

    Returns
    -------
    dict with keys:
        'best_fidelity'       : float, best fidelity achieved
        'fidelity_history'    : list[float], fidelity at each generation
        'infidelity_history'  : list[float], 1 - fidelity at each generation
        'n_evaluations_used'  : int, total oracle calls made
        'best_params'         : np.ndarray, best parameter vector
        'wall_time_s'         : float, total wall time
        'converged'           : bool, whether differential_evolution converged
        'message'             : str, optimizer message
    """
    n_params = n_steps * 2  # I and Q channels
    time_grid = PiecewiseConstantTimeGrid.uniform(steps=n_steps, dt_s=dt_s)

    channel_specs = [
        ModelControlChannelSpec(
            name="qubit",
            target="qubit",
            quadratures=("I", "Q"),
            amplitude_bounds=(-float(amp_bound), float(amp_bound)),
        ),
    ]

    objective = UnitaryObjective(
        target_operator=np.asarray(target_matrix, dtype=np.complex128),
        subspace=grape_subspace,
        ignore_global_phase=True,
        allow_diagonal_phase=False,
        probe_strategy="basis_plus_uniform",
        weight=1.0,
        name="X_gate_bb",
    )

    # Build problem from truth model (for evaluation only)
    problem = build_control_problem_from_model(
        truth_model,
        frame=frame,
        time_grid=time_grid,
        channel_specs=channel_specs,
        objectives=[objective],
        penalties=[],
    )

    eval_case = ControlEvaluationCase(
        model=truth_model,
        label="black_box_eval",
        frame=frame,
        noise=None,
        weight=1.0,
    )

    fidelity_history: list[float] = []
    eval_count = [0]

    def objective_function(params: np.ndarray) -> float:
        """
        Black-box objective: return negative fidelity (for minimization).

        Converts the flat parameter vector to a (n_controls, n_slices) array
        and evaluates using evaluate_control_with_simulator.
        """
        eval_count[0] += 1
        try:
            eval_result = evaluate_control_with_simulator(
                problem,
                params.reshape(problem.parameterization.parameter_shape),
                cases=[eval_case],
            )
            member = eval_result.member_reports[0]
            fidelity = float(member.aggregate_fidelity)
        except Exception:
            fidelity = 0.0

        return -fidelity   # minimization

    def callback_fn(xk, convergence):
        """Record best fidelity at each generation."""
        f = -objective_function(xk)
        fidelity_history.append(float(f))

    # Parameter bounds: all parameters in [-amp_bound, +amp_bound]
    bounds = [(-float(amp_bound), float(amp_bound))] * n_params

    t_start = time.perf_counter()

    de_result = differential_evolution(
        objective_function,
        bounds,
        maxiter=int(maxiter_de),
        popsize=int(popsize),
        seed=int(seed),
        tol=1e-5,
        mutation=(0.5, 1.0),
        recombination=0.7,
        callback=callback_fn,
        init="sobol",
        updating="deferred",
        workers=1,
    )

    wall_time = time.perf_counter() - t_start

    best_fidelity = float(-de_result.fun)
    best_params = np.asarray(de_result.x, dtype=float)

    # Ensure the final best is recorded in history
    if not fidelity_history or abs(fidelity_history[-1] - best_fidelity) > 1e-8:
        fidelity_history.append(best_fidelity)

    return {
        "best_fidelity": best_fidelity,
        "fidelity_history": fidelity_history,
        "infidelity_history": [1.0 - f for f in fidelity_history],
        "n_evaluations_used": int(eval_count[0]),
        "best_params": best_params,
        "wall_time_s": float(wall_time),
        "converged": bool(de_result.success),
        "message": str(de_result.message),
        "n_steps": int(n_steps),
        "dt_s": float(dt_s),
        "gate_duration_s": float(n_steps * dt_s),
    }


if __name__ == "__main__":
    import sys
    sys.path.insert(0, ".")
    from models import (
        make_truth_model, make_frame, make_grape_subspace, make_target_matrix
    )

    print("Running black-box optimizer (quick test, maxiter_de=5)...")
    truth_model = make_truth_model()
    frame = make_frame(truth_model)
    subspace = make_grape_subspace()
    target = make_target_matrix()

    result = run_black_box(
        truth_model=truth_model,
        frame=frame,
        grape_subspace=subspace,
        target_matrix=target,
        n_steps=8,
        dt_s=10e-9,
        seed=42,
        maxiter_de=5,
        popsize=8,
    )

    print(f"Best fidelity: {result['best_fidelity']:.6f}")
    print(f"Evaluations used: {result['n_evaluations_used']}")
    print(f"Wall time: {result['wall_time_s']:.1f} s")
    print(f"Fidelity history: {[f'{f:.4f}' for f in result['fidelity_history']]}")
    print(f"Converged: {result['converged']}")
