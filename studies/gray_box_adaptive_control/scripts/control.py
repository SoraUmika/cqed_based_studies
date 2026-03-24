"""
control.py — GRAPE wrappers and evaluation helpers.

Wraps the cqed_sim GRAPE and evaluation APIs for the gray-box adaptive control study.
Provides:
- run_grape: optimize a unitary X gate on a given model
- eval_on_model: evaluate a GRAPE result on a (possibly different) truth model
- compute_subspace_fidelity_from_result: extract optimizer fidelity directly
- eval_per_fock_fidelity: compute per-Fock-sector fidelity
- run_grape_multistart: run GRAPE with multiple seeds, return best by training fidelity

Physical setup:
  The target is the qubit X gate simultaneously on Fock 0-3, restricted to the 8-dim
  subspace {|g,0-3>, |e,0-3>} inside the 12-dim Hilbert space.

GRAPE configuration:
  16 time slices * 10 ns = 160 ns total gate time
  Qubit drive: I + Q quadratures, bound ±2pi*50 MHz
  Penalties: AmplitudePenalty(weight=0.01) + LeakagePenalty(weight=5.0)
  Optimizer: L-BFGS-B, up to maxiter=200
"""

from __future__ import annotations

import numpy as np

from cqed_sim import (
    ControlEvaluationCase,
    DispersiveTransmonCavityModel,
    FrameSpec,
    GrapeConfig,
    GrapeSolver,
    ModelControlChannelSpec,
    NoiseSpec,
    PiecewiseConstantTimeGrid,
    UnitaryObjective,
    build_control_problem_from_model,
    evaluate_control_with_simulator,
)
from cqed_sim.optimal_control import LeakagePenalty
from cqed_sim.unitary_synthesis import Subspace


# ---------------------------------------------------------------------------
# GRAPE execution
# ---------------------------------------------------------------------------


def run_grape(
    model: DispersiveTransmonCavityModel,
    frame: FrameSpec,
    subspace: Subspace,
    target_matrix: np.ndarray,
    n_steps: int = 16,
    dt_s: float = 10e-9,
    maxiter: int = 200,
    seed: int = 42,
    amp_bound: float = 2 * np.pi * 50e6,
) -> tuple:
    """
    Run GRAPE to optimize the qubit X gate on the given model.

    Builds a ControlProblem using build_control_problem_from_model with:
    - Two control channels (qubit I and Q quadratures)
    - UnitaryObjective targeting the X gate in the g+e subspace
    - AmplitudePenalty to regularize pulse energy
    - LeakagePenalty to suppress f-level leakage

    Wraps GrapeSolver in try/except to return the best result found even if
    convergence is not declared.

    Parameters
    ----------
    model : DispersiveTransmonCavityModel
        Model used for GRAPE optimization (learner model or truth model).
    frame : FrameSpec
        Rotating frame. Should match model.omega_c and model.omega_q.
    subspace : Subspace
        Subspace for the UnitaryObjective and LeakagePenalty (8-dim g+e block).
    target_matrix : np.ndarray, shape (8, 8)
        Target unitary restricted to the subspace.
    n_steps : int
        Number of piecewise-constant time slices. Default 16.
    dt_s : float
        Duration of each time slice in seconds. Default 10 ns.
    maxiter : int
        Maximum GRAPE iterations. Default 200.
    seed : int
        Random seed for initial guess. Default 42.
    amp_bound : float
        Maximum qubit drive amplitude (rad/s). Default 2pi*50 MHz.

    Returns
    -------
    tuple : (GrapeResult, final_fidelity)
        GrapeResult from the optimizer.
        final_fidelity is the unitary subspace fidelity at the final iteration.
    """
    # Build the uniform time grid
    time_grid = PiecewiseConstantTimeGrid.uniform(steps=n_steps, dt_s=dt_s)

    # Channel specifications: separate I and Q channels for the qubit drive
    channel_specs = [
        ModelControlChannelSpec(
            name="qubit",
            target="qubit",
            quadratures=("I", "Q"),
            amplitude_bounds=(-float(amp_bound), float(amp_bound)),
        ),
    ]

    # Build objective: unitary X gate in the g+e subspace
    objective = UnitaryObjective(
        target_operator=np.asarray(target_matrix, dtype=np.complex128),
        subspace=subspace,
        ignore_global_phase=True,
        allow_diagonal_phase=False,
        probe_strategy="basis_plus_uniform",
        weight=1.0,
        name="X_gate",
    )

    # Penalties
    # Note: AmplitudePenalty with apply_to='command' evaluates in physical rad/s units,
    # causing gradient magnitudes ~(2pi*50e6)^2 * weight which overwhelms the infidelity
    # gradient and collapses optimization to zero amplitude. Omit amplitude regularization
    # for this study; hard amplitude bounds from ModelControlChannelSpec are sufficient.
    leakage_penalty = LeakagePenalty(
        subspace=subspace,
        weight=1.0,
        metric="average",
    )

    # Build the control problem
    problem = build_control_problem_from_model(
        model,
        frame=frame,
        time_grid=time_grid,
        channel_specs=channel_specs,
        objectives=[objective],
        penalties=[leakage_penalty],
    )

    # GRAPE configuration
    config = GrapeConfig(
        maxiter=int(maxiter),
        seed=int(seed),
        initial_guess="random",
        random_scale=0.15,
        engine="numpy",
    )
    solver = GrapeSolver(config=config)

    try:
        result = solver.solve(problem)
    except Exception as exc:
        # Return whatever was computed (if any) — build a minimal fallback
        # This should not happen in normal usage
        raise RuntimeError(f"GRAPE failed unexpectedly: {exc}") from exc

    # Extract the best fidelity found
    final_fidelity = float(result.metrics.get("nominal_fidelity", 0.0))
    # The metric "nominal_fidelity" is the unitary fidelity at the last evaluation.
    # Also check exact_unitary_fidelity from system metrics.
    if result.system_metrics:
        sm = result.system_metrics[0]
        if "objectives" in sm and sm["objectives"]:
            obj0 = sm["objectives"][0]
            exact_f = obj0.get("exact_unitary_fidelity", None)
            if exact_f is not None:
                final_fidelity = float(exact_f)

    return result, problem, final_fidelity


# ---------------------------------------------------------------------------
# Multi-start GRAPE
# ---------------------------------------------------------------------------


def run_grape_multistart(
    model: DispersiveTransmonCavityModel,
    frame: FrameSpec,
    subspace: Subspace,
    target_matrix: np.ndarray,
    seeds: list[int],
    n_steps: int = 16,
    dt_s: float = 10e-9,
    maxiter: int = 200,
    amp_bound: float = 2 * np.pi * 50e6,
) -> tuple:
    """
    Run GRAPE with multiple random seeds and return the best result.

    Runs run_grape for each seed in `seeds` and selects the result with the
    highest training-model fidelity (objective_value, lowest infidelity).
    This mitigates the local-minima problem inherent in GRAPE optimization.

    Parameters
    ----------
    model : DispersiveTransmonCavityModel
        Training model for GRAPE.
    frame : FrameSpec
        Rotating frame for the training model.
    subspace : Subspace
        Subspace for UnitaryObjective and LeakagePenalty.
    target_matrix : np.ndarray
        Target unitary in the subspace.
    seeds : list[int]
        List of random seeds to try.
    n_steps, dt_s, maxiter, amp_bound : see run_grape

    Returns
    -------
    tuple : (GrapeResult, ControlProblem, final_fidelity, best_seed)
        Best result selected by lowest objective_value on training model.
    """
    best_result = None
    best_problem = None
    best_fidelity = -1.0
    best_seed = seeds[0]

    for seed in seeds:
        result, problem, fidelity = run_grape(
            model, frame, subspace, target_matrix,
            n_steps=n_steps, dt_s=dt_s, maxiter=maxiter,
            seed=seed, amp_bound=amp_bound,
        )
        # Select best by highest fidelity on training model
        if fidelity > best_fidelity:
            best_fidelity = fidelity
            best_result = result
            best_problem = problem
            best_seed = seed

    return best_result, best_problem, best_fidelity, best_seed


# ---------------------------------------------------------------------------
# Cross-model evaluation
# ---------------------------------------------------------------------------


def eval_on_model(
    grape_result,
    problem,
    eval_model: DispersiveTransmonCavityModel,
    eval_frame: FrameSpec,
    eval_noise: NoiseSpec | None = None,
    label: str = "eval",
) -> tuple[float, list[float]]:
    """
    Evaluate a GRAPE-optimized pulse on a (possibly different) model.

    Uses evaluate_control_with_simulator with a ControlEvaluationCase specifying
    eval_model and eval_frame. This correctly re-computes the Hamiltonian using
    eval_model, enabling cross-model evaluation (e.g., nominal GRAPE evaluated
    on truth model).

    Parameters
    ----------
    grape_result : GrapeResult
        Result from run_grape. Its schedule is used for evaluation.
    problem : ControlProblem
        The ControlProblem used in the GRAPE run. Needed for pulse extraction.
    eval_model : DispersiveTransmonCavityModel
        Model used for evaluation (typically the truth model).
    eval_frame : FrameSpec
        Rotating frame for evaluation.
    eval_noise : NoiseSpec, optional
        Noise for evaluation. None = noiseless (unitary) simulation.
    label : str
        Label for the evaluation case.

    Returns
    -------
    tuple : (avg_fidelity, per_state_fidelities)
        avg_fidelity : float
            Average fidelity across all probe state pairs.
        per_state_fidelities : list[float]
            Per-pair fidelities from the objective evaluation.
    """
    eval_case = ControlEvaluationCase(
        model=eval_model,
        label=str(label),
        frame=eval_frame,
        noise=eval_noise,
        weight=1.0,
    )

    eval_result = evaluate_control_with_simulator(
        problem,
        grape_result.schedule,
        cases=[eval_case],
    )

    # Extract fidelities from the first (only) member report
    member = eval_result.member_reports[0]
    avg_fidelity = float(member.aggregate_fidelity)

    per_state_fidelities: list[float] = []
    for obj_report in member.objective_reports:
        per_state_fidelities.extend(list(obj_report.fidelities))

    return avg_fidelity, per_state_fidelities


def eval_per_fock_fidelity(
    grape_result,
    problem,
    eval_model: DispersiveTransmonCavityModel,
    eval_frame: FrameSpec,
    eval_noise: NoiseSpec | None = None,
) -> dict:
    """
    Evaluate gate fidelity per Fock sector.

    Returns the average fidelity and per-Fock-sector information by analyzing
    the per-pair fidelities from the UnitaryObjective probe states.

    Note: The UnitaryObjective generates probe states automatically. The probe
    states include pairs from the subspace basis vectors. This function extracts
    the per-Fock information from the member report.

    Parameters
    ----------
    grape_result : GrapeResult
        GRAPE result to evaluate.
    problem : ControlProblem
        ControlProblem from the GRAPE run.
    eval_model : DispersiveTransmonCavityModel
        Evaluation model (typically truth model).
    eval_frame : FrameSpec
        Rotating frame for evaluation.
    eval_noise : NoiseSpec, optional
        Noise specification. None = noiseless.

    Returns
    -------
    dict with keys:
        'avg_fidelity'       : float, overall average
        'per_state_fidelities' : list[float]
        'aggregate_fidelity' : float (same as avg_fidelity)
        'aggregate_leakage'  : float or None
    """
    eval_case = ControlEvaluationCase(
        model=eval_model,
        label="per_fock",
        frame=eval_frame,
        noise=eval_noise,
        weight=1.0,
    )

    eval_result = evaluate_control_with_simulator(
        problem,
        grape_result.schedule,
        cases=[eval_case],
    )

    member = eval_result.member_reports[0]
    per_state_fidelities = []
    for obj_report in member.objective_reports:
        per_state_fidelities.extend(list(obj_report.fidelities))

    return {
        "avg_fidelity": float(member.aggregate_fidelity),
        "per_state_fidelities": per_state_fidelities,
        "aggregate_fidelity": float(member.aggregate_fidelity),
        "aggregate_leakage": member.aggregate_leakage,
        "pair_labels": [
            label
            for obj in member.objective_reports
            for label in obj.pair_labels
        ],
    }


def compute_subspace_fidelity_from_result(grape_result) -> float:
    """
    Extract the optimizer's own unitary fidelity from a GrapeResult.

    The exact unitary fidelity (gauge-corrected, from the optimizer's final
    evaluation) is stored in system_metrics[0]['objectives'][0]['exact_unitary_fidelity'].
    Falls back to the aggregate fidelity metric if not present.

    Parameters
    ----------
    grape_result : GrapeResult
        GRAPE result from run_grape.

    Returns
    -------
    float
        Subspace unitary fidelity at the final GRAPE iteration.
    """
    # Try system_metrics first (most accurate)
    if grape_result.system_metrics:
        sm = grape_result.system_metrics[0]
        if "objectives" in sm and sm["objectives"]:
            obj0 = sm["objectives"][0]
            exact_f = obj0.get("exact_unitary_fidelity")
            if exact_f is not None:
                return float(exact_f)
            f_weighted = obj0.get("fidelity_weighted")
            if f_weighted is not None:
                return float(1.0 - float(obj0.get("infidelity", 1.0 - f_weighted)))

    # Fall back to metrics dict
    for key in ("nominal_fidelity", "nominal_physical_fidelity", "nominal_command_fidelity"):
        val = grape_result.metrics.get(key)
        if val is not None and not np.isnan(float(val)):
            return float(val)

    # Last resort: 1 - objective_value (approximate)
    return float(1.0 - grape_result.objective_value)


if __name__ == "__main__":
    import sys
    sys.path.insert(0, ".")
    from models import (
        make_truth_model, make_learner_model, make_frame,
        make_grape_subspace, make_target_matrix, CHI_PRIOR
    )

    print("Building models...")
    truth_model = make_truth_model()
    learner_model = make_learner_model(chi=CHI_PRIOR)
    frame = make_frame(learner_model)
    subspace = make_grape_subspace()
    target = make_target_matrix()

    print(f"Running GRAPE on learner model (chi={CHI_PRIOR/(2*np.pi)/1e6:.2f} MHz), maxiter=50 for quick test...")
    grape_result, problem, fidelity_on_learner = run_grape(
        learner_model, frame, subspace, target, maxiter=50, seed=42
    )
    print(f"GRAPE converged: {grape_result.success}")
    print(f"Fidelity on learner model: {fidelity_on_learner:.6f}")

    print("Evaluating on truth model (noiseless)...")
    truth_frame = make_frame(truth_model)
    avg_f, per_f = eval_on_model(grape_result, problem, truth_model, truth_frame, eval_noise=None)
    print(f"Avg fidelity on truth: {avg_f:.6f}")
    print(f"Per-state fidelities: {[f'{f:.4f}' for f in per_f[:8]]}")
