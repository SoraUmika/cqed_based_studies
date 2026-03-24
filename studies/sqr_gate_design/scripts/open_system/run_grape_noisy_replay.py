"""Solve representative GRAPE SQR controls and replay them under realistic noise."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.linalg import block_diag

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import runtime_compat  # noqa: F401

from common import (
    DATA_DIR,
    KAPPA_STORAGE_DEFAULT,
    N_CAV_TWO_MODE,
    N_STORAGE_LOGICAL,
    N_TR,
    PHI_TARGET,
    TARGET_STORAGE_LEVEL,
    THETA_TARGET,
    build_basis_initial_states,
    build_frame,
    build_multilevel_noise_spec,
    build_session,
    build_two_mode_model,
    duration_from_chi_t,
    target_qubit_unitary,
)

from cqed_sim import (
    ControlEvaluationCase,
    ControlProblem,
    GrapeConfig,
    GrapeSolver,
    ModelControlChannelSpec,
    PiecewiseConstantTimeGrid,
    StateTransferObjective,
    StateTransferPair,
    UnitaryObjective,
    build_control_problem_from_model,
    evaluate_control_with_simulator,
)
from cqed_sim.unitary_synthesis import Subspace

GRAPE_CHI_T = np.array([1.0, 2.0, 3.0, 5.0], dtype=float)
N_SLICES = 48
AMP_BOUND = 2.0 * np.pi * 50.0e6
MAXITER = 200
REPLAY_SUBSTEPS_PER_SLICE = 512
CONTROL_SHAPE = (2, N_SLICES)
OUTPUT_PATH = DATA_DIR / "grape_noisy_replay.npz"
CHECKPOINT_PATH = DATA_DIR / "grape_noisy_replay.partial.npz"


def replay_dt_for_duration(
    duration: float,
    *,
    n_slices: int = N_SLICES,
    substeps_per_slice: int = REPLAY_SUBSTEPS_PER_SLICE,
) -> float:
    return float(duration) / (int(n_slices) * int(substeps_per_slice))


def build_grape_problem(duration: float, *, n_slices: int = N_SLICES):
    model = build_two_mode_model()
    frame = build_frame(model)
    target, subspace, phase_blocks = build_target_subspace()
    problem = build_control_problem_from_model(
        model,
        frame=frame,
        time_grid=PiecewiseConstantTimeGrid.uniform(steps=int(n_slices), dt_s=float(duration) / int(n_slices)),
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
        objectives=(
            UnitaryObjective(
                target_operator=target,
                subspace=subspace,
                ignore_global_phase=True,
                phase_blocks=phase_blocks,
            ),
        ),
    )
    return model, frame, problem


def command_values_from_result(result) -> np.ndarray:
    if result.command_values is not None:
        return np.asarray(result.command_values, dtype=float)
    return np.asarray(result.schedule.command_values(), dtype=float)


def representative_noise_spec():
    return build_multilevel_noise_spec(
        transmon_t1=(30.0e-6, 10.0e-6),
        kappa_storage=KAPPA_STORAGE_DEFAULT,
        nth_storage=0.02,
        for_three_mode=False,
    )


def build_replay_case(
    model,
    frame,
    duration: float,
    *,
    replay_substeps_per_slice: int = REPLAY_SUBSTEPS_PER_SLICE,
) -> ControlEvaluationCase:
    replay_dt = replay_dt_for_duration(
        duration,
        n_slices=N_SLICES,
        substeps_per_slice=replay_substeps_per_slice,
    )
    return ControlEvaluationCase(
        model=model,
        label="representative_noisy",
        frame=frame,
        noise=representative_noise_spec(),
        compiler_dt_s=replay_dt,
        max_step_s=replay_dt,
        metadata={
            "replay_substeps_per_slice": int(replay_substeps_per_slice),
        },
    )


def evaluate_schedule_target_fidelity(
    problem,
    model,
    frame,
    schedule,
    duration: float,
    *,
    replay_substeps_per_slice: int = REPLAY_SUBSTEPS_PER_SLICE,
) -> float:
    evaluation = evaluate_control_with_simulator(
        problem,
        schedule,
        cases=(build_replay_case(model, frame, duration, replay_substeps_per_slice=replay_substeps_per_slice),),
        waveform_mode="command",
    )
    return float(evaluation.metrics["aggregate_fidelity"])


def _ideal_target_states(problem, result, model) -> tuple[object, ...]:
    _, initial_states = build_basis_initial_states(model, n_storage_levels=N_STORAGE_LOGICAL)
    nominal_unitary = None if result.nominal_final_unitary is None else np.asarray(result.nominal_final_unitary, dtype=np.complex128)
    if nominal_unitary is None:
        replay_dt = replay_dt_for_duration(float(problem.time_grid.duration_s))
        pulses, drive_ops, _ = result.to_pulses()
        session = build_session(model, build_frame(model), list(pulses), drive_ops, duration=problem.time_grid.duration_s, noise=None, dt=replay_dt)
        return tuple(session.run(state).final_state for state in initial_states)
    return tuple(
        nominal_unitary @ np.asarray(initial_state.full(), dtype=np.complex128).reshape(-1)
        for initial_state in initial_states
    )


def build_ideal_replay_problem(problem, result, model) -> ControlProblem:
    labels, initial_states = build_basis_initial_states(model, n_storage_levels=N_STORAGE_LOGICAL)
    target_states = _ideal_target_states(problem, result, model)
    objective = StateTransferObjective(
        pairs=tuple(
            StateTransferPair(
                initial_state=initial_state,
                target_state=target_state,
                label=f"q{int(qubit_level)}_s{int(storage_level)}",
            )
            for (qubit_level, storage_level), initial_state, target_state in zip(labels, initial_states, target_states, strict=True)
        ),
        name="ideal_replay_match",
    )
    return ControlProblem(
        parameterization=problem.parameterization,
        systems=problem.systems,
        objectives=(objective,),
        penalties=problem.penalties,
        ensemble_aggregate=problem.ensemble_aggregate,
        hardware_model=problem.hardware_model,
        metadata=dict(problem.metadata),
    )


def noisy_replay_ideal_fidelity(
    problem,
    model,
    frame,
    result,
    schedule,
    duration: float,
    *,
    replay_substeps_per_slice: int = REPLAY_SUBSTEPS_PER_SLICE,
) -> float:
    ideal_problem = build_ideal_replay_problem(problem, result, model)
    evaluation = evaluate_control_with_simulator(
        ideal_problem,
        schedule,
        cases=(build_replay_case(model, frame, duration, replay_substeps_per_slice=replay_substeps_per_slice),),
        waveform_mode="command",
    )
    return float(evaluation.metrics["aggregate_fidelity"])


def replay_archived_command_values(
    chi_t_value: float,
    command_values: np.ndarray,
    *,
    n_slices: int = N_SLICES,
    replay_substeps_per_slice: int = REPLAY_SUBSTEPS_PER_SLICE,
) -> float:
    duration = duration_from_chi_t(float(chi_t_value))
    model, frame, problem = build_grape_problem(duration, n_slices=n_slices)
    return evaluate_schedule_target_fidelity(
        problem,
        model,
        frame,
        np.asarray(command_values, dtype=float),
        duration,
        replay_substeps_per_slice=replay_substeps_per_slice,
    )


def save_checkpoint(
    noisy_fidelity_to_ideal: np.ndarray,
    noisy_fidelity_to_target: np.ndarray,
    objective_fidelity: np.ndarray,
    converged: np.ndarray,
    done_mask: np.ndarray,
    command_values: np.ndarray,
) -> None:
    np.savez(
        CHECKPOINT_PATH,
        chi_t_values=GRAPE_CHI_T,
        objective_fidelity=objective_fidelity,
        noisy_fidelity_to_ideal=noisy_fidelity_to_ideal,
        noisy_fidelity_to_target=noisy_fidelity_to_target,
        converged=converged,
        done_mask=done_mask,
        command_values=command_values,
        n_slices=N_SLICES,
        amp_bound=AMP_BOUND,
        maxiter=MAXITER,
        replay_substeps_per_slice=REPLAY_SUBSTEPS_PER_SLICE,
        replay_dt_s=np.array([replay_dt_for_duration(duration_from_chi_t(float(value))) for value in GRAPE_CHI_T], dtype=float),
        kappa_storage_rad_s=KAPPA_STORAGE_DEFAULT,
    )


def load_checkpoint() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    noisy_fidelity_to_ideal = np.zeros(len(GRAPE_CHI_T), dtype=float)
    noisy_fidelity_to_target = np.zeros(len(GRAPE_CHI_T), dtype=float)
    objective_fidelity = np.zeros(len(GRAPE_CHI_T), dtype=float)
    converged = np.zeros(len(GRAPE_CHI_T), dtype=bool)
    done_mask = np.zeros(len(GRAPE_CHI_T), dtype=bool)
    command_values = np.zeros((len(GRAPE_CHI_T),) + CONTROL_SHAPE, dtype=float)
    if not CHECKPOINT_PATH.exists():
        return noisy_fidelity_to_ideal, noisy_fidelity_to_target, objective_fidelity, converged, done_mask, command_values

    payload = np.load(CHECKPOINT_PATH, allow_pickle=True)
    if (
        tuple(payload["objective_fidelity"].shape) != (len(GRAPE_CHI_T),)
        or "command_values" not in payload.files
        or tuple(payload["command_values"].shape) != ((len(GRAPE_CHI_T),) + CONTROL_SHAPE)
    ):
        return noisy_fidelity_to_ideal, noisy_fidelity_to_target, objective_fidelity, converged, done_mask, command_values

    noisy_fidelity_to_ideal[...] = payload["noisy_fidelity_to_ideal"]
    noisy_fidelity_to_target[...] = payload["noisy_fidelity_to_target"]
    objective_fidelity[...] = payload["objective_fidelity"]
    converged[...] = payload["converged"]
    done_mask[...] = payload["done_mask"]
    command_values[...] = payload["command_values"]
    return noisy_fidelity_to_ideal, noisy_fidelity_to_target, objective_fidelity, converged, done_mask, command_values


def build_target_subspace() -> tuple[np.ndarray, Subspace, tuple[tuple[int, int], ...]]:
    """Return the SQR target operator, the working subspace, and phase blocks."""
    indices: list[int] = []
    labels: list[str] = []
    for storage_level in range(N_STORAGE_LOGICAL):
        indices.append(0 * N_CAV_TWO_MODE + storage_level)
        indices.append(1 * N_CAV_TWO_MODE + storage_level)
        labels.append(f"|g,{storage_level}>")
        labels.append(f"|e,{storage_level}>")
    subspace = Subspace(
        full_dim=N_TR * N_CAV_TWO_MODE,
        indices=tuple(indices),
        labels=tuple(labels),
    )
    identity_2 = np.eye(2, dtype=np.complex128)
    target_block = target_qubit_unitary(THETA_TARGET, PHI_TARGET)
    target = block_diag(
        *[
            target_block if storage_level == TARGET_STORAGE_LEVEL else identity_2
            for storage_level in range(N_STORAGE_LOGICAL)
        ]
    )
    phase_blocks = tuple((2 * storage_level, 2 * storage_level + 1) for storage_level in range(N_STORAGE_LOGICAL))
    return target, subspace, phase_blocks


def solve_grape_for_duration(duration: float, *, n_slices: int = N_SLICES, maxiter: int = MAXITER):
    """Solve the cphase-SQR GRAPE problem for one gate duration."""
    model, frame, problem = build_grape_problem(duration, n_slices=n_slices)
    result = GrapeSolver(GrapeConfig(maxiter=int(maxiter), seed=17)).solve(problem)
    return model, frame, problem, result


def extract_reported_fidelity(result) -> float:
    """Return the main fidelity metric from a GRAPE result."""
    if "nominal_fidelity" in result.metrics:
        return float(result.metrics["nominal_fidelity"])
    if "fidelity" in result.metrics:
        return float(result.metrics["fidelity"])
    objective_value = float(result.objective_value)
    return float(1.0 - objective_value) if objective_value <= 1.0 else 0.0


def main() -> None:
    (
        noisy_fidelity_to_ideal,
        noisy_fidelity_to_target,
        objective_fidelity,
        converged,
        done_mask,
        command_values,
    ) = load_checkpoint()
    if CHECKPOINT_PATH.exists():
        print(f"Resuming from {CHECKPOINT_PATH}")

    print("=" * 68)
    print("A2 GRAPE noisy replay")
    print("=" * 68)

    for idx, chi_t_value in enumerate(GRAPE_CHI_T):
        if done_mask[idx]:
            print(f"chiT/2pi={chi_t_value:3.1f}  resumed")
            continue
        duration = duration_from_chi_t(float(chi_t_value))
        model, frame, problem, result = solve_grape_for_duration(duration)
        command_values[idx] = command_values_from_result(result)

        archived_noisy_target_fidelity = replay_archived_command_values(
            float(chi_t_value),
            command_values[idx],
            n_slices=N_SLICES,
            replay_substeps_per_slice=REPLAY_SUBSTEPS_PER_SLICE,
        )
        archived_noisy_ideal_fidelity = noisy_replay_ideal_fidelity(
            problem,
            model,
            frame,
            result,
            command_values[idx],
            duration,
            replay_substeps_per_slice=REPLAY_SUBSTEPS_PER_SLICE,
        )

        noisy_fidelity_to_ideal[idx] = archived_noisy_ideal_fidelity
        noisy_fidelity_to_target[idx] = archived_noisy_target_fidelity
        objective_fidelity[idx] = extract_reported_fidelity(result)
        converged[idx] = bool(result.success)
        done_mask[idx] = True
        save_checkpoint(
            noisy_fidelity_to_ideal,
            noisy_fidelity_to_target,
            objective_fidelity,
            converged,
            done_mask,
            command_values,
        )

        print(
            f"chiT/2pi={chi_t_value:3.1f}  objective_F={objective_fidelity[idx]:.6f}  "
            f"noisy_to_ideal={noisy_fidelity_to_ideal[idx]:.6f}  "
            f"noisy_to_target={noisy_fidelity_to_target[idx]:.6f}  "
            f"success={converged[idx]}"
        )

    np.savez(
        OUTPUT_PATH,
        chi_t_values=GRAPE_CHI_T,
        objective_fidelity=objective_fidelity,
        noisy_fidelity_to_ideal=noisy_fidelity_to_ideal,
        noisy_fidelity_to_target=noisy_fidelity_to_target,
        converged=converged,
        command_values=command_values,
        n_slices=N_SLICES,
        amp_bound=AMP_BOUND,
        maxiter=MAXITER,
        replay_substeps_per_slice=REPLAY_SUBSTEPS_PER_SLICE,
        replay_dt_s=np.array([replay_dt_for_duration(duration_from_chi_t(float(value))) for value in GRAPE_CHI_T], dtype=float),
        kappa_storage_rad_s=KAPPA_STORAGE_DEFAULT,
    )
    CHECKPOINT_PATH.unlink(missing_ok=True)
    print(f"Saved {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
