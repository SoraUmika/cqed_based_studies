"""Shared helpers for the simultaneous multitone SQR study."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

import runtime_compat  # noqa: F401


SCRIPT_DIR = Path(__file__).resolve().parent
STUDY_DIR = SCRIPT_DIR.parent
DATA_DIR = STUDY_DIR / "data"
FIGURES_DIR = STUDY_DIR / "figures"
REPORT_DIR = STUDY_DIR / "report"
REPO_ROOT = STUDY_DIR.parent.parent
USER_ROOT = REPO_ROOT.parent
CQED_SIM_ROOT = USER_ROOT / "cQED_simulation"


def ensure_cqed_sim_on_path() -> None:
    """Add the local cqed_sim checkout to sys.path when needed."""
    candidate = str(CQED_SIM_ROOT)
    if candidate not in sys.path:
        sys.path.insert(0, candidate)


ensure_cqed_sim_on_path()

from cqed_sim.calibration import (  # noqa: E402
    ConditionedMultitoneRunConfig,
    ConditionedQubitTargets,
    TargetedSubspaceObjectiveWeights,
)
from cqed_sim.core.frame import FrameSpec  # noqa: E402
from cqed_sim.core.frequencies import manifold_transition_frequency  # noqa: E402
from cqed_sim.core.model import DispersiveTransmonCavityModel  # noqa: E402


OMEGA_Q = 2.0 * np.pi * 6.150e9
OMEGA_C = 2.0 * np.pi * 5.241e9
ALPHA = 2.0 * np.pi * (-255.0e6)
CHI = 2.0 * np.pi * (-2.84e6)

LOGICAL_LEVELS = (0, 1, 2, 3)
N_CAV = 6
N_TR_REDUCED = 2
N_TR_QUTRIT = 3

DT_S = 4.0e-9
SIGMA_FRACTION = 1.0 / 6.0
CHI_T_VALUES = np.asarray([1.0, 2.0, 3.0, 5.0], dtype=float)
THETA_VALUES = np.pi * np.asarray([0.125, 0.25, 0.5, 1.0], dtype=float)

CASE_SPECS: dict[str, tuple[int, ...]] = {
    "pair_adjacent": (0, 1),
    "pair_separated": (0, 2),
    "triple_low": (0, 1, 2),
    "all_four": (0, 1, 2, 3),
}

CASE_LABELS: dict[str, str] = {
    "pair_adjacent": "Pair (0,1)",
    "pair_separated": "Pair (0,2)",
    "triple_low": "Triple (0,1,2)",
    "all_four": "All four",
}

OBJECTIVE_WEIGHTS = TargetedSubspaceObjectiveWeights(
    qubit_weight=0.3,
    subspace_weight=1.0,
    preservation_weight=0.25,
    leakage_weight=0.25,
)


def duration_from_chi_t(chi_t_2pi: float, chi: float = CHI) -> float:
    """Return the gate duration associated with the chosen chi-period count."""
    f_chi_hz = abs(float(chi)) / (2.0 * np.pi)
    return float(chi_t_2pi) / f_chi_hz


def build_model(*, n_tr: int = N_TR_REDUCED, n_cav: int = N_CAV) -> DispersiveTransmonCavityModel:
    """Build the dispersive transmon-cavity model used throughout the study."""
    return DispersiveTransmonCavityModel(
        omega_c=OMEGA_C,
        omega_q=OMEGA_Q,
        alpha=ALPHA,
        chi=CHI,
        n_cav=int(n_cav),
        n_tr=int(n_tr),
    )


def build_frame(model: DispersiveTransmonCavityModel) -> FrameSpec:
    """Use a frame rotating at the bare qubit and cavity frequencies."""
    return FrameSpec(omega_c_frame=float(model.omega_c), omega_q_frame=float(model.omega_q))


def run_config_for_chi_t(
    model: DispersiveTransmonCavityModel,
    chi_t_2pi: float,
    *,
    dt_s: float = DT_S,
) -> ConditionedMultitoneRunConfig:
    """Build the standard multitone run configuration for one duration."""
    frame = build_frame(model)
    fock_fqs_hz = tuple(
        float(manifold_transition_frequency(model, n, frame=frame) / (2.0 * np.pi))
        for n in LOGICAL_LEVELS
    )
    return ConditionedMultitoneRunConfig(
        frame=frame,
        duration_s=duration_from_chi_t(float(chi_t_2pi)),
        dt_s=float(dt_s),
        sigma_fraction=float(SIGMA_FRACTION),
        tone_cutoff=1.0e-12,
        include_all_levels=False,
        max_step_s=float(dt_s),
        fock_fqs_hz=fock_fqs_hz,
    )


def build_targets(target_levels: tuple[int, ...], theta: float) -> ConditionedQubitTargets:
    """Build a logical target with X-axis rotations on selected branches."""
    target_set = {int(level) for level in target_levels}
    spec = {
        int(level): (float(theta) if int(level) in target_set else 0.0, 0.0)
        for level in LOGICAL_LEVELS
    }
    return ConditionedQubitTargets.from_spec(spec, n_levels=len(LOGICAL_LEVELS))


def process_fidelity(target_operator: np.ndarray, actual_operator: np.ndarray) -> float:
    """Return the standard restricted process fidelity."""
    target = np.asarray(target_operator, dtype=np.complex128)
    actual = np.asarray(actual_operator, dtype=np.complex128)
    dim = float(target.shape[0])
    overlap = np.trace(target.conj().T @ actual)
    return float(np.clip(abs(overlap) ** 2 / (dim * dim), 0.0, 1.0))


def target_and_spectator_means(
    sector_metrics,
    target_levels: tuple[int, ...],
) -> tuple[float, float]:
    """Average conditioned fidelities over target and spectator branches."""
    target_set = {int(level) for level in target_levels}
    target_values = [float(metric.fidelity) for metric in sector_metrics if int(metric.n) in target_set]
    spectator_values = [float(metric.fidelity) for metric in sector_metrics if int(metric.n) not in target_set]
    target_mean = float(np.mean(target_values)) if target_values else float("nan")
    spectator_mean = float(np.mean(spectator_values)) if spectator_values else float("nan")
    return target_mean, spectator_mean


def logical_indices(model: DispersiveTransmonCavityModel, logical_levels: tuple[int, ...] = LOGICAL_LEVELS) -> tuple[int, ...]:
    """Indices of the selected logical qubit-cavity basis states in the full Hilbert space."""
    indices: list[int] = []
    for level in logical_levels:
        indices.extend([int(level), int(model.n_cav) + int(level)])
    return tuple(indices)
