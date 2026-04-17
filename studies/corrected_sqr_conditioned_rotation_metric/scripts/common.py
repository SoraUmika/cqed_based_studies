"""Shared study helpers for the corrected SQR conditioned-rotation study."""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

import runtime_compat  # noqa: F401


SCRIPT_DIR = Path(__file__).resolve().parent
STUDY_DIR = SCRIPT_DIR.parent
DATA_DIR = STUDY_DIR / "data"
FIGURES_DIR = STUDY_DIR / "figures"
ARTIFACTS_DIR = STUDY_DIR / "artifacts"
REPORT_DIR = STUDY_DIR / "report"
REPO_ROOT = STUDY_DIR.parent.parent
USER_ROOT = REPO_ROOT.parent
CQED_SIM_ROOT = USER_ROOT / "cQED_simulation"


def ensure_cqed_sim_on_path() -> None:
    candidate = str(CQED_SIM_ROOT)
    if candidate not in sys.path:
        sys.path.insert(0, candidate)


ensure_cqed_sim_on_path()

from cqed_sim.calibration import (  # noqa: E402
    ConditionedMultitoneCorrections,
    ConditionedMultitoneRunConfig,
    ConditionedOptimizationConfig,
    ConditionedQubitTargets,
)
from cqed_sim.core.frame import FrameSpec  # noqa: E402
from cqed_sim.core.frequencies import manifold_transition_frequency  # noqa: E402
from cqed_sim.core.model import DispersiveTransmonCavityModel  # noqa: E402


OMEGA_Q = 2.0 * np.pi * 6.150e9
OMEGA_C = 2.0 * np.pi * 5.241e9
ALPHA = 2.0 * np.pi * (-255.0e6)
CHI = 2.0 * np.pi * (-2.84e6)

N_CAV = 8
N_TR = 2
DT_S = 4.0e-9
SIGMA_FRACTION = 1.0 / 6.0
ACTIVE_WINDOWS = (1, 2, 3, 4)
CHI_T_VALUES = (1.0, 3.0, 5.0)

# Smooth, nontrivial corrected-SQR profile used throughout the study.
THETA_PROFILE_RAD = np.pi * np.asarray([0.20, 0.35, 0.50, 0.65], dtype=float)
PHI_PROFILE_RAD = np.pi * np.asarray([0.00, 0.25, 0.50, 0.75], dtype=float)


def wrap_pi(value: float) -> float:
    return float((float(value) + np.pi) % (2.0 * np.pi) - np.pi)


def duration_from_chi_t(chi_t_over_2pi: float, chi_rad_s: float = CHI) -> float:
    f_chi_hz = abs(float(chi_rad_s)) / (2.0 * np.pi)
    return float(chi_t_over_2pi) / f_chi_hz


def build_model(*, n_cav: int = N_CAV, n_tr: int = N_TR) -> DispersiveTransmonCavityModel:
    return DispersiveTransmonCavityModel(
        omega_c=OMEGA_C,
        omega_q=OMEGA_Q,
        alpha=ALPHA,
        chi=CHI,
        n_cav=int(n_cav),
        n_tr=int(n_tr),
    )


def build_frame(model: DispersiveTransmonCavityModel) -> FrameSpec:
    return FrameSpec(omega_c_frame=float(model.omega_c), omega_q_frame=float(model.omega_q))


def build_run_config(
    model: DispersiveTransmonCavityModel,
    *,
    duration_s: float,
    dt_s: float = DT_S,
    sigma_fraction: float = SIGMA_FRACTION,
) -> ConditionedMultitoneRunConfig:
    frame = build_frame(model)
    return ConditionedMultitoneRunConfig(
        frame=frame,
        duration_s=float(duration_s),
        dt_s=float(dt_s),
        sigma_fraction=float(sigma_fraction),
        tone_cutoff=1.0e-12,
        include_all_levels=False,
        max_step_s=float(dt_s),
        # Let cqed_sim compute the in-frame manifold frequencies internally.
        # The optional fock_fqs_hz override is interpreted as an absolute qubit
        # transition-frequency list, not as an already frame-shifted detuning list.
        fock_fqs_hz=None,
    )


def run_config_for_chi_t(
    model: DispersiveTransmonCavityModel,
    chi_t_over_2pi: float,
    *,
    dt_s: float = DT_S,
) -> ConditionedMultitoneRunConfig:
    return build_run_config(
        model,
        duration_s=duration_from_chi_t(float(chi_t_over_2pi)),
        dt_s=float(dt_s),
    )


def build_targets(n_active: int) -> ConditionedQubitTargets:
    count = int(n_active)
    theta = tuple(float(value) for value in THETA_PROFILE_RAD[:count])
    phi = tuple(float(value) for value in PHI_PROFILE_RAD[:count])
    weights = tuple(float(1.0 / count) for _ in range(count))
    return ConditionedQubitTargets(theta=theta, phi=phi, weights=weights)


def target_rows(n_active: int) -> list[dict[str, float]]:
    targets = build_targets(n_active)
    return targets.as_rows()


def as_serializable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(key): as_serializable(inner) for key, inner in value.items()}
    if isinstance(value, (list, tuple)):
        return [as_serializable(inner) for inner in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def save_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(as_serializable(payload), indent=2), encoding="utf-8")
