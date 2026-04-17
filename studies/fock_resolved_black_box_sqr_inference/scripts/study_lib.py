"""Inference and case-generation utilities for the Fock-resolved SQR study."""

from __future__ import annotations

import math
import time
from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np
import qutip as qt
from scipy.optimize import minimize
from scipy.special import logsumexp

import runtime_compat  # noqa: F401

from cqed_sim.calibration.conditioned_multitone import (
    ConditionedOptimizationConfig,
    optimize_conditioned_multitone,
)
from cqed_sim.core import displacement_op
from cqed_sim.core.ideal_gates import qubit_rotation_axis, qubit_rotation_xy

from common import (
    CHI,
    CHI_PRIME,
    DEFAULT_CPSQR_PHASES,
    DEFAULT_PHI,
    DEFAULT_POPULATIONS,
    DEFAULT_THETA,
    IMPERFECT_DURATION_S,
    N_ACTIVE,
    NEAR_IDEAL_DURATION_S,
    PAULI_Z,
    SectorSummary,
    as_dm,
    bloch_from_density,
    build_frame,
    build_model,
    build_multitone_waveform,
    cavity_mixture_state,
    cavity_superposition_state,
    coherent_state,
    compile_waveform,
    conditioned_multitone_targets,
    conditioned_run_config,
    cpsqr_like_operator,
    density_matrix_fidelity,
    displacement_probability_matrix,
    ideal_sqr_operator,
    json_ready,
    normalize_probabilities,
    qubit_density_from_bloch,
    qubit_excited_dm,
    qubit_ground_dm,
    qubit_plus_dm,
    qubit_plus_y_dm,
    sector_summaries_from_state,
    simulate_compiled_on_states,
    trace_distance,
)


MEASUREMENT_AXES = ("X", "Y", "Z")
DEFAULT_SHOT_GRID = (100, 300, 1_000, 3_000, 10_000)
DEFAULT_WAIT_GRID_S = tuple(float(x) for x in np.linspace(0.0, 0.35e-6, 11))
DEFAULT_COMBINED_WAIT_GRID_S = tuple(float(x) for x in np.linspace(0.0, 0.35e-6, 9))
DEFAULT_ALPHA_GRID = (
    0.0 + 0.0j,
    0.35 + 0.0j,
    0.55 + 0.0j,
    0.35j,
    0.55j,
    0.35 + 0.35j,
    0.55 + 0.25j,
)
DEFAULT_COMBINED_ALPHA_GRID = (
    0.0 + 0.0j,
    0.35 + 0.0j,
    0.55 + 0.0j,
    0.35j,
    0.35 + 0.35j,
)
KERNEL_DIM = 10
EPS = 1.0e-12


@dataclass(frozen=True)
class MeasurementSetting:
    protocol: str
    alpha: complex
    wait_s: float
    label: str

    def as_dict(self) -> dict[str, Any]:
        return {
            "protocol": str(self.protocol),
            "alpha_real": float(np.real(self.alpha)),
            "alpha_imag": float(np.imag(self.alpha)),
            "wait_s": float(self.wait_s),
            "label": str(self.label),
        }


@dataclass(frozen=True)
class NoiseModel:
    shots: int
    rotation_sigma_rad: float = 0.0
    displacement_rel_sigma: float = 0.0
    displacement_phase_sigma_rad: float = 0.0
    chi_rel_sigma: float = 0.0
    chi_prime_rel_sigma: float = 0.0
    t1_s: float | None = None
    t2_s: float | None = None
    sample_mode: str = "binomial"

    def as_dict(self) -> dict[str, Any]:
        return {
            "shots": int(self.shots),
            "rotation_sigma_rad": float(self.rotation_sigma_rad),
            "displacement_rel_sigma": float(self.displacement_rel_sigma),
            "displacement_phase_sigma_rad": float(self.displacement_phase_sigma_rad),
            "chi_rel_sigma": float(self.chi_rel_sigma),
            "chi_prime_rel_sigma": float(self.chi_prime_rel_sigma),
            "t1_s": None if self.t1_s is None else float(self.t1_s),
            "t2_s": None if self.t2_s is None else float(self.t2_s),
            "sample_mode": str(self.sample_mode),
        }


@dataclass(frozen=True)
class StudyCase:
    case_id: str
    family: str
    description: str
    model_class: str
    state: qt.Qobj
    model: Any
    frame: Any
    truth_sectors: tuple[SectorSummary, ...]
    metadata: dict[str, Any]

    @property
    def weighted_transverse_truth(self) -> np.ndarray:
        return np.asarray(
            [summary.population * (summary.x + 1.0j * summary.y) for summary in self.truth_sectors],
            dtype=np.complex128,
        )

    @property
    def weighted_longitudinal_truth(self) -> np.ndarray:
        return np.asarray(
            [summary.population * summary.z for summary in self.truth_sectors],
            dtype=float,
        )

    @property
    def z_total_truth(self) -> float:
        return float(np.sum(self.weighted_longitudinal_truth))

    @property
    def populations_truth(self) -> np.ndarray:
        return np.asarray([summary.population for summary in self.truth_sectors], dtype=float)


def complex_to_dict(value: complex) -> dict[str, float]:
    return {"real": float(np.real(value)), "imag": float(np.imag(value))}


def qubit_input_state(label: str) -> qt.Qobj:
    key = str(label).lower()
    if key == "g":
        return qubit_ground_dm()
    if key == "e":
        return qubit_excited_dm()
    if key in {"+", "plus", "x+"}:
        return qubit_plus_dm()
    if key in {"+y", "plus_y", "y+"}:
        return qubit_plus_y_dm()
    raise ValueError(f"Unsupported qubit input label: {label}")


def embed_qubit_operator(op_2x2: qt.Qobj, n_tr: int) -> qt.Qobj:
    n_tr = int(n_tr)
    mat = np.eye(n_tr, dtype=np.complex128)
    mat[:2, :2] = np.asarray(op_2x2.full(), dtype=np.complex128)
    return qt.Qobj(mat, dims=[[n_tr], [n_tr]])


def embed_qubit_density(rho_q: qt.Qobj, n_tr: int) -> qt.Qobj:
    return embed_qubit_operator(rho_q, n_tr)


def embedded_sigma_z(n_tr: int) -> qt.Qobj:
    diag = np.zeros(int(n_tr), dtype=np.complex128)
    diag[0] = 1.0
    if int(n_tr) > 1:
        diag[1] = -1.0
    return qt.Qobj(np.diag(diag), dims=[[int(n_tr)], [int(n_tr)]])


def tomography_prerotation(axis: str, angle_error: float = 0.0) -> qt.Qobj:
    key = str(axis).upper()
    delta = float(angle_error)
    if key == "Z":
        return qt.qeye(2)
    if key == "X":
        return qubit_rotation_xy(-0.5 * np.pi + delta, 0.5 * np.pi)
    if key == "Y":
        return qubit_rotation_xy(0.5 * np.pi + delta, 0.0)
    raise ValueError(f"Unsupported tomography axis: {axis}")


def branch_expectations_from_density(rho_q: qt.Qobj, *, angle_errors: dict[str, float] | None = None) -> dict[str, float]:
    errors = {} if angle_errors is None else {str(key).upper(): float(value) for key, value in angle_errors.items()}
    out: dict[str, float] = {}
    for axis in MEASUREMENT_AXES:
        pre = tomography_prerotation(axis, errors.get(axis, 0.0))
        out[axis] = float(np.real(qt.expect(PAULI_Z, pre * rho_q * pre.dag())))
    return out


def branch_expectations_from_joint_state(
    state: qt.Qobj,
    *,
    n_tr: int,
    angle_errors: dict[str, float] | None = None,
) -> dict[str, float]:
    errors = {} if angle_errors is None else {str(key).upper(): float(value) for key, value in angle_errors.items()}
    reduced = state.ptrace(0)
    sigma_z = embedded_sigma_z(n_tr)
    out: dict[str, float] = {}
    for axis in MEASUREMENT_AXES:
        pre = embed_qubit_operator(tomography_prerotation(axis, errors.get(axis, 0.0)), n_tr=n_tr)
        out[axis] = float(np.real(qt.expect(sigma_z, pre * reduced * pre.dag())))
    return out


def project_density_to_physical(matrix: np.ndarray) -> np.ndarray:
    herm = 0.5 * (matrix + matrix.conj().T)
    evals, evecs = np.linalg.eigh(herm)
    evals = np.clip(np.real(evals), 0.0, None)
    if float(np.sum(evals)) <= EPS:
        return np.asarray(qubit_ground_dm().full(), dtype=np.complex128)
    clipped = (evecs * evals) @ evecs.conj().T
    return clipped / float(np.real(np.trace(clipped)))


def cholesky_raw_to_density(raw: Sequence[float]) -> qt.Qobj:
    values = np.asarray(raw, dtype=float).reshape(4)
    t1 = math.exp(float(np.clip(values[0], -12.0, 12.0)))
    t4 = math.exp(float(np.clip(values[3], -12.0, 12.0)))
    t2 = complex(float(values[1]), float(values[2]))
    tri = np.asarray([[t1, 0.0], [t2, t4]], dtype=np.complex128)
    rho = tri.conj().T @ tri
    rho = rho / float(np.real(np.trace(rho)))
    return qt.Qobj(rho, dims=[[2], [2]])


def density_to_cholesky_raw(rho_q: qt.Qobj) -> np.ndarray:
    rho = project_density_to_physical(np.asarray(rho_q.full(), dtype=np.complex128))
    rho11 = float(max(np.real(rho[1, 1]), EPS))
    t4 = math.sqrt(rho11)
    if t4 <= math.sqrt(EPS):
        t2 = 0.0 + 0.0j
        t1 = math.sqrt(float(max(np.real(rho[0, 0]), EPS)))
    else:
        t2 = np.conj(rho[0, 1]) / t4
        t1_sq = float(max(np.real(rho[0, 0] - abs(t2) ** 2), EPS))
        t1 = math.sqrt(t1_sq)
    return np.asarray([math.log(t1), float(np.real(t2)), float(np.imag(t2)), math.log(t4)], dtype=float)


def linear_inversion_density(expectations: dict[str, float]) -> qt.Qobj:
    rho = qubit_density_from_bloch(
        float(expectations["X"]),
        float(expectations["Y"]),
        float(expectations["Z"]),
    )
    clipped = project_density_to_physical(np.asarray(rho.full(), dtype=np.complex128))
    return qt.Qobj(clipped, dims=[[2], [2]])


def single_qubit_target_state() -> qt.Qobj:
    u_target = qubit_rotation_axis(np.pi / 3.0, "x") * qubit_rotation_axis(np.pi / 4.0, "z")
    return u_target * qubit_ground_dm() * u_target.dag()


def sample_branch_measurement(
    expectation: float,
    *,
    shots: int,
    rng: np.random.Generator,
    mode: str,
) -> tuple[float, int | None]:
    value = float(np.clip(expectation, -1.0, 1.0))
    if str(mode) == "gaussian":
        sigma = 1.0 / math.sqrt(max(int(shots), 1))
        return float(np.clip(rng.normal(value, sigma), -1.0, 1.0)), None
    p_plus = float(np.clip(0.5 * (1.0 + value), EPS, 1.0 - EPS))
    counts_plus = int(rng.binomial(int(shots), p_plus))
    measured = (2.0 * counts_plus / float(shots)) - 1.0
    return float(measured), counts_plus


def single_qubit_dataset(
    rho_true: qt.Qobj,
    *,
    shots: int,
    rng: np.random.Generator,
    sample_mode: str = "binomial",
) -> dict[str, Any]:
    exact = branch_expectations_from_density(rho_true)
    measured: dict[str, float] = {}
    counts_plus: dict[str, int | None] = {}
    for axis in MEASUREMENT_AXES:
        measured_value, counts = sample_branch_measurement(
            exact[axis],
            shots=int(shots),
            rng=rng,
            mode=sample_mode,
        )
        measured[axis] = float(measured_value)
        counts_plus[axis] = counts
    return {
        "shots": int(shots),
        "exact_expectations": exact,
        "measured_expectations": measured,
        "counts_plus": counts_plus,
        "sample_mode": str(sample_mode),
    }


def fit_single_qubit_ls(dataset: dict[str, Any]) -> dict[str, Any]:
    measured = {axis: float(dataset["measured_expectations"][axis]) for axis in MEASUREMENT_AXES}
    x0 = density_to_cholesky_raw(linear_inversion_density(measured))

    def objective(raw: np.ndarray) -> float:
        pred = branch_expectations_from_density(cholesky_raw_to_density(raw))
        return float(sum((pred[axis] - measured[axis]) ** 2 for axis in MEASUREMENT_AXES))

    bounds = [(-8.0, 8.0), (-2.5, 2.5), (-2.5, 2.5), (-8.0, 8.0)]
    result = minimize(objective, x0, method="L-BFGS-B", bounds=bounds)
    if not result.success:
        result = minimize(objective, x0, method="Powell", options={"maxiter": 400})
    rho_opt = cholesky_raw_to_density(result.x)
    return {
        "rho": rho_opt,
        "objective": float(objective(result.x)),
        "expectations": branch_expectations_from_density(rho_opt),
        "optimizer_success": bool(result.success),
        "optimizer_message": str(result.message),
        "raw_params": np.asarray(result.x, dtype=float),
    }


def fit_single_qubit_mle(dataset: dict[str, Any]) -> dict[str, Any]:
    measured = {axis: float(dataset["measured_expectations"][axis]) for axis in MEASUREMENT_AXES}
    counts_plus = dataset.get("counts_plus", {})
    shots = int(dataset["shots"])
    x0 = density_to_cholesky_raw(linear_inversion_density(measured))

    def nll(raw: np.ndarray) -> float:
        pred = branch_expectations_from_density(cholesky_raw_to_density(raw))
        total = 0.0
        for axis in MEASUREMENT_AXES:
            count = counts_plus.get(axis)
            if count is None:
                continue
            prob = float(np.clip(0.5 * (1.0 + pred[axis]), EPS, 1.0 - EPS))
            total -= float(count) * math.log(prob)
            total -= float(shots - int(count)) * math.log(1.0 - prob)
        if total == 0.0:
            return float(sum((pred[axis] - measured[axis]) ** 2 for axis in MEASUREMENT_AXES))
        return total

    bounds = [(-8.0, 8.0), (-2.5, 2.5), (-2.5, 2.5), (-8.0, 8.0)]
    result = minimize(nll, x0, method="L-BFGS-B", bounds=bounds)
    if not result.success:
        result = minimize(nll, x0, method="Powell", options={"maxiter": 400})
    rho_opt = cholesky_raw_to_density(result.x)
    return {
        "rho": rho_opt,
        "objective": float(nll(result.x)),
        "expectations": branch_expectations_from_density(rho_opt),
        "optimizer_success": bool(result.success),
        "optimizer_message": str(result.message),
        "raw_params": np.asarray(result.x, dtype=float),
    }


def run_single_qubit_baseline(
    *,
    shot_grid: Sequence[int] = DEFAULT_SHOT_GRID,
    repeats: int = 40,
    seed: int = 1234,
    sample_mode: str = "binomial",
) -> dict[str, Any]:
    rho_true = single_qubit_target_state()
    rng = np.random.default_rng(int(seed))
    rows: list[dict[str, Any]] = []
    for shots in shot_grid:
        for trial in range(int(repeats)):
            dataset = single_qubit_dataset(
                rho_true,
                shots=int(shots),
                rng=rng,
                sample_mode=sample_mode,
            )
            fit_ls = fit_single_qubit_ls(dataset)
            fit_mle = fit_single_qubit_mle(dataset)
            rows.append(
                {
                    "shots": int(shots),
                    "trial": int(trial),
                    "sample_mode": str(sample_mode),
                    "fidelity_ls": density_matrix_fidelity(rho_true, fit_ls["rho"]),
                    "fidelity_mle": density_matrix_fidelity(rho_true, fit_mle["rho"]),
                    "trace_distance_ls": trace_distance(rho_true, fit_ls["rho"]),
                    "trace_distance_mle": trace_distance(rho_true, fit_mle["rho"]),
                }
            )
    summary_rows: list[dict[str, Any]] = []
    for shots in shot_grid:
        shot_rows = [row for row in rows if int(row["shots"]) == int(shots)]
        for method in ("ls", "mle"):
            values = np.asarray([float(row[f"fidelity_{method}"]) for row in shot_rows], dtype=float)
            summary_rows.append(
                {
                    "shots": int(shots),
                    "method": str(method),
                    "mean_fidelity": float(np.mean(values)),
                    "std_fidelity": float(np.std(values)),
                    "median_fidelity": float(np.median(values)),
                    "min_fidelity": float(np.min(values)),
                    "max_fidelity": float(np.max(values)),
                }
            )
    return {
        "true_state_bloch": dict(zip(("x", "y", "z"), bloch_from_density(rho_true), strict=True)),
        "rows": rows,
        "summary_rows": summary_rows,
    }


def measurement_settings_wait_only(wait_grid_s: Sequence[float] = DEFAULT_WAIT_GRID_S) -> list[MeasurementSetting]:
    return [
        MeasurementSetting(protocol="wait_only", alpha=0.0 + 0.0j, wait_s=float(wait_s), label=f"wait_{index:02d}")
        for index, wait_s in enumerate(wait_grid_s)
    ]


def measurement_settings_displacement_only(alpha_grid: Sequence[complex] = DEFAULT_ALPHA_GRID) -> list[MeasurementSetting]:
    return [
        MeasurementSetting(protocol="displacement_only", alpha=complex(alpha), wait_s=0.0, label=f"disp_{index:02d}")
        for index, alpha in enumerate(alpha_grid)
    ]


def measurement_settings_combined(
    alpha_grid: Sequence[complex] = DEFAULT_COMBINED_ALPHA_GRID,
    wait_grid_s: Sequence[float] = DEFAULT_COMBINED_WAIT_GRID_S,
) -> list[MeasurementSetting]:
    settings: list[MeasurementSetting] = []
    counter = 0
    for alpha in alpha_grid:
        for wait_s in wait_grid_s:
            settings.append(
                MeasurementSetting(
                    protocol="combined",
                    alpha=complex(alpha),
                    wait_s=float(wait_s),
                    label=f"comb_{counter:03d}",
                )
            )
            counter += 1
    return settings


def measurement_settings_by_protocol() -> dict[str, list[MeasurementSetting]]:
    return {
        "wait_only": measurement_settings_wait_only(),
        "displacement_only": measurement_settings_displacement_only(),
        "combined": measurement_settings_combined(),
    }


def dispersive_phase_array(wait_s: float, *, n_cav: int, chi: float = CHI, chi_prime: float = CHI_PRIME) -> np.ndarray:
    levels = np.arange(int(n_cav), dtype=float)
    return (float(chi) * levels + float(chi_prime) * levels * (levels - 1.0)) * float(wait_s)


def diagonal_kernel_vector(
    setting: MeasurementSetting,
    *,
    n_active: int = N_ACTIVE,
    kernel_dim: int = KERNEL_DIM,
    chi: float = CHI,
    chi_prime: float = CHI_PRIME,
) -> np.ndarray:
    if int(n_active) > int(kernel_dim):
        raise ValueError("kernel_dim must be at least n_active.")
    matrix = displacement_probability_matrix(int(kernel_dim), complex(setting.alpha))
    phases = np.exp(-1.0j * dispersive_phase_array(setting.wait_s, n_cav=kernel_dim, chi=chi, chi_prime=chi_prime))
    return np.asarray(
        [np.sum(matrix[:, level] * phases) for level in range(int(n_active))],
        dtype=np.complex128,
    )


def diagonal_predict_observables(
    settings: Sequence[MeasurementSetting],
    *,
    weighted_transverse: Sequence[complex],
    z_total: float,
    kernel_dim: int = KERNEL_DIM,
    chi: float = CHI,
    chi_prime: float = CHI_PRIME,
) -> dict[str, np.ndarray]:
    u = np.asarray(weighted_transverse, dtype=np.complex128)
    x_values: list[float] = []
    y_values: list[float] = []
    z_values: list[float] = []
    for setting in settings:
        kernel = diagonal_kernel_vector(
            setting,
            n_active=u.size,
            kernel_dim=kernel_dim,
            chi=chi,
            chi_prime=chi_prime,
        )
        transverse = np.sum(u * kernel)
        x_values.append(float(np.real(transverse)))
        y_values.append(float(np.imag(transverse)))
        z_values.append(float(z_total))
    return {
        "X": np.asarray(x_values, dtype=float),
        "Y": np.asarray(y_values, dtype=float),
        "Z": np.asarray(z_values, dtype=float),
    }


def diagonal_protocol_design_matrix(
    settings: Sequence[MeasurementSetting],
    *,
    n_active: int = N_ACTIVE,
    kernel_dim: int = KERNEL_DIM,
    chi: float = CHI,
    chi_prime: float = CHI_PRIME,
) -> np.ndarray:
    rows: list[np.ndarray] = []
    for setting in settings:
        kernel = diagonal_kernel_vector(
            setting,
            n_active=n_active,
            kernel_dim=kernel_dim,
            chi=chi,
            chi_prime=chi_prime,
        )
        row_x = np.zeros(2 * int(n_active) + 1, dtype=float)
        row_y = np.zeros(2 * int(n_active) + 1, dtype=float)
        row_z = np.zeros(2 * int(n_active) + 1, dtype=float)
        row_x[:n_active] = np.real(kernel)
        row_x[n_active : 2 * n_active] = -np.imag(kernel)
        row_y[:n_active] = np.imag(kernel)
        row_y[n_active : 2 * n_active] = np.real(kernel)
        row_z[-1] = 1.0
        rows.extend((row_x, row_y, row_z))
    return np.asarray(rows, dtype=float)


def protocol_identifiability_summary(
    settings: Sequence[MeasurementSetting],
    *,
    n_active: int = N_ACTIVE,
    kernel_dim: int = KERNEL_DIM,
) -> dict[str, Any]:
    matrix = diagonal_protocol_design_matrix(settings, n_active=n_active, kernel_dim=kernel_dim)
    transverse = matrix[:, : 2 * int(n_active)]
    svals = np.linalg.svd(matrix, compute_uv=False)
    svals_t = np.linalg.svd(transverse, compute_uv=False)
    cond_full = None if np.min(svals) <= EPS else float(np.max(svals) / np.min(svals))
    cond_transverse = None if np.min(svals_t) <= EPS else float(np.max(svals_t) / np.min(svals_t))
    return {
        "n_settings": int(len(settings)),
        "full_rank": int(np.linalg.matrix_rank(matrix)),
        "transverse_rank": int(np.linalg.matrix_rank(transverse)),
        "full_condition_number": cond_full,
        "transverse_condition_number": cond_transverse,
        "smallest_singular_value": float(np.min(svals)),
        "smallest_transverse_singular_value": float(np.min(svals_t)),
    }


def diagonal_dataset_from_weighted_truth(
    settings: Sequence[MeasurementSetting],
    *,
    weighted_transverse: Sequence[complex],
    z_total: float,
    noise: NoiseModel,
    rng: np.random.Generator,
    kernel_dim: int = KERNEL_DIM,
) -> dict[str, Any]:
    alpha_scale = 1.0 + rng.normal(0.0, float(noise.displacement_rel_sigma))
    alpha_phase = rng.normal(0.0, float(noise.displacement_phase_sigma_rad))
    chi_actual = float(CHI) * (1.0 + rng.normal(0.0, float(noise.chi_rel_sigma)))
    chi_prime_actual = float(CHI_PRIME) * (1.0 + rng.normal(0.0, float(noise.chi_prime_rel_sigma)))
    angle_errors = {
        axis: float(rng.normal(0.0, float(noise.rotation_sigma_rad)))
        for axis in MEASUREMENT_AXES
    }

    rows: list[dict[str, Any]] = []
    for setting in settings:
        alpha_actual = complex(setting.alpha) * complex(alpha_scale * np.exp(1.0j * alpha_phase))
        setting_actual = MeasurementSetting(
            protocol=setting.protocol,
            alpha=alpha_actual,
            wait_s=float(setting.wait_s),
            label=setting.label,
        )
        exact = diagonal_predict_observables(
            [setting_actual],
            weighted_transverse=weighted_transverse,
            z_total=float(z_total),
            kernel_dim=kernel_dim,
            chi=chi_actual,
            chi_prime=chi_prime_actual,
        )
        bloch = {axis: float(exact[axis][0]) for axis in MEASUREMENT_AXES}
        if noise.t2_s is not None and float(noise.t2_s) > 0.0:
            decay = math.exp(-float(setting.wait_s) / float(noise.t2_s))
            bloch["X"] *= decay
            bloch["Y"] *= decay
        if noise.t1_s is not None and float(noise.t1_s) > 0.0:
            bloch["Z"] = 1.0 - (1.0 - bloch["Z"]) * math.exp(-float(setting.wait_s) / float(noise.t1_s))
        rho_meas = linear_inversion_density(bloch)
        branch_exact = branch_expectations_from_density(rho_meas, angle_errors=angle_errors)
        for axis in MEASUREMENT_AXES:
            measured, counts_plus = sample_branch_measurement(
                branch_exact[axis],
                shots=int(noise.shots),
                rng=rng,
                mode=noise.sample_mode,
            )
            rows.append(
                {
                    **setting.as_dict(),
                    "axis": str(axis),
                    "exact_expectation": float(branch_exact[axis]),
                    "measured_expectation": float(measured),
                    "counts_plus": counts_plus,
                    "shots": int(noise.shots),
                }
            )
    return {
        "rows": rows,
        "systematics": {
            "alpha_scale": float(alpha_scale),
            "alpha_phase_rad": float(alpha_phase),
            "chi_actual": float(chi_actual),
            "chi_prime_actual": float(chi_prime_actual),
            "rotation_angle_errors_rad": dict(angle_errors),
        },
        "noise_model": noise.as_dict(),
    }


def _ordered_row_vector(rows: Sequence[dict[str, Any]], settings: Sequence[MeasurementSetting], field: str) -> np.ndarray:
    mapping = {(str(row["label"]), str(row["axis"]).upper()): row for row in rows}
    out: list[float] = []
    for setting in settings:
        for axis in MEASUREMENT_AXES:
            out.append(float(mapping[(setting.label, axis)][field]))
    return np.asarray(out, dtype=float)


def _ordered_counts_vector(rows: Sequence[dict[str, Any]], settings: Sequence[MeasurementSetting]) -> tuple[np.ndarray, np.ndarray]:
    mapping = {(str(row["label"]), str(row["axis"]).upper()): row for row in rows}
    counts: list[float] = []
    shots: list[float] = []
    for setting in settings:
        for axis in MEASUREMENT_AXES:
            row = mapping[(setting.label, axis)]
            count = row.get("counts_plus")
            if count is None:
                meas = float(row["measured_expectation"])
                shot = int(row["shots"])
                count = int(round(0.5 * (1.0 + np.clip(meas, -1.0, 1.0)) * shot))
            counts.append(float(count))
            shots.append(float(row["shots"]))
    return np.asarray(counts, dtype=float), np.asarray(shots, dtype=float)


def unpack_weighted_parameter_vector(theta: Sequence[float], *, n_active: int) -> tuple[np.ndarray, float]:
    arr = np.asarray(theta, dtype=float).reshape(2 * int(n_active) + 1)
    real = arr[:n_active]
    imag = arr[n_active : 2 * n_active]
    return real + 1.0j * imag, float(arr[-1])


def weighted_parameter_vector(weighted_transverse: Sequence[complex], z_total: float) -> np.ndarray:
    u = np.asarray(weighted_transverse, dtype=np.complex128)
    return np.concatenate([np.real(u), np.imag(u), np.asarray([float(z_total)], dtype=float)])


def infer_weighted_ls(
    rows: Sequence[dict[str, Any]],
    settings: Sequence[MeasurementSetting],
    *,
    n_active: int = N_ACTIVE,
    kernel_dim: int = KERNEL_DIM,
) -> dict[str, Any]:
    matrix = diagonal_protocol_design_matrix(settings, n_active=n_active, kernel_dim=kernel_dim)
    y = _ordered_row_vector(rows, settings, "measured_expectation")
    theta, *_ = np.linalg.lstsq(matrix, y, rcond=None)
    weighted_transverse, z_total = unpack_weighted_parameter_vector(theta, n_active=n_active)
    pred = matrix @ theta
    residual = y - pred
    return {
        "weighted_transverse": weighted_transverse,
        "z_total": float(z_total),
        "theta": np.asarray(theta, dtype=float),
        "predicted_vector": pred,
        "residual_rms": float(np.sqrt(np.mean(residual ** 2))),
        "residual_max_abs": float(np.max(np.abs(residual))),
    }


def infer_weighted_mle(
    rows: Sequence[dict[str, Any]],
    settings: Sequence[MeasurementSetting],
    *,
    n_active: int = N_ACTIVE,
    kernel_dim: int = KERNEL_DIM,
) -> dict[str, Any]:
    matrix = diagonal_protocol_design_matrix(settings, n_active=n_active, kernel_dim=kernel_dim)
    counts_plus, shots = _ordered_counts_vector(rows, settings)
    ls = infer_weighted_ls(rows, settings, n_active=n_active, kernel_dim=kernel_dim)
    x0 = np.asarray(ls["theta"], dtype=float)
    bounds = [(-1.0, 1.0)] * (2 * int(n_active) + 1)

    def objective(theta: np.ndarray) -> float:
        pred = np.clip(matrix @ np.asarray(theta, dtype=float), -0.999999, 0.999999)
        prob = np.clip(0.5 * (1.0 + pred), EPS, 1.0 - EPS)
        return float(-np.sum(counts_plus * np.log(prob) + (shots - counts_plus) * np.log(1.0 - prob)))

    result = minimize(objective, x0, method="L-BFGS-B", bounds=bounds)
    if not result.success:
        result = minimize(objective, x0, method="Powell", options={"maxiter": 800})
    weighted_transverse, z_total = unpack_weighted_parameter_vector(result.x, n_active=n_active)
    pred = matrix @ np.asarray(result.x, dtype=float)
    residual = _ordered_row_vector(rows, settings, "measured_expectation") - pred
    return {
        "weighted_transverse": weighted_transverse,
        "z_total": float(z_total),
        "theta": np.asarray(result.x, dtype=float),
        "predicted_vector": pred,
        "objective": float(objective(result.x)),
        "residual_rms": float(np.sqrt(np.mean(residual ** 2))),
        "residual_max_abs": float(np.max(np.abs(residual))),
        "optimizer_success": bool(result.success),
        "optimizer_message": str(result.message),
    }


def recoverable_error_summary(
    inferred_weighted: Sequence[complex],
    true_weighted: Sequence[complex],
    *,
    inferred_z_total: float,
    true_z_total: float,
    true_populations: Sequence[float] | None = None,
) -> dict[str, Any]:
    inferred = np.asarray(inferred_weighted, dtype=np.complex128)
    truth = np.asarray(true_weighted, dtype=np.complex128)
    errors = inferred - truth
    payload: dict[str, Any] = {
        "per_sector_abs_error": np.abs(errors).tolist(),
        "per_sector_weighted_x_error": np.real(errors).tolist(),
        "per_sector_weighted_y_error": np.imag(errors).tolist(),
        "mean_abs_error": float(np.mean(np.abs(errors))),
        "max_abs_error": float(np.max(np.abs(errors))),
        "z_total_error": float(float(inferred_z_total) - float(true_z_total)),
        "weighted_rmse": float(np.sqrt(np.mean(np.real(errors) ** 2 + np.imag(errors) ** 2))),
    }
    if true_populations is not None:
        probs = np.asarray(true_populations, dtype=float)
        assist = np.zeros_like(inferred)
        truth_norm = np.zeros_like(truth)
        mask = probs > 1.0e-6
        assist[mask] = inferred[mask] / probs[mask]
        truth_norm[mask] = truth[mask] / probs[mask]
        payload["oracle_transverse_component_error"] = np.abs(assist - truth_norm).tolist()
        payload["oracle_transverse_component_rmse"] = float(np.sqrt(np.mean(np.abs(assist[mask] - truth_norm[mask]) ** 2))) if np.any(mask) else 0.0
    return payload


def construct_exact_gauge_family(
    weighted_transverse: Sequence[complex],
    z_total: float,
    *,
    patterns: Sequence[Sequence[float]] | None = None,
) -> list[dict[str, Any]]:
    u = np.asarray(weighted_transverse, dtype=np.complex128)
    min_prob = np.abs(u)
    slack_total = float(1.0 - np.sum(min_prob))
    if slack_total <= 1.0e-8:
        return []
    if patterns is None:
        patterns = (
            (4.0, 1.0, 3.0, 2.0),
            (1.0, 3.0, 2.0, 4.0),
        )
    solutions: list[dict[str, Any]] = []
    for pattern in patterns:
        weights = normalize_probabilities(pattern[: u.size])
        probs = min_prob + slack_total * weights
        caps = np.sqrt(np.maximum(probs ** 2 - np.abs(u) ** 2, 0.0))
        cap_sum = float(np.sum(caps))
        if cap_sum <= EPS or abs(float(z_total)) > cap_sum + 1.0e-9:
            continue
        weighted_z = float(z_total) * caps / cap_sum
        z_values = np.divide(weighted_z, probs, out=np.zeros_like(weighted_z), where=probs > EPS)
        bloch_rows: list[dict[str, float]] = []
        for index in range(u.size):
            xy = u[index] / probs[index]
            bloch_rows.append(
                {
                    "x": float(np.real(xy)),
                    "y": float(np.imag(xy)),
                    "z": float(z_values[index]),
                }
            )
        solutions.append(
            {
                "pattern": [float(value) for value in pattern[: u.size]],
                "probabilities": probs.tolist(),
                "bloch_rows": bloch_rows,
                "reconstructed_z_total": float(np.sum(probs * z_values)),
                "reconstructed_weighted_transverse": [complex_to_dict(value) for value in u],
            }
        )
    return solutions


def full_state_random_initial(n_active: int, rng: np.random.Generator) -> np.ndarray:
    logits = rng.normal(0.0, 0.8, size=int(n_active))
    raw_states = rng.normal(0.0, 0.5, size=4 * int(n_active))
    return np.concatenate([logits, raw_states])


def decode_full_state_params(raw: Sequence[float], *, n_active: int) -> tuple[np.ndarray, list[qt.Qobj]]:
    arr = np.asarray(raw, dtype=float).reshape(int(n_active) + 4 * int(n_active))
    logits = arr[:n_active]
    probs = np.exp(logits - logsumexp(logits))
    states: list[qt.Qobj] = []
    offset = n_active
    for _ in range(int(n_active)):
        states.append(cholesky_raw_to_density(arr[offset : offset + 4]))
        offset += 4
    return probs, states


def full_state_prediction_vector(
    theta: Sequence[float],
    settings: Sequence[MeasurementSetting],
    *,
    n_active: int = N_ACTIVE,
    kernel_dim: int = KERNEL_DIM,
) -> np.ndarray:
    probs, states = decode_full_state_params(theta, n_active=n_active)
    weighted_transverse = np.asarray(
        [
            complex(float(probs[index]) * float(bloch_from_density(state)[0]), float(probs[index]) * float(bloch_from_density(state)[1]))
            for index, state in enumerate(states)
        ],
        dtype=np.complex128,
    )
    z_total = float(sum(float(probs[index]) * float(bloch_from_density(state)[2]) for index, state in enumerate(states)))
    matrix = diagonal_protocol_design_matrix(settings, n_active=n_active, kernel_dim=kernel_dim)
    return matrix @ weighted_parameter_vector(weighted_transverse, z_total)


def run_full_state_restart_diagnostic(
    rows: Sequence[dict[str, Any]],
    settings: Sequence[MeasurementSetting],
    *,
    n_active: int = N_ACTIVE,
    kernel_dim: int = KERNEL_DIM,
    n_restarts: int = 8,
    seed: int = 4321,
) -> dict[str, Any]:
    rng = np.random.default_rng(int(seed))
    counts_plus, shots = _ordered_counts_vector(rows, settings)
    measured = _ordered_row_vector(rows, settings, "measured_expectation")
    results: list[dict[str, Any]] = []
    bounds = [(-5.0, 5.0)] * int(n_active) + [(-8.0, 8.0), (-2.5, 2.5), (-2.5, 2.5), (-8.0, 8.0)] * int(n_active)

    def objective(theta: np.ndarray) -> float:
        pred = np.clip(
            full_state_prediction_vector(theta, settings, n_active=n_active, kernel_dim=kernel_dim),
            -0.999999,
            0.999999,
        )
        prob = np.clip(0.5 * (1.0 + pred), EPS, 1.0 - EPS)
        return float(-np.sum(counts_plus * np.log(prob) + (shots - counts_plus) * np.log(1.0 - prob)))

    for restart in range(int(n_restarts)):
        x0 = full_state_random_initial(n_active, rng)
        result = minimize(objective, x0, method="L-BFGS-B", bounds=bounds)
        if not result.success:
            result = minimize(objective, x0, method="Powell", options={"maxiter": 800})
        probs, states = decode_full_state_params(result.x, n_active=n_active)
        bloch_rows = []
        weighted = []
        for prob, state in zip(probs, states, strict=True):
            x, y, z = bloch_from_density(state)
            bloch_rows.append({"x": float(x), "y": float(y), "z": float(z)})
            weighted.append(complex(float(prob) * float(x), float(prob) * float(y)))
        pred = full_state_prediction_vector(result.x, settings, n_active=n_active, kernel_dim=kernel_dim)
        residual = measured - pred
        results.append(
            {
                "restart": int(restart),
                "objective": float(objective(result.x)),
                "residual_rms": float(np.sqrt(np.mean(residual ** 2))),
                "optimizer_success": bool(result.success),
                "optimizer_message": str(result.message),
                "probabilities": probs.tolist(),
                "bloch_rows": bloch_rows,
                "weighted_transverse": [complex_to_dict(value) for value in weighted],
                "z_total": float(sum(float(probs[index]) * float(bloch_rows[index]["z"]) for index in range(int(n_active)))),
            }
        )
    results.sort(key=lambda item: float(item["objective"]))
    probability_matrix = np.asarray([item["probabilities"] for item in results], dtype=float)
    z_matrix = np.asarray([[row["z"] for row in item["bloch_rows"]] for item in results], dtype=float)
    return {
        "rows": results,
        "objective_span": float(max(item["objective"] for item in results) - min(item["objective"] for item in results)),
        "probability_std_by_sector": np.std(probability_matrix, axis=0).tolist(),
        "z_std_by_sector": np.std(z_matrix, axis=0).tolist(),
    }


def make_wait_unitary(model: Any, frame: Any, wait_s: float) -> qt.Qobj:
    if float(wait_s) <= 0.0:
        return qt.tensor(qt.qeye(int(model.n_tr)), qt.qeye(int(model.n_cav)))
    h0 = model.static_hamiltonian(frame=frame)
    return (-1.0j * float(wait_s) * h0).expm()


def analysis_evolved_state(
    state: qt.Qobj,
    *,
    model: Any,
    frame: Any,
    alpha: complex,
    wait_s: float,
) -> qt.Qobj:
    disp = qt.tensor(qt.qeye(int(model.n_tr)), displacement_op(int(model.n_cav), complex(alpha)))
    wait = make_wait_unitary(model, frame, float(wait_s))
    return wait * (disp * state * disp.dag()) * wait.dag()


def exact_joint_dataset(
    state: qt.Qobj,
    settings: Sequence[MeasurementSetting],
    *,
    model: Any,
    frame: Any,
    noise: NoiseModel,
    rng: np.random.Generator,
) -> dict[str, Any]:
    alpha_scale = 1.0 + rng.normal(0.0, float(noise.displacement_rel_sigma))
    alpha_phase = rng.normal(0.0, float(noise.displacement_phase_sigma_rad))
    angle_errors = {axis: float(rng.normal(0.0, float(noise.rotation_sigma_rad))) for axis in MEASUREMENT_AXES}
    rows: list[dict[str, Any]] = []
    for setting in settings:
        alpha_actual = complex(setting.alpha) * complex(alpha_scale * np.exp(1.0j * alpha_phase))
        evolved = analysis_evolved_state(
            state,
            model=model,
            frame=frame,
            alpha=alpha_actual,
            wait_s=float(setting.wait_s),
        )
        exact = branch_expectations_from_joint_state(
            evolved,
            n_tr=int(model.n_tr),
            angle_errors=angle_errors,
        )
        for axis in MEASUREMENT_AXES:
            mode = "gaussian" if int(model.n_tr) > 2 or str(noise.sample_mode) == "gaussian" else "binomial"
            measured, counts_plus = sample_branch_measurement(
                exact[axis],
                shots=int(noise.shots),
                rng=rng,
                mode=mode,
            )
            rows.append(
                {
                    **setting.as_dict(),
                    "axis": str(axis),
                    "exact_expectation": float(exact[axis]),
                    "measured_expectation": float(measured),
                    "counts_plus": counts_plus,
                    "shots": int(noise.shots),
                }
            )
    return {
        "rows": rows,
        "noise_model": noise.as_dict(),
        "systematics": {
            "alpha_scale": float(alpha_scale),
            "alpha_phase_rad": float(alpha_phase),
            "rotation_angle_errors_rad": dict(angle_errors),
        },
    }


def padded_rotation_arrays(*, n_cav: int) -> tuple[np.ndarray, np.ndarray]:
    theta = np.zeros(int(n_cav), dtype=float)
    phi = np.zeros(int(n_cav), dtype=float)
    theta[: N_ACTIVE] = np.asarray(DEFAULT_THETA, dtype=float)
    phi[: N_ACTIVE] = np.asarray(DEFAULT_PHI, dtype=float)
    return theta, phi


def padded_phase_array(*, n_cav: int) -> np.ndarray:
    phases = np.zeros(int(n_cav), dtype=float)
    phases[: N_ACTIVE] = np.asarray(DEFAULT_CPSQR_PHASES, dtype=float)
    return phases


def cavity_population_vector(n_cav: int, probabilities: Sequence[float] = DEFAULT_POPULATIONS) -> np.ndarray:
    probs = normalize_probabilities(probabilities)
    padded = np.zeros(int(n_cav), dtype=float)
    padded[: probs.size] = probs
    return padded


def build_ideal_reference_cases(*, n_cav: int = 6) -> list[StudyCase]:
    model = build_model(n_cav=n_cav, n_tr=2)
    frame = build_frame(model)
    theta, phi = padded_rotation_arrays(n_cav=n_cav)
    sqr = ideal_sqr_operator(theta, phi)
    phases = padded_phase_array(n_cav=n_cav)
    cpsqr = cpsqr_like_operator(phases)
    populations = cavity_population_vector(n_cav)
    cavity_diag = cavity_mixture_state(populations)
    cavity_coherent = cavity_superposition_state(np.sqrt(populations.astype(np.complex128)))
    cases: list[StudyCase] = []
    for qubit_label in ("g", "e", "+", "+y"):
        rho_in = qt.tensor(qubit_input_state(qubit_label), cavity_diag)
        rho_out = sqr * rho_in * sqr.dag()
        cases.append(
            StudyCase(
                case_id=f"ideal_diag_{qubit_label}",
                family="analytic",
                description=f"Ideal block-diagonal SQR on a cavity-diagonal input with qubit {qubit_label}.",
                model_class="fock_diagonal",
                state=rho_out,
                model=model,
                frame=frame,
                truth_sectors=tuple(sector_summaries_from_state(rho_out, range(N_ACTIVE))),
                metadata={"gate": "ideal_sqr", "qubit_input": qubit_label, "cavity_input": "mixture"},
            )
        )
    rho_coherent = qt.tensor(qubit_ground_dm(), as_dm(cavity_coherent))
    rho_coherent_out = sqr * rho_coherent * sqr.dag()
    cases.append(
        StudyCase(
            case_id="ideal_coherent_ground",
            family="analytic",
            description="Ideal block-diagonal SQR on a cavity superposition; same sector marginals as the diagonal benchmark but with cavity coherences.",
            model_class="coherent_block",
            state=rho_coherent_out,
            model=model,
            frame=frame,
            truth_sectors=tuple(sector_summaries_from_state(rho_coherent_out, range(N_ACTIVE))),
            metadata={"gate": "ideal_sqr", "qubit_input": "g", "cavity_input": "superposition"},
        )
    )
    rho_cpsqr = qt.tensor(qubit_plus_dm(), cavity_diag)
    rho_cpsqr_out = cpsqr * rho_cpsqr * cpsqr.dag()
    cases.append(
        StudyCase(
            case_id="cpsqr_like_plus",
            family="analytic",
            description="CPSQR-like conditional phase gate on a cavity-diagonal input with qubit |+>.",
            model_class="fock_diagonal",
            state=rho_cpsqr_out,
            model=model,
            frame=frame,
            truth_sectors=tuple(sector_summaries_from_state(rho_cpsqr_out, range(N_ACTIVE))),
            metadata={"gate": "cpsqr_like", "qubit_input": "+", "cavity_input": "mixture"},
        )
    )
    return cases


def qubit_density_for_model_input(model: Any, qubit_label: str) -> qt.Qobj:
    return embed_qubit_density(qubit_input_state(qubit_label), int(model.n_tr))


def build_multitone_gate_payload(
    *,
    n_cav: int = 6,
    duration_s: float = NEAR_IDEAL_DURATION_S,
    optimize: bool = True,
    simulation_mode: str = "reduced",
    label: str = "study_multitone",
    maxiter_stage1: int = 8,
    maxiter_stage2: int = 12,
) -> dict[str, Any]:
    model = build_model(n_cav=n_cav, n_tr=2)
    targets = conditioned_multitone_targets(DEFAULT_THETA, DEFAULT_PHI)
    run_config = conditioned_run_config(model, duration_s=duration_s)
    waveform, _, _ = build_multitone_waveform(model, duration_s=duration_s)
    payload: dict[str, Any] = {
        "model": model,
        "frame": run_config.frame,
        "targets": targets,
        "run_config": run_config,
        "compiled": compile_waveform(waveform, run_config),
        "waveform": waveform,
        "drive_ops": {"qubit": "qubit"},
        "optimization": None,
        "construction": "seed",
    }
    if not optimize:
        return payload
    opt_cfg = ConditionedOptimizationConfig(
        active_levels=tuple(range(N_ACTIVE)),
        parameters=("d_lambda", "d_alpha", "d_omega"),
        maxiter_stage1=int(maxiter_stage1),
        maxiter_stage2=int(maxiter_stage2),
        d_lambda_bounds=(-0.75, 0.75),
        d_alpha_bounds=(-np.pi, np.pi),
        d_omega_hz_bounds=(-3.0e6, 3.0e6),
        regularization_lambda=5.0e-4,
        regularization_alpha=5.0e-4,
        regularization_omega=5.0e-4,
    )
    start = time.perf_counter()
    optimization = optimize_conditioned_multitone(
        model,
        targets,
        run_config,
        optimization_config=opt_cfg,
        simulation_mode=simulation_mode,
        label=label,
    )
    runtime_s = time.perf_counter() - start
    optimized = optimization.optimized_result
    payload.update(
        {
            "compiled": optimized.compiled,
            "waveform": optimized.waveform,
            "optimization": optimization,
            "optimization_runtime_s": float(runtime_s),
            "drive_ops": {str(optimized.waveform.drive_channel): "qubit"},
            "construction": "optimized",
        }
    )
    return payload


def simulate_pulse_case_state(
    *,
    compiled: Any,
    model: Any,
    frame: Any,
    drive_ops: dict[str, str],
    qubit_label: str,
    cavity_label: str,
) -> qt.Qobj:
    if cavity_label == "mixture":
        cavity_state = cavity_mixture_state(cavity_population_vector(int(model.n_cav)))
    elif cavity_label == "superposition":
        cavity_state = as_dm(cavity_superposition_state(np.sqrt(cavity_population_vector(int(model.n_cav)).astype(np.complex128))))
    elif cavity_label == "coherent":
        cavity_state = as_dm(coherent_state(0.55 + 0.15j, n_cav=int(model.n_cav)))
    else:
        raise ValueError(f"Unsupported cavity label: {cavity_label}")
    rho0 = qt.tensor(qubit_density_for_model_input(model, qubit_label), cavity_state)
    return simulate_compiled_on_states(
        model,
        compiled,
        frame=frame,
        drive_ops=drive_ops,
        initial_states=[rho0],
    )[0]


def build_pulse_level_cases() -> tuple[list[StudyCase], dict[str, Any]]:
    optimized_payload = build_multitone_gate_payload(
        n_cav=6,
        duration_s=NEAR_IDEAL_DURATION_S,
        optimize=True,
        simulation_mode="reduced",
        label="near_ideal_black_box",
        maxiter_stage1=6,
        maxiter_stage2=8,
    )
    seed_payload = build_multitone_gate_payload(
        n_cav=6,
        duration_s=IMPERFECT_DURATION_S,
        optimize=False,
        label="seed_imperfect_black_box",
    )
    cases: list[StudyCase] = []
    for qubit_label in ("g", "+"):
        state = simulate_pulse_case_state(
            compiled=optimized_payload["compiled"],
            model=optimized_payload["model"],
            frame=optimized_payload["frame"],
            drive_ops=optimized_payload["drive_ops"],
            qubit_label=qubit_label,
            cavity_label="mixture",
        )
        cases.append(
            StudyCase(
                case_id=f"pulse_optimized_mix_{qubit_label}",
                family="pulse_level",
                description=f"Optimized multitone black-box SQR on a cavity mixture with qubit {qubit_label}.",
                model_class="pulse_level",
                state=state,
                model=optimized_payload["model"],
                frame=optimized_payload["frame"],
                truth_sectors=tuple(sector_summaries_from_state(state, range(N_ACTIVE))),
                metadata={"construction": "optimized_multitone", "qubit_input": qubit_label, "cavity_input": "mixture"},
            )
        )
    coherent_state_out = simulate_pulse_case_state(
        compiled=optimized_payload["compiled"],
        model=optimized_payload["model"],
        frame=optimized_payload["frame"],
        drive_ops=optimized_payload["drive_ops"],
        qubit_label="g",
        cavity_label="superposition",
    )
    cases.append(
        StudyCase(
            case_id="pulse_optimized_superposition_g",
            family="pulse_level",
            description="Optimized multitone black-box SQR on a cavity superposition with qubit |g>.",
            model_class="pulse_level",
            state=coherent_state_out,
            model=optimized_payload["model"],
            frame=optimized_payload["frame"],
            truth_sectors=tuple(sector_summaries_from_state(coherent_state_out, range(N_ACTIVE))),
            metadata={"construction": "optimized_multitone", "qubit_input": "g", "cavity_input": "superposition"},
        )
    )
    imperfect_state = simulate_pulse_case_state(
        compiled=seed_payload["compiled"],
        model=seed_payload["model"],
        frame=seed_payload["frame"],
        drive_ops=seed_payload["drive_ops"],
        qubit_label="g",
        cavity_label="superposition",
    )
    cases.append(
        StudyCase(
            case_id="pulse_seed_superposition_g",
            family="pulse_level",
            description="Short, unoptimized multitone seed waveform on a cavity superposition with qubit |g>.",
            model_class="pulse_level",
            state=imperfect_state,
            model=seed_payload["model"],
            frame=seed_payload["frame"],
            truth_sectors=tuple(sector_summaries_from_state(imperfect_state, range(N_ACTIVE))),
            metadata={"construction": "seed_multitone", "qubit_input": "g", "cavity_input": "superposition"},
        )
    )
    leakage_model = build_model(n_cav=6, n_tr=3)
    leakage_frame = build_frame(leakage_model)
    leakage_init = qt.tensor(
        embed_qubit_density(qubit_ground_dm(), int(leakage_model.n_tr)),
        as_dm(cavity_superposition_state(np.sqrt(cavity_population_vector(int(leakage_model.n_cav)).astype(np.complex128)))),
    )
    leakage_state = simulate_compiled_on_states(
        leakage_model,
        optimized_payload["compiled"],
        frame=leakage_frame,
        drive_ops=optimized_payload["drive_ops"],
        initial_states=[leakage_init],
    )[0]
    cases.append(
        StudyCase(
            case_id="pulse_optimized_leakage_replay",
            family="pulse_level",
            description="Optimized multitone waveform replayed on a three-level transmon model to expose leakage sensitivity.",
            model_class="pulse_level_leakage",
            state=leakage_state,
            model=leakage_model,
            frame=leakage_frame,
            truth_sectors=tuple(sector_summaries_from_state(leakage_state, range(N_ACTIVE))),
            metadata={"construction": "optimized_multitone_replay_ntr3", "qubit_input": "g", "cavity_input": "superposition"},
        )
    )
    return cases, {
        "optimized_payload": optimized_payload,
        "seed_payload": seed_payload,
    }


def serialize_case(case: StudyCase) -> dict[str, Any]:
    return {
        "case_id": str(case.case_id),
        "family": str(case.family),
        "description": str(case.description),
        "model_class": str(case.model_class),
        "truth_sectors": [sector.as_dict() for sector in case.truth_sectors],
        "metadata": json_ready(case.metadata),
    }


def protocol_fit_for_case(
    case: StudyCase,
    *,
    settings: Sequence[MeasurementSetting],
    noise: NoiseModel,
    rng: np.random.Generator,
    fit_method: str = "ls",
    kernel_dim: int | None = None,
) -> dict[str, Any]:
    use_kernel_dim = int(case.model.n_cav) if kernel_dim is None else int(kernel_dim)
    dataset = exact_joint_dataset(
        case.state,
        settings,
        model=case.model,
        frame=case.frame,
        noise=noise,
        rng=rng,
    )
    if str(fit_method) == "mle":
        fit = infer_weighted_mle(dataset["rows"], settings, n_active=N_ACTIVE, kernel_dim=use_kernel_dim)
    else:
        fit = infer_weighted_ls(dataset["rows"], settings, n_active=N_ACTIVE, kernel_dim=use_kernel_dim)
    return {
        "case": serialize_case(case),
        "protocol": str(settings[0].protocol if settings else "unknown"),
        "fit_method": str(fit_method),
        "noise": noise.as_dict(),
        "dataset": dataset,
        "fit": {
            "weighted_transverse": [complex_to_dict(value) for value in fit["weighted_transverse"]],
            "z_total": float(fit["z_total"]),
            "residual_rms": float(fit["residual_rms"]),
            "residual_max_abs": float(fit["residual_max_abs"]),
        },
        "fit_errors": recoverable_error_summary(
            fit["weighted_transverse"],
            case.weighted_transverse_truth,
            inferred_z_total=float(fit["z_total"]),
            true_z_total=float(case.z_total_truth),
            true_populations=case.populations_truth,
        ),
    }
