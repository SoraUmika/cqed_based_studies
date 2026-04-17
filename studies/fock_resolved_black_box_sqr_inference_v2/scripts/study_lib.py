"""Inference and case-generation library for the Fock-resolved SQR study (v2).

Structure mirrors the open study plan exactly:
  Part 1  — single-qubit MLE building block (forward sim + LS + MLE + F-vs-N)
  Part 2  — full joint-system study
    2A  forward model and identifiability
    2B  full {p_n, rho_q^(n)} MLE — attempted first, non-uniqueness characterised
    2C  recoverable-subspace inference (weighted-transverse LS and MLE)
    2D  per-sector metrics (fidelity, trace distance, Bloch-vector error)
    2E  black-box case library (ideal analytic + pulse-level)
    2F  Model-B coherence sweep
    2G  robustness / noise sweep
  Part 3  — comparison questions helpers
"""

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
from cqed_sim.core import FrameSpec, displacement_op
from cqed_sim.core.ideal_gates import qubit_rotation_axis, qubit_rotation_xy

from common import (
    CHI,
    CHI_PRIME,
    DEFAULT_CPSQR_PHASES,
    DEFAULT_PHI,
    DEFAULT_POPULATIONS,
    DEFAULT_THETA,
    LONG_DURATION_S,
    SHORT_DURATION_S,
    N_ACTIVE,
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

# Coherence-sweep mixing fractions: 0 = purely diagonal, 1 = fully coherent.
DEFAULT_COHERENCE_FRACTIONS = (0.0, 0.1, 0.25, 0.5, 0.75, 1.0)


# ---------------------------------------------------------------------------
# Shared data structures
# ---------------------------------------------------------------------------

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
            [s.population * (s.x + 1.0j * s.y) for s in self.truth_sectors],
            dtype=np.complex128,
        )

    @property
    def z_total_truth(self) -> float:
        return float(sum(s.population * s.z for s in self.truth_sectors))

    @property
    def populations_truth(self) -> np.ndarray:
        return np.asarray([s.population for s in self.truth_sectors], dtype=float)

    @property
    def normalized_qubit_states_truth(self) -> list[qt.Qobj]:
        """Normalized rho_q^(n) from ground truth, for per-sector fidelity."""
        states = []
        for s in self.truth_sectors:
            states.append(qubit_density_from_bloch(s.x, s.y, s.z))
        return states


def complex_to_dict(value: complex) -> dict[str, float]:
    return {"real": float(np.real(value)), "imag": float(np.imag(value))}


# ---------------------------------------------------------------------------
# Part 1 — Single-qubit tomography building block
# ---------------------------------------------------------------------------

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


def embedded_sigma_z(n_tr: int) -> qt.Qobj:
    diag = np.zeros(int(n_tr), dtype=np.complex128)
    diag[0] = 1.0
    if int(n_tr) > 1:
        diag[1] = -1.0
    return qt.Qobj(np.diag(diag), dims=[[int(n_tr)], [int(n_tr)]])


def tomography_prerotation(axis: str, angle_error: float = 0.0) -> qt.Qobj:
    """Pre-measurement unitary for tomography axis.

    Branch   Pre-rotation
    Z        I
    X        Ry(-pi/2)   rotates X-eigenstate into Z-eigenstate
    Y        Rx(+pi/2)   rotates Y-eigenstate into Z-eigenstate
    """
    key = str(axis).upper()
    delta = float(angle_error)
    if key == "Z":
        return qt.qeye(2)
    if key == "X":
        return qubit_rotation_xy(-0.5 * np.pi + delta, 0.5 * np.pi)
    if key == "Y":
        return qubit_rotation_xy(0.5 * np.pi + delta, 0.0)
    raise ValueError(f"Unsupported tomography axis: {axis}")


def branch_expectations_from_density(
    rho_q: qt.Qobj,
    *,
    angle_errors: dict[str, float] | None = None,
) -> dict[str, float]:
    """E_k = Tr(U_k rho_q U_k† sigma_z) for k in {X, Y, Z}."""
    errors = {} if angle_errors is None else {str(k).upper(): float(v) for k, v in angle_errors.items()}
    return {
        axis: float(np.real(qt.expect(PAULI_Z, tomography_prerotation(axis, errors.get(axis, 0.0)) * rho_q * tomography_prerotation(axis, errors.get(axis, 0.0)).dag())))
        for axis in MEASUREMENT_AXES
    }


def branch_expectations_from_joint_state(
    state: qt.Qobj,
    *,
    n_tr: int,
    angle_errors: dict[str, float] | None = None,
) -> dict[str, float]:
    errors = {} if angle_errors is None else {str(k).upper(): float(v) for k, v in angle_errors.items()}
    reduced = state.ptrace(0)
    sigma_z = embedded_sigma_z(n_tr)
    return {
        axis: float(np.real(qt.expect(sigma_z, embed_qubit_operator(tomography_prerotation(axis, errors.get(axis, 0.0)), n_tr=n_tr) * reduced * embed_qubit_operator(tomography_prerotation(axis, errors.get(axis, 0.0)), n_tr=n_tr).dag())))
        for axis in MEASUREMENT_AXES
    }


def project_density_to_physical(matrix: np.ndarray) -> np.ndarray:
    herm = 0.5 * (matrix + matrix.conj().T)
    evals, evecs = np.linalg.eigh(herm)
    evals = np.clip(np.real(evals), 0.0, None)
    if float(np.sum(evals)) <= EPS:
        return np.asarray(qubit_ground_dm().full(), dtype=np.complex128)
    clipped = (evecs * evals) @ evecs.conj().T
    return clipped / float(np.real(np.trace(clipped)))


def cholesky_raw_to_density(raw: Sequence[float]) -> qt.Qobj:
    """T†T / Tr(T†T) where T is lower-triangular with positive diagonal."""
    values = np.asarray(raw, dtype=float).reshape(4)
    t1 = math.exp(float(np.clip(values[0], -12.0, 12.0)))
    t4 = math.exp(float(np.clip(values[3], -12.0, 12.0)))
    t2 = complex(float(values[1]), float(values[2]))
    tri = np.asarray([[t1, 0.0], [t2, t4]], dtype=np.complex128)
    rho = tri.conj().T @ tri
    return qt.Qobj(rho / float(np.real(np.trace(rho))), dims=[[2], [2]])


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
    """U = Rx(pi/3) * Rz(pi/4), applied to |g><g|."""
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
        # Plan Part 1 specifies m_k ~ N(E_k, 1/sqrt(N))
        sigma = 1.0 / math.sqrt(max(int(shots), 1))
        return float(np.clip(rng.normal(value, sigma), -1.0, 1.0)), None
    # Binomial mode (used for Part 2 / MLE)
    p_plus = float(np.clip(0.5 * (1.0 + value), EPS, 1.0 - EPS))
    counts_plus = int(rng.binomial(int(shots), p_plus))
    measured = (2.0 * counts_plus / float(shots)) - 1.0
    return float(measured), counts_plus


def single_qubit_dataset(
    rho_true: qt.Qobj,
    *,
    shots: int,
    rng: np.random.Generator,
    sample_mode: str = "gaussian",
) -> dict[str, Any]:
    """Gaussian noise as specified by plan Part 1: m_k ~ N(E_k, 1/sqrt(N))."""
    exact = branch_expectations_from_density(rho_true)
    measured: dict[str, float] = {}
    counts_plus: dict[str, int | None] = {}
    for axis in MEASUREMENT_AXES:
        m_val, cnt = sample_branch_measurement(
            exact[axis],
            shots=int(shots),
            rng=rng,
            mode=sample_mode,
        )
        measured[axis] = float(m_val)
        counts_plus[axis] = cnt
    return {
        "shots": int(shots),
        "exact_expectations": exact,
        "measured_expectations": measured,
        "counts_plus": counts_plus,
        "sample_mode": str(sample_mode),
    }


def fit_single_qubit_ls(dataset: dict[str, Any]) -> dict[str, Any]:
    """Least-squares fit: minimize sum_k (Tr(U_k rho U_k† sz) - m_k)^2."""
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
    }


def fit_single_qubit_mle(dataset: dict[str, Any]) -> dict[str, Any]:
    """Binomial MLE fit using Cholesky parameterization."""
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
                # Fall back to LS cost when binomial counts not available
                total += (pred[axis] - measured[axis]) ** 2
                continue
            prob = float(np.clip(0.5 * (1.0 + pred[axis]), EPS, 1.0 - EPS))
            total -= float(count) * math.log(prob)
            total -= float(shots - int(count)) * math.log(1.0 - prob)
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
    }


def run_single_qubit_baseline(
    *,
    shot_grid: Sequence[int] = DEFAULT_SHOT_GRID,
    repeats: int = 40,
    seed: int = 1234,
) -> dict[str, Any]:
    """Part 1: F vs N scaling curve. Noise is Gaussian per plan spec."""
    rho_true = single_qubit_target_state()
    rng = np.random.default_rng(int(seed))
    rows: list[dict[str, Any]] = []
    for shots in shot_grid:
        for trial in range(int(repeats)):
            dataset = single_qubit_dataset(rho_true, shots=int(shots), rng=rng, sample_mode="gaussian")
            fit_ls = fit_single_qubit_ls(dataset)
            fit_mle = fit_single_qubit_mle(dataset)
            rows.append({
                "shots": int(shots),
                "trial": int(trial),
                "fidelity_ls": density_matrix_fidelity(rho_true, fit_ls["rho"]),
                "fidelity_mle": density_matrix_fidelity(rho_true, fit_mle["rho"]),
                "trace_distance_ls": trace_distance(rho_true, fit_ls["rho"]),
                "trace_distance_mle": trace_distance(rho_true, fit_mle["rho"]),
            })
    summary_rows: list[dict[str, Any]] = []
    for shots in shot_grid:
        shot_rows = [r for r in rows if r["shots"] == shots]
        for method in ("ls", "mle"):
            vals = np.asarray([r[f"fidelity_{method}"] for r in shot_rows], dtype=float)
            summary_rows.append({
                "shots": int(shots),
                "method": str(method),
                "mean_fidelity": float(np.mean(vals)),
                "std_fidelity": float(np.std(vals)),
                "median_fidelity": float(np.median(vals)),
                "min_fidelity": float(np.min(vals)),
            })
    best_mle = max(
        (r for r in summary_rows if r["method"] == "mle"),
        key=lambda r: r["mean_fidelity"],
    )
    return {
        "true_state_bloch": dict(zip(("x", "y", "z"), bloch_from_density(rho_true))),
        "rows": rows,
        "summary_rows": summary_rows,
        "best_mle_mean_fidelity": float(best_mle["mean_fidelity"]),
    }


# ---------------------------------------------------------------------------
# Part 2A — Forward model and identifiability
# ---------------------------------------------------------------------------

def measurement_settings_wait_only(wait_grid_s: Sequence[float] = DEFAULT_WAIT_GRID_S) -> list[MeasurementSetting]:
    return [
        MeasurementSetting(protocol="wait_only", alpha=0.0 + 0.0j, wait_s=float(w), label=f"wait_{i:02d}")
        for i, w in enumerate(wait_grid_s)
    ]


def measurement_settings_displacement_only(alpha_grid: Sequence[complex] = DEFAULT_ALPHA_GRID) -> list[MeasurementSetting]:
    return [
        MeasurementSetting(protocol="displacement_only", alpha=complex(a), wait_s=0.0, label=f"disp_{i:02d}")
        for i, a in enumerate(alpha_grid)
    ]


def measurement_settings_combined(
    alpha_grid: Sequence[complex] = DEFAULT_COMBINED_ALPHA_GRID,
    wait_grid_s: Sequence[float] = DEFAULT_COMBINED_WAIT_GRID_S,
) -> list[MeasurementSetting]:
    settings = []
    counter = 0
    for alpha in alpha_grid:
        for wait_s in wait_grid_s:
            settings.append(MeasurementSetting(
                protocol="combined",
                alpha=complex(alpha),
                wait_s=float(wait_s),
                label=f"comb_{counter:03d}",
            ))
            counter += 1
    return settings


def measurement_settings_by_protocol() -> dict[str, list[MeasurementSetting]]:
    return {
        "wait_only": measurement_settings_wait_only(),
        "displacement_only": measurement_settings_displacement_only(),
        "combined": measurement_settings_combined(),
    }


def dispersive_phase_array(
    wait_s: float,
    *,
    n_cav: int,
    chi: float = CHI,
    chi_prime: float = CHI_PRIME,
) -> np.ndarray:
    """phi_m(t) = (chi*m + chi'*m*(m-1)) * t  for m = 0..n_cav-1."""
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
    """K_n(alpha, t) = sum_m |<m|D(alpha)|n>|^2 * exp(-i phi_m(t))."""
    if int(n_active) > int(kernel_dim):
        raise ValueError("kernel_dim must be >= n_active.")
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
    """X(a,t)+iY(a,t) = sum_n u_n K_n(a,t);  Z(a,t) = z_total."""
    u = np.asarray(weighted_transverse, dtype=np.complex128)
    x_vals, y_vals, z_vals = [], [], []
    for setting in settings:
        kernel = diagonal_kernel_vector(setting, n_active=u.size, kernel_dim=kernel_dim, chi=chi, chi_prime=chi_prime)
        transverse = np.sum(u * kernel)
        x_vals.append(float(np.real(transverse)))
        y_vals.append(float(np.imag(transverse)))
        z_vals.append(float(z_total))
    return {
        "X": np.asarray(x_vals, dtype=float),
        "Y": np.asarray(y_vals, dtype=float),
        "Z": np.asarray(z_vals, dtype=float),
    }


def diagonal_protocol_design_matrix(
    settings: Sequence[MeasurementSetting],
    *,
    n_active: int = N_ACTIVE,
    kernel_dim: int = KERNEL_DIM,
    chi: float = CHI,
    chi_prime: float = CHI_PRIME,
) -> np.ndarray:
    """3*len(settings) x (2*n_active+1) matrix mapping theta=[Re u, Im u, z_total] to [X, Y, Z]."""
    rows = []
    for setting in settings:
        kernel = diagonal_kernel_vector(setting, n_active=n_active, kernel_dim=kernel_dim, chi=chi, chi_prime=chi_prime)
        row_x = np.zeros(2 * int(n_active) + 1, dtype=float)
        row_y = np.zeros(2 * int(n_active) + 1, dtype=float)
        row_z = np.zeros(2 * int(n_active) + 1, dtype=float)
        row_x[:n_active] = np.real(kernel)
        row_x[n_active:2 * n_active] = -np.imag(kernel)
        row_y[:n_active] = np.imag(kernel)
        row_y[n_active:2 * n_active] = np.real(kernel)
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
    transverse = matrix[:, :2 * int(n_active)]
    svals = np.linalg.svd(matrix, compute_uv=False)
    svals_t = np.linalg.svd(transverse, compute_uv=False)
    cond_full = None if float(np.min(svals)) <= EPS else float(np.max(svals) / np.min(svals))
    cond_t = None if float(np.min(svals_t)) <= EPS else float(np.max(svals_t) / np.min(svals_t))
    return {
        "n_settings": int(len(settings)),
        "full_rank": int(np.linalg.matrix_rank(matrix)),
        "transverse_rank": int(np.linalg.matrix_rank(transverse)),
        "full_condition_number": cond_full,
        "transverse_condition_number": cond_t,
        "smallest_singular_value": float(np.min(svals)),
        "smallest_transverse_singular_value": float(np.min(svals_t)),
    }


# ---------------------------------------------------------------------------
# Part 2B — Full {p_n, rho_q^(n)} MLE — attempted first
# ---------------------------------------------------------------------------

def full_state_random_initial(n_active: int, rng: np.random.Generator) -> np.ndarray:
    logits = rng.normal(0.0, 0.8, size=int(n_active))
    raw_states = rng.normal(0.0, 0.5, size=4 * int(n_active))
    return np.concatenate([logits, raw_states])


def decode_full_state_params(raw: Sequence[float], *, n_active: int) -> tuple[np.ndarray, list[qt.Qobj]]:
    """Softmax populations + Cholesky qubit states."""
    arr = np.asarray(raw, dtype=float).reshape(5 * int(n_active))
    logits = arr[:n_active]
    probs = np.exp(logits - logsumexp(logits))
    states = []
    offset = n_active
    for _ in range(int(n_active)):
        states.append(cholesky_raw_to_density(arr[offset:offset + 4]))
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
        [complex(float(probs[i]) * float(bloch_from_density(states[i])[0]),
                 float(probs[i]) * float(bloch_from_density(states[i])[1]))
         for i in range(int(n_active))],
        dtype=np.complex128,
    )
    z_total = float(sum(float(probs[i]) * float(bloch_from_density(states[i])[2]) for i in range(int(n_active))))
    matrix = diagonal_protocol_design_matrix(settings, n_active=n_active, kernel_dim=kernel_dim)
    theta_lin = np.concatenate([np.real(weighted_transverse), np.imag(weighted_transverse), [z_total]])
    return matrix @ theta_lin


def run_full_state_mle_attempt(
    rows: Sequence[dict[str, Any]],
    settings: Sequence[MeasurementSetting],
    *,
    n_active: int = N_ACTIVE,
    kernel_dim: int = KERNEL_DIM,
    n_restarts: int = 8,
    seed: int = 4321,
) -> dict[str, Any]:
    """Attempt full {p_n, rho_q^(n)} recovery via MLE with random restarts.

    Returns results from all restarts, exposing the spread in recovered
    populations and qubit states. This characterises the non-uniqueness
    failure mode before the study pivots to the recoverable subspace.
    """
    rng = np.random.default_rng(int(seed))
    counts_plus, shots = _ordered_counts_vector(rows, settings)
    measured = _ordered_row_vector(rows, settings, "measured_expectation")
    bounds = [(-5.0, 5.0)] * int(n_active) + [(-8.0, 8.0), (-2.5, 2.5), (-2.5, 2.5), (-8.0, 8.0)] * int(n_active)

    def objective(theta: np.ndarray) -> float:
        pred = np.clip(
            full_state_prediction_vector(theta, settings, n_active=n_active, kernel_dim=kernel_dim),
            -0.999999, 0.999999,
        )
        prob = np.clip(0.5 * (1.0 + pred), EPS, 1.0 - EPS)
        return float(-np.sum(counts_plus * np.log(prob) + (shots - counts_plus) * np.log(1.0 - prob)))

    results = []
    for restart in range(int(n_restarts)):
        x0 = full_state_random_initial(n_active, rng)
        result = minimize(objective, x0, method="L-BFGS-B", bounds=bounds)
        if not result.success:
            result = minimize(objective, x0, method="Powell", options={"maxiter": 800})
        probs, states = decode_full_state_params(result.x, n_active=n_active)
        bloch_rows = []
        weighted = []
        per_sector_fidelity = []
        for i, (prob, state) in enumerate(zip(probs, states)):
            x, y, z = bloch_from_density(state)
            bloch_rows.append({"x": float(x), "y": float(y), "z": float(z)})
            weighted.append(complex(float(prob) * float(x), float(prob) * float(y)))
            per_sector_fidelity.append(None)  # No ground truth available here; filled by caller
        pred = full_state_prediction_vector(result.x, settings, n_active=n_active, kernel_dim=kernel_dim)
        residual = measured - pred
        results.append({
            "restart": int(restart),
            "objective": float(objective(result.x)),
            "residual_rms": float(np.sqrt(np.mean(residual ** 2))),
            "optimizer_success": bool(result.success),
            "probabilities": probs.tolist(),
            "bloch_rows": bloch_rows,
            "weighted_transverse": [complex_to_dict(v) for v in weighted],
            "z_total": float(sum(float(probs[i]) * bloch_rows[i]["z"] for i in range(int(n_active)))),
        })
    results.sort(key=lambda r: float(r["objective"]))
    probability_matrix = np.asarray([r["probabilities"] for r in results], dtype=float)
    z_matrix = np.asarray([[row["z"] for row in r["bloch_rows"]] for r in results], dtype=float)
    objectives = np.asarray([r["objective"] for r in results], dtype=float)
    return {
        "rows": results,
        "objective_span": float(np.max(objectives) - np.min(objectives)),
        "probability_std_by_sector": np.std(probability_matrix, axis=0).tolist(),
        "z_std_by_sector": np.std(z_matrix, axis=0).tolist(),
        "conclusion": (
            "Non-unique: multiple restarts yield different {p_n, Z_n} at similar objective values."
            if float(np.std(probability_matrix)) > 0.05
            else "Converged: restarts agree on {p_n, Z_n}."
        ),
    }


def construct_exact_gauge_family(
    weighted_transverse: Sequence[complex],
    z_total: float,
    *,
    patterns: Sequence[Sequence[float]] | None = None,
) -> list[dict[str, Any]]:
    """Algebraically construct distinct {p_n, Z_n} decompositions with identical observables."""
    u = np.asarray(weighted_transverse, dtype=np.complex128)
    min_prob = np.abs(u)
    slack_total = float(1.0 - np.sum(min_prob))
    if slack_total <= 1.0e-8:
        return []
    if patterns is None:
        patterns = ((4.0, 1.0, 3.0, 2.0), (1.0, 3.0, 2.0, 4.0))
    solutions = []
    for pattern in patterns:
        weights = normalize_probabilities(pattern[:u.size])
        probs = min_prob + slack_total * weights
        caps = np.sqrt(np.maximum(probs ** 2 - np.abs(u) ** 2, 0.0))
        cap_sum = float(np.sum(caps))
        if cap_sum <= EPS or abs(float(z_total)) > cap_sum + 1.0e-9:
            continue
        weighted_z = float(z_total) * caps / cap_sum
        z_values = np.divide(weighted_z, probs, out=np.zeros_like(weighted_z), where=probs > EPS)
        bloch_rows = []
        for i in range(u.size):
            xy = u[i] / probs[i]
            bloch_rows.append({"x": float(np.real(xy)), "y": float(np.imag(xy)), "z": float(z_values[i])})
        solutions.append({
            "pattern": [float(v) for v in pattern[:u.size]],
            "probabilities": probs.tolist(),
            "bloch_rows": bloch_rows,
            "reconstructed_z_total": float(np.sum(probs * z_values)),
        })
    return solutions


# ---------------------------------------------------------------------------
# Part 2C — Recoverable-subspace inference (weighted-transverse LS and MLE)
# ---------------------------------------------------------------------------

def _ordered_row_vector(
    rows: Sequence[dict[str, Any]],
    settings: Sequence[MeasurementSetting],
    field: str,
) -> np.ndarray:
    mapping = {(str(r["label"]), str(r["axis"]).upper()): r for r in rows}
    return np.asarray(
        [float(mapping[(s.label, axis)][field]) for s in settings for axis in MEASUREMENT_AXES],
        dtype=float,
    )


def _ordered_counts_vector(
    rows: Sequence[dict[str, Any]],
    settings: Sequence[MeasurementSetting],
) -> tuple[np.ndarray, np.ndarray]:
    mapping = {(str(r["label"]), str(r["axis"]).upper()): r for r in rows}
    counts, shots = [], []
    for s in settings:
        for axis in MEASUREMENT_AXES:
            row = mapping[(s.label, axis)]
            cnt = row.get("counts_plus")
            if cnt is None:
                meas = float(row["measured_expectation"])
                shot = int(row["shots"])
                cnt = int(round(0.5 * (1.0 + float(np.clip(meas, -1.0, 1.0))) * shot))
            counts.append(float(cnt))
            shots.append(float(row["shots"]))
    return np.asarray(counts, dtype=float), np.asarray(shots, dtype=float)


def unpack_weighted_parameter_vector(theta: Sequence[float], *, n_active: int) -> tuple[np.ndarray, float]:
    arr = np.asarray(theta, dtype=float).reshape(2 * int(n_active) + 1)
    return arr[:n_active] + 1.0j * arr[n_active:2 * n_active], float(arr[-1])


def weighted_parameter_vector(weighted_transverse: Sequence[complex], z_total: float) -> np.ndarray:
    u = np.asarray(weighted_transverse, dtype=np.complex128)
    return np.concatenate([np.real(u), np.imag(u), [float(z_total)]])


def infer_weighted_ls(
    rows: Sequence[dict[str, Any]],
    settings: Sequence[MeasurementSetting],
    *,
    n_active: int = N_ACTIVE,
    kernel_dim: int = KERNEL_DIM,
) -> dict[str, Any]:
    """Pseudoinverse LS solve for theta = [Re u_n, Im u_n, z_total]."""
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
    """Binomial MLE over the recoverable parameter vector theta."""
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
        "objective": float(objective(result.x)),
        "residual_rms": float(np.sqrt(np.mean(residual ** 2))),
        "residual_max_abs": float(np.max(np.abs(residual))),
        "optimizer_success": bool(result.success),
    }


# ---------------------------------------------------------------------------
# Part 2D — Per-sector metrics
# ---------------------------------------------------------------------------

def per_sector_metrics(
    inferred_weighted: Sequence[complex],
    true_sectors: Sequence[SectorSummary],
) -> list[dict[str, Any]]:
    """Per-Fock-sector accuracy metrics.

    For each sector n:
    - weighted transverse error: |u_n^inf - u_n^true|
    - oracle bloch vector error: ||v_n^inf - v_n^true|| (using true p_n to normalise)
    - oracle fidelity: F(rho_q^(n)_true, rho_q^(n)_oracle_inferred)
    - oracle trace distance

    'oracle' means we use the true p_n to normalise the inferred u_n, giving
    the best estimate of rho_q^(n) achievable from the inferred u_n.
    This is an oracle quantity — experimentally p_n is unknown.
    """
    u_inf = np.asarray(inferred_weighted, dtype=np.complex128)
    out = []
    for n, sector in enumerate(true_sectors):
        p_true = float(sector.population)
        u_true = complex(p_true * sector.x, p_true * sector.y)
        u_diff = u_inf[n] - u_true
        weighted_abs_error = float(abs(u_diff))

        entry: dict[str, Any] = {
            "level": int(sector.level),
            "p_true": float(p_true),
            "weighted_transverse_error": float(weighted_abs_error),
            "weighted_transverse_abs_inferred": float(abs(u_inf[n])),
            "weighted_transverse_abs_true": float(abs(u_true)),
        }

        if p_true > 1.0e-6:
            # Oracle: normalise by true population
            x_inf = float(np.real(u_inf[n])) / p_true
            y_inf = float(np.imag(u_inf[n])) / p_true
            # Z is not individually recoverable; use true Z_n as oracle stand-in
            z_inf = float(sector.z)
            bloch_inf = np.asarray([x_inf, y_inf, z_inf], dtype=float)
            bloch_true = sector.bloch_vector

            # Clip to Bloch sphere
            norm_inf = float(np.linalg.norm(bloch_inf[:2]))  # transverse norm
            if norm_inf > 1.0:
                bloch_inf[:2] /= norm_inf

            bloch_error = float(np.linalg.norm(bloch_inf - bloch_true))
            rho_inf = qubit_density_from_bloch(float(bloch_inf[0]), float(bloch_inf[1]), float(bloch_inf[2]))
            rho_true = qubit_density_from_bloch(float(bloch_true[0]), float(bloch_true[1]), float(bloch_true[2]))
            entry["oracle_bloch_error"] = float(bloch_error)
            entry["oracle_fidelity"] = density_matrix_fidelity(rho_true, rho_inf)
            entry["oracle_trace_distance"] = trace_distance(rho_true, rho_inf)
            entry["oracle_x_error"] = float(x_inf - sector.x)
            entry["oracle_y_error"] = float(y_inf - sector.y)
        else:
            entry["oracle_bloch_error"] = None
            entry["oracle_fidelity"] = None
            entry["oracle_trace_distance"] = None
            entry["oracle_x_error"] = None
            entry["oracle_y_error"] = None

        out.append(entry)
    return out


def recoverable_error_summary(
    inferred_weighted: Sequence[complex],
    true_weighted: Sequence[complex],
    *,
    inferred_z_total: float,
    true_z_total: float,
    true_sectors: Sequence[SectorSummary] | None = None,
) -> dict[str, Any]:
    u_inf = np.asarray(inferred_weighted, dtype=np.complex128)
    u_true = np.asarray(true_weighted, dtype=np.complex128)
    errors = u_inf - u_true
    payload: dict[str, Any] = {
        "per_sector_abs_error": np.abs(errors).tolist(),
        "mean_abs_error": float(np.mean(np.abs(errors))),
        "max_abs_error": float(np.max(np.abs(errors))),
        "weighted_rmse": float(np.sqrt(np.mean(np.real(errors) ** 2 + np.imag(errors) ** 2))),
        "z_total_error": float(float(inferred_z_total) - float(true_z_total)),
    }
    if true_sectors is not None:
        sector_metrics = per_sector_metrics(inferred_weighted, true_sectors)
        payload["per_sector_metrics"] = sector_metrics
        oracle_fids = [m["oracle_fidelity"] for m in sector_metrics if m["oracle_fidelity"] is not None]
        oracle_tds = [m["oracle_trace_distance"] for m in sector_metrics if m["oracle_trace_distance"] is not None]
        oracle_bloch = [m["oracle_bloch_error"] for m in sector_metrics if m["oracle_bloch_error"] is not None]
        payload["mean_oracle_fidelity"] = float(np.mean(oracle_fids)) if oracle_fids else None
        payload["mean_oracle_trace_distance"] = float(np.mean(oracle_tds)) if oracle_tds else None
        payload["mean_oracle_bloch_error"] = float(np.mean(oracle_bloch)) if oracle_bloch else None
    return payload


# ---------------------------------------------------------------------------
# Part 2E — Black-box case library
# ---------------------------------------------------------------------------

def padded_rotation_arrays(*, n_cav: int) -> tuple[np.ndarray, np.ndarray]:
    theta = np.zeros(int(n_cav), dtype=float)
    phi = np.zeros(int(n_cav), dtype=float)
    theta[:N_ACTIVE] = np.asarray(DEFAULT_THETA, dtype=float)
    phi[:N_ACTIVE] = np.asarray(DEFAULT_PHI, dtype=float)
    return theta, phi


def padded_phase_array(*, n_cav: int) -> np.ndarray:
    phases = np.zeros(int(n_cav), dtype=float)
    phases[:N_ACTIVE] = np.asarray(DEFAULT_CPSQR_PHASES, dtype=float)
    return phases


def cavity_population_vector(n_cav: int, probabilities: Sequence[float] = DEFAULT_POPULATIONS) -> np.ndarray:
    probs = normalize_probabilities(probabilities)
    padded = np.zeros(int(n_cav), dtype=float)
    padded[:probs.size] = probs
    return padded


def build_ideal_reference_cases(*, n_cav: int = 6) -> list[StudyCase]:
    """Analytic reference cases covering all planned model types.

    Includes:
    - Ideal SQR on diagonal cavity, all 4 qubit inputs (g, e, +, +y) -> Model A
    - Ideal SQR on cavity |0> and |1> separately -> Model A, single-sector
    - Ideal SQR on D(beta)|0> coherent input -> Model A (approx diagonal)
    - Ideal SQR on superposition cavity -> Model B (off-diagonal blocks)
    - CPSQR-like conditional phase gate on |+> qubit -> Model A
    """
    model = build_model(n_cav=n_cav, n_tr=2)
    frame = build_frame(model)
    theta, phi = padded_rotation_arrays(n_cav=n_cav)
    sqr = ideal_sqr_operator(theta, phi)
    phases = padded_phase_array(n_cav=n_cav)
    cpsqr = cpsqr_like_operator(phases)
    populations = cavity_population_vector(n_cav)
    cavity_diag = cavity_mixture_state(populations)
    cavity_coherent = cavity_superposition_state(np.sqrt(populations[:n_cav].astype(np.complex128)))
    cases: list[StudyCase] = []

    # Model A: all four qubit inputs on diagonal cavity
    for qubit_label in ("g", "e", "+", "+y"):
        rho_in = qt.tensor(qubit_input_state(qubit_label), cavity_diag)
        rho_out = sqr * rho_in * sqr.dag()
        cases.append(StudyCase(
            case_id=f"ideal_diag_{qubit_label}",
            family="analytic",
            description=f"Ideal block-diagonal SQR on diagonal cavity, qubit |{qubit_label}>.",
            model_class="fock_diagonal",
            state=rho_out,
            model=model,
            frame=frame,
            truth_sectors=tuple(sector_summaries_from_state(rho_out, range(N_ACTIVE))),
            metadata={"gate": "ideal_sqr", "qubit_input": qubit_label, "cavity_input": "mixture"},
        ))

    # Model A: cavity |0> only (trivial single-sector case)
    cavity_0 = cavity_mixture_state([1.0] + [0.0] * (n_cav - 1))
    rho_in_0 = qt.tensor(qubit_ground_dm(), cavity_0)
    rho_out_0 = sqr * rho_in_0 * sqr.dag()
    cases.append(StudyCase(
        case_id="ideal_diag_g_cav0",
        family="analytic",
        description="Ideal SQR on cavity |0>, qubit |g>. Single active sector.",
        model_class="fock_diagonal",
        state=rho_out_0,
        model=model,
        frame=frame,
        truth_sectors=tuple(sector_summaries_from_state(rho_out_0, range(N_ACTIVE))),
        metadata={"gate": "ideal_sqr", "qubit_input": "g", "cavity_input": "fock_0"},
    ))

    # Model A: cavity |1> only
    cavity_1 = cavity_mixture_state([0.0, 1.0] + [0.0] * (n_cav - 2))
    rho_in_1 = qt.tensor(qubit_ground_dm(), cavity_1)
    rho_out_1 = sqr * rho_in_1 * sqr.dag()
    cases.append(StudyCase(
        case_id="ideal_diag_g_cav1",
        family="analytic",
        description="Ideal SQR on cavity |1>, qubit |g>. Single active sector.",
        model_class="fock_diagonal",
        state=rho_out_1,
        model=model,
        frame=frame,
        truth_sectors=tuple(sector_summaries_from_state(rho_out_1, range(N_ACTIVE))),
        metadata={"gate": "ideal_sqr", "qubit_input": "g", "cavity_input": "fock_1"},
    ))

    # Model A: displaced vacuum D(beta)|0> as initial cavity state
    beta = 0.8 + 0.0j
    cav_coherent_input = as_dm(coherent_state(beta, n_cav=n_cav))
    rho_in_coh = qt.tensor(qubit_ground_dm(), cav_coherent_input)
    rho_out_coh = sqr * rho_in_coh * sqr.dag()
    cases.append(StudyCase(
        case_id="ideal_diag_g_coherent_input",
        family="analytic",
        description="Ideal SQR on D(0.8)|0> coherent cavity input, qubit |g>.",
        model_class="fock_diagonal",
        state=rho_out_coh,
        model=model,
        frame=frame,
        truth_sectors=tuple(sector_summaries_from_state(rho_out_coh, range(N_ACTIVE))),
        metadata={"gate": "ideal_sqr", "qubit_input": "g", "cavity_input": "coherent_0.8"},
    ))

    # Model B: cavity superposition -> off-diagonal blocks (coherence witness target)
    rho_coherent_in = qt.tensor(qubit_ground_dm(), as_dm(cavity_coherent))
    rho_coherent_out = sqr * rho_coherent_in * sqr.dag()
    cases.append(StudyCase(
        case_id="ideal_coherent_g",
        family="analytic",
        description="Ideal SQR on cavity superposition, qubit |g>. Same marginals as mixture but with off-diagonal cavity blocks.",
        model_class="coherent_block",
        state=rho_coherent_out,
        model=model,
        frame=frame,
        truth_sectors=tuple(sector_summaries_from_state(rho_coherent_out, range(N_ACTIVE))),
        metadata={"gate": "ideal_sqr", "qubit_input": "g", "cavity_input": "superposition"},
    ))

    # Model A: CPSQR-like conditional phase gate
    rho_cpsqr_in = qt.tensor(qubit_plus_dm(), cavity_diag)
    rho_cpsqr_out = cpsqr * rho_cpsqr_in * cpsqr.dag()
    cases.append(StudyCase(
        case_id="cpsqr_like_plus",
        family="analytic",
        description="CPSQR-like conditional phase gate on diagonal cavity, qubit |+>.",
        model_class="fock_diagonal",
        state=rho_cpsqr_out,
        model=model,
        frame=frame,
        truth_sectors=tuple(sector_summaries_from_state(rho_cpsqr_out, range(N_ACTIVE))),
        metadata={"gate": "cpsqr_like", "qubit_input": "+", "cavity_input": "mixture"},
    ))

    return cases


def _make_joint_state(
    qubit_label: str,
    cavity_state: qt.Qobj,
    gate_op: qt.Qobj,
) -> qt.Qobj:
    rho_in = qt.tensor(qubit_input_state(qubit_label), as_dm(cavity_state))
    return gate_op * rho_in * gate_op.dag()


def build_pulse_level_cases(
    *,
    n_cav: int = 6,
    n_tr_standard: int = 2,
    n_tr_leakage: int = 3,
    waveform_artifacts: dict[str, Any] | None = None,
    seed: int = 999,
) -> tuple[list[StudyCase], dict[str, Any]]:
    """Build all planned pulse-level black-box cases.

    Returns (cases, waveform_metadata_dict) where the metadata contains the
    optimized corrections for artifact saving.

    Cases:
    1. long_pulse_optimized_mix_g   : near-ideal SQR, qubit |g>, cavity mixture
    2. long_pulse_optimized_mix_+   : near-ideal SQR, qubit |+>, cavity mixture
    3. short_pulse_seed_mix_g       : shorter / distorted SQR (seed only, no optimize)
    4. short_pulse_seed_superpos_g  : shorter pulse on superposition cavity -> coherent output
    5. cpsqr_pulse_mix_+            : analytic CPSQR-like gate (waveform-free stand-in)
    6. leakage_replay_mix_g         : long pulse optimized on n_tr=2 but replayed on n_tr=3 -> leakage
    7. spectator_distorted_mix_g    : long pulse with extra detuning (non-ideal spectator simulation)
    """
    model_std = build_model(n_cav=n_cav, n_tr=n_tr_standard)
    model_leak = build_model(n_cav=n_cav, n_tr=n_tr_leakage)
    frame_std = build_frame(model_std)
    frame_leak = build_frame(model_leak)

    populations = cavity_population_vector(n_cav)
    cavity_diag = cavity_mixture_state(populations)
    cavity_super = cavity_superposition_state(np.sqrt(populations[:n_cav].astype(np.complex128)))
    phases = padded_phase_array(n_cav=n_cav)
    cpsqr_op = cpsqr_like_operator(phases)

    wf_meta: dict[str, Any] = {}
    cases: list[StudyCase] = []

    # --- LONG PULSE (optimized) ---
    opt_cfg = ConditionedOptimizationConfig(
        active_levels=tuple(range(N_ACTIVE)),
        parameters=("d_lambda", "d_alpha", "d_omega"),
        maxiter_stage1=6,
        maxiter_stage2=8,
        d_lambda_bounds=(-0.75, 0.75),
        d_alpha_bounds=(-np.pi, np.pi),
        d_omega_hz_bounds=(-3.0e6, 3.0e6),
        regularization_lambda=5.0e-4,
        regularization_alpha=5.0e-4,
        regularization_omega=5.0e-4,
    )
    targets = conditioned_multitone_targets()
    run_config_long = conditioned_run_config(model_std, duration_s=LONG_DURATION_S)

    optimization = optimize_conditioned_multitone(
        model_std,
        targets,
        run_config_long,
        optimization_config=opt_cfg,
        simulation_mode="reduced",
        label="long_pulse_near_ideal",
    )
    optimized = optimization.optimized_result
    compiled_long = optimized.compiled
    drive_ops_long = {str(optimized.waveform.drive_channel): "qubit"}
    wf_meta["long_optimized"] = True

    for qubit_label in ("g", "+"):
        initial_state = qt.tensor(qubit_input_state(qubit_label), cavity_diag)
        out_states = simulate_compiled_on_states(
            model_std, compiled_long,
            frame=frame_std,
            drive_ops=drive_ops_long,
            initial_states=[initial_state],
        )
        rho_out = out_states[0]
        cases.append(StudyCase(
            case_id=f"long_pulse_optimized_mix_{qubit_label}",
            family="pulse_level",
            description=f"Near-ideal optimized multitone SQR ({LONG_DURATION_S*1e6:.1f} µs), qubit |{qubit_label}>, diagonal cavity.",
            model_class="pulse_level",
            state=rho_out,
            model=model_std,
            frame=frame_std,
            truth_sectors=tuple(sector_summaries_from_state(rho_out, range(N_ACTIVE))),
            metadata={"duration_s": LONG_DURATION_S, "qubit_input": qubit_label, "pulse_type": "near_ideal"},
        ))

    # --- SHORT PULSE (seed, no optimization = distorted SQR) ---
    wf_short, _, run_config_short = build_multitone_waveform(model_std, duration_s=SHORT_DURATION_S)
    compiled_short = compile_waveform(wf_short, run_config_short)
    drive_ops_short = {str(wf_short.drive_channel): "qubit"}
    wf_meta["short_seed_waveform"] = True  # flag for artifact saving

    for cavity_label, cavity_state in [("mixture", cavity_diag), ("superpos", as_dm(cavity_super))]:
        initial_state = qt.tensor(qubit_ground_dm(), cavity_state)
        out_states = simulate_compiled_on_states(
            model_std, compiled_short,
            frame=frame_std,
            drive_ops=drive_ops_short,
            initial_states=[initial_state],
        )
        rho_out = out_states[0]
        model_class = "pulse_level" if cavity_label == "mixture" else "coherent_block"
        cases.append(StudyCase(
            case_id=f"short_pulse_seed_{cavity_label}_g",
            family="pulse_level",
            description=f"Shorter unoptimized multitone SQR ({SHORT_DURATION_S*1e6:.2f} µs), qubit |g>, cavity {cavity_label}.",
            model_class=model_class,
            state=rho_out,
            model=model_std,
            frame=frame_std,
            truth_sectors=tuple(sector_summaries_from_state(rho_out, range(N_ACTIVE))),
            metadata={"duration_s": SHORT_DURATION_S, "qubit_input": "g", "pulse_type": "distorted"},
        ))

    # --- CPSQR-like (analytic operator, stand-in for waveform case) ---
    rho_cpsqr = _make_joint_state("+", cavity_diag, cpsqr_op)
    cases.append(StudyCase(
        case_id="cpsqr_pulse_mix_+",
        family="analytic_pulse_standin",
        description="CPSQR-like conditional phase gate on diagonal cavity, qubit |+>. Analytic operator stand-in.",
        model_class="fock_diagonal",
        state=rho_cpsqr,
        model=model_std,
        frame=frame_std,
        truth_sectors=tuple(sector_summaries_from_state(rho_cpsqr, range(N_ACTIVE))),
        metadata={"gate": "cpsqr_analytic", "qubit_input": "+", "pulse_type": "cpsqr"},
    ))

    # --- LEAKAGE: long optimized pulse replayed on n_tr=3 model ---
    wf_long_leak, _, run_config_long_leak = build_multitone_waveform(model_leak, duration_s=LONG_DURATION_S)
    compiled_long_leak = compile_waveform(wf_long_leak, run_config_long_leak)
    drive_ops_leak = {str(wf_long_leak.drive_channel): "qubit"}
    initial_g_leak = qt.tensor(qt.basis(n_tr_leakage, 0).proj(), cavity_diag)
    out_leak = simulate_compiled_on_states(
        model_leak, compiled_long_leak,
        frame=frame_leak,
        drive_ops=drive_ops_leak,
        initial_states=[initial_g_leak],
    )
    rho_out_leak = out_leak[0]
    cases.append(StudyCase(
        case_id="leakage_replay_mix_g",
        family="pulse_level",
        description="Near-ideal pulse replayed on n_tr=3 model; leakage into |f> state.",
        model_class="pulse_level_leakage",
        state=rho_out_leak,
        model=model_leak,
        frame=frame_leak,
        truth_sectors=tuple(sector_summaries_from_state(rho_out_leak, range(N_ACTIVE))),
        metadata={"duration_s": LONG_DURATION_S, "qubit_input": "g", "pulse_type": "leakage", "n_tr": n_tr_leakage},
    ))

    # --- NON-IDEAL SPECTATOR: long pulse with frequency offset (detuned drive) ---
    # Simulate by running the long waveform with a deliberately wrong frame
    # (offset omega_q_frame by chi/4 to mimic an off-resonance spectator qubit)
    spectator_frame = FrameSpec(
        omega_q_frame=float(model_std.omega_q) + float(abs(CHI)) * 0.25,
        omega_c_frame=float(model_std.omega_c),
    )
    out_spectator = simulate_compiled_on_states(
        model_std, compiled_long,
        frame=spectator_frame,
        drive_ops=drive_ops_long,
        initial_states=[qt.tensor(qubit_ground_dm(), cavity_diag)],
    )
    rho_spectator = out_spectator[0]
    cases.append(StudyCase(
        case_id="spectator_distorted_mix_g",
        family="pulse_level",
        description="Long pulse run with 0.25 chi frequency offset to simulate non-ideal spectator.",
        model_class="pulse_level",
        state=rho_spectator,
        model=model_std,
        frame=frame_std,
        truth_sectors=tuple(sector_summaries_from_state(rho_spectator, range(N_ACTIVE))),
        metadata={"duration_s": LONG_DURATION_S, "qubit_input": "g", "pulse_type": "spectator"},
    ))

    return cases, wf_meta


# ---------------------------------------------------------------------------
# Part 2F — Model-B coherence sweep
# ---------------------------------------------------------------------------

def build_coherence_sweep_states(
    *,
    n_cav: int = 6,
    fractions: Sequence[float] = DEFAULT_COHERENCE_FRACTIONS,
) -> list[dict[str, Any]]:
    """Mix a Fock-diagonal state with a coherent-cavity state at varying fractions.

    rho_mix(f) = (1-f) * rho_diagonal + f * rho_coherent_block

    Both rho_diagonal and rho_coherent_block share the same cavity sector
    marginals so that the qubit-level effect is purely from cavity coherences.
    Returns a list of StudyCase-like dicts plus the mixing fraction.
    """
    model = build_model(n_cav=n_cav, n_tr=2)
    frame = build_frame(model)
    theta, phi = padded_rotation_arrays(n_cav=n_cav)
    sqr = ideal_sqr_operator(theta, phi)
    populations = cavity_population_vector(n_cav)
    cavity_diag = cavity_mixture_state(populations)
    cavity_coherent = as_dm(cavity_superposition_state(np.sqrt(populations[:n_cav].astype(np.complex128))))

    rho_in_diag = qt.tensor(qubit_ground_dm(), cavity_diag)
    rho_in_coh = qt.tensor(qubit_ground_dm(), cavity_coherent)
    rho_out_diag = sqr * rho_in_diag * sqr.dag()
    rho_out_coh = sqr * rho_in_coh * sqr.dag()

    sweep_entries = []
    for f in fractions:
        f = float(f)
        rho_mix = (1.0 - f) * rho_out_diag + f * rho_out_coh
        # Re-normalise (should already be normalised but guard against float errors)
        rho_mix = rho_mix / float(rho_mix.tr())
        sweep_entries.append({
            "coherence_fraction": f,
            "state": rho_mix,
            "model": model,
            "frame": frame,
            "truth_sectors": tuple(sector_summaries_from_state(rho_mix, range(N_ACTIVE))),
        })
    return sweep_entries


# ---------------------------------------------------------------------------
# Measurement data generation (shared)
# ---------------------------------------------------------------------------

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
    rows = []
    for setting in settings:
        alpha_actual = complex(setting.alpha) * complex(alpha_scale * np.exp(1.0j * alpha_phase))
        evolved = analysis_evolved_state(state, model=model, frame=frame, alpha=alpha_actual, wait_s=float(setting.wait_s))
        exact = branch_expectations_from_joint_state(evolved, n_tr=int(model.n_tr), angle_errors=angle_errors)
        for axis in MEASUREMENT_AXES:
            mode = "gaussian" if int(model.n_tr) > 2 or str(noise.sample_mode) == "gaussian" else "binomial"
            measured, counts_plus = sample_branch_measurement(
                exact[axis], shots=int(noise.shots), rng=rng, mode=mode,
            )
            rows.append({
                **setting.as_dict(),
                "axis": str(axis),
                "exact_expectation": float(exact[axis]),
                "measured_expectation": float(measured),
                "counts_plus": counts_plus,
                "shots": int(noise.shots),
            })
    return {"rows": rows, "noise_model": noise.as_dict()}


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
    angle_errors = {axis: float(rng.normal(0.0, float(noise.rotation_sigma_rad))) for axis in MEASUREMENT_AXES}
    rows = []
    for setting in settings:
        alpha_actual = complex(setting.alpha) * complex(alpha_scale * np.exp(1.0j * alpha_phase))
        setting_actual = MeasurementSetting(
            protocol=setting.protocol, alpha=alpha_actual, wait_s=float(setting.wait_s), label=setting.label,
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
        rho_meas = qt.Qobj(
            project_density_to_physical(
                np.asarray(qubit_density_from_bloch(bloch["X"], bloch["Y"], bloch["Z"]).full(), dtype=np.complex128)
            ),
            dims=[[2], [2]],
        )
        branch_exact = branch_expectations_from_density(rho_meas, angle_errors=angle_errors)
        for axis in MEASUREMENT_AXES:
            measured, counts_plus = sample_branch_measurement(
                branch_exact[axis], shots=int(noise.shots), rng=rng, mode=noise.sample_mode,
            )
            rows.append({
                **setting.as_dict(),
                "axis": str(axis),
                "exact_expectation": float(branch_exact[axis]),
                "measured_expectation": float(measured),
                "counts_plus": counts_plus,
                "shots": int(noise.shots),
            })
    return {"rows": rows, "noise_model": noise.as_dict()}
