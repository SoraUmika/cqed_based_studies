"""Local reduced multitone effective-unitary metric for the corrected SQR study."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np
import qutip as qt
from scipy.optimize import Bounds, minimize

from common import wrap_pi
from cqed_sim.calibration import (
    ConditionedMultitoneCorrections,
    ConditionedOptimizationConfig,
    ConditionedQubitTargets,
)
from cqed_sim.calibration.conditioned_multitone import (
    _classify_error,
    build_conditioned_multitone_tones,
    build_conditioned_multitone_waveform,
    compile_conditioned_multitone_waveform,
    qubit_density_matrix_from_angles,
)
from cqed_sim.core.frequencies import manifold_transition_frequency
from cqed_sim.pulses.calibration import sqr_lambda0_rad_s
from cqed_sim.pulses.envelopes import normalized_gaussian


_SIGMA_PLUS = qt.create(2)
_SIGMA_MINUS = qt.destroy(2)
_SIGMA_X = qt.sigmax()
_SIGMA_Y = qt.sigmay()
_SIGMA_Z = qt.sigmaz()
_N_Q = qt.num(2)


def _normalize_unitary(matrix: np.ndarray) -> np.ndarray:
    out = np.asarray(matrix, dtype=np.complex128)
    det = np.linalg.det(out)
    if abs(det) > 1.0e-15:
        out = out * np.exp(-0.5j * np.angle(det))
    return out


def _bloch_from_ket(ket: np.ndarray) -> tuple[float, float, float]:
    vec = np.asarray(ket, dtype=np.complex128).reshape((2, 1))
    rho = vec @ vec.conj().T
    return (
        2.0 * float(np.real(rho[0, 1])),
        2.0 * float(np.imag(rho[0, 1])),
        float(np.real(rho[0, 0] - rho[1, 1])),
    )


def _unitary_rotation_parameters(unitary: np.ndarray) -> tuple[float, float, float]:
    u = _normalize_unitary(np.asarray(unitary, dtype=np.complex128))
    cos_half = float(np.clip(np.real(np.trace(u) / 2.0), -1.0, 1.0))
    theta = float(2.0 * np.arccos(cos_half))
    if theta < 1.0e-12:
        return 0.0, 0.0, 0.0
    sin_half = float(np.sin(theta / 2.0))
    nx = float(np.real(1.0j * np.trace(_SIGMA_X.full() @ u) / (2.0 * sin_half)))
    ny = float(np.real(1.0j * np.trace(_SIGMA_Y.full() @ u) / (2.0 * sin_half)))
    nz = float(np.real(1.0j * np.trace(_SIGMA_Z.full() @ u) / (2.0 * sin_half)))
    norm = float(np.sqrt(nx * nx + ny * ny + nz * nz))
    if norm > 1.0e-12:
        nx /= norm
        ny /= norm
        nz /= norm
    phi = float(np.mod(np.arctan2(ny, nx), 2.0 * np.pi))
    return theta, phi, nz


def _process_fidelity(target_unitary: np.ndarray, simulated_unitary: np.ndarray) -> float:
    target = np.asarray(target_unitary, dtype=np.complex128)
    simulated = np.asarray(simulated_unitary, dtype=np.complex128)
    overlap = np.trace(target.conj().T @ simulated)
    return float(np.clip(abs(overlap) ** 2 / 4.0, 0.0, 1.0))


def _target_unitary(theta: float, phi: float) -> np.ndarray:
    return np.asarray(
        (
            np.cos(theta / 2.0) * np.eye(2, dtype=np.complex128)
            - 1j
            * np.sin(theta / 2.0)
            * (np.cos(phi) * _SIGMA_X.full() + np.sin(phi) * _SIGMA_Y.full())
        ),
        dtype=np.complex128,
    )


def _reduced_sector_unitary(model, compiled, waveform, run_config, n: int) -> np.ndarray:
    coeff = np.asarray(compiled.channels[waveform.drive_channel].distorted, dtype=np.complex128)
    detuning = float(manifold_transition_frequency(model, int(n), frame=run_config.frame))
    hamiltonian = [
        detuning * _N_Q,
        [_SIGMA_PLUS, coeff],
        [_SIGMA_MINUS, np.conj(coeff)],
    ]
    options: dict[str, Any] = {
        "atol": 1.0e-8,
        "rtol": 1.0e-7,
        "nsteps": 200000,
    }
    if run_config.max_step_s is not None:
        options["max_step"] = float(run_config.max_step_s)
    propagators = qt.propagator(
        hamiltonian,
        compiled.tlist,
        options=options,
        tlist=compiled.tlist,
    )
    final = propagators[-1] if isinstance(propagators, list) else propagators
    return _normalize_unitary(np.asarray(final.full(), dtype=np.complex128))


@dataclass(frozen=True)
class ReducedUnitarySectorMetrics:
    n: int
    weight: float
    process_fidelity: float
    state_fidelity: float
    target_theta_rad: float
    target_phi_rad: float
    achieved_theta_rad: float
    achieved_phi_rad: float
    achieved_axis_z: float
    theta_error_rad: float
    phi_error_rad: float
    bloch_x: float
    bloch_y: float
    bloch_z: float
    dominant_error: str

    def as_dict(self) -> dict[str, float | int | str]:
        return {
            "n": int(self.n),
            "weight": float(self.weight),
            "process_fidelity": float(self.process_fidelity),
            "state_fidelity": float(self.state_fidelity),
            "target_theta_rad": float(self.target_theta_rad),
            "target_phi_rad": float(self.target_phi_rad),
            "achieved_theta_rad": float(self.achieved_theta_rad),
            "achieved_phi_rad": float(self.achieved_phi_rad),
            "achieved_axis_z": float(self.achieved_axis_z),
            "theta_error_rad": float(self.theta_error_rad),
            "phi_error_rad": float(self.phi_error_rad),
            "bloch_x": float(self.bloch_x),
            "bloch_y": float(self.bloch_y),
            "bloch_z": float(self.bloch_z),
            "dominant_error": str(self.dominant_error),
        }


@dataclass(frozen=True)
class ReducedUnitaryValidationResult:
    targets: ConditionedQubitTargets
    corrections: ConditionedMultitoneCorrections
    waveform: Any
    compiled: Any
    sector_metrics: tuple[ReducedUnitarySectorMetrics, ...]
    aggregate_cost: float
    weighted_mean_process_fidelity: float
    weighted_mean_state_fidelity: float
    metadata: dict[str, Any]

    def as_dict(self) -> dict[str, Any]:
        return {
            "aggregate_cost": float(self.aggregate_cost),
            "weighted_mean_process_fidelity": float(self.weighted_mean_process_fidelity),
            "weighted_mean_state_fidelity": float(self.weighted_mean_state_fidelity),
            "sector_metrics": [metric.as_dict() for metric in self.sector_metrics],
            "corrections": {
                "d_lambda": [float(x) for x in self.corrections.d_lambda],
                "d_alpha": [float(x) for x in self.corrections.d_alpha],
                "d_omega_rad_s": [float(x) for x in self.corrections.d_omega_rad_s],
            },
            "metadata": self.metadata,
        }


@dataclass(frozen=True)
class ReducedUnitaryOptimizationResult:
    initial_result: ReducedUnitaryValidationResult
    optimized_result: ReducedUnitaryValidationResult
    optimized_corrections: ConditionedMultitoneCorrections
    active_levels: tuple[int, ...]
    parameters: tuple[str, ...]
    history: tuple[dict[str, float], ...]
    success_stage1: bool
    success_stage2: bool
    message_stage1: str
    message_stage2: str
    start_label: str

    def improvement_summary(self) -> dict[str, Any]:
        return {
            "start_label": str(self.start_label),
            "initial_cost": float(self.initial_result.aggregate_cost),
            "optimized_cost": float(self.optimized_result.aggregate_cost),
            "cost_reduction": float(self.initial_result.aggregate_cost - self.optimized_result.aggregate_cost),
            "initial_weighted_mean_process_fidelity": float(self.initial_result.weighted_mean_process_fidelity),
            "optimized_weighted_mean_process_fidelity": float(self.optimized_result.weighted_mean_process_fidelity),
            "initial_weighted_mean_state_fidelity": float(self.initial_result.weighted_mean_state_fidelity),
            "optimized_weighted_mean_state_fidelity": float(self.optimized_result.weighted_mean_state_fidelity),
            "active_levels": [int(level) for level in self.active_levels],
            "parameters": [str(name) for name in self.parameters],
        }


def run_reduced_unitary_validation(
    model,
    targets: ConditionedQubitTargets,
    run_config,
    *,
    corrections: ConditionedMultitoneCorrections | None = None,
) -> ReducedUnitaryValidationResult:
    corr = ConditionedMultitoneCorrections.zeros(targets.n_levels) if corrections is None else corrections.padded(targets.n_levels)
    tones = build_conditioned_multitone_tones(model, targets, run_config, corrections=corr)
    waveform = build_conditioned_multitone_waveform(tones, run_config, label="reduced_unitary")
    compiled = compile_conditioned_multitone_waveform(waveform, run_config)
    sector_metrics: list[ReducedUnitarySectorMetrics] = []
    for n in range(targets.n_levels):
        simulated = _reduced_sector_unitary(model, compiled, waveform, run_config, n)
        target = _target_unitary(float(targets.theta[n]), float(targets.phi[n]))
        process_fidelity = _process_fidelity(target, simulated)
        achieved_theta, achieved_phi, achieved_axis_z = _unitary_rotation_parameters(simulated)
        theta_error = wrap_pi(achieved_theta - float(targets.theta[n]))
        phi_error = wrap_pi(achieved_phi - float(np.mod(targets.phi[n], 2.0 * np.pi)))
        ket = np.asarray(simulated[:, 0], dtype=np.complex128).reshape((2, 1))
        x, y, z = _bloch_from_ket(ket)
        target_dm = qubit_density_matrix_from_angles(float(targets.theta[n]), float(targets.phi[n]))
        rho = qt.Qobj(ket @ ket.conj().T, dims=[[2], [2]])
        state_fidelity = float(np.clip(np.real((target_dm * rho).tr()), 0.0, 1.0))
        dominant = _classify_error(theta_error, phi_error, 1.0)
        sector_metrics.append(
            ReducedUnitarySectorMetrics(
                n=int(n),
                weight=float(targets.weights[n]),
                process_fidelity=float(process_fidelity),
                state_fidelity=float(state_fidelity),
                target_theta_rad=float(targets.theta[n]),
                target_phi_rad=float(np.mod(targets.phi[n], 2.0 * np.pi)),
                achieved_theta_rad=float(achieved_theta),
                achieved_phi_rad=float(achieved_phi),
                achieved_axis_z=float(achieved_axis_z),
                theta_error_rad=float(theta_error),
                phi_error_rad=float(phi_error),
                bloch_x=float(x),
                bloch_y=float(y),
                bloch_z=float(z),
                dominant_error=str(dominant),
            )
        )
    weights = np.asarray(targets.weights, dtype=float)
    process_values = np.asarray([metric.process_fidelity for metric in sector_metrics], dtype=float)
    state_values = np.asarray([metric.state_fidelity for metric in sector_metrics], dtype=float)
    aggregate_cost = float(np.sum(weights * (1.0 - process_values)))
    metadata = {
        "t_s": np.asarray(compiled.tlist, dtype=float).tolist(),
        "tone_specs": [tone.as_dict() for tone in tones],
        "waveform_metadata": dict(waveform.metadata),
    }
    return ReducedUnitaryValidationResult(
        targets=targets,
        corrections=corr,
        waveform=waveform,
        compiled=compiled,
        sector_metrics=tuple(sector_metrics),
        aggregate_cost=float(aggregate_cost),
        weighted_mean_process_fidelity=float(1.0 - aggregate_cost),
        weighted_mean_state_fidelity=float(np.sum(weights * state_values)),
        metadata=metadata,
    )


def analytic_kernel_matrix(model, run_config, n_active: int, *, n_pts: int = 4097) -> np.ndarray:
    levels = range(int(n_active))
    detunings = np.asarray(
        [float(manifold_transition_frequency(model, n, frame=run_config.frame)) for n in levels],
        dtype=float,
    )
    grid = np.linspace(0.0, 1.0, int(n_pts))
    envelope = np.asarray(normalized_gaussian(grid, sigma_fraction=float(run_config.sigma_fraction)), dtype=np.complex128)
    trapezoid = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    kernel = np.zeros((int(n_active), int(n_active)), dtype=np.complex128)
    for row, delta_n in enumerate(detunings):
        for col, delta_m in enumerate(detunings):
            phase = np.exp(1j * (delta_n - delta_m) * float(run_config.duration_s) * grid)
            kernel[row, col] = trapezoid(envelope * phase, grid)
    return kernel


def analytic_warm_start_corrections(
    targets: ConditionedQubitTargets,
    model,
    run_config,
    *,
    clip_bounds: tuple[float, float] | None = None,
) -> tuple[ConditionedMultitoneCorrections, dict[str, Any]]:
    n_active = int(targets.n_levels)
    kernel = analytic_kernel_matrix(model, run_config, n_active)
    lam0 = sqr_lambda0_rad_s(float(run_config.duration_s))
    target_beta = np.asarray(
        [np.sin(float(theta) / 2.0) * np.exp(1j * float(phi)) for theta, phi in zip(targets.theta, targets.phi)],
        dtype=np.complex128,
    )
    amps = np.linalg.solve(kernel, target_beta / float(run_config.duration_s))
    d_lambda: list[float] = []
    d_alpha: list[float] = []
    for idx, amp in enumerate(amps):
        base_amp = float(targets.theta[idx]) / (2.0 * float(run_config.duration_s))
        abs_amp = float(abs(amp))
        delta_lambda = 0.0 if abs(lam0) <= 1.0e-15 else float((abs_amp - base_amp) / lam0)
        if clip_bounds is not None:
            delta_lambda = float(np.clip(delta_lambda, float(clip_bounds[0]), float(clip_bounds[1])))
        d_lambda.append(delta_lambda)
        d_alpha.append(wrap_pi(float(np.angle(amp)) - float(targets.phi[idx])))
    corr = ConditionedMultitoneCorrections(
        d_lambda=tuple(float(x) for x in d_lambda),
        d_alpha=tuple(float(x) for x in d_alpha),
        d_omega_rad_s=tuple(0.0 for _ in range(n_active)),
    )
    payload = {
        "kernel_matrix": kernel,
        "kernel_condition_number": float(np.linalg.cond(kernel)),
        "analytic_complex_amplitudes_rad_s": amps,
    }
    return corr, payload


def _active_levels(targets: ConditionedQubitTargets, config: ConditionedOptimizationConfig) -> tuple[int, ...]:
    if config.active_levels:
        return tuple(int(level) for level in config.active_levels)
    return tuple(int(level) for level in range(targets.n_levels) if float(targets.weights[level]) > 0.0)


def _vector_from_corrections(
    corrections: ConditionedMultitoneCorrections,
    active_levels: Sequence[int],
    parameters: Sequence[str],
) -> np.ndarray:
    max_level = int(max(active_levels, default=-1)) + 1
    corr = corrections.padded(max(max_level, len(corrections.d_lambda), len(corrections.d_alpha), len(corrections.d_omega_rad_s)))
    vector: list[float] = []
    for level in active_levels:
        for name in parameters:
            if name == "d_lambda":
                vector.append(float(corr.d_lambda[int(level)]))
            elif name == "d_alpha":
                vector.append(float(corr.d_alpha[int(level)]))
            elif name == "d_omega":
                vector.append(float(corr.d_omega_rad_s[int(level)]))
    return np.asarray(vector, dtype=float)


def _corrections_from_vector(
    base: ConditionedMultitoneCorrections,
    vector: np.ndarray,
    n_levels: int,
    active_levels: Sequence[int],
    parameters: Sequence[str],
) -> ConditionedMultitoneCorrections:
    corr = base.padded(n_levels)
    d_lambda = np.asarray(corr.d_lambda, dtype=float)
    d_alpha = np.asarray(corr.d_alpha, dtype=float)
    d_omega = np.asarray(corr.d_omega_rad_s, dtype=float)
    data = np.asarray(vector, dtype=float).reshape(-1)
    expected = len(active_levels) * len(parameters)
    if data.size != expected:
        raise ValueError(f"Expected optimization vector of length {expected}, received {data.size}.")
    offset = 0
    for level in active_levels:
        for name in parameters:
            if name == "d_lambda":
                d_lambda[int(level)] = float(data[offset])
            elif name == "d_alpha":
                d_alpha[int(level)] = float(data[offset])
            elif name == "d_omega":
                d_omega[int(level)] = float(data[offset])
            offset += 1
    return ConditionedMultitoneCorrections(
        d_lambda=tuple(float(x) for x in d_lambda),
        d_alpha=tuple(float(x) for x in d_alpha),
        d_omega_rad_s=tuple(float(x) for x in d_omega),
    )


def _optimization_bounds(config: ConditionedOptimizationConfig, active_levels: Sequence[int]) -> Bounds:
    lower: list[float] = []
    upper: list[float] = []
    for _level in active_levels:
        for name in config.parameters:
            if name == "d_lambda":
                lower.append(float(config.d_lambda_bounds[0]))
                upper.append(float(config.d_lambda_bounds[1]))
            elif name == "d_alpha":
                lower.append(float(config.d_alpha_bounds[0]))
                upper.append(float(config.d_alpha_bounds[1]))
            elif name == "d_omega":
                lower.append(float(2.0 * np.pi * config.d_omega_hz_bounds[0]))
                upper.append(float(2.0 * np.pi * config.d_omega_hz_bounds[1]))
    return Bounds(np.asarray(lower, dtype=float), np.asarray(upper, dtype=float))


def _regularization_cost(
    corrections: ConditionedMultitoneCorrections,
    active_levels: Sequence[int],
    parameters: Sequence[str],
    config: ConditionedOptimizationConfig,
) -> float:
    value = 0.0
    for level in active_levels:
        d_lambda, d_alpha, d_omega = corrections.correction_for_n(int(level))
        for name in parameters:
            if name == "d_lambda":
                value += float(config.regularization_lambda) * float(d_lambda**2)
            elif name == "d_alpha":
                value += float(config.regularization_alpha) * float(d_alpha**2)
            elif name == "d_omega":
                value += float(config.regularization_omega) * float(d_omega**2)
    return float(value)


def optimize_reduced_unitary_multitone(
    model,
    targets: ConditionedQubitTargets,
    run_config,
    *,
    initial_corrections: ConditionedMultitoneCorrections | None = None,
    optimization_config: ConditionedOptimizationConfig | None = None,
    start_label: str = "zero",
) -> ReducedUnitaryOptimizationResult:
    opt_cfg = ConditionedOptimizationConfig() if optimization_config is None else optimization_config
    base_corr = ConditionedMultitoneCorrections.zeros(targets.n_levels) if initial_corrections is None else initial_corrections.padded(targets.n_levels)
    active_levels = _active_levels(targets, opt_cfg)
    x0 = _vector_from_corrections(base_corr, active_levels, opt_cfg.parameters)
    bounds = _optimization_bounds(opt_cfg, active_levels)
    history: list[dict[str, float]] = []

    initial_result = run_reduced_unitary_validation(model, targets, run_config, corrections=base_corr)

    def objective(vector: np.ndarray) -> float:
        corr = _corrections_from_vector(base_corr, vector, targets.n_levels, active_levels, opt_cfg.parameters)
        validation = run_reduced_unitary_validation(model, targets, run_config, corrections=corr)
        reg = _regularization_cost(corr, active_levels, opt_cfg.parameters, opt_cfg)
        objective_value = float(validation.aggregate_cost + reg)
        history.append(
            {
                "evaluation": float(len(history)),
                "aggregate_cost": float(validation.aggregate_cost),
                "regularization": float(reg),
                "objective": float(objective_value),
                "weighted_mean_process_fidelity": float(validation.weighted_mean_process_fidelity),
                "weighted_mean_state_fidelity": float(validation.weighted_mean_state_fidelity),
            }
        )
        return objective_value

    stage1 = minimize(
        objective,
        x0=x0,
        method=str(opt_cfg.method_stage1),
        bounds=bounds,
        options={"maxiter": int(opt_cfg.maxiter_stage1), "disp": False},
    )

    stage2 = None
    if opt_cfg.method_stage2:
        stage2 = minimize(
            objective,
            x0=np.asarray(stage1.x, dtype=float),
            method=str(opt_cfg.method_stage2),
            bounds=bounds,
            options={"maxiter": int(opt_cfg.maxiter_stage2)},
        )

    candidates = [_vector_from_corrections(base_corr, active_levels, opt_cfg.parameters), np.asarray(stage1.x, dtype=float)]
    scores = [float(initial_result.aggregate_cost), float(stage1.fun)]
    if stage2 is not None:
        candidates.append(np.asarray(stage2.x, dtype=float))
        scores.append(float(stage2.fun))
    best = candidates[int(np.argmin(scores))]
    optimized_corrections = _corrections_from_vector(base_corr, best, targets.n_levels, active_levels, opt_cfg.parameters)
    optimized_result = run_reduced_unitary_validation(model, targets, run_config, corrections=optimized_corrections)
    return ReducedUnitaryOptimizationResult(
        initial_result=initial_result,
        optimized_result=optimized_result,
        optimized_corrections=optimized_corrections,
        active_levels=tuple(int(level) for level in active_levels),
        parameters=tuple(str(name) for name in opt_cfg.parameters),
        history=tuple(history),
        success_stage1=bool(stage1.success),
        success_stage2=False if stage2 is None else bool(stage2.success),
        message_stage1=str(stage1.message),
        message_stage2="" if stage2 is None else str(stage2.message),
        start_label=str(start_label),
    )
