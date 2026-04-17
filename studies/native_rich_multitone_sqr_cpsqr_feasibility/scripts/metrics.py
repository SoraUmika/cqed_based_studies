from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np
from scipy.linalg import norm, polar
from scipy.optimize import minimize_scalar


TWO_PI = 2.0 * np.pi

PAULI_X = np.asarray([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
PAULI_Y = np.asarray([[0.0, -1.0j], [1.0j, 0.0]], dtype=np.complex128)
PAULI_Z = np.asarray([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)

PROBE_QUBIT_STATES: dict[str, np.ndarray] = {
    "g": np.asarray([1.0, 0.0], dtype=np.complex128),
    "e": np.asarray([0.0, 1.0], dtype=np.complex128),
    "plus_x": np.asarray([1.0, 1.0], dtype=np.complex128) / np.sqrt(2.0),
    "plus_y": np.asarray([1.0, 1.0j], dtype=np.complex128) / np.sqrt(2.0),
}

PROBE_TIERS: dict[str, tuple[str, ...]] = {
    "single_ground": ("g",),
    "selected_pair": ("g", "plus_x"),
    "spanning_quartet": ("g", "e", "plus_x", "plus_y"),
}


@dataclass(frozen=True)
class CPSQRFit:
    delta_rad: float
    block_phase_rad: float
    process_fidelity: float
    target_block: np.ndarray


def wrap_pi(value: float) -> float:
    return float((float(value) + np.pi) % (2.0 * np.pi) - np.pi)


def qubit_rx(theta: float) -> np.ndarray:
    half = 0.5 * float(theta)
    return np.asarray(
        [[np.cos(half), -1.0j * np.sin(half)], [-1.0j * np.sin(half), np.cos(half)]],
        dtype=np.complex128,
    )


def qubit_rz(theta: float) -> np.ndarray:
    half = 0.5 * float(theta)
    return np.asarray(
        [[np.exp(-1.0j * half), 0.0], [0.0, np.exp(1.0j * half)]],
        dtype=np.complex128,
    )


def cpsqr_block(theta: float, delta_rad: float, block_phase_rad: float = 0.0) -> np.ndarray:
    return np.exp(1.0j * float(block_phase_rad)) * qubit_rz(delta_rad) @ qubit_rx(theta)


def process_fidelity(target_operator: np.ndarray, actual_operator: np.ndarray) -> float:
    target = np.asarray(target_operator, dtype=np.complex128)
    actual = np.asarray(actual_operator, dtype=np.complex128)
    dim = float(target.shape[0])
    overlap = np.trace(target.conj().T @ actual)
    return float(np.clip(abs(overlap) ** 2 / (dim * dim), 0.0, 1.0))


def average_gate_fidelity(target_operator: np.ndarray, actual_operator: np.ndarray) -> float:
    dim = float(np.asarray(target_operator).shape[0])
    proc = process_fidelity(target_operator, actual_operator)
    return float((dim * proc + 1.0) / (dim + 1.0))


def frobenius_error(target_operator: np.ndarray, actual_operator: np.ndarray) -> float:
    target = np.asarray(target_operator, dtype=np.complex128)
    actual = np.asarray(actual_operator, dtype=np.complex128)
    denom = float(norm(target, ord="fro"))
    raw = float(norm(actual - target, ord="fro"))
    return raw if denom <= 0.0 else raw / denom


def operator_2norm_error(target_operator: np.ndarray, actual_operator: np.ndarray) -> float:
    target = np.asarray(target_operator, dtype=np.complex128)
    actual = np.asarray(actual_operator, dtype=np.complex128)
    return float(norm(actual - target, ord=2))


def nearest_unitary(matrix: np.ndarray) -> np.ndarray:
    unitary, _ = polar(np.asarray(matrix, dtype=np.complex128))
    det = np.linalg.det(unitary)
    if abs(det) > 1.0e-15:
        unitary = unitary * np.exp(-0.5j * np.angle(det))
    return np.asarray(unitary, dtype=np.complex128)


def unitary_rotation_parameters(unitary: np.ndarray) -> tuple[float, float, float]:
    u = nearest_unitary(unitary)
    trace = np.trace(u)
    cos_half = float(np.clip(np.real(trace / 2.0), -1.0, 1.0))
    theta = float(2.0 * np.arccos(cos_half))
    if theta < 1.0e-12:
        return 0.0, 0.0, 0.0
    sin_half = float(np.sin(theta / 2.0))
    nx = float(np.real(1.0j * np.trace(PAULI_X @ u) / (2.0 * sin_half)))
    ny = float(np.real(1.0j * np.trace(PAULI_Y @ u) / (2.0 * sin_half)))
    nz = float(np.real(1.0j * np.trace(PAULI_Z @ u) / (2.0 * sin_half)))
    axis_norm = float(np.sqrt(nx * nx + ny * ny + nz * nz))
    if axis_norm > 1.0e-12:
        nx /= axis_norm
        ny /= axis_norm
        nz /= axis_norm
    phi = float(np.mod(np.arctan2(ny, nx), TWO_PI))
    return theta, phi, nz


def bloch_vector_from_density_matrix(rho: np.ndarray) -> tuple[float, float, float]:
    matrix = np.asarray(rho, dtype=np.complex128)
    return (
        float(np.real(np.trace(matrix @ PAULI_X))),
        float(np.real(np.trace(matrix @ PAULI_Y))),
        float(np.real(np.trace(matrix @ PAULI_Z))),
    )


def state_density_matrix(qubit_state: Sequence[complex]) -> np.ndarray:
    vec = np.asarray(qubit_state, dtype=np.complex128).reshape((2, 1))
    return np.asarray(vec @ vec.conj().T, dtype=np.complex128)


def state_fidelity_from_dm(actual_rho: np.ndarray, target_rho: np.ndarray) -> float:
    actual = np.asarray(actual_rho, dtype=np.complex128)
    target = np.asarray(target_rho, dtype=np.complex128)
    return float(np.clip(np.real(np.trace(actual @ target)), 0.0, 1.0))


def full_state_fidelity(actual_state: np.ndarray, target_state: np.ndarray) -> float:
    actual = np.asarray(actual_state, dtype=np.complex128).reshape(-1)
    target = np.asarray(target_state, dtype=np.complex128).reshape(-1)
    return float(np.clip(abs(np.vdot(target, actual)) ** 2, 0.0, 1.0))


def qubit_channel_kraus_from_full(full_operator: np.ndarray, n_cav: int, input_level: int) -> tuple[np.ndarray, ...]:
    full = np.asarray(full_operator, dtype=np.complex128)
    kraus: list[np.ndarray] = []
    for output_level in range(int(n_cav)):
        block = np.zeros((2, 2), dtype=np.complex128)
        for q_out in range(2):
            for q_in in range(2):
                row = q_out * int(n_cav) + int(output_level)
                col = q_in * int(n_cav) + int(input_level)
                block[q_out, q_in] = full[row, col]
        kraus.append(block)
    return tuple(kraus)


def same_manifold_block(full_operator: np.ndarray, n_cav: int, level: int) -> np.ndarray:
    return np.asarray(qubit_channel_kraus_from_full(full_operator, n_cav, level)[int(level)], dtype=np.complex128)


def apply_channel(kraus_ops: Sequence[np.ndarray], rho: np.ndarray) -> np.ndarray:
    state = np.asarray(rho, dtype=np.complex128)
    out = np.zeros((2, 2), dtype=np.complex128)
    for op in kraus_ops:
        mat = np.asarray(op, dtype=np.complex128)
        out = out + mat @ state @ mat.conj().T
    return out


def channel_process_fidelity_to_unitary(kraus_ops: Sequence[np.ndarray], target_unitary: np.ndarray) -> float:
    target = np.asarray(target_unitary, dtype=np.complex128)
    accum = 0.0
    for op in kraus_ops:
        accum += abs(np.trace(target.conj().T @ np.asarray(op, dtype=np.complex128))) ** 2
    return float(np.clip(accum / 4.0, 0.0, 1.0))


def fit_cpsqr_channel(kraus_ops: Sequence[np.ndarray], theta_target: float) -> CPSQRFit:
    def objective(delta_rad: float) -> float:
        return -channel_process_fidelity_to_unitary(kraus_ops, qubit_rz(delta_rad) @ qubit_rx(theta_target))

    result = minimize_scalar(objective, bounds=(-np.pi, np.pi), method="bounded", options={"xatol": 1.0e-4})
    delta = float(result.x)
    target = qubit_rz(delta) @ qubit_rx(theta_target)
    return CPSQRFit(
        delta_rad=delta,
        block_phase_rad=0.0,
        process_fidelity=channel_process_fidelity_to_unitary(kraus_ops, target),
        target_block=target,
    )


def fit_cpsqr_block(actual_block: np.ndarray, theta_target: float) -> CPSQRFit:
    actual = np.asarray(actual_block, dtype=np.complex128)

    def overlap_abs(delta_rad: float) -> float:
        target = qubit_rz(delta_rad) @ qubit_rx(theta_target)
        return float(abs(np.trace(target.conj().T @ actual)))

    result = minimize_scalar(lambda delta: -overlap_abs(delta), bounds=(-np.pi, np.pi), method="bounded", options={"xatol": 1.0e-4})
    delta = float(result.x)
    base_target = qubit_rz(delta) @ qubit_rx(theta_target)
    term = np.trace(base_target.conj().T @ actual)
    block_phase = float(np.angle(term))
    target = np.exp(1.0j * block_phase) * base_target
    return CPSQRFit(
        delta_rad=delta,
        block_phase_rad=block_phase,
        process_fidelity=process_fidelity(target, actual),
        target_block=target,
    )


def build_cpsqr_joint_target(restricted_operator: np.ndarray, theta_values: Sequence[float]) -> tuple[np.ndarray, list[dict[str, float]]]:
    restricted = np.asarray(restricted_operator, dtype=np.complex128)
    n_blocks = restricted.shape[0] // 2
    target = np.zeros_like(restricted)
    rows: list[dict[str, float]] = []
    for idx in range(n_blocks):
        block = restricted[2 * idx : 2 * idx + 2, 2 * idx : 2 * idx + 2]
        fit = fit_cpsqr_block(block, float(theta_values[idx]))
        target[2 * idx : 2 * idx + 2, 2 * idx : 2 * idx + 2] = fit.target_block
        rows.append(
            {
                "level": int(idx),
                "delta_rad": float(fit.delta_rad),
                "block_phase_rad": float(fit.block_phase_rad),
                "process_fidelity": float(fit.process_fidelity),
            }
        )
    return target, rows


def qubit_probe_fidelity_rows_for_channel(
    kraus_ops: Sequence[np.ndarray],
    *,
    strict_unitary: np.ndarray,
    cpsqr_unitary: np.ndarray,
) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    for name, state in PROBE_QUBIT_STATES.items():
        rho_in = state_density_matrix(state)
        actual = apply_channel(kraus_ops, rho_in)
        strict_target = strict_unitary @ np.asarray(state, dtype=np.complex128).reshape((2, 1))
        cpsqr_target = cpsqr_unitary @ np.asarray(state, dtype=np.complex128).reshape((2, 1))
        strict_rho = strict_target @ strict_target.conj().T
        cpsqr_rho = cpsqr_target @ cpsqr_target.conj().T
        rows.append(
            {
                "probe_label": str(name),
                "strict_fidelity": state_fidelity_from_dm(actual, strict_rho),
                "cpsqr_fidelity": state_fidelity_from_dm(actual, cpsqr_rho),
                "bloch_x": bloch_vector_from_density_matrix(actual)[0],
                "bloch_y": bloch_vector_from_density_matrix(actual)[1],
                "bloch_z": bloch_vector_from_density_matrix(actual)[2],
            }
        )
    return rows


def probe_tier_summary(probe_rows: Sequence[dict[str, Any]], metric_key: str) -> dict[str, dict[str, float]]:
    summary: dict[str, dict[str, float]] = {}
    for tier_name, labels in PROBE_TIERS.items():
        values = [float(row[metric_key]) for row in probe_rows if str(row["probe_label"]) in labels]
        if values:
            summary[tier_name] = {
                "mean_fidelity": float(np.mean(values)),
                "min_fidelity": float(np.min(values)),
            }
        else:
            summary[tier_name] = {"mean_fidelity": float("nan"), "min_fidelity": float("nan")}
    return summary


def coherent_error_decomposition(actual_block: np.ndarray, target_block: np.ndarray) -> dict[str, float]:
    actual = nearest_unitary(actual_block)
    target = nearest_unitary(target_block)
    error = nearest_unitary(target.conj().T @ actual)
    trace = np.trace(error)
    cos_half = float(np.clip(np.real(trace / 2.0), -1.0, 1.0))
    theta = float(2.0 * np.arccos(cos_half))
    if theta < 1.0e-12:
        rotvec = np.zeros(3, dtype=float)
    else:
        sin_half = float(np.sin(theta / 2.0))
        nx = float(np.real(1.0j * np.trace(PAULI_X @ error) / (2.0 * sin_half)))
        ny = float(np.real(1.0j * np.trace(PAULI_Y @ error) / (2.0 * sin_half)))
        nz = float(np.real(1.0j * np.trace(PAULI_Z @ error) / (2.0 * sin_half)))
        axis = np.asarray([nx, ny, nz], dtype=float)
        norm_axis = float(np.linalg.norm(axis))
        if norm_axis > 1.0e-12:
            axis = axis / norm_axis
        rotvec = axis * theta
    return {
        "error_rotation_angle_rad": float(np.linalg.norm(rotvec)),
        "residual_z_error_rad": float(abs(rotvec[2])),
        "transverse_error_rad": float(np.linalg.norm(rotvec[:2])),
        "error_rotvec_x_rad": float(rotvec[0]),
        "error_rotvec_y_rad": float(rotvec[1]),
        "error_rotvec_z_rad": float(rotvec[2]),
    }


def addressed_indices(n_cav: int, levels: Sequence[int]) -> np.ndarray:
    rows: list[int] = []
    for level in levels:
        rows.extend([int(level), int(n_cav) + int(level)])
    return np.asarray(rows, dtype=int)


def leakage_outside_indices(state_vector: np.ndarray, keep_indices: Sequence[int]) -> float:
    vec = np.asarray(state_vector, dtype=np.complex128).reshape(-1)
    keep = np.asarray(tuple(int(index) for index in keep_indices), dtype=int)
    keep_prob = float(np.sum(np.abs(vec[keep]) ** 2))
    return float(max(0.0, 1.0 - keep_prob))
