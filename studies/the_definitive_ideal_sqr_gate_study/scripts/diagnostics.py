from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.linalg import polar


PAULI_X = np.asarray([[0.0, 1.0], [1.0, 0.0]], dtype=np.complex128)
PAULI_Y = np.asarray([[0.0, -1.0j], [1.0j, 0.0]], dtype=np.complex128)
PAULI_Z = np.asarray([[1.0, 0.0], [0.0, -1.0]], dtype=np.complex128)


@dataclass(frozen=True)
class DiagnosticTask:
    study: str
    case_id: str
    construction: str
    artifact_path: str
    target_family: str
    n_active: int
    random_seed: int | None
    target_kind: str
    strict_process_fidelity: float | None


def complex_matrix_from_json(payload: Any) -> np.ndarray:
    if isinstance(payload, dict) and "real" in payload and "imag" in payload:
        return np.asarray(payload["real"], dtype=float) + 1.0j * np.asarray(payload["imag"], dtype=float)
    return np.asarray(payload, dtype=np.complex128)


def nearest_unitary(matrix: np.ndarray) -> np.ndarray:
    unitary, _ = polar(np.asarray(matrix, dtype=np.complex128))
    det = np.linalg.det(unitary)
    if abs(det) > 1.0e-15:
        unitary = unitary * np.exp(-0.5j * np.angle(det))
    return np.asarray(unitary, dtype=np.complex128)


def rx(theta_rad: float) -> np.ndarray:
    half = 0.5 * float(theta_rad)
    return np.asarray(
        [[np.cos(half), -1.0j * np.sin(half)], [-1.0j * np.sin(half), np.cos(half)]],
        dtype=np.complex128,
    )


def ideal_theta_values(target_family: str, n_active: int, random_seed: int | None = None) -> tuple[float, ...]:
    n = int(n_active)
    if target_family == "smooth_x":
        return tuple(float(min(0.30 * np.pi + 0.18 * np.pi * level, 0.92 * np.pi)) for level in range(n))
    if target_family == "staggered_x":
        rows = [0.22 * np.pi + 0.58 * np.pi * (level % 2) + 0.06 * np.pi * (level // 2) for level in range(n)]
        return tuple(float(min(value, 0.94 * np.pi)) for value in rows)
    if target_family == "random_x":
        rng = np.random.default_rng(None if random_seed is None else int(random_seed))
        return tuple(float(rng.uniform(0.18 * np.pi, 0.92 * np.pi)) for _ in range(n))
    raise ValueError(f"Unsupported ideal-SQR family '{target_family}'.")


def rotation_vector_from_error_unitary(error_unitary: np.ndarray) -> np.ndarray:
    error = nearest_unitary(error_unitary)
    trace = np.trace(error)
    cos_half = float(np.clip(np.real(trace / 2.0), -1.0, 1.0))
    theta = float(2.0 * np.arccos(cos_half))
    if theta < 1.0e-12:
        return np.zeros(3, dtype=float)
    sin_half = float(np.sin(theta / 2.0))
    if abs(sin_half) < 1.0e-12:
        return np.zeros(3, dtype=float)
    nx = float(np.real(1.0j * np.trace(PAULI_X @ error) / (2.0 * sin_half)))
    ny = float(np.real(1.0j * np.trace(PAULI_Y @ error) / (2.0 * sin_half)))
    nz = float(np.real(1.0j * np.trace(PAULI_Z @ error) / (2.0 * sin_half)))
    axis = np.asarray([nx, ny, nz], dtype=float)
    axis_norm = float(np.linalg.norm(axis))
    if axis_norm > 1.0e-12:
        axis = axis / axis_norm
    return axis * theta


def same_manifold_block(full_operator: np.ndarray, n_cav: int, level: int) -> np.ndarray:
    rows = [int(level), int(n_cav) + int(level)]
    cols = [int(level), int(n_cav) + int(level)]
    return np.asarray(full_operator[np.ix_(rows, cols)], dtype=np.complex128)


def unitary_rotation_angle(unitary: np.ndarray) -> float:
    matrix = nearest_unitary(unitary)
    trace = np.trace(matrix)
    cos_half = float(np.clip(np.real(trace / 2.0), -1.0, 1.0))
    return float(2.0 * np.arccos(cos_half))


def _rows_from_native_rich_artifact(task: DiagnosticTask, artifact: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for item in artifact.get("reduced_level_rows", []):
        eps_x = float(item["strict_error_rotvec_x_rad"])
        eps_y = float(item["strict_error_rotvec_y_rad"])
        eps_z = float(item["strict_error_rotvec_z_rad"])
        rows.append(
            {
                "study": task.study,
                "case_id": task.case_id,
                "construction": task.construction,
                "artifact_path": task.artifact_path,
                "target_family": task.target_family,
                "n_active": int(task.n_active),
                "n": int(item["level"]),
                "eps_x": eps_x,
                "eps_y": eps_y,
                "eps_z": eps_z,
                "eps_norm": float(np.linalg.norm([eps_x, eps_y, eps_z])),
                "target_theta_rad": float(item["target_theta_rad"]),
                "actual_theta_rad": float(item["achieved_theta_rad"]),
                "strict_process_fidelity": task.strict_process_fidelity,
            }
        )
    return rows


def _rows_from_legacy_artifact(task: DiagnosticTask, artifact: dict[str, Any]) -> list[dict[str, Any]]:
    full_operator = complex_matrix_from_json(artifact["full_operator_columns_on_logical_inputs"])
    n_cav = int(full_operator.shape[0] // 2)
    thetas = ideal_theta_values(task.target_family, task.n_active, task.random_seed)
    rows: list[dict[str, Any]] = []
    for level, theta in enumerate(thetas):
        actual_block = nearest_unitary(same_manifold_block(full_operator, n_cav, level))
        ideal_block = rx(theta)
        error_operator = actual_block @ ideal_block.conj().T
        rotvec = rotation_vector_from_error_unitary(error_operator)
        rows.append(
            {
                "study": task.study,
                "case_id": task.case_id,
                "construction": task.construction,
                "artifact_path": task.artifact_path,
                "target_family": task.target_family,
                "n_active": int(task.n_active),
                "n": int(level),
                "eps_x": float(rotvec[0]),
                "eps_y": float(rotvec[1]),
                "eps_z": float(rotvec[2]),
                "eps_norm": float(np.linalg.norm(rotvec)),
                "target_theta_rad": float(theta),
                "actual_theta_rad": unitary_rotation_angle(actual_block),
                "strict_process_fidelity": task.strict_process_fidelity,
            }
        )
    return rows


def run_task(task: DiagnosticTask) -> list[dict[str, Any]]:
    artifact = json.loads(Path(task.artifact_path).read_text(encoding="utf-8"))
    if "reduced_level_rows" in artifact:
        return _rows_from_native_rich_artifact(task, artifact)
    return _rows_from_legacy_artifact(task, artifact)


def build_tasks(records: list[dict[str, Any]]) -> list[DiagnosticTask]:
    tasks: list[DiagnosticTask] = []
    for row in records:
        if row.get("target_kind") != "ideal_sqr":
            continue
        artifact_path = row.get("artifact_path")
        if not artifact_path:
            continue
        path = Path(str(artifact_path))
        if not path.exists():
            continue
        random_seed = row.get("random_seed")
        if random_seed in (None, ""):
            parsed_seed = None
        else:
            try:
                parsed_seed = int(random_seed)
            except (TypeError, ValueError):
                parsed_seed = None
        tasks.append(
            DiagnosticTask(
                study=str(row["study"]),
                case_id=str(row["case_id"]),
                construction=str(row["construction"]),
                artifact_path=str(path),
                target_family=str(row["target_family"]),
                n_active=int(row["n_active"]),
                random_seed=parsed_seed,
                target_kind=str(row["target_kind"]),
                strict_process_fidelity=None if row.get("strict_process_fidelity") in (None, "") else float(row["strict_process_fidelity"]),
            )
        )
    return tasks
