"""
Reusable phase-compilation diagnostics for the SQR follow-up study.

This module keeps the simulation layer inside ``cqed_sim`` and provides a
consistent post-processing stack for:

- restricted logical-operator extraction,
- cavity-only block-phase fitting and polynomial compression,
- branch-local virtual-Z diagnostics,
- same-block / leakage accounting,
- coherent superposition transfer checks.
"""

from __future__ import annotations

from itertools import combinations

import numpy as np
from scipy.linalg import block_diag

from common import (
    extract_branch_unitaries,
    extract_leakage,
    target_qubit_unitary,
    z_corrected_target_fidelity,
)
from cqed_sim.core.ideal_gates import logical_block_phase_op
from cqed_sim.unitary_synthesis.metrics import logical_block_phase_diagnostics


def logical_basis_indices(model, logical_n: int) -> tuple[int, ...]:
    """Indices of the logical |g,n>, |e,n> basis in the full Hilbert space."""
    indices: list[int] = []
    for n in range(int(logical_n)):
        indices.extend([n, int(model.n_cav) + n])
    return tuple(indices)


def block_slices(logical_n: int) -> tuple[slice, ...]:
    """Return 2x2 logical block slices ordered by Fock level."""
    return tuple(slice(2 * n, 2 * n + 2) for n in range(int(logical_n)))


def restricted_operator(full_operator: np.ndarray, model, logical_n: int) -> np.ndarray:
    """Restrict a full operator to the logical |g/e> x Fock subspace."""
    indices = np.asarray(logical_basis_indices(model, logical_n), dtype=int)
    return np.asarray(full_operator[np.ix_(indices, indices)], dtype=np.complex128)


def full_operator_from_basis_outputs(final_states, model, logical_n: int) -> np.ndarray:
    """Reconstruct the logical columns of the full propagator from basis outputs."""
    full_dim = int(model.n_tr) * int(model.n_cav)
    operator = np.zeros((full_dim, full_dim), dtype=np.complex128)
    for n in range(int(logical_n)):
        for qubit_level in (0, 1):
            column_index = qubit_level * int(model.n_cav) + n
            operator[:, column_index] = np.asarray(
                final_states[2 * n + qubit_level].full(),
                dtype=np.complex128,
            ).reshape(-1)
    return operator


def single_target_blocks(logical_n: int, target_branch: int, theta: float, phi: float) -> list[np.ndarray]:
    """Target blocks for an SQR acting only on ``target_branch``."""
    target_block = target_qubit_unitary(theta, phi)
    identity = np.eye(2, dtype=np.complex128)
    return [
        target_block if n == int(target_branch) else identity
        for n in range(int(logical_n))
    ]


def allbranch_blocks(logical_n: int, theta: float, phi: float) -> list[np.ndarray]:
    """Target blocks for simultaneous rotation on every logical Fock branch."""
    target_block = target_qubit_unitary(theta, phi)
    return [target_block.copy() for _ in range(int(logical_n))]


def restricted_target_from_blocks(target_blocks: list[np.ndarray]) -> np.ndarray:
    """Block-diagonal restricted target operator."""
    return block_diag(*target_blocks)


def full_target_from_blocks(model, logical_n: int, target_blocks: list[np.ndarray]) -> np.ndarray:
    """Embed 2x2 logical target blocks into the full qutrit-cavity space."""
    full_dim = int(model.n_tr) * int(model.n_cav)
    target = np.eye(full_dim, dtype=np.complex128)
    for n, block in enumerate(target_blocks[: int(logical_n)]):
        rows = [n, int(model.n_cav) + n]
        target[np.ix_(rows, rows)] = np.asarray(block, dtype=np.complex128)
    return target


def process_fidelity(actual: np.ndarray, target: np.ndarray) -> float:
    """Standard process fidelity |Tr(V^dagger U)|^2 / d^2."""
    d = float(target.shape[0])
    return float(np.clip(abs(np.trace(target.conj().T @ actual)) ** 2 / (d * d), 0.0, 1.0))


def best_global_qubit_z_fidelity(
    actual: np.ndarray,
    target: np.ndarray,
    logical_n: int,
    *,
    n_grid: int = 1441,
) -> tuple[float, float]:
    """Maximize restricted process fidelity over one global virtual-Z on the qubit."""
    best_fid = -1.0
    best_alpha = 0.0
    phase_grid = np.linspace(0.0, 2.0 * np.pi, int(n_grid))
    for alpha in phase_grid:
        z_layer = np.diag(np.tile([1.0, np.exp(1j * alpha)], int(logical_n)))
        fid = process_fidelity(z_layer @ actual, target)
        if fid > best_fid:
            best_fid = float(fid)
            best_alpha = float(alpha)
    return best_fid, best_alpha


def global_qubit_z_full_matrix(model, alpha: float) -> np.ndarray:
    """Full-space global qubit Z on the g/e manifold; |f> is left unchanged."""
    diag = np.ones(int(model.n_tr), dtype=np.complex128)
    if int(model.n_tr) >= 2:
        diag[1] = np.exp(1j * float(alpha))
    return np.kron(np.diag(diag), np.eye(int(model.n_cav), dtype=np.complex128))


def gauge_fixed_phase_profile(phases_rad: np.ndarray | list[float]) -> np.ndarray:
    """Continuously unwrap and reference the phase profile to level n=0."""
    profile = np.unwrap(np.asarray(phases_rad, dtype=float).reshape(-1))
    return profile - float(profile[0]) if profile.size else profile


def polynomial_phase_profile(phases_rad: np.ndarray | list[float], order: int) -> tuple[np.ndarray, np.ndarray, float]:
    """Fit a low-order polynomial to a gauge-fixed phase profile."""
    profile = gauge_fixed_phase_profile(phases_rad)
    n = np.arange(profile.size, dtype=float)
    if profile.size == 0:
        coeffs = np.zeros(int(order) + 1, dtype=float)
        return profile.copy(), coeffs, 0.0
    coeffs = np.polyfit(n, profile, int(order))
    fitted = np.polyval(coeffs, n)
    fitted = fitted - float(fitted[0])
    residual = float(np.sqrt(np.mean((fitted - profile) ** 2)))
    return fitted, coeffs, residual


def apply_cavity_phase_compilation(
    full_operator: np.ndarray,
    model,
    phases_rad: np.ndarray | list[float],
) -> np.ndarray:
    """Apply a cavity-only phase layer to the full qutrit-cavity operator."""
    layer = logical_block_phase_op(
        np.asarray(phases_rad, dtype=float),
        cavity_dim=int(model.n_cav),
        qubit_dim=int(model.n_tr),
    )
    return np.asarray(layer.full(), dtype=np.complex128) @ np.asarray(full_operator, dtype=np.complex128)


def local_z_relaxed_fidelity(
    branch_unitaries: list[np.ndarray],
    target_blocks: list[np.ndarray],
) -> tuple[float, np.ndarray, np.ndarray]:
    """Per-block virtual-Z-relaxed fidelity and fitted phases."""
    branch_fids = np.zeros(len(branch_unitaries), dtype=float)
    branch_phases = np.zeros(len(branch_unitaries), dtype=float)
    accum = 0.0
    for index, (actual_block, target_block) in enumerate(zip(branch_unitaries, target_blocks, strict=True)):
        fid, alpha = z_corrected_target_fidelity(actual_block, target_block)
        branch_fids[index] = float(fid)
        branch_phases[index] = float(alpha)
        accum += np.sqrt(float(fid))
    relaxed = (accum / max(len(branch_unitaries), 1)) ** 2
    return float(np.clip(relaxed, 0.0, 1.0)), branch_fids, branch_phases


def same_block_population_metrics(full_operator: np.ndarray, model, logical_n: int) -> tuple[float, float]:
    """Average and minimum population retained inside the intended 2x2 logical block."""
    values: list[float] = []
    for n in range(int(logical_n)):
        for qubit_level in (0, 1):
            logical_column = qubit_level * int(model.n_cav) + n
            column = np.asarray(full_operator[:, logical_column], dtype=np.complex128).reshape(-1)
            same_block = abs(column[n]) ** 2 + abs(column[int(model.n_cav) + n]) ** 2
            values.append(float(same_block))
    arr = np.asarray(values, dtype=float)
    return float(np.mean(arr)), float(np.min(arr))


def pair_superposition_fidelity_stats(
    corrected_full_operator: np.ndarray,
    target_full_operator: np.ndarray,
    model,
    logical_n: int,
) -> dict[str, float]:
    """Mean/min fidelity over explicit inter-Fock coherent superposition probes."""
    fidelities: list[float] = []
    full_dim = int(model.n_tr) * int(model.n_cav)
    for left, right in combinations(range(int(logical_n)), 2):
        for qubit_level in (0, 1):
            psi = np.zeros(full_dim, dtype=np.complex128)
            psi[qubit_level * int(model.n_cav) + left] = 1.0 / np.sqrt(2.0)
            psi[qubit_level * int(model.n_cav) + right] = 1.0 / np.sqrt(2.0)
            out = corrected_full_operator @ psi
            target = target_full_operator @ psi
            fidelities.append(float(np.clip(abs(np.vdot(target, out)) ** 2, 0.0, 1.0)))

        psi_cross = np.zeros(full_dim, dtype=np.complex128)
        psi_cross[0 * int(model.n_cav) + left] = 1.0 / np.sqrt(2.0)
        psi_cross[1 * int(model.n_cav) + right] = 1.0 / np.sqrt(2.0)
        out_cross = corrected_full_operator @ psi_cross
        target_cross = target_full_operator @ psi_cross
        fidelities.append(float(np.clip(abs(np.vdot(target_cross, out_cross)) ** 2, 0.0, 1.0)))

    arr = np.asarray(fidelities, dtype=float)
    return {
        "mean": float(np.mean(arr)) if arr.size else 0.0,
        "min": float(np.min(arr)) if arr.size else 0.0,
        "count": int(arr.size),
    }


def compiled_phase_metrics(
    full_operator: np.ndarray,
    final_states,
    model,
    logical_n: int,
    target_blocks: list[np.ndarray],
    *,
    coherence_stats: bool = True,
) -> dict[str, float | np.ndarray]:
    """Compute strict, compiled, and branch-local phase diagnostics for one case."""
    logical_n = int(logical_n)
    target_restricted = restricted_target_from_blocks(target_blocks)
    target_full = full_target_from_blocks(model, logical_n, target_blocks)
    restricted = restricted_operator(full_operator, model, logical_n)
    branches = extract_branch_unitaries(final_states, model, logical_n)
    leakage = extract_leakage(final_states, model, logical_n)

    raw_strict_fid = process_fidelity(restricted, target_restricted)
    raw_global_z_fid, raw_global_z_alpha = best_global_qubit_z_fidelity(
        restricted, target_restricted, logical_n
    )

    relaxed_fid, branch_local_z_fids, branch_local_z_phases = local_z_relaxed_fidelity(
        branches, target_blocks
    )

    diagnostics = logical_block_phase_diagnostics(
        restricted,
        target_restricted,
        block_slices=block_slices(logical_n),
    )
    exact_cavity_phases = gauge_fixed_phase_profile(diagnostics.best_fit_correction_phases_rad)
    exact_compiled_full = apply_cavity_phase_compilation(full_operator, model, exact_cavity_phases)
    exact_compiled_restricted = restricted_operator(exact_compiled_full, model, logical_n)
    exact_compiled_fid, exact_compiled_alpha = best_global_qubit_z_fidelity(
        exact_compiled_restricted, target_restricted, logical_n
    )

    fit_profiles: dict[str, np.ndarray] = {}
    fit_coeffs: dict[str, np.ndarray] = {}
    fit_rms: dict[str, float] = {}
    fit_fids: dict[str, float] = {}
    for label, order in (("linear", 1), ("quadratic", 2), ("cubic", 3)):
        fitted_profile, coeffs, residual = polynomial_phase_profile(exact_cavity_phases, order)
        fit_profiles[label] = fitted_profile
        fit_coeffs[label] = coeffs
        fit_rms[label] = residual
        compiled_full = apply_cavity_phase_compilation(full_operator, model, fitted_profile)
        compiled_restricted = restricted_operator(compiled_full, model, logical_n)
        fit_fids[label], _ = best_global_qubit_z_fidelity(compiled_restricted, target_restricted, logical_n)

    same_block_mean, same_block_min = same_block_population_metrics(full_operator, model, logical_n)

    raw_corrected_full = global_qubit_z_full_matrix(model, raw_global_z_alpha) @ np.asarray(
        full_operator, dtype=np.complex128
    )
    exact_corrected_full = global_qubit_z_full_matrix(model, exact_compiled_alpha) @ exact_compiled_full

    metrics: dict[str, float | np.ndarray] = {
        "raw_strict_fid": float(raw_strict_fid),
        "raw_global_z_fid": float(raw_global_z_fid),
        "raw_global_z_alpha": float(raw_global_z_alpha),
        "branch_local_z_relaxed_fid": float(relaxed_fid),
        "branch_local_z_fids": branch_local_z_fids,
        "branch_local_z_phases": branch_local_z_phases,
        "exact_cavity_compiled_fid": float(exact_compiled_fid),
        "exact_cavity_global_z_alpha": float(exact_compiled_alpha),
        "exact_cavity_phases": exact_cavity_phases,
        "linear_cavity_compiled_fid": float(fit_fids["linear"]),
        "quadratic_cavity_compiled_fid": float(fit_fids["quadratic"]),
        "cubic_cavity_compiled_fid": float(fit_fids["cubic"]),
        "linear_phase_fit_rms": float(fit_rms["linear"]),
        "quadratic_phase_fit_rms": float(fit_rms["quadratic"]),
        "cubic_phase_fit_rms": float(fit_rms["cubic"]),
        "linear_phase_profile": fit_profiles["linear"],
        "quadratic_phase_profile": fit_profiles["quadratic"],
        "cubic_phase_profile": fit_profiles["cubic"],
        "linear_phase_coeffs": fit_coeffs["linear"],
        "quadratic_phase_coeffs": fit_coeffs["quadratic"],
        "cubic_phase_coeffs": fit_coeffs["cubic"],
        "same_block_population_mean": float(same_block_mean),
        "same_block_population_min": float(same_block_min),
        "leakage_mean": float(np.mean(leakage)),
        "leakage_max": float(np.max(leakage)),
    }

    if coherence_stats:
        raw_pairs = pair_superposition_fidelity_stats(raw_corrected_full, target_full, model, logical_n)
        exact_pairs = pair_superposition_fidelity_stats(exact_corrected_full, target_full, model, logical_n)
        metrics.update(
            {
                "pair_superposition_count": float(raw_pairs["count"]),
                "pair_superposition_raw_mean": float(raw_pairs["mean"]),
                "pair_superposition_raw_min": float(raw_pairs["min"]),
                "pair_superposition_compiled_mean": float(exact_pairs["mean"]),
                "pair_superposition_compiled_min": float(exact_pairs["min"]),
            }
        )

    return metrics
