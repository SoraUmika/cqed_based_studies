"""Linear-dispersive bounds and matched-filter calculations."""

from __future__ import annotations

from dataclasses import dataclass
from math import erf, sqrt

import numpy as np


@dataclass(frozen=True)
class LinearResponse:
    """State-conditioned linear resonator response under an arbitrary waveform."""

    tlist: np.ndarray
    epsilon: np.ndarray
    alpha_g: np.ndarray
    alpha_e: np.ndarray
    output_g: np.ndarray
    output_e: np.ndarray
    delta_g: float
    delta_e: float
    kappa: float


def solve_linear_response(
    epsilon: np.ndarray,
    tlist: np.ndarray,
    *,
    kappa: float,
    chi: float,
    delta_g: float,
) -> LinearResponse:
    """Exact piecewise-constant propagation of the linear dispersive ODE."""
    eps = np.asarray(epsilon, dtype=np.complex128)
    t = np.asarray(tlist, dtype=float)
    if len(t) != len(eps) + 1:
        raise ValueError("tlist must have length len(epsilon) + 1.")

    alpha_g = propagate_piecewise_constant(eps, np.diff(t), kappa=kappa, delta=delta_g)
    alpha_e = propagate_piecewise_constant(eps, np.diff(t), kappa=kappa, delta=delta_g + chi)
    root_kappa = sqrt(max(float(kappa), 0.0))
    return LinearResponse(
        tlist=t,
        epsilon=eps,
        alpha_g=alpha_g,
        alpha_e=alpha_e,
        output_g=root_kappa * alpha_g,
        output_e=root_kappa * alpha_e,
        delta_g=float(delta_g),
        delta_e=float(delta_g + chi),
        kappa=float(kappa),
    )


def propagate_piecewise_constant(
    epsilon: np.ndarray,
    durations: np.ndarray,
    *,
    kappa: float,
    delta: float,
    alpha0: complex = 0.0,
) -> np.ndarray:
    """Exact scalar propagation for `dot(alpha) = -(kappa/2 + i delta) alpha - i epsilon`."""
    eps = np.asarray(epsilon, dtype=np.complex128)
    tau = np.asarray(durations, dtype=float)
    lam = complex(0.5 * kappa + 1j * delta)
    alpha = np.empty(len(eps) + 1, dtype=np.complex128)
    alpha[0] = complex(alpha0)
    for idx, (eps_k, tau_k) in enumerate(zip(eps, tau, strict=True)):
        a_k = np.exp(-lam * tau_k)
        if abs(lam) > 1.0e-30:
            b_k = -1j * (a_k - 1.0) / lam
        else:
            b_k = -1j * tau_k
        alpha[idx + 1] = a_k * alpha[idx] + b_k * eps_k
    return alpha


def segment_coefficients(duration: float, *, kappa: float, delta: float) -> tuple[complex, complex]:
    """Return `(A, B)` for a single piecewise-constant segment."""
    lam = complex(0.5 * kappa + 1j * delta)
    a_seg = np.exp(-lam * float(duration))
    if abs(lam) > 1.0e-30:
        b_seg = -1j * (a_seg - 1.0) / lam
    else:
        b_seg = -1j * float(duration)
    return complex(a_seg), complex(b_seg)


def matched_filter_snr2(
    output_g: np.ndarray,
    output_e: np.ndarray,
    tlist: np.ndarray,
    *,
    eta: float = 1.0,
) -> float:
    """Quantum-limited matched-filter SNR² for complex output traces."""
    diff = np.asarray(output_e) - np.asarray(output_g)
    return float(max(eta, 0.0) * np.trapezoid(np.abs(diff) ** 2, np.asarray(tlist, dtype=float)))


def assignment_fidelity_from_snr2(snr2: float) -> float:
    """Binary Gaussian discrimination fidelity for matched-filter readout."""
    return float(0.5 * (1.0 + erf(sqrt(max(float(snr2), 0.0) / 2.0))))


def t1_weighted_snr2(
    output_g: np.ndarray,
    output_e: np.ndarray,
    tlist: np.ndarray,
    *,
    eta: float,
    t1: float,
) -> float:
    """Upper-bound SNR² after weighting information by qubit survival probability."""
    if not np.isfinite(t1) or t1 <= 0.0:
        return matched_filter_snr2(output_g, output_e, tlist, eta=eta)
    times = np.asarray(tlist, dtype=float)
    diff = np.asarray(output_e) - np.asarray(output_g)
    weight = np.exp(-times / float(t1))
    return float(max(eta, 0.0) * np.trapezoid(weight * np.abs(diff) ** 2, times))


def t1_limited_assignment_bound(
    output_g: np.ndarray,
    output_e: np.ndarray,
    tlist: np.ndarray,
    *,
    eta: float,
    t1: float,
) -> float:
    """Upper bound that combines matched-filter discrimination with excited-state decay."""
    times = np.asarray(tlist, dtype=float)
    snr2_eff = t1_weighted_snr2(output_g, output_e, times, eta=eta, t1=t1)
    survive = 1.0 if t1 <= 0.0 or not np.isfinite(t1) else float(np.exp(-times[-1] / t1))
    return float(survive * assignment_fidelity_from_snr2(snr2_eff) + (1.0 - survive) * 0.5)


def signal_projection(
    signal: np.ndarray,
    template_g: np.ndarray,
    template_e: np.ndarray,
    tlist: np.ndarray,
) -> tuple[float, float]:
    """Return `(projection, norm)` onto the optimal matched-filter decision axis."""
    diff = np.asarray(template_e) - np.asarray(template_g)
    centered = np.asarray(signal) - 0.5 * (np.asarray(template_e) + np.asarray(template_g))
    norm = float(np.trapezoid(np.abs(diff) ** 2, np.asarray(tlist, dtype=float)))
    proj = float(np.trapezoid(np.real(np.conjugate(diff) * centered), np.asarray(tlist, dtype=float)))
    return proj, norm


def classification_probability(
    signal: np.ndarray,
    template_g: np.ndarray,
    template_e: np.ndarray,
    tlist: np.ndarray,
    *,
    eta: float,
    target: str,
) -> float:
    """Probability that `signal` is classified as `target` under matched filtering."""
    proj, norm = signal_projection(signal, template_g, template_e, tlist)
    if norm <= 1.0e-30:
        return 0.5
    z_score = 2.0 * sqrt(max(eta, 0.0)) * proj / sqrt(norm)
    p_e = 0.5 * (1.0 + erf(z_score / sqrt(2.0)))
    return float(p_e if target == "e" else 1.0 - p_e)


def solve_two_segment_nulling(
    alpha_g: complex,
    alpha_e: complex,
    *,
    tau_3: float,
    tau_4: float,
    kappa: float,
    chi: float,
    delta_g: float,
) -> tuple[complex, complex]:
    """Solve the last two complex segment amplitudes that zero both final states."""
    a3_g, b3_g = segment_coefficients(tau_3, kappa=kappa, delta=delta_g)
    a4_g, b4_g = segment_coefficients(tau_4, kappa=kappa, delta=delta_g)
    a3_e, b3_e = segment_coefficients(tau_3, kappa=kappa, delta=delta_g + chi)
    a4_e, b4_e = segment_coefficients(tau_4, kappa=kappa, delta=delta_g + chi)

    system = np.array(
        [
            [a4_g * b3_g, b4_g],
            [a4_e * b3_e, b4_e],
        ],
        dtype=np.complex128,
    )
    rhs = -np.array(
        [
            a4_g * a3_g * alpha_g,
            a4_e * a3_e * alpha_e,
        ],
        dtype=np.complex128,
    )
    solution = np.linalg.solve(system, rhs)
    return complex(solution[0]), complex(solution[1])
