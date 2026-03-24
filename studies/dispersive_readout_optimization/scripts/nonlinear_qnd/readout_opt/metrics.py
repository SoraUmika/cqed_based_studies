"""Metrics for readout discrimination, cavity emptying, and QND-aware behavior."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .bounds import (
    assignment_fidelity_from_snr2,
    classification_probability,
    matched_filter_snr2,
    t1_limited_assignment_bound,
)


@dataclass(frozen=True)
class ReadoutMetrics:
    """Common metric surface shared by the linear and richer-model evaluations."""

    snr2_ideal: float
    fidelity_ideal: float
    fidelity_eta: float
    fidelity_t1_bound: float
    residual_photons: float
    residual_after_wait: float
    peak_photons: float
    qnd_preservation: float
    qnd_defect: float
    leakage: float
    transition_error: float
    baseline_transition: float
    measurement_induced_transition: float
    measurement_induced_qnd_defect: float
    repeat_consistency: float
    score_balanced: float
    score_info: float
    score_emptying: float


def cavity_expectation_from_quadratures(x_trace: np.ndarray, p_trace: np.ndarray) -> np.ndarray:
    """Recover `<a>` from `x = a + a^†` and `p = -i(a - a^†)` traces."""
    return 0.5 * (np.asarray(x_trace, dtype=np.complex128) + 1j * np.asarray(p_trace, dtype=np.complex128))


def residual_photons(alpha_g: np.ndarray, alpha_e: np.ndarray) -> float:
    """Mean end-of-pulse residual photons over the two basis preparations."""
    return float(0.5 * (abs(alpha_g[-1]) ** 2 + abs(alpha_e[-1]) ** 2))


def peak_photons(alpha_g: np.ndarray, alpha_e: np.ndarray) -> float:
    """Maximum conditioned photon occupancy during the pulse."""
    return float(max(np.max(np.abs(alpha_g) ** 2), np.max(np.abs(alpha_e) ** 2)))


def balanced_score(
    fidelity_eta: float,
    *,
    residual: float,
    peak: float,
    measurement_induced_transition: float,
    measurement_induced_qnd_defect: float,
    leakage: float,
    repeat_consistency: float,
    n_crit: float,
) -> float:
    """Primary scalarization used for richer-model optimization."""
    peak_ratio = 0.0 if not np.isfinite(n_crit) or n_crit <= 0.0 else max(0.0, peak / n_crit - 0.45)
    return float(
        fidelity_eta
        - 0.08 * residual
        - 0.20 * peak_ratio
        - 0.85 * measurement_induced_qnd_defect
        - 1.15 * measurement_induced_transition
        - 1.20 * leakage
        - 0.50 * (1.0 - repeat_consistency)
    )


def info_only_score(fidelity_eta: float) -> float:
    return float(fidelity_eta)


def emptying_score(
    fidelity_eta: float,
    *,
    residual: float,
    measurement_induced_qnd_defect: float,
    leakage: float,
) -> float:
    return float(
        fidelity_eta
        - 0.20 * residual
        - 0.65 * measurement_induced_qnd_defect
        - 1.00 * leakage
    )


def build_metrics(
    *,
    output_g: np.ndarray,
    output_e: np.ndarray,
    tlist: np.ndarray,
    eta: float,
    t1: float,
    alpha_g: np.ndarray,
    alpha_e: np.ndarray,
    residual_after_wait_value: float,
    qnd_preservation: float,
    leakage: float,
    repeat_consistency: float,
    n_crit: float,
    baseline_transition: float = 0.0,
) -> ReadoutMetrics:
    """Construct the common metric bundle from signal and preservation data."""
    snr2_ideal = matched_filter_snr2(output_g, output_e, tlist, eta=1.0)
    fidelity_ideal = assignment_fidelity_from_snr2(snr2_ideal)
    fidelity_eta = assignment_fidelity_from_snr2(matched_filter_snr2(output_g, output_e, tlist, eta=eta))
    fidelity_t1 = t1_limited_assignment_bound(output_g, output_e, tlist, eta=eta, t1=t1)
    residual_value = residual_photons(alpha_g, alpha_e)
    peak_value = peak_photons(alpha_g, alpha_e)
    qnd_defect = max(0.0, 1.0 - qnd_preservation)
    transition_error = max(0.0, 1.0 - qnd_preservation - leakage)
    measurement_induced_transition = max(0.0, transition_error - baseline_transition)
    measurement_induced_qnd_defect = max(0.0, qnd_defect - baseline_transition)
    return ReadoutMetrics(
        snr2_ideal=snr2_ideal,
        fidelity_ideal=fidelity_ideal,
        fidelity_eta=fidelity_eta,
        fidelity_t1_bound=fidelity_t1,
        residual_photons=residual_value,
        residual_after_wait=residual_after_wait_value,
        peak_photons=peak_value,
        qnd_preservation=qnd_preservation,
        qnd_defect=qnd_defect,
        leakage=leakage,
        transition_error=transition_error,
        baseline_transition=float(max(baseline_transition, 0.0)),
        measurement_induced_transition=measurement_induced_transition,
        measurement_induced_qnd_defect=measurement_induced_qnd_defect,
        repeat_consistency=repeat_consistency,
        score_balanced=balanced_score(
            fidelity_eta,
            residual=residual_value,
            peak=peak_value,
            measurement_induced_transition=measurement_induced_transition,
            measurement_induced_qnd_defect=measurement_induced_qnd_defect,
            leakage=leakage,
            repeat_consistency=repeat_consistency,
            n_crit=n_crit,
        ),
        score_info=info_only_score(fidelity_eta),
        score_emptying=emptying_score(
            fidelity_eta,
            residual=residual_value,
            measurement_induced_qnd_defect=measurement_induced_qnd_defect,
            leakage=leakage,
        ),
    )


def repeated_consistency_from_signals(
    *,
    first_g: np.ndarray,
    first_e: np.ndarray,
    second_g: np.ndarray,
    second_e: np.ndarray,
    tlist: np.ndarray,
    eta: float,
) -> float:
    """Average probability that the second readout agrees with the prepared basis state."""
    p_g = classification_probability(second_g, first_g, first_e, tlist, eta=eta, target="g")
    p_e = classification_probability(second_e, first_g, first_e, tlist, eta=eta, target="e")
    return float(0.5 * (p_g + p_e))
