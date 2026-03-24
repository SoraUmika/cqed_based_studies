"""Estimate Purcell-limited T1 and readout-induced backaction for the A2 study."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import qutip as qt

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import runtime_compat  # noqa: F401

from common import (
    DATA_DIR,
    KAPPA_R,
    OMEGA_Q,
    OMEGA_R,
    TWO_PI,
    build_nominal_readout_chain,
)

DETUNING_GHZ = np.array([0.5, 1.0, 2.0, 3.0, 5.0], dtype=float)
READOUT_EPSILON_MHZ = np.array([0.5, 1.0, 2.0, 5.0], dtype=float)
BACKACTION_DURATIONS_NS = np.array([50.0, 100.0, 200.0, 400.0, 800.0], dtype=float)


def main() -> None:
    omega_q_grid = OMEGA_R + TWO_PI * DETUNING_GHZ * 1.0e9

    chain_no_filter = build_nominal_readout_chain(
        epsilon=TWO_PI * 1.0e6,
        include_filter=False,
    )
    chain_with_filter = build_nominal_readout_chain(
        epsilon=TWO_PI * 1.0e6,
        include_filter=True,
        filter_bandwidth=KAPPA_R,
    )

    purcell_rate_no_filter = np.array([chain_no_filter.purcell_rate(omega_q) for omega_q in omega_q_grid], dtype=float)
    purcell_rate_with_filter = np.array([chain_with_filter.purcell_rate(omega_q) for omega_q in omega_q_grid], dtype=float)
    purcell_t1_no_filter = np.array([chain_no_filter.purcell_limited_t1(omega_q) for omega_q in omega_q_grid], dtype=float)
    purcell_t1_with_filter = np.array([chain_with_filter.purcell_limited_t1(omega_q) for omega_q in omega_q_grid], dtype=float)

    rho_plus = ((qt.basis(2, 0) + qt.basis(2, 1)).unit()).proj()
    coherence_matrix = np.zeros((len(READOUT_EPSILON_MHZ), len(BACKACTION_DURATIONS_NS)), dtype=float)
    excited_population_matrix = np.zeros_like(coherence_matrix)
    gamma_meas_values = np.zeros(len(READOUT_EPSILON_MHZ), dtype=float)

    for eps_index, epsilon_mhz in enumerate(READOUT_EPSILON_MHZ):
        chain = build_nominal_readout_chain(
            epsilon=TWO_PI * epsilon_mhz * 1.0e6,
            include_filter=True,
            filter_bandwidth=KAPPA_R,
        )
        gamma_meas_values[eps_index] = chain.gamma_meas()
        for dur_index, duration_ns in enumerate(BACKACTION_DURATIONS_NS):
            rho_after = chain.apply_backaction(
                rho_plus,
                omega_q=OMEGA_Q,
                duration=float(duration_ns) * 1.0e-9,
                include_measurement_dephasing=True,
                include_purcell_relaxation=True,
            )
            rho_arr = rho_after.full()
            coherence_matrix[eps_index, dur_index] = float(abs(rho_arr[0, 1]))
            excited_population_matrix[eps_index, dur_index] = float(np.real(rho_arr[1, 1]))

    output_path = DATA_DIR / "purcell_and_backaction.npz"
    np.savez(
        output_path,
        detuning_ghz=DETUNING_GHZ,
        readout_epsilon_mhz=READOUT_EPSILON_MHZ,
        backaction_durations_ns=BACKACTION_DURATIONS_NS,
        purcell_rate_no_filter=purcell_rate_no_filter,
        purcell_rate_with_filter=purcell_rate_with_filter,
        purcell_t1_no_filter_s=purcell_t1_no_filter,
        purcell_t1_with_filter_s=purcell_t1_with_filter,
        gamma_meas_rad_s=gamma_meas_values,
        coherence_after_backaction=coherence_matrix,
        excited_population_after_backaction=excited_population_matrix,
        kappa_r_rad_s=KAPPA_R,
    )
    print(f"Saved {output_path}")


if __name__ == "__main__":
    main()