from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import qutip as qt

from common import (
    ARTIFACTS_DIR,
    MultimodeStudyConfig,
    SingleModeStudyConfig,
    bose_occupation,
    fit_decay_rate,
    make_dispersive_model,
    make_multimode_dispersive_model,
    make_multimode_dressed_model,
    make_multimode_noise,
    make_single_mode_noise,
    mode_occupation_metrics,
    pure_dephasing_ratio,
    qubit_populations,
    save_json,
    steady_state,
    simulate_idle,
)


def run_validation(
    *,
    single: SingleModeStudyConfig | None = None,
    multi: MultimodeStudyConfig | None = None,
    save_output: bool = True,
) -> dict:
    single = SingleModeStudyConfig() if single is None else single
    multi = MultimodeStudyConfig() if multi is None else multi
    checks: list[dict[str, object]] = []

    model_zero, frame_zero = make_dispersive_model(single, n_cav=12)
    rho_zero = steady_state(model_zero, noise=make_single_mode_noise(single, 0.0), frame=frame_zero)
    zero_n = mode_occupation_metrics(rho_zero, alias="cavity")["mean_n"]
    zero_pe = qubit_populations(rho_zero, n_tr=single.n_tr)["p_e"]
    checks.append(
        {
            "name": "zero_temperature_limit",
            "value": {"cavity_occupation": zero_n, "qubit_excited_population": zero_pe},
            "criterion": "n_c < 1e-8 and p_e < 1e-8",
            "passed": bool(zero_n < 1.0e-8 and zero_pe < 1.0e-8),
        }
    )

    target_temperature = 0.20
    target_nth = float(bose_occupation(target_temperature, single.omega_r))
    model_bose, frame_bose = make_dispersive_model(single, n_cav=18)
    rho_bose = steady_state(model_bose, noise=make_single_mode_noise(single, target_nth), frame=frame_bose)
    bose_n = mode_occupation_metrics(rho_bose, alias="cavity")["mean_n"]
    bose_rel_err = abs(bose_n - target_nth) / target_nth
    checks.append(
        {
            "name": "bose_consistency_0p20K",
            "value": {"target_nth": target_nth, "simulated_nth": bose_n, "relative_error": bose_rel_err},
            "criterion": "relative error < 1e-3",
            "passed": bool(bose_rel_err < 1.0e-3),
        }
    )

    hot_temperature = 2.0
    hot_nth = float(bose_occupation(hot_temperature, single.omega_r))
    model_30, frame_30 = make_dispersive_model(single, n_cav=30)
    rho_30 = steady_state(model_30, noise=make_single_mode_noise(single, hot_nth), frame=frame_30)
    model_36, frame_36 = make_dispersive_model(single, n_cav=36)
    rho_36 = steady_state(model_36, noise=make_single_mode_noise(single, hot_nth), frame=frame_36)
    n_30 = mode_occupation_metrics(rho_30, alias="cavity")["mean_n"]
    n_36 = mode_occupation_metrics(rho_36, alias="cavity")["mean_n"]
    thermometry_diff = abs(n_36 - n_30) / n_36
    checks.append(
        {
            "name": "thermometry_truncation_2K",
            "value": {"n_cav_30": n_30, "n_cav_36": n_36, "relative_difference": thermometry_diff},
            "criterion": "relative difference < 0.03",
            "passed": bool(thermometry_diff < 0.03),
        }
    )

    dephasing_temperature = 0.20
    dephasing_nth = float(bose_occupation(dephasing_temperature, single.omega_r))
    plus12 = qt.tensor((qt.basis(single.n_tr, 0) + qt.basis(single.n_tr, 1)).unit(), qt.basis(12, 0))
    plus14 = qt.tensor((qt.basis(single.n_tr, 0) + qt.basis(single.n_tr, 1)).unit(), qt.basis(14, 0))
    model_12, frame_12 = make_dispersive_model(single, n_cav=12)
    times_12, states_12 = simulate_idle(
        model_12,
        duration=30.0e-6,
        sample_dt=50.0e-9,
        frame=frame_12,
        initial_state=plus12,
        noise=make_single_mode_noise(single, dephasing_nth),
        max_step=10.0e-9,
    )
    pure_12 = np.array([pure_dephasing_ratio(state) for state in states_12], dtype=float)
    gamma_12 = fit_decay_rate(times_12, pure_12 / pure_12[0])["gamma"]

    model_14, frame_14 = make_dispersive_model(single, n_cav=14)
    times_14, states_14 = simulate_idle(
        model_14,
        duration=30.0e-6,
        sample_dt=50.0e-9,
        frame=frame_14,
        initial_state=plus14,
        noise=make_single_mode_noise(single, dephasing_nth),
        max_step=10.0e-9,
    )
    pure_14 = np.array([pure_dephasing_ratio(state) for state in states_14], dtype=float)
    gamma_14 = fit_decay_rate(times_14, pure_14 / pure_14[0])["gamma"]
    dephasing_diff = abs(gamma_14 - gamma_12) / gamma_14
    checks.append(
        {
            "name": "dephasing_truncation_0p20K",
            "value": {"gamma_12": gamma_12, "gamma_14": gamma_14, "relative_difference": dephasing_diff},
            "criterion": "relative difference < 1e-3",
            "passed": bool(dephasing_diff < 1.0e-3),
        }
    )

    low_temps = np.asarray([0.05, 0.07, 0.10], dtype=float)
    analytic_errors = []
    for temperature in low_temps:
        nth = float(bose_occupation(temperature, single.omega_r))
        model, frame = make_dispersive_model(single, n_cav=12)
        plus = qt.tensor((qt.basis(single.n_tr, 0) + qt.basis(single.n_tr, 1)).unit(), qt.basis(12, 0))
        times, states = simulate_idle(
            model,
            duration=30.0e-6,
            sample_dt=50.0e-9,
            frame=frame,
            initial_state=plus,
            noise=make_single_mode_noise(single, nth),
            max_step=10.0e-9,
        )
        pure = np.array([pure_dephasing_ratio(state) for state in states], dtype=float)
        gamma_num = fit_decay_rate(times, pure / pure[0])["gamma"]
        gamma_analytic = 4.0 * (abs(single.chi) ** 2) * nth * (nth + 1.0) / single.kappa_readout
        analytic_errors.append(abs(gamma_num - gamma_analytic) / max(gamma_num, 1.0e-12))
    max_analytic_error = float(max(analytic_errors))
    checks.append(
        {
            "name": "analytic_dephasing_scaling_low_temperature",
            "value": {"temperatures_K": low_temps.tolist(), "max_relative_error": max_analytic_error},
            "criterion": "relative error < 3.0",
            "passed": bool(max_analytic_error < 3.0),
        }
    )

    model_weak, frame_weak, omega_weak = make_multimode_dispersive_model(multi, detuning_mhz=0.0, coupling_mhz=0.0)
    weak_rho = steady_state(
        model_weak,
        noise=make_multimode_noise(
            multi,
            nth_storage=float(bose_occupation(multi.hot_storage_temperature, omega_weak)),
            nth_readout=float(bose_occupation(multi.cold_readout_temperature, multi.omega_r)),
        ),
        frame=frame_weak,
    )
    weak_readout = mode_occupation_metrics(weak_rho, alias="readout")["mean_n"]
    cold_readout = float(bose_occupation(multi.cold_readout_temperature, multi.omega_r))
    weak_limit_diff = abs(weak_readout - cold_readout) / max(cold_readout, 1.0e-12)
    checks.append(
        {
            "name": "multimode_weak_coupling_limit",
            "value": {"simulated_readout": weak_readout, "cold_readout_baseline": cold_readout, "relative_difference": weak_limit_diff},
            "criterion": "relative difference < 0.1",
            "passed": bool(weak_limit_diff < 0.1),
        }
    )

    model_m6, frame_m6, omega_m6 = make_multimode_dispersive_model(multi, detuning_mhz=-60.0, coupling_mhz=6.0, n_storage=6, n_readout=6)
    rho_m6 = steady_state(
        model_m6,
        noise=make_multimode_noise(
            multi,
            nth_storage=float(bose_occupation(multi.hot_storage_temperature, omega_m6)),
            nth_readout=float(bose_occupation(multi.cold_readout_temperature, multi.omega_r)),
        ),
        frame=frame_m6,
    )
    readout_6 = mode_occupation_metrics(rho_m6, alias="readout")["mean_n"]
    model_m8, frame_m8, omega_m8 = make_multimode_dispersive_model(multi, detuning_mhz=-60.0, coupling_mhz=6.0, n_storage=8, n_readout=8)
    rho_m8 = steady_state(
        model_m8,
        noise=make_multimode_noise(
            multi,
            nth_storage=float(bose_occupation(multi.hot_storage_temperature, omega_m8)),
            nth_readout=float(bose_occupation(multi.cold_readout_temperature, multi.omega_r)),
        ),
        frame=frame_m8,
    )
    readout_8 = mode_occupation_metrics(rho_m8, alias="readout")["mean_n"]

    model_d6, frame_d6, omega_d6 = make_multimode_dressed_model(multi, detuning_mhz=-60.0, coupling_mhz=6.0, n_storage=6, n_readout=6)
    rho_d6 = steady_state(
        model_d6,
        noise=make_multimode_noise(
            multi,
            nth_storage=float(bose_occupation(multi.hot_storage_temperature, omega_d6)),
            nth_readout=float(bose_occupation(multi.cold_readout_temperature, multi.omega_r)),
        ),
        frame=frame_d6,
    )
    pe_6 = qubit_populations(rho_d6, n_tr=2)["p_e"]
    model_d8, frame_d8, omega_d8 = make_multimode_dressed_model(multi, detuning_mhz=-60.0, coupling_mhz=6.0, n_storage=8, n_readout=8)
    rho_d8 = steady_state(
        model_d8,
        noise=make_multimode_noise(
            multi,
            nth_storage=float(bose_occupation(multi.hot_storage_temperature, omega_d8)),
            nth_readout=float(bose_occupation(multi.cold_readout_temperature, multi.omega_r)),
        ),
        frame=frame_d8,
    )
    pe_8 = qubit_populations(rho_d8, n_tr=2)["p_e"]
    readout_diff = abs(readout_8 - readout_6) / readout_8
    excitation_diff = abs(pe_8 - pe_6) / pe_8
    checks.append(
        {
            "name": "multimode_truncation_representative_point",
            "value": {
                "readout_6": readout_6,
                "readout_8": readout_8,
                "readout_relative_difference": readout_diff,
                "pe_6": pe_6,
                "pe_8": pe_8,
                "excitation_relative_difference": excitation_diff,
            },
            "criterion": "both relative differences < 0.03",
            "passed": bool(readout_diff < 0.03 and excitation_diff < 0.03),
        }
    )

    payload = {
        "checks": checks,
        "passed_all": bool(all(bool(item["passed"]) for item in checks)),
    }
    if save_output:
        save_json(ARTIFACTS_DIR / "validation_summary.json", payload)

    return payload


def main() -> None:
    payload = run_validation(save_output=True)

    for item in payload["checks"]:
        label = "PASS" if item["passed"] else "FAIL"
        print(f"[{label}] {item['name']}: {item['criterion']}")


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    main()
