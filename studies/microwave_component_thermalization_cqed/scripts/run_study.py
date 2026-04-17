from __future__ import annotations

import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import qutip as qt
from joblib import Parallel, delayed

from common import (
    ARTIFACTS_DIR,
    DATA_DIR,
    FIGURES_DIR,
    REPORT_DIR,
    SingleModeStudyConfig,
    MultimodeStudyConfig,
    bose_occupation,
    cavity_temperature_resolution_proxy,
    configure_matplotlib,
    dataclass_payload,
    elapsed_s,
    ensure_directories,
    fit_decay_rate,
    fit_step_response,
    format_seconds,
    make_dispersive_model,
    make_dressed_single_mode_model,
    make_multimode_dispersive_model,
    make_multimode_dressed_model,
    make_multimode_noise,
    make_single_mode_noise,
    mode_occupation_metrics,
    payload_with_runtime,
    pure_dephasing_ratio,
    qubit_coherence,
    qubit_populations,
    qubit_temperature_resolution_proxy,
    readout_noise_penalty,
    save_figure,
    save_json,
    spectroscopy_broadening_mhz,
    spectroscopy_shift_mhz,
    steady_state,
    simulate_idle,
)


def thermometer_payload(config: SingleModeStudyConfig | None = None) -> dict:
    start = time.perf_counter()
    config = SingleModeStudyConfig() if config is None else config
    temperatures = np.asarray(config.temperature_grid, dtype=float)
    nth_values = bose_occupation(temperatures, config.omega_r)

    dispersive_model, dispersive_frame = make_dispersive_model(config, n_cav=config.n_cav_steady)
    dressed_model, dressed_frame = make_dressed_single_mode_model(config, n_cav=config.n_cav_steady)

    cavity_n = []
    cavity_var = []
    qubit_excited = []
    qubit_f = []
    dressed_cavity_n = []
    spectroscopy_shift = []
    spectroscopy_width = []
    readout_penalty = []

    for nth in nth_values:
        rho_disp = steady_state(dispersive_model, noise=make_single_mode_noise(config, nth), frame=dispersive_frame)
        mode_metrics = mode_occupation_metrics(rho_disp, alias="cavity")
        cavity_n.append(mode_metrics["mean_n"])
        cavity_var.append(mode_metrics["variance_n"])
        spectroscopy_shift.append(spectroscopy_shift_mhz(config.chi, mode_metrics["mean_n"]))
        spectroscopy_width.append(spectroscopy_broadening_mhz(config.chi, mode_metrics["variance_n"]))
        readout_penalty.append(readout_noise_penalty(mode_metrics["mean_n"]))

        rho_dressed = steady_state(dressed_model, noise=make_single_mode_noise(config, nth, dressed=True), frame=dressed_frame)
        qubit_metrics = qubit_populations(rho_dressed, n_tr=config.n_tr)
        dressed_metrics = mode_occupation_metrics(rho_dressed, alias="readout")
        qubit_excited.append(qubit_metrics["p_e"])
        qubit_f.append(qubit_metrics["p_f"])
        dressed_cavity_n.append(dressed_metrics["mean_n"])

    cavity_n = np.asarray(cavity_n, dtype=float)
    cavity_var = np.asarray(cavity_var, dtype=float)
    qubit_excited = np.asarray(qubit_excited, dtype=float)
    qubit_f = np.asarray(qubit_f, dtype=float)
    dressed_cavity_n = np.asarray(dressed_cavity_n, dtype=float)
    spectroscopy_shift = np.asarray(spectroscopy_shift, dtype=float)
    spectroscopy_width = np.asarray(spectroscopy_width, dtype=float)
    readout_penalty = np.asarray(readout_penalty, dtype=float)

    cavity_temp_resolution = cavity_temperature_resolution_proxy(temperatures, cavity_n)
    qubit_temp_resolution = qubit_temperature_resolution_proxy(temperatures, qubit_excited)
    cavity_dT = np.gradient(cavity_n, temperatures, edge_order=1)
    qubit_dT = np.gradient(qubit_excited, temperatures, edge_order=1)

    thermometer_preference = np.where(qubit_temp_resolution < cavity_temp_resolution, "qubit_excited_population", "cavity_occupation")
    crossover_index = next(
        (index for index in range(1, len(temperatures)) if thermometer_preference[index] != thermometer_preference[index - 1]),
        None,
    )
    crossover_temperature = None if crossover_index is None else float(temperatures[crossover_index])

    payload = {
        "config": dataclass_payload(config),
        "temperatures_K": temperatures.tolist(),
        "nth_values": nth_values.tolist(),
        "dispersive_cavity_occupation": cavity_n.tolist(),
        "dispersive_cavity_variance": cavity_var.tolist(),
        "dressed_qubit_excited_population": qubit_excited.tolist(),
        "dressed_qubit_f_population": qubit_f.tolist(),
        "dressed_readout_occupation": dressed_cavity_n.tolist(),
        "spectroscopy_shift_MHz": spectroscopy_shift.tolist(),
        "spectroscopy_width_MHz": spectroscopy_width.tolist(),
        "readout_visibility_penalty": readout_penalty.tolist(),
        "cavity_temperature_resolution_proxy_K": cavity_temp_resolution.tolist(),
        "qubit_temperature_resolution_proxy_K": qubit_temp_resolution.tolist(),
        "cavity_responsivity_per_K": cavity_dT.tolist(),
        "qubit_responsivity_per_K": qubit_dT.tolist(),
        "preferred_thermometer": thermometer_preference.tolist(),
        "crossover_temperature_K": crossover_temperature,
        "summary": {
            "lowest_temperature_with_qubit_excitation_above_1e-4": float(
                next((temp for temp, p_e in zip(temperatures, qubit_excited) if p_e >= 1.0e-4), temperatures[-1])
            ),
            "temperature_where_readout_penalty_drops_below_0p9": float(
                next((temp for temp, penalty in zip(temperatures, readout_penalty) if penalty <= 0.9), temperatures[-1])
            ),
            "temperature_where_spectroscopy_width_exceeds_0p1_MHz": float(
                next((temp for temp, width in zip(temperatures, spectroscopy_width) if width >= 0.1), temperatures[-1])
            ),
        },
    }

    figure, axes = plt.subplots(2, 2, figsize=(10.5, 8.0), constrained_layout=True)

    axes[0, 0].plot(temperatures, cavity_n, marker="o", label="Dispersive cavity occupation")
    axes[0, 0].plot(temperatures, qubit_excited, marker="s", label="Dressed qubit excited population")
    axes[0, 0].set_xlabel("Bath temperature (K)")
    axes[0, 0].set_ylabel("Observable")
    axes[0, 0].set_title("Calibration curves")
    axes[0, 0].legend(loc="upper left")

    axes[0, 1].plot(temperatures, spectroscopy_shift, marker="o", label="Mean shift")
    axes[0, 1].plot(temperatures, spectroscopy_width, marker="s", label="Thermal broadening")
    axes[0, 1].set_xlabel("Bath temperature (K)")
    axes[0, 1].set_ylabel("Spectroscopy signature (MHz)")
    axes[0, 1].set_title("Spectroscopy proxy")
    axes[0, 1].legend(loc="upper left")

    axes[1, 0].plot(temperatures, cavity_temp_resolution, marker="o", label="Cavity occupation proxy")
    axes[1, 0].plot(temperatures, qubit_temp_resolution, marker="s", label="Qubit population proxy")
    axes[1, 0].set_yscale("log")
    axes[1, 0].set_xlabel("Bath temperature (K)")
    axes[1, 0].set_ylabel("Temperature resolution proxy (K)")
    axes[1, 0].set_title("Lower is better")
    axes[1, 0].legend(loc="upper right")

    axes[1, 1].plot(temperatures, cavity_dT, marker="o", label="d<n>/dT")
    axes[1, 1].plot(temperatures, qubit_dT, marker="s", label="dPe/dT")
    axes[1, 1].plot(temperatures, 1.0 - readout_penalty, marker="^", label="Readout noise increase")
    axes[1, 1].set_xlabel("Bath temperature (K)")
    axes[1, 1].set_ylabel("Responsivity / penalty")
    axes[1, 1].set_title("Sensitivity tradeoffs")
    axes[1, 1].legend(loc="upper left")

    save_figure(figure, "thermometer_calibration_curves")
    return payload_with_runtime(payload, runtime_s=elapsed_s(start))


def dephasing_payload(config: SingleModeStudyConfig | None = None) -> dict:
    start = time.perf_counter()
    config = SingleModeStudyConfig() if config is None else config
    temperatures = np.asarray(config.dephasing_temperatures, dtype=float)
    nth_values = bose_occupation(temperatures, config.omega_r)

    dispersive_model, dispersive_frame = make_dispersive_model(config, n_cav=config.n_cav_dynamic)
    dressed_model, dressed_frame = make_dressed_single_mode_model(config, n_cav=config.n_cav_dynamic)

    coherence_traces: dict[str, list[float]] = {}
    heating_traces: dict[str, list[float]] = {}

    gamma_total = []
    gamma_pure = []
    gamma_analytic = []
    gamma_up = []
    p_ss = []
    t2_total = []

    duration_coherence = 30.0e-6
    duration_heating = 12.0e-6

    qubit_plus = qt.tensor((qt.basis(config.n_tr, 0) + qt.basis(config.n_tr, 1)).unit(), qt.basis(config.n_cav_dynamic, 0))
    qubit_ground = qt.tensor(qt.basis(config.n_tr, 0), qt.basis(config.n_cav_dynamic, 0))

    for temperature, nth in zip(temperatures, nth_values):
        disp_noise = make_single_mode_noise(config, nth)
        times, states = simulate_idle(
            dispersive_model,
            duration=duration_coherence,
            sample_dt=50.0e-9,
            frame=dispersive_frame,
            initial_state=qubit_plus,
            noise=disp_noise,
            max_step=10.0e-9,
        )
        coherence = np.array([2.0 * qubit_coherence(state) for state in states], dtype=float)
        pure = np.array([pure_dephasing_ratio(state) for state in states], dtype=float)
        coherence_traces[f"{temperature:.2f} K"] = coherence.tolist()

        total_fit = fit_decay_rate(times, coherence / max(coherence[0], 1.0e-12))
        pure_fit = fit_decay_rate(times, pure / max(pure[0], 1.0e-12))
        gamma_total.append(total_fit["gamma"])
        gamma_pure.append(pure_fit["gamma"])
        gamma_analytic.append(4.0 * (abs(config.chi) ** 2) * float(nth) * float(nth + 1.0) / config.kappa_readout)

        dressed_noise = make_single_mode_noise(config, nth, dressed=True)
        heat_times, heat_states = simulate_idle(
            dressed_model,
            duration=duration_heating,
            sample_dt=25.0e-9,
            frame=dressed_frame,
            initial_state=qubit_ground,
            noise=dressed_noise,
            max_step=10.0e-9,
        )
        p_e_trace = np.array([qubit_populations(state, n_tr=config.n_tr)["p_e"] for state in heat_states], dtype=float)
        heating_traces[f"{temperature:.2f} K"] = p_e_trace.tolist()
        step_fit = fit_step_response(heat_times, p_e_trace)
        p_ss_value = max(step_fit["y_inf"], 0.0)
        gamma_up_value = 0.0 if not np.isfinite(step_fit["tau"]) or step_fit["tau"] <= 0.0 else p_ss_value / step_fit["tau"]
        gamma_up.append(gamma_up_value)
        p_ss.append(p_ss_value)
        t2_total.append(float(np.inf if gamma_up_value + pure_fit["gamma"] <= 0.0 else 1.0 / (pure_fit["gamma"] + 0.5 * gamma_up_value)))

    gamma_total = np.asarray(gamma_total, dtype=float)
    gamma_pure = np.asarray(gamma_pure, dtype=float)
    gamma_analytic = np.asarray(gamma_analytic, dtype=float)
    gamma_up = np.asarray(gamma_up, dtype=float)
    p_ss = np.asarray(p_ss, dtype=float)
    t2_total = np.asarray(t2_total, dtype=float)

    tolerable_temperature = float(
        next((temp for temp, t2 in zip(temperatures, t2_total) if t2 < config.coherence_target), temperatures[-1])
    )

    payload = {
        "config": dataclass_payload(config),
        "temperatures_K": temperatures.tolist(),
        "nth_values": nth_values.tolist(),
        "gamma_total_per_s": gamma_total.tolist(),
        "gamma_pure_per_s": gamma_pure.tolist(),
        "gamma_analytic_per_s": gamma_analytic.tolist(),
        "gamma_up_per_s": gamma_up.tolist(),
        "steady_state_qubit_excited_population": p_ss.tolist(),
        "t2_total_s": t2_total.tolist(),
        "coherence_traces": coherence_traces,
        "heating_traces": heating_traces,
        "coherence_target_s": config.coherence_target,
        "tolerable_temperature_K": tolerable_temperature,
        "summary": {
            "max_relative_error_low_occupation": float(
                np.max(np.abs(gamma_pure[:4] - gamma_analytic[:4]) / np.clip(gamma_pure[:4], 1.0e-12, None))
            ),
            "dominant_mechanism_above_0p15K": "pure_dephasing" if gamma_pure[4] > 2.0 * gamma_up[4] else "mixed",
        },
    }

    figure, axes = plt.subplots(2, 2, figsize=(11.0, 8.0), constrained_layout=True)
    representative_keys = [f"{temperatures[index]:.2f} K" for index in (0, 3, 5, 7)]
    time_axis = np.linspace(0.0, duration_coherence * 1.0e6, len(next(iter(coherence_traces.values()))))
    heat_axis = np.linspace(0.0, duration_heating * 1.0e6, len(next(iter(heating_traces.values()))))

    for key in representative_keys:
        axes[0, 0].plot(time_axis, coherence_traces[key], label=key)
        axes[0, 1].plot(heat_axis, heating_traces[key], label=key)

    axes[0, 0].set_xlabel("Idle time (us)")
    axes[0, 0].set_ylabel("Normalized coherence")
    axes[0, 0].set_title("Ramsey-like decay")
    axes[0, 0].legend(loc="upper right")

    axes[0, 1].set_xlabel("Idle time (us)")
    axes[0, 1].set_ylabel("Qubit excited population")
    axes[0, 1].set_title("Thermal population heating")
    axes[0, 1].legend(loc="upper left")

    axes[1, 0].plot(nth_values, gamma_pure / 1.0e6, marker="o", label="Numerical pure dephasing")
    axes[1, 0].plot(nth_values, gamma_analytic / 1.0e6, marker="s", linestyle="--", label="Analytic low-occupation estimate")
    axes[1, 0].plot(nth_values, gamma_up / 1.0e6, marker="^", label="Heating rate")
    axes[1, 0].set_xlabel("Thermal occupation")
    axes[1, 0].set_ylabel("Rate (MHz)")
    axes[1, 0].set_title("Dephasing and heating rates")
    axes[1, 0].legend(loc="upper left")

    axes[1, 1].plot(temperatures, 1.0e6 * t2_total, marker="o", label="Thermal-limited T2")
    axes[1, 1].axhline(1.0e6 * config.coherence_target, color="black", linestyle="--", label="20 us target")
    axes[1, 1].set_xlabel("Bath temperature (K)")
    axes[1, 1].set_ylabel("Coherence time (us)")
    axes[1, 1].set_title("Tolerable temperature window")
    axes[1, 1].legend(loc="upper right")

    save_figure(figure, "dephasing_and_coherence_limits")
    return payload_with_runtime(payload, runtime_s=elapsed_s(start))


def multimode_dispersive_point(config: MultimodeStudyConfig, detuning_mhz: float, coupling_mhz: float) -> tuple[float, float, float]:
    model, frame, storage_omega = make_multimode_dispersive_model(config, detuning_mhz=detuning_mhz, coupling_mhz=coupling_mhz)
    nth_storage = float(bose_occupation(config.hot_storage_temperature, storage_omega))
    nth_readout = float(bose_occupation(config.cold_readout_temperature, config.omega_r))
    rho = steady_state(model, noise=make_multimode_noise(config, nth_storage=nth_storage, nth_readout=nth_readout), frame=frame)
    storage_metrics = mode_occupation_metrics(rho, alias="storage")
    readout_metrics = mode_occupation_metrics(rho, alias="readout")
    gamma_proxy = 4.0 * (abs(config.chi_r) ** 2) * readout_metrics["variance_n"] / config.kappa_readout
    return readout_metrics["mean_n"], gamma_proxy, storage_metrics["mean_n"]


def multimode_dressed_point(config: MultimodeStudyConfig, detuning_mhz: float, coupling_mhz: float) -> tuple[float, float]:
    model, frame, storage_omega = make_multimode_dressed_model(config, detuning_mhz=detuning_mhz, coupling_mhz=coupling_mhz)
    nth_storage = float(bose_occupation(config.hot_storage_temperature, storage_omega))
    nth_readout = float(bose_occupation(config.cold_readout_temperature, config.omega_r))
    rho = steady_state(model, noise=make_multimode_noise(config, nth_storage=nth_storage, nth_readout=nth_readout), frame=frame)
    qubit_metrics = qubit_populations(rho, n_tr=2)
    readout_metrics = mode_occupation_metrics(rho, alias="readout")
    return qubit_metrics["p_e"], readout_metrics["mean_n"]


def multimode_payload(config: MultimodeStudyConfig | None = None) -> dict:
    start = time.perf_counter()
    config = MultimodeStudyConfig() if config is None else config
    detuning_grid = np.asarray(config.detuning_grid_mhz, dtype=float)
    coupling_grid = np.asarray(config.coupling_grid_mhz, dtype=float)

    grid_points = [(detuning, coupling) for detuning in detuning_grid for coupling in coupling_grid]
    dispersive_results = Parallel(n_jobs=-1, prefer="threads")(
        delayed(multimode_dispersive_point)(config, detuning, coupling) for detuning, coupling in grid_points
    )
    dressed_results = Parallel(n_jobs=-1, prefer="threads")(
        delayed(multimode_dressed_point)(config, detuning, coupling) for detuning, coupling in grid_points
    )

    readout_heating = np.asarray([item[0] for item in dispersive_results], dtype=float).reshape(len(detuning_grid), len(coupling_grid))
    gamma_proxy = np.asarray([item[1] for item in dispersive_results], dtype=float).reshape(len(detuning_grid), len(coupling_grid))
    storage_occupation = np.asarray([item[2] for item in dispersive_results], dtype=float).reshape(len(detuning_grid), len(coupling_grid))
    dressed_excitation = np.asarray([item[0] for item in dressed_results], dtype=float).reshape(len(detuning_grid), len(coupling_grid))
    dressed_readout = np.asarray([item[1] for item in dressed_results], dtype=float).reshape(len(detuning_grid), len(coupling_grid))
    safe_mask = (
        (readout_heating <= config.safe_readout_threshold)
        & (dressed_excitation <= config.safe_excitation_threshold)
        & (gamma_proxy <= 1.0 / 20.0e-6)
    )

    sensitivity_scan = []
    for linewidth_ns in config.kappa_storage_sensitivity_ns:
        linewidth = 1.0 / (linewidth_ns * 1.0e-9)
        model, frame, storage_omega = make_multimode_dispersive_model(config, detuning_mhz=60.0, coupling_mhz=6.0)
        rho = steady_state(
            model,
            noise=make_multimode_noise(
                config,
                nth_storage=float(bose_occupation(config.hot_storage_temperature, storage_omega)),
                nth_readout=float(bose_occupation(config.cold_readout_temperature, config.omega_r)),
                kappa_storage=linewidth,
            ),
            frame=frame,
        )
        metrics = mode_occupation_metrics(rho, alias="readout")
        sensitivity_scan.append(
            {
                "storage_linewidth_ns": float(linewidth_ns),
                "readout_occupation": metrics["mean_n"],
                "dephasing_proxy_per_s": 4.0 * (abs(config.chi_r) ** 2) * metrics["variance_n"] / config.kappa_readout,
            }
        )

    payload = {
        "config": dataclass_payload(config),
        "detuning_grid_MHz": detuning_grid.tolist(),
        "coupling_grid_MHz": coupling_grid.tolist(),
        "induced_readout_occupation": readout_heating.tolist(),
        "dephasing_proxy_per_s": gamma_proxy.tolist(),
        "storage_occupation": storage_occupation.tolist(),
        "dressed_qubit_excited_population": dressed_excitation.tolist(),
        "dressed_readout_occupation": dressed_readout.tolist(),
        "safe_mask": safe_mask.astype(int).tolist(),
        "safe_fraction": float(np.mean(safe_mask)),
        "linewidth_sensitivity": sensitivity_scan,
        "summary": {
            "most_dangerous_point": {
                "detuning_MHz": float(detuning_grid[np.unravel_index(np.argmax(dressed_excitation), dressed_excitation.shape)[0]]),
                "coupling_MHz": float(coupling_grid[np.unravel_index(np.argmax(dressed_excitation), dressed_excitation.shape)[1]]),
                "qubit_excited_population": float(np.max(dressed_excitation)),
            }
        },
    }

    figure, axes = plt.subplots(1, 3, figsize=(13.5, 4.0), constrained_layout=True)
    mesh_x, mesh_y = np.meshgrid(coupling_grid, detuning_grid)

    image0 = axes[0].pcolormesh(mesh_x, mesh_y, readout_heating, shading="nearest")
    axes[0].set_title("Readout heating from hot auxiliary mode")
    axes[0].set_xlabel("Storage-readout exchange (MHz)")
    axes[0].set_ylabel("Storage detuning from readout (MHz)")
    figure.colorbar(image0, ax=axes[0], label="Mean readout occupation")

    image1 = axes[1].pcolormesh(mesh_x, mesh_y, gamma_proxy / 1.0e6, shading="nearest")
    axes[1].set_title("Dispersive dephasing proxy")
    axes[1].set_xlabel("Storage-readout exchange (MHz)")
    axes[1].set_ylabel("Storage detuning from readout (MHz)")
    figure.colorbar(image1, ax=axes[1], label="Proxy rate (MHz)")

    image2 = axes[2].pcolormesh(mesh_x, mesh_y, dressed_excitation, shading="nearest")
    axes[2].contour(mesh_x, mesh_y, safe_mask.astype(float), levels=[0.5], colors="white", linewidths=1.5)
    axes[2].set_title("Dressed qubit heating with safe boundary")
    axes[2].set_xlabel("Storage-readout exchange (MHz)")
    axes[2].set_ylabel("Storage detuning from readout (MHz)")
    figure.colorbar(image2, ax=axes[2], label="Qubit excited population")

    save_figure(figure, "multimode_heating_maps")
    return payload_with_runtime(payload, runtime_s=elapsed_s(start))


def transient_payload(config: SingleModeStudyConfig | None = None) -> dict:
    start = time.perf_counter()
    config = SingleModeStudyConfig() if config is None else config
    dispersive_model, dispersive_frame = make_dispersive_model(config, n_cav=config.n_cav_dynamic)
    dressed_model, dressed_frame = make_dressed_single_mode_model(config, n_cav=config.n_cav_dynamic)

    cold_nth = float(bose_occupation(config.cold_temperature, config.omega_r))
    rho_cold_disp = steady_state(dispersive_model, noise=make_single_mode_noise(config, cold_nth), frame=dispersive_frame)
    rho_cold_dressed = steady_state(dressed_model, noise=make_single_mode_noise(config, cold_nth, dressed=True), frame=dressed_frame)

    transient_results = {}
    figure, axes = plt.subplots(1, 2, figsize=(10.5, 4.0), constrained_layout=True)

    for hot_temperature in config.transient_temperatures:
        hot_nth = float(bose_occupation(hot_temperature, config.omega_r))
        times_disp, states_disp = simulate_idle(
            dispersive_model,
            duration=6.0e-6,
            sample_dt=20.0e-9,
            frame=dispersive_frame,
            initial_state=rho_cold_disp,
            noise=make_single_mode_noise(config, hot_nth),
            max_step=5.0e-9,
        )
        n_c_trace = np.array([mode_occupation_metrics(state, alias="cavity")["mean_n"] for state in states_disp], dtype=float)
        fit_n_c = fit_step_response(times_disp, n_c_trace)

        times_dressed, states_dressed = simulate_idle(
            dressed_model,
            duration=6.0e-6,
            sample_dt=20.0e-9,
            frame=dressed_frame,
            initial_state=rho_cold_dressed,
            noise=make_single_mode_noise(config, hot_nth, dressed=True),
            max_step=5.0e-9,
        )
        p_e_trace = np.array([qubit_populations(state, n_tr=config.n_tr)["p_e"] for state in states_dressed], dtype=float)
        fit_p_e = fit_step_response(times_dressed, p_e_trace)

        transient_results[f"{hot_temperature:.2f} K"] = {
            "bath_occupation": hot_nth,
            "times_s": times_disp.tolist(),
            "cavity_occupation_trace": n_c_trace.tolist(),
            "qubit_excited_trace": p_e_trace.tolist(),
            "cavity_response_tau_s": fit_n_c["tau"],
            "qubit_response_tau_s": fit_p_e["tau"],
        }

        axes[0].plot(times_disp * 1.0e6, n_c_trace, label=f"{hot_temperature:.2f} K")
        axes[1].plot(times_dressed * 1.0e6, p_e_trace, label=f"{hot_temperature:.2f} K")

    axes[0].set_xlabel("Time after bath step (us)")
    axes[0].set_ylabel("Cavity occupation")
    axes[0].set_title("Intrinsic quantum response of cavity")
    axes[0].legend(loc="lower right")

    axes[1].set_xlabel("Time after bath step (us)")
    axes[1].set_ylabel("Qubit excited population")
    axes[1].set_title("Intrinsic quantum response of qubit thermometer")
    axes[1].legend(loc="lower right")

    save_figure(figure, "transient_thermal_step_response")
    payload = {
        "config": dataclass_payload(config),
        "cold_temperature_K": config.cold_temperature,
        "results": transient_results,
        "summary": {
            "captured_time_scale_statement": "The simulated response times reflect only the internal cQED Lindblad dynamics after the bath occupation is changed instantaneously.",
            "missing_time_scale_statement": "Any slower VTS or hardware thermalization time must come from macroscopic heat flow and boundary resistance outside the quantum model.",
        },
    }
    return payload_with_runtime(payload, runtime_s=elapsed_s(start))


def run_full_study(
    *,
    single_config: SingleModeStudyConfig | None = None,
    multimode_config: MultimodeStudyConfig | None = None,
    save_outputs: bool = True,
) -> dict:
    ensure_directories()
    configure_matplotlib()
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    single = SingleModeStudyConfig() if single_config is None else single_config
    multimode = MultimodeStudyConfig() if multimode_config is None else multimode_config

    thermometer = thermometer_payload(single)
    if save_outputs:
        save_json(ARTIFACTS_DIR / "thermometer_summary.json", thermometer)

    dephasing = dephasing_payload(single)
    if save_outputs:
        save_json(ARTIFACTS_DIR / "dephasing_summary.json", dephasing)

    multimode_payload_result = multimode_payload(multimode)
    if save_outputs:
        save_json(ARTIFACTS_DIR / "multimode_summary.json", multimode_payload_result)

    transient = transient_payload(single)
    if save_outputs:
        save_json(ARTIFACTS_DIR / "transient_summary.json", transient)

    combined = {
        "thermometer": thermometer["summary"],
        "dephasing": {
            "tolerable_temperature_K": dephasing["tolerable_temperature_K"],
            "lowest_t2_us": float(min(np.asarray(dephasing["t2_total_s"], dtype=float)) * 1.0e6),
        },
        "multimode": {
            "safe_fraction": multimode_payload_result["safe_fraction"],
            "most_dangerous_point": multimode_payload_result["summary"]["most_dangerous_point"],
        },
        "transient": {
            key: {
                "cavity_response_tau_us": float(value["cavity_response_tau_s"] * 1.0e6),
                "qubit_response_tau_us": float(value["qubit_response_tau_s"] * 1.0e6),
            }
            for key, value in transient["results"].items()
        },
        "figure_stems": [
            "thermometer_calibration_curves",
            "dephasing_and_coherence_limits",
            "multimode_heating_maps",
            "transient_thermal_step_response",
        ],
    }
    if save_outputs:
        save_json(DATA_DIR / "study_summary.json", combined)

    return {
        "thermometer": thermometer,
        "dephasing": dephasing,
        "multimode": multimode_payload_result,
        "transient": transient,
        "summary": combined,
    }


def main() -> None:
    results = run_full_study(save_outputs=True)

    print("Thermometer runtime:", format_seconds(results["thermometer"]["runtime_s"]))
    print("Dephasing runtime:", format_seconds(results["dephasing"]["runtime_s"]))
    print("Multimode runtime:", format_seconds(results["multimode"]["runtime_s"]))
    print("Transient runtime:", format_seconds(results["transient"]["runtime_s"]))
    print("Saved summary to", DATA_DIR / "study_summary.json")


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    main()
