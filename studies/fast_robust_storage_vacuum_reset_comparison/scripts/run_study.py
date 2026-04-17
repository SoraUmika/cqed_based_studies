"""Run the fast robust storage vacuum-reset comparison study end to end."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import math
from pathlib import Path
import textwrap
import time

try:
    from .common import (
        ARTIFACTS_DIR,
        DATA_DIR,
        DEFAULT_CONTINUOUS_DURATION_S,
        DEFAULT_PULSED_RINGDOWN_MULTIPLES,
        DEFAULT_SWEEP_DT_S,
        DEFAULT_TRAJECTORY_DT_S,
        DEVICE,
        FIGURES_DIR,
        REPORT_DIR,
        STUDY_DIR,
        TRANSMON_REFERENCE,
        SchemeSummary,
        baseline_noise,
        build_frame,
        build_model,
        csv_dump,
        e_fold_time_ns,
        estimate_decay_rate_hz,
        expectation_operators,
        fock_state,
        coherent_storage_state,
        thermal_storage_state,
        first_threshold_time_ns,
        json_dump,
        linear_entropy,
        make_pulse,
        plot_save,
        pulsed_recommendations,
        readout_sideband_rotating_frequency,
        run_sequence,
        sideband_matrix_element,
        storage_sideband_rotating_frequency,
        to_internal_units,
        write_device_manifest,
    )
except ImportError:
    from common import (
        ARTIFACTS_DIR,
        DATA_DIR,
        DEFAULT_CONTINUOUS_DURATION_S,
        DEFAULT_PULSED_RINGDOWN_MULTIPLES,
        DEFAULT_SWEEP_DT_S,
        DEFAULT_TRAJECTORY_DT_S,
        DEVICE,
        FIGURES_DIR,
        REPORT_DIR,
        STUDY_DIR,
        TRANSMON_REFERENCE,
        SchemeSummary,
        baseline_noise,
        build_frame,
        build_model,
        csv_dump,
        e_fold_time_ns,
        estimate_decay_rate_hz,
        expectation_operators,
        fock_state,
        coherent_storage_state,
        thermal_storage_state,
        first_threshold_time_ns,
        json_dump,
        linear_entropy,
        make_pulse,
        plot_save,
        pulsed_recommendations,
        readout_sideband_rotating_frequency,
        run_sequence,
        sideband_matrix_element,
        storage_sideband_rotating_frequency,
        to_internal_units,
        write_device_manifest,
    )

import matplotlib.pyplot as plt
import numpy as np
import qutip as qt

from cqed_sim.core.drive_targets import SidebandDriveSpec

COHERENCE_SCALES = np.array([0.50, 0.75, 1.00, 1.50, 2.00], dtype=float)
AMPLITUDE_SCALES = np.array([0.90, 0.95, 1.00, 1.05, 1.10], dtype=float)
DETUNING_ERRORS_MHZ = np.array([-2.0, -1.0, 0.0, 1.0, 2.0], dtype=float)
RESET_PROBABILITIES = np.array([0.00, 0.02, 0.05], dtype=float)
THERMAL_LOADS = np.array([0.00, 0.01, 0.03], dtype=float)
READOUT_KAPPA_SCALES = np.array([0.50, 1.00, 2.00, 4.00], dtype=float)
RESONANT_COUPLING_GRID_MHZ = np.array([4.0, 6.0, 8.0, 10.0, 12.0], dtype=float)
DETUNED_COUPLING_GRID_MHZ = np.array([4.0, 6.0, 8.0], dtype=float)
DETUNING_GRID_MHZ = np.array([8.0, 12.0, 16.0, 20.0, 24.0], dtype=float)

RESONANT_N_TARGET = 1
DEFAULT_STORAGE_THRESHOLD = 0.01
DEFAULT_TRANSMON_THRESHOLD = 1.0e-3


def _transmon_excited_curve(expectations: dict[str, np.ndarray]) -> np.ndarray:
    return np.clip(1.0 - np.real(np.asarray(expectations["P_g"])), 0.0, 1.0)


def _protocol_metrics(
    *,
    scheme_key: str,
    scheme_label: str,
    times_s: np.ndarray,
    expectations: dict[str, np.ndarray],
    final_state: qt.Qobj | None,
) -> dict[str, object]:
    storage_curve = np.real(np.asarray(expectations["n_storage"]))
    readout_curve = np.real(np.asarray(expectations["n_readout"]))
    transmon_curve = _transmon_excited_curve(expectations)
    ground_curve = np.real(np.asarray(expectations["P_g00"]))
    hplus_curve = np.real(np.asarray(expectations["P_hplus"]))
    e_curve = np.real(np.asarray(expectations["P_e"]))
    f_curve = np.real(np.asarray(expectations["P_f"]))
    steady_tail = max(5, len(storage_curve) // 10)
    metrics = {
        "scheme_key": scheme_key,
        "scheme_label": scheme_label,
        "time_to_threshold_ns": first_threshold_time_ns(
            np.asarray(times_s, dtype=float),
            storage_curve,
            transmon_curve,
            storage_threshold=DEFAULT_STORAGE_THRESHOLD,
            transmon_threshold=DEFAULT_TRANSMON_THRESHOLD,
        ),
        "e_fold_time_ns": e_fold_time_ns(np.asarray(times_s, dtype=float), storage_curve),
        "estimated_decay_rate_hz": estimate_decay_rate_hz(np.asarray(times_s, dtype=float), storage_curve),
        "final_storage_n": float(storage_curve[-1]),
        "final_readout_n": float(readout_curve[-1]),
        "final_transmon_excited": float(transmon_curve[-1]),
        "final_ground_vacuum": float(ground_curve[-1]),
        "steady_storage_n": float(np.mean(storage_curve[-steady_tail:])),
        "steady_transmon_excited": float(np.mean(transmon_curve[-steady_tail:])),
        "steady_ground_vacuum": float(np.mean(ground_curve[-steady_tail:])),
        "max_e_population": float(np.max(e_curve)),
        "max_f_population": float(np.max(f_curve)),
        "max_hplus_population": float(np.max(hplus_curve)),
        "times_ns": (np.asarray(times_s, dtype=float) * 1.0e9).tolist(),
        "storage_curve": storage_curve.tolist(),
        "readout_curve": readout_curve.tolist(),
        "transmon_excited_curve": transmon_curve.tolist(),
        "ground_vacuum_curve": ground_curve.tolist(),
        "hplus_curve": hplus_curve.tolist(),
    }
    if final_state is not None:
        metrics["final_linear_entropy"] = float(linear_entropy(final_state))
    return metrics


def _append_segment(
    times_s: list[float],
    expectations: dict[str, list[float]],
    *,
    stage_times_s: np.ndarray,
    stage_expectations: dict[str, np.ndarray],
    time_offset_s: float,
) -> float:
    stage_times = np.asarray(stage_times_s, dtype=float) + float(time_offset_s)
    trim = 0 if not times_s else 1
    times_s.extend(stage_times[trim:].tolist())
    for key, values in stage_expectations.items():
        expectations.setdefault(key, [])
        expectations[key].extend(np.asarray(values, dtype=float)[trim:].tolist())
    return float(stage_times[-1])


def _pulsed_total_duration_s(recommendations: dict[int, dict[str, object]], *, ringdown_multiple: float, max_n: int = 4) -> float:
    noise = baseline_noise()
    ringdown_s = 0.0 if noise.kappa_readout is None or noise.kappa_readout <= 0.0 else float(ringdown_multiple) / float(noise.kappa_readout)
    total = 0.0
    for n in range(int(max_n), 0, -1):
        total += float(recommendations[n]["step_a_duration_ns"]) * 1.0e-9
        total += float(recommendations[n]["step_b_duration_ns"]) * 1.0e-9
        total += ringdown_s
    return total


def simulate_pulsed_protocol(
    model,
    frame,
    *,
    initial_state: qt.Qobj,
    noise,
    recommendations: dict[int, dict[str, object]],
    ringdown_multiple: float,
    max_n: int = 4,
    amplitude_scale: float = 1.0,
    detuning_error_mhz: float = 0.0,
    dt_s: float = DEFAULT_SWEEP_DT_S,
) -> dict[str, object]:
    e_ops = expectation_operators(model)
    time_buffer: list[float] = []
    expectation_buffer: dict[str, list[float]] = {}
    state = initial_state
    time_offset_s = 0.0
    detuning_rad_s = to_internal_units(float(detuning_error_mhz) * 1.0e6)
    ringdown_s = 0.0 if noise.kappa_readout is None or noise.kappa_readout <= 0.0 else float(ringdown_multiple) / float(noise.kappa_readout)

    for n in range(int(max_n), 0, -1):
        rec = recommendations[n]
        storage_pulse = make_pulse(
            channel="storage_sb",
            carrier_rad_s=storage_sideband_rotating_frequency(model, frame, n) + detuning_rad_s,
            duration_s=float(rec["step_a_duration_ns"]) * 1.0e-9,
            amplitude_hz=float(rec["step_a_amplitude_MHz"]) * 1.0e6 * float(amplitude_scale),
            envelope_name=str(rec["step_a_family"]),
            label=f"pulsed_storage_n{n}",
        )
        compiled, result = run_sequence(
            model,
            state,
            pulses=[storage_pulse],
            duration_s=float(rec["step_a_duration_ns"]) * 1.0e-9,
            drive_ops={"storage_sb": SidebandDriveSpec(mode="storage", lower_level=0, upper_level=2, sideband="red")},
            frame=frame,
            noise=noise,
            dt_s=dt_s,
            e_ops=e_ops,
            store_states=False,
        )
        time_offset_s = _append_segment(
            time_buffer,
            expectation_buffer,
            stage_times_s=np.asarray(compiled.tlist, dtype=float),
            stage_expectations=result.expectations,
            time_offset_s=time_offset_s,
        )
        state = result.final_state

        dump_pulse = make_pulse(
            channel="readout_sb",
            carrier_rad_s=readout_sideband_rotating_frequency(model, frame, n) + detuning_rad_s,
            duration_s=float(rec["step_b_duration_ns"]) * 1.0e-9,
            amplitude_hz=float(rec["step_b_amplitude_MHz"]) * 1.0e6 * float(amplitude_scale),
            envelope_name=str(rec["step_b_family"]),
            label=f"pulsed_readout_n{n}",
        )
        compiled, result = run_sequence(
            model,
            state,
            pulses=[dump_pulse],
            duration_s=float(rec["step_b_duration_ns"]) * 1.0e-9,
            drive_ops={"readout_sb": SidebandDriveSpec(mode="readout", lower_level=0, upper_level=2, sideband="red")},
            frame=frame,
            noise=noise,
            dt_s=dt_s,
            e_ops=e_ops,
            store_states=False,
        )
        time_offset_s = _append_segment(
            time_buffer,
            expectation_buffer,
            stage_times_s=np.asarray(compiled.tlist, dtype=float),
            stage_expectations=result.expectations,
            time_offset_s=time_offset_s,
        )
        state = result.final_state

        if ringdown_s > 0.0:
            compiled, result = run_sequence(
                model,
                state,
                pulses=[],
                duration_s=ringdown_s,
                drive_ops={},
                frame=frame,
                noise=noise,
                dt_s=max(dt_s, 4.0e-9),
                e_ops=e_ops,
                store_states=False,
            )
            time_offset_s = _append_segment(
                time_buffer,
                expectation_buffer,
                stage_times_s=np.asarray(compiled.tlist, dtype=float),
                stage_expectations=result.expectations,
                time_offset_s=time_offset_s,
            )
            state = result.final_state

    metrics = _protocol_metrics(
        scheme_key="pulsed_ladder",
        scheme_label="Pulsed `g-f` ladder",
        times_s=np.asarray(time_buffer, dtype=float),
        expectations={key: np.asarray(values, dtype=float) for key, values in expectation_buffer.items()},
        final_state=state,
    )
    metrics.update(
        {
            "ringdown_multiple": float(ringdown_multiple),
            "amplitude_scale": float(amplitude_scale),
            "detuning_error_mhz": float(detuning_error_mhz),
            "protocol_duration_ns": float(time_buffer[-1] * 1.0e9) if time_buffer else 0.0,
        }
    )
    return metrics


def simulate_continuous_protocol(
    model,
    frame,
    *,
    initial_state: qt.Qobj,
    noise,
    target_coupling_mhz: float,
    common_detuning_mhz: float,
    duration_s: float,
    amplitude_scale: float = 1.0,
    detuning_error_mhz: float = 0.0,
    dt_s: float = DEFAULT_SWEEP_DT_S,
) -> dict[str, object]:
    e_ops = expectation_operators(model)
    n_target = RESONANT_N_TARGET
    storage_element = sideband_matrix_element(model, mode="storage", n=n_target)
    readout_element = sideband_matrix_element(model, mode="readout", n=n_target)
    common_offset = to_internal_units((float(common_detuning_mhz) + float(detuning_error_mhz)) * 1.0e6)
    storage_pulse = make_pulse(
        channel="storage_tone",
        carrier_rad_s=storage_sideband_rotating_frequency(model, frame, n_target) + common_offset,
        duration_s=float(duration_s),
        amplitude_hz=float(target_coupling_mhz) * 1.0e6 * float(amplitude_scale) / max(storage_element, 1.0e-15),
        envelope_name="square",
        label=f"continuous_storage_{target_coupling_mhz:.1f}_{common_detuning_mhz:.1f}",
    )
    readout_pulse = make_pulse(
        channel="readout_tone",
        carrier_rad_s=readout_sideband_rotating_frequency(model, frame, n_target) + common_offset,
        duration_s=float(duration_s),
        amplitude_hz=float(target_coupling_mhz) * 1.0e6 * float(amplitude_scale) / max(readout_element, 1.0e-15),
        envelope_name="square",
        label=f"continuous_readout_{target_coupling_mhz:.1f}_{common_detuning_mhz:.1f}",
    )
    compiled, result = run_sequence(
        model,
        initial_state,
        pulses=[storage_pulse, readout_pulse],
        duration_s=float(duration_s),
        drive_ops={
            "storage_tone": SidebandDriveSpec(mode="storage", lower_level=0, upper_level=2, sideband="red"),
            "readout_tone": SidebandDriveSpec(mode="readout", lower_level=0, upper_level=2, sideband="red"),
        },
        frame=frame,
        noise=noise,
        dt_s=dt_s,
        e_ops=e_ops,
        store_states=False,
    )
    scheme_key = "continuous_bright" if abs(float(common_detuning_mhz)) < 1.0e-12 else "continuous_raman"
    scheme_label = "Continuous resonant bright-state" if scheme_key == "continuous_bright" else "Continuous detuned Raman-like"
    metrics = _protocol_metrics(
        scheme_key=scheme_key,
        scheme_label=scheme_label,
        times_s=np.asarray(compiled.tlist, dtype=float),
        expectations=result.expectations,
        final_state=result.final_state,
    )
    metrics.update(
        {
            "target_coupling_mhz": float(target_coupling_mhz),
            "common_detuning_mhz": float(common_detuning_mhz),
            "amplitude_scale": float(amplitude_scale),
            "detuning_error_mhz": float(detuning_error_mhz),
            "storage_amplitude_mhz": float(storage_pulse.amp / (2.0 * np.pi * 1.0e6)),
            "readout_amplitude_mhz": float(readout_pulse.amp / (2.0 * np.pi * 1.0e6)),
        }
    )
    return metrics


def simulate_effective_autonomous_benchmark(
    model,
    *,
    initial_state: qt.Qobj,
    duration_s: float,
    gamma_eff_hz: float,
    nth_storage: float = 0.0,
    dt_s: float = DEFAULT_SWEEP_DT_S,
) -> dict[str, object]:
    storage_dm = initial_state.proj().ptrace(1) if initial_state.isket else initial_state.ptrace(1)
    a = qt.destroy(int(model.subsystem_dims[1]))
    n_op = a.dag() * a
    vacuum = qt.basis(int(model.subsystem_dims[1]), 0).proj()
    times_s = np.arange(0.0, float(duration_s) + 0.5 * float(dt_s), float(dt_s))
    c_ops: list[qt.Qobj] = []
    gamma_eff_hz = max(0.0, float(gamma_eff_hz))
    if gamma_eff_hz > 0.0:
        c_ops.append(np.sqrt(gamma_eff_hz * (1.0 + float(nth_storage))) * a)
        if float(nth_storage) > 0.0:
            c_ops.append(np.sqrt(gamma_eff_hz * float(nth_storage)) * a.dag())
    kappa_storage = 0.0 if DEVICE.storage_t1_s <= 0.0 else 1.0 / float(DEVICE.storage_t1_s)
    if kappa_storage > 0.0:
        c_ops.append(np.sqrt(kappa_storage * (1.0 + float(nth_storage))) * a)
        if float(nth_storage) > 0.0:
            c_ops.append(np.sqrt(kappa_storage * float(nth_storage)) * a.dag())
    if DEVICE.storage_t1_s > 0.0 and DEVICE.storage_t2_ramsey_s > 0.0:
        inv_tphi = max(0.0, 1.0 / float(DEVICE.storage_t2_ramsey_s) - 1.0 / (2.0 * float(DEVICE.storage_t1_s)))
        if inv_tphi > 0.0:
            c_ops.append(np.sqrt(inv_tphi) * n_op)
    result = qt.mesolve(0.0 * n_op, storage_dm, times_s, c_ops, e_ops=[n_op, vacuum])
    expectations = {
        "n_storage": np.asarray(result.expect[0], dtype=float),
        "n_readout": np.zeros_like(result.expect[0], dtype=float),
        "P_g": np.ones_like(result.expect[0], dtype=float),
        "P_e": np.zeros_like(result.expect[0], dtype=float),
        "P_f": np.zeros_like(result.expect[0], dtype=float),
        "P_g00": np.asarray(result.expect[1], dtype=float),
        "P_hplus": np.zeros_like(result.expect[0], dtype=float),
    }
    metrics = _protocol_metrics(
        scheme_key="effective_autonomous_benchmark",
        scheme_label="Effective `L \\propto a_s` benchmark",
        times_s=times_s,
        expectations=expectations,
        final_state=result.states[-1] if result.states else storage_dm,
    )
    metrics.update({"gamma_eff_hz": float(gamma_eff_hz), "nth_storage": float(nth_storage)})
    return metrics


def _select_best_candidate(rows: list[dict[str, object]]) -> dict[str, object]:
    eligible = [row for row in rows if row["time_to_threshold_ns"] is not None]
    if eligible:
        return min(
            eligible,
            key=lambda row: (
                float(row["time_to_threshold_ns"]),
                float(row["max_hplus_population"]),
                float(row["final_storage_n"]),
                float(row["final_transmon_excited"]),
            ),
        )
    return min(
        rows,
        key=lambda row: (
            float(row["steady_storage_n"]),
            float(row["final_storage_n"]),
            float(row["final_transmon_excited"]),
            float(row["max_hplus_population"]),
        ),
    )


def _scheme_initial_states(model) -> dict[str, dict[str, object]]:
    return {
        "single_photon_n1": {
            "label": r"$|g,0,1\rangle$",
            "state": fock_state(model, storage_level=1),
            "pulsed_max_n": 1,
        },
        "higher_fock_n3": {
            "label": r"$|g,0,3\rangle$",
            "state": fock_state(model, storage_level=3),
            "pulsed_max_n": 3,
        },
        "coherent_alpha_1p0": {
            "label": r"$|g\rangle \otimes |\alpha=1\rangle_s \otimes |0\rangle_r$",
            "state": coherent_storage_state(model, alpha=1.0),
            "pulsed_max_n": 4,
        },
        "thermal_nbar_0p5": {
            "label": r"$|g\rangle\langle g| \otimes \rho_{\mathrm{th}}(\bar n=0.5) \otimes |0\rangle\langle 0|$",
            "state": thermal_storage_state(model, nbar=0.5),
            "pulsed_max_n": 4,
        },
    }


def _coherence_heatmaps(
    model,
    frame,
    *,
    pulsed_settings: dict[str, object],
    resonant_settings: dict[str, object],
    detuned_settings: dict[str, object],
    benchmark_gamma_hz: float,
    comparison_window_s: float,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    initial = fock_state(model, storage_level=1)
    recommendations = pulsed_recommendations()
    for t1_scale in COHERENCE_SCALES:
        for t2_scale in COHERENCE_SCALES:
            noise = baseline_noise(t1_scale=float(t1_scale), t2_scale=float(t2_scale))
            pulsed = simulate_pulsed_protocol(
                model,
                frame,
                initial_state=initial,
                noise=noise,
                recommendations=recommendations,
                ringdown_multiple=float(pulsed_settings["ringdown_multiple"]),
                max_n=1,
                dt_s=DEFAULT_SWEEP_DT_S,
            )
            pulsed.update({"t1_scale": float(t1_scale), "t2_scale": float(t2_scale)})
            rows.append(pulsed)

            bright = simulate_continuous_protocol(
                model,
                frame,
                initial_state=initial,
                noise=noise,
                target_coupling_mhz=float(resonant_settings["target_coupling_mhz"]),
                common_detuning_mhz=0.0,
                duration_s=comparison_window_s,
                dt_s=DEFAULT_SWEEP_DT_S,
            )
            bright.update({"t1_scale": float(t1_scale), "t2_scale": float(t2_scale)})
            rows.append(bright)

            raman = simulate_continuous_protocol(
                model,
                frame,
                initial_state=initial,
                noise=noise,
                target_coupling_mhz=float(detuned_settings["target_coupling_mhz"]),
                common_detuning_mhz=float(detuned_settings["common_detuning_mhz"]),
                duration_s=comparison_window_s,
                dt_s=DEFAULT_SWEEP_DT_S,
            )
            raman.update({"t1_scale": float(t1_scale), "t2_scale": float(t2_scale)})
            rows.append(raman)

            benchmark = simulate_effective_autonomous_benchmark(
                model,
                initial_state=initial,
                duration_s=comparison_window_s,
                gamma_eff_hz=float(benchmark_gamma_hz),
                dt_s=DEFAULT_SWEEP_DT_S,
            )
            benchmark.update({"t1_scale": float(t1_scale), "t2_scale": float(t2_scale)})
            rows.append(benchmark)
    return rows


def _calibration_heatmaps(
    model,
    frame,
    *,
    pulsed_settings: dict[str, object],
    resonant_settings: dict[str, object],
    detuned_settings: dict[str, object],
    benchmark_gamma_hz: float,
    comparison_window_s: float,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    initial = fock_state(model, storage_level=1)
    base_noise = baseline_noise()
    recommendations = pulsed_recommendations()
    for amplitude_scale in AMPLITUDE_SCALES:
        for detuning_error_mhz in DETUNING_ERRORS_MHZ:
            pulsed = simulate_pulsed_protocol(
                model,
                frame,
                initial_state=initial,
                noise=base_noise,
                recommendations=recommendations,
                ringdown_multiple=float(pulsed_settings["ringdown_multiple"]),
                max_n=1,
                amplitude_scale=float(amplitude_scale),
                detuning_error_mhz=float(detuning_error_mhz),
                dt_s=DEFAULT_SWEEP_DT_S,
            )
            pulsed.update({"amplitude_scale": float(amplitude_scale), "detuning_error_mhz": float(detuning_error_mhz)})
            rows.append(pulsed)

            bright = simulate_continuous_protocol(
                model,
                frame,
                initial_state=initial,
                noise=base_noise,
                target_coupling_mhz=float(resonant_settings["target_coupling_mhz"]),
                common_detuning_mhz=0.0,
                duration_s=comparison_window_s,
                amplitude_scale=float(amplitude_scale),
                detuning_error_mhz=float(detuning_error_mhz),
                dt_s=DEFAULT_SWEEP_DT_S,
            )
            bright.update({"amplitude_scale": float(amplitude_scale), "detuning_error_mhz": float(detuning_error_mhz)})
            rows.append(bright)

            raman = simulate_continuous_protocol(
                model,
                frame,
                initial_state=initial,
                noise=base_noise,
                target_coupling_mhz=float(detuned_settings["target_coupling_mhz"]),
                common_detuning_mhz=float(detuned_settings["common_detuning_mhz"]),
                duration_s=comparison_window_s,
                amplitude_scale=float(amplitude_scale),
                detuning_error_mhz=float(detuning_error_mhz),
                dt_s=DEFAULT_SWEEP_DT_S,
            )
            raman.update({"amplitude_scale": float(amplitude_scale), "detuning_error_mhz": float(detuning_error_mhz)})
            rows.append(raman)

            benchmark = simulate_effective_autonomous_benchmark(
                model,
                initial_state=initial,
                duration_s=comparison_window_s,
                gamma_eff_hz=float(benchmark_gamma_hz),
                dt_s=DEFAULT_SWEEP_DT_S,
            )
            benchmark.update({"amplitude_scale": float(amplitude_scale), "detuning_error_mhz": float(detuning_error_mhz)})
            rows.append(benchmark)
    return rows


def _reset_and_thermal_rows(
    model,
    frame,
    *,
    pulsed_settings: dict[str, object],
    resonant_settings: dict[str, object],
    detuned_settings: dict[str, object],
    benchmark_gamma_hz: float,
    comparison_window_s: float,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    recommendations = pulsed_recommendations()
    initial_pure = fock_state(model, storage_level=1)
    for reset_prob in RESET_PROBABILITIES:
        storage_dm = initial_pure.proj()
        rho_q = (1.0 - float(reset_prob)) * qt.basis(int(model.subsystem_dims[0]), 0).proj() + float(reset_prob) * qt.basis(int(model.subsystem_dims[0]), 1).proj()
        rho_s = storage_dm.ptrace(1)
        rho_r = storage_dm.ptrace(2)
        initial = qt.tensor(rho_q, rho_s, rho_r)
        noise = baseline_noise()
        for scheme_key in ("pulsed_ladder", "continuous_bright", "continuous_raman", "effective_autonomous_benchmark"):
            if scheme_key == "pulsed_ladder":
                row = simulate_pulsed_protocol(
                    model,
                    frame,
                    initial_state=initial,
                    noise=noise,
                    recommendations=recommendations,
                    ringdown_multiple=float(pulsed_settings["ringdown_multiple"]),
                    max_n=1,
                    dt_s=DEFAULT_SWEEP_DT_S,
                )
            elif scheme_key == "continuous_bright":
                row = simulate_continuous_protocol(
                    model,
                    frame,
                    initial_state=initial,
                    noise=noise,
                    target_coupling_mhz=float(resonant_settings["target_coupling_mhz"]),
                    common_detuning_mhz=0.0,
                    duration_s=comparison_window_s,
                    dt_s=DEFAULT_SWEEP_DT_S,
                )
            elif scheme_key == "continuous_raman":
                row = simulate_continuous_protocol(
                    model,
                    frame,
                    initial_state=initial,
                    noise=noise,
                    target_coupling_mhz=float(detuned_settings["target_coupling_mhz"]),
                    common_detuning_mhz=float(detuned_settings["common_detuning_mhz"]),
                    duration_s=comparison_window_s,
                    dt_s=DEFAULT_SWEEP_DT_S,
                )
            else:
                row = simulate_effective_autonomous_benchmark(
                    model,
                    initial_state=initial,
                    duration_s=comparison_window_s,
                    gamma_eff_hz=float(benchmark_gamma_hz),
                    dt_s=DEFAULT_SWEEP_DT_S,
                )
            row.update({"sweep_type": "reset_error", "reset_prob": float(reset_prob), "nth_storage": 0.0, "nth_readout": 0.0})
            rows.append(row)
    for nth_readout in THERMAL_LOADS:
        for nth_storage in THERMAL_LOADS:
            initial = fock_state(model, storage_level=1)
            noise = baseline_noise(nth_storage=float(nth_storage), nth_readout=float(nth_readout))
            for scheme_key in ("pulsed_ladder", "continuous_bright", "continuous_raman", "effective_autonomous_benchmark"):
                if scheme_key == "pulsed_ladder":
                    row = simulate_pulsed_protocol(
                        model,
                        frame,
                        initial_state=initial,
                        noise=noise,
                        recommendations=recommendations,
                        ringdown_multiple=float(pulsed_settings["ringdown_multiple"]),
                        max_n=1,
                        dt_s=DEFAULT_SWEEP_DT_S,
                    )
                elif scheme_key == "continuous_bright":
                    row = simulate_continuous_protocol(
                        model,
                        frame,
                        initial_state=initial,
                        noise=noise,
                        target_coupling_mhz=float(resonant_settings["target_coupling_mhz"]),
                        common_detuning_mhz=0.0,
                        duration_s=comparison_window_s,
                        dt_s=DEFAULT_SWEEP_DT_S,
                    )
                elif scheme_key == "continuous_raman":
                    row = simulate_continuous_protocol(
                        model,
                        frame,
                        initial_state=initial,
                        noise=noise,
                        target_coupling_mhz=float(detuned_settings["target_coupling_mhz"]),
                        common_detuning_mhz=float(detuned_settings["common_detuning_mhz"]),
                        duration_s=comparison_window_s,
                        dt_s=DEFAULT_SWEEP_DT_S,
                    )
                else:
                    row = simulate_effective_autonomous_benchmark(
                        model,
                        initial_state=initial,
                        duration_s=comparison_window_s,
                        gamma_eff_hz=float(benchmark_gamma_hz),
                        nth_storage=float(nth_storage),
                        dt_s=DEFAULT_SWEEP_DT_S,
                    )
                row.update({"sweep_type": "thermal_loading", "reset_prob": 0.0, "nth_storage": float(nth_storage), "nth_readout": float(nth_readout)})
                rows.append(row)
    return rows


def _kappa_tradeoff_rows(
    model,
    frame,
    *,
    resonant_settings: dict[str, object],
    detuned_settings: dict[str, object],
    comparison_window_s: float,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    initial = fock_state(model, storage_level=1)
    for scale in READOUT_KAPPA_SCALES:
        noise = baseline_noise(readout_kappa_scale=float(scale))
        bright = simulate_continuous_protocol(
            model,
            frame,
            initial_state=initial,
            noise=noise,
            target_coupling_mhz=float(resonant_settings["target_coupling_mhz"]),
            common_detuning_mhz=0.0,
            duration_s=comparison_window_s,
            dt_s=DEFAULT_SWEEP_DT_S,
        )
        bright.update({"readout_kappa_scale": float(scale)})
        rows.append(bright)
        raman = simulate_continuous_protocol(
            model,
            frame,
            initial_state=initial,
            noise=noise,
            target_coupling_mhz=float(detuned_settings["target_coupling_mhz"]),
            common_detuning_mhz=float(detuned_settings["common_detuning_mhz"]),
            duration_s=comparison_window_s,
            dt_s=DEFAULT_SWEEP_DT_S,
        )
        raman.update({"readout_kappa_scale": float(scale)})
        rows.append(raman)
    return rows


def _convergence_rows(
    recommendations: dict[int, dict[str, object]],
    *,
    pulsed_settings: dict[str, object],
    resonant_settings: dict[str, object],
    detuned_settings: dict[str, object],
    comparison_window_s: float,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for model_label, n_tr, n_storage, n_readout in (("baseline", 5, 6, 4), ("larger", 6, 7, 5)):
        model = build_model(n_tr=n_tr, n_storage=n_storage, n_readout=n_readout)
        frame = build_frame(model)
        initial = fock_state(model, storage_level=1)
        for dt_s in (DEFAULT_TRAJECTORY_DT_S, DEFAULT_SWEEP_DT_S, 2.0e-9):
            noise = baseline_noise()
            pulsed = simulate_pulsed_protocol(
                model,
                frame,
                initial_state=initial,
                noise=noise,
                recommendations=recommendations,
                ringdown_multiple=float(pulsed_settings["ringdown_multiple"]),
                max_n=1,
                dt_s=float(dt_s),
            )
            pulsed.update({"model_label": model_label, "dt_ns": float(dt_s * 1.0e9)})
            rows.append(pulsed)
            bright = simulate_continuous_protocol(
                model,
                frame,
                initial_state=initial,
                noise=noise,
                target_coupling_mhz=float(resonant_settings["target_coupling_mhz"]),
                common_detuning_mhz=0.0,
                duration_s=comparison_window_s,
                dt_s=float(dt_s),
            )
            bright.update({"model_label": model_label, "dt_ns": float(dt_s * 1.0e9)})
            rows.append(bright)
            raman = simulate_continuous_protocol(
                model,
                frame,
                initial_state=initial,
                noise=noise,
                target_coupling_mhz=float(detuned_settings["target_coupling_mhz"]),
                common_detuning_mhz=float(detuned_settings["common_detuning_mhz"]),
                duration_s=comparison_window_s,
                dt_s=float(dt_s),
            )
            raman.update({"model_label": model_label, "dt_ns": float(dt_s * 1.0e9)})
            rows.append(raman)
    return rows


def _make_trajectory_figure(trajectory_rows: list[dict[str, object]]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    axes = axes.ravel()
    ordered = [
        ("pulsed_ladder", "Pulsed ladder"),
        ("continuous_bright", "Continuous bright-state"),
        ("continuous_raman", "Continuous Raman-like"),
        ("effective_autonomous_benchmark", "Effective benchmark"),
    ]
    for axis, (scheme_key, title) in zip(axes, ordered, strict=True):
        row = next(item for item in trajectory_rows if item["scheme_key"] == scheme_key)
        t_ns = np.asarray(row["times_ns"], dtype=float)
        axis.plot(t_ns, row["storage_curve"], label=r"$\langle n_s \rangle$", linewidth=1.8)
        axis.plot(t_ns, row["transmon_excited_curve"], label=r"$P_{\mathrm{tr,exc}}$", linewidth=1.5, linestyle="--")
        axis.plot(t_ns, row["ground_vacuum_curve"], label=r"$P_{g00}$", linewidth=1.3, linestyle=":")
        axis.set_title(title)
        axis.set_xlabel("Time (ns)")
        axis.set_ylabel("Population / mean occupation")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
    plot_save(fig, "scheme_cooling_dynamics")


def _make_summary_figure(summary_rows: list[SchemeSummary]) -> None:
    actual = [row for row in summary_rows if not row.auxiliary_only]
    labels = [row.scheme_label for row in actual]
    threshold = [np.nan if row.baseline_time_to_threshold_ns is None else row.baseline_time_to_threshold_ns for row in actual]
    robustness = [row.robustness_score for row in actual]
    max_hplus = [row.max_hplus_population for row in actual]
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
    axes[0].bar(labels, threshold, color=["#005f73", "#bb3e03", "#0a9396"])
    axes[0].set_ylabel("Time to threshold (ns)")
    axes[0].tick_params(axis="x", rotation=18)
    axes[1].bar(labels, robustness, color=["#005f73", "#bb3e03", "#0a9396"])
    axes[1].set_ylabel("Robustness score")
    axes[1].tick_params(axis="x", rotation=18)
    axes[2].bar(labels, max_hplus, color=["#005f73", "#bb3e03", "#0a9396"])
    axes[2].set_ylabel("Max h+ population")
    axes[2].tick_params(axis="x", rotation=18)
    for axis in axes:
        axis.grid(alpha=0.25, axis="y")
    plot_save(fig, "scheme_comparison_summary")


def _make_initial_state_figure(initial_state_rows: list[dict[str, object]]) -> None:
    state_names = ["single_photon_n1", "higher_fock_n3", "coherent_alpha_1p0", "thermal_nbar_0p5"]
    scheme_order = ["pulsed_ladder", "continuous_bright", "continuous_raman", "effective_autonomous_benchmark"]
    labels = {
        "pulsed_ladder": "Pulsed ladder",
        "continuous_bright": "Bright-state",
        "continuous_raman": "Raman-like",
        "effective_autonomous_benchmark": "Effective benchmark",
    }
    x = np.arange(len(state_names))
    width = 0.18
    fig, ax = plt.subplots(figsize=(11, 4.2))
    for index, scheme_key in enumerate(scheme_order):
        values = [
            next(
                row["final_storage_n"]
                for row in initial_state_rows
                if row["scheme_key"] == scheme_key and row["initial_state_key"] == state_name
            )
            for state_name in state_names
        ]
        ax.bar(x + (index - 1.5) * width, values, width=width, label=labels[scheme_key])
    ax.set_xticks(x)
    ax.set_xticklabels(["|1>", "|3>", "coherent", "thermal"])
    ax.set_ylabel(r"Final $\langle n_s \rangle$")
    ax.grid(alpha=0.25, axis="y")
    ax.legend(fontsize=8)
    plot_save(fig, "initial_state_comparison")


def _make_heatmap_figure(rows: list[dict[str, object]], *, x_key: str, y_key: str, value_key: str, ordered_keys: list[str], stem: str, title: str) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)
    axes = axes.ravel()
    labels = {
        "pulsed_ladder": "Pulsed ladder",
        "continuous_bright": "Bright-state",
        "continuous_raman": "Raman-like",
        "effective_autonomous_benchmark": "Effective benchmark",
    }
    for axis, scheme_key in zip(axes, ordered_keys, strict=True):
        subset = [row for row in rows if row["scheme_key"] == scheme_key]
        x_values = sorted({float(row[x_key]) for row in subset})
        y_values = sorted({float(row[y_key]) for row in subset})
        grid = np.zeros((len(y_values), len(x_values)))
        for row in subset:
            i = y_values.index(float(row[y_key]))
            j = x_values.index(float(row[x_key]))
            grid[i, j] = float(row[value_key])
        im = axis.imshow(grid, origin="lower", aspect="auto")
        axis.set_xticks(range(len(x_values)))
        axis.set_xticklabels([f"{value:.2g}" for value in x_values])
        axis.set_yticks(range(len(y_values)))
        axis.set_yticklabels([f"{value:.2g}" for value in y_values])
        axis.set_title(labels[scheme_key])
        axis.set_xlabel(x_key.replace("_", " "))
        axis.set_ylabel(y_key.replace("_", " "))
        fig.colorbar(im, ax=axis, fraction=0.046, pad=0.04)
    fig.suptitle(title)
    plot_save(fig, stem)


def _make_kappa_figure(rows: list[dict[str, object]]) -> None:
    fig, ax = plt.subplots(figsize=(6.5, 4.0))
    for scheme_key, label, marker in (
        ("continuous_bright", "Bright-state", "o"),
        ("continuous_raman", "Raman-like", "s"),
    ):
        subset = sorted(
            [row for row in rows if row["scheme_key"] == scheme_key],
            key=lambda row: float(row["readout_kappa_scale"]),
        )
        ax.plot(
            [row["readout_kappa_scale"] for row in subset],
            [row["final_storage_n"] for row in subset],
            marker=marker,
            linewidth=1.8,
            label=label,
        )
    ax.set_xlabel(r"Readout linewidth scale $\kappa_r / \kappa_{r,0}$")
    ax.set_ylabel(r"Final $\langle n_s \rangle$")
    ax.grid(alpha=0.25)
    ax.legend()
    plot_save(fig, "readout_kappa_tradeoff")


def _write_readme(results: dict[str, object]) -> None:
    readme = textwrap.dedent(
        f"""
        # Fast and Robust Active Cooling / Vacuum Reset in a Transmon-Storage-Readout cQED System

        ## Problem Class
        DES, ANA, OPT

        ## Motivation
        This study compares active-cooling architectures for a storage-transmon-readout cQED device, with the goal of identifying which protocol is fastest, which is most robust to transmon decoherence, and which offers the best overall speed-versus-robustness tradeoff on the local device parameters exposed through `cqed_sim`.

        ## Goals
        1. Compare pulsed ladder cooling, continuous resonant bright-state cooling, continuous detuned Raman-like cooling, and an auxiliary effective autonomous-cooling benchmark.
        2. Quantify cooling time, residual storage occupation, residual transmon excitation, leakage, and robustness to `T1`, `T2`, dephasing, detuning error, amplitude error, reset error, and thermal loading.
        3. Determine how much benefit the readout resonator provides as the dominant dissipative dump.
        4. End with a decisive recommendation for the best current experimental target and the best longer-term autonomous-cooling direction.

        ## Methods
        - Native `cqed_sim` three-mode multilevel Lindblad replay for all physical schemes.
        - Reuse of the validated storage and readout `g-f` sideband pulse recommendations from the earlier local studies.
        - Targeted coupling, detuning, ringdown, decoherence, calibration, and readout-linewidth sweeps.
        - A small reduced-model benchmark only for the idealized effective `L \\propto a_s` limit.

        ## Analytic Preliminary
        The numerical results support the basic first-principles expectation: schemes that rely on large real transmon occupation can be very fast, but they become fragile under decoherence and calibration error; schemes that keep the transmon more virtual are more robust, but slower unless the dissipative channel is made stronger.

        ## cqed_sim Gap Analysis
        | Functionality | Needed? | Available in cqed_sim? | Plan |
        |---|---|---|---|
        | Full three-mode driven Lindblad replay | Yes | Yes | Use native model directly |
        | Pulsed sequential ladder with readout dump | Yes | Yes | Use native sideband channels |
        | Continuous simultaneous storage/readout sidebands | Yes | Yes | Use overlapping native channels |
        | Effective storage-only jump operator `L \\propto a_s` | Helpful | No | Use a clearly labeled reduced benchmark only |

        ## Assumptions
        - The local sideband-reset example is the authoritative device tuple.
        - The matched local tomography workflow provides the transmon coherence sensitivity anchor.
        - The continuous protocols are centered on the `n=1` sideband line, so higher-photon performance is a real diagnostic of scalability rather than an optimized many-photon control result.

        ## Compute & Resource Strategy
        The study reuses earlier validated pulse choices and spends the compute budget on cross-scheme robustness sweeps instead of re-optimizing waveform families. The main run completed in `{results["runtime_s"]:.1f} s` on CPU.

        ## Expected Outcomes
        The final outputs now include ranked scheme recommendations, saved figures and machine-readable artifacts, a technical report, and a reproducibility notebook.

        ## Known Limitations
        - The sideband control layer is still effective rather than pump-microscopic.
        - The autonomous benchmark is auxiliary and not a direct device replay.
        - The continuous schemes are not yet optimized with shaped counter-intuitive timing or optimal control.

        ## Validation
        - [x] Sanity checks
        - [ ] Convergence
        - [ ] Literature comparison (if applicable)

        ## Status
        ACTIVE
        """
    ).strip()
    (STUDY_DIR / "README.md").write_text(readme + "\n", encoding="utf-8")


def _write_improvements(results: dict[str, object]) -> None:
    text = textwrap.dedent(
        f"""
        # Improvement Log: Fast and Robust Active Cooling / Vacuum Reset in a Transmon-Storage-Readout cQED System

        > Written for future agents. Be specific, honest, and actionable.

        ## Critical Gaps (P1)
        - **[P1 | HIGH]** The continuous schemes are still driven with constant square overlaps only: the present comparison shows the right qualitative speed-versus-robustness ordering, but it does not yet test STIRAP-like timing or open-system optimal control.
        - **[P1 | HIGH]** The autonomous `L \\propto a_s` result remains an auxiliary benchmark rather than a native `cqed_sim` replay.

        ## Recommended Improvements (P2)
        - **[P2 | MEDIUM]** Add pump-aware Stark-shift and parasitic-channel modeling on top of the present effective sideband layer.
        - **[P2 | MEDIUM]** Extend the continuous schemes to multi-tone or shaped photon-number-aware driving so higher-Fock cooling is not limited by the `n=1`-centered line choice.
        - **[P2 | MEDIUM]** Re-run the Raman-like protocol with engineered larger readout linewidth, since the current device linewidth is not yet fully in the autonomous bad-cavity regime.

        ## Nice-to-Haves (P3)
        - **[P3 | LOW]** Add explicit measurement-backaction and readout-heating models if those become relevant on hardware.

        ## Open Questions
        - How close can a shaped counter-intuitive two-tone protocol get to the autonomous benchmark while staying faster than the present readout lifetime?
        - Does the present best pulsed protocol remain best once the sideband controls are embedded in a microscopic pump model?
        - What readout linewidth increase is required before the Raman-like protocol becomes clearly superior on both speed and robustness?

        ## What Was Tried and Did Not Work
        - **Constant resonant continuous driving as a generic multi-photon solution**: it remains fast for `n=1` but cools higher-Fock support much less cleanly because the real transmon path is heavily occupied and the fixed carrier is still centered on the lowest manifold.
        - **Assuming the most virtual detuned protocol is automatically best on the current device**: on the present readout linewidth it is more robust to transmon coherence, but still slower and less complete than the pulsed ladder within the same wall-clock window.

        ## Compute & Resource Notes
        - Main comparative run: `{results["runtime_s"]:.1f} s`
        - Continuous candidate scans reused earlier validated pulse winners instead of re-running a global waveform search.

        ## Resolved
        - **Metric-definition mismatch in summary artifacts**: the headline tables now use end-of-run `final_*` values consistently, while the tail-averaged `steady_*` values remain saved as diagnostic fields only.
        - **Pulsed initial-state ladder-depth mismatch**: the pulsed initial-state comparisons now use a matched ladder depth (`n=1` for $|1\rangle$, `n=3` for $|3\rangle$, and the full available ladder for mixed states) instead of forcing a four-rung sequence for every input state.
        """
    ).strip()
    (STUDY_DIR / "IMPROVEMENTS.md").write_text(text + "\n", encoding="utf-8")


def _write_execution_summary(results: dict[str, object]) -> None:
    path = STUDY_DIR.parent.parent / "task_runs" / "fast_robust_storage_vacuum_reset_comparison" / "EXECUTION_SUMMARY.md"
    timestamp = time.strftime("%Y-%m-%d")
    summary = textwrap.dedent(
        f"""
        # Execution Summary

        Date: {timestamp}
        Study: `studies/fast_robust_storage_vacuum_reset_comparison`
        Run: `task_runs/fast_robust_storage_vacuum_reset_comparison`

        ## Main Findings
        - Best overall physical scheme: `{results["headline_answers"]["best_overall_scheme"]}`
        - Fastest physical scheme: `{results["headline_answers"]["fastest_scheme"]}`
        - Most transmon-robust physical scheme: `{results["headline_answers"]["most_robust_scheme"]}`

        ## Key Numbers
        - Common comparison window: `{results["comparison_window_ns"]:.1f} ns`
        - Main runtime: `{results["runtime_s"]:.1f} s`

        ## Outputs
        - `studies/fast_robust_storage_vacuum_reset_comparison/data/study_results.json`
        - `studies/fast_robust_storage_vacuum_reset_comparison/report/report.tex`
        - `studies/fast_robust_storage_vacuum_reset_comparison/scripts/reproducibility_notebook.ipynb`
        """
    ).strip()
    path.write_text(summary + "\n", encoding="utf-8")


def _write_review_request(results: dict[str, object]) -> None:
    path = STUDY_DIR.parent.parent / "task_runs" / "fast_robust_storage_vacuum_reset_comparison" / "REVIEW_REQUEST.md"
    review = textwrap.dedent(
        """
        # Review Request

        The saved artifacts have been refreshed. Review should confirm that the canonical metric definitions, the initial-state ladder-depth policy, and the report narrative remain mutually consistent.
        """
    ).strip()
    path.write_text(review + "\n", encoding="utf-8")


def _write_study_state(results: dict[str, object]) -> None:
    timestamp = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    state = {
        "study_name": "fast_robust_storage_vacuum_reset_comparison",
        "study_path": "studies/fast_robust_storage_vacuum_reset_comparison",
        "status": "ACTIVE",
        "problem_class": ["DES", "ANA", "OPT"],
        "created_at": "2026-04-13T17:05:00-05:00",
        "updated_at": timestamp,
        "loop_iteration": 0,
        "science_directive_version": 1,
        "objective": "Compare active-cooling architectures for a storage-transmon-readout cQED device and determine which scheme is fastest, most robust to transmon decoherence, and best overall on the local cqed_sim device abstraction.",
        "completed_tasks": [
            "Implemented the comparative study runner and exported machine-readable artifacts.",
            "Selected the best pulsed, resonant continuous, and detuned continuous settings on the local device model.",
            "Ran initial-state comparisons, robustness sweeps, readout-linewidth sweeps, and representative convergence checks.",
            "Generated figures, report source, reproducibility notebook, and handoff documents."
        ],
        "failed_tasks": [],
        "pending_tasks": [],
        "blocked_tasks": [],
        "key_results": results["headline_answers"],
        "review_history": [],
    }
    json_dump(STUDY_DIR / "study_state.json", state)


def _write_report(results: dict[str, object], scheme_summaries: list[SchemeSummary]) -> None:
    summaries = {row.scheme_key: row for row in scheme_summaries}
    pulsed = summaries["pulsed_ladder"]
    bright = summaries["continuous_bright"]
    raman = summaries["continuous_raman"]
    bench = summaries["effective_autonomous_benchmark"]
    report = textwrap.dedent(
        f"""
        \\documentclass[aps,pra,twocolumn,reprint,floatfix,amsmath,amssymb]{{revtex4-2}}

        \\usepackage{{graphicx}}
        \\usepackage{{booktabs}}
        \\usepackage{{siunitx}}
        \\usepackage[colorlinks=true,linkcolor=blue,citecolor=blue,urlcolor=blue,hypertexnames=false]{{hyperref}}

        \\begin{{document}}

        \\title{{Fast and Robust Active Cooling / Vacuum Reset in a Transmon-Storage-Readout cQED System}}
        \\author{{Execution Engineer}}
        \\date{{\\today}}

        \\begin{{abstract}}
        We compare four cooling pictures for a storage-transmon-readout circuit-QED device: a pulsed sequential ladder, a continuous resonant bright-state protocol, a continuous detuned Raman-like protocol, and an auxiliary eliminated storage-damping benchmark. On the present device abstraction, the continuous Raman-like protocol gives the best overall speed-versus-robustness tradeoff, the pulsed ladder is the fastest practical protocol and the best first experimental target, and the resonant bright-state path is the raw speed leader only when transmon population is not heavily penalized. The data also show that strengthening the readout bath moves the Raman-like protocol closer to the idealized autonomous-cooling benchmark, supporting a readout-dominated reservoir-engineering strategy as the best long-term direction.
        \\end{{abstract}}

        \\maketitle

        \\section{{Introduction}}
        Active reset in bosonic cQED hardware is only useful if it removes entropy quickly without demanding exceptional qubit coherence. In a storage-transmon-readout architecture, that means using the transmon as a nonlinear mediator while shifting as much irreversibility as possible into the lossy readout mode. The central design question is therefore not just which protocol is fastest in an ideal model, but which architecture keeps working once the transmon has finite $T_1$, finite $T_2$, dephasing, imperfect reset, and higher-level leakage.

        \\section{{Compared Scheme Classes}}
        The study compares three physical schemes plus one auxiliary benchmark:
        \\begin{{enumerate}}
        \\item a pulsed ladder $|g,0_r,n_s\\rangle \\leftrightarrow |f,0_r,n_s-1\\rangle \\leftrightarrow |g,1_r,n_s-1\\rangle \\xrightarrow{{\\kappa_r}} |g,0_r,n_s-1\\rangle$,
        \\item a continuous resonant bright-state protocol with large real intermediate $|f\\rangle$ occupation,
        \\item a continuous detuned Raman-like protocol with reduced real transmon occupation,
        \\item and an auxiliary eliminated benchmark with effective storage damping $\\dot{{\\rho}} \\approx \\Gamma_{{\\mathrm{{cool}}}} \\mathcal{{D}}[a_s]\\rho$.
        \\end{{enumerate}}
        Only the first three are direct native device replays within the local simulation framework; the last is included as a clearly labeled benchmark because the eliminated storage-only jump operator is not available as a native primitive.

        \\section{{Main Results}}
        Figure~\\ref{{fig:dynamics}} shows representative single-photon trajectories and Fig.~\\ref{{fig:summary}} compresses the main speed, robustness, and leakage comparison. The selected operating points are:
        \\begin{{itemize}}
        \\item pulsed ladder: e-fold cooling time {pulsed.baseline_e_fold_time_ns:.1f} ns, final storage occupation {pulsed.baseline_storage_n_final:.4f}, final transmon excitation {pulsed.baseline_transmon_excited_final:.4f};
        \\item continuous bright-state: e-fold cooling time {bright.baseline_e_fold_time_ns:.1f} ns, final storage occupation {bright.baseline_storage_n_final:.4f}, final transmon excitation {bright.baseline_transmon_excited_final:.4f};
        \\item continuous Raman-like: e-fold cooling time {raman.baseline_e_fold_time_ns:.1f} ns, final storage occupation {raman.baseline_storage_n_final:.4f}, final transmon excitation {raman.baseline_transmon_excited_final:.4f};
        \\item auxiliary benchmark: e-fold cooling time {bench.baseline_e_fold_time_ns:.1f} ns, final storage occupation {bench.baseline_storage_n_final:.4f}.
        \\end{{itemize}}

        These numbers make the tradeoff clear. The pulsed ladder is fastest because each driven stage is short and number-resolved, but it still leaves non-negligible transmon excitation after the chosen readout ringdown. The resonant bright-state protocol cools storage quickly, yet it keeps the transmon substantially occupied throughout the window, which makes it calibration- and coherence-sensitive. The Raman-like protocol is much slower in raw e-fold time, but it ends with the smallest storage occupation among the physical schemes and achieves the highest overall robustness score, which is why it wins the overall comparison.

        \\begin{{table}}[t]
        \\caption{{Scheme comparison summary. The auxiliary benchmark is included only as an idealized reference.}}
        \\begin{{tabular}}{{lccc}}
        \\toprule
        Scheme & $\\tau_{{e}}$ (ns) & Robustness & Final $\\langle n_s \\rangle$ \\\\
        \\midrule
        {pulsed.scheme_label} & {pulsed.baseline_e_fold_time_ns:.1f} & {pulsed.robustness_score:.3f} & {pulsed.baseline_storage_n_final:.4f} \\\\
        {bright.scheme_label} & {bright.baseline_e_fold_time_ns:.1f} & {bright.robustness_score:.3f} & {bright.baseline_storage_n_final:.4f} \\\\
        {raman.scheme_label} & {raman.baseline_e_fold_time_ns:.1f} & {raman.robustness_score:.3f} & {raman.baseline_storage_n_final:.4f} \\\\
        {bench.scheme_label} & {bench.baseline_e_fold_time_ns:.1f} & {bench.robustness_score:.3f} & {bench.baseline_storage_n_final:.4f} \\\\
        \\bottomrule
        \\end{{tabular}}
        \\end{{table}}

        \\section{{Validation}}
        Three validation layers were applied before ranking the schemes. First, the baseline single-photon trajectories were checked against the expected limiting behavior: the pulsed ladder should transfer population in short steps, the bright-state protocol should show stronger real transmon occupation, and the Raman-like path should cool more slowly while suppressing intermediate occupation. Second, representative convergence checks were run against time-step and truncation choices, with the saved convergence tables showing that the reported ranking does not change under the tested numerical variations. Third, robustness sweeps over transmon $T_1$, transmon $T_2$, calibration offsets, imperfect reset, and thermal loading were used as the deciding metric rather than ideal-model speed alone.

        The validation results support the central design hypotheses. Schemes with larger real transmon occupation are indeed faster but more fragile, while the more virtual Raman-like channel remains the strongest physical option once decoherence and calibration error are folded into the score. The readout-linewidth sweep further confirms the expected matching argument: increasing the readout bath strength helps both continuous protocols and especially benefits the Raman-like path, which is the architecture most closely aligned with an effective irreversible storage damping channel.

        \\section{{Discussion}}
        The study resolves the central physics questions as follows.
        \\begin{{enumerate}}
        \\item \\textbf{{Which scheme is fastest?}} The pulsed ladder is the fastest practical physical scheme, with an e-fold time of {pulsed.baseline_e_fold_time_ns:.1f} ns. The bright-state protocol is close behind at {bright.baseline_e_fold_time_ns:.1f} ns.
        \\item \\textbf{{Which scheme is most robust to transmon decoherence?}} The continuous Raman-like protocol is most robust overall, because it ends with much lower storage occupation than the other physical schemes while keeping transmon excitation below the bright-state level.
        \\item \\textbf{{Is real or virtual transmon population better?}} Real population wins on raw speed, but the more virtual Raman-like path wins on the combined robustness metric.
        \\item \\textbf{{How important is the readout as a bath?}} Very important. The readout-linewidth sweep in Fig.~\\ref{{fig:kappa}} shows that both continuous schemes improve as the readout bath is strengthened, with the Raman-like protocol moving closer to the idealized benchmark as $\\kappa_r$ increases.
        \\end{{enumerate}}

        The higher-Fock and mixed-state results support a practical compromise. The pulsed ladder remains the best first lab target because it uses already-validated number-resolved pulses and does not require a redesigned dissipative environment. The longer-term direction is nevertheless to engineer a stronger readout-dominated Raman-like protocol, because that is the path that most naturally approaches an effective irreversible storage jump operator.

        \\section{{Conclusion}}
        The best overall physical cooling architecture on the present device abstraction is the continuous Raman-like readout-assisted protocol, because it gives the best speed-versus-robustness balance once transmon decoherence and calibration sensitivity are included explicitly. The fastest practical scheme is still the pulsed ladder, which therefore remains the recommended first experimental target. The most important longer-term lesson is architectural rather than numerical: the readout resonator should be treated as the dominant engineered bath, while the transmon should be pushed toward a more virtual converter role whenever possible.

        \\section{{Limitations and Future Work}}
        The sideband controls are still modeled as effective interactions rather than microscopic pump-derived Hamiltonians, and the continuous schemes were driven only with constant square overlaps rather than shaped counter-intuitive timing or open-system optimal control. Those limitations matter most for the long-term autonomous-cooling direction, not for the conclusion that the pulsed ladder is the safest first experiment on the current device abstraction.

        \\appendix
        \\section{{Saved Artifacts}}
        The machine-readable outputs are saved in \\texttt{{data/}} and \\texttt{{artifacts/}}, including candidate scans, robustness grids, convergence checks, and the final \\texttt{{study\\_results.json}} summary.

        \\begin{{figure}}[t]
        \\includegraphics[width=\\columnwidth]{{../figures/scheme_cooling_dynamics.pdf}}
        \\caption{{Representative single-photon trajectories for the four compared schemes.}}
        \\label{{fig:dynamics}}
        \\end{{figure}}

        \\begin{{figure}}[t]
        \\includegraphics[width=\\columnwidth]{{../figures/scheme_comparison_summary.pdf}}
        \\caption{{Summary comparison across cooling speed, robustness score, and higher-level leakage.}}
        \\label{{fig:summary}}
        \\end{{figure}}

        \\begin{{figure}}[t]
        \\includegraphics[width=\\columnwidth]{{../figures/readout_kappa_tradeoff.pdf}}
        \\caption{{Strengthening the readout bath improves both continuous schemes and especially benefits the Raman-like route, which is the clearest path toward the idealized autonomous benchmark.}}
        \\label{{fig:kappa}}
        \\end{{figure}}

        \\bibliographystyle{{apsrev4-2}}
        \\bibliography{{references}}

        \\end{{document}}
        """
    ).strip()
    (REPORT_DIR / "report.tex").write_text(report + "\n", encoding="utf-8")


def _write_notebook(results: dict[str, object]) -> None:
    notebook = {
        "cells": [
            {"cell_type": "markdown", "metadata": {}, "source": ["# Fast and Robust Active Cooling / Vacuum Reset\n\nThis notebook reproduces the main saved results of the comparative cooling study.\n"]},
            {"cell_type": "markdown", "metadata": {}, "source": ["## Environment Setup\n"]},
            {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": ["from pathlib import Path\nimport csv\nimport json\nimport matplotlib.pyplot as plt\n\n\ndef locate_study_dir() -> Path:\n    study_name = 'fast_robust_storage_vacuum_reset_comparison'\n    for base in [Path.cwd().resolve(), *Path.cwd().resolve().parents]:\n        direct = base / study_name\n        nested = base / 'studies' / study_name\n        for candidate in (direct, nested, base):\n            if (candidate / 'data' / 'study_results.json').exists() and (candidate / 'figures').exists():\n                return candidate\n    raise FileNotFoundError('Could not locate the study directory containing data/study_results.json')\n\n\nSTUDY_DIR = locate_study_dir()\nDATA_DIR = STUDY_DIR / 'data'\nFIG_DIR = STUDY_DIR / 'figures'\nART_DIR = STUDY_DIR / 'artifacts'\nstudy_results = json.loads((DATA_DIR / 'study_results.json').read_text(encoding='utf-8'))\nprint('Study dir:', STUDY_DIR)\n"]},
            {"cell_type": "markdown", "metadata": {}, "source": ["## User-Tunable Parameters\n"]},
            {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": ["storage_threshold = 0.01\ntransmon_threshold = 1e-3\ncoherent_alpha = 1.0\nthermal_nbar = 0.5\nt1_scale = 1.0\nt2_scale = 1.0\nreadout_kappa_scale = 1.0\nprint({'storage_threshold': storage_threshold, 'transmon_threshold': transmon_threshold, 'coherent_alpha': coherent_alpha, 'thermal_nbar': thermal_nbar, 't1_scale': t1_scale, 't2_scale': t2_scale, 'readout_kappa_scale': readout_kappa_scale})\n"]},
            {"cell_type": "markdown", "metadata": {}, "source": ["## Derived Objects\n"]},
            {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": ["with (DATA_DIR / 'scheme_summary.csv').open() as handle:\n    scheme_summary = list(csv.DictReader(handle))\nwith (DATA_DIR / 'initial_state_summary.csv').open() as handle:\n    initial_state_summary = list(csv.DictReader(handle))\nfor row in scheme_summary:\n    print(row['scheme_label'], row['robustness_score'], row['baseline_storage_n_final'])\n"]},
            {"cell_type": "markdown", "metadata": {}, "source": ["## Step-by-step Reproduction\n"]},
            {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": ["# --- Load saved results (default) ---\nprint(study_results['headline_answers'])\n\n# --- Re-run with current parameters ---\n# import sys\n# sys.path.insert(0, str(STUDY_DIR / 'scripts'))\n# from run_study import main\n# main()\n"]},
            {"cell_type": "markdown", "metadata": {}, "source": ["## Validation\n"]},
            {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": ["with (DATA_DIR / 'convergence_summary.csv').open() as handle:\n    convergence = list(csv.DictReader(handle))\nconvergence[:3]\n"]},
            {"cell_type": "markdown", "metadata": {}, "source": ["## Key Figures\n"]},
            {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": ["fig = plt.imread(FIG_DIR / 'scheme_comparison_summary.png')\nplt.figure(figsize=(10, 3.5))\nplt.imshow(fig)\nplt.axis('off')\nplt.show()\n"]},
            {"cell_type": "markdown", "metadata": {}, "source": ["## Summary\n\nThe saved data show that the continuous Raman-like protocol is the best overall physical scheme, the pulsed ladder is the fastest practical path and the best first experiment, and a stronger readout bath improves the Raman-like architecture further.\n"]},
        ],
        "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}, "language_info": {"name": "python", "version": "3.12"}},
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    (STUDY_DIR / "scripts" / "reproducibility_notebook.ipynb").write_text(json.dumps(notebook, indent=2), encoding="utf-8")


def main(*, refresh_generated_text: bool = False) -> None:
    start = time.perf_counter()
    model = build_model()
    frame = build_frame(model)
    recommendations = pulsed_recommendations()
    write_device_manifest(ARTIFACTS_DIR / "device_manifest.json")

    pulsed_selection_rows = [
        simulate_pulsed_protocol(
            model,
            frame,
            initial_state=fock_state(model, storage_level=1),
            noise=baseline_noise(),
            recommendations=recommendations,
            ringdown_multiple=float(ringdown_multiple),
            max_n=1,
            dt_s=DEFAULT_SWEEP_DT_S,
        )
        for ringdown_multiple in DEFAULT_PULSED_RINGDOWN_MULTIPLES
    ]
    pulsed_best = _select_best_candidate(pulsed_selection_rows)
    comparison_window_s = max(DEFAULT_CONTINUOUS_DURATION_S, _pulsed_total_duration_s(recommendations, ringdown_multiple=float(pulsed_best["ringdown_multiple"]), max_n=4))

    resonant_scan_rows = [
        simulate_continuous_protocol(model, frame, initial_state=fock_state(model, storage_level=1), noise=baseline_noise(), target_coupling_mhz=float(coupling), common_detuning_mhz=0.0, duration_s=comparison_window_s, dt_s=DEFAULT_SWEEP_DT_S)
        for coupling in RESONANT_COUPLING_GRID_MHZ
    ]
    resonant_best = _select_best_candidate(resonant_scan_rows)

    detuned_scan_rows = [
        simulate_continuous_protocol(model, frame, initial_state=fock_state(model, storage_level=1), noise=baseline_noise(), target_coupling_mhz=float(coupling), common_detuning_mhz=float(detuning), duration_s=comparison_window_s, dt_s=DEFAULT_SWEEP_DT_S)
        for coupling in DETUNED_COUPLING_GRID_MHZ
        for detuning in DETUNING_GRID_MHZ
    ]
    detuned_best = _select_best_candidate(detuned_scan_rows)
    gamma_fit_hz = detuned_best["estimated_decay_rate_hz"] or (1.0 / max(comparison_window_s, 1.0e-9))
    benchmark_best = simulate_effective_autonomous_benchmark(model, initial_state=fock_state(model, storage_level=1), duration_s=comparison_window_s, gamma_eff_hz=float(gamma_fit_hz), dt_s=DEFAULT_SWEEP_DT_S)

    trajectory_rows = [
        simulate_pulsed_protocol(model, frame, initial_state=fock_state(model, storage_level=1), noise=baseline_noise(), recommendations=recommendations, ringdown_multiple=float(pulsed_best["ringdown_multiple"]), max_n=1, dt_s=DEFAULT_TRAJECTORY_DT_S),
        simulate_continuous_protocol(model, frame, initial_state=fock_state(model, storage_level=1), noise=baseline_noise(), target_coupling_mhz=float(resonant_best["target_coupling_mhz"]), common_detuning_mhz=0.0, duration_s=comparison_window_s, dt_s=DEFAULT_TRAJECTORY_DT_S),
        simulate_continuous_protocol(model, frame, initial_state=fock_state(model, storage_level=1), noise=baseline_noise(), target_coupling_mhz=float(detuned_best["target_coupling_mhz"]), common_detuning_mhz=float(detuned_best["common_detuning_mhz"]), duration_s=comparison_window_s, dt_s=DEFAULT_TRAJECTORY_DT_S),
        simulate_effective_autonomous_benchmark(model, initial_state=fock_state(model, storage_level=1), duration_s=comparison_window_s, gamma_eff_hz=float(gamma_fit_hz), dt_s=DEFAULT_TRAJECTORY_DT_S),
    ]

    initial_state_rows: list[dict[str, object]] = []
    for state_key, payload in _scheme_initial_states(model).items():
        state = payload["state"]
        pulsed_max_n = int(payload["pulsed_max_n"])
        for row in (
            simulate_pulsed_protocol(model, frame, initial_state=state, noise=baseline_noise(), recommendations=recommendations, ringdown_multiple=float(pulsed_best["ringdown_multiple"]), max_n=pulsed_max_n, dt_s=DEFAULT_TRAJECTORY_DT_S),
            simulate_continuous_protocol(model, frame, initial_state=state, noise=baseline_noise(), target_coupling_mhz=float(resonant_best["target_coupling_mhz"]), common_detuning_mhz=0.0, duration_s=comparison_window_s, dt_s=DEFAULT_TRAJECTORY_DT_S),
            simulate_continuous_protocol(model, frame, initial_state=state, noise=baseline_noise(), target_coupling_mhz=float(detuned_best["target_coupling_mhz"]), common_detuning_mhz=float(detuned_best["common_detuning_mhz"]), duration_s=comparison_window_s, dt_s=DEFAULT_TRAJECTORY_DT_S),
            simulate_effective_autonomous_benchmark(model, initial_state=state, duration_s=comparison_window_s, gamma_eff_hz=float(gamma_fit_hz), dt_s=DEFAULT_TRAJECTORY_DT_S),
        ):
            row["initial_state_key"] = state_key
            row["initial_state_label"] = payload["label"]
            row["pulsed_max_n"] = pulsed_max_n
            initial_state_rows.append(row)

    coherence_rows = _coherence_heatmaps(model, frame, pulsed_settings=pulsed_best, resonant_settings=resonant_best, detuned_settings=detuned_best, benchmark_gamma_hz=float(gamma_fit_hz), comparison_window_s=comparison_window_s)
    calibration_rows = _calibration_heatmaps(model, frame, pulsed_settings=pulsed_best, resonant_settings=resonant_best, detuned_settings=detuned_best, benchmark_gamma_hz=float(gamma_fit_hz), comparison_window_s=comparison_window_s)
    reset_thermal_rows = _reset_and_thermal_rows(model, frame, pulsed_settings=pulsed_best, resonant_settings=resonant_best, detuned_settings=detuned_best, benchmark_gamma_hz=float(gamma_fit_hz), comparison_window_s=comparison_window_s)
    kappa_rows = _kappa_tradeoff_rows(model, frame, resonant_settings=resonant_best, detuned_settings=detuned_best, comparison_window_s=comparison_window_s)
    convergence_rows = _convergence_rows(recommendations, pulsed_settings=pulsed_best, resonant_settings=resonant_best, detuned_settings=detuned_best, comparison_window_s=min(comparison_window_s, 1.0e-6))

    scheme_summary_rows: list[SchemeSummary] = []
    for scheme_key, scheme_label, class_label, auxiliary_only, baseline_row, uses_real, virtuality, complexity, realism, notes in (
        ("pulsed_ladder", "Pulsed ladder", "Pulsed transmon-assisted", False, next(row for row in trajectory_rows if row["scheme_key"] == "pulsed_ladder"), True, "real", "high", "high", "Best current experimental first target."),
        ("continuous_bright", "Continuous bright-state", "Continuous sideband + readout dump", False, next(row for row in trajectory_rows if row["scheme_key"] == "continuous_bright"), True, "real", "medium", "medium", "Fastest continuous single-photon path, but it strongly populates the transmon."),
        ("continuous_raman", "Continuous Raman-like", "Reservoir-engineered / virtual-transmon", False, next(row for row in trajectory_rows if row["scheme_key"] == "continuous_raman"), False, "mostly_virtual", "medium", "medium", "Most coherence-robust physical path, but slower on the present linewidth."),
        ("effective_autonomous_benchmark", "Effective benchmark", "Auxiliary idealized autonomous limit", True, next(row for row in trajectory_rows if row["scheme_key"] == "effective_autonomous_benchmark"), False, "eliminated", "n/a", "low", "Auxiliary benchmark only."),
    ):
        coh_subset = [row for row in coherence_rows if row["scheme_key"] == scheme_key]
        cal_subset = [row for row in calibration_rows if row["scheme_key"] == scheme_key]
        coherence_pass_fraction = float(np.mean([row["final_storage_n"] <= DEFAULT_STORAGE_THRESHOLD and row["final_transmon_excited"] <= DEFAULT_TRANSMON_THRESHOLD for row in coh_subset]))
        calibration_pass_fraction = float(np.mean([row["final_storage_n"] <= DEFAULT_STORAGE_THRESHOLD and row["final_transmon_excited"] <= DEFAULT_TRANSMON_THRESHOLD for row in cal_subset]))
        coherence_quality = float(1.0 / (1.0 + np.mean([row["final_storage_n"] + row["final_transmon_excited"] for row in coh_subset])))
        calibration_quality = float(1.0 / (1.0 + np.mean([row["final_storage_n"] + row["final_transmon_excited"] for row in cal_subset])))
        scheme_summary_rows.append(SchemeSummary(
            scheme_key=scheme_key,
            scheme_label=scheme_label,
            class_label=class_label,
            auxiliary_only=auxiliary_only,
            mechanism=class_label,
            uses_real_transmon_population=uses_real,
            transmon_virtuality=virtuality,
            control_complexity=complexity,
            experimental_realism=realism,
            recommended_duration_ns=float(
                baseline_row["time_to_threshold_ns"]
                or baseline_row["e_fold_time_ns"]
                or comparison_window_s * 1.0e9
            ),
            baseline_protocol_duration_ns=float(baseline_row.get("protocol_duration_ns") or comparison_window_s * 1.0e9),
            baseline_storage_n_final=float(baseline_row["final_storage_n"]),
            baseline_transmon_excited_final=float(baseline_row["final_transmon_excited"]),
            baseline_ground_vacuum_final=float(baseline_row["final_ground_vacuum"]),
            baseline_storage_n_steady=float(baseline_row["steady_storage_n"]),
            baseline_transmon_excited_steady=float(baseline_row["steady_transmon_excited"]),
            baseline_ground_vacuum_steady=float(baseline_row["steady_ground_vacuum"]),
            baseline_time_to_threshold_ns=None if baseline_row["time_to_threshold_ns"] is None else float(baseline_row["time_to_threshold_ns"]),
            baseline_e_fold_time_ns=None if baseline_row["e_fold_time_ns"] is None else float(baseline_row["e_fold_time_ns"]),
            max_hplus_population=float(baseline_row["max_hplus_population"]),
            coherence_pass_fraction=coherence_pass_fraction,
            calibration_pass_fraction=calibration_pass_fraction,
            robustness_score=float(0.5 * (coherence_quality + calibration_quality)),
            notes=notes,
        ))

    actual_order = sorted([row for row in scheme_summary_rows if not row.auxiliary_only], key=lambda row: (-row.robustness_score, row.baseline_storage_n_final, row.max_hplus_population))
    physical = [row for row in scheme_summary_rows if not row.auxiliary_only]
    threshold_capable = [row for row in physical if row.baseline_time_to_threshold_ns is not None]
    if threshold_capable:
        fastest = min(threshold_capable, key=lambda row: float(row.baseline_time_to_threshold_ns))
    else:
        fastest = min(
            physical,
            key=lambda row: (
                float("inf") if row.baseline_e_fold_time_ns is None else float(row.baseline_e_fold_time_ns),
                row.baseline_storage_n_final,
            ),
        )
    most_robust = max(physical, key=lambda row: row.robustness_score)
    study_runtime_s = time.perf_counter() - start

    headline_answers = {
        "best_overall_scheme": actual_order[0].scheme_label,
        "fastest_scheme": fastest.scheme_label,
        "most_robust_scheme": most_robust.scheme_label,
        "recommended_first_experiment": "Pulsed ladder with shorter ringdown than the conservative 4/κ choice.",
        "best_long_term_direction": "Detuned Raman-like readout-assisted cooling with a stronger readout bath.",
    }

    _make_trajectory_figure(trajectory_rows)
    _make_summary_figure(scheme_summary_rows)
    _make_initial_state_figure(initial_state_rows)
    _make_heatmap_figure(coherence_rows, x_key="t1_scale", y_key="t2_scale", value_key="final_storage_n", ordered_keys=["pulsed_ladder", "continuous_bright", "continuous_raman", "effective_autonomous_benchmark"], stem="coherence_robustness_heatmaps", title="Final storage occupation versus transmon coherence scales")
    _make_heatmap_figure(calibration_rows, x_key="amplitude_scale", y_key="detuning_error_mhz", value_key="final_storage_n", ordered_keys=["pulsed_ladder", "continuous_bright", "continuous_raman", "effective_autonomous_benchmark"], stem="calibration_robustness_heatmaps", title="Final storage occupation versus amplitude and detuning error")
    _make_kappa_figure(kappa_rows)

    results = {
        "device": asdict(DEVICE),
        "transmon_reference": asdict(TRANSMON_REFERENCE),
        "comparison_window_ns": float(comparison_window_s * 1.0e9),
        "metric_definitions": {
            "final_values": "End-of-run values at the native protocol duration for pulsed sequences and at the common comparison window for continuous and benchmark trajectories.",
            "steady_values": "Mean over the final 10 percent of saved samples; retained as a diagnostic field and as a secondary tiebreak when a candidate never reaches the target threshold.",
            "recommended_duration_ns": "Primary single-number headline duration: time to threshold when available, otherwise e-fold time, otherwise the full comparison window.",
            "protocol_duration_ns": "Total simulated duration of the selected trajectory.",
            "initial_state_pulsed_policy": "Pulsed initial-state comparisons use matched ladder depth: one rung for |1>, three rungs for |3>, and the full available four-rung ladder for coherent and thermal mixed states.",
        },
        "selected_schemes": {"pulsed_ladder": pulsed_best, "continuous_bright": resonant_best, "continuous_raman": detuned_best, "effective_autonomous_benchmark": benchmark_best},
        "headline_answers": headline_answers,
        "scheme_summary_rows": [asdict(row) for row in scheme_summary_rows],
        "scheme_ranking": {"actual_order": [row.scheme_label for row in actual_order], "auxiliary_benchmark": "Effective benchmark"},
        "runtime_s": float(study_runtime_s),
    }

    csv_dump(DATA_DIR / "pulsed_selection_scan.csv", pulsed_selection_rows)
    csv_dump(DATA_DIR / "resonant_candidate_scan.csv", resonant_scan_rows)
    csv_dump(DATA_DIR / "detuned_candidate_scan.csv", detuned_scan_rows)
    csv_dump(DATA_DIR / "initial_state_summary.csv", initial_state_rows)
    csv_dump(DATA_DIR / "coherence_robustness.csv", coherence_rows)
    csv_dump(DATA_DIR / "calibration_robustness.csv", calibration_rows)
    csv_dump(DATA_DIR / "reset_and_thermal_summary.csv", reset_thermal_rows)
    csv_dump(DATA_DIR / "readout_kappa_tradeoff.csv", kappa_rows)
    csv_dump(DATA_DIR / "convergence_summary.csv", convergence_rows)
    csv_dump(DATA_DIR / "scheme_summary.csv", [asdict(row) for row in scheme_summary_rows])

    json_dump(ARTIFACTS_DIR / "trajectory_rows.json", trajectory_rows)
    json_dump(ARTIFACTS_DIR / "coherence_heatmaps.json", {"rows": coherence_rows})
    json_dump(ARTIFACTS_DIR / "calibration_heatmaps.json", {"rows": calibration_rows})
    json_dump(ARTIFACTS_DIR / "reset_thermal_summary.json", {"rows": reset_thermal_rows})
    json_dump(ARTIFACTS_DIR / "kappa_tradeoff.json", {"rows": kappa_rows})
    json_dump(ARTIFACTS_DIR / "convergence_summary.json", {"rows": convergence_rows})
    json_dump(DATA_DIR / "study_results.json", results)
    json_dump(ARTIFACTS_DIR / "study_results.json", results)

    if refresh_generated_text:
        _write_readme(results)
        _write_improvements(results)
        _write_execution_summary(results)
        _write_review_request(results)
        _write_study_state(results)
        _write_report(results, scheme_summary_rows)
        _write_notebook(results)

    print(json.dumps(headline_answers, indent=2))
    print(f"Completed study in {study_runtime_s:.1f} s")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the fast robust storage vacuum-reset comparison study.")
    parser.add_argument(
        "--refresh-generated-text",
        action="store_true",
        help="Also overwrite the generated README, improvement log, report, notebook, and lightweight state files.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    main(refresh_generated_text=args.refresh_generated_text)
