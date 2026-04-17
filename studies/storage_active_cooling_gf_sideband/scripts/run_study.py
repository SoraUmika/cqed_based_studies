"""Run the storage gf-sideband active-cooling study end to end."""

from __future__ import annotations

from dataclasses import asdict
import json
import math
import time

import numpy as np
import qutip as qt

from cqed_sim import FloquetConfig, FloquetProblem
from cqed_sim.core.drive_targets import SidebandDriveSpec, TransmonTransitionDriveSpec
from cqed_sim.floquet import build_target_drive_term, solve_floquet

from common import (
    ARTIFACTS_DIR,
    DATA_DIR,
    DEFAULT_DT_S,
    DEVICE,
    build_frame,
    build_model,
    build_noise,
    basis_label,
    cavity_wigner,
    csv_dump,
    expected_pi_duration_s,
    hz,
    json_dump,
    ladder_sequence_times,
    make_pulse,
    plot_save,
    readout_dump_lab_frequency,
    readout_dump_rotating_frequency,
    reduced_storage_state,
    sideband_matrix_element,
    simulate_single_stage,
    sorted_basis_populations,
    state_population,
    storage_photon_number,
    readout_photon_number,
    storage_sideband_lab_frequency,
    storage_sideband_rotating_frequency,
    transmon_excited_population,
    write_device_manifest,
)

PULSE_FAMILY_INFO = {
    "square": {
        "envelope": "square",
        "drag": 0.0,
        "paper_role": "baseline constant-amplitude control",
    },
    "gaussian": {
        "envelope": "gaussian",
        "drag": 0.0,
        "paper_role": "baseline smooth finite-bandwidth control",
    },
    "gaussian_drag": {
        "envelope": "gaussian",
        "drag": 0.35,
        "paper_role": "DRAG-like derivative quadrature correction",
    },
    "cosine_squared": {
        "envelope": "cosine_squared",
        "drag": 0.0,
        "paper_role": "smooth ramp alternative to a square pulse",
    },
    "bump": {
        "envelope": "bump",
        "drag": 0.0,
        "paper_role": "paper-motivated adiabatic-to-Floquet ramp",
    },
    "phase_modulated_bump": {
        "envelope": "phase_modulated_bump",
        "drag": 0.0,
        "paper_role": "parametrically modulated sideband-style ramp",
    },
    "bb1_like": {
        "envelope": "bb1_like",
        "drag": 0.0,
        "paper_role": "composite-pulse robustness screen",
    },
}
PULSE_FAMILIES = tuple(PULSE_FAMILY_INFO.keys())
STORAGE_AMPLITUDES_MHZ = np.array([4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0])
STORAGE_DURATIONS_NS = np.array([10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 120.0, 160.0])
DUMP_AMPLITUDES_MHZ = np.linspace(4.0, 14.0, 6)
DUMP_DURATIONS_NS = np.array([20.0, 30.0, 40.0, 50.0, 60.0, 80.0])
DETUNING_OFFSETS_MHZ = np.linspace(-4.0, 4.0, 17)
CALIBRATION_DETUNING_OFFSETS_MHZ = np.linspace(-2.0, 2.0, 17)
AMPLITUDE_SCALE_GRID = np.linspace(0.9, 1.1, 11)
FLOQUET_TRUNCATION = {"n_tr": 6, "n_storage": 8, "n_readout": 4}


def nearest_unwanted_detuning_hz(values_hz: list[float], index: int) -> float:
    target = values_hz[index]
    other = [abs(target - value) for j, value in enumerate(values_hz) if j != index]
    return float(min(other)) if other else float("nan")


def make_family_pulse(
    *,
    channel: str,
    carrier_rad_s: float,
    duration_s: float,
    amplitude_hz: float,
    family: str,
    label: str | None = None,
):
    family_spec = PULSE_FAMILY_INFO[str(family)]
    return make_pulse(
        channel=channel,
        carrier_rad_s=carrier_rad_s,
        duration_s=duration_s,
        amplitude_hz=amplitude_hz,
        envelope_name=str(family_spec["envelope"]),
        drag=float(family_spec["drag"]),
        label=label,
    )


def spectrum_and_frequency_table(model, frame):
    storage_lab_values = [hz(storage_sideband_lab_frequency(model, n)) for n in range(1, 5)]
    readout_lab_values = [hz(readout_dump_lab_frequency(model, n)) for n in range(1, 5)]

    rows = []
    for n in range(1, 5):
        storage_lab_hz = hz(storage_sideband_lab_frequency(model, n))
        storage_rot_hz = hz(storage_sideband_rotating_frequency(model, frame, n))
        readout_lab_hz = hz(readout_dump_lab_frequency(model, n))
        readout_rot_hz = hz(readout_dump_rotating_frequency(model, frame, n))
        rows.append(
            {
                "n": n,
                "initial_state": basis_label(0, n, 0),
                "step_a_target_state": basis_label(2, n - 1, 0),
                "step_b_target_state": basis_label(0, n - 1, 1),
                "storage_sideband_lab_GHz": round(storage_lab_hz / 1.0e9, 9),
                "storage_sideband_rot_MHz": round(storage_rot_hz / 1.0e6, 6),
                "readout_dump_lab_GHz": round(readout_lab_hz / 1.0e9, 9),
                "readout_dump_rot_MHz": round(readout_rot_hz / 1.0e6, 6),
                "nearest_storage_sideband_detuning_MHz": round(
                    nearest_unwanted_detuning_hz(storage_lab_values, n - 1) / 1.0e6,
                    6,
                ),
                "nearest_readout_dump_detuning_MHz": round(
                    nearest_unwanted_detuning_hz(readout_lab_values, n - 1) / 1.0e6,
                    6,
                ),
                "gf_carrier_lab_GHz": round(
                    hz(
                        model.transmon_transition_frequency(
                            storage_level=n,
                            readout_level=0,
                            lower_level=0,
                            upper_level=2,
                        )
                    )
                    / 1.0e9,
                    9,
                ),
                "storage_sideband_matrix_element": round(sideband_matrix_element(model, mode="storage", n=n), 6),
                "readout_dump_matrix_element": round(sideband_matrix_element(model, mode="readout", n=n), 6),
                "storage_nominal_minus_model_kHz": round(
                    (DEVICE.storage_gf_sideband_nominal_hz - storage_lab_hz) / 1.0e3 if n == 1 else 0.0,
                    6,
                ),
                "spectral_crowding_comment": (
                    "Adjacent storage-sideband lines are only a few MHz apart; selectivity is set by the pulse bandwidth."
                ),
            }
        )

    h_static = model.static_hamiltonian(frame=frame)
    eigenvalues, eigenstates = h_static.eigenstates()
    labels = {}
    tracked_targets = {
        basis_label(0, n, 0): model.basis_state(0, n, 0) for n in range(0, 5)
    } | {
        basis_label(2, n, 0): model.basis_state(2, n, 0) for n in range(0, 4)
    } | {
        basis_label(0, n, 1): model.basis_state(0, n, 1) for n in range(0, 4)
    }
    for label, ket in tracked_targets.items():
        overlaps = [abs(ket.overlap(eig)) ** 2 for eig in eigenstates]
        idx = int(np.argmax(overlaps))
        labels[label] = {
            "eigen_index": idx,
            "eigen_energy_rot_MHz": round(hz(eigenvalues[idx]) / 1.0e6, 6),
            "max_overlap": float(overlaps[idx]),
        }

    for row in rows:
        n = int(row["n"])
        g_storage = labels[basis_label(0, n, 0)]
        f_storage = labels[basis_label(2, n - 1, 0)]
        g_readout = labels[basis_label(0, n - 1, 1)]
        row["energy_g_0r_ns_rot_MHz"] = g_storage["eigen_energy_rot_MHz"]
        row["energy_f_0r_nsm1_rot_MHz"] = f_storage["eigen_energy_rot_MHz"]
        row["energy_g_1r_nsm1_rot_MHz"] = g_readout["eigen_energy_rot_MHz"]
        row["overlap_g_0r_ns"] = round(float(g_storage["max_overlap"]), 6)
        row["overlap_f_0r_nsm1"] = round(float(f_storage["max_overlap"]), 6)
        row["overlap_g_1r_nsm1"] = round(float(g_readout["max_overlap"]), 6)

    return {"frequency_rows": rows, "spectrum_labels": labels}


def simulate_transfer(model, frame, *, n: int, family: str, amplitude_mhz: float, duration_ns: float, mode: str):
    duration_s = float(duration_ns) * 1.0e-9
    if mode == "storage":
        initial = model.basis_state(0, n, 0)
        target = model.basis_state(2, n - 1, 0)
        target_frequency = storage_sideband_rotating_frequency(model, frame, n)
        drive_target = SidebandDriveSpec(mode="storage", lower_level=0, upper_level=2, sideband="red")
        channel = "storage_sb"
    elif mode == "readout":
        initial = model.basis_state(2, n - 1, 0)
        target = model.basis_state(0, n - 1, 1)
        target_frequency = readout_dump_rotating_frequency(model, frame, n)
        drive_target = SidebandDriveSpec(mode="readout", lower_level=0, upper_level=2, sideband="red")
        channel = "readout_sb"
    else:
        raise ValueError(f"Unsupported mode '{mode}'.")

    pulse = make_family_pulse(
        channel=channel,
        carrier_rad_s=target_frequency,
        duration_s=duration_s,
        amplitude_hz=float(amplitude_mhz) * 1.0e6,
        family=family,
        label=f"{mode}_{family}_n{n}",
    )
    solver_failure = None
    try:
        _, result = simulate_single_stage(
            model,
            initial,
            pulse=pulse,
            duration_s=duration_s,
            drive_ops={channel: drive_target},
            frame=frame,
            noise=None,
            dt_s=DEFAULT_DT_S,
            store_states=False,
        )
        final_state = result.final_state
        p_target = state_population(final_state, target)
        p_initial = state_population(final_state, initial)
        leakage = max(0.0, 1.0 - p_target - p_initial)
        dominant_basis_states = sorted_basis_populations(final_state, model, cutoff=3)
    except Exception as exc:
        p_target = 0.0
        p_initial = 0.0
        leakage = 1.0
        dominant_basis_states = []
        solver_failure = str(exc)
    return {
        "mode": mode,
        "n": n,
        "family": family,
        "paper_role": str(PULSE_FAMILY_INFO[str(family)]["paper_role"]),
        "amplitude_MHz": float(amplitude_mhz),
        "duration_ns": float(duration_ns),
        "target_probability": float(p_target),
        "return_probability": float(p_initial),
        "leakage_probability": float(leakage),
        "dominant_basis_states": dominant_basis_states,
        "solver_failure": solver_failure,
    }


def choose_practical_best(cases: list[dict]) -> dict:
    """Prefer the shortest high-fidelity pulse over long many-cycle revivals."""

    good_cases = [case for case in cases if float(case["target_probability"]) >= 0.98]
    candidates = good_cases if good_cases else cases
    return min(
        candidates,
        key=lambda case: (
            float(case["duration_ns"]),
            -float(case["target_probability"]),
            float(case["leakage_probability"]),
            float(case["amplitude_MHz"]),
        ),
    )


def best_case_per_family(cases: list[dict]) -> dict[str, dict]:
    summary: dict[str, dict] = {}
    for family in PULSE_FAMILIES:
        family_cases = [case for case in cases if str(case["family"]) == str(family)]
        if family_cases:
            summary[str(family)] = choose_practical_best(family_cases)
    return summary


def optimize_pulses(model, frame):
    storage_scan_rows = []
    dump_scan_rows = []
    best_storage = {}
    best_dump = {}

    for n in range(1, 5):
        storage_cases_n = []
        dump_cases_n = []
        for family in PULSE_FAMILIES:
            for amplitude_mhz in STORAGE_AMPLITUDES_MHZ:
                for duration_ns in STORAGE_DURATIONS_NS:
                    case = simulate_transfer(
                        model,
                        frame,
                        n=n,
                        family=family,
                        amplitude_mhz=float(amplitude_mhz),
                        duration_ns=float(duration_ns),
                        mode="storage",
                    )
                    storage_scan_rows.append(case)
                    storage_cases_n.append(case)
            for amplitude_mhz in DUMP_AMPLITUDES_MHZ:
                for duration_ns in DUMP_DURATIONS_NS:
                    case = simulate_transfer(
                        model,
                        frame,
                        n=n,
                        family=family,
                        amplitude_mhz=float(amplitude_mhz),
                        duration_ns=float(duration_ns),
                        mode="readout",
                    )
                    dump_scan_rows.append(case)
                    dump_cases_n.append(case)
        best_storage[n] = choose_practical_best(storage_cases_n)
        best_dump[n] = choose_practical_best(dump_cases_n)
    return {
        "storage_scan_rows": storage_scan_rows,
        "dump_scan_rows": dump_scan_rows,
        "best_storage": best_storage,
        "best_dump": best_dump,
        "best_storage_by_family": {
            n: best_case_per_family([case for case in storage_scan_rows if int(case["n"]) == n]) for n in range(1, 5)
        },
        "best_dump_by_family": {
            n: best_case_per_family([case for case in dump_scan_rows if int(case["n"]) == n]) for n in range(1, 5)
        },
    }


def simulate_direct_carrier_control(model, frame):
    rows = []
    for n in range(1, 5):
        frequency = model.transmon_transition_frequency(
            storage_level=n,
            readout_level=0,
            lower_level=0,
            upper_level=2,
            frame=frame,
        )
        pulse = make_family_pulse(
            channel="gf_carrier",
            carrier_rad_s=float(frequency),
            duration_s=50.0e-9,
            amplitude_hz=5.0e6,
            family="gaussian",
            label=f"direct_gf_n{n}",
        )
        _, result = simulate_single_stage(
            model,
            model.basis_state(0, n, 0),
            pulse=pulse,
            duration_s=50.0e-9,
            drive_ops={"gf_carrier": TransmonTransitionDriveSpec(lower_level=0, upper_level=2)},
            frame=frame,
            noise=None,
            dt_s=DEFAULT_DT_S,
            store_states=False,
        )
        final_state = result.final_state
        rows.append(
            {
                "n": n,
                "carrier_rot_MHz": round(hz(frequency) / 1.0e6, 6),
                "p_f_n": float(state_population(final_state, model.basis_state(2, n, 0))),
                "p_f_n_minus_1": float(state_population(final_state, model.basis_state(2, n - 1, 0))),
                "p_g_n_minus_1": float(state_population(final_state, model.basis_state(0, n - 1, 0))),
                "comment": "Direct transmon g-f carrier excitation does not remove a storage photon.",
            }
        )
    return rows


def sensitivity_scan(model, frame, best_storage):
    rows = []
    for n in range(1, 5):
        best = best_storage[n]
        base_frequency = storage_sideband_rotating_frequency(model, frame, n)
        for detuning_mhz in DETUNING_OFFSETS_MHZ:
            for amp_scale in AMPLITUDE_SCALE_GRID:
                pulse = make_family_pulse(
                    channel="storage_sb",
                    carrier_rad_s=base_frequency + 2.0 * np.pi * float(detuning_mhz) * 1.0e6,
                    duration_s=float(best["duration_ns"]) * 1.0e-9,
                    amplitude_hz=float(best["amplitude_MHz"]) * 1.0e6 * float(amp_scale),
                    family=str(best["family"]),
                    label=f"sensitivity_n{n}",
                )
                _, result = simulate_single_stage(
                    model,
                    model.basis_state(0, n, 0),
                    pulse=pulse,
                    duration_s=float(best["duration_ns"]) * 1.0e-9,
                    drive_ops={"storage_sb": SidebandDriveSpec(mode="storage", lower_level=0, upper_level=2, sideband="red")},
                    frame=frame,
                    noise=None,
                    dt_s=DEFAULT_DT_S,
                    store_states=False,
                )
                rows.append(
                    {
                        "n": n,
                        "detuning_MHz": float(detuning_mhz),
                        "amplitude_scale": float(amp_scale),
                        "target_probability": float(state_population(result.final_state, model.basis_state(2, n - 1, 0))),
                    }
                )
    return rows


def calibration_frequency_scan(model, frame, best_storage, best_dump):
    rows = []
    summaries = []
    for mode, best_cases, initial_fn, target_fn, frequency_fn, lab_frequency_fn, channel in (
        (
            "storage",
            best_storage,
            lambda n: model.basis_state(0, n, 0),
            lambda n: model.basis_state(2, n - 1, 0),
            storage_sideband_rotating_frequency,
            storage_sideband_lab_frequency,
            "storage_sb",
        ),
        (
            "readout",
            best_dump,
            lambda n: model.basis_state(2, n - 1, 0),
            lambda n: model.basis_state(0, n - 1, 1),
            readout_dump_rotating_frequency,
            readout_dump_lab_frequency,
            "readout_sb",
        ),
    ):
        for n in range(1, 5):
            best = best_cases[n]
            base_frequency = frequency_fn(model, frame, n)
            drive_target = SidebandDriveSpec(mode=mode, lower_level=0, upper_level=2, sideband="red")
            best_row = None
            for detuning_mhz in CALIBRATION_DETUNING_OFFSETS_MHZ:
                pulse = make_family_pulse(
                    channel=channel,
                    carrier_rad_s=base_frequency + 2.0 * np.pi * float(detuning_mhz) * 1.0e6,
                    duration_s=float(best["duration_ns"]) * 1.0e-9,
                    amplitude_hz=float(best["amplitude_MHz"]) * 1.0e6,
                    family=str(best["family"]),
                    label=f"{mode}_freq_cal_n{n}",
                )
                _, result = simulate_single_stage(
                    model,
                    initial_fn(n),
                    pulse=pulse,
                    duration_s=float(best["duration_ns"]) * 1.0e-9,
                    drive_ops={channel: drive_target},
                    frame=frame,
                    noise=None,
                    dt_s=DEFAULT_DT_S,
                    store_states=False,
                )
                row = {
                    "mode": mode,
                    "n": n,
                    "family": str(best["family"]),
                    "amplitude_MHz": float(best["amplitude_MHz"]),
                    "duration_ns": float(best["duration_ns"]),
                    "detuning_MHz": float(detuning_mhz),
                    "target_probability": float(state_population(result.final_state, target_fn(n))),
                }
                rows.append(row)
                if best_row is None or row["target_probability"] > best_row["target_probability"]:
                    best_row = row
            assert best_row is not None
            summaries.append(
                {
                    "mode": mode,
                    "n": n,
                    "family": str(best["family"]),
                    "static_frequency_lab_GHz": round(hz(lab_frequency_fn(model, n)) / 1.0e9, 9),
                    "optimal_detuning_MHz": float(best_row["detuning_MHz"]),
                    "optimal_target_probability": float(best_row["target_probability"]),
                }
            )
    return {"rows": rows, "summary": summaries}


def floquet_sideband_summary(best_storage, best_dump):
    floquet_model = build_model(**FLOQUET_TRUNCATION)
    floquet_frame = build_frame(floquet_model)
    rows = []
    for mode, best_cases, frequency_fn, bare_pair_fn in (
        (
            "storage",
            best_storage,
            storage_sideband_rotating_frequency,
            lambda n: [floquet_model.basis_state(0, n, 0), floquet_model.basis_state(2, n - 1, 0)],
        ),
        (
            "readout",
            best_dump,
            readout_dump_rotating_frequency,
            lambda n: [floquet_model.basis_state(2, n - 1, 0), floquet_model.basis_state(0, n - 1, 1)],
        ),
    ):
        for n in range(1, 5):
            case = best_cases[n]
            drive_frequency = abs(frequency_fn(floquet_model, floquet_frame, n))
            drive = build_target_drive_term(
                floquet_model,
                SidebandDriveSpec(mode=mode, lower_level=0, upper_level=2, sideband="red"),
                amplitude=2.0 * np.pi * float(case["amplitude_MHz"]) * 1.0e6,
                frequency=drive_frequency,
                waveform="cos",
                label=f"{mode}_floquet_n{n}",
            )
            result = solve_floquet(
                FloquetProblem(
                    model=floquet_model,
                    frame=floquet_frame,
                    periodic_terms=[drive],
                    period=2.0 * np.pi / drive_frequency,
                    label=f"{mode}_floquet_problem_n{n}",
                ),
                FloquetConfig(n_time_samples=96),
            )
            modes = result.modes(0.0)
            target_pair = bare_pair_fn(n)
            supports = []
            for idx, mode_state in enumerate(modes):
                support = sum(abs(ket.overlap(mode_state)) ** 2 for ket in target_pair)
                supports.append((idx, float(support)))
            supports.sort(key=lambda item: item[1], reverse=True)
            dominant = supports[:2]
            i0, i1 = dominant[0][0], dominant[1][0]
            quasienergy_split_hz = hz(abs(result.quasienergies[i0] - result.quasienergies[i1]))
            expected_split_hz = float(case["amplitude_MHz"]) * 1.0e6 * (math.sqrt(float(n)) if mode == "storage" else 1.0)
            rows.append(
                {
                    "mode": mode,
                    "n": n,
                    "drive_amplitude_MHz": float(case["amplitude_MHz"]),
                    "dominant_pair_support_1": float(dominant[0][1]),
                    "dominant_pair_support_2": float(dominant[1][1]),
                    "quasienergy_split_MHz": round(quasienergy_split_hz / 1.0e6, 6),
                    "expected_bosonic_split_MHz": round(expected_split_hz / 1.0e6, 6),
                    "split_minus_expected_kHz": round((quasienergy_split_hz - expected_split_hz) / 1.0e3, 6),
                    "warnings": list(result.warnings),
                }
            )
    return {
        "rows": rows,
        "floquet_truncation": dict(FLOQUET_TRUNCATION),
    }


def simulate_cooling_primitive(model, frame, noise, *, n: int, best_storage_case: dict, best_dump_case: dict):
    storage_duration_s = float(best_storage_case["duration_ns"]) * 1.0e-9
    dump_duration_s = float(best_dump_case["duration_ns"]) * 1.0e-9
    ringdown_s = ladder_sequence_times(noise)

    storage_pulse = make_family_pulse(
        channel="storage_sb",
        carrier_rad_s=storage_sideband_rotating_frequency(model, frame, n),
        duration_s=storage_duration_s,
        amplitude_hz=float(best_storage_case["amplitude_MHz"]) * 1.0e6,
        family=str(best_storage_case["family"]),
        label=f"storage_sb_n{n}",
    )
    dump_pulse = make_family_pulse(
        channel="readout_sb",
        carrier_rad_s=readout_dump_rotating_frequency(model, frame, n),
        duration_s=dump_duration_s,
        amplitude_hz=float(best_dump_case["amplitude_MHz"]) * 1.0e6,
        family=str(best_dump_case["family"]),
        label=f"readout_sb_n{n}",
    )

    initial = model.basis_state(0, n, 0)
    storage_compiled, storage_result = simulate_single_stage(
        model,
        initial,
        pulse=storage_pulse,
        duration_s=storage_duration_s,
        drive_ops={"storage_sb": SidebandDriveSpec(mode="storage", lower_level=0, upper_level=2, sideband="red")},
        frame=frame,
        noise=noise,
        dt_s=DEFAULT_DT_S,
        store_states=True,
    )
    dump_compiled, dump_result = simulate_single_stage(
        model,
        storage_result.final_state,
        pulse=dump_pulse,
        duration_s=dump_duration_s,
        drive_ops={"readout_sb": SidebandDriveSpec(mode="readout", lower_level=0, upper_level=2, sideband="red")},
        frame=frame,
        noise=noise,
        dt_s=DEFAULT_DT_S,
        store_states=True,
    )
    idle_compiled, idle_result = simulate_single_stage(
        model,
        dump_result.final_state,
        pulse=None,
        duration_s=ringdown_s,
        drive_ops={},
        frame=frame,
        noise=noise,
        dt_s=4.0e-9,
        store_states=True,
    )
    final_state = idle_result.final_state
    return {
        "n": n,
        "initial_mean_storage_n": float(storage_photon_number(initial)),
        "after_storage_mean_storage_n": float(storage_photon_number(storage_result.final_state)),
        "after_dump_mean_storage_n": float(storage_photon_number(dump_result.final_state)),
        "final_mean_storage_n": float(storage_photon_number(final_state)),
        "final_readout_n": float(readout_photon_number(final_state)),
        "final_transmon_excited_population": float(transmon_excited_population(final_state)),
        "success_probability": float(state_population(final_state, model.basis_state(0, n - 1, 0))),
        "residual_same_n_probability": float(state_population(final_state, model.basis_state(0, n, 0))),
        "storage_stage_top": sorted_basis_populations(storage_result.final_state, model, cutoff=5),
        "dump_stage_top": sorted_basis_populations(dump_result.final_state, model, cutoff=5),
        "final_stage_top": sorted_basis_populations(final_state, model, cutoff=5),
        "storage_t_ns": list((storage_compiled.tlist * 1.0e9).tolist()),
        "storage_target_p": [float(state_population(state, model.basis_state(2, n - 1, 0))) for state in storage_result.states],
        "dump_t_ns": list((dump_compiled.tlist * 1.0e9).tolist()),
        "dump_target_p": [float(state_population(state, model.basis_state(0, n - 1, 1))) for state in dump_result.states],
        "ringdown_t_ns": list((idle_compiled.tlist * 1.0e9).tolist()),
        "ringdown_storage_n": [float(storage_photon_number(state)) for state in idle_result.states],
        "ringdown_readout_n": [float(readout_photon_number(state)) for state in idle_result.states],
    }


def simulate_ladder_protocol(model, frame, noise, best_storage, best_dump, initial_state):
    state = initial_state
    cycle_rows = []
    for n in range(4, 0, -1):
        storage_pulse = make_family_pulse(
            channel="storage_sb",
            carrier_rad_s=storage_sideband_rotating_frequency(model, frame, n),
            duration_s=float(best_storage[n]["duration_ns"]) * 1.0e-9,
            amplitude_hz=float(best_storage[n]["amplitude_MHz"]) * 1.0e6,
            family=str(best_storage[n]["family"]),
            label=f"ladder_storage_n{n}",
        )
        _, storage_result = simulate_single_stage(
            model,
            state,
            pulse=storage_pulse,
            duration_s=float(best_storage[n]["duration_ns"]) * 1.0e-9,
            drive_ops={"storage_sb": SidebandDriveSpec(mode="storage", lower_level=0, upper_level=2, sideband="red")},
            frame=frame,
            noise=noise,
            dt_s=DEFAULT_DT_S,
            store_states=False,
        )
        dump_pulse = make_family_pulse(
            channel="readout_sb",
            carrier_rad_s=readout_dump_rotating_frequency(model, frame, n),
            duration_s=float(best_dump[n]["duration_ns"]) * 1.0e-9,
            amplitude_hz=float(best_dump[n]["amplitude_MHz"]) * 1.0e6,
            family=str(best_dump[n]["family"]),
            label=f"ladder_dump_n{n}",
        )
        _, dump_result = simulate_single_stage(
            model,
            storage_result.final_state,
            pulse=dump_pulse,
            duration_s=float(best_dump[n]["duration_ns"]) * 1.0e-9,
            drive_ops={"readout_sb": SidebandDriveSpec(mode="readout", lower_level=0, upper_level=2, sideband="red")},
            frame=frame,
            noise=noise,
            dt_s=DEFAULT_DT_S,
            store_states=False,
        )
        _, idle_result = simulate_single_stage(
            model,
            dump_result.final_state,
            pulse=None,
            duration_s=ladder_sequence_times(noise),
            drive_ops={},
            frame=frame,
            noise=noise,
            dt_s=4.0e-9,
            store_states=False,
        )
        state = idle_result.final_state
        cycle_rows.append(
            {
                "target_n": n,
                "final_mean_storage_n": float(storage_photon_number(state)),
                "final_readout_n": float(readout_photon_number(state)),
                "final_transmon_excited_population": float(transmon_excited_population(state)),
            }
        )
    return {"final_state": state, "cycles": cycle_rows}


def validation_summary(model, frame, best_storage):
    rows = []
    for n in range(1, 5):
        predicted_rot = hz(storage_sideband_rotating_frequency(model, frame, n))
        analytic_rot = DEVICE.qubit_anharmonicity_hz + 2.0 * DEVICE.chi_storage_hz * (n - 1) - DEVICE.storage_kerr_hz * (n - 1)
        rows.append(
            {
                "n": n,
                "predicted_rot_MHz": round(predicted_rot / 1.0e6, 6),
                "analytic_rot_MHz": round(analytic_rot / 1.0e6, 6),
                "difference_kHz": round((predicted_rot - analytic_rot) / 1.0e3, 9),
                "expected_pi_ns": round(
                    expected_pi_duration_s(
                        amplitude_hz=float(best_storage[n]["amplitude_MHz"]) * 1.0e6,
                        matrix_element=math.sqrt(float(n)),
                    )
                    * 1.0e9,
                    6,
                ),
                "optimized_pi_ns": round(float(best_storage[n]["duration_ns"]), 6),
            }
        )
    return rows


def combined_recommendation_table(frequency_rows, best_storage, best_dump, cooling_rows, calibration_summary):
    calibration_lookup = {
        (str(row["mode"]), int(row["n"])): row for row in calibration_summary
    }
    rows = []
    for row in frequency_rows:
        n = int(row["n"])
        storage_case = best_storage[n]
        dump_case = best_dump[n]
        cooling = cooling_rows[n]
        storage_cal = calibration_lookup[("storage", n)]
        dump_cal = calibration_lookup[("readout", n)]
        rows.append(
            {
                "n": n,
                "step_a_frequency_GHz": row["storage_sideband_lab_GHz"],
                "step_b_frequency_GHz": row["readout_dump_lab_GHz"],
                "step_a_family": str(storage_case["family"]),
                "step_b_family": str(dump_case["family"]),
                "step_a_duration_ns": float(storage_case["duration_ns"]),
                "step_b_duration_ns": float(dump_case["duration_ns"]),
                "step_a_amplitude_MHz": float(storage_case["amplitude_MHz"]),
                "step_b_amplitude_MHz": float(dump_case["amplitude_MHz"]),
                "step_a_transfer": float(storage_case["target_probability"]),
                "step_b_transfer": float(dump_case["target_probability"]),
                "single_cycle_success": float(cooling["success_probability"]),
                "step_a_optimal_detuning_MHz": float(storage_cal["optimal_detuning_MHz"]),
                "step_b_optimal_detuning_MHz": float(dump_cal["optimal_detuning_MHz"]),
                "dominant_leakage_channel": str(cooling["final_stage_top"][1][0]) if len(cooling["final_stage_top"]) > 1 else "none",
                "expected_experimental_difficulty": "moderate" if n <= 2 else "high",
            }
        )
    return rows


def make_figures(
    frequency_data,
    pulse_data,
    sensitivity_rows,
    calibration_data,
    floquet_data,
    cooling_rows,
    ladder_data,
    coherent_data,
    thermal_data,
):
    frequency_rows = frequency_data["frequency_rows"]
    best_storage = pulse_data["best_storage"]

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    ns = [row["n"] for row in frequency_rows]
    ax.plot(ns, [row["storage_sideband_lab_GHz"] for row in frequency_rows], marker="o", label="Storage conversion")
    ax.plot(ns, [row["readout_dump_lab_GHz"] for row in frequency_rows], marker="s", label="Readout dump")
    ax.set_xlabel("Storage Fock number n")
    ax.set_ylabel("Lab-frame frequency (GHz)")
    ax.set_title("gf active-cooling transition map")
    ax.legend()
    plot_save(fig, "transition_map")

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.6), sharex=True, sharey=True)
    for ax, n in zip(axes.flat, range(1, 5), strict=True):
        subset = [row for row in pulse_data["storage_scan_rows"] if row["n"] == n and row["family"] == best_storage[n]["family"]]
        amp_vals = sorted({row["amplitude_MHz"] for row in subset})
        dur_vals = sorted({row["duration_ns"] for row in subset})
        heat = np.zeros((len(amp_vals), len(dur_vals)))
        for row in subset:
            i = amp_vals.index(row["amplitude_MHz"])
            j = dur_vals.index(row["duration_ns"])
            heat[i, j] = row["target_probability"]
        im = ax.imshow(
            heat,
            aspect="auto",
            origin="lower",
            extent=(dur_vals[0], dur_vals[-1], amp_vals[0], amp_vals[-1]),
            vmin=0.0,
            vmax=1.0,
            cmap="viridis",
        )
        ax.set_title(f"n={n}, {best_storage[n]['family']}")
        ax.scatter([best_storage[n]["duration_ns"]], [best_storage[n]["amplitude_MHz"]], c="r", s=18)
    axes[1, 0].set_xlabel("Duration (ns)")
    axes[1, 1].set_xlabel("Duration (ns)")
    axes[0, 0].set_ylabel("Amplitude (MHz)")
    axes[1, 0].set_ylabel("Amplitude (MHz)")
    fig.colorbar(im, ax=axes.ravel().tolist(), label="Target transfer probability")
    plot_save(fig, "storage_scan_heatmaps")

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.6), sharex=True, sharey=True)
    for ax, n in zip(axes.flat, range(1, 5), strict=True):
        family_cases = pulse_data["best_storage_by_family"][n]
        families = list(family_cases.keys())
        transfers = [float(family_cases[family]["target_probability"]) for family in families]
        ax.bar(range(len(families)), transfers, color="tab:blue")
        ax.set_title(f"Step A family screen, n={n}")
        ax.set_xticks(range(len(families)))
        ax.set_xticklabels(families, rotation=45, ha="right")
        ax.set_ylim(0.0, 1.02)
    axes[0, 0].set_ylabel("Best transfer probability")
    axes[1, 0].set_ylabel("Best transfer probability")
    plot_save(fig, "storage_family_comparison")

    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    primitive = cooling_rows[4]
    t_storage = np.asarray(primitive["storage_t_ns"], dtype=float)
    t_dump = np.asarray(primitive["dump_t_ns"], dtype=float) + (t_storage[-1] if len(t_storage) else 0.0)
    t_ring = np.asarray(primitive["ringdown_t_ns"], dtype=float) + (t_dump[-1] if len(t_dump) else 0.0)
    ax.plot(t_storage, primitive["storage_target_p"], label="P(|f,0,3>) during storage conversion")
    ax.plot(t_dump, primitive["dump_target_p"], label="P(|g,1,3>) during readout dump")
    ax.plot(t_ring, primitive["ringdown_storage_n"], label="<n_s> during ringdown")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Population / mean occupation")
    ax.set_title("Representative n=4 cooling primitive")
    ax.legend()
    plot_save(fig, "cooling_primitive_trajectories")

    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    ax.plot(range(1, 5), [row["final_mean_storage_n"] for row in ladder_data["cycles"]], marker="o", label="Initial |g,0,4>")
    ax.plot(range(1, 5), [row["final_mean_storage_n"] for row in coherent_data["cycles"]], marker="s", label="Initial coherent state")
    ax.plot(range(1, 5), [row["final_mean_storage_n"] for row in thermal_data["cycles"]], marker="^", label="Initial thermal mixture")
    ax.set_xlabel("Cooling stage index")
    ax.set_ylabel("Final mean storage occupation")
    ax.set_title("Repeated ladder cooling performance")
    ax.legend()
    plot_save(fig, "cooling_per_cycle")

    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.6), sharex=True, sharey=True)
    for ax, n in zip(axes.flat, range(1, 5), strict=True):
        subset = [row for row in sensitivity_rows if row["n"] == n]
        amp_vals = sorted({row["amplitude_scale"] for row in subset})
        det_vals = sorted({row["detuning_MHz"] for row in subset})
        heat = np.zeros((len(amp_vals), len(det_vals)))
        for row in subset:
            i = amp_vals.index(row["amplitude_scale"])
            j = det_vals.index(row["detuning_MHz"])
            heat[i, j] = row["target_probability"]
        im = ax.imshow(
            heat,
            aspect="auto",
            origin="lower",
            extent=(det_vals[0], det_vals[-1], amp_vals[0], amp_vals[-1]),
            vmin=0.0,
            vmax=1.0,
            cmap="magma",
        )
        ax.set_title(f"n={n}")
    axes[1, 0].set_xlabel("Detuning error (MHz)")
    axes[1, 1].set_xlabel("Detuning error (MHz)")
    axes[0, 0].set_ylabel("Amplitude scale")
    axes[1, 0].set_ylabel("Amplitude scale")
    fig.colorbar(im, ax=axes.ravel().tolist(), label="Target transfer probability")
    plot_save(fig, "sensitivity_heatmaps")

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2), sharey=True)
    for ax, mode in zip(axes, ("storage", "readout"), strict=True):
        for n in range(1, 5):
            subset = [row for row in calibration_data["rows"] if row["mode"] == mode and int(row["n"]) == n]
            subset.sort(key=lambda row: float(row["detuning_MHz"]))
            ax.plot(
                [float(row["detuning_MHz"]) for row in subset],
                [float(row["target_probability"]) for row in subset],
                marker="o",
                label=f"n={n}",
            )
        ax.axvline(0.0, color="black", linestyle="--", linewidth=1.0)
        ax.set_title(f"{mode.capitalize()} frequency calibration")
        ax.set_xlabel("Detuning from static resonance (MHz)")
    axes[0].set_ylabel("Target transfer probability")
    axes[1].legend()
    plot_save(fig, "frequency_calibration_curves")

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2), sharey=True)
    for ax, mode in zip(axes, ("storage", "readout"), strict=True):
        subset = [row for row in floquet_data["rows"] if row["mode"] == mode]
        subset.sort(key=lambda row: int(row["n"]))
        ax.plot(
            [int(row["n"]) for row in subset],
            [float(row["quasienergy_split_MHz"]) for row in subset],
            marker="o",
            label="Floquet split",
        )
        ax.plot(
            [int(row["n"]) for row in subset],
            [float(row["expected_bosonic_split_MHz"]) for row in subset],
            marker="s",
            linestyle="--",
            label="Expected split",
        )
        ax.set_title(f"{mode.capitalize()} Floquet doublet")
        ax.set_xlabel("Storage Fock number n")
    axes[0].set_ylabel("Split (MHz)")
    axes[1].legend()
    plot_save(fig, "floquet_doublet_splitting")

    storage_rho = reduced_storage_state(coherent_data["final_state"])
    xvec, yvec, wigner = cavity_wigner(storage_rho, n_points=81, extent=3.5, coordinate="quadrature")
    fig, ax = plt.subplots(figsize=(5.5, 4.8))
    pcm = ax.pcolormesh(xvec, yvec, wigner, shading="auto", cmap="RdBu_r")
    ax.set_xlabel("x")
    ax.set_ylabel("p")
    ax.set_title("Storage Wigner after ladder cooling")
    fig.colorbar(pcm, ax=ax, label="Wigner")
    plot_save(fig, "coherent_final_wigner")


def main() -> None:
    start = time.time()
    model = build_model()
    frame = build_frame(model)
    noise = build_noise()

    write_device_manifest(ARTIFACTS_DIR / "device_manifest.json")

    frequency_data = spectrum_and_frequency_table(model, frame)
    csv_dump(DATA_DIR / "frequency_table.csv", frequency_data["frequency_rows"])
    json_dump(ARTIFACTS_DIR / "frequency_analysis.json", frequency_data)

    direct_rows = simulate_direct_carrier_control(model, frame)
    json_dump(ARTIFACTS_DIR / "direct_carrier_comparison.json", direct_rows)

    pulse_data = optimize_pulses(model, frame)
    csv_dump(DATA_DIR / "storage_scan_summary.csv", pulse_data["storage_scan_rows"])
    csv_dump(DATA_DIR / "dump_scan_summary.csv", pulse_data["dump_scan_rows"])
    json_dump(ARTIFACTS_DIR / "pulse_optimization.json", pulse_data)

    sensitivity_rows = sensitivity_scan(model, frame, pulse_data["best_storage"])
    csv_dump(DATA_DIR / "sensitivity_scan.csv", sensitivity_rows)
    json_dump(ARTIFACTS_DIR / "sensitivity_analysis.json", sensitivity_rows)

    calibration_data = calibration_frequency_scan(model, frame, pulse_data["best_storage"], pulse_data["best_dump"])
    csv_dump(DATA_DIR / "frequency_calibration_scan.csv", calibration_data["rows"])
    json_dump(ARTIFACTS_DIR / "frequency_calibration_scan.json", calibration_data)

    floquet_data = floquet_sideband_summary(pulse_data["best_storage"], pulse_data["best_dump"])
    csv_dump(DATA_DIR / "floquet_summary.csv", floquet_data["rows"])
    json_dump(ARTIFACTS_DIR / "floquet_summary.json", floquet_data)

    cooling_rows = {
        n: simulate_cooling_primitive(
            model,
            frame,
            noise,
            n=n,
            best_storage_case=pulse_data["best_storage"][n],
            best_dump_case=pulse_data["best_dump"][n],
        )
        for n in range(1, 5)
    }
    json_dump(ARTIFACTS_DIR / "cooling_primitive_metrics.json", cooling_rows)

    ladder_data = simulate_ladder_protocol(
        model,
        frame,
        noise,
        pulse_data["best_storage"],
        pulse_data["best_dump"],
        model.basis_state(0, 4, 0),
    )
    coherent_storage = qt.coherent(model.n_storage, 1.1)
    coherent_initial = qt.tensor(qt.basis(model.n_tr, 0), coherent_storage, qt.basis(model.n_readout, 0))
    coherent_data = simulate_ladder_protocol(
        model,
        frame,
        noise,
        pulse_data["best_storage"],
        pulse_data["best_dump"],
        coherent_initial,
    )
    probs = np.array([0.44, 0.28, 0.16, 0.08, 0.04], dtype=float)
    probs = probs / probs.sum()
    thermal_storage = 0.0 * qt.basis(model.n_storage, 0).proj()
    for n, prob in enumerate(probs):
        thermal_storage = thermal_storage + prob * qt.basis(model.n_storage, n).proj()
    thermal_initial = qt.tensor(qt.basis(model.n_tr, 0).proj(), thermal_storage, qt.basis(model.n_readout, 0).proj())
    thermal_data = simulate_ladder_protocol(
        model,
        frame,
        noise,
        pulse_data["best_storage"],
        pulse_data["best_dump"],
        thermal_initial,
    )
    json_dump(
        ARTIFACTS_DIR / "ladder_protocols.json",
        {
            "basis_g4": {"cycles": ladder_data["cycles"]},
            "coherent": {"cycles": coherent_data["cycles"]},
            "thermal": {"cycles": thermal_data["cycles"]},
        },
    )

    validation_rows = validation_summary(model, frame, pulse_data["best_storage"])
    json_dump(ARTIFACTS_DIR / "validation_summary.json", validation_rows)

    recommendation_rows = combined_recommendation_table(
        frequency_data["frequency_rows"],
        pulse_data["best_storage"],
        pulse_data["best_dump"],
        cooling_rows,
        calibration_data["summary"],
    )
    csv_dump(DATA_DIR / "recommendation_table.csv", recommendation_rows)
    json_dump(ARTIFACTS_DIR / "recommendation_table.json", recommendation_rows)

    make_figures(
        frequency_data,
        pulse_data,
        sensitivity_rows,
        calibration_data,
        floquet_data,
        cooling_rows,
        ladder_data,
        coherent_data,
        thermal_data,
    )

    results = {
        "device": asdict(DEVICE),
        "best_storage": pulse_data["best_storage"],
        "best_storage_by_family": pulse_data["best_storage_by_family"],
        "best_dump": pulse_data["best_dump"],
        "best_dump_by_family": pulse_data["best_dump_by_family"],
        "frequency_calibration": calibration_data["summary"],
        "floquet_summary": floquet_data["rows"],
        "recommendation_table": recommendation_rows,
        "cooling_primitives": cooling_rows,
        "ladder_basis_g4": ladder_data["cycles"],
        "ladder_coherent": coherent_data["cycles"],
        "ladder_thermal": thermal_data["cycles"],
        "validation": validation_rows,
        "runtime_s": time.time() - start,
    }
    json_dump(DATA_DIR / "study_results.json", results)
    print(json.dumps(results, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
