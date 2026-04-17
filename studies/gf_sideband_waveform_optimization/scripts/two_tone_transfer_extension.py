"""Simultaneous two-tone storage-to-readout transfer extension for the gf-sideband study."""

from __future__ import annotations

from dataclasses import asdict
import math
import time

import numpy as np

from cqed_sim.core.drive_targets import SidebandDriveSpec
from cqed_sim.sequence.scheduler import SequenceCompiler
from cqed_sim.sim.runner import SimulationConfig, simulate_sequence

from common import (
    ARTIFACTS_DIR,
    DATA_DIR,
    DEFAULT_FINAL_DT_S,
    DEFAULT_SWEEP_DT_S,
    DEVICE,
    FAMILY_VARIANTS,
    TRANSMON_REFERENCE,
    basis_label,
    build_frame,
    build_model,
    build_noise,
    csv_dump,
    ghz,
    json_dump,
    make_pulse,
    mhz,
    plot_save,
    sanitize_for_json,
    sorted_basis_populations,
    state_population,
    to_internal_units,
)

TWO_PI = 2.0 * np.pi
TRANSFER_N_VALUES = (1, 2, 3)
TARGET_COUPLING_GRID_MHZ = np.array([4.0, 6.0, 8.0, 10.0, 12.0, 14.0])
COMMON_DETUNING_GRID_MHZ = np.array([0.0, 4.0, 8.0, 12.0, 16.0, 20.0])
DETUNED_MIN_MHZ = 8.0
RESONANT_TARGET_THRESHOLD = 0.99
DETUNED_TARGET_THRESHOLD = 0.95
MAX_DETUNED_INTERMEDIATE_THRESHOLD = 0.20
MAX_SCAN_DURATION_S = 0.6e-6

SQUARE_VARIANT = next(variant for variant in FAMILY_VARIANTS if variant.family == "square")


def two_tone_noise_scenarios() -> list[dict[str, object]]:
    return [
        {
            "noise_scenario": "mode_only",
            "description": "Storage relaxation/dephasing and readout linewidth only.",
            "transmon_t1_s": None,
            "transmon_t2_ramsey_s": None,
            "transmon_tphi_s": None,
            "source": str(TRANSMON_REFERENCE.parameter_source),
        },
        {
            "noise_scenario": "transmon_reference",
            "description": "Mode noise plus the matched local transmon coherence reference.",
            "transmon_t1_s": TRANSMON_REFERENCE.qubit_t1_s,
            "transmon_t2_ramsey_s": TRANSMON_REFERENCE.qubit_t2_ramsey_s,
            "transmon_tphi_s": TRANSMON_REFERENCE.qubit_tphi_ramsey_s,
            "source": str(TRANSMON_REFERENCE.parameter_source),
        },
    ]


def transfer_states(model, n: int):
    source = model.basis_state(0, int(n), 0)
    intermediate = model.basis_state(2, int(n) - 1, 0)
    target = model.basis_state(0, int(n) - 1, 1)
    return source, intermediate, target


def storage_leg_rotating_frequency(model, frame, n: int) -> float:
    return float(
        model.sideband_transition_frequency(
            mode="storage",
            storage_level=int(n) - 1,
            readout_level=0,
            lower_level=0,
            upper_level=2,
            sideband="red",
            frame=frame,
        )
    )


def readout_leg_rotating_frequency(model, frame, n: int) -> float:
    return float(
        model.sideband_transition_frequency(
            mode="readout",
            storage_level=int(n) - 1,
            readout_level=0,
            lower_level=0,
            upper_level=2,
            sideband="red",
            frame=frame,
        )
    )


def storage_leg_lab_frequency(model, n: int) -> float:
    source, intermediate, _ = transfer_states(model, n)
    del source
    return float(model.basis_energy(2, int(n) - 1, 0) - model.basis_energy(0, int(n), 0))


def readout_leg_lab_frequency(model, n: int) -> float:
    _, intermediate, target = transfer_states(model, n)
    del intermediate, target
    return float(model.basis_energy(2, int(n) - 1, 0) - model.basis_energy(0, int(n) - 1, 1))


def analytic_storage_leg_frequency_hz(n: int) -> float:
    return float(DEVICE.qubit_anharmonicity_hz + 2.0 * DEVICE.chi_storage_hz * (int(n) - 1) - DEVICE.storage_kerr_hz * (int(n) - 1))


def analytic_readout_leg_frequency_hz(n: int) -> float:
    return float(DEVICE.qubit_anharmonicity_hz + (2.0 * DEVICE.chi_storage_hz - DEVICE.chi_storage_readout_hz) * (int(n) - 1))


def storage_matrix_element(model, n: int) -> float:
    _, storage_operator = model.sideband_drive_operators(mode="storage", lower_level=0, upper_level=2, sideband="red")
    source, intermediate, _ = transfer_states(model, n)
    return float(abs(source.overlap(storage_operator * intermediate)))


def readout_matrix_element(model, n: int) -> float:
    _, readout_operator = model.sideband_drive_operators(mode="readout", lower_level=0, upper_level=2, sideband="red")
    _, intermediate, target = transfer_states(model, n)
    return float(abs(target.overlap(readout_operator * intermediate)))


def reduced_model_dynamics(times_s: np.ndarray, *, g_storage_rad_s: float, g_readout_rad_s: float, detuning_rad_s: float) -> dict[str, np.ndarray]:
    hamiltonian = np.array(
        [
            [0.0, g_storage_rad_s, 0.0],
            [g_storage_rad_s, detuning_rad_s, g_readout_rad_s],
            [0.0, g_readout_rad_s, 0.0],
        ],
        dtype=np.complex128,
    )
    eigenvalues, eigenvectors = np.linalg.eigh(hamiltonian)
    initial = np.array([1.0, 0.0, 0.0], dtype=np.complex128)
    coefficients = eigenvectors.conj().T @ initial
    phases = np.exp(-1.0j * np.outer(eigenvalues, np.asarray(times_s, dtype=float)))
    states = eigenvectors @ (phases * coefficients[:, None])
    probabilities = np.abs(states) ** 2
    return {
        "source_curve": probabilities[0],
        "intermediate_curve": probabilities[1],
        "target_curve": probabilities[2],
    }


def resonant_peak_time_s(g_storage_rad_s: float, g_readout_rad_s: float) -> float:
    return float(np.pi / np.sqrt(g_storage_rad_s**2 + g_readout_rad_s**2))


def raman_peak_time_s(g_storage_rad_s: float, g_readout_rad_s: float, detuning_rad_s: float) -> float | None:
    if abs(detuning_rad_s) < 1.0e-15:
        return None
    return float(np.pi * abs(detuning_rad_s) / (2.0 * g_storage_rad_s * g_readout_rad_s))


def suggested_duration_s(*, g_storage_rad_s: float, g_readout_rad_s: float, detuning_rad_s: float) -> float:
    resonant_guess = resonant_peak_time_s(g_storage_rad_s, g_readout_rad_s)
    raman_guess = raman_peak_time_s(g_storage_rad_s, g_readout_rad_s, detuning_rad_s)
    if raman_guess is None:
        estimate = max(80.0e-9, 1.6 * resonant_guess)
    else:
        estimate = max(120.0e-9, 1.35 * max(resonant_guess, raman_guess))
    return float(min(estimate, MAX_SCAN_DURATION_S))


def two_tone_drive_targets() -> dict[str, SidebandDriveSpec]:
    return {
        "storage_two_tone": SidebandDriveSpec(mode="storage", lower_level=0, upper_level=2, sideband="red"),
        "readout_two_tone": SidebandDriveSpec(mode="readout", lower_level=0, upper_level=2, sideband="red"),
    }


def simulate_two_tone_case(
    model,
    frame,
    *,
    n: int,
    target_coupling_mhz: float,
    common_detuning_mhz: float,
    noise=None,
    dt_s: float = DEFAULT_SWEEP_DT_S,
    duration_s: float | None = None,
    include_curves: bool = False,
) -> dict[str, object]:
    source, intermediate, target = transfer_states(model, n)
    storage_frequency = storage_leg_rotating_frequency(model, frame, n)
    readout_frequency = readout_leg_rotating_frequency(model, frame, n)
    storage_element = storage_matrix_element(model, n)
    readout_element = readout_matrix_element(model, n)
    target_coupling_hz = float(target_coupling_mhz) * 1.0e6
    detuning_rad_s = TWO_PI * float(common_detuning_mhz) * 1.0e6
    g_storage_rad_s = to_internal_units(target_coupling_hz)
    g_readout_rad_s = to_internal_units(target_coupling_hz)
    storage_amplitude_hz = target_coupling_hz / max(storage_element, 1.0e-15)
    readout_amplitude_hz = target_coupling_hz / max(readout_element, 1.0e-15)
    if duration_s is None:
        duration_s = suggested_duration_s(
            g_storage_rad_s=g_storage_rad_s,
            g_readout_rad_s=g_readout_rad_s,
            detuning_rad_s=detuning_rad_s,
        )

    storage_pulse = make_pulse(
        channel="storage_two_tone",
        carrier_rad_s=storage_frequency + detuning_rad_s,
        duration_s=float(duration_s),
        amplitude_hz=storage_amplitude_hz,
        variant=SQUARE_VARIANT,
        label=f"two_tone_storage_n{n}_{target_coupling_mhz:.1f}MHz_{common_detuning_mhz:.1f}MHz",
    )
    readout_pulse = make_pulse(
        channel="readout_two_tone",
        carrier_rad_s=readout_frequency + detuning_rad_s,
        duration_s=float(duration_s),
        amplitude_hz=readout_amplitude_hz,
        variant=SQUARE_VARIANT,
        label=f"two_tone_readout_n{n}_{target_coupling_mhz:.1f}MHz_{common_detuning_mhz:.1f}MHz",
    )
    compiled = SequenceCompiler(dt=float(dt_s)).compile([storage_pulse, readout_pulse], t_end=float(max(duration_s, dt_s)))
    result = simulate_sequence(
        model,
        compiled,
        source,
        two_tone_drive_targets(),
        SimulationConfig(frame=frame, store_states=True, max_step=float(dt_s)),
        noise=noise,
    )

    times_s = np.asarray(compiled.tlist, dtype=float)
    source_curve = np.asarray([state_population(state, source) for state in result.states], dtype=float)
    intermediate_curve = np.asarray([state_population(state, intermediate) for state in result.states], dtype=float)
    target_curve = np.asarray([state_population(state, target) for state in result.states], dtype=float)
    leakage_curve = np.clip(1.0 - source_curve - intermediate_curve - target_curve, 0.0, None)
    peak_index = int(np.argmax(target_curve))
    reduced_curves = reduced_model_dynamics(
        times_s,
        g_storage_rad_s=g_storage_rad_s,
        g_readout_rad_s=g_readout_rad_s,
        detuning_rad_s=detuning_rad_s,
    )
    reduced_peak_index = int(np.argmax(reduced_curves["target_curve"]))
    payload: dict[str, object] = {
        "n": int(n),
        "mechanism": "resonant_bright_state" if abs(common_detuning_mhz) < 1.0e-12 else "detuned_raman_like",
        "target_coupling_MHz": float(target_coupling_mhz),
        "common_detuning_MHz": float(common_detuning_mhz),
        "storage_leg_lab_frequency_GHz": float(ghz(storage_leg_lab_frequency(model, n))),
        "readout_leg_lab_frequency_GHz": float(ghz(readout_leg_lab_frequency(model, n))),
        "storage_leg_rotating_frequency_MHz": float(mhz(storage_frequency)),
        "readout_leg_rotating_frequency_MHz": float(mhz(readout_frequency)),
        "storage_matrix_element": float(storage_element),
        "readout_matrix_element": float(readout_element),
        "storage_amplitude_MHz": float(storage_amplitude_hz / 1.0e6),
        "readout_amplitude_MHz": float(readout_amplitude_hz / 1.0e6),
        "simulation_duration_ns": float(times_s[-1] * 1.0e9),
        "peak_target_probability": float(target_curve[peak_index]),
        "time_to_peak_target_ns": float(times_s[peak_index] * 1.0e9),
        "max_intermediate_probability": float(np.max(intermediate_curve)),
        "intermediate_at_target_peak": float(intermediate_curve[peak_index]),
        "leakage_at_target_peak": float(leakage_curve[peak_index]),
        "final_target_probability": float(target_curve[-1]),
        "reduced_peak_target_probability": float(reduced_curves["target_curve"][reduced_peak_index]),
        "reduced_peak_time_ns": float(times_s[reduced_peak_index] * 1.0e9),
        "reduced_max_intermediate_probability": float(np.max(reduced_curves["intermediate_curve"])),
        "resonant_peak_time_ns": float(resonant_peak_time_s(g_storage_rad_s, g_readout_rad_s) * 1.0e9),
        "raman_peak_time_ns": None if common_detuning_mhz == 0.0 else float(raman_peak_time_s(g_storage_rad_s, g_readout_rad_s, detuning_rad_s) * 1.0e9),
        "peak_target_error_vs_reduced": float(target_curve[peak_index] - reduced_curves["target_curve"][reduced_peak_index]),
        "peak_time_error_vs_reduced_ns": float(times_s[peak_index] * 1.0e9 - times_s[reduced_peak_index] * 1.0e9),
        "dominant_peak_states": sorted_basis_populations(result.states[peak_index], model, cutoff=6),
        "dominant_final_states": sorted_basis_populations(result.states[-1], model, cutoff=6),
        "source_label": basis_label(0, int(n), 0),
        "intermediate_label": basis_label(2, int(n) - 1, 0),
        "target_label": basis_label(0, int(n) - 1, 1),
    }
    if include_curves:
        payload.update(
            {
                "times_ns": (times_s * 1.0e9).tolist(),
                "source_curve": source_curve.tolist(),
                "intermediate_curve": intermediate_curve.tolist(),
                "target_curve": target_curve.tolist(),
                "leakage_curve": leakage_curve.tolist(),
                "reduced_source_curve": reduced_curves["source_curve"].tolist(),
                "reduced_intermediate_curve": reduced_curves["intermediate_curve"].tolist(),
                "reduced_target_curve": reduced_curves["target_curve"].tolist(),
            }
        )
    return payload


def frequency_rows(model, frame) -> list[dict[str, object]]:
    rows = []
    for n in TRANSFER_N_VALUES:
        storage_exact_mhz = mhz(storage_leg_rotating_frequency(model, frame, n))
        readout_exact_mhz = mhz(readout_leg_rotating_frequency(model, frame, n))
        storage_analytic_mhz = analytic_storage_leg_frequency_hz(n) / 1.0e6
        readout_analytic_mhz = analytic_readout_leg_frequency_hz(n) / 1.0e6
        rows.append(
            {
                "n": int(n),
                "storage_leg_lab_frequency_GHz": float(ghz(storage_leg_lab_frequency(model, n))),
                "readout_leg_lab_frequency_GHz": float(ghz(readout_leg_lab_frequency(model, n))),
                "storage_leg_exact_rotating_frequency_MHz": float(storage_exact_mhz),
                "storage_leg_analytic_rotating_frequency_MHz": float(storage_analytic_mhz),
                "storage_leg_difference_kHz": float((storage_exact_mhz - storage_analytic_mhz) * 1.0e3),
                "readout_leg_exact_rotating_frequency_MHz": float(readout_exact_mhz),
                "readout_leg_analytic_rotating_frequency_MHz": float(readout_analytic_mhz),
                "readout_leg_difference_kHz": float((readout_exact_mhz - readout_analytic_mhz) * 1.0e3),
                "storage_matrix_element": float(storage_matrix_element(model, n)),
                "readout_matrix_element": float(readout_matrix_element(model, n)),
            }
        )
    return rows


def scan_rows(model, frame) -> list[dict[str, object]]:
    rows = []
    for n in TRANSFER_N_VALUES:
        for target_coupling_mhz in TARGET_COUPLING_GRID_MHZ:
            for common_detuning_mhz in COMMON_DETUNING_GRID_MHZ:
                rows.append(
                    simulate_two_tone_case(
                        model,
                        frame,
                        n=int(n),
                        target_coupling_mhz=float(target_coupling_mhz),
                        common_detuning_mhz=float(common_detuning_mhz),
                        dt_s=DEFAULT_SWEEP_DT_S,
                        include_curves=False,
                    )
                )
    return rows


def select_case_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    selected: list[dict[str, object]] = []
    for n in TRANSFER_N_VALUES:
        resonant_candidates = [row for row in rows if int(row["n"]) == int(n) and abs(float(row["common_detuning_MHz"])) < 1.0e-12]
        resonant_eligible = [
            row
            for row in resonant_candidates
            if float(row["peak_target_probability"]) >= RESONANT_TARGET_THRESHOLD and float(row["leakage_at_target_peak"]) <= 0.02
        ]
        if resonant_eligible:
            resonant_best = min(
                resonant_eligible,
                key=lambda row: (
                    float(row["time_to_peak_target_ns"]),
                    -float(row["peak_target_probability"]),
                    float(row["max_intermediate_probability"]),
                ),
            )
            resonant_status = "meets_threshold"
        else:
            resonant_best = max(
                resonant_candidates,
                key=lambda row: (float(row["peak_target_probability"]), -float(row["time_to_peak_target_ns"])),
            )
            resonant_status = "best_available"
        resonant_record = dict(resonant_best)
        resonant_record["case_role"] = "best_resonant"
        resonant_record["selection_status"] = resonant_status
        resonant_record["selection_rule"] = "Minimize transfer time subject to peak target >= 0.99 in the resonant manifold."
        selected.append(resonant_record)

        detuned_candidates = [
            row for row in rows if int(row["n"]) == int(n) and float(row["common_detuning_MHz"]) >= DETUNED_MIN_MHZ
        ]
        detuned_primary = [row for row in detuned_candidates if float(row["peak_target_probability"]) >= DETUNED_TARGET_THRESHOLD]
        detuned_low_f = [
            row
            for row in detuned_primary
            if float(row["max_intermediate_probability"]) <= MAX_DETUNED_INTERMEDIATE_THRESHOLD
        ]
        if detuned_low_f:
            detuned_best = min(
                detuned_low_f,
                key=lambda row: (
                    float(row["max_intermediate_probability"]),
                    float(row["time_to_peak_target_ns"]),
                    -float(row["peak_target_probability"]),
                ),
            )
            detuned_status = "meets_threshold"
        elif detuned_primary:
            detuned_best = min(
                detuned_primary,
                key=lambda row: (
                    float(row["max_intermediate_probability"]),
                    float(row["time_to_peak_target_ns"]),
                    -float(row["peak_target_probability"]),
                ),
            )
            detuned_status = "peak_only"
        else:
            detuned_best = max(
                detuned_candidates,
                key=lambda row: (
                    float(row["peak_target_probability"]) - 0.5 * float(row["max_intermediate_probability"]),
                    -float(row["time_to_peak_target_ns"]),
                ),
            )
            detuned_status = "best_available"
        detuned_record = dict(detuned_best)
        detuned_record["case_role"] = "best_detuned"
        detuned_record["selection_status"] = detuned_status
        detuned_record["selection_rule"] = "Minimize intermediate-state population subject to peak target >= 0.95 in the detuned Raman-like regime."
        selected.append(detuned_record)
    return selected


def trajectory_rows(model, frame, selected_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for row in selected_rows:
        rows.append(
            simulate_two_tone_case(
                model,
                frame,
                n=int(row["n"]),
                target_coupling_mhz=float(row["target_coupling_MHz"]),
                common_detuning_mhz=float(row["common_detuning_MHz"]),
                dt_s=DEFAULT_FINAL_DT_S,
                duration_s=float(row["simulation_duration_ns"]) * 1.0e-9,
                include_curves=True,
            )
            | {
                "case_role": str(row["case_role"]),
                "selection_status": str(row["selection_status"]),
            }
        )
    return rows


def open_system_rows(model, frame, selected_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for selected in selected_rows:
        for scenario in two_tone_noise_scenarios():
            noise = build_noise(
                transmon_t1_s=scenario["transmon_t1_s"],
                transmon_t2_ramsey_s=scenario["transmon_t2_ramsey_s"],
                transmon_tphi_s=scenario["transmon_tphi_s"],
            )
            replay = simulate_two_tone_case(
                model,
                frame,
                n=int(selected["n"]),
                target_coupling_mhz=float(selected["target_coupling_MHz"]),
                common_detuning_mhz=float(selected["common_detuning_MHz"]),
                noise=noise,
                dt_s=DEFAULT_FINAL_DT_S,
                duration_s=float(selected["simulation_duration_ns"]) * 1.0e-9,
                include_curves=False,
            )
            rows.append(
                {
                    "n": int(selected["n"]),
                    "case_role": str(selected["case_role"]),
                    "mechanism": str(selected["mechanism"]),
                    "target_coupling_MHz": float(selected["target_coupling_MHz"]),
                    "common_detuning_MHz": float(selected["common_detuning_MHz"]),
                    "noise_scenario": str(scenario["noise_scenario"]),
                    "noise_description": str(scenario["description"]),
                    "storage_amplitude_MHz": float(selected["storage_amplitude_MHz"]),
                    "readout_amplitude_MHz": float(selected["readout_amplitude_MHz"]),
                    "peak_target_probability": float(replay["peak_target_probability"]),
                    "time_to_peak_target_ns": float(replay["time_to_peak_target_ns"]),
                    "max_intermediate_probability": float(replay["max_intermediate_probability"]),
                    "leakage_at_target_peak": float(replay["leakage_at_target_peak"]),
                    "final_target_probability": float(replay["final_target_probability"]),
                    "closed_peak_target_probability": float(selected["peak_target_probability"]),
                    "closed_time_to_peak_target_ns": float(selected["time_to_peak_target_ns"]),
                    "closed_max_intermediate_probability": float(selected["max_intermediate_probability"]),
                }
            )
    return rows


def validation_rows(model, frame, selected_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    representatives = [row for row in selected_rows if int(row["n"]) in (1, 3)]
    rows = []
    for row in representatives:
        baseline = simulate_two_tone_case(
            model,
            frame,
            n=int(row["n"]),
            target_coupling_mhz=float(row["target_coupling_MHz"]),
            common_detuning_mhz=float(row["common_detuning_MHz"]),
            dt_s=DEFAULT_FINAL_DT_S,
            duration_s=float(row["simulation_duration_ns"]) * 1.0e-9,
            include_curves=False,
        )
        coarse_dt = simulate_two_tone_case(
            model,
            frame,
            n=int(row["n"]),
            target_coupling_mhz=float(row["target_coupling_MHz"]),
            common_detuning_mhz=float(row["common_detuning_MHz"]),
            dt_s=1.0e-9,
            duration_s=float(row["simulation_duration_ns"]) * 1.0e-9,
            include_curves=False,
        )
        medium_dt = simulate_two_tone_case(
            model,
            frame,
            n=int(row["n"]),
            target_coupling_mhz=float(row["target_coupling_MHz"]),
            common_detuning_mhz=float(row["common_detuning_MHz"]),
            dt_s=DEFAULT_SWEEP_DT_S,
            duration_s=float(row["simulation_duration_ns"]) * 1.0e-9,
            include_curves=False,
        )
        larger_truncation = simulate_two_tone_case(
            build_model(n_storage=6, n_readout=6, n_tr=5),
            build_frame(build_model(n_storage=6, n_readout=6, n_tr=5)),
            n=int(row["n"]),
            target_coupling_mhz=float(row["target_coupling_MHz"]),
            common_detuning_mhz=float(row["common_detuning_MHz"]),
            dt_s=DEFAULT_FINAL_DT_S,
            duration_s=float(row["simulation_duration_ns"]) * 1.0e-9,
            include_curves=False,
        )
        rows.append(
            {
                "n": int(row["n"]),
                "case_role": str(row["case_role"]),
                "baseline_peak_target_probability": float(baseline["peak_target_probability"]),
                "baseline_peak_time_ns": float(baseline["time_to_peak_target_ns"]),
                "dt_1p0ns_peak_target_delta": float(coarse_dt["peak_target_probability"] - baseline["peak_target_probability"]),
                "dt_0p5ns_peak_target_delta": float(medium_dt["peak_target_probability"] - baseline["peak_target_probability"]),
                "dt_1p0ns_peak_time_delta_ns": float(coarse_dt["time_to_peak_target_ns"] - baseline["time_to_peak_target_ns"]),
                "larger_truncation_peak_target_delta": float(larger_truncation["peak_target_probability"] - baseline["peak_target_probability"]),
                "larger_truncation_peak_time_delta_ns": float(larger_truncation["time_to_peak_target_ns"] - baseline["time_to_peak_target_ns"]),
                "reduced_model_peak_target_error": float(baseline["peak_target_error_vs_reduced"]),
                "reduced_model_peak_time_error_ns": float(baseline["peak_time_error_vs_reduced_ns"]),
            }
        )
    return rows


def make_scan_summary_figure(rows: list[dict[str, object]]) -> None:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, len(TRANSFER_N_VALUES), figsize=(13.5, 9.5), sharex=True)
    coupling_values = list(TARGET_COUPLING_GRID_MHZ)
    detuning_values = list(COMMON_DETUNING_GRID_MHZ)
    value_specs = (
        ("peak_target_probability", "Peak target probability", 0.0, 1.0, "viridis"),
        ("max_intermediate_probability", "Max intermediate probability", 0.0, 1.0, "magma_r"),
        ("time_to_peak_target_ns", "Time to peak target (ns)", None, None, "cividis"),
    )
    for column, n in enumerate(TRANSFER_N_VALUES):
        subset = [row for row in rows if int(row["n"]) == int(n)]
        for row_index, (field, title, vmin, vmax, cmap) in enumerate(value_specs):
            heat = np.zeros((len(coupling_values), len(detuning_values)))
            for entry in subset:
                i = coupling_values.index(float(entry["target_coupling_MHz"]))
                j = detuning_values.index(float(entry["common_detuning_MHz"]))
                heat[i, j] = float(entry[field])
            ax = axes[row_index, column]
            image = ax.imshow(
                heat,
                origin="lower",
                aspect="auto",
                extent=(detuning_values[0], detuning_values[-1], coupling_values[0], coupling_values[-1]),
                vmin=vmin,
                vmax=vmax,
                cmap=cmap,
            )
            ax.set_title(f"n={n}: {title}")
            if row_index == len(value_specs) - 1:
                ax.set_xlabel("Common detuning (MHz)")
            if column == 0:
                ax.set_ylabel("Equal leg coupling (MHz)")
            fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    plot_save(fig, "two_tone_scan_summary")


def make_dynamics_figure(trajectories: list[dict[str, object]]) -> None:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(len(TRANSFER_N_VALUES), 2, figsize=(12.5, 10.0), sharex=True, sharey=True)
    axis_lookup = {}
    for row_index, n in enumerate(TRANSFER_N_VALUES):
        axis_lookup[(n, "best_resonant")] = axes[row_index, 0]
        axis_lookup[(n, "best_detuned")] = axes[row_index, 1]
    for trajectory in trajectories:
        ax = axis_lookup[(int(trajectory["n"]), str(trajectory["case_role"]))]
        ax.plot(trajectory["times_ns"], trajectory["source_curve"], linewidth=1.8, label="Source")
        ax.plot(trajectory["times_ns"], trajectory["intermediate_curve"], linewidth=1.8, label="Intermediate")
        ax.plot(trajectory["times_ns"], trajectory["target_curve"], linewidth=1.8, label="Target")
        ax.plot(trajectory["times_ns"], trajectory["leakage_curve"], linewidth=1.6, label="Leakage")
        ax.plot(trajectory["times_ns"], trajectory["reduced_target_curve"], linestyle="--", linewidth=1.2, color="black", label="Reduced target")
        mechanism = "resonant" if str(trajectory["case_role"]) == "best_resonant" else "detuned"
        ax.set_title(
            f"n={trajectory['n']} {mechanism}: g={trajectory['target_coupling_MHz']:.0f} MHz, "
            f"detuning={trajectory['common_detuning_MHz']:.0f} MHz"
        )
        ax.set_ylim(0.0, 1.02)
    axes[-1, 0].set_xlabel("Time (ns)")
    axes[-1, 1].set_xlabel("Time (ns)")
    for row in range(len(TRANSFER_N_VALUES)):
        axes[row, 0].set_ylabel("Population")
    axes[0, 1].legend(loc="upper right")
    plot_save(fig, "two_tone_population_dynamics")


def make_open_system_figure(open_rows: list[dict[str, object]]) -> None:
    import matplotlib.pyplot as plt

    labels = []
    closed_peak = []
    mode_only_peak = []
    transmon_peak = []
    closed_f = []
    mode_only_f = []
    transmon_f = []
    for n in TRANSFER_N_VALUES:
        for case_role in ("best_resonant", "best_detuned"):
            subset = [row for row in open_rows if int(row["n"]) == int(n) and str(row["case_role"]) == case_role]
            subset_by_scenario = {str(row["noise_scenario"]): row for row in subset}
            representative = subset_by_scenario["mode_only"]
            labels.append(f"n={n}\n{'res' if case_role == 'best_resonant' else 'det'}")
            closed_peak.append(float(representative["closed_peak_target_probability"]))
            closed_f.append(float(representative["closed_max_intermediate_probability"]))
            mode_only_peak.append(float(subset_by_scenario["mode_only"]["peak_target_probability"]))
            mode_only_f.append(float(subset_by_scenario["mode_only"]["max_intermediate_probability"]))
            transmon_peak.append(float(subset_by_scenario["transmon_reference"]["peak_target_probability"]))
            transmon_f.append(float(subset_by_scenario["transmon_reference"]["max_intermediate_probability"]))

    x = np.arange(len(labels), dtype=float)
    width = 0.24
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6), sharex=True)
    axes[0].bar(x - width, closed_peak, width=width, label="Closed")
    axes[0].bar(x, mode_only_peak, width=width, label="Mode-only noise")
    axes[0].bar(x + width, transmon_peak, width=width, label="Transmon reference")
    axes[0].set_ylabel("Peak target probability")
    axes[0].set_ylim(0.0, 1.02)
    axes[0].set_title("Transfer efficiency under noise")
    axes[0].legend()

    axes[1].bar(x - width, closed_f, width=width, label="Closed")
    axes[1].bar(x, mode_only_f, width=width, label="Mode-only noise")
    axes[1].bar(x + width, transmon_f, width=width, label="Transmon reference")
    axes[1].set_ylabel("Max intermediate probability")
    axes[1].set_ylim(0.0, 1.02)
    axes[1].set_title("Intermediate-state cost under noise")

    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
    plot_save(fig, "two_tone_open_system_summary")


def save_trajectory_csvs(trajectories: list[dict[str, object]]) -> None:
    for trajectory in trajectories:
        csv_dump(
            DATA_DIR / f"two_tone_trajectory_n{trajectory['n']}_{trajectory['case_role']}.csv",
            [
                {
                    "time_ns": float(trajectory["times_ns"][index]),
                    "source_probability": float(trajectory["source_curve"][index]),
                    "intermediate_probability": float(trajectory["intermediate_curve"][index]),
                    "target_probability": float(trajectory["target_curve"][index]),
                    "leakage_probability": float(trajectory["leakage_curve"][index]),
                    "reduced_target_probability": float(trajectory["reduced_target_curve"][index]),
                }
                for index in range(len(trajectory["times_ns"]))
            ],
        )


def run_two_tone_transfer_extension(model=None, frame=None) -> dict[str, object]:
    start = time.time()
    if model is None:
        model = build_model()
    if frame is None:
        frame = build_frame(model)

    frequency_table = frequency_rows(model, frame)
    sweep = scan_rows(model, frame)
    selected = select_case_rows(sweep)
    trajectories = trajectory_rows(model, frame, selected)
    noisy = open_system_rows(model, frame, selected)
    validation = validation_rows(model, frame, selected)

    csv_dump(DATA_DIR / "two_tone_frequency_table.csv", frequency_table)
    csv_dump(DATA_DIR / "two_tone_scan.csv", sweep)
    csv_dump(DATA_DIR / "two_tone_selected_cases.csv", selected)
    csv_dump(DATA_DIR / "two_tone_open_system.csv", noisy)
    csv_dump(DATA_DIR / "two_tone_validation.csv", validation)
    save_trajectory_csvs(trajectories)

    make_scan_summary_figure(sweep)
    make_dynamics_figure(trajectories)
    make_open_system_figure(noisy)

    payload = {
        "runtime_s": float(time.time() - start),
        "device": asdict(DEVICE),
        "transmon_reference": asdict(TRANSMON_REFERENCE),
        "scan_grid": {
            "n_values": list(TRANSFER_N_VALUES),
            "target_coupling_grid_MHz": TARGET_COUPLING_GRID_MHZ.tolist(),
            "common_detuning_grid_MHz": COMMON_DETUNING_GRID_MHZ.tolist(),
        },
        "frequency_table": frequency_table,
        "scan_rows": sweep,
        "selected_cases": selected,
        "trajectories": trajectories,
        "open_system_rows": noisy,
        "validation_rows": validation,
    }
    json_dump(DATA_DIR / "two_tone_results.json", payload)
    json_dump(ARTIFACTS_DIR / "two_tone_results.json", payload)
    return payload


def main() -> None:
    payload = run_two_tone_transfer_extension()
    print(
        "Two-tone extension complete: "
        f"{len(payload['scan_rows'])} closed-system scan points, "
        f"{len(payload['selected_cases'])} selected cases, "
        f"{len(payload['open_system_rows'])} open-system replays."
    )


if __name__ == "__main__":
    main()