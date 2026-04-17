"""Run the gf-sideband waveform optimization study end to end."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import asdict
import math
import time

import numpy as np

from common import (
    ARTIFACTS_DIR,
    DATA_DIR,
    DEFAULT_FINAL_DT_S,
    DEFAULT_SWEEP_DT_S,
    DEVICE,
    FAMILY_DESCRIPTIONS,
    FAMILY_ORDER,
    FAMILY_VARIANTS,
    MODES,
    N_VALUES,
    TRANSMON_REFERENCE,
    analytic_rotating_sideband_frequency,
    basis_label,
    basis_state_for_mode,
    build_frame,
    build_model,
    build_noise,
    candidate_target_transfer,
    csv_dump,
    device_parameter_rows,
    expected_pi_amplitude_hz,
    expected_pi_duration_s,
    ghz,
    json_dump,
    make_pulse,
    mhz,
    plot_save,
    projected_swap_metrics,
    pulse_timeseries,
    rank_cases,
    sideband_drive_target,
    sideband_lab_frequency,
    sideband_matrix_element,
    sideband_rotating_frequency,
    simulate_single_pulse,
    sorted_basis_populations,
    spectrum_magnitude,
    state_population,
    write_device_manifest,
)
from two_tone_transfer_extension import run_two_tone_transfer_extension

DURATION_GRID_NS = np.array([12.0, 16.0, 20.0, 24.0, 32.0, 40.0, 60.0, 80.0, 120.0, 160.0, 220.0, 300.0])
AMPLITUDE_SCALE_GRID = np.array([0.65, 0.80, 0.95, 1.05, 1.20, 1.35, 1.50])
FINAL_DETUNING_GRID_MHZ = np.array([-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5])
ROBUSTNESS_AMPLITUDE_GRID = np.linspace(0.90, 1.10, 9)
ROBUSTNESS_DETUNING_GRID_MHZ = np.linspace(-1.5, 1.5, 13)
TRANSMON_REFERENCE_SCALE_FACTORS = (0.5, 1.0, 2.0)

UNSELECTIVE_THRESHOLDS = {
    "target_probability": 0.985,
    "leakage_probability": 0.030,
}
SELECTIVE_THRESHOLDS = {
    "target_probability": 0.990,
    "leakage_probability": 0.020,
    "max_neighbor_transfer": 0.010,
}
GATE_THRESHOLDS = {
    "target_probability": 0.985,
    "leakage_probability": 0.030,
    "max_neighbor_transfer": 0.020,
}


def mode_frequency_spacing_mhz(device=DEVICE) -> dict[str, float]:
    return {
        "storage": abs(2.0 * device.chi_storage_hz - device.storage_kerr_hz) / 1.0e6,
        "readout": abs(2.0 * device.chi_readout_hz - device.readout_kerr_hz) / 1.0e6,
    }


def transmon_noise_scenarios() -> list[dict[str, object]]:
    scenarios = [
        {
            "noise_scenario": "mode_only",
            "description": "Storage/readout dissipation only; no transmon relaxation or dephasing.",
            "transmon_t1_s": None,
            "transmon_t2_ramsey_s": None,
            "transmon_tphi_s": None,
            "scale_factor": None,
            "source": str(TRANSMON_REFERENCE.parameter_source),
        }
    ]
    for scale in TRANSMON_REFERENCE_SCALE_FACTORS:
        t1_s = None if TRANSMON_REFERENCE.qubit_t1_s is None else float(TRANSMON_REFERENCE.qubit_t1_s) * float(scale)
        t2_s = None if TRANSMON_REFERENCE.qubit_t2_ramsey_s is None else float(TRANSMON_REFERENCE.qubit_t2_ramsey_s) * float(scale)
        tphi_s = None if TRANSMON_REFERENCE.qubit_tphi_ramsey_s is None else float(TRANSMON_REFERENCE.qubit_tphi_ramsey_s) * float(scale)
        label = "transmon_reference" if math.isclose(float(scale), 1.0) else f"transmon_reference_x{scale:.1f}"
        scenarios.append(
            {
                "noise_scenario": label,
                "description": (
                    "Storage/readout dissipation plus transmon relaxation and Ramsey-derived pure dephasing "
                    f"scaled by {scale:.1f} relative to the matched local tomography reference."
                ),
                "transmon_t1_s": t1_s,
                "transmon_t2_ramsey_s": t2_s,
                "transmon_tphi_s": tphi_s,
                "scale_factor": float(scale),
                "source": str(TRANSMON_REFERENCE.parameter_source),
            }
        )
    return scenarios


def build_frequency_table(model, frame) -> list[dict[str, object]]:
    rows = []
    for mode in MODES:
        for n in N_VALUES:
            row = {
                "mode": mode,
                "n": n,
                "source_state": basis_label(
                    2,
                    int(n) - 1 if mode == "storage" else 0,
                    0 if mode == "storage" else int(n) - 1,
                ),
                "target_state": basis_label(
                    0,
                    int(n) if mode == "storage" else 0,
                    0 if mode == "storage" else int(n),
                ),
                "exact_lab_frequency_GHz": round(ghz(sideband_lab_frequency(model, mode, n)), 9),
                "exact_rotating_frequency_MHz": round(mhz(sideband_rotating_frequency(model, frame, mode, n)), 6),
                "analytic_rotating_frequency_MHz": round(analytic_rotating_sideband_frequency(mode, n) / 1.0e6, 6),
                "rotating_frequency_difference_kHz": round(
                    (
                        sideband_rotating_frequency(model, frame, mode, n) / (2.0 * np.pi)
                        - analytic_rotating_sideband_frequency(mode, n)
                    )
                    / 1.0e3,
                    9,
                ),
                "matrix_element": round(sideband_matrix_element(model, mode, n), 6),
                "expected_sqrt_n": round(math.sqrt(float(n)), 6),
                "nearest_same_mode_spacing_MHz": round(mode_frequency_spacing_mhz()[mode], 6),
                "expected_pi_time_ns_at_8MHz": round(
                    expected_pi_duration_s(amplitude_hz=8.0e6, matrix_element=sideband_matrix_element(model, mode, n))
                    * 1.0e9,
                    6,
                ),
            }
            rows.append(row)
    return rows


def simulate_candidate(model, frame, *, mode: str, n: int, variant, duration_ns: float, amplitude_scale: float, dt_s: float):
    matrix_element = sideband_matrix_element(model, mode, n)
    duration_s = float(duration_ns) * 1.0e-9
    guess_hz = expected_pi_amplitude_hz(duration_s, matrix_element)
    actual_amp_hz = guess_hz * float(amplitude_scale)
    source, target = basis_state_for_mode(model, mode, n)
    pulse = make_pulse(
        channel=f"{mode}_sb",
        carrier_rad_s=sideband_rotating_frequency(model, frame, mode, n),
        duration_s=duration_s,
        amplitude_hz=actual_amp_hz,
        variant=variant,
        label=f"{mode}_{variant.family}_{variant.variant}_n{n}_{duration_ns:.1f}ns_{amplitude_scale:.2f}",
    )
    _, result = simulate_single_pulse(
        model,
        source,
        pulse=pulse,
        duration_s=duration_s,
        drive_target=sideband_drive_target(mode),
        frame=frame,
        dt_s=dt_s,
        noise=None,
        store_states=False,
    )
    metrics = candidate_target_transfer(result.final_state, source, target)
    return {
        "mode": mode,
        "n": n,
        "family": variant.family,
        "variant": variant.variant,
        "description": variant.description,
        "duration_ns": float(duration_ns),
        "amplitude_scale": float(amplitude_scale),
        "expected_pi_amplitude_MHz": float(guess_hz / 1.0e6),
        "amplitude_MHz": float(actual_amp_hz / 1.0e6),
        "carrier_rotating_MHz": float(mhz(sideband_rotating_frequency(model, frame, mode, n))),
        **metrics,
    }


def coarse_sweep(model, frame) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for mode in MODES:
        for n in N_VALUES:
            for variant in FAMILY_VARIANTS:
                for duration_ns in DURATION_GRID_NS:
                    for amplitude_scale in AMPLITUDE_SCALE_GRID:
                        rows.append(
                            simulate_candidate(
                                model,
                                frame,
                                mode=mode,
                                n=n,
                                variant=variant,
                                duration_ns=float(duration_ns),
                                amplitude_scale=float(amplitude_scale),
                                dt_s=DEFAULT_SWEEP_DT_S,
                            )
                        )
    return rows


def shortlist_by_duration(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, int, str, float], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["mode"]), int(row["n"]), str(row["family"]), float(row["duration_ns"]))].append(row)
    shortlist = [rank_cases(group)[0] for group in grouped.values()]
    shortlist.sort(key=lambda row: (str(row["mode"]), int(row["n"]), FAMILY_ORDER.index(str(row["family"])), float(row["duration_ns"])))
    return shortlist


def neighbor_transfer_metrics(
    model,
    frame,
    *,
    pulse,
    duration_s: float,
    mode: str,
    target_n: int,
    noise=None,
    dt_s: float = DEFAULT_SWEEP_DT_S,
) -> dict[str, object]:
    rows = []
    for n_other in N_VALUES:
        if int(n_other) == int(target_n):
            continue
        source, target = basis_state_for_mode(model, mode, n_other)
        _, result = simulate_single_pulse(
            model,
            source,
            pulse=pulse,
            duration_s=duration_s,
            drive_target=sideband_drive_target(mode),
            frame=frame,
            dt_s=dt_s,
            noise=noise,
            store_states=False,
        )
        rows.append({"n_other": int(n_other), "target_probability": float(state_population(result.final_state, target))})
    return {
        "neighbor_transfer_rows": rows,
        "max_neighbor_transfer": float(max((float(row["target_probability"]) for row in rows), default=0.0)),
    }


def enrich_shortlist(model, frame, shortlist: list[dict[str, object]]) -> list[dict[str, object]]:
    enriched = []
    for row in shortlist:
        variant = next(candidate for candidate in FAMILY_VARIANTS if candidate.family == row["family"] and candidate.variant == row["variant"])
        mode = str(row["mode"])
        n = int(row["n"])
        duration_s = float(row["duration_ns"]) * 1.0e-9
        pulse = make_pulse(
            channel=f"{mode}_sb",
            carrier_rad_s=sideband_rotating_frequency(model, frame, mode, n),
            duration_s=duration_s,
            amplitude_hz=float(row["amplitude_MHz"]) * 1.0e6,
            variant=variant,
            label=f"enriched_{mode}_{row['family']}_{row['variant']}_n{n}",
        )
        source, _ = basis_state_for_mode(model, mode, n)
        _, source_result = simulate_single_pulse(
            model,
            source,
            pulse=pulse,
            duration_s=duration_s,
            drive_target=sideband_drive_target(mode),
            frame=frame,
            dt_s=DEFAULT_SWEEP_DT_S,
            noise=None,
            store_states=False,
        )
        projected = projected_swap_metrics(
            model,
            frame,
            mode=mode,
            n=n,
            pulse=pulse,
            duration_s=duration_s,
            dt_s=DEFAULT_SWEEP_DT_S,
        )
        neighbor = neighbor_transfer_metrics(model, frame, pulse=pulse, duration_s=duration_s, mode=mode, target_n=n)
        enriched_row = dict(row)
        enriched_row.update(projected)
        enriched_row.update(neighbor)
        enriched_row["selectivity_ratio"] = float(enriched_row["target_probability"]) / max(float(enriched_row["max_neighbor_transfer"]), 1.0e-6)
        enriched_row["dominant_basis_states"] = sorted_basis_populations(source_result.final_state, model, cutoff=5)
        enriched.append(enriched_row)
    return enriched


def open_system_duration_metrics(model, frame, duration_rows: list[dict[str, object]], *, scenario: dict[str, object]) -> list[dict[str, object]]:
    noise = build_noise(
        transmon_t1_s=scenario["transmon_t1_s"],
        transmon_t2_ramsey_s=scenario["transmon_t2_ramsey_s"],
        transmon_tphi_s=scenario["transmon_tphi_s"],
    )
    enriched = []
    for row in duration_rows:
        variant = next(candidate for candidate in FAMILY_VARIANTS if candidate.family == row["family"] and candidate.variant == row["variant"])
        mode = str(row["mode"])
        n = int(row["n"])
        duration_s = float(row["duration_ns"]) * 1.0e-9
        pulse = make_pulse(
            channel=f"{mode}_sb",
            carrier_rad_s=sideband_rotating_frequency(model, frame, mode, n),
            duration_s=duration_s,
            amplitude_hz=float(row["amplitude_MHz"]) * 1.0e6,
            variant=variant,
            label=f"open_system_{scenario['noise_scenario']}_{mode}_{row['family']}_{row['variant']}_n{n}",
        )
        source, target = basis_state_for_mode(model, mode, n)
        _, result = simulate_single_pulse(
            model,
            source,
            pulse=pulse,
            duration_s=duration_s,
            drive_target=sideband_drive_target(mode),
            frame=frame,
            dt_s=DEFAULT_SWEEP_DT_S,
            noise=noise,
            store_states=False,
        )
        metrics = candidate_target_transfer(result.final_state, source, target)
        neighbor = neighbor_transfer_metrics(
            model,
            frame,
            pulse=pulse,
            duration_s=duration_s,
            mode=mode,
            target_n=n,
            noise=noise,
            dt_s=DEFAULT_SWEEP_DT_S,
        )
        updated = dict(row)
        updated.update(metrics)
        updated.update(neighbor)
        updated["noise_scenario"] = str(scenario["noise_scenario"])
        updated["noise_description"] = str(scenario["description"])
        enriched.append(updated)
    return enriched


def passes_thresholds(row: dict[str, object], thresholds: dict[str, float]) -> bool:
    for key, threshold in thresholds.items():
        value = float(row[key])
        if key in {"leakage_probability", "max_neighbor_transfer"}:
            if value > float(threshold):
                return False
        else:
            if value < float(threshold):
                return False
    return True


def duration_level_table(enriched: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for row in enriched:
        new_row = dict(row)
        new_row["passes_unselective"] = passes_thresholds(new_row, UNSELECTIVE_THRESHOLDS)
        new_row["passes_selective"] = passes_thresholds(new_row, SELECTIVE_THRESHOLDS)
        rows.append(new_row)
    return rows


def family_threshold_cases(duration_rows: list[dict[str, object]], criterion: str) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, int], list[dict[str, object]]] = defaultdict(list)
    for row in duration_rows:
        if not bool(row[f"passes_{criterion}"]):
            continue
        grouped[(str(row["mode"]), str(row["family"]), int(row["n"]))].append(row)
    selected = []
    for cases in grouped.values():
        best = min(
            cases,
            key=lambda row: (
                float(row["duration_ns"]),
                -float(row["target_probability"]),
                float(row["leakage_probability"]),
                -float(row["projected_swap_fidelity"]),
                float(row["max_neighbor_transfer"]),
            ),
        )
        record = dict(best)
        record["criterion"] = criterion
        selected.append(record)
    selected.sort(key=lambda row: (str(row["mode"]), str(row["criterion"]), FAMILY_ORDER.index(str(row["family"])), int(row["n"])))
    return selected


def mode_level_winners(selected_cases: list[dict[str, object]]) -> list[dict[str, object]]:
    by_key: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in selected_cases:
        by_key[(str(row["mode"]), str(row["criterion"]), str(row["family"]))].append(row)

    winners = []
    for mode in MODES:
        for criterion in ("selective", "unselective"):
            candidates = []
            for family in FAMILY_ORDER:
                rows = by_key.get((mode, criterion, family), [])
                if {int(row["n"]) for row in rows} != set(N_VALUES):
                    continue
                rows_sorted = sorted(rows, key=lambda row: int(row["n"]))
                hardest = max(rows_sorted, key=lambda row: float(row["duration_ns"]))
                candidates.append(
                    {
                        "mode": mode,
                        "criterion": criterion,
                        "family": family,
                        "family_description": FAMILY_DESCRIPTIONS[family],
                        "conservative_duration_ns": max(float(row["duration_ns"]) for row in rows_sorted),
                        "mean_duration_ns": float(np.mean([float(row["duration_ns"]) for row in rows_sorted])),
                        "mean_target_probability": float(np.mean([float(row["target_probability"]) for row in rows_sorted])),
                        "mean_leakage_probability": float(np.mean([float(row["leakage_probability"]) for row in rows_sorted])),
                        "mean_max_neighbor_transfer": float(np.mean([float(row["max_neighbor_transfer"]) for row in rows_sorted])),
                        "hardest_n": int(hardest["n"]),
                        "per_n_cases": rows_sorted,
                    }
                )
            if not candidates:
                continue
            winners.append(
                min(
                    candidates,
                    key=lambda row: (
                        float(row["conservative_duration_ns"]),
                        float(row["mean_duration_ns"]),
                        -float(row["mean_target_probability"]),
                        float(row["mean_leakage_probability"]),
                        float(row["mean_max_neighbor_transfer"]),
                    ),
                )
            )
    return winners


def gate_selected_cases(duration_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, int], list[dict[str, object]]] = defaultdict(list)
    for row in duration_rows:
        if not passes_thresholds(row, GATE_THRESHOLDS):
            continue
        grouped[(str(row["mode"]), str(row["family"]), int(row["n"]))].append(row)
    selected = []
    for cases in grouped.values():
        best = max(
            cases,
            key=lambda row: (
                float(row["projected_swap_fidelity"]),
                float(row["target_probability"]),
                -float(row["leakage_probability"]),
                -float(row["duration_ns"]),
                -abs(float(row["offdiag_phase_asymmetry_rad"])),
            ),
        )
        record = dict(best)
        record["criterion"] = "gate_oriented"
        selected.append(record)
    selected.sort(key=lambda row: (str(row["mode"]), FAMILY_ORDER.index(str(row["family"])), int(row["n"])))
    return selected


def gate_family_summary_rows(selected_cases: list[dict[str, object]]) -> list[dict[str, object]]:
    by_key: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in selected_cases:
        by_key[(str(row["mode"]), str(row["family"]))].append(row)

    winners = []
    for mode in MODES:
        candidates = []
        for family in FAMILY_ORDER:
            rows = by_key.get((mode, family), [])
            if {int(row["n"]) for row in rows} != set(N_VALUES):
                continue
            rows_sorted = sorted(rows, key=lambda row: int(row["n"]))
            hardest = max(rows_sorted, key=lambda row: float(row["duration_ns"]))
            candidates.append(
                {
                    "mode": mode,
                    "family": family,
                    "family_description": FAMILY_DESCRIPTIONS[family],
                    "conservative_duration_ns": max(float(row["duration_ns"]) for row in rows_sorted),
                    "mean_duration_ns": float(np.mean([float(row["duration_ns"]) for row in rows_sorted])),
                    "mean_target_probability": float(np.mean([float(row["target_probability"]) for row in rows_sorted])),
                    "mean_leakage_probability": float(np.mean([float(row["leakage_probability"]) for row in rows_sorted])),
                    "mean_max_neighbor_transfer": float(np.mean([float(row["max_neighbor_transfer"]) for row in rows_sorted])),
                    "mean_projected_swap_fidelity": float(np.mean([float(row["projected_swap_fidelity"]) for row in rows_sorted])),
                    "mean_abs_phase_asymmetry_rad": float(np.mean([abs(float(row["offdiag_phase_asymmetry_rad"])) for row in rows_sorted])),
                    "hardest_n": int(hardest["n"]),
                    "per_n_cases": rows_sorted,
                }
            )
        if candidates:
            winners.append(
                max(
                    candidates,
                    key=lambda row: (
                        float(row["mean_projected_swap_fidelity"]),
                        float(row["mean_target_probability"]),
                        -float(row["mean_leakage_probability"]),
                        -float(row["conservative_duration_ns"]),
                    ),
                )
            )
    return winners


def refine_detuning(model, frame, *, case: dict[str, object], criterion: str) -> dict[str, object]:
    variant = next(candidate for candidate in FAMILY_VARIANTS if candidate.family == case["family"] and candidate.variant == case["variant"])
    mode = str(case["mode"])
    n = int(case["n"])
    duration_s = float(case["duration_ns"]) * 1.0e-9
    base_frequency = sideband_rotating_frequency(model, frame, mode, n)
    source, target = basis_state_for_mode(model, mode, n)
    best = None
    for detuning_mhz in FINAL_DETUNING_GRID_MHZ:
        pulse = make_pulse(
            channel=f"{mode}_sb",
            carrier_rad_s=base_frequency + 2.0 * np.pi * float(detuning_mhz) * 1.0e6,
            duration_s=duration_s,
            amplitude_hz=float(case["amplitude_MHz"]) * 1.0e6,
            variant=variant,
            label=f"refine_{mode}_{criterion}_n{n}",
        )
        _, result = simulate_single_pulse(
            model,
            source,
            pulse=pulse,
            duration_s=duration_s,
            drive_target=sideband_drive_target(mode),
            frame=frame,
            dt_s=DEFAULT_FINAL_DT_S,
            noise=None,
            store_states=False,
        )
        metrics = candidate_target_transfer(result.final_state, source, target)
        current = {"detuning_MHz": float(detuning_mhz), **metrics}
        if best is None or (
            current["target_probability"],
            -current["leakage_probability"],
            -abs(current["detuning_MHz"]),
        ) > (
            best["target_probability"],
            -best["leakage_probability"],
            -abs(best["detuning_MHz"]),
        ):
            best = current
    assert best is not None
    updated = dict(case)
    updated["criterion"] = criterion
    updated["optimal_detuning_MHz"] = float(best["detuning_MHz"])
    updated["refined_target_probability"] = float(best["target_probability"])
    updated["refined_leakage_probability"] = float(best["leakage_probability"])
    return updated


def open_system_followup(model, frame, noise, *, case: dict[str, object], scenario: dict[str, object]) -> dict[str, object]:
    variant = next(candidate for candidate in FAMILY_VARIANTS if candidate.family == case["family"] and candidate.variant == case["variant"])
    mode = str(case["mode"])
    n = int(case["n"])
    duration_s = float(case["duration_ns"]) * 1.0e-9
    source, target = basis_state_for_mode(model, mode, n)
    pulse = make_pulse(
        channel=f"{mode}_sb",
        carrier_rad_s=sideband_rotating_frequency(model, frame, mode, n) + 2.0 * np.pi * float(case.get("optimal_detuning_MHz", 0.0)) * 1.0e6,
        duration_s=duration_s,
        amplitude_hz=float(case["amplitude_MHz"]) * 1.0e6,
        variant=variant,
        label=f"noisy_{mode}_{case['criterion']}_n{n}",
    )
    _, result = simulate_single_pulse(
        model,
        source,
        pulse=pulse,
        duration_s=duration_s,
        drive_target=sideband_drive_target(mode),
        frame=frame,
        dt_s=DEFAULT_FINAL_DT_S,
        noise=noise,
        store_states=False,
    )
    metrics = candidate_target_transfer(result.final_state, source, target)
    return {
        "mode": mode,
        "criterion": str(case["criterion"]),
        "family": str(case["family"]),
        "n": n,
        "noise_scenario": str(scenario["noise_scenario"]),
        "noise_description": str(scenario["description"]),
        "transmon_t1_us": None if scenario["transmon_t1_s"] is None else float(scenario["transmon_t1_s"]) * 1.0e6,
        "transmon_t2_ramsey_us": None if scenario["transmon_t2_ramsey_s"] is None else float(scenario["transmon_t2_ramsey_s"]) * 1.0e6,
        "transmon_tphi_us": None if scenario["transmon_tphi_s"] is None else float(scenario["transmon_tphi_s"]) * 1.0e6,
        "noise_source": str(scenario["source"]),
        "duration_ns": float(case["duration_ns"]),
        "amplitude_MHz": float(case["amplitude_MHz"]),
        "optimal_detuning_MHz": float(case.get("optimal_detuning_MHz", 0.0)),
        "noisy_target_probability": float(metrics["target_probability"]),
        "noisy_leakage_probability": float(metrics["leakage_probability"]),
        "closed_target_probability": float(case.get("refined_target_probability", case["target_probability"])),
        "closed_leakage_probability": float(case.get("refined_leakage_probability", case["leakage_probability"])),
    }


def open_system_summary_rows(noisy_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in noisy_rows:
        grouped[(str(row["noise_scenario"]), str(row["mode"]), str(row["criterion"]))].append(row)
    summary = []
    for (noise_scenario, mode, criterion), rows in sorted(grouped.items()):
        summary.append(
            {
                "noise_scenario": noise_scenario,
                "mode": mode,
                "criterion": criterion,
                "family": str(rows[0]["family"]),
                "mean_noisy_target_probability": float(np.mean([float(row["noisy_target_probability"]) for row in rows])),
                "mean_noisy_leakage_probability": float(np.mean([float(row["noisy_leakage_probability"]) for row in rows])),
                "worst_noisy_target_probability": float(min(float(row["noisy_target_probability"]) for row in rows)),
                "worst_noisy_leakage_probability": float(max(float(row["noisy_leakage_probability"]) for row in rows)),
            }
        )
    return summary


def robustness_map(model, frame, *, case: dict[str, object]) -> dict[str, object]:
    variant = next(candidate for candidate in FAMILY_VARIANTS if candidate.family == case["family"] and candidate.variant == case["variant"])
    mode = str(case["mode"])
    n = int(case["n"])
    duration_s = float(case["duration_ns"]) * 1.0e-9
    source, target = basis_state_for_mode(model, mode, n)
    base_frequency = sideband_rotating_frequency(model, frame, mode, n) + 2.0 * np.pi * float(case.get("optimal_detuning_MHz", 0.0)) * 1.0e6
    rows = []
    for amplitude_scale in ROBUSTNESS_AMPLITUDE_GRID:
        for detuning_mhz in ROBUSTNESS_DETUNING_GRID_MHZ:
            pulse = make_pulse(
                channel=f"{mode}_sb",
                carrier_rad_s=base_frequency + 2.0 * np.pi * float(detuning_mhz) * 1.0e6,
                duration_s=duration_s,
                amplitude_hz=float(case["amplitude_MHz"]) * 1.0e6 * float(amplitude_scale),
                variant=variant,
                label=f"robust_{mode}_{case['criterion']}_n{n}",
            )
            _, result = simulate_single_pulse(
                model,
                source,
                pulse=pulse,
                duration_s=duration_s,
                drive_target=sideband_drive_target(mode),
                frame=frame,
                dt_s=DEFAULT_FINAL_DT_S,
                noise=None,
                store_states=False,
            )
            metrics = candidate_target_transfer(result.final_state, source, target)
            rows.append(
                {
                    "mode": mode,
                    "criterion": str(case["criterion"]),
                    "n": n,
                    "amplitude_scale": float(amplitude_scale),
                    "detuning_MHz": float(detuning_mhz),
                    "target_probability": float(metrics["target_probability"]),
                    "leakage_probability": float(metrics["leakage_probability"]),
                }
            )
    return {"mode": mode, "criterion": str(case["criterion"]), "family": str(case["family"]), "n": n, "rows": rows}


def trajectory_data(model, frame, *, case: dict[str, object]) -> dict[str, object]:
    variant = next(candidate for candidate in FAMILY_VARIANTS if candidate.family == case["family"] and candidate.variant == case["variant"])
    mode = str(case["mode"])
    n = int(case["n"])
    duration_s = float(case["duration_ns"]) * 1.0e-9
    source, target = basis_state_for_mode(model, mode, n)
    pulse = make_pulse(
        channel=f"{mode}_sb",
        carrier_rad_s=sideband_rotating_frequency(model, frame, mode, n) + 2.0 * np.pi * float(case.get("optimal_detuning_MHz", 0.0)) * 1.0e6,
        duration_s=duration_s,
        amplitude_hz=float(case["amplitude_MHz"]) * 1.0e6,
        variant=variant,
        label=f"traj_{mode}_{case['criterion']}_n{n}",
    )
    compiled, result = simulate_single_pulse(
        model,
        source,
        pulse=pulse,
        duration_s=duration_s,
        drive_target=sideband_drive_target(mode),
        frame=frame,
        dt_s=DEFAULT_FINAL_DT_S,
        noise=None,
        store_states=True,
    )
    times_ns = list((compiled.tlist * 1.0e9).tolist())
    source_curve = [float(state_population(state, source)) for state in result.states]
    target_curve = [float(state_population(state, target)) for state in result.states]
    leakage_curve = [max(0.0, 1.0 - source_curve[i] - target_curve[i]) for i in range(len(times_ns))]
    return {
        "mode": mode,
        "criterion": str(case["criterion"]),
        "family": str(case["family"]),
        "n": n,
        "times_ns": times_ns,
        "source_curve": source_curve,
        "target_curve": target_curve,
        "leakage_curve": leakage_curve,
        "final_top_states": sorted_basis_populations(result.final_state, model, cutoff=6),
    }


def winner_summary_rows(winners: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for winner in winners:
        per_n = {int(row["n"]): row for row in winner["per_n_cases"]}
        rows.append(
            {
                "mode": winner["mode"],
                "criterion": winner["criterion"],
                "winning_family": winner["family"],
                "conservative_duration_ns": float(winner["conservative_duration_ns"]),
                "mean_duration_ns": float(winner["mean_duration_ns"]),
                "hardest_n": int(winner["hardest_n"]),
                "n1_duration_ns": float(per_n[1]["duration_ns"]),
                "n2_duration_ns": float(per_n[2]["duration_ns"]),
                "n3_duration_ns": float(per_n[3]["duration_ns"]),
                "mean_target_probability": float(winner["mean_target_probability"]),
                "mean_leakage_probability": float(winner["mean_leakage_probability"]),
                "mean_max_neighbor_transfer": float(winner["mean_max_neighbor_transfer"]),
            }
        )
    return rows


def gate_winner_summary_rows(winners: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for winner in winners:
        per_n = {int(row["n"]): row for row in winner["per_n_cases"]}
        rows.append(
            {
                "mode": winner["mode"],
                "winning_family": winner["family"],
                "conservative_duration_ns": float(winner["conservative_duration_ns"]),
                "mean_duration_ns": float(winner["mean_duration_ns"]),
                "hardest_n": int(winner["hardest_n"]),
                "n1_projected_swap_fidelity": float(per_n[1]["projected_swap_fidelity"]),
                "n2_projected_swap_fidelity": float(per_n[2]["projected_swap_fidelity"]),
                "n3_projected_swap_fidelity": float(per_n[3]["projected_swap_fidelity"]),
                "mean_projected_swap_fidelity": float(winner["mean_projected_swap_fidelity"]),
                "mean_target_probability": float(winner["mean_target_probability"]),
                "mean_abs_phase_asymmetry_rad": float(winner["mean_abs_phase_asymmetry_rad"]),
            }
        )
    return rows


def make_figures(
    frequency_rows: list[dict[str, object]],
    duration_rows: list[dict[str, object]],
    winners: list[dict[str, object]],
    open_system_summary: list[dict[str, object]],
    robustness_results: list[dict[str, object]],
    trajectory_results: list[dict[str, object]],
) -> None:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.0))
    for ax, mode in zip(axes, MODES, strict=True):
        subset = [row for row in frequency_rows if row["mode"] == mode]
        ax.plot([row["n"] for row in subset], [row["exact_lab_frequency_GHz"] for row in subset], marker="o", linewidth=2.0)
        ax.set_title(f"{mode.capitalize()} sideband frequency map")
        ax.set_xlabel("Target bosonic manifold n")
        ax.set_ylabel("Lab-frame frequency (GHz)")
    plot_save(fig, "sideband_frequency_map")

    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.0), sharex=True, sharey="row")
    for column, mode in enumerate(MODES):
        n_subset = [row for row in duration_rows if row["mode"] == mode and int(row["n"]) == 1]
        for family in FAMILY_ORDER:
            family_rows = [row for row in n_subset if row["family"] == family]
            family_rows.sort(key=lambda row: float(row["duration_ns"]))
            axes[0, column].plot(
                [row["duration_ns"] for row in family_rows],
                [row["target_probability"] for row in family_rows],
                marker="o",
                linewidth=1.4,
                label=family,
            )
            axes[1, column].plot(
                [row["duration_ns"] for row in family_rows],
                [row["max_neighbor_transfer"] for row in family_rows],
                marker="o",
                linewidth=1.4,
                label=family,
            )
        axes[0, column].set_title(f"{mode.capitalize()} n=1 transfer")
        axes[1, column].set_title(f"{mode.capitalize()} n=1 neighbor response")
        axes[1, column].set_xlabel("Duration (ns)")
    axes[0, 0].set_ylabel("Target probability")
    axes[1, 0].set_ylabel("Max neighboring-manifold transfer")
    axes[0, 1].legend(loc="center left", bbox_to_anchor=(1.02, 0.5))
    plot_save(fig, "duration_tradeoff_n1")

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2), sharey=True)
    for ax, mode in zip(axes, MODES, strict=True):
        for criterion, marker in (("selective", "o"), ("unselective", "s")):
            winner = next((row for row in winners if row["mode"] == mode and row["criterion"] == criterion), None)
            if winner is None:
                continue
            x_values = [float(case["duration_ns"]) for case in winner["per_n_cases"]]
            y_values = [float(case["max_neighbor_transfer"]) for case in winner["per_n_cases"]]
            ax.plot(x_values, y_values, marker=marker, linewidth=2.0, label=f"{criterion}: {winner['family']}")
            for case in winner["per_n_cases"]:
                ax.annotate(f"n={int(case['n'])}", (float(case["duration_ns"]), float(case["max_neighbor_transfer"])), textcoords="offset points", xytext=(4, 4), fontsize=8)
        ax.set_title(f"{mode.capitalize()} winner frontier")
        ax.set_xlabel("Duration (ns)")
        ax.set_ylabel("Max neighboring-manifold transfer")
        ax.set_yscale("log")
        ax.legend()
    plot_save(fig, "winner_frontier")

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.5))
    for ax, result in zip(axes.flat, trajectory_results, strict=False):
        ax.plot(result["times_ns"], result["source_curve"], linewidth=1.8, label="Source")
        ax.plot(result["times_ns"], result["target_curve"], linewidth=1.8, label="Target")
        ax.plot(result["times_ns"], result["leakage_curve"], linewidth=1.8, label="Leakage")
        ax.set_title(f"{result['mode'].capitalize()} {result['criterion']} (n={result['n']}, {result['family']})")
        ax.set_xlabel("Time (ns)")
        ax.set_ylabel("Population")
        ax.set_ylim(0.0, 1.02)
        ax.legend()
    plot_save(fig, "population_dynamics_representative")

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.5), sharex=True, sharey=True)
    for ax, result in zip(axes.flat, robustness_results, strict=False):
        rows = result["rows"]
        amp_vals = sorted({float(row["amplitude_scale"]) for row in rows})
        det_vals = sorted({float(row["detuning_MHz"]) for row in rows})
        heat = np.zeros((len(amp_vals), len(det_vals)))
        for row in rows:
            i = amp_vals.index(float(row["amplitude_scale"]))
            j = det_vals.index(float(row["detuning_MHz"]))
            heat[i, j] = float(row["target_probability"])
        im = ax.imshow(
            heat,
            aspect="auto",
            origin="lower",
            extent=(det_vals[0], det_vals[-1], amp_vals[0], amp_vals[-1]),
            vmin=0.0,
            vmax=1.0,
            cmap="magma",
        )
        ax.set_title(f"{result['mode'].capitalize()} {result['criterion']} (n={result['n']})")
        ax.set_xlabel("Detuning error (MHz)")
        ax.set_ylabel("Amplitude scale")
    fig.colorbar(im, ax=axes.ravel().tolist(), label="Target probability")
    plot_save(fig, "robustness_maps")

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.0), sharey=True)
    compare_variants = {
        "square": next(variant for variant in FAMILY_VARIANTS if variant.family == "square"),
        "flat_top_cosine": next(variant for variant in FAMILY_VARIANTS if variant.family == "flat_top_cosine" and variant.variant == "ramp_0p25"),
        "smooth_bump": next(variant for variant in FAMILY_VARIANTS if variant.family == "smooth_bump"),
        "blackman": next(variant for variant in FAMILY_VARIANTS if variant.family == "blackman"),
    }
    for family, variant in compare_variants.items():
        freqs_mhz, magnitude = spectrum_magnitude(variant, duration_s=80.0e-9)
        mask = np.abs(freqs_mhz) <= 40.0
        axes[0].plot(freqs_mhz[mask], magnitude[mask], linewidth=1.7, label=family)
        times_ns, pulse_mhz = pulse_timeseries(variant, duration_s=80.0e-9, amplitude_hz=8.0e6)
        axes[1].plot(times_ns, np.real(pulse_mhz), linewidth=1.7, label=family)
    axes[0].set_title("Spectral intuition at 80 ns")
    axes[0].set_xlabel("Baseband detuning (MHz)")
    axes[0].set_ylabel("Normalized magnitude")
    axes[0].legend()
    axes[1].set_title("Envelope comparison at 80 ns")
    axes[1].set_xlabel("Time (ns)")
    axes[1].set_ylabel("Amplitude (MHz)")
    plot_save(fig, "spectral_intuition")

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.5), sharex=True, sharey=True)
    axis_lookup = {
        ("storage", "selective"): axes[0, 0],
        ("storage", "unselective"): axes[0, 1],
        ("readout", "selective"): axes[1, 0],
        ("readout", "unselective"): axes[1, 1],
    }
    for (mode, criterion), ax in axis_lookup.items():
        subset = [row for row in open_system_summary if row["mode"] == mode and row["criterion"] == criterion]
        subset.sort(key=lambda row: str(row["noise_scenario"]))
        labels = [str(row["noise_scenario"]).replace("transmon_reference", "transmon") for row in subset]
        values = [float(row["mean_noisy_target_probability"]) for row in subset]
        ax.bar(range(len(values)), values, color="#4c78a8")
        ax.set_title(f"{mode.capitalize()} {criterion}")
        ax.set_ylabel("Mean noisy target probability")
        ax.set_ylim(0.0, 1.0)
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=25, ha="right")
    plot_save(fig, "open_system_scenario_summary")


def results_summary_payload(**kwargs) -> dict[str, object]:
    return {
        "device": asdict(DEVICE),
        "transmon_reference": asdict(TRANSMON_REFERENCE),
        "thresholds": {
            "selective": SELECTIVE_THRESHOLDS,
            "unselective": UNSELECTIVE_THRESHOLDS,
            "gate_oriented": GATE_THRESHOLDS,
        },
        **kwargs,
    }


def print_executive_summary(winners: list[dict[str, object]], refined_cases: list[dict[str, object]]) -> None:
    refined_lookup = {(str(row["mode"]), str(row["criterion"]), int(row["n"])): row for row in refined_cases}
    for winner in winners:
        refined = refined_lookup[(str(winner["mode"]), str(winner["criterion"]), int(winner["hardest_n"]))]
        print(
            f"{winner['mode'].capitalize()} {winner['criterion']}: family={winner['family']}, "
            f"conservative_duration_ns={winner['conservative_duration_ns']:.1f}, "
            f"hardest_n={winner['hardest_n']}, "
            f"refined_target={refined['refined_target_probability']:.4f}, "
            f"detuning_MHz={refined['optimal_detuning_MHz']:.2f}"
        )


def main() -> None:
    start = time.time()
    model = build_model()
    frame = build_frame(model)
    parameter_rows = device_parameter_rows()
    noise_scenarios = transmon_noise_scenarios()

    write_device_manifest(ARTIFACTS_DIR / "device_manifest.json")
    csv_dump(DATA_DIR / "device_parameter_table.csv", parameter_rows)
    json_dump(ARTIFACTS_DIR / "device_parameter_table.json", parameter_rows)
    csv_dump(DATA_DIR / "noise_scenarios.csv", noise_scenarios)
    json_dump(ARTIFACTS_DIR / "noise_scenarios.json", noise_scenarios)

    frequency_rows = build_frequency_table(model, frame)
    csv_dump(DATA_DIR / "frequency_table.csv", frequency_rows)
    json_dump(ARTIFACTS_DIR / "frequency_table.json", frequency_rows)

    sweep_rows = coarse_sweep(model, frame)
    csv_dump(DATA_DIR / "coarse_sweep.csv", sweep_rows)

    shortlist = shortlist_by_duration(sweep_rows)
    csv_dump(DATA_DIR / "duration_shortlist.csv", shortlist)

    enriched = enrich_shortlist(model, frame, shortlist)
    duration_rows = duration_level_table(enriched)
    csv_dump(DATA_DIR / "duration_level_metrics.csv", duration_rows)

    selected_cases = family_threshold_cases(duration_rows, "selective") + family_threshold_cases(duration_rows, "unselective")
    csv_dump(DATA_DIR / "selected_cases.csv", selected_cases)

    winners = mode_level_winners(selected_cases)
    csv_dump(DATA_DIR / "winner_table.csv", winner_summary_rows(winners))
    gate_cases = gate_selected_cases(duration_rows)
    gate_winners = gate_family_summary_rows(gate_cases)
    csv_dump(DATA_DIR / "gate_selected_cases.csv", gate_cases)
    csv_dump(DATA_DIR / "gate_winner_table.csv", gate_winner_summary_rows(gate_winners))
    transmon_reference_scenario = next(scenario for scenario in noise_scenarios if scenario["noise_scenario"] == "transmon_reference")
    open_system_duration_rows = duration_level_table(
        open_system_duration_metrics(model, frame, duration_rows, scenario=transmon_reference_scenario)
    )
    open_system_selected_cases = family_threshold_cases(open_system_duration_rows, "selective") + family_threshold_cases(open_system_duration_rows, "unselective")
    open_system_winners = mode_level_winners(open_system_selected_cases)
    csv_dump(DATA_DIR / "open_system_transmon_reference_duration_metrics.csv", open_system_duration_rows)
    csv_dump(DATA_DIR / "open_system_transmon_reference_selected_cases.csv", open_system_selected_cases)
    csv_dump(DATA_DIR / "open_system_transmon_reference_winner_table.csv", winner_summary_rows(open_system_winners))

    refined_cases = []
    noisy_rows = []
    robustness_results = []
    trajectory_results = []
    for winner in winners:
        for case in winner["per_n_cases"]:
            refined = refine_detuning(model, frame, case=case, criterion=str(winner["criterion"]))
            refined_cases.append(refined)
            for scenario in noise_scenarios:
                scenario_noise = build_noise(
                    transmon_t1_s=scenario["transmon_t1_s"],
                    transmon_t2_ramsey_s=scenario["transmon_t2_ramsey_s"],
                    transmon_tphi_s=scenario["transmon_tphi_s"],
                )
                noisy_rows.append(open_system_followup(model, frame, scenario_noise, case=refined, scenario=scenario))
        hardest_case = max(winner["per_n_cases"], key=lambda row: float(row["duration_ns"]))
        refined_hardest = refine_detuning(model, frame, case=hardest_case, criterion=str(winner["criterion"]))
        robustness_results.append(robustness_map(model, frame, case=refined_hardest))
        trajectory_results.append(trajectory_data(model, frame, case=refined_hardest))

    csv_dump(DATA_DIR / "refined_cases.csv", refined_cases)
    csv_dump(DATA_DIR / "open_system_followup.csv", noisy_rows)
    open_system_summary = open_system_summary_rows(noisy_rows)
    csv_dump(DATA_DIR / "open_system_summary.csv", open_system_summary)
    for result in robustness_results:
        csv_dump(DATA_DIR / f"robustness_{result['mode']}_{result['criterion']}_n{result['n']}.csv", result["rows"])
    for result in trajectory_results:
        csv_dump(
            DATA_DIR / f"trajectory_{result['mode']}_{result['criterion']}_n{result['n']}.csv",
            [
                {
                    "time_ns": result["times_ns"][i],
                    "source_probability": result["source_curve"][i],
                    "target_probability": result["target_curve"][i],
                    "leakage_probability": result["leakage_curve"][i],
                }
                for i in range(len(result["times_ns"]))
            ],
        )

    two_tone_results = run_two_tone_transfer_extension(model=model, frame=frame)

    make_figures(frequency_rows, duration_rows, winners, open_system_summary, robustness_results, trajectory_results)

    runtime_s = time.time() - start
    results = results_summary_payload(
        runtime_s=runtime_s,
        device_parameter_rows=parameter_rows,
        noise_scenarios=noise_scenarios,
        frequency_table=frequency_rows,
        duration_level_rows=duration_rows,
        selected_cases=selected_cases,
        winners=winners,
        gate_selected_cases=gate_cases,
        gate_winners=gate_winners,
        open_system_transmon_reference_duration_rows=open_system_duration_rows,
        open_system_transmon_reference_selected_cases=open_system_selected_cases,
        open_system_transmon_reference_winners=open_system_winners,
        refined_cases=refined_cases,
        open_system_followup=noisy_rows,
        open_system_summary=open_system_summary,
        robustness_results=robustness_results,
        trajectory_results=trajectory_results,
        two_tone_extension=two_tone_results,
    )
    json_dump(DATA_DIR / "study_results.json", results)
    json_dump(ARTIFACTS_DIR / "study_results.json", results)
    print_executive_summary(winners, refined_cases)


if __name__ == "__main__":
    main()
