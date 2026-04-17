"""Validation pass for the gf-sideband waveform optimization study."""

from __future__ import annotations

import csv
from pathlib import Path

from common import (
    ARTIFACTS_DIR,
    DATA_DIR,
    DEFAULT_FINAL_DT_S,
    DEFAULT_SWEEP_DT_S,
    FAMILY_VARIANTS,
    MODES,
    N_VALUES,
    analytic_rotating_sideband_frequency,
    basis_state_for_mode,
    build_frame,
    build_model,
    candidate_target_transfer,
    csv_dump,
    json_dump,
    make_pulse,
    mhz,
    sideband_drive_target,
    sideband_rotating_frequency,
    simulate_single_pulse,
)
from two_tone_transfer_extension import frequency_rows as two_tone_frequency_rows
from two_tone_transfer_extension import validation_rows as two_tone_validation_rows

RESULTS_PATH = DATA_DIR / "study_results.json"
REFINED_CASES_PATH = DATA_DIR / "refined_cases.csv"
TWO_TONE_SELECTED_CASES_PATH = DATA_DIR / "two_tone_selected_cases.csv"


def load_refined_cases() -> list[dict[str, object]]:
    rows = []
    with REFINED_CASES_PATH.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            rows.append(
                {
                    "mode": row["mode"],
                    "criterion": row["criterion"],
                    "family": row["family"],
                    "variant": row["variant"],
                    "n": int(row["n"]),
                    "duration_ns": float(row["duration_ns"]),
                    "amplitude_MHz": float(row["amplitude_MHz"]),
                    "optimal_detuning_MHz": float(row["optimal_detuning_MHz"]),
                }
            )
    return rows


def load_two_tone_selected_cases() -> list[dict[str, object]]:
    rows = []
    with TWO_TONE_SELECTED_CASES_PATH.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            rows.append(
                {
                    "n": int(row["n"]),
                    "mechanism": row["mechanism"],
                    "target_coupling_MHz": float(row["target_coupling_MHz"]),
                    "common_detuning_MHz": float(row["common_detuning_MHz"]),
                    "simulation_duration_ns": float(row["simulation_duration_ns"]),
                    "case_role": row["case_role"],
                    "selection_status": row["selection_status"],
                }
            )
    return rows


def find_case(refined_cases: list[dict[str, object]], *, mode: str, criterion: str, n: int) -> dict[str, object]:
    return next(row for row in refined_cases if row["mode"] == mode and row["criterion"] == criterion and int(row["n"]) == int(n))


def simulate_case(case: dict[str, object], *, dt_s: float, n_storage: int = 5, n_readout: int = 5, n_tr: int = 4) -> dict[str, float]:
    model = build_model(n_storage=n_storage, n_readout=n_readout, n_tr=n_tr)
    frame = build_frame(model)
    variant = next(candidate for candidate in FAMILY_VARIANTS if candidate.family == case["family"] and candidate.variant == case["variant"])
    mode = str(case["mode"])
    n = int(case["n"])
    duration_s = float(case["duration_ns"]) * 1.0e-9
    pulse = make_pulse(
        channel=f"{mode}_sb",
        carrier_rad_s=sideband_rotating_frequency(model, frame, mode, n) + 2.0 * 3.141592653589793 * float(case["optimal_detuning_MHz"]) * 1.0e6,
        duration_s=duration_s,
        amplitude_hz=float(case["amplitude_MHz"]) * 1.0e6,
        variant=variant,
        label=f"validate_{mode}_{case['criterion']}_n{n}",
    )
    source, target = basis_state_for_mode(model, mode, n)
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
    return candidate_target_transfer(result.final_state, source, target)


def analytic_frequency_rows() -> list[dict[str, object]]:
    model = build_model()
    frame = build_frame(model)
    rows = []
    for mode in MODES:
        for n in N_VALUES:
            exact_mhz = mhz(sideband_rotating_frequency(model, frame, mode, n))
            analytic_mhz = analytic_rotating_sideband_frequency(mode, n) / 1.0e6
            rows.append(
                {
                    "mode": mode,
                    "n": n,
                    "exact_rotating_frequency_MHz": float(exact_mhz),
                    "analytic_rotating_frequency_MHz": float(analytic_mhz),
                    "difference_kHz": float((exact_mhz - analytic_mhz) * 1.0e3),
                }
            )
    return rows


def convergence_rows(refined_cases: list[dict[str, object]]) -> list[dict[str, object]]:
    representatives = [
        find_case(refined_cases, mode="storage", criterion="selective", n=2),
        find_case(refined_cases, mode="storage", criterion="unselective", n=1),
        find_case(refined_cases, mode="readout", criterion="selective", n=1),
        find_case(refined_cases, mode="readout", criterion="unselective", n=3),
    ]
    rows = []
    for case in representatives:
        baseline = simulate_case(case, dt_s=DEFAULT_FINAL_DT_S)
        coarse_dt = simulate_case(case, dt_s=1.0e-9)
        medium_dt = simulate_case(case, dt_s=DEFAULT_SWEEP_DT_S)
        larger_truncation = simulate_case(case, dt_s=DEFAULT_FINAL_DT_S, n_storage=6, n_readout=6, n_tr=5)
        rows.append(
            {
                "mode": case["mode"],
                "criterion": case["criterion"],
                "n": case["n"],
                "family": case["family"],
                "baseline_target_probability": float(baseline["target_probability"]),
                "dt_1p0ns_delta": float(coarse_dt["target_probability"] - baseline["target_probability"]),
                "dt_0p5ns_delta": float(medium_dt["target_probability"] - baseline["target_probability"]),
                "larger_truncation_delta": float(larger_truncation["target_probability"] - baseline["target_probability"]),
                "baseline_leakage_probability": float(baseline["leakage_probability"]),
            }
        )
    return rows


def literature_alignment_note() -> dict[str, str]:
    return {
        "selective_regime": "The selective winners are smooth Gaussian-family pulses, consistent with the usual sideband-control expectation that smoother edges suppress off-resonant spectral weight.",
        "fast_regime": "The fastest winners are square pulses, consistent with the expectation that hard-edged pulses maximize transfer rate when selectivity is not required.",
        "device_specific_caveat": "The readout-mode sideband becomes strongly limited by the available readout decay channel once long selective pulses are used, so the closed-system family ranking does not translate directly into practical selective readout control.",
    }


def main() -> None:
    refined_cases = load_refined_cases()
    frequency_rows = analytic_frequency_rows()
    convergence = convergence_rows(refined_cases)
    model = build_model()
    frame = build_frame(model)
    two_tone_frequency = two_tone_frequency_rows(model, frame)
    two_tone_selected_cases = load_two_tone_selected_cases()
    two_tone_convergence = two_tone_validation_rows(model, frame, two_tone_selected_cases)
    payload = {
        "analytic_frequency_validation": frequency_rows,
        "convergence_summary": convergence,
        "two_tone_frequency_validation": two_tone_frequency,
        "two_tone_convergence_summary": two_tone_convergence,
        "literature_alignment": literature_alignment_note(),
    }
    csv_dump(DATA_DIR / "validation_frequency.csv", frequency_rows)
    csv_dump(DATA_DIR / "validation_convergence.csv", convergence)
    csv_dump(DATA_DIR / "validation_two_tone_frequency.csv", two_tone_frequency)
    csv_dump(DATA_DIR / "validation_two_tone_convergence.csv", two_tone_convergence)
    json_dump(ARTIFACTS_DIR / "validation_summary.json", payload)
    print(
        "Validation written to data/validation_frequency.csv, data/validation_convergence.csv, "
        "data/validation_two_tone_frequency.csv, data/validation_two_tone_convergence.csv, and artifacts/validation_summary.json"
    )


if __name__ == "__main__":
    main()
