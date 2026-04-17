"""Targeted truncation and timestep convergence checks for the gf-cooling study."""

from __future__ import annotations

import json
import time
from pathlib import Path

from cqed_sim.core.drive_targets import SidebandDriveSpec

from common import (
    ARTIFACTS_DIR,
    DATA_DIR,
    DEFAULT_DT_S,
    build_frame,
    build_model,
    build_noise,
    csv_dump,
    json_dump,
    ladder_sequence_times,
    readout_dump_rotating_frequency,
    readout_photon_number,
    simulate_single_stage,
    state_population,
    storage_photon_number,
    storage_sideband_lab_frequency,
    storage_sideband_rotating_frequency,
    transmon_excited_population,
)
from run_study import make_family_pulse


RESULTS_PATH = DATA_DIR / "study_results.json"
N_TARGET = 4


def load_best_cases() -> tuple[dict, dict]:
    payload = json.loads(RESULTS_PATH.read_text(encoding="utf-8"))
    return payload["best_storage"][str(N_TARGET)], payload["best_dump"][str(N_TARGET)]


def simulate_variant(
    *,
    n_tr: int,
    n_storage: int,
    n_readout: int,
    dt_s: float,
    storage_case: dict,
    dump_case: dict,
) -> dict:
    model = build_model(n_tr=n_tr, n_storage=n_storage, n_readout=n_readout)
    frame = build_frame(model)
    noise = build_noise()

    storage_duration_s = float(storage_case["duration_ns"]) * 1.0e-9
    dump_duration_s = float(dump_case["duration_ns"]) * 1.0e-9
    ringdown_s = ladder_sequence_times(noise)

    storage_pulse = make_family_pulse(
        channel="storage_sb",
        carrier_rad_s=storage_sideband_rotating_frequency(model, frame, N_TARGET),
        duration_s=storage_duration_s,
        amplitude_hz=float(storage_case["amplitude_MHz"]) * 1.0e6,
        family=str(storage_case["family"]),
        label="conv_storage",
    )
    dump_pulse = make_family_pulse(
        channel="readout_sb",
        carrier_rad_s=readout_dump_rotating_frequency(model, frame, N_TARGET),
        duration_s=dump_duration_s,
        amplitude_hz=float(dump_case["amplitude_MHz"]) * 1.0e6,
        family=str(dump_case["family"]),
        label="conv_dump",
    )

    initial = model.basis_state(0, N_TARGET, 0)
    _, storage_result = simulate_single_stage(
        model,
        initial,
        pulse=storage_pulse,
        duration_s=storage_duration_s,
        drive_ops={"storage_sb": SidebandDriveSpec(mode="storage", lower_level=0, upper_level=2, sideband="red")},
        frame=frame,
        noise=noise,
        dt_s=dt_s,
        store_states=False,
    )
    _, dump_result = simulate_single_stage(
        model,
        storage_result.final_state,
        pulse=dump_pulse,
        duration_s=dump_duration_s,
        drive_ops={"readout_sb": SidebandDriveSpec(mode="readout", lower_level=0, upper_level=2, sideband="red")},
        frame=frame,
        noise=noise,
        dt_s=dt_s,
        store_states=False,
    )
    _, idle_result = simulate_single_stage(
        model,
        dump_result.final_state,
        pulse=None,
        duration_s=ringdown_s,
        drive_ops={},
        frame=frame,
        noise=noise,
        dt_s=4.0e-9,
        store_states=False,
    )
    final_state = idle_result.final_state

    return {
        "n_tr": n_tr,
        "n_storage": n_storage,
        "n_readout": n_readout,
        "dt_ns": dt_s * 1.0e9,
        "storage_sideband_lab_GHz": storage_sideband_lab_frequency(model, N_TARGET) / (2.0 * 3.141592653589793 * 1.0e9),
        "success_probability": state_population(final_state, model.basis_state(0, N_TARGET - 1, 0)),
        "final_mean_storage_n": storage_photon_number(final_state),
        "final_readout_n": readout_photon_number(final_state),
        "final_transmon_excited_population": transmon_excited_population(final_state),
    }


def main() -> None:
    start = time.time()
    storage_case, dump_case = load_best_cases()
    baseline_case = {"label": "baseline", "n_tr": 4, "n_storage": 7, "n_readout": 3, "dt_s": DEFAULT_DT_S}
    cases = [
        baseline_case,
        {"label": "transmon_5", "n_tr": 5, "n_storage": 7, "n_readout": 3, "dt_s": DEFAULT_DT_S},
        {"label": "storage_6", "n_tr": 4, "n_storage": 6, "n_readout": 3, "dt_s": DEFAULT_DT_S},
        {"label": "storage_8", "n_tr": 4, "n_storage": 8, "n_readout": 3, "dt_s": DEFAULT_DT_S},
        {"label": "readout_2", "n_tr": 4, "n_storage": 7, "n_readout": 2, "dt_s": DEFAULT_DT_S},
        {"label": "readout_4", "n_tr": 4, "n_storage": 7, "n_readout": 4, "dt_s": DEFAULT_DT_S},
        {"label": "dt_0p125", "n_tr": 4, "n_storage": 7, "n_readout": 3, "dt_s": 0.125e-9},
        {"label": "dt_0p5", "n_tr": 4, "n_storage": 7, "n_readout": 3, "dt_s": 0.5e-9},
    ]

    rows = []
    for case in cases:
        row = simulate_variant(
            n_tr=int(case["n_tr"]),
            n_storage=int(case["n_storage"]),
            n_readout=int(case["n_readout"]),
            dt_s=float(case["dt_s"]),
            storage_case=storage_case,
            dump_case=dump_case,
        )
        row["label"] = str(case["label"])
        rows.append(row)

    baseline = next(row for row in rows if row["label"] == "baseline")
    for row in rows:
        row["delta_frequency_kHz"] = (row["storage_sideband_lab_GHz"] - baseline["storage_sideband_lab_GHz"]) * 1.0e6
        row["delta_success_probability"] = row["success_probability"] - baseline["success_probability"]
        row["delta_final_mean_storage_n"] = row["final_mean_storage_n"] - baseline["final_mean_storage_n"]
        row["delta_final_readout_n"] = row["final_readout_n"] - baseline["final_readout_n"]
        row["delta_final_transmon_excited_population"] = (
            row["final_transmon_excited_population"] - baseline["final_transmon_excited_population"]
        )

    payload = {
        "n_target": N_TARGET,
        "baseline": baseline,
        "rows": rows,
        "runtime_s": time.time() - start,
    }
    csv_dump(DATA_DIR / "convergence_checks.csv", rows)
    json_dump(ARTIFACTS_DIR / "convergence_checks.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
