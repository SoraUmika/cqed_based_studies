"""Experiment-facing spectroscopy follow-up for the gf-sideband cooling study."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from cqed_sim.core.drive_targets import SidebandDriveSpec

from common import (
    ARTIFACTS_DIR,
    DATA_DIR,
    FIGURES_DIR,
    basis_label,
    build_frame,
    build_model,
    csv_dump,
    hz,
    json_dump,
    plot_save,
    readout_dump_lab_frequency,
    readout_dump_rotating_frequency,
    readout_photon_number,
    sideband_matrix_element,
    simulate_single_stage,
    state_population,
    storage_sideband_lab_frequency,
    storage_sideband_rotating_frequency,
    transmon_level_populations,
)
from run_study import make_family_pulse


DETUNING_OFFSETS_MHZ = np.linspace(-10.0, 10.0, 41)
SPECTROSCOPY_DURATION_S = 400.0e-9
SPECTROSCOPY_DT_S = 0.5e-9
STORAGE_PROBE_BASE_AMPLITUDE_MHZ = 0.30
READOUT_PROBE_AMPLITUDE_MHZ = 0.30


def load_previous_artifacts() -> tuple[dict[str, object], list[dict[str, object]]]:
    study_results = json.loads((DATA_DIR / "study_results.json").read_text(encoding="utf-8"))
    direct_carrier = json.loads((ARTIFACTS_DIR / "direct_carrier_comparison.json").read_text(encoding="utf-8"))
    return study_results, direct_carrier


def storage_probe_amplitude_mhz(model, n: int) -> float:
    matrix_element = sideband_matrix_element(model, mode="storage", n=n)
    return float(STORAGE_PROBE_BASE_AMPLITUDE_MHZ / matrix_element)


def step_a_followup_scan(model, frame, best_dump: dict[int, dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for n in range(1, 5):
        probe_amplitude = storage_probe_amplitude_mhz(model, n)
        base_frequency = storage_sideband_rotating_frequency(model, frame, n)
        witness_dump = best_dump[n]
        dump_pulse = make_family_pulse(
            channel="readout_sb",
            carrier_rad_s=readout_dump_rotating_frequency(model, frame, n),
            duration_s=float(witness_dump["duration_ns"]) * 1.0e-9,
            amplitude_hz=float(witness_dump["amplitude_MHz"]) * 1.0e6,
            family=str(witness_dump["family"]),
            label=f"step_a_dump_witness_n{n}",
        )
        for detuning_mhz in DETUNING_OFFSETS_MHZ:
            storage_pulse = make_family_pulse(
                channel="storage_sb",
                carrier_rad_s=base_frequency + 2.0 * np.pi * float(detuning_mhz) * 1.0e6,
                duration_s=SPECTROSCOPY_DURATION_S,
                amplitude_hz=probe_amplitude * 1.0e6,
                family="square",
                label=f"step_a_probe_n{n}",
            )
            _, storage_result = simulate_single_stage(
                model,
                model.basis_state(0, n, 0),
                pulse=storage_pulse,
                duration_s=SPECTROSCOPY_DURATION_S,
                drive_ops={"storage_sb": SidebandDriveSpec(mode="storage", lower_level=0, upper_level=2, sideband="red")},
                frame=frame,
                noise=None,
                dt_s=SPECTROSCOPY_DT_S,
                store_states=False,
            )
            storage_state = storage_result.final_state
            _, dump_result = simulate_single_stage(
                model,
                storage_state,
                pulse=dump_pulse,
                duration_s=float(witness_dump["duration_ns"]) * 1.0e-9,
                drive_ops={"readout_sb": SidebandDriveSpec(mode="readout", lower_level=0, upper_level=2, sideband="red")},
                frame=frame,
                noise=None,
                dt_s=SPECTROSCOPY_DT_S,
                store_states=False,
            )
            dump_state = dump_result.final_state
            transmon_populations = transmon_level_populations(storage_state)
            rows.append(
                {
                    "n": n,
                    "detuning_MHz": float(detuning_mhz),
                    "probe_duration_ns": SPECTROSCOPY_DURATION_S * 1.0e9,
                    "probe_amplitude_MHz": probe_amplitude,
                    "storage_target_probability": float(state_population(storage_state, model.basis_state(2, n - 1, 0))),
                    "storage_total_f_probability": float(transmon_populations.get(2, 0.0)),
                    "dump_witness_target_probability": float(state_population(dump_state, model.basis_state(0, n - 1, 1))),
                    "dump_witness_readout_n": float(readout_photon_number(dump_state)),
                }
            )
    return rows


def step_b_followup_scan(model, frame) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for n in range(1, 5):
        base_frequency = readout_dump_rotating_frequency(model, frame, n)
        for detuning_mhz in DETUNING_OFFSETS_MHZ:
            readout_pulse = make_family_pulse(
                channel="readout_sb",
                carrier_rad_s=base_frequency + 2.0 * np.pi * float(detuning_mhz) * 1.0e6,
                duration_s=SPECTROSCOPY_DURATION_S,
                amplitude_hz=READOUT_PROBE_AMPLITUDE_MHZ * 1.0e6,
                family="square",
                label=f"step_b_probe_n{n}",
            )
            _, result = simulate_single_stage(
                model,
                model.basis_state(2, n - 1, 0),
                pulse=readout_pulse,
                duration_s=SPECTROSCOPY_DURATION_S,
                drive_ops={"readout_sb": SidebandDriveSpec(mode="readout", lower_level=0, upper_level=2, sideband="red")},
                frame=frame,
                noise=None,
                dt_s=SPECTROSCOPY_DT_S,
                store_states=False,
            )
            final_state = result.final_state
            transmon_populations = transmon_level_populations(final_state)
            rows.append(
                {
                    "n": n,
                    "detuning_MHz": float(detuning_mhz),
                    "probe_duration_ns": SPECTROSCOPY_DURATION_S * 1.0e9,
                    "probe_amplitude_MHz": READOUT_PROBE_AMPLITUDE_MHZ,
                    "dump_target_probability": float(state_population(final_state, model.basis_state(0, n - 1, 1))),
                    "readout_n": float(readout_photon_number(final_state)),
                    "final_f_probability": float(transmon_populations.get(2, 0.0)),
                }
            )
    return rows


def exact_transition_summary(model, frame, direct_carrier_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    direct_lookup = {int(row["n"]): row for row in direct_carrier_rows}
    storage_lines = {n: hz(storage_sideband_lab_frequency(model, n)) / 1.0e6 for n in range(1, 5)}
    readout_lines = {n: hz(readout_dump_lab_frequency(model, n)) / 1.0e6 for n in range(1, 5)}
    rows: list[dict[str, object]] = []
    for n in range(1, 5):
        direct_lab_frequency_hz = hz(model.basis_energy(2, n, 0) - model.basis_energy(0, n, 0))
        nearest_storage_spacing_mhz = min(abs(storage_lines[n] - storage_lines[other]) for other in storage_lines if other != n)
        nearest_readout_spacing_mhz = min(abs(readout_lines[n] - readout_lines[other]) for other in readout_lines if other != n)
        rows.append(
            {
                "n": n,
                "step_a_initial_state": basis_label(0, n, 0),
                "step_a_target_state": basis_label(2, n - 1, 0),
                "step_b_initial_state": basis_label(2, n - 1, 0),
                "step_b_target_state": basis_label(0, n - 1, 1),
                "direct_gf_carrier_lab_GHz": round(direct_lab_frequency_hz / 1.0e9, 9),
                "direct_gf_comment": str(direct_lookup[n]["comment"]),
                "step_a_lab_GHz": round(hz(storage_sideband_lab_frequency(model, n)) / 1.0e9, 9),
                "step_b_lab_GHz": round(hz(readout_dump_lab_frequency(model, n)) / 1.0e9, 9),
                "step_a_rot_MHz": round(hz(storage_sideband_rotating_frequency(model, frame, n)) / 1.0e6, 6),
                "step_b_rot_MHz": round(hz(readout_dump_rotating_frequency(model, frame, n)) / 1.0e6, 6),
                "step_a_nearest_spacing_MHz": round(nearest_storage_spacing_mhz, 6),
                "step_b_nearest_spacing_MHz": round(nearest_readout_spacing_mhz, 6),
                "step_a_probe_amplitude_MHz": round(storage_probe_amplitude_mhz(model, n), 6),
                "step_b_probe_amplitude_MHz": READOUT_PROBE_AMPLITUDE_MHZ,
                "recommended_step_a_observable": "P_f immediately after the storage-sideband probe",
                "fallback_step_a_observable": "fixed Step B dump witness: P(|g,1_r,n_s-1>) or integrated readout ringdown",
                "recommended_step_b_observable": "readout occupation or integrated homodyne ringdown",
            }
        )
    return rows


def peak_summary(rows: list[dict[str, object]], *, target_key: str) -> list[dict[str, object]]:
    summaries: list[dict[str, object]] = []
    for n in range(1, 5):
        subset = [row for row in rows if int(row["n"]) == n]
        best = max(subset, key=lambda row: float(row[target_key]))
        summaries.append(
            {
                "n": n,
                "peak_detuning_MHz": float(best["detuning_MHz"]),
                "peak_target_probability": float(best[target_key]),
            }
        )
    return summaries


def make_step_a_figure(rows: list[dict[str, object]]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True, sharey=True)
    for ax, n in zip(axes.ravel(), range(1, 5), strict=True):
        subset = [row for row in rows if int(row["n"]) == n]
        detuning = [float(row["detuning_MHz"]) for row in subset]
        target = [float(row["storage_target_probability"]) for row in subset]
        witness = [float(row["dump_witness_target_probability"]) for row in subset]
        ax.plot(detuning, target, label=r"$P_{f,n-1}$")
        ax.plot(detuning, witness, "--", label=r"dump witness $P_{g,1_r,n-1}$")
        ax.set_title(f"Step A spectroscopy, n={n}")
        ax.set_xlabel("Detuning (MHz)")
        ax.set_ylabel("Signal")
        ax.axvline(0.0, color="k", alpha=0.2, linewidth=1.0)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False)
    fig.suptitle("Storage-sideband spectroscopy observables")
    plot_save(fig, "spectroscopy_followup_step_a")


def make_step_b_figure(rows: list[dict[str, object]]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True, sharey=True)
    for ax, n in zip(axes.ravel(), range(1, 5), strict=True):
        subset = [row for row in rows if int(row["n"]) == n]
        detuning = [float(row["detuning_MHz"]) for row in subset]
        target = [float(row["dump_target_probability"]) for row in subset]
        readout_n = [float(row["readout_n"]) for row in subset]
        ax.plot(detuning, target, label=r"$P_{g,1_r,n-1}$")
        ax.plot(detuning, readout_n, "--", label=r"$\langle n_r \rangle$")
        ax.set_title(f"Step B spectroscopy, n={n}")
        ax.set_xlabel("Detuning (MHz)")
        ax.set_ylabel("Signal")
        ax.axvline(0.0, color="k", alpha=0.2, linewidth=1.0)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False)
    fig.suptitle("Readout-dump spectroscopy observables")
    plot_save(fig, "spectroscopy_followup_step_b")


def write_followup_note(
    summary_rows: list[dict[str, object]],
    step_a_peaks: list[dict[str, object]],
    step_b_peaks: list[dict[str, object]],
) -> None:
    step_a_peak_lookup = {int(row["n"]): row for row in step_a_peaks}
    step_b_peak_lookup = {int(row["n"]): row for row in step_b_peaks}
    lines = [
        "# Experiment-Facing Spectroscopy Follow-Up",
        "",
        "## Main answers",
        "- To locate the Step A storage-cooling transition, the primary observable should be the population in `|f,0_r,n_s-1>` after a weak long storage-sideband probe. In practice this means read out `P_f`, not an `X/Y` projection, during the first spectroscopy pass.",
        "- A direct transmon `g-f` carrier probe is **not** the cooling transition because it keeps the storage photon number fixed. It is still useful as an auxiliary calibration for qutrit readout and for preparing `|f,0_r,n_s-1>` when calibrating Step B.",
        "- If direct `f` discrimination is weak, the best fallback is a dump-assisted witness: scan Step A, then apply the calibrated Step B dump pulse and read out the resulting readout occupation or integrated ringdown.",
        "- `X/Y` projections in the `g-f` manifold are needed only after the resonance is found, for Ramsey/tomography-style checks of phase, Stark shift, and axis control.",
        "- Neighboring `n_s`-resolved Step A and Step B lines are only about `5.680842 MHz` apart, so the first search should use weak, long probes rather than the short high-fidelity control pulses.",
        "",
        "## Exact transition search workflow",
        "1. Prepare `|g,0_r,n_s>` for the desired storage manifold.",
        "2. Coarse scan the Step A sideband around the predicted line with a weak square probe of duration `400 ns` and amplitude `0.30/sqrt(n_s) MHz`.",
        "3. Measure either direct `P_f` or the dump-assisted readout witness.",
        "4. Fit the peak and then repeat a fine scan over `+-2 MHz` with smaller steps.",
        "5. At the fitted frequency, perform a duration sweep to calibrate the Rabi period and then switch to the faster high-fidelity pulse family from the main study.",
        "6. For Step B, prepare `|f,0_r,n_s-1>` using either the calibrated Step A pulse or the direct transmon `g-f` carrier, then scan the readout sideband and monitor readout occupation.",
        "7. Only after those `Z`-axis calibrations are stable should you run `g-f` Ramsey or analysis-pulse tomography to project onto `X_{gf}` and `Y_{gf}`.",
        "",
        "## Predicted exact lines and recommended observables",
        "| n_s | Direct g-f carrier (GHz) | Step A sideband (GHz) | Step B dump (GHz) | Preferred Step A readout | Preferred Step B readout |",
        "|---|---:|---:|---:|---|---|",
    ]
    for row in summary_rows:
        n = int(row["n"])
        lines.append(
            "| "
            f"{n} | {float(row['direct_gf_carrier_lab_GHz']):.9f} | {float(row['step_a_lab_GHz']):.9f} | "
            f"{float(row['step_b_lab_GHz']):.9f} | `P_f` or dump witness | `\\langle n_r \\rangle` or integrated ringdown |"
        )
    lines.extend(
        [
            "",
            "## Simulated weak-drive spectroscopy peaks",
            "| n_s | Step A peak detuning (MHz) | Step A peak signal | Step B peak detuning (MHz) | Step B peak signal |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for row in summary_rows:
        n = int(row["n"])
        lines.append(
            "| "
            f"{n} | {float(step_a_peak_lookup[n]['peak_detuning_MHz']):+.2f} | "
            f"{float(step_a_peak_lookup[n]['peak_target_probability']):.4f} | "
            f"{float(step_b_peak_lookup[n]['peak_detuning_MHz']):+.2f} | "
            f"{float(step_b_peak_lookup[n]['peak_target_probability']):.4f} |"
        )
    lines.extend(
        [
            "",
            "## What to measure in each experiment",
            "- **Step A coarse spectroscopy:** measure `P_f`. This is a `Z_{gf}` measurement, not an `X/Y` projection.",
            "- **Step A fallback when `f` readout is weak:** append the calibrated Step B pulse and detect the induced readout photon as a dump witness.",
            "- **Step A Rabi calibration:** stay with `P_f` versus pulse duration.",
            "- **Step A Stark-shift and axis check:** switch to `g-f` Ramsey with a final `pi/2` analysis pulse so that `X_{gf}` or `Y_{gf}` is mapped onto readout.",
            "- **Step B spectroscopy:** measure readout occupation or integrated homodyne ringdown, because the target state already lives in `|g,1_r,n_s-1>`.",
            "- **Full cooling validation:** compare storage number-splitting or Wigner/tomography before and after repeated cycles, not just transmon population.",
            "",
            "## Most useful experiment sequence",
            "The cleanest lab sequence is:",
            "1. Calibrate qutrit transmon readout (`g/e/f`) and the direct `g-f` carrier first.",
            "2. Use the storage-sideband probe to find Step A at each `n_s`.",
            "3. Use the direct `g-f` carrier from `|g,0_r,n_s-1>` to prepare `|f,0_r,n_s-1>` and calibrate Step B independently.",
            "4. Reconnect the two calibrated pieces into the full cooling primitive.",
            "",
            "The central practical point is that you should not begin with full `g-f` tomography. Begin with population spectroscopy in the `g/f` manifold, then add `X/Y` projections only after the line center and pulse area are known.",
        ]
    )
    (ARTIFACTS_DIR / "spectroscopy_measurement_followup.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    model = build_model()
    frame = build_frame(model)
    study_results, direct_carrier_rows = load_previous_artifacts()
    best_dump = {int(key): value for key, value in study_results["best_dump"].items()}

    step_a_rows = step_a_followup_scan(model, frame, best_dump)
    step_b_rows = step_b_followup_scan(model, frame)
    summary_rows = exact_transition_summary(model, frame, direct_carrier_rows)
    step_a_peaks = peak_summary(step_a_rows, target_key="storage_target_probability")
    step_b_peaks = peak_summary(step_b_rows, target_key="dump_target_probability")

    csv_dump(DATA_DIR / "spectroscopy_followup_step_a.csv", step_a_rows)
    csv_dump(DATA_DIR / "spectroscopy_followup_step_b.csv", step_b_rows)
    csv_dump(DATA_DIR / "spectroscopy_followup_summary.csv", summary_rows)
    json_dump(
        ARTIFACTS_DIR / "spectroscopy_followup.json",
        {
            "summary_rows": summary_rows,
            "step_a_peak_summary": step_a_peaks,
            "step_b_peak_summary": step_b_peaks,
            "detuning_offsets_MHz": [float(value) for value in DETUNING_OFFSETS_MHZ],
            "probe_duration_ns": SPECTROSCOPY_DURATION_S * 1.0e9,
            "storage_probe_base_amplitude_MHz": STORAGE_PROBE_BASE_AMPLITUDE_MHZ,
            "readout_probe_amplitude_MHz": READOUT_PROBE_AMPLITUDE_MHZ,
        },
    )
    write_followup_note(summary_rows, step_a_peaks, step_b_peaks)
    make_step_a_figure(step_a_rows)
    make_step_b_figure(step_b_rows)


if __name__ == "__main__":
    main()
