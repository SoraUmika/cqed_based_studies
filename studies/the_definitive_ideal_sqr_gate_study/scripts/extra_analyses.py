from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[3]
NATIVE_RICH_SCRIPTS = REPO_ROOT / "studies" / "native_rich_multitone_sqr_cpsqr_feasibility" / "scripts"

if str(NATIVE_RICH_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(NATIVE_RICH_SCRIPTS))

# cqed_sim does not expose a one-call helper for replaying a stored pulse and
# extracting manifold-resolved diagnostics, so we reuse the already-audited
# helper layer from the patched native-rich study instead of re-implementing it.
from common import (  # type: ignore[import-not-found]
    CHI,
    DEFAULT_DT,
    DEFAULT_SIGMA_FRACTION,
    IDEAL_X_PI,
    PI_PULSE_DURATION_S,
    PI_PULSE_SIGMA_FRACTION,
    build_frame,
    build_model,
    compile_pulse_sequence,
    duration_from_chi_t,
    logical_levels,
    make_gaussian_qubit_rotation_pulse,
    manifold_transition_frequencies_hz,
    simulate_full_operator_on_logical_inputs,
)
from metrics import (  # type: ignore[import-not-found]
    channel_process_fidelity_to_unitary,
    nearest_unitary,
    qubit_channel_kraus_from_full,
    same_manifold_block,
    unitary_rotation_parameters,
)


XPI_DURATION_GRID_NS = (20.0, 40.0, 60.0, 80.0, 100.0, 120.0)
XPI_SIGMA_GRID = (0.18, 0.25, 0.33)
XPI_AMPLITUDE_SCALE_GRID = (0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15)
XPI_CARRIER_LEVEL_GRID = (0, 1, 2)
XPI_EVAL_LEVELS = tuple(range(5))

SPECTRAL_DURATION_GRID = (3.0, 5.0, 8.0, 12.0, 20.0)


def _model_rows() -> tuple[tuple[str, bool], ...]:
    return (("chi_only", False), ("chi_plus_chiprime", True))


def _axis_components(phi_rad: float, axis_z: float) -> tuple[float, float, float]:
    clipped_axis_z = float(np.clip(axis_z, -1.0, 1.0))
    radial = float(np.sqrt(max(0.0, 1.0 - clipped_axis_z * clipped_axis_z)))
    return (
        float(radial * np.cos(phi_rad)),
        float(radial * np.sin(phi_rad)),
        clipped_axis_z,
    )


def _simulate_refocusing_config(
    *,
    model_variant: str,
    include_chi_prime: bool,
    duration_ns: float,
    sigma_fraction: float,
    amplitude_scale: float,
    carrier_level: int,
    variant: str,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    model = build_model(include_chi_prime=include_chi_prime, n_active=max(XPI_EVAL_LEVELS) + 1)
    frame = build_frame(model)
    pulse = make_gaussian_qubit_rotation_pulse(
        model,
        frame,
        theta=np.pi * float(amplitude_scale),
        phase=0.0,
        duration_s=float(duration_ns) * 1.0e-9,
        manifold_level=int(carrier_level),
        sigma_fraction=float(sigma_fraction),
        t0=0.0,
        label=f"{variant}_{model_variant}_xpi",
    )
    compiled = compile_pulse_sequence(
        [pulse],
        dt_s=float(DEFAULT_DT),
        total_duration_s=float(duration_ns) * 1.0e-9,
    )
    full_operator = simulate_full_operator_on_logical_inputs(
        model,
        compiled,
        frame=frame,
        drive_ops={str(pulse.channel): "qubit"},
        levels=logical_levels(max(XPI_EVAL_LEVELS) + 1),
    )

    rows: list[dict[str, Any]] = []
    fidelity_values: list[float] = []
    for level in XPI_EVAL_LEVELS:
        same_block = nearest_unitary(same_manifold_block(full_operator, int(model.n_cav), int(level)))
        theta_rad, phi_rad, axis_z = unitary_rotation_parameters(same_block)
        axis_x, axis_y, axis_z_val = _axis_components(phi_rad, axis_z)
        kraus_ops = qubit_channel_kraus_from_full(full_operator, int(model.n_cav), int(level))
        process_fidelity = float(channel_process_fidelity_to_unitary(kraus_ops, IDEAL_X_PI))
        fidelity_values.append(process_fidelity)
        rows.append(
            {
                "model_variant": str(model_variant),
                "include_chi_prime": bool(include_chi_prime),
                "variant": str(variant),
                "level": int(level),
                "duration_ns": float(duration_ns),
                "sigma_fraction": float(sigma_fraction),
                "amplitude_scale": float(amplitude_scale),
                "carrier_level": int(carrier_level),
                "theta_rad": float(theta_rad),
                "axis_x": float(axis_x),
                "axis_y": float(axis_y),
                "axis_z": float(axis_z_val),
                "process_fidelity": process_fidelity,
            }
        )

    summary = {
        "model_variant": str(model_variant),
        "include_chi_prime": bool(include_chi_prime),
        "variant": str(variant),
        "duration_ns": float(duration_ns),
        "sigma_fraction": float(sigma_fraction),
        "amplitude_scale": float(amplitude_scale),
        "carrier_level": int(carrier_level),
        "worst_process_fidelity": float(min(fidelity_values)),
        "mean_process_fidelity": float(np.mean(fidelity_values)),
    }
    return rows, summary


def compute_xpi_characterization() -> dict[str, Any]:
    all_rows: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    best_configs: list[dict[str, Any]] = []

    for model_variant, include_chi_prime in _model_rows():
        baseline_rows, baseline_summary = _simulate_refocusing_config(
            model_variant=model_variant,
            include_chi_prime=include_chi_prime,
            duration_ns=PI_PULSE_DURATION_S * 1.0e9,
            sigma_fraction=PI_PULSE_SIGMA_FRACTION,
            amplitude_scale=1.0,
            carrier_level=0,
            variant="baseline",
        )
        all_rows.extend(baseline_rows)
        summaries.append(baseline_summary)

        candidate_summaries: list[dict[str, Any]] = []
        best_summary: dict[str, Any] | None = None
        best_rows: list[dict[str, Any]] = []
        for duration_ns in XPI_DURATION_GRID_NS:
            for sigma_fraction in XPI_SIGMA_GRID:
                for amplitude_scale in XPI_AMPLITUDE_SCALE_GRID:
                    for carrier_level in XPI_CARRIER_LEVEL_GRID:
                        robust_rows, robust_summary = _simulate_refocusing_config(
                            model_variant=model_variant,
                            include_chi_prime=include_chi_prime,
                            duration_ns=duration_ns,
                            sigma_fraction=sigma_fraction,
                            amplitude_scale=amplitude_scale,
                            carrier_level=carrier_level,
                            variant="robust_candidate",
                        )
                        candidate_summaries.append(robust_summary)
                        if best_summary is None:
                            best_summary = robust_summary
                            best_rows = robust_rows
                            continue
                        incumbent = (best_summary["worst_process_fidelity"], best_summary["mean_process_fidelity"])
                        challenger = (robust_summary["worst_process_fidelity"], robust_summary["mean_process_fidelity"])
                        if challenger > incumbent:
                            best_summary = robust_summary
                            best_rows = robust_rows

        if best_summary is None:
            continue

        final_rows = []
        for row in best_rows:
            copied = dict(row)
            copied["variant"] = "robust"
            final_rows.append(copied)
        all_rows.extend(final_rows)
        robust_summary = dict(best_summary)
        robust_summary["variant"] = "robust"
        summaries.append(robust_summary)
        best_configs.append(robust_summary)

        summaries.extend(
            [
                {
                    "model_variant": str(model_variant),
                    "include_chi_prime": bool(include_chi_prime),
                    "variant": "search_space",
                    "candidate_count": len(candidate_summaries),
                    "best_worst_process_fidelity": float(best_summary["worst_process_fidelity"]),
                    "best_mean_process_fidelity": float(best_summary["mean_process_fidelity"]),
                }
            ]
        )

    return {
        "rows": all_rows,
        "summaries": summaries,
        "best_configs": best_configs,
        "notes": {
            "purpose": "Standalone manifold-resolved characterization of the finite Gaussian refocusing pulse used in the echoed source studies, plus a lightweight compromise pulse selected by worst-manifold process fidelity.",
            "eval_levels": list(XPI_EVAL_LEVELS),
            "duration_grid_ns": list(XPI_DURATION_GRID_NS),
            "sigma_fraction_grid": list(XPI_SIGMA_GRID),
            "amplitude_scale_grid": list(XPI_AMPLITUDE_SCALE_GRID),
            "carrier_level_grid": list(XPI_CARRIER_LEVEL_GRID),
        },
    }


def compute_spectral_crowding() -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    plot_rows: list[dict[str, Any]] = []
    for model_variant, include_chi_prime in _model_rows():
        model = build_model(include_chi_prime=include_chi_prime, n_active=8)
        frame = build_frame(model)
        levels = tuple(range(8))
        frequencies_hz = manifold_transition_frequencies_hz(model, levels, frame)
        gaps_hz = np.diff(frequencies_hz)

        for level, frequency_hz in zip(levels, frequencies_hz, strict=True):
            plot_rows.append(
                {
                    "model_variant": str(model_variant),
                    "level": int(level),
                    "transition_frequency_hz": float(frequency_hz),
                }
            )

        for chi_t in SPECTRAL_DURATION_GRID:
            duration_s = duration_from_chi_t(float(chi_t), chi=float(CHI))
            sigma_s = float(DEFAULT_SIGMA_FRACTION) * float(duration_s)
            bandwidth_hz = 1.0 / max(sigma_s, 1.0e-15)
            for idx, gap_hz in enumerate(gaps_hz):
                rows.append(
                    {
                        "model_variant": str(model_variant),
                        "lower_level": int(idx),
                        "upper_level": int(idx + 1),
                        "chi_t_over_2pi": float(chi_t),
                        "duration_ns": float(duration_s * 1.0e9),
                        "sigma_ns": float(sigma_s * 1.0e9),
                        "transition_gap_hz": float(abs(gap_hz)),
                        "tone_bandwidth_hz": float(bandwidth_hz),
                        "crowding_ratio": float(abs(bandwidth_hz / gap_hz)),
                    }
                )

    return {
        "rows": rows,
        "plot_rows": plot_rows,
        "notes": {
            "purpose": "Analytic spectral-crowding estimate for Gaussian direct-waveform tones.",
            "duration_grid_chi_t_over_2pi": list(SPECTRAL_DURATION_GRID),
            "bandwidth_definition": "The Gaussian tone bandwidth proxy is 1/sigma, where sigma = sigma_fraction * T and sigma_fraction = 1/6 from the direct-waveform seed.",
        },
    }


def compute_scaling_summary(df: pd.DataFrame) -> dict[str, Any]:
    native = df[
        (df["study"] == "native_rich_multitone_sqr_cpsqr_feasibility")
        & (df["target_kind"] == "ideal_sqr")
        & df["strict_process_fidelity"].notna()
        & df["target_family"].isin({"smooth_x", "staggered_x"})
    ].copy()
    if native.empty:
        return {"rows": [], "fits": [], "notes": {"status": "no_native_rows"}}

    best_rows = (
        native.sort_values("strict_process_fidelity", ascending=False)
        .drop_duplicates(subset=["construction", "model_variant", "n_active", "target_family", "case_id"])
        .copy()
    )
    grouped = (
        best_rows.groupby(["construction", "construction_display", "construction_family", "model_variant", "n_active"], as_index=False)[
            "strict_process_fidelity"
        ]
        .max()
        .sort_values(["construction", "model_variant", "n_active"])
    )

    fit_rows: list[dict[str, Any]] = []
    for (construction, model_variant), group in grouped.groupby(["construction", "model_variant"]):
        ordered = group.sort_values("n_active")
        if len(ordered) < 2:
            continue
        x = ordered["n_active"].to_numpy(dtype=float)
        residual = np.maximum(1.0e-6, 1.0 - ordered["strict_process_fidelity"].to_numpy(dtype=float))
        coeffs = np.polyfit(x, np.log(residual), deg=1)
        slope = float(coeffs[0])
        intercept = float(coeffs[1])
        predicted_n8 = float(np.clip(1.0 - np.exp(intercept + slope * 8.0), 0.0, 1.0))
        fit_rows.append(
            {
                "construction": str(construction),
                "model_variant": str(model_variant),
                "fit_kind": "log-linear residual vs N_active",
                "slope": slope,
                "intercept": intercept,
                "predicted_strict_process_at_n8": predicted_n8,
            }
        )

    return {
        "rows": grouped.to_dict(orient="records"),
        "fits": fit_rows,
        "notes": {
            "purpose": "Strict ideal-SQR scaling summary across the native-rich corpus.",
            "target_families": ["smooth_x", "staggered_x"],
        },
    }


def xpi_dataframe(payload: dict[str, Any]) -> pd.DataFrame:
    return pd.DataFrame(payload.get("rows", []))


def spectral_crowding_dataframe(payload: dict[str, Any]) -> pd.DataFrame:
    return pd.DataFrame(payload.get("rows", []))


def scaling_dataframe(payload: dict[str, Any]) -> pd.DataFrame:
    return pd.DataFrame(payload.get("rows", []))
