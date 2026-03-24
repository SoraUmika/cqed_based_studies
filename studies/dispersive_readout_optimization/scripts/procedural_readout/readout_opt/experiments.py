"""Orchestration for the procedural readout optimization study."""

from __future__ import annotations

from dataclasses import asdict, replace
import json
from pathlib import Path
from typing import Any

import numpy as np

from .config import DATA_DIR, DEFAULT_CONFIG, ReadoutStudyConfig
from .optimize import OptimizationOutcome, optimize_family
from .pulse_families import PulseDesign, get_family, set_nulling_tail_kappa
from .simulate import ReplayEvaluation, evaluate_full_design

LINEAR_FAMILIES = (
    "square",
    "smooth_square",
    "ring_hold",
    "procedural_segments",
    "nulling_tail",
    "fourier_basis",
    "piecewise_reference",
)

FULL_FAMILIES = (
    "square",
    "smooth_square",
    "ring_hold",
    "procedural_segments",
    "nulling_tail",
    "fourier_basis",
)

REPRESENTATIVE_FAMILIES = (
    "square",
    "procedural_segments",
    "nulling_tail",
    "fourier_basis",
    "piecewise_reference",
)


def _to_serializable(obj: Any) -> Any:
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    if isinstance(obj, complex):
        return {"real": float(obj.real), "imag": float(obj.imag)}
    if isinstance(obj, dict):
        return {str(key): _to_serializable(value) for key, value in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_serializable(value) for value in obj]
    return obj


def design_to_record(design: PulseDesign) -> dict[str, Any]:
    return {
        "family": design.family,
        "params": _to_serializable(design.params),
        "duration_ns": float(design.duration * 1.0e9),
        "dt_ns": float(design.dt * 1.0e9),
        "delta_g_mhz": float(design.delta_g / (2.0 * np.pi * 1.0e6)),
        "metadata": _to_serializable(design.metadata),
    }


def evaluation_to_record(evaluation: ReplayEvaluation) -> dict[str, Any]:
    return {
        "regime": evaluation.regime,
        "design": design_to_record(evaluation.design),
        "metrics": _to_serializable(asdict(evaluation.metrics)),
        "metadata": _to_serializable(evaluation.metadata),
    }


def outcome_to_record(outcome: OptimizationOutcome) -> dict[str, Any]:
    return {
        "family": outcome.family,
        "regime": outcome.regime,
        "objective_name": outcome.objective_name,
        "duration_ns": float(outcome.duration * 1.0e9),
        "x_best": _to_serializable(outcome.x_best),
        "evaluation": evaluation_to_record(outcome.evaluation),
    }


def run_linear_sweep(cfg: ReadoutStudyConfig) -> dict[str, dict[float, OptimizationOutcome]]:
    results: dict[str, dict[float, OptimizationOutcome]] = {family: {} for family in LINEAR_FAMILIES}
    for duration in cfg.duration_grid:
        for family in LINEAR_FAMILIES:
            objective = "info" if family == "piecewise_reference" else "balanced"
            seed = cfg.seed + int(round(duration * 1.0e9))
            if family == "piecewise_reference":
                outcome = optimize_family(
                    family=family,
                    regime="linear",
                    objective_name=objective,
                    duration=float(duration),
                    cfg=cfg,
                    seed=seed,
                    n_random=24,
                    n_local=2,
                    maxiter=36,
                )
            else:
                outcome = optimize_family(
                    family=family,
                    regime="linear",
                    objective_name=objective,
                    duration=float(duration),
                    cfg=cfg,
                    seed=seed,
                    n_random=16,
                    n_local=2,
                    maxiter=30,
                )
            results[family][float(duration)] = outcome
    return results


def run_full_sweep(
    cfg: ReadoutStudyConfig,
    linear_results: dict[str, dict[float, OptimizationOutcome]],
) -> dict[str, dict[float, OptimizationOutcome]]:
    results: dict[str, dict[float, OptimizationOutcome]] = {family: {} for family in FULL_FAMILIES}
    for duration in cfg.duration_grid:
        for family in FULL_FAMILIES:
            warm_start = linear_results[family][float(duration)].x_best
            outcome = optimize_family(
                family=family,
                regime="full",
                objective_name="balanced",
                duration=float(duration),
                cfg=cfg,
                warm_start=warm_start,
                seed=cfg.seed + 1000 + int(round(duration * 1.0e9)),
                n_random=10,
                n_local=2,
                maxiter=22,
            )
            results[family][float(duration)] = outcome
    return results


def compute_bound_hierarchy(
    cfg: ReadoutStudyConfig,
    linear_results: dict[str, dict[float, OptimizationOutcome]],
    full_results: dict[str, dict[float, OptimizationOutcome]],
) -> list[dict[str, float | str]]:
    hierarchy: list[dict[str, float | str]] = []
    for duration in cfg.duration_grid:
        ideal = linear_results["piecewise_reference"][float(duration)].evaluation.metrics.fidelity_ideal
        detector = linear_results["piecewise_reference"][float(duration)].evaluation.metrics.fidelity_eta
        t1_bound = linear_results["piecewise_reference"][float(duration)].evaluation.metrics.fidelity_t1_bound
        family_best = max(
            FULL_FAMILIES,
            key=lambda family: full_results[family][float(duration)].evaluation.metrics.score_balanced,
        )
        full_best = full_results[family_best][float(duration)].evaluation.metrics.fidelity_eta
        hierarchy.append(
            {
                "duration_ns": float(duration * 1.0e9),
                "ideal_bound": float(ideal),
                "detector_bound": float(detector),
                "t1_bound": float(t1_bound),
                "realistic_best": float(full_best),
                "best_family": family_best,
            }
        )
    return hierarchy


def representative_catalog(
    cfg: ReadoutStudyConfig,
    linear_results: dict[str, dict[float, OptimizationOutcome]],
    full_results: dict[str, dict[float, OptimizationOutcome]],
) -> dict[str, ReplayEvaluation]:
    duration = float(cfg.representative_duration)
    rep: dict[str, ReplayEvaluation] = {}
    for family in REPRESENTATIVE_FAMILIES:
        if family == "piecewise_reference":
            rep[family] = linear_results[family][duration].evaluation
        elif family in full_results:
            rep[family] = full_results[family][duration].evaluation
    return rep


def run_tradeoff_slice(
    cfg: ReadoutStudyConfig,
    linear_results: dict[str, dict[float, OptimizationOutcome]],
    full_results: dict[str, dict[float, OptimizationOutcome]],
) -> dict[str, dict[str, dict[str, Any]]]:
    duration = float(cfg.representative_duration)
    slice_results: dict[str, dict[str, dict[str, Any]]] = {}
    for family in ("square", "procedural_segments", "nulling_tail", "fourier_basis"):
        warm_start = linear_results[family][duration].x_best
        slice_results[family] = {}
        for objective in ("info", "balanced", "emptying"):
            if objective == "balanced":
                outcome = full_results[family][duration]
            else:
                outcome = optimize_family(
                    family=family,
                    regime="full",
                    objective_name=objective,
                    duration=duration,
                    cfg=cfg,
                    warm_start=warm_start,
                    seed=cfg.seed + 2000 + len(slice_results) * 100,
                    n_random=8,
                    n_local=2,
                    maxiter=18,
                )
            slice_results[family][objective] = outcome_to_record(outcome)
    return slice_results


def _scaled_design(design: PulseDesign, *, amp_scale: float = 1.0, phase_shift: float = 0.0, delta_shift: float = 0.0) -> PulseDesign:
    return PulseDesign(
        family=design.family,
        params=np.asarray(design.params, dtype=float).copy(),
        waveform=np.asarray(design.waveform, dtype=np.complex128) * amp_scale * np.exp(1j * phase_shift),
        dt=design.dt,
        duration=design.duration,
        delta_g=float(design.delta_g + delta_shift),
        metadata=dict(design.metadata),
    )


def run_robustness_suite(
    cfg: ReadoutStudyConfig,
    representative: dict[str, ReplayEvaluation],
) -> dict[str, dict[str, Any]]:
    robustness: dict[str, dict[str, Any]] = {}
    for family in ("square", "procedural_segments", "nulling_tail", "fourier_basis"):
        nominal = representative[family]
        family_records: list[dict[str, Any]] = []
        for scale in cfg.chi_probe_scale:
            cfg_shift = replace(cfg, chi=cfg.chi * scale)
            ev = evaluate_full_design(nominal.design, cfg_shift)
            family_records.append({"kind": "chi_scale", "value": float(scale), "fidelity_eta": float(ev.metrics.fidelity_eta)})
        for scale in cfg.kappa_probe_scale:
            cfg_shift = replace(cfg, kappa=cfg.kappa * scale)
            ev = evaluate_full_design(nominal.design, cfg_shift)
            family_records.append({"kind": "kappa_scale", "value": float(scale), "fidelity_eta": float(ev.metrics.fidelity_eta)})
        for scale in cfg.amp_probe_scale:
            ev = evaluate_full_design(_scaled_design(nominal.design, amp_scale=scale), cfg)
            family_records.append({"kind": "amp_scale", "value": float(scale), "fidelity_eta": float(ev.metrics.fidelity_eta)})
        for shift_mhz in cfg.detuning_probe_mhz:
            ev = evaluate_full_design(_scaled_design(nominal.design, delta_shift=2.0 * np.pi * shift_mhz * 1.0e6), cfg)
            family_records.append({"kind": "detuning_mhz", "value": float(shift_mhz), "fidelity_eta": float(ev.metrics.fidelity_eta)})
        for phase_deg in cfg.phase_probe_deg:
            ev = evaluate_full_design(_scaled_design(nominal.design, phase_shift=np.deg2rad(phase_deg)), cfg)
            family_records.append({"kind": "phase_deg", "value": float(phase_deg), "fidelity_eta": float(ev.metrics.fidelity_eta)})
        fidelity_values = np.array([record["fidelity_eta"] for record in family_records], dtype=float)
        robustness[family] = {
            "nominal_fidelity_eta": float(nominal.metrics.fidelity_eta),
            "mean_fidelity_eta": float(np.mean(fidelity_values)),
            "worst_fidelity_eta": float(np.min(fidelity_values)),
            "records": family_records,
        }
    return robustness


def run_convergence_checks(
    cfg: ReadoutStudyConfig,
    representative: dict[str, ReplayEvaluation],
) -> dict[str, Any]:
    family = "procedural_segments"
    nominal = representative[family]
    cfg_fine = replace(cfg, dt=0.5 * cfg.dt)
    set_nulling_tail_kappa(cfg_fine.kappa)
    fine_design = get_family(nominal.design.family).builder(
        np.asarray(nominal.design.params, dtype=float),
        nominal.design.duration,
        cfg_fine.dt,
        cfg_fine.amp_max,
        cfg_fine.chi,
    )
    fine_eval = evaluate_full_design(
        fine_design,
        cfg_fine,
    )
    cfg_trunc = replace(cfg, n_cav=cfg.truncation_probe[-1])
    trunc_eval = evaluate_full_design(nominal.design, cfg_trunc)
    return {
        "dt_ref_ns": float(cfg.dt * 1.0e9),
        "dt_fine_ns": float(cfg_fine.dt * 1.0e9),
        "fidelity_dt_delta": float(abs(nominal.metrics.fidelity_eta - fine_eval.metrics.fidelity_eta)),
        "qnd_dt_delta": float(abs(nominal.metrics.qnd_preservation - fine_eval.metrics.qnd_preservation)),
        "residual_dt_delta": float(abs(nominal.metrics.residual_photons - fine_eval.metrics.residual_photons)),
        "n_cav_ref": int(cfg.n_cav),
        "n_cav_probe": int(cfg_trunc.n_cav),
        "fidelity_trunc_delta": float(abs(nominal.metrics.fidelity_eta - trunc_eval.metrics.fidelity_eta)),
        "qnd_trunc_delta": float(abs(nominal.metrics.qnd_preservation - trunc_eval.metrics.qnd_preservation)),
        "residual_trunc_delta": float(abs(nominal.metrics.residual_photons - trunc_eval.metrics.residual_photons)),
    }


def export_representative_traces(
    representative: dict[str, ReplayEvaluation],
    *,
    output_path: Path,
) -> None:
    payload: dict[str, Any] = {}
    for family, evaluation in representative.items():
        payload[f"{family}_t_ns"] = 1.0e9 * evaluation.t_pulse
        payload[f"{family}_alpha_g"] = evaluation.alpha_g
        payload[f"{family}_alpha_e"] = evaluation.alpha_e
        payload[f"{family}_signal_g"] = evaluation.first_signal_g
        payload[f"{family}_signal_e"] = evaluation.first_signal_e
        payload[f"{family}_signal_g_second"] = evaluation.second_signal_g
        payload[f"{family}_signal_e_second"] = evaluation.second_signal_e
    np.savez_compressed(output_path, **payload)


def save_results(
    results: dict[str, Any],
    *,
    output_json: Path,
) -> None:
    output_json.write_text(json.dumps(_to_serializable(results), indent=2))


def run_all_experiments(cfg: ReadoutStudyConfig = DEFAULT_CONFIG) -> dict[str, Any]:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    linear_results = run_linear_sweep(cfg)
    full_results = run_full_sweep(cfg, linear_results)
    hierarchy = compute_bound_hierarchy(cfg, linear_results, full_results)
    representative = representative_catalog(cfg, linear_results, full_results)
    tradeoff_slice = run_tradeoff_slice(cfg, linear_results, full_results)
    robustness = run_robustness_suite(cfg, representative)
    convergence = run_convergence_checks(cfg, representative)
    export_representative_traces(representative, output_path=DATA_DIR / "representative_traces.npz")

    summary = {
        "config": cfg.as_dict(),
        "linear_results": {
            family: [outcome_to_record(outcomes[float(duration)]) for duration in cfg.duration_grid]
            for family, outcomes in linear_results.items()
        },
        "full_results": {
            family: [outcome_to_record(outcomes[float(duration)]) for duration in cfg.duration_grid]
            for family, outcomes in full_results.items()
        },
        "hierarchy": hierarchy,
        "representative": {family: evaluation_to_record(evaluation) for family, evaluation in representative.items()},
        "tradeoff_slice": tradeoff_slice,
        "robustness": robustness,
        "convergence": convergence,
    }
    save_results(summary, output_json=DATA_DIR / "study_summary.json")
    return summary
