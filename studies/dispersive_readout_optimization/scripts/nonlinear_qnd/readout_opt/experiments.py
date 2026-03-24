"""Orchestration for the nonlinear-QND and hardware-realistic readout follow-up study."""

from __future__ import annotations

from dataclasses import asdict, replace
import json
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .config import DATA_DIR, DEFAULT_CONFIG, ReadoutStudyConfig
from .optimize import OptimizationOutcome, optimize_family
from .pulse_families import PulseDesign, get_family, set_nulling_tail_kappa
from .simulate import (
    ReplayEvaluation,
    evaluate_full_design,
    evaluate_hardware_design,
    evaluate_multilevel_design,
    evaluate_nonlinear_design,
    evaluate_rich_design,
    transport_analysis,
)

LINEAR_FAMILIES = (
    "square",
    "smooth_square",
    "ring_hold",
    "procedural_segments",
    "nulling_tail",
    "fourier_basis",
    "piecewise_reference",
)

PRACTICAL_FAMILIES = (
    "square",
    "smooth_square",
    "ring_hold",
    "procedural_segments",
    "nulling_tail",
    "fourier_basis",
)

KEY_FAMILIES = (
    "square",
    "procedural_segments",
    "nulling_tail",
    "fourier_basis",
)

REPRESENTATIVE_FAMILIES = KEY_FAMILIES + ("piecewise_reference",)


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


def run_linear_sweep(cfg: ReadoutStudyConfig) -> dict[str, dict[float, OptimizationOutcome]]:
    results: dict[str, dict[float, OptimizationOutcome]] = {family: {} for family in LINEAR_FAMILIES}
    for duration in cfg.duration_grid:
        duration_ns = int(round(duration * 1.0e9))
        for family in LINEAR_FAMILIES:
            print(f"[linear] {family} @ {duration_ns} ns", flush=True)
            outcome = optimize_family(
                family=family,
                regime="linear",
                objective_name="info" if family == "piecewise_reference" else "balanced",
                duration=float(duration),
                cfg=cfg,
                seed=cfg.seed + duration_ns + len(family),
                n_random=18 if family == "piecewise_reference" else 12,
                n_local=2 if family == "piecewise_reference" else 1,
                maxiter=28 if family == "piecewise_reference" else 22,
            )
            results[family][float(duration)] = outcome
    return results


def run_physical_sweep(
    cfg: ReadoutStudyConfig,
    linear_results: dict[str, dict[float, OptimizationOutcome]],
    *,
    regime: str,
) -> dict[str, dict[float, OptimizationOutcome]]:
    results: dict[str, dict[float, OptimizationOutcome]] = {family: {} for family in PRACTICAL_FAMILIES}
    regime_offset = {"full": 1000, "rich": 2000}[regime]
    for duration in cfg.duration_grid:
        duration_ns = int(round(duration * 1.0e9))
        for family in PRACTICAL_FAMILIES:
            print(f"[{regime}] {family} @ {duration_ns} ns", flush=True)
            warm_start = linear_results[family][float(duration)].x_best
            outcome = optimize_family(
                family=family,
                regime=regime,
                objective_name="balanced",
                duration=float(duration),
                cfg=cfg,
                warm_start=warm_start,
                seed=cfg.seed + regime_offset + duration_ns + len(family),
                n_random=6,
                n_local=0,
                maxiter=14,
            )
            results[family][float(duration)] = outcome
    return results


def run_nominal_rich_replay(
    cfg: ReadoutStudyConfig,
    full_results: dict[str, dict[float, OptimizationOutcome]],
) -> dict[str, dict[float, ReplayEvaluation]]:
    replays: dict[str, dict[float, ReplayEvaluation]] = {family: {} for family in PRACTICAL_FAMILIES}
    for duration in cfg.duration_grid:
        duration_ns = int(round(duration * 1.0e9))
        for family in PRACTICAL_FAMILIES:
            print(f"[nominal-rich] {family} @ {duration_ns} ns", flush=True)
            design = full_results[family][float(duration)].design
            replays[family][float(duration)] = evaluate_rich_design(design, cfg)
    return replays


def compute_bound_hierarchy(
    cfg: ReadoutStudyConfig,
    linear_results: dict[str, dict[float, OptimizationOutcome]],
    full_results: dict[str, dict[float, OptimizationOutcome]],
    nominal_rich_replay: dict[str, dict[float, ReplayEvaluation]],
    rich_results: dict[str, dict[float, OptimizationOutcome]],
) -> list[dict[str, float | str]]:
    hierarchy: list[dict[str, float | str]] = []
    for duration in cfg.duration_grid:
        duration_key = float(duration)
        ideal = linear_results["piecewise_reference"][duration_key].evaluation.metrics.fidelity_ideal
        detector = linear_results["piecewise_reference"][duration_key].evaluation.metrics.fidelity_eta
        t1_bound = linear_results["piecewise_reference"][duration_key].evaluation.metrics.fidelity_t1_bound

        full_best_family = max(
            PRACTICAL_FAMILIES,
            key=lambda family: full_results[family][duration_key].evaluation.metrics.score_balanced,
        )
        nominal_rich_best_family = max(
            PRACTICAL_FAMILIES,
            key=lambda family: nominal_rich_replay[family][duration_key].metrics.score_balanced,
        )
        rich_best_family = max(
            PRACTICAL_FAMILIES,
            key=lambda family: rich_results[family][duration_key].evaluation.metrics.score_balanced,
        )

        qnd_candidates = [
            rich_results[family][duration_key].evaluation
            for family in PRACTICAL_FAMILIES
            if rich_results[family][duration_key].evaluation.metrics.qnd_preservation >= cfg.benchmark_qnd_min
        ]
        qnd_best = max((candidate.metrics.fidelity_eta for candidate in qnd_candidates), default=float("nan"))

        hierarchy.append(
            {
                "duration_ns": float(duration * 1.0e9),
                "ideal_bound": float(ideal),
                "detector_bound": float(detector),
                "t1_bound": float(t1_bound),
                "legacy_best": float(full_results[full_best_family][duration_key].evaluation.metrics.fidelity_eta),
                "legacy_best_family": full_best_family,
                "nominal_rich_best": float(nominal_rich_replay[nominal_rich_best_family][duration_key].metrics.fidelity_eta),
                "nominal_rich_best_family": nominal_rich_best_family,
                "rich_best": float(rich_results[rich_best_family][duration_key].evaluation.metrics.fidelity_eta),
                "rich_best_family": rich_best_family,
                "rich_qnd_constrained_best": float(qnd_best),
            }
        )
    return hierarchy


def run_representative_regime_breakdown(
    cfg: ReadoutStudyConfig,
    full_results: dict[str, dict[float, OptimizationOutcome]],
    linear_results: dict[str, dict[float, OptimizationOutcome]],
) -> dict[str, dict[str, ReplayEvaluation]]:
    duration = float(cfg.representative_duration)
    evaluators: dict[str, Callable[[PulseDesign, ReadoutStudyConfig], ReplayEvaluation]] = {
        "full": evaluate_full_design,
        "multilevel": evaluate_multilevel_design,
        "hardware": evaluate_hardware_design,
        "nonlinear": evaluate_nonlinear_design,
        "rich": evaluate_rich_design,
    }
    breakdown: dict[str, dict[str, ReplayEvaluation]] = {}
    for family in REPRESENTATIVE_FAMILIES:
        print(f"[breakdown] {family}", flush=True)
        if family == "piecewise_reference":
            design = linear_results[family][duration].design
        else:
            design = full_results[family][duration].design
        breakdown[family] = {name: evaluator(design, cfg) for name, evaluator in evaluators.items()}
    return breakdown


def run_representative_reference_slice(
    cfg: ReadoutStudyConfig,
    linear_results: dict[str, dict[float, OptimizationOutcome]],
) -> dict[str, OptimizationOutcome]:
    duration = float(cfg.representative_duration)
    warm_start = linear_results["piecewise_reference"][duration].x_best
    outcomes: dict[str, OptimizationOutcome] = {}
    for offset, objective in enumerate(("info", "balanced", "emptying"), start=1):
        print(f"[reference] {objective}", flush=True)
        outcomes[objective] = optimize_family(
            family="piecewise_reference",
            regime="rich",
            objective_name=objective,
            duration=duration,
            cfg=cfg,
            warm_start=warm_start,
            seed=cfg.seed + 3000 + offset,
            n_random=8,
            n_local=0,
            maxiter=16,
        )
    return outcomes


def run_tradeoff_slice(
    cfg: ReadoutStudyConfig,
    rich_results: dict[str, dict[float, OptimizationOutcome]],
    reference_slice: dict[str, OptimizationOutcome],
) -> dict[str, dict[str, dict[str, Any]]]:
    duration = float(cfg.representative_duration)
    slice_results: dict[str, dict[str, dict[str, Any]]] = {}
    for family in KEY_FAMILIES:
        print(f"[tradeoff] {family}", flush=True)
        warm_start = rich_results[family][duration].x_best
        slice_results[family] = {}
        for objective in ("info", "balanced", "emptying"):
            if objective == "balanced":
                outcome = rich_results[family][duration]
            else:
                outcome = optimize_family(
                    family=family,
                    regime="rich",
                    objective_name=objective,
                    duration=duration,
                    cfg=cfg,
                    warm_start=warm_start,
                    seed=cfg.seed + 3200 + len(slice_results) * 50,
                    n_random=6,
                    n_local=0,
                    maxiter=14,
                )
            slice_results[family][objective] = outcome_to_record(outcome)
    slice_results["piecewise_reference"] = {objective: outcome_to_record(outcome) for objective, outcome in reference_slice.items()}
    return slice_results


def run_qnd_stress_test(
    cfg: ReadoutStudyConfig,
    representative_rich: dict[str, ReplayEvaluation],
) -> dict[str, list[dict[str, Any]]]:
    stress: dict[str, list[dict[str, Any]]] = {}
    for family in KEY_FAMILIES:
        print(f"[qnd-stress] {family}", flush=True)
        design = representative_rich[family].design
        records: list[dict[str, Any]] = []
        for scale in cfg.amp_probe_scale:
            ev = evaluate_rich_design(_scaled_design(design, amp_scale=scale), cfg)
            records.append(
                {
                    "amp_scale": float(scale),
                    "fidelity_eta": float(ev.metrics.fidelity_eta),
                    "qnd_preservation": float(ev.metrics.qnd_preservation),
                    "measurement_induced_transition": float(ev.metrics.measurement_induced_transition),
                    "leakage": float(ev.metrics.leakage),
                    "peak_photons": float(ev.metrics.peak_photons),
                }
            )
        stress[family] = records
    return stress


def run_robustness_suite(
    cfg: ReadoutStudyConfig,
    representative_rich: dict[str, ReplayEvaluation],
) -> dict[str, dict[str, Any]]:
    robustness: dict[str, dict[str, Any]] = {}
    for family in KEY_FAMILIES:
        print(f"[robustness] {family}", flush=True)
        nominal = representative_rich[family]
        family_records: list[dict[str, Any]] = []
        for scale in cfg.chi_probe_scale:
            cfg_shift = replace(cfg, chi=cfg.chi * scale)
            ev = evaluate_rich_design(nominal.design, cfg_shift)
            family_records.append({"kind": "chi_scale", "value": float(scale), "fidelity_eta": float(ev.metrics.fidelity_eta)})
        for scale in cfg.kappa_probe_scale:
            cfg_shift = replace(cfg, kappa=cfg.kappa * scale)
            ev = evaluate_rich_design(nominal.design, cfg_shift)
            family_records.append({"kind": "kappa_scale", "value": float(scale), "fidelity_eta": float(ev.metrics.fidelity_eta)})
        for scale in cfg.amp_probe_scale:
            ev = evaluate_rich_design(_scaled_design(nominal.design, amp_scale=scale), cfg)
            family_records.append({"kind": "amp_scale", "value": float(scale), "fidelity_eta": float(ev.metrics.fidelity_eta)})
        for shift_mhz in cfg.detuning_probe_mhz:
            ev = evaluate_rich_design(_scaled_design(nominal.design, delta_shift=2.0 * np.pi * shift_mhz * 1.0e6), cfg)
            family_records.append({"kind": "detuning_mhz", "value": float(shift_mhz), "fidelity_eta": float(ev.metrics.fidelity_eta)})
        for phase_deg in cfg.phase_probe_deg:
            ev = evaluate_rich_design(_scaled_design(nominal.design, phase_shift=np.deg2rad(phase_deg)), cfg)
            family_records.append({"kind": "phase_deg", "value": float(phase_deg), "fidelity_eta": float(ev.metrics.fidelity_eta)})
        for bw in cfg.hardware_bandwidth_probe_hz:
            cfg_shift = replace(cfg, hardware=replace(cfg.hardware, lowpass_bw_hz=bw))
            ev = evaluate_rich_design(nominal.design, cfg_shift)
            family_records.append({"kind": "hardware_bw_hz", "value": float(bw), "fidelity_eta": float(ev.metrics.fidelity_eta)})
        for skew_deg in cfg.hardware_skew_probe_deg:
            cfg_shift = replace(cfg, hardware=replace(cfg.hardware, quadrature_skew_deg=skew_deg))
            ev = evaluate_rich_design(nominal.design, cfg_shift)
            family_records.append({"kind": "hardware_skew_deg", "value": float(skew_deg), "fidelity_eta": float(ev.metrics.fidelity_eta)})
        for gain_q in cfg.hardware_gain_q_probe:
            cfg_shift = replace(cfg, hardware=replace(cfg.hardware, gain_q=gain_q))
            ev = evaluate_rich_design(nominal.design, cfg_shift)
            family_records.append({"kind": "hardware_gain_q", "value": float(gain_q), "fidelity_eta": float(ev.metrics.fidelity_eta)})
        for bits in cfg.hardware_bits_probe:
            cfg_shift = replace(cfg, hardware=replace(cfg.hardware, amplitude_bits=int(bits)))
            ev = evaluate_rich_design(nominal.design, cfg_shift)
            family_records.append({"kind": "hardware_bits", "value": int(bits), "fidelity_eta": float(ev.metrics.fidelity_eta)})
        combined_cfg = replace(
            cfg,
            chi=cfg.chi * cfg.chi_probe_scale[0],
            hardware=replace(
                cfg.hardware,
                lowpass_bw_hz=min(cfg.hardware_bandwidth_probe_hz),
                gain_q=min(cfg.hardware_gain_q_probe),
                quadrature_skew_deg=max(cfg.hardware_skew_probe_deg),
            ),
        )
        combined_ev = evaluate_rich_design(_scaled_design(nominal.design, amp_scale=min(cfg.amp_probe_scale)), combined_cfg)
        family_records.append({"kind": "combined_worst_case", "value": 1.0, "fidelity_eta": float(combined_ev.metrics.fidelity_eta)})
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
    representative_rich: dict[str, ReplayEvaluation],
) -> dict[str, Any]:
    print("[convergence] procedural_segments", flush=True)
    family = "procedural_segments"
    nominal = representative_rich[family]
    cfg_fine = replace(cfg, dt=0.5 * cfg.dt)
    set_nulling_tail_kappa(cfg_fine.kappa)
    fine_design = get_family(nominal.design.family).builder(
        np.asarray(nominal.design.params, dtype=float),
        nominal.design.duration,
        cfg_fine.dt,
        cfg_fine.amp_max,
        cfg_fine.chi,
    )
    fine_eval = evaluate_rich_design(fine_design, cfg_fine)
    cfg_trunc = replace(cfg, n_cav=cfg.truncation_probe[-1])
    trunc_eval = evaluate_rich_design(nominal.design, cfg_trunc)
    return {
        "dt_ref_ns": float(cfg.dt * 1.0e9),
        "dt_fine_ns": float(cfg_fine.dt * 1.0e9),
        "fidelity_dt_delta": float(abs(nominal.metrics.fidelity_eta - fine_eval.metrics.fidelity_eta)),
        "qnd_dt_delta": float(abs(nominal.metrics.qnd_preservation - fine_eval.metrics.qnd_preservation)),
        "transition_dt_delta": float(abs(nominal.metrics.measurement_induced_transition - fine_eval.metrics.measurement_induced_transition)),
        "residual_dt_delta": float(abs(nominal.metrics.residual_photons - fine_eval.metrics.residual_photons)),
        "n_cav_ref": int(cfg.n_cav),
        "n_cav_probe": int(cfg_trunc.n_cav),
        "fidelity_trunc_delta": float(abs(nominal.metrics.fidelity_eta - trunc_eval.metrics.fidelity_eta)),
        "qnd_trunc_delta": float(abs(nominal.metrics.qnd_preservation - trunc_eval.metrics.qnd_preservation)),
        "transition_trunc_delta": float(abs(nominal.metrics.measurement_induced_transition - trunc_eval.metrics.measurement_induced_transition)),
        "residual_trunc_delta": float(abs(nominal.metrics.residual_photons - trunc_eval.metrics.residual_photons)),
    }


def export_representative_traces(
    representative_rich: dict[str, ReplayEvaluation],
    representative_breakdown: dict[str, dict[str, ReplayEvaluation]],
    representative_reference: dict[str, OptimizationOutcome],
    *,
    cfg: ReadoutStudyConfig,
    output_path: Path,
) -> None:
    payload: dict[str, Any] = {}
    for family, evaluation in representative_rich.items():
        payload[f"{family}_t_ns"] = 1.0e9 * evaluation.t_pulse
        payload[f"{family}_alpha_g"] = evaluation.alpha_g
        payload[f"{family}_alpha_e"] = evaluation.alpha_e
        payload[f"{family}_signal_g"] = evaluation.first_signal_g
        payload[f"{family}_signal_e"] = evaluation.first_signal_e
        payload[f"{family}_signal_g_second"] = evaluation.second_signal_g
        payload[f"{family}_signal_e_second"] = evaluation.second_signal_e
        transport = transport_analysis(evaluation.design, cfg, regime="rich")
        payload[f"{family}_program_waveform"] = transport["program_waveform"]
        payload[f"{family}_transport_waveform"] = transport["transport_waveform"]
        payload[f"{family}_distorted_waveform"] = transport["distorted_waveform"]
        nominal_rich = representative_breakdown[family]["rich"]
        payload[f"{family}_nominal_rich_signal_g"] = nominal_rich.first_signal_g
        payload[f"{family}_nominal_rich_signal_e"] = nominal_rich.first_signal_e
    for objective, outcome in representative_reference.items():
        payload[f"piecewise_reference_{objective}_signal_g"] = outcome.evaluation.first_signal_g
        payload[f"piecewise_reference_{objective}_signal_e"] = outcome.evaluation.first_signal_e
        payload[f"piecewise_reference_{objective}_alpha_g"] = outcome.evaluation.alpha_g
        payload[f"piecewise_reference_{objective}_alpha_e"] = outcome.evaluation.alpha_e
        transport = transport_analysis(outcome.design, cfg, regime="rich")
        payload[f"piecewise_reference_{objective}_program_waveform"] = transport["program_waveform"]
        payload[f"piecewise_reference_{objective}_distorted_waveform"] = transport["distorted_waveform"]
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
    full_results = run_physical_sweep(cfg, linear_results, regime="full")
    rich_results = run_physical_sweep(cfg, linear_results, regime="rich")
    nominal_rich_replay = run_nominal_rich_replay(cfg, full_results)
    hierarchy = compute_bound_hierarchy(cfg, linear_results, full_results, nominal_rich_replay, rich_results)
    representative_breakdown = run_representative_regime_breakdown(cfg, full_results, linear_results)
    representative_reference = run_representative_reference_slice(cfg, linear_results)
    representative_rich = {
        family: rich_results[family][float(cfg.representative_duration)].evaluation for family in KEY_FAMILIES
    }
    representative_rich["piecewise_reference"] = representative_reference["balanced"].evaluation
    tradeoff_slice = run_tradeoff_slice(cfg, rich_results, representative_reference)
    qnd_stress = run_qnd_stress_test(cfg, representative_rich)
    robustness = run_robustness_suite(cfg, representative_rich)
    convergence = run_convergence_checks(cfg, representative_rich)
    export_representative_traces(
        representative_rich,
        representative_breakdown,
        representative_reference,
        cfg=cfg,
        output_path=DATA_DIR / "representative_traces.npz",
    )

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
        "rich_results": {
            family: [outcome_to_record(outcomes[float(duration)]) for duration in cfg.duration_grid]
            for family, outcomes in rich_results.items()
        },
        "nominal_rich_replay": {
            family: [evaluation_to_record(outcomes[float(duration)]) for duration in cfg.duration_grid]
            for family, outcomes in nominal_rich_replay.items()
        },
        "hierarchy": hierarchy,
        "representative_breakdown": {
            family: {regime: evaluation_to_record(evaluation) for regime, evaluation in regimes.items()}
            for family, regimes in representative_breakdown.items()
        },
        "representative_reference": {objective: outcome_to_record(outcome) for objective, outcome in representative_reference.items()},
        "representative_rich": {family: evaluation_to_record(evaluation) for family, evaluation in representative_rich.items()},
        "tradeoff_slice": tradeoff_slice,
        "qnd_stress": qnd_stress,
        "robustness": robustness,
        "convergence": convergence,
    }
    save_results(summary, output_json=DATA_DIR / "study_summary.json")
    return summary
