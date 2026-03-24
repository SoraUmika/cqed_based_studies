"""Run the follow-up study with stage markers for long executions."""

from __future__ import annotations

from dataclasses import replace

from readout_opt import DEFAULT_CONFIG
from readout_opt.config import DATA_DIR
from readout_opt.experiments import (
    compute_bound_hierarchy,
    evaluation_to_record,
    export_representative_traces,
    outcome_to_record,
    run_convergence_checks,
    run_linear_sweep,
    run_nominal_rich_replay,
    run_physical_sweep,
    run_qnd_stress_test,
    run_representative_reference_slice,
    run_representative_regime_breakdown,
    run_robustness_suite,
    run_tradeoff_slice,
    save_results,
)


def main() -> None:
    cfg = replace(DEFAULT_CONFIG, duration_grid_ns=(96.0, 240.0, 496.0), representative_duration_ns=240.0)
    print("stage: linear", flush=True)
    linear = run_linear_sweep(cfg)
    print("stage: full", flush=True)
    full = run_physical_sweep(cfg, linear, regime="full")
    print("stage: rich", flush=True)
    rich = run_physical_sweep(cfg, linear, regime="rich")
    print("stage: nominal-rich", flush=True)
    nominal = run_nominal_rich_replay(cfg, full)
    print("stage: hierarchy", flush=True)
    hierarchy = compute_bound_hierarchy(cfg, linear, full, nominal, rich)
    print("stage: breakdown", flush=True)
    breakdown = run_representative_regime_breakdown(cfg, full, linear)
    print("stage: reference", flush=True)
    reference = run_representative_reference_slice(cfg, linear)
    representative_rich = {
        family: rich[family][float(cfg.representative_duration)].evaluation
        for family in ("square", "procedural_segments", "nulling_tail", "fourier_basis")
    }
    representative_rich["piecewise_reference"] = reference["balanced"].evaluation
    print("stage: tradeoff", flush=True)
    tradeoff = run_tradeoff_slice(cfg, rich, reference)
    print("stage: qnd-stress", flush=True)
    qnd_stress = run_qnd_stress_test(cfg, representative_rich)
    print("stage: robustness", flush=True)
    robustness = run_robustness_suite(cfg, representative_rich)
    print("stage: convergence", flush=True)
    convergence = run_convergence_checks(cfg, representative_rich)
    print("stage: export", flush=True)
    export_representative_traces(
        representative_rich,
        breakdown,
        reference,
        cfg=cfg,
        output_path=DATA_DIR / "representative_traces.npz",
    )
    print("stage: save", flush=True)
    summary = {
        "config": cfg.as_dict(),
        "linear_results": {
            family: [outcome_to_record(outcomes[float(duration)]) for duration in cfg.duration_grid]
            for family, outcomes in linear.items()
        },
        "full_results": {
            family: [outcome_to_record(outcomes[float(duration)]) for duration in cfg.duration_grid]
            for family, outcomes in full.items()
        },
        "rich_results": {
            family: [outcome_to_record(outcomes[float(duration)]) for duration in cfg.duration_grid]
            for family, outcomes in rich.items()
        },
        "nominal_rich_replay": {
            family: [evaluation_to_record(outcomes[float(duration)]) for duration in cfg.duration_grid]
            for family, outcomes in nominal.items()
        },
        "hierarchy": hierarchy,
        "representative_breakdown": {
            family: {regime: evaluation_to_record(evaluation) for regime, evaluation in regimes.items()}
            for family, regimes in breakdown.items()
        },
        "representative_reference": {objective: outcome_to_record(outcome) for objective, outcome in reference.items()},
        "representative_rich": {family: evaluation_to_record(evaluation) for family, evaluation in representative_rich.items()},
        "tradeoff_slice": tradeoff,
        "qnd_stress": qnd_stress,
        "robustness": robustness,
        "convergence": convergence,
    }
    save_results(summary, output_json=DATA_DIR / "study_summary.json")
    print("stage: done", flush=True)


if __name__ == "__main__":
    main()