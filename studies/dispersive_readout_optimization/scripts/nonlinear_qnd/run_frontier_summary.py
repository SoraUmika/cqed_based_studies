"""Run the stable multi-duration frontier subset of the follow-up study."""

from __future__ import annotations

from dataclasses import replace

from readout_opt import DEFAULT_CONFIG
from readout_opt.config import DATA_DIR
from readout_opt.experiments import (
    compute_bound_hierarchy,
    evaluation_to_record,
    outcome_to_record,
    run_linear_sweep,
    run_nominal_rich_replay,
    run_physical_sweep,
    save_results,
)


def main() -> None:
    cfg = replace(DEFAULT_CONFIG, duration_grid_ns=(96.0, 240.0, 496.0), representative_duration_ns=240.0)
    print("frontier: linear", flush=True)
    linear = run_linear_sweep(cfg)
    print("frontier: full", flush=True)
    full = run_physical_sweep(cfg, linear, regime="full")
    print("frontier: rich", flush=True)
    rich = run_physical_sweep(cfg, linear, regime="rich")
    print("frontier: nominal-rich", flush=True)
    nominal = run_nominal_rich_replay(cfg, full)
    print("frontier: hierarchy", flush=True)
    hierarchy = compute_bound_hierarchy(cfg, linear, full, nominal, rich)
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
    }
    save_results(summary, output_json=DATA_DIR / "frontier_summary.json")
    print("frontier: done", flush=True)


if __name__ == "__main__":
    main()