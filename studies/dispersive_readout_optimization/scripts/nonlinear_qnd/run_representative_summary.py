"""Run the representative-duration diagnostics subset of the follow-up study."""

from __future__ import annotations

from dataclasses import replace

from readout_opt import DEFAULT_CONFIG
from readout_opt.config import DATA_DIR
from readout_opt.experiments import run_all_experiments, save_results


def main() -> None:
    cfg = replace(DEFAULT_CONFIG, duration_grid_ns=(240.0,), representative_duration_ns=240.0)
    print("representative: run", flush=True)
    summary = run_all_experiments(cfg)
    save_results(summary, output_json=DATA_DIR / "representative_summary.json")
    print("representative: done", flush=True)


if __name__ == "__main__":
    main()