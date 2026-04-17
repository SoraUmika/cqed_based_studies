"""Regenerate the initial-state comparison artifacts from the saved selected settings."""

from __future__ import annotations

import json
from pathlib import Path

try:
    from . import run_study as rs
except ImportError:
    import run_study as rs


def main() -> None:
    study_results_path = rs.DATA_DIR / "study_results.json"
    results = json.loads(study_results_path.read_text(encoding="utf-8"))

    model = rs.build_model()
    frame = rs.build_frame(model)
    recommendations = rs.pulsed_recommendations()
    pulsed_best = results["selected_schemes"]["pulsed_ladder"]
    resonant_best = results["selected_schemes"]["continuous_bright"]
    detuned_best = results["selected_schemes"]["continuous_raman"]
    benchmark_best = results["selected_schemes"]["effective_autonomous_benchmark"]
    comparison_window_s = float(results["comparison_window_ns"]) * 1.0e-9

    initial_state_rows: list[dict[str, object]] = []
    for state_key, payload in rs._scheme_initial_states(model).items():
        state = payload["state"]
        pulsed_max_n = int(payload["pulsed_max_n"])
        rows = [
            rs.simulate_pulsed_protocol(
                model,
                frame,
                initial_state=state,
                noise=rs.baseline_noise(),
                recommendations=recommendations,
                ringdown_multiple=float(pulsed_best["ringdown_multiple"]),
                max_n=pulsed_max_n,
                dt_s=rs.DEFAULT_TRAJECTORY_DT_S,
            ),
            rs.simulate_continuous_protocol(
                model,
                frame,
                initial_state=state,
                noise=rs.baseline_noise(),
                target_coupling_mhz=float(resonant_best["target_coupling_mhz"]),
                common_detuning_mhz=0.0,
                duration_s=comparison_window_s,
                dt_s=rs.DEFAULT_TRAJECTORY_DT_S,
            ),
            rs.simulate_continuous_protocol(
                model,
                frame,
                initial_state=state,
                noise=rs.baseline_noise(),
                target_coupling_mhz=float(detuned_best["target_coupling_mhz"]),
                common_detuning_mhz=float(detuned_best["common_detuning_mhz"]),
                duration_s=comparison_window_s,
                dt_s=rs.DEFAULT_TRAJECTORY_DT_S,
            ),
            rs.simulate_effective_autonomous_benchmark(
                model,
                initial_state=state,
                duration_s=comparison_window_s,
                gamma_eff_hz=float(benchmark_best["gamma_eff_hz"]),
                dt_s=rs.DEFAULT_TRAJECTORY_DT_S,
            ),
        ]
        for row in rows:
            row["initial_state_key"] = state_key
            row["initial_state_label"] = payload["label"]
            row["pulsed_max_n"] = pulsed_max_n
            initial_state_rows.append(row)

    rs.csv_dump(rs.DATA_DIR / "initial_state_summary.csv", initial_state_rows)
    rs.json_dump(rs.ARTIFACTS_DIR / "initial_state_summary.json", {"rows": initial_state_rows})
    rs._make_initial_state_figure(initial_state_rows)

    print(f"Wrote {rs.DATA_DIR / 'initial_state_summary.csv'}")
    print(f"Wrote {rs.ARTIFACTS_DIR / 'initial_state_summary.json'}")


if __name__ == "__main__":
    main()