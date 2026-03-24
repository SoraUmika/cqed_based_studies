"""Run the nonlinear-QND readout follow-up study and save the serialized outputs."""

from __future__ import annotations

from readout_opt.config import DATA_DIR, DEFAULT_CONFIG
from readout_opt.experiments import run_all_experiments
from readout_opt.plots import generate_all_figures


def main() -> None:
    summary = run_all_experiments(DEFAULT_CONFIG)
    generate_all_figures(
        summary_path=DATA_DIR / "study_summary.json",
        traces_path=DATA_DIR / "representative_traces.npz",
    )

    hierarchy = summary["hierarchy"]
    best_row = max(hierarchy, key=lambda row: row["rich_best"])
    print("Nonlinear-QND readout follow-up study complete.")
    print(
        "Best rich-model balanced fidelity: "
        f"{best_row['rich_best']:.4f} at {best_row['duration_ns']:.0f} ns "
        f"with {best_row['rich_best_family']}."
    )
    print(
        "Best nominal-rich replay fidelity at the same durations reached "
        f"{max(row['nominal_rich_best'] for row in hierarchy):.4f}."
    )


if __name__ == "__main__":
    main()
