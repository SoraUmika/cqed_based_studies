"""Run the full procedural readout study and save the serialized outputs."""

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
    best_row = max(hierarchy, key=lambda row: row["realistic_best"])
    print("Procedural readout study complete.")
    print(
        "Best realistic balanced fidelity: "
        f"{best_row['realistic_best']:.4f} at {best_row['duration_ns']:.0f} ns "
        f"with {best_row['best_family']}."
    )


if __name__ == "__main__":
    main()
