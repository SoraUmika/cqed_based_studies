"""Build the reproducibility notebook for the strong-validation study."""

from __future__ import annotations

import json
from pathlib import Path


STUDY_DIR = Path(__file__).resolve().parent.parent
NOTEBOOK_PATH = STUDY_DIR / "scripts" / "reproducibility_notebook.ipynb"
SUMMARY_PATH = STUDY_DIR / "data" / "study_summary.json"


def md_cell(text: str) -> dict:
    return {"cell_type": "markdown", "metadata": {}, "source": [line + "\n" for line in text.strip().splitlines()]}


def code_cell(text: str) -> dict:
    return {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": [line + "\n" for line in text.strip().splitlines()],
    }


def main() -> None:
    summary = json.loads(SUMMARY_PATH.read_text(encoding="utf-8"))
    best_strict = summary["best_strict_case"]
    best_relaxed = summary["best_relaxed_case"]

    cells = [
        md_cell(
            f"""
            # Strong Validation of SQR / CPSQR for Arbitrary Fock-Conditional Qubit Rotations

            This notebook reproduces the saved results from the study at `{STUDY_DIR}`.

            Main goals:
            1. inspect the strict and relaxed highlight cases,
            2. reload saved artifacts and summary tables,
            3. regenerate the key comparison plots from the saved CSV,
            4. provide commented live-rerun cells for selected cases.
            """
        ),
        md_cell(
            """
            ## Environment Setup

            This cell configures imports and makes the study `scripts/` directory importable.
            The notebook is artifact-first: loading saved JSON/CSV is the default path, while expensive recomputation is left in commented cells below.
            """
        ),
        code_cell(
            f"""
            from pathlib import Path
            import json
            import sys

            import matplotlib.pyplot as plt
            import numpy as np
            import pandas as pd

            STUDY_DIR = Path(r\"{STUDY_DIR}\")
            SCRIPTS_DIR = STUDY_DIR / "scripts"
            if str(SCRIPTS_DIR) not in sys.path:
                sys.path.insert(0, str(SCRIPTS_DIR))

            import run_study
            import study_lib as lib

            SUMMARY_PATH = STUDY_DIR / "data" / "study_summary.json"
            RESULTS_CSV = STUDY_DIR / "data" / "study_results.csv"
            ARTIFACT_DIR = STUDY_DIR / "artifacts" / "cases"

            print("Study directory:", STUDY_DIR)
            """
        ),
        md_cell(
            """
            ## User-Tunable Parameters

            Edit this single cell if you want to switch the representative case, compare a different target family, or rerun a selected optimization.
            The values here are the only knobs that downstream cells assume.
            """
        ),
        code_cell(
            f"""
            # Highlight artifacts loaded by default.
            STRICT_CASE_ID = "{best_strict['case_id']}"
            STRICT_FAMILY = "{best_strict['family_name']}"

            RELAXED_CASE_ID = "{best_relaxed['case_id']}"
            RELAXED_FAMILY = "{best_relaxed['family_name']}"

            # Optional live-rerun request.
            SELECT_STAGE = "main"
            SELECT_TARGET_FAMILY = "structured_zyz"
            SELECT_MODEL_VARIANT = "chi_plus_chiprime"
            SELECT_INCLUDE_CHI_PRIME = True
            SELECT_N_ACTIVE = 3
            SELECT_CHI_T_OVER_2PI = 5.0
            SELECT_FAMILY_NAMES = ("single_pulse_gaussian",)
            SELECT_RANDOM_SEED = None

            print("Strict default:", STRICT_CASE_ID, STRICT_FAMILY)
            print("Relaxed default:", RELAXED_CASE_ID, RELAXED_FAMILY)
            print("Live request:", SELECT_TARGET_FAMILY, SELECT_MODEL_VARIANT, SELECT_N_ACTIVE, SELECT_CHI_T_OVER_2PI)
            """
        ),
        md_cell(
            """
            ## Derived Objects

            This cell derives the artifact paths and a `CaseRequest` object from the tunable parameters above.
            If you change the tunable parameters, rerun this cell before the live-rerun cells later in the notebook.
            """
        ),
        code_cell(
            """
            summary = json.loads(SUMMARY_PATH.read_text(encoding="utf-8"))
            results_df = pd.read_csv(RESULTS_CSV)

            strict_artifact_path = ARTIFACT_DIR / f"{STRICT_CASE_ID}_{STRICT_FAMILY}.json"
            relaxed_artifact_path = ARTIFACT_DIR / f"{RELAXED_CASE_ID}_{RELAXED_FAMILY}.json"

            selected_request = run_study.CaseRequest(
                stage=SELECT_STAGE,
                target_family=SELECT_TARGET_FAMILY,
                model_variant=SELECT_MODEL_VARIANT,
                include_chi_prime=SELECT_INCLUDE_CHI_PRIME,
                n_active=SELECT_N_ACTIVE,
                chi_t_over_2pi=SELECT_CHI_T_OVER_2PI,
                family_names=tuple(SELECT_FAMILY_NAMES),
                random_seed=SELECT_RANDOM_SEED,
            )

            print("Strict artifact:", strict_artifact_path.name)
            print("Relaxed artifact:", relaxed_artifact_path.name)
            print("Selected request:", selected_request)
            """
        ),
        md_cell(
            """
            ## Load Saved Study Summary

            This is the fastest overview path. It loads the saved machine-readable summary and the full result table.
            Tunable parameters do not affect this cell unless you point the notebook at a different study directory.
            """
        ),
        code_cell(
            """
            display(pd.DataFrame(summary["family_summary"]))
            display(results_df.head())
            """
        ),
        md_cell(
            """
            ## Load Saved Strict Highlight Case

            This cell loads the saved artifact for the best strict case and shows the headline metrics, target metadata, and reduced six-state validation rows.
            """
        ),
        code_cell(
            """
            strict_artifact = json.loads(strict_artifact_path.read_text(encoding="utf-8"))
            strict_row = strict_artifact["summary_row"]
            print("Strict headline metrics:")
            for key in ("strict_joint_process_fidelity", "relaxed_joint_process_fidelity", "classification"):
                print(f"  {key}: {strict_row[key]}")
            display(pd.DataFrame(strict_artifact["target_spec"]["block_rows"]))
            display(pd.DataFrame(strict_artifact["reduced_probe_rows"]).head(12))
            """
        ),
        md_cell(
            """
            ## Re-run Strict Highlight Case

            This is the slow path. Uncomment and run the cell below only if you want to rebuild the selected request live with the current package.
            The tunable parameters above control which case is rebuilt.
            """
        ),
        code_cell(
            """
            # selected_context = run_study.build_case_context(selected_request)
            # live_row, live_artifact, _ = run_study.evaluate_family(selected_context, SELECT_FAMILY_NAMES[0])
            # print(live_row["strict_joint_process_fidelity"], live_row["relaxed_joint_process_fidelity"], live_row["classification"])
            """
        ),
        md_cell(
            """
            ## Load Saved Relaxed Highlight Case

            This cell loads the saved artifact for the best relaxed case and shows the same main diagnostics.
            """
        ),
        code_cell(
            """
            relaxed_artifact = json.loads(relaxed_artifact_path.read_text(encoding="utf-8"))
            relaxed_row = relaxed_artifact["summary_row"]
            print("Relaxed headline metrics:")
            for key in ("strict_joint_process_fidelity", "relaxed_joint_process_fidelity", "classification"):
                print(f"  {key}: {relaxed_row[key]}")
            display(pd.DataFrame(relaxed_artifact["relaxed_fit_rows"]))
            display(pd.DataFrame(relaxed_artifact["cross_block_rows"]))
            """
        ),
        md_cell(
            """
            ## Re-run Relaxed Highlight Case

            This is the live recomputation path for a relaxed-target family. Uncomment only when you want to repeat the expensive optimization.
            """
        ),
        code_cell(
            """
            # relaxed_request = run_study.CaseRequest(
            #     stage="stress",
            #     target_family="stress_zyz",
            #     model_variant="chi_plus_chiprime",
            #     include_chi_prime=True,
            #     n_active=3,
            #     chi_t_over_2pi=5.0,
            #     family_names=("segmented_relaxed",),
            # )
            # relaxed_context = run_study.build_case_context(relaxed_request)
            # live_row, live_artifact, _ = run_study.evaluate_family(relaxed_context, "segmented_relaxed")
            # print(live_row["strict_joint_process_fidelity"], live_row["relaxed_joint_process_fidelity"], live_row["classification"])
            """
        ),
        md_cell(
            """
            ## Validation and Key Figures

            These cells reproduce the main family-comparison and strict-vs-relaxed plots from the saved CSV, not from a fresh optimization run.
            """
        ),
        code_cell(
            """
            family_summary = (
                results_df.groupby("family_name", as_index=False)[
                    ["strict_joint_process_fidelity", "relaxed_joint_process_fidelity"]
                ]
                .mean()
                .sort_values("relaxed_joint_process_fidelity", ascending=False)
            )

            fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
            axes[0].bar(family_summary["family_name"], family_summary["strict_joint_process_fidelity"])
            axes[0].set_title("Mean strict joint process")
            axes[0].tick_params(axis="x", rotation=25)

            axes[1].scatter(results_df["strict_joint_process_fidelity"], results_df["relaxed_joint_process_fidelity"], s=16, alpha=0.7)
            axes[1].set_xlabel("Strict joint")
            axes[1].set_ylabel("Relaxed joint")
            axes[1].set_title("Strict vs relaxed")
            plt.tight_layout()
            plt.show()
            """
        ),
        md_cell(
            """
            ## Summary

            This final cell prints a compact parameter map. If you change the tunable parameters, rerun the parameter and derived-object cells first.
            """
        ),
        code_cell(
            """
            param_rows = pd.DataFrame(
                [
                    {"parameter": "STRICT_CASE_ID", "value": STRICT_CASE_ID, "effect": "Loads the default strict highlight artifact."},
                    {"parameter": "STRICT_FAMILY", "value": STRICT_FAMILY, "effect": "Chooses the strict-highlight control family."},
                    {"parameter": "RELAXED_CASE_ID", "value": RELAXED_CASE_ID, "effect": "Loads the default relaxed highlight artifact."},
                    {"parameter": "RELAXED_FAMILY", "value": RELAXED_FAMILY, "effect": "Chooses the relaxed-highlight control family."},
                    {"parameter": "SELECT_TARGET_FAMILY", "value": SELECT_TARGET_FAMILY, "effect": "Sets the live rerun target family."},
                    {"parameter": "SELECT_MODEL_VARIANT", "value": SELECT_MODEL_VARIANT, "effect": "Selects chi-only or chi+chi' for live reruns."},
                    {"parameter": "SELECT_N_ACTIVE", "value": SELECT_N_ACTIVE, "effect": "Changes the number of addressed Fock manifolds in the live request."},
                    {"parameter": "SELECT_CHI_T_OVER_2PI", "value": SELECT_CHI_T_OVER_2PI, "effect": "Changes the live pulse-duration setting."},
                    {"parameter": "SELECT_FAMILY_NAMES", "value": SELECT_FAMILY_NAMES, "effect": "Chooses the live control family to rerun."},
                    {"parameter": "SELECT_RANDOM_SEED", "value": SELECT_RANDOM_SEED, "effect": "Sets the random target seed when applicable."},
                ]
            )
            display(param_rows)
            """
        ),
    ]

    notebook = {
        "cells": cells,
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python", "version": "3.12"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    NOTEBOOK_PATH.write_text(json.dumps(notebook, indent=2), encoding="utf-8")
    print(f"Wrote {NOTEBOOK_PATH}")


if __name__ == "__main__":
    main()
