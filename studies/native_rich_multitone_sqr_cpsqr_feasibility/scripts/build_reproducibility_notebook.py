"""Build the reproducibility notebook for the final multitone SQR/CPSQR study."""

from __future__ import annotations

import json
from pathlib import Path


NB_PATH = Path(__file__).resolve().parent / "reproducibility_notebook.ipynb"


def md(text: str) -> dict:
    return {
        "cell_type": "markdown",
        "metadata": {},
        "source": text.strip("\n").splitlines(keepends=True),
    }


def code(text: str) -> dict:
    return {
        "cell_type": "code",
        "metadata": {},
        "execution_count": None,
        "outputs": [],
        "source": text.strip("\n").splitlines(keepends=True),
    }


def main() -> None:
    cells = [
        md(
            """
            # Native / Rich Multitone Feasibility for Ideal SQR and CPSQR

            This notebook reproduces the main saved results of the final patched-package study. It is organized to answer five practical questions:

            1. Does direct multitone realize strict ideal SQR?
            2. Does echoed/composite multitone realize CPSQR?
            3. What changes once the full joint qubit-cavity action is enforced?
            4. What are the shortest convincing durations on the refined grids?
            5. Which saved artifacts should future users inspect first?

            The notebook defaults to loading saved artifacts and figures. Re-run cells are included but commented out because some optimizations are intentionally expensive.
            """
        ),
        md(
            """
            ## Environment Setup

            This cell finds the study root, makes the shared `scripts/` directory importable, and loads the utilities needed to inspect saved study outputs. The tunable parameters in the next section control both which saved artifact we inspect and which live re-run command would be launched if you uncomment the expensive cells.
            """
        ),
        code(
            """
            from __future__ import annotations

            import json
            import sys
            from pathlib import Path

            import matplotlib.pyplot as plt
            import pandas as pd
            from IPython.display import Markdown, Image, display

            candidate_roots = [Path.cwd().resolve(), Path.cwd().resolve().parent, Path.cwd().resolve().parent.parent]
            study_dir = None
            for candidate in candidate_roots:
                if (candidate / "data" / "study_summary.json").exists():
                    study_dir = candidate
                    break
            if study_dir is None:
                raise FileNotFoundError("Could not locate the study root from the current working directory.")

            scripts_dir = study_dir / "scripts"
            data_dir = study_dir / "data"
            figures_dir = study_dir / "figures"
            artifacts_dir = study_dir / "artifacts"

            if str(scripts_dir) not in sys.path:
                sys.path.insert(0, str(scripts_dir))

            from common import CHI, DEFAULT_DT, DEFAULT_SIGMA_FRACTION, PI_PULSE_DURATION_S, build_frame, build_model, duration_from_chi_t, logical_levels, make_run_config, manifold_transition_frequencies_hz, target_spec

            def load_json(path: Path):
                return json.loads(path.read_text(encoding="utf-8-sig"))

            def load_rows(path: Path) -> pd.DataFrame:
                payload = load_json(path)
                rows = payload["rows"] if isinstance(payload, dict) and "rows" in payload else payload
                return pd.DataFrame(rows)

            print(f"Study root: {study_dir}")
            """
        ),
        md(
            """
            ## User-Tunable Parameters

            Every adjustable knob used by this notebook lives in the next code cell. Update these values and then re-run the following derived-object cell to propagate the change through the rest of the notebook. The defaults point to the two most informative representative cases:

            - the best strict direct-SQR case,
            - the clearest echoed CPSQR case whose strict joint fidelity is still poor.
            """
        ),
        code(
            """
            # Core selection knobs
            strict_model_variant = "chi_only"
            strict_target_family = "smooth_x"
            strict_family_name = "reduced_unitary_direct"
            strict_n_active = 2
            strict_chi_t_over_2pi = 5.0

            cpsqr_model_variant = "chi_plus_chiprime"
            cpsqr_target_family = "staggered_x"
            cpsqr_family_name = "echoed_independent"
            cpsqr_n_active = 2
            cpsqr_chi_t_over_2pi = 5.0

            # Model / simulation knobs
            n_tr = 2
            n_cav_padding = 2
            dt_s = DEFAULT_DT
            sigma_fraction = DEFAULT_SIGMA_FRACTION
            inserted_x_pi_duration_s = PI_PULSE_DURATION_S

            # Diagnostic knobs
            reduced_threshold = 0.99
            joint_threshold = 0.99
            duration_grid = [0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0]

            active_config = {
                "strict_case": {
                    "model_variant": strict_model_variant,
                    "target_family": strict_target_family,
                    "family_name": strict_family_name,
                    "n_active": strict_n_active,
                    "chi_t_over_2pi": strict_chi_t_over_2pi,
                },
                "cpsqr_case": {
                    "model_variant": cpsqr_model_variant,
                    "target_family": cpsqr_target_family,
                    "family_name": cpsqr_family_name,
                    "n_active": cpsqr_n_active,
                    "chi_t_over_2pi": cpsqr_chi_t_over_2pi,
                },
                "n_tr": n_tr,
                "n_cav_padding": n_cav_padding,
                "dt_s": dt_s,
                "sigma_fraction": sigma_fraction,
                "inserted_x_pi_duration_s": inserted_x_pi_duration_s,
                "reduced_threshold": reduced_threshold,
                "joint_threshold": joint_threshold,
                "duration_grid": duration_grid,
            }

            print(json.dumps(active_config, indent=2))
            """
        ),
        md(
            """
            ## Derived Objects

            This cell constructs the two representative models, target specifications, run configurations, and manifold transition-frequency summaries directly from the tunable parameters above. Re-running this cell after changing the previous one is sufficient to update the downstream notebook logic.
            """
        ),
        code(
            """
            def include_chi_prime(model_variant: str) -> bool:
                return str(model_variant) == "chi_plus_chiprime"

            strict_model = build_model(
                include_chi_prime=include_chi_prime(strict_model_variant),
                n_active=strict_n_active,
                n_tr=n_tr,
                n_cav_padding=n_cav_padding,
            )
            strict_frame = build_frame(strict_model)
            strict_levels = logical_levels(strict_n_active)
            strict_duration_s = duration_from_chi_t(strict_chi_t_over_2pi)
            strict_run_config = make_run_config(strict_model, n_active=strict_n_active, duration_s=strict_duration_s, dt_s=dt_s, sigma_fraction=sigma_fraction)
            strict_spec = target_spec(strict_target_family, strict_n_active)
            strict_freqs_hz = manifold_transition_frequencies_hz(strict_model, strict_levels, strict_frame)

            cpsqr_model = build_model(
                include_chi_prime=include_chi_prime(cpsqr_model_variant),
                n_active=cpsqr_n_active,
                n_tr=n_tr,
                n_cav_padding=n_cav_padding,
            )
            cpsqr_frame = build_frame(cpsqr_model)
            cpsqr_levels = logical_levels(cpsqr_n_active)
            cpsqr_duration_s = duration_from_chi_t(cpsqr_chi_t_over_2pi)
            cpsqr_run_config = make_run_config(cpsqr_model, n_active=cpsqr_n_active, duration_s=cpsqr_duration_s, dt_s=dt_s, sigma_fraction=sigma_fraction)
            cpsqr_spec = target_spec(cpsqr_target_family, cpsqr_n_active)
            cpsqr_freqs_hz = manifold_transition_frequencies_hz(cpsqr_model, cpsqr_levels, cpsqr_frame)

            display(Markdown("### Strict representative manifold frequencies (Hz)"))
            display(pd.DataFrame({"level": strict_levels, "frequency_hz": strict_freqs_hz}))
            display(Markdown("### CPSQR representative manifold frequencies (Hz)"))
            display(pd.DataFrame({"level": cpsqr_levels, "frequency_hz": cpsqr_freqs_hz}))
            """
        ),
        md(
            """
            ## Step 1: Load Saved Study Summary

            This default cell loads the saved study summary and the main comparison tables. It is the fastest way to verify the headline claims about selected families, threshold behavior, and the patched negative control. The tunable parameters above do not change the saved summary itself, but they control which rows we inspect downstream.
            """
        ),
        code(
            """
            # --- Load saved results (default) ---
            summary = load_json(data_dir / "study_summary.json")
            comparison_df = load_rows(data_dir / "comparison_results.json")
            duration_df = load_rows(data_dir / "duration_results.json")
            duration_echo_df = load_rows(data_dir / "duration_echoed_results.json")
            negative_control = load_json(data_dir / "negative_controls.json")

            display(Markdown("### Executive Summary"))
            for line in summary["executive_summary"]:
                display(Markdown(f"- {line}"))

            display(Markdown("### Selected Families"))
            display(pd.DataFrame({"family_name": summary["selected_families"]}))
            """
        ),
        md(
            """
            ## Step 2: Inspect the Saved Comparison Grid

            This cell filters the saved comparison grid using the currently selected strict and CPSQR settings. It is the main bridge from the study-wide summary to the specific artifact files inspected below. If you change the family, target, duration, or model in the tunable-parameter cell, re-run this cell to update the filtered leaderboards.
            """
        ),
        code(
            """
            strict_filter = (
                (comparison_df["model_variant"] == strict_model_variant)
                & (comparison_df["target_family"] == strict_target_family)
                & (comparison_df["n_active"] == strict_n_active)
                & (comparison_df["chi_t_over_2pi"] == strict_chi_t_over_2pi)
            )
            cpsqr_filter = (
                (comparison_df["model_variant"] == cpsqr_model_variant)
                & (comparison_df["target_family"] == cpsqr_target_family)
                & (comparison_df["n_active"] == cpsqr_n_active)
                & (comparison_df["chi_t_over_2pi"] == cpsqr_chi_t_over_2pi)
            )

            strict_cols = [
                "case_id",
                "family_name",
                "strict_reduced_quartet_mean",
                "strict_full_quartet_mean",
                "strict_joint_process_fidelity",
                "cpsqr_joint_process_fidelity",
            ]

            display(Markdown("### Strict-case leaderboard"))
            display(comparison_df.loc[strict_filter, strict_cols].sort_values("strict_joint_process_fidelity", ascending=False).reset_index(drop=True))

            display(Markdown("### CPSQR-case leaderboard"))
            display(comparison_df.loc[cpsqr_filter, strict_cols].sort_values("cpsqr_joint_process_fidelity", ascending=False).reset_index(drop=True))
            """
        ),
        md(
            """
            ## Step 3: Load the Representative Strict-SQR Artifact

            This cell loads the saved case artifact for the strongest direct strict-SQR solution. It reports the saved summary row, the optimized tone parameters, and the target angles. The current tunable parameters affect which comparison row is selected here.
            """
        ),
        code(
            """
            strict_row = (
                comparison_df.loc[strict_filter & (comparison_df["family_name"] == strict_family_name)]
                .sort_values("strict_joint_process_fidelity", ascending=False)
                .iloc[0]
            )
            strict_artifact_path = artifacts_dir / "cases" / f"{strict_row.case_id}_{strict_row.family_name}.json"
            strict_artifact = load_json(strict_artifact_path)

            display(Markdown(f"### Loaded strict artifact: `{strict_artifact_path.name}`"))
            display(pd.DataFrame([strict_artifact["summary_row"]]).T.rename(columns={0: "value"}))
            display(Markdown("### Direct tone rows"))
            display(pd.DataFrame(strict_artifact["metadata"]["tone_rows"]))
            """
        ),
        md(
            """
            ### Re-run with Current Parameters

            The next cell is intentionally commented out. Uncomment it only if you want to re-run the selected direct case live with the current parameter choices from the tunable-parameter cell. This step can take time because it launches the single-case optimization path.
            """
        ),
        code(
            """
            # --- Re-run with current parameters ---
            # import subprocess
            #
            # cmd = [
            #     sys.executable,
            #     str(scripts_dir / "run_study.py"),
            #     "--single-family",
            #     strict_family_name,
            #     "--single-target",
            #     strict_target_family,
            #     "--single-model",
            #     strict_model_variant,
            #     "--single-na",
            #     str(strict_n_active),
            #     "--single-chiT",
            #     str(strict_chi_t_over_2pi),
            # ]
            # print(" ".join(cmd))
            # subprocess.run(cmd, check=True, cwd=study_dir)
            """
        ),
        md(
            """
            ## Step 4: Load the Representative CPSQR / Echoed Artifact

            This cell loads the saved echoed solution that most cleanly demonstrates the difference between strict SQR and CPSQR. It reports the saved summary row, the segment-level correction parameters, and the fitted per-manifold CPSQR phase corrections.
            """
        ),
        code(
            """
            cpsqr_row = (
                comparison_df.loc[cpsqr_filter & (comparison_df["family_name"] == cpsqr_family_name)]
                .sort_values("cpsqr_joint_process_fidelity", ascending=False)
                .iloc[0]
            )
            cpsqr_artifact_path = artifacts_dir / "cases" / f"{cpsqr_row.case_id}_{cpsqr_row.family_name}.json"
            cpsqr_artifact = load_json(cpsqr_artifact_path)

            display(Markdown(f"### Loaded CPSQR artifact: `{cpsqr_artifact_path.name}`"))
            display(pd.DataFrame([cpsqr_artifact["summary_row"]]).T.rename(columns={0: "value"}))
            display(Markdown("### Echo segment 1 tone rows"))
            display(pd.DataFrame(cpsqr_artifact["metadata"]["tone_rows_segment_1"]))
            display(Markdown("### Echo segment 2 tone rows"))
            display(pd.DataFrame(cpsqr_artifact["metadata"]["tone_rows_segment_2"]))
            display(Markdown("### Fitted CPSQR block phases"))
            display(pd.DataFrame(cpsqr_artifact["cpsqr_fit_rows"]))
            """
        ),
        md(
            """
            ### Re-run with Current Parameters

            The next cell is also commented out. Uncomment it only if you want to re-run the selected composite case live. This is the most expensive path in the notebook because the echoed family optimizes two multitone halves.
            """
        ),
        code(
            """
            # --- Re-run with current parameters ---
            # import subprocess
            #
            # cmd = [
            #     sys.executable,
            #     str(scripts_dir / "run_study.py"),
            #     "--single-family",
            #     cpsqr_family_name,
            #     "--single-target",
            #     cpsqr_target_family,
            #     "--single-model",
            #     cpsqr_model_variant,
            #     "--single-na",
            #     str(cpsqr_n_active),
            #     "--single-chiT",
            #     str(cpsqr_chi_t_over_2pi),
            # ]
            # print(" ".join(cmd))
            # subprocess.run(cmd, check=True, cwd=study_dir)
            """
        ),
        md(
            """
            ## Step 5: Validation Checks

            This section reproduces the main validation arguments from the report:

            - the patched negative control,
            - the refined duration thresholds,
            - the fact that reduced or full-state quartet success does not certify the strict joint operator.

            The tunable model and family choices affect which threshold rows are shown in the filtered tables below.
            """
        ),
        code(
            """
            display(Markdown("### Patched negative control"))
            display(pd.DataFrame([negative_control]).T.rename(columns={0: "value"}))

            duration_all = pd.concat([duration_df, duration_echo_df], ignore_index=True)
            threshold_cols = ["family_name", "model_variant", "n_active", "metric_name", "threshold", "minimum_chi_t_over_2pi", "best_value"]

            focus_thresholds = duration_all[
                ((duration_all["model_variant"] == strict_model_variant) & (duration_all["n_active"] == strict_n_active))
                | ((duration_all["model_variant"] == cpsqr_model_variant) & (duration_all["n_active"] == cpsqr_n_active))
            ]

            display(Markdown("### Filtered duration thresholds"))
            display(focus_thresholds[threshold_cols].sort_values(["family_name", "metric_name", "model_variant", "n_active"]).reset_index(drop=True))

            quartet_gap = pd.DataFrame(
                [
                    {
                        "case": strict_artifact_path.name,
                        "strict_reduced_quartet_mean": strict_artifact["summary_row"]["strict_reduced_quartet_mean"],
                        "strict_full_quartet_mean": strict_artifact["summary_row"]["strict_full_quartet_mean"],
                        "strict_joint_process_fidelity": strict_artifact["summary_row"]["strict_joint_process_fidelity"],
                    },
                    {
                        "case": cpsqr_artifact_path.name,
                        "strict_reduced_quartet_mean": cpsqr_artifact["summary_row"]["strict_reduced_quartet_mean"],
                        "strict_full_quartet_mean": cpsqr_artifact["summary_row"]["strict_full_quartet_mean"],
                        "strict_joint_process_fidelity": cpsqr_artifact["summary_row"]["strict_joint_process_fidelity"],
                    },
                ]
            )
            display(Markdown("### Quartet-versus-joint comparison for the two representative artifacts"))
            display(quartet_gap)
            """
        ),
        md(
            """
            ## Step 6: Key Figures

            This section re-displays the main report figures from the saved image files. It is the fastest way to visually confirm the main claims without re-running any optimization. The selected tunable parameters do not change which figure files are loaded here; the figures correspond to the saved report artifacts.
            """
        ),
        code(
            """
            for figure_name in [
                "reduced_and_joint_fidelity_vs_duration.png",
                "direct_vs_echoed_comparison.png",
                "sqr_vs_cpsqr_joint_comparison.png",
                "residual_z_and_transverse_vs_duration.png",
                "strict_best_waveform_and_spectrum.png",
                "cpsqr_best_waveform_and_spectrum.png",
            ]:
                display(Markdown(f"### {figure_name}"))
                display(Image(filename=str(figures_dir / figure_name)))
            """
        ),
        md(
            """
            ## Step 7: Summary and Parameter Map

            The final cell builds a compact parameter map so users and future agents can see how each tunable parameter affects what is loaded or re-run. This is also the place to sanity-check that the notebook is still pointed at the intended representative artifacts.
            """
        ),
        code(
            """
            parameter_map = pd.DataFrame(
                [
                    {"parameter": "strict_model_variant", "default_value": strict_model_variant, "effect": "Chooses the model for the representative strict-SQR artifact and re-run command."},
                    {"parameter": "strict_target_family", "default_value": strict_target_family, "effect": "Chooses the target profile for the representative strict-SQR artifact."},
                    {"parameter": "strict_family_name", "default_value": strict_family_name, "effect": "Chooses which direct family is loaded or re-run for strict SQR."},
                    {"parameter": "strict_n_active", "default_value": strict_n_active, "effect": "Sets the addressed manifold count for the strict representative case."},
                    {"parameter": "strict_chi_t_over_2pi", "default_value": strict_chi_t_over_2pi, "effect": "Sets the dimensionless duration for the strict representative case."},
                    {"parameter": "cpsqr_model_variant", "default_value": cpsqr_model_variant, "effect": "Chooses the model for the representative echoed CPSQR artifact and re-run command."},
                    {"parameter": "cpsqr_target_family", "default_value": cpsqr_target_family, "effect": "Chooses the target profile for the representative CPSQR artifact."},
                    {"parameter": "cpsqr_family_name", "default_value": cpsqr_family_name, "effect": "Chooses which composite family is loaded or re-run for CPSQR."},
                    {"parameter": "cpsqr_n_active", "default_value": cpsqr_n_active, "effect": "Sets the addressed manifold count for the representative CPSQR case."},
                    {"parameter": "cpsqr_chi_t_over_2pi", "default_value": cpsqr_chi_t_over_2pi, "effect": "Sets the dimensionless duration for the representative CPSQR case."},
                    {"parameter": "n_tr", "default_value": n_tr, "effect": "Sets the transmon truncation used when derived objects are rebuilt in this notebook."},
                    {"parameter": "n_cav_padding", "default_value": n_cav_padding, "effect": "Controls the cavity truncation via N_cav = N_active + padding."},
                    {"parameter": "dt_s", "default_value": dt_s, "effect": "Controls the simulation time step in derived run configurations and re-run commands."},
                    {"parameter": "sigma_fraction", "default_value": sigma_fraction, "effect": "Controls the Gaussian width used in the derived run configurations."},
                    {"parameter": "inserted_x_pi_duration_s", "default_value": inserted_x_pi_duration_s, "effect": "Documents the echo-pulse duration assumed by the saved echoed construction."},
                    {"parameter": "reduced_threshold", "default_value": reduced_threshold, "effect": "Sets the reduced-fidelity benchmark used when reading threshold tables."},
                    {"parameter": "joint_threshold", "default_value": joint_threshold, "effect": "Sets the joint-fidelity benchmark used when reading threshold tables."},
                    {"parameter": "duration_grid", "default_value": duration_grid, "effect": "Documents the refined dimensionless duration sweep used in the saved study."},
                ]
            )

            display(parameter_map)
            display(
                Markdown(
                    \"\"\"\n### Final caveat\nThe notebook is intentionally artifact-first. Strong reduced or even per-input full-state fidelities should not be interpreted as full strict-SQR success unless the corresponding joint operator metrics also agree.\n\"\"\"\n                )
            )
            """
        ),
    ]

    notebook = {
        "cells": cells,
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3",
            },
            "language_info": {
                "name": "python",
                "version": "3.12",
            },
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }

    NB_PATH.write_text(json.dumps(notebook, indent=2), encoding="utf-8")
    print(f"Wrote notebook to {NB_PATH}")


if __name__ == "__main__":
    main()
