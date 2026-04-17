from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import subprocess
from time import perf_counter
from typing import Any

import numpy as np
import pandas as pd

from diagnostics import build_tasks, run_task
from extra_analyses import compute_scaling_summary, compute_spectral_crowding, compute_xpi_characterization, xpi_dataframe
from load_prior_studies import DEFINITIVE_STUDY, STUDY_PATHS, combined_dataframe, load_normalized_results, save_snapshots
from parallel_utils import add_parallel_cli, default_worker_count, run_tasks
from plotting import generate_all_figures
from validate_results import write_validation_outputs


DATA_DIR = DEFINITIVE_STUDY / "data"
NEW_RESULTS_DIR = DATA_DIR / "new_results"
VALIDATION_DIR = DATA_DIR / "validation"
FIGURE_ROOT = DEFINITIVE_STUDY / "figures"
REPORT_DIR = DEFINITIVE_STUDY / "report"
SCRIPT_DIR = DEFINITIVE_STUDY / "scripts"


def _save_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _latex_escape(value: Any) -> str:
    text = str(value)
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
    }
    for src, dst in replacements.items():
        text = text.replace(src, dst)
    return text


def _deduplicate_rows(df: pd.DataFrame) -> pd.DataFrame:
    scored = df.copy()
    scored["_score"] = scored["strict_process_fidelity"].fillna(scored["avg_gate_fidelity"]).fillna(-1.0)
    deduped = (
        scored.sort_values(["study", "case_id", "construction", "_score"], ascending=[True, True, True, False])
        .drop_duplicates(subset=["study", "case_id", "construction"], keep="first")
        .drop(columns="_score")
        .reset_index(drop=True)
    )
    return deduped


def _fundamental_limits(df: pd.DataFrame) -> dict[str, Any]:
    ideal = df[(df["target_kind"] == "ideal_sqr") & df["strict_process_fidelity"].notna()].copy()
    ratios = ideal.copy()
    ratios["parameter_ratio"] = ratios["parameter_count"] / ratios["target_parameter_count"]
    duration_native = ideal[ideal["study"] == "native_rich_multitone_sqr_cpsqr_feasibility"].copy()
    if not duration_native.empty:
        duration_native["chi_t_over_2pi"] = duration_native["case_id"].str.extract(r"chiT([0-9p]+)")[0].str.replace("p", ".", regex=False).astype(float)
    extrapolations = []
    for construction, group in duration_native.groupby("construction"):
        grouped = group.groupby("chi_t_over_2pi")["strict_process_fidelity"].max().reset_index().sort_values("chi_t_over_2pi")
        valid = grouped[grouped["strict_process_fidelity"] < 0.999]
        if len(valid) < 2:
            continue
        residual = np.maximum(1.0e-6, 1.0 - valid["strict_process_fidelity"].to_numpy())
        coeffs = np.polyfit(valid["chi_t_over_2pi"].to_numpy(), np.log(residual), deg=1)
        slope = float(coeffs[0])
        intercept = float(coeffs[1])
        target_residual = 1.0 - 0.99
        chi_t_needed = None
        if slope < 0.0:
            chi_t_needed = float((math.log(target_residual) - intercept) / slope)
        extrapolations.append(
            {
                "construction": construction,
                "fit_kind": "log-linear residual",
                "slope": slope,
                "intercept": intercept,
                "estimated_chi_t_for_0p99": chi_t_needed,
            }
        )
    return {
        "target_parameter_rule": "An ideal x-axis SQR on N_active addressed manifolds needs 2 * N_active real parameters: one rotation angle and one phase-lock constraint per manifold.",
        "parameter_count_rows": ratios[
            [
                "study",
                "case_id",
                "construction",
                "n_active",
                "parameter_count",
                "target_parameter_count",
                "parameter_ratio",
                "strict_process_fidelity",
            ]
        ].to_dict(orient="records"),
        "duration_extrapolations": extrapolations,
        "qsl_note": "This definitive pass does not derive a formal quantum speed limit from a control-norm bound. The limit analysis instead uses ansatz parameter counts plus empirical duration scaling from the saved native-rich runs.",
    }


def _construction_summary(df: pd.DataFrame, error_df: pd.DataFrame) -> dict[str, Any]:
    summary_rows = (
        df[df["strict_process_fidelity"].notna()]
        .groupby(["study", "construction", "construction_display", "construction_family"], as_index=False)[
            ["strict_process_fidelity", "avg_gate_fidelity", "cpsqr_process_fidelity", "parameter_count"]
        ]
        .max()
        .sort_values("strict_process_fidelity", ascending=False)
    )
    if error_df.empty:
        error_summary = []
    else:
        temp = error_df.copy()
        temp["abs_eps_x"] = np.abs(temp["eps_x"])
        temp["abs_eps_y"] = np.abs(temp["eps_y"])
        temp["abs_eps_z"] = np.abs(temp["eps_z"])
        temp["abs_eps_norm"] = np.abs(temp["eps_norm"])
        error_summary = (
            temp.groupby(["construction"], as_index=False)[["abs_eps_x", "abs_eps_y", "abs_eps_z", "abs_eps_norm"]]
            .mean()
            .to_dict(orient="records")
        )
    return {
        "top_strict_rows": summary_rows.head(24).to_dict(orient="records"),
        "error_summary": error_summary,
    }


def _master_results_table(df: pd.DataFrame, error_df: pd.DataFrame) -> pd.DataFrame:
    merged = df.copy()
    if not error_df.empty:
        temp = error_df.copy()
        temp["abs_eps_z"] = np.abs(temp["eps_z"])
        temp["abs_eps_perp"] = np.sqrt(np.asarray(temp["eps_x"], dtype=float) ** 2 + np.asarray(temp["eps_y"], dtype=float) ** 2)
        error_summary = (
            temp.groupby(["case_id", "construction"], as_index=False)
            .agg(mean_abs_eps_z=("abs_eps_z", "mean"), mean_abs_eps_perp=("abs_eps_perp", "mean"))
        )
        merged = merged.merge(error_summary, on=["case_id", "construction"], how="left")
    merged["f_unitary"] = merged["strict_process_fidelity"].fillna(merged["avg_gate_fidelity"])
    total_duration_s = merged["total_gate_duration_ns"].fillna(merged["duration_ns"]).fillna(0.0) * 1.0e-9
    merged["f_decoherence"] = np.exp(-total_duration_s / 50.0e-6 - total_duration_s / 30.0e-6)
    merged["f_practical"] = merged["f_unitary"] * merged["f_decoherence"]
    merged["mean_abs_eps_z"] = merged.get("mean_abs_eps_z", pd.Series(np.nan, index=merged.index))
    merged["mean_abs_eps_perp"] = merged.get("mean_abs_eps_perp", pd.Series(np.nan, index=merged.index))
    merged = merged.sort_values("f_practical", ascending=False).reset_index(drop=True)
    merged["rank"] = np.arange(1, len(merged) + 1)
    return merged[
        [
            "rank",
            "construction_display",
            "ansatz",
            "envelope",
            "duration_ns",
            "n_active",
            "f_unitary",
            "f_decoherence",
            "f_practical",
            "mean_abs_eps_z",
            "mean_abs_eps_perp",
            "parameter_count",
            "study",
            "case_id",
        ]
    ].rename(columns={"construction_display": "construction"})


def _objective_summaries(
    *,
    combined: pd.DataFrame,
    new_results: pd.DataFrame,
    xpi_payload: dict[str, Any],
    scaling_payload: dict[str, Any],
    figure_stems: list[str],
    wall_clock_seconds: float,
) -> list[dict[str, Any]]:
    strict_best = float(new_results["strict_process_fidelity"].max())
    cpsqr_best = float(new_results["cpsqr_process_fidelity"].max())
    baseline_direct = combined[
        (combined["study"] == "ideal_sqr_direct_vs_echoed_multitone") & (combined["construction"] == "direct_multitone")
    ]
    baseline_avg = float(baseline_direct["avg_gate_fidelity"].max())
    xpi_best_rows = xpi_payload.get("best_configs", [])
    return [
        {
            "objective": 1,
            "tasks_completed": [
                "Loaded and normalized all prior SQR studies plus the native-rich extension.",
                "Computed reusable per-manifold Pauli error decompositions for the ideal-SQR artifacts with saved propagators.",
                "Generated unified prior-work and full-grid comparison figures.",
            ],
            "headline_result": "The original baseline negative result persists, but the broader repository corpus already contains later direct-waveform strict-SQR successes above 0.99.",
            "best_fidelity": strict_best,
            "comparison_to_baseline": f"Baseline direct average gate fidelity remained {baseline_avg:.4f}, while the native-rich strict joint fidelity reached {strict_best:.5f}.",
            "figures_generated": [stem for stem in figure_stems if stem in {"unified_prior_work_comparison", "per_manifold_error_budget", "full_case_construction_heatmap"}],
            "wall_clock_seconds": wall_clock_seconds,
        },
        {
            "objective": 2,
            "tasks_completed": [
                "Characterized the standalone finite Gaussian refocusing pulse across addressed manifolds.",
                "Scanned a lightweight compromise pulse over duration, amplitude scale, sigma fraction, and carrier manifold.",
            ],
            "headline_result": "The inserted finite refocusing pulse is explicitly manifold dependent; a simple compromise pulse improves the worst-manifold fidelity but does not by itself certify echoed strict-SQR success.",
            "best_fidelity": None if not xpi_best_rows else float(max(row["worst_process_fidelity"] for row in xpi_best_rows)),
            "comparison_to_baseline": "This objective remains a standalone pulse audit only; a full echoed-grid rerun with the compromise pulse is still future work.",
            "figures_generated": [stem for stem in figure_stems if stem == "refocusing_pulse_manifold_dependence"],
            "wall_clock_seconds": wall_clock_seconds,
        },
        {
            "objective": 3,
            "tasks_completed": [
                "Aggregated richer direct families from the native-rich study and compared them across duration and model variant.",
                "Extracted a parameter-count-versus-fidelity tradeoff and a fixed-duration ansatz ranking.",
            ],
            "headline_result": "The reduced-unitary direct family is the strongest strict-SQR ansatz in the current repository corpus.",
            "best_fidelity": strict_best,
            "comparison_to_baseline": f"The strongest direct family reaches {strict_best:.5f} strict joint fidelity, far above the original 0.7245 average-gate baseline.",
            "figures_generated": [stem for stem in figure_stems if stem in {"fidelity_vs_duration_scaling", "ansatz_comparison_fixed_duration", "parameter_count_vs_fidelity"}],
            "wall_clock_seconds": wall_clock_seconds,
        },
        {
            "objective": 4,
            "tasks_completed": [
                "Reviewed the repository corpus for hybrid direct-plus-echo constructions.",
            ],
            "headline_result": "No dedicated hybrid SQR family was present in the saved native-rich or baseline studies, so the definitive study cannot make a data-backed hybrid claim yet.",
            "best_fidelity": None,
            "comparison_to_baseline": "Hybrid direct-plus-echo optimization remains an unexecuted follow-up rather than a validated result.",
            "figures_generated": [],
            "wall_clock_seconds": wall_clock_seconds,
        },
        {
            "objective": 5,
            "tasks_completed": [
                "Constructed a practical ranking table with the decoherence-only fallback factor allowed by the prompt.",
            ],
            "headline_result": "The practical ranking remains provisional because finite-difference drift sweeps and open-system replays were not rerun in this definitive pass.",
            "best_fidelity": strict_best,
            "comparison_to_baseline": "The report now ranks constructions by a transparent fallback figure of merit instead of by unitary fidelity alone.",
            "figures_generated": [],
            "wall_clock_seconds": wall_clock_seconds,
        },
        {
            "objective": 6,
            "tasks_completed": [
                "Computed an explicit spectral-crowding analysis from the patched Hamiltonian parameters.",
                "Extracted best strict fidelities versus addressed-manifold count and fit simple extrapolations toward N_active = 8.",
            ],
            "headline_result": "Strict ideal-SQR performance degrades systematically with addressed-manifold count, especially in the harder model variant.",
            "best_fidelity": strict_best,
            "comparison_to_baseline": "The original baseline stopped at small addressed windows; the definitive pass now exposes the observed scaling trend and crowding mechanism.",
            "figures_generated": [stem for stem in figure_stems if stem in {"n_active_scaling", "spectral_crowding_diagram"}],
            "wall_clock_seconds": wall_clock_seconds,
        },
        {
            "objective": 7,
            "tasks_completed": [
                "Computed ansatz parameter-count ratios for all loaded ideal-SQR rows.",
                "Derived empirical duration extrapolations toward the 0.99 threshold from the saved native-rich runs.",
            ],
            "headline_result": "The current failure cases are not explained by a simple parameter-count shortage alone; the ansatz-rich direct family still loses ground as the manifold ladder crowds.",
            "best_fidelity": strict_best,
            "comparison_to_baseline": f"The later native-rich study already reaches {cpsqr_best:.9f} CPSQR and {strict_best:.5f} strict SQR, so the present limit analysis is about generality, not existence.",
            "figures_generated": [stem for stem in figure_stems if stem == "parameter_count_vs_fidelity"],
            "wall_clock_seconds": wall_clock_seconds,
        },
    ]


def _write_notebook(summary: dict[str, Any]) -> None:
    notebook = {
        "cells": [
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "# The Definitive Ideal SQR Gate Study\n",
                    "\n",
                    "This notebook loads the definitive-study outputs generated from the combined SQR artifact corpus and reproduces the headline tables, figures, and objective summaries.\n",
                ],
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "from pathlib import Path\n",
                    "import json\n",
                    "import pandas as pd\n",
                    "ROOT = Path.cwd().parents[1]\n",
                    "STUDY = ROOT / 'studies' / 'the_definitive_ideal_sqr_gate_study'\n",
                    "NEW = STUDY / 'data' / 'new_results'\n",
                ],
            },
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": ["## User-Tunable Parameters\n", "Adjust these selectors before re-running the cells below.\n"],
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "selected_study = 'native_rich_multitone_sqr_cpsqr_feasibility'\n",
                    "selected_case_substring = 'chi_plus_chiprime_smooth_x'\n",
                    "construction_filter = 'direct'\n",
                    "print({'selected_study': selected_study, 'selected_case_substring': selected_case_substring, 'construction_filter': construction_filter})\n",
                ],
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "full_grid = pd.read_csv(NEW / 'full_grid_results.csv')\n",
                    "best_per_case = pd.read_csv(NEW / 'best_per_case.csv')\n",
                    "master = pd.read_csv(NEW / 'master_results_table.csv')\n",
                    "objectives = json.loads((NEW / 'objective_summaries.json').read_text())\n",
                    "full_grid.head()\n",
                ],
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": ["master.head(15)\n"],
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "subset = full_grid[(full_grid['study'] == selected_study) & (full_grid['case_id'].str.contains(selected_case_substring)) & (full_grid['construction_family'] == construction_filter)]\n",
                    "subset[['case_id','construction_display','strict_process_fidelity','cpsqr_process_fidelity','mean_residual_z_error_rad','mean_transverse_error_rad']].sort_values('strict_process_fidelity', ascending=False).head(20)\n",
                ],
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": ["xpi = json.loads((NEW / 'xpi_characterization.json').read_text())\n", "xpi['best_configs']\n"],
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [f"summary = {json.dumps(summary, indent=2)}\n", "summary['headline_findings']\n"],
            },
        ],
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python", "version": "3.12"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    (SCRIPT_DIR / "reproducibility_notebook.ipynb").write_text(json.dumps(notebook, indent=2), encoding="utf-8")


def _write_report(
    *,
    combined: pd.DataFrame,
    new_results: pd.DataFrame,
    master_table: pd.DataFrame,
    xpi_payload: dict[str, Any],
    scaling_payload: dict[str, Any],
) -> None:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    strict_best_row = combined.sort_values("strict_process_fidelity", ascending=False).iloc[0]
    cpsqr_best_row = combined.sort_values("cpsqr_process_fidelity", ascending=False, na_position="last").iloc[0]
    baseline_direct = combined[
        (combined["study"] == "ideal_sqr_direct_vs_echoed_multitone") & (combined["construction"] == "direct_multitone")
    ]
    baseline_best_avg = float(baseline_direct["avg_gate_fidelity"].max())
    baseline_best_strict = float(baseline_direct["strict_process_fidelity"].max())
    structured_native = new_results[new_results["target_family"].isin({"smooth_x", "staggered_x"})].copy()
    structured_native_best_by_case = structured_native.sort_values("strict_process_fidelity", ascending=False).drop_duplicates(subset=["case_id"])
    structured_best_mean = float(structured_native_best_by_case["strict_process_fidelity"].mean())
    hard_direct = new_results[
        (new_results["construction"] == "reduced_unitary_direct")
        & (new_results["model_variant"] == "chi_plus_chiprime")
        & (new_results["n_active"] == 3)
    ]
    hard_direct_best = float(hard_direct["strict_process_fidelity"].max()) if not hard_direct.empty else float("nan")
    prior_rows = (
        combined.groupby(["study_short", "study_role"], as_index=False)[["strict_process_fidelity", "cpsqr_process_fidelity"]]
        .max()
        .sort_values("strict_process_fidelity", ascending=False)
    )
    prior_table_rows = "\n".join(
        f"{_latex_escape(row.study_short)} & {_latex_escape(row.study_role)} & {row.strict_process_fidelity:.4f} & "
        f"{('--' if pd.isna(row.cpsqr_process_fidelity) else f'{row.cpsqr_process_fidelity:.4f}')} \\\\"
        for row in prior_rows.itertuples(index=False)
    )
    ansatz_short = {
        "reduced unitary direct": "reduced direct",
        "native direct strict": "native direct",
        "symmetric two segment": "two-segment",
        "complex envelope": "rich envelope",
        "basis expanded": "basis-expanded",
    }
    ranking_rows = "\n".join(
        f"{int(row.rank)} & {_latex_escape(row.construction)} & {_latex_escape(ansatz_short.get(str(row.ansatz), str(row.ansatz)))} & {int(row.n_active)} & "
        f"{row.f_unitary:.4f} & {row.f_practical:.4f} \\\\"
        for row in master_table.head(10).itertuples(index=False)
    )
    xpi_best_rows = xpi_payload.get("best_configs", [])
    xpi_text = "A lightweight compromise pulse was selected by maximizing the worst addressed-manifold channel fidelity to an ideal $X_\\pi$ rotation across $n=0,\\ldots,4$."
    if xpi_best_rows:
        formatted = []
        for row in xpi_best_rows:
            label = "simpler ladder" if row["model_variant"] == "chi_only" else "curved ladder"
            formatted.append(
                f"In the {label}, the best compromise pulse uses $T_\\pi={row['duration_ns']:.0f}$ ns and reaches worst-manifold fidelity {row['worst_process_fidelity']:.4f}"
            )
        xpi_text = ". ".join(formatted) + "."
    scaling_fits = scaling_payload.get("fits", [])
    scaling_text = "The saved native-rich corpus already shows a monotone degradation with addressed-manifold count."
    if scaling_fits:
        first_fit = scaling_fits[0]
        scaling_text = (
            f"A simple log-linear residual fit for {str(first_fit['construction']).replace('_', ' ')} in "
            f"{str(first_fit['model_variant']).replace('_', ' ')} predicts a strict process near "
            f"{float(first_fit['predicted_strict_process_at_n8']):.3f} at $N_{{\\mathrm{{active}}}}=8$."
        )
    tex = rf"""
\documentclass[aps,pra,twocolumn,reprint,floatfix,amsmath,amssymb]{{revtex4-2}}
\usepackage{{graphicx}}
\usepackage{{booktabs}}
\usepackage{{siunitx}}
\usepackage{{mathtools}}
\usepackage[colorlinks=true,linkcolor=blue,citecolor=blue,urlcolor=blue]{{hyperref}}
\emergencystretch=2em
\begin{{document}}
\sloppy
\title{{The Definitive Ideal SQR Gate Study}}
\author{{Codex}}
\affiliation{{Autonomous cQED Study Workspace}}
\date{{\today}}
\begin{{abstract}}
This definitive study unifies the repository's selective-qubit-rotation line into one self-contained answer. It re-presents the three historical SQR baselines, normalizes their saved outputs into a single comparison frame, imports the later patched native-rich extension as the strongest new-results layer, and adds fresh cross-study diagnostics: per-manifold ideal-SQR error generators, a standalone finite-refocusing-pulse audit, explicit spectral-crowding estimates, addressed-manifold scaling fits, and a ranked practical table. The central conclusion is conditional rather than universal. Parameterized direct multitone waveforms do realize strict ideal selective qubit rotations above 0.99 on physically credible addressed windows, with the best loaded strict joint fidelity reaching {float(strict_best_row.strict_process_fidelity):.5f}. That success does not generalize across the broader grid: the mean best strict fidelity on the structured native-rich case set is only {structured_best_mean:.4f}, far below the user's 0.95 grid-average landmark, and the harder three-manifold model with higher-order dispersive curvature still saturates around {hard_direct_best:.4f}. Echoed constructions remain compelling mainly for the relaxed conditional-phase SQR target, whose best loaded joint fidelity reaches {float(cpsqr_best_row.cpsqr_process_fidelity):.9f}. The decisive obstruction is therefore not the existence of strict ideal SQR, but its lack of broad robustness across the tested ansatz families and addressed-manifold counts.
\end{{abstract}}
\maketitle

\section{{Introduction}}
The ideal selective qubit rotation problem asks whether a single qubit-drive waveform can realize
\begin{{align}}
U_{{\mathrm{{SQR}}}}^{{\mathrm{{ideal}}}}
=
\sum_{{n\in\mathcal{{A}}}} |n\rangle\langle n| \otimes R_x(\theta_n)
\end{{align}}
on an addressed cavity window $\mathcal{{A}}=\{{0,\ldots,N_{{\mathrm{{active}}}}-1\}}$. The strict target locks both the rotation angle and the rotation axis on every addressed manifold while also preserving the relative inter-manifold phase structure of the full joint qubit-cavity operator. The guiding question of this report is therefore sharper than the original baseline: can a parameterized multitone waveform realize $U_{{\mathrm{{SQR}}}}^{{\mathrm{{ideal}}}}$ with $F_{{\mathrm{{strict}}}}>0.99$, and if not broadly, where does the remaining error floor live?

\section{{System Model}}
All numerical source studies use the patched public \texttt{{cqed\_sim}} dispersive transmon-cavity model \cite{{cqedsim2026,blais2021cqed,koch2007transmon}}. The addressed manifold transition frequencies are
\begin{{align}}
\Delta_n = \chi n + \chi' n(n-1),
\end{{align}}
with a lighter model that retains only $\chi$ and a harder model that keeps the higher-order $\chi'$ term. The study corpus spans addressed windows with $N_{{\mathrm{{active}}}}=2,3,4$ and durations reported as the dimensionless product $|\chi|T/2\pi$. The two principal Hamiltonian variants are referred to below as the simpler dispersive ladder and the curved dispersive ladder.

\section{{Prior Results, Unified}}
Figure~\ref{{fig:prior_work}} and Table~\ref{{tab:prior_work}} place the full repository history in one frame. The arbitrary-target and residual-phase studies broadened the ansatz space and clarified the role of coherent phase structure, but they did not answer the strict ideal-SQR question directly. The direct-versus-echoed baseline then answered that strict question negatively on a $16$-case Gaussian grid, reaching only a best direct average gate fidelity of {baseline_best_avg:.4f} and a best strict process fidelity of {baseline_best_strict:.4f}. The later native-rich extension changed the answer qualitatively by finding strict joint fidelities above $0.99$ on restricted direct-waveform cases while also demonstrating near-unit conditional-phase SQR on echoed cases.
\begin{{figure}}[t]
  \centering
  \includegraphics[width=\columnwidth]{{../figures/pdf/unified_prior_work_comparison.pdf}}
  \caption{{Best strict process fidelities from the three historical baselines and the native-rich extension, replotted in one consistent visual frame.}}
  \label{{fig:prior_work}}
\end{{figure}}
\begin{{table}}[t]
\scriptsize
\caption{{Unified repository priors plus the current new-results layer.}}
\label{{tab:prior_work}}
\begin{{ruledtabular}}
\begin{{tabular}}{{lccc}}
Study & Role & Best strict & Best CPSQR \\
\midrule
{prior_table_rows}
\end{{tabular}}
\end{{ruledtabular}}
\end{{table}}

\section{{New Methods in This Definitive Pass}}
The new work in this folder is not a fresh multi-hour optimization campaign. It is a strict unification-and-diagnostics pass over the saved repository corpus. The definitive scripts now perform four new tasks that did not exist in one place before. They normalize every prior SQR study into one common schema. They recompute per-manifold ideal-SQR error generators from saved propagators. They audit the finite Gaussian refocusing pulse used by the echoed constructions. They also add explicit spectral-crowding and addressed-manifold scaling analyses. The refocusing-pulse audit is intentionally lightweight but informative. {xpi_text} The scaling analysis is likewise empirical rather than axiomatic. {scaling_text}

\section{{Per-Manifold Error Anatomy}}
Figure~\ref{{fig:error_budget}} shows the missing diagnostic that the original baseline lacked: the per-manifold error-generator decomposition
\begin{{align}}
U_n^{{\mathrm{{real}}}} &= \exp\!\left(-\frac{{i}}{{2}} E_n\right) U_n^{{\mathrm{{ideal}}}},
E_n &= \epsilon_x^{{(n)}} X \notag\\
&\quad + \epsilon_y^{{(n)}} Y + \epsilon_z^{{(n)}} Z .
\end{{align}}
The strong direct solutions suppress all three components simultaneously, but the hard cases are not pure residual-$Z$ failures. The remaining error floor is mixed: transverse under-rotation and axis tilt remain comparable to or larger than the residual conditional phase on the hardest structured windows. That observation is consistent with the broader direct-waveform story from the native-rich study: more parameters help, but they do not erase the crowding-induced control tension once the addressed ladder bends.
\begin{{figure}}[t]
  \centering
  \includegraphics[width=\columnwidth]{{../figures/pdf/per_manifold_error_budget.pdf}}
  \caption{{Representative per-manifold ideal-SQR error budgets for the strongest saved constructions.}}
  \label{{fig:error_budget}}
\end{{figure}}

\section{{Head-to-Head Construction Comparison}}
The strict comparison between direct and non-direct constructions is summarized in Fig.~\ref{{fig:direct_scatter}}. The direct waveform remains the clear winner for strict ideal SQR. Echoed and composite constructions help primarily when the target is relaxed to a conditional-phase SQR, not when the full strict joint operator is enforced. This is the central conceptual correction to the original echoed narrative: the echoed pulses are not worthless, but they are solving a weaker problem than strict ideal SQR.
\begin{{figure}}[t]
  \centering
  \includegraphics[width=\columnwidth]{{../figures/pdf/direct_vs_echoed_fidelity_scatter.pdf}}
  \caption{{Best non-direct strict fidelity versus best direct strict fidelity for every ideal-SQR case that has both.}}
  \label{{fig:direct_scatter}}
\end{{figure}}
\begin{{figure}}[t]
  \centering
  \includegraphics[width=\columnwidth]{{../figures/pdf/strict_vs_cpsqr_joint_comparison.pdf}}
  \caption{{Family-averaged strict versus conditional-phase joint fidelity across the native-rich corpus. Composite families shift upward because they repair conditional phase more effectively than they repair strict inter-manifold phase structure.}}
  \label{{fig:strict_vs_cpsqr}}
\end{{figure}}

\section{{The Push to Fidelity Above 0.99}}
The repository now contains a genuine positive answer to the existence question. The best loaded strict direct case reaches {float(strict_best_row.strict_process_fidelity):.5f}, well above the user's 0.99 threshold, and several other direct cases also cross that line on the easier addressed windows. Figure~\ref{{fig:duration_scaling}} shows how that success depends on duration: strict direct fidelities improve systematically with $|\chi|T/2\pi$, whereas echoed strict fidelities remain bounded well below the same threshold. Figure~\ref{{fig:ansatz_count}} shows the price of that improvement. Richer direct ansatz families move the fidelity frontier upward, but only up to the point allowed by the addressed-manifold crowding.
\begin{{figure}}[t]
  \centering
  \includegraphics[width=\columnwidth]{{../figures/pdf/fidelity_vs_duration_scaling.pdf}}
  \caption{{Strict fidelity versus dimensionless duration across the native-rich construction families.}}
  \label{{fig:duration_scaling}}
\end{{figure}}
\begin{{figure}}[t]
  \centering
  \includegraphics[width=\columnwidth]{{../figures/pdf/ansatz_comparison_fixed_duration.pdf}}
  \caption{{Longest-duration strict fidelity comparison across the saved native-rich ansatz families.}}
  \label{{fig:ansatz_long}}
\end{{figure}}
\begin{{figure}}[t]
  \centering
  \includegraphics[width=\columnwidth]{{../figures/pdf/parameter_count_vs_fidelity.pdf}}
  \caption{{Approximate ansatz parameter count versus best achieved strict fidelity.}}
  \label{{fig:ansatz_count}}
\end{{figure}}

\section{{Practical Viability, Scaling, and Crowding}}
The user's strongest grid-average landmark is not met in the current corpus: the mean best strict fidelity over the structured native-rich cases is only {structured_best_mean:.4f}. Figures~\ref{{fig:xpi}}--\ref{{fig:crowding}} explain why. The finite inserted $X_\pi$ pulse used by the echoed constructions is itself manifold dependent, so the clean first-order echo logic \cite{{hahn1950echo}} is physically obstructed before any higher-order optimization issue is reached. At the same time, the manifold ladder crowds as $N_{{\mathrm{{active}}}}$ grows, especially in the curved dispersive ladder, so direct ansatz richness buys only partial relief. The practical ranking table in Table~\ref{{tab:ranking}} therefore uses a transparent decoherence-only fallback factor and is explicitly provisional until full drift and open-system sweeps are replayed.
\begin{{figure}}[t]
  \centering
  \includegraphics[width=\columnwidth]{{../figures/pdf/refocusing_pulse_manifold_dependence.pdf}}
  \caption{{Standalone manifold-resolved characterization of the baseline finite Gaussian refocusing pulse and a lightweight compromise pulse.}}
  \label{{fig:xpi}}
\end{{figure}}
\begin{{figure}}[t]
  \centering
  \includegraphics[width=\columnwidth]{{../figures/pdf/n_active_scaling.pdf}}
  \caption{{Best strict fidelity versus addressed-manifold count across the native-rich corpus, with simple residual extrapolations toward $N_{{\mathrm{{active}}}}=8$.}}
  \label{{fig:na_scaling}}
\end{{figure}}
\begin{{figure}}[t]
  \centering
  \includegraphics[width=\columnwidth]{{../figures/pdf/spectral_crowding_diagram.pdf}}
  \caption{{Transition ladders and Gaussian bandwidth bands for a representative direct-waveform duration. The band overlap grows rapidly with addressed-manifold count in the curved dispersive ladder.}}
  \label{{fig:crowding}}
\end{{figure}}
\begin{{table*}}[t]
\scriptsize
\caption{{Top provisional practical rankings from the combined corpus.}}
\label{{tab:ranking}}
\begin{{ruledtabular}}
\begin{{tabular}}{{lccccc}}
Rank & Construction & Ansatz & $N_{{\mathrm{{active}}}}$ & $F_{{\mathrm{{unitary}}}}$ & $F_{{\mathrm{{practical}}}}$ \\
\midrule
{ranking_rows}
\end{{tabular}}
\end{{ruledtabular}}
\end{{table*}}

\section{{Discussion}}
The definitive answer is therefore tiered. On the positive side, strict ideal SQR is already demonstrated in the repository: direct multitone waveforms can exceed $0.99$ on physically credible cases, so the problem is not impossible within the dispersive model. On the negative side, the saved ansatz families do not achieve that success broadly. The grid-average landmark of $0.95$ is missed by a wide margin, the hard three-manifold curved-ladder cases saturate below the same threshold, and no echoed construction rescues strict joint SQR once finite refocusing pulses and full inter-manifold phase structure are enforced. The best interpretation is that the remaining error floor is partly ansatz-limited and partly crowding-limited. Richer direct parameterizations and longer durations help, but the improvement weakens as the addressed ladder curves. Echo helps strongly for conditional-phase SQR because it mainly removes the phase structure that the relaxed target tolerates, not because it produces a broad strict-SQR fix.

\section{{Conclusion and Roadmap}}
The central question can now be answered cleanly. Yes: parameterized direct multitone waveforms can realize strict ideal selective qubit rotations above $0.99$ on selected addressed windows, and the combined repository corpus already contains that Tier-1 landmark result. No: the tested ansatz families do not make that success broad across the primary structured grid, and the saved echoed constructions do not change that answer for the strict full joint operator. If a user needs a strict SQR tomorrow, the best available construction is the strongest direct reduced-unitary family on the easier addressed windows. If the application tolerates a conditional post-rotation phase, the echoed conditional-phase family is the more robust practical recommendation. The most valuable next step is not another qualitative echoed claim, but a quantitative follow-up that combines a truly robust refocusing pulse, an inter-manifold-phase-aware composite objective, and explicit open-system plus drift replays for the top direct and echoed candidates.

\section{{Validation}}
The definitive study inherits the source-study validation files and also checks four repository-level facts directly: the baseline negative result is preserved, the native-rich extension contains strict ideal-SQR cases above $0.99$, the same extension contains near-unit conditional-phase SQR cases, and the full corpus loads into one normalized table without silent omissions.

\section{{Limitations and Future Work}}
This definitive pass is still an analysis-and-unification study rather than a fresh from-scratch optimization campaign. The finite-difference drift sweeps, full Lindblad replays, hybrid direct-plus-echo family, and true unconstrained GRAPE upper bound remain future work. Those omissions matter, so the practical ranking should be read as a provisional engineering guide, not as the final last word on hardware robustness.

\appendix
\section{{Detailed Results and Data}}
Figure~\ref{{fig:heatmap}} shows the full ideal-SQR case-by-construction heatmap for the loaded corpus. Figure~\ref{{fig:validation_ladder}} summarizes the gap between reduced-manifold diagnostics, full-manifold quartet diagnostics, and the decisive joint-process metric.
\begin{{figure*}}[t]
  \centering
  \includegraphics[width=0.98\textwidth,height=0.80\textheight,keepaspectratio]{{../figures/pdf/full_case_construction_heatmap.pdf}}
  \caption{{Full ideal-SQR case-by-construction heatmap across the loaded studies.}}
  \label{{fig:heatmap}}
\end{{figure*}}
\begin{{figure}}[t]
  \centering
  \includegraphics[width=\columnwidth]{{../figures/pdf/one_state_vs_quartet_validation.pdf}}
  \caption{{Mean reduced-quartet, full-quartet, and joint-process fidelities across the native-rich families.}}
  \label{{fig:validation_ladder}}
\end{{figure}}

\section{{Reproducibility}}
The study folder contains machine-readable prior-study snapshots, normalized strict-result tables, per-manifold error decompositions, standalone refocusing-pulse diagnostics, spectral-crowding rows, figure JSON sidecars, a ranked master table, and a reproducibility notebook in \texttt{{scripts/reproducibility\_notebook.ipynb}}.

\bibliographystyle{{apsrev4-2}}
\bibliography{{references}}
\end{{document}}
"""
    (REPORT_DIR / "main.tex").write_text(tex.strip() + "\n", encoding="utf-8")
    (REPORT_DIR / "report.tex").write_text("\\input{main}\n", encoding="utf-8")
    bib = r"""
@article{blais2021cqed,
  author = {Blais, Alexandre and Grimsmo, Arne L. and Girvin, Steven M. and Wallraff, Andreas},
  title = {Circuit Quantum Electrodynamics},
  journal = {Reviews of Modern Physics},
  volume = {93},
  number = {2},
  pages = {025005},
  year = {2021},
  doi = {10.1103/RevModPhys.93.025005}
}

@article{koch2007transmon,
  author = {Koch, Jens and Yu, Terri M. and Gambetta, Jay and Houck, A. A. and Schuster, D. I. and Majer, J. and Blais, Alexandre and Devoret, M. H. and Girvin, S. M. and Schoelkopf, R. J.},
  title = {Charge-Insensitive Qubit Design Derived from the Cooper Pair Box},
  journal = {Physical Review A},
  volume = {76},
  number = {4},
  pages = {042319},
  year = {2007},
  doi = {10.1103/PhysRevA.76.042319}
}

@article{hahn1950echo,
  author = {Hahn, E. L.},
  title = {Spin Echoes},
  journal = {Physical Review},
  volume = {80},
  number = {4},
  pages = {580--594},
  year = {1950},
  doi = {10.1103/PhysRev.80.580}
}

@misc{cqedsim2026,
  author = {{Shankar Quantum Circuits Group}},
  title = {cqed\_sim: Circuit Quantum Electrodynamics Simulation Framework},
  year = {2026},
  howpublished = {\url{https://github.com/SoraUmika/qubox_cQEDsim}},
  note = {Accessed March 28, 2026}
}
"""
    (REPORT_DIR / "references.bib").write_text(bib.strip() + "\n", encoding="utf-8")


def _compile_report() -> bool:
    command_sets = [
        [
            ["latexmk", "-pdf", "-interaction=nonstopmode", "main.tex"],
        ],
        [
            ["pdflatex", "-interaction=nonstopmode", "main.tex"],
            ["bibtex", "main"],
            ["pdflatex", "-interaction=nonstopmode", "main.tex"],
            ["pdflatex", "-interaction=nonstopmode", "main.tex"],
        ],
    ]
    for commands in command_sets:
        ok = True
        for command in commands:
            try:
                completed = subprocess.run(command, cwd=REPORT_DIR, capture_output=True, text=True, check=False)
            except FileNotFoundError:
                ok = False
                break
            if completed.returncode != 0:
                ok = False
                break
        if ok:
            return True
    return False


def main() -> None:
    parser = argparse.ArgumentParser(description="Build the definitive ideal-SQR study from the saved repository corpus.")
    add_parallel_cli(parser)
    parser.add_argument("--quick", action="store_true", help="Limit the unified table to the ideal-SQR baseline and native-rich extension.")
    parser.add_argument("--full", action="store_true", help="Load the full configured SQR corpus. This is the default behavior.")
    args = parser.parse_args()

    started = perf_counter()
    print(f"Detected worker count: {default_worker_count()} (using {args.n_workers} workers)")
    print("Parallelization plan: artifact-level error decomposition tasks are distributed across workers; cross-study aggregation, standalone pulse audits, and figure generation remain local.")

    normalized = load_normalized_results()
    save_snapshots(normalized)
    combined = _deduplicate_rows(combined_dataframe(normalized))
    if args.quick and not args.full:
        combined = combined[
            combined["study"].isin({"ideal_sqr_direct_vs_echoed_multitone", "native_rich_multitone_sqr_cpsqr_feasibility"})
        ].copy()

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    NEW_RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    VALIDATION_DIR.mkdir(parents=True, exist_ok=True)

    combined.to_csv(DATA_DIR / "combined_results.csv", index=False)
    new_results = combined[combined["study"] == "native_rich_multitone_sqr_cpsqr_feasibility"].copy()
    new_results.to_csv(NEW_RESULTS_DIR / "full_grid_results.csv", index=False)
    best_per_case = new_results.sort_values("strict_process_fidelity", ascending=False).drop_duplicates(subset=["case_id", "construction"])
    best_per_case.to_csv(NEW_RESULTS_DIR / "best_per_case.csv", index=False)

    tasks = build_tasks(combined.to_dict(orient="records"))
    task_results = run_tasks(tasks, run_task, n_workers=args.n_workers, sequential=args.sequential, desc="error diagnostics")
    error_rows = []
    worker_failures = []
    for result in task_results:
        if result.ok:
            error_rows.extend(result.value)
        else:
            worker_failures.append({"artifact": result.item.artifact_path, "error": result.error})
    error_df = pd.DataFrame(error_rows)
    error_df.to_csv(NEW_RESULTS_DIR / "error_decomposition.csv", index=False)
    _save_json(NEW_RESULTS_DIR / "error_decomposition.json", {"rows": error_rows, "worker_failures": worker_failures})

    xpi_payload = compute_xpi_characterization()
    _save_json(NEW_RESULTS_DIR / "xpi_characterization.json", xpi_payload)
    xpi_df = xpi_dataframe(xpi_payload)

    spectral_payload = compute_spectral_crowding()
    _save_json(NEW_RESULTS_DIR / "spectral_crowding.json", spectral_payload)

    scaling_payload = compute_scaling_summary(combined)
    _save_json(NEW_RESULTS_DIR / "scaling_analysis.json", scaling_payload)

    figure_stems = generate_all_figures(
        combined,
        error_df,
        figure_root=FIGURE_ROOT,
        xpi_df=xpi_df,
        spectral_payload=spectral_payload,
        scaling_payload=scaling_payload,
    )
    limits = _fundamental_limits(combined)
    _save_json(NEW_RESULTS_DIR / "fundamental_limits.json", limits)
    construction_summary = _construction_summary(combined, error_df)
    _save_json(NEW_RESULTS_DIR / "construction_summary.json", construction_summary)
    _save_json(
        NEW_RESULTS_DIR / "sensitivity_analysis.json",
        {
            "status": "not_run",
            "reason": "This definitive pass builds the practical ranking from the saved artifact corpus and does not yet replay the best constructions under parameter drifts or open-system dynamics.",
        },
    )

    master_table = _master_results_table(combined, error_df)
    master_table.to_csv(NEW_RESULTS_DIR / "master_results_table.csv", index=False)

    validation = write_validation_outputs(combined, output_dir=VALIDATION_DIR, source_root=STUDY_PATHS["study1"].parent)
    wall_clock_seconds = perf_counter() - started
    objective_summaries = _objective_summaries(
        combined=combined,
        new_results=new_results,
        xpi_payload=xpi_payload,
        scaling_payload=scaling_payload,
        figure_stems=figure_stems,
        wall_clock_seconds=wall_clock_seconds,
    )
    _save_json(NEW_RESULTS_DIR / "objective_summaries.json", objective_summaries)

    structured_native = new_results[new_results["target_family"].isin({"smooth_x", "staggered_x"})].copy()
    structured_best_mean = float(
        structured_native.sort_values("strict_process_fidelity", ascending=False).drop_duplicates(subset=["case_id"])["strict_process_fidelity"].mean()
    )
    summary = {
        "headline_findings": [
            f"Best loaded strict ideal-SQR process fidelity: {combined['strict_process_fidelity'].max():.5f}.",
            f"Original baseline best direct average gate fidelity remained {combined[(combined['study'] == 'ideal_sqr_direct_vs_echoed_multitone') & (combined['construction'] == 'direct_multitone')]['avg_gate_fidelity'].max():.4f}.",
            f"Native-rich best strict process reached {new_results['strict_process_fidelity'].max():.5f}, confirming that strict ideal SQR is achievable on selected addressed windows.",
            f"Native-rich best CPSQR process reached {new_results['cpsqr_process_fidelity'].max():.9f}, confirming that echoed gains are strongest for the relaxed conditional-phase target.",
            f"The structured native-rich best-by-case mean strict fidelity is {structured_best_mean:.4f}, so the current corpus still falls short of the user's 0.95 grid-average landmark.",
        ],
        "figure_stems": figure_stems,
        "validation": validation,
        "worker_failures": worker_failures,
        "wall_clock_seconds": wall_clock_seconds,
        "objective_summaries_path": str((NEW_RESULTS_DIR / 'objective_summaries.json').resolve()),
    }
    _save_json(NEW_RESULTS_DIR / "study_summary.json", summary)
    _write_notebook(summary)
    _write_report(
        combined=combined,
        new_results=new_results,
        master_table=master_table,
        xpi_payload=xpi_payload,
        scaling_payload=scaling_payload,
    )
    compiled = _compile_report()
    print(f"Report compiled: {compiled}")
    print(f"Wall clock seconds: {summary['wall_clock_seconds']:.2f}")


if __name__ == "__main__":
    main()
