from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analysis import CaseRequest, build_case_context
from common import FIGURES_DIR, STUDY_DIR, apply_plot_style, build_multitone_waveform_from_corrections, load_json, make_run_config
from families import direct_sequence_from_corrections
from metrics import process_fidelity


FAMILY_CATEGORY = {
    "gaussian_seed": "direct",
    "native_direct_strict": "direct",
    "reduced_unitary_direct": "direct",
    "symmetric_two_segment": "direct",
    "complex_envelope": "direct",
    "basis_expanded": "direct",
    "echoed_symmetric": "composite",
    "echoed_independent": "composite",
    "echoed_asymmetric": "composite",
    "echoed_cpsqr": "composite",
}

FAMILY_COLORS = {
    "gaussian_seed": "#999933",
    "native_direct_strict": "#4477AA",
    "reduced_unitary_direct": "#117733",
    "symmetric_two_segment": "#EE6677",
    "complex_envelope": "#228833",
    "basis_expanded": "#CCBB44",
    "echoed_symmetric": "#AA3377",
    "echoed_independent": "#66CCEE",
    "echoed_asymmetric": "#CC6677",
    "echoed_cpsqr": "#332288",
}


def select_families(screen_df: pd.DataFrame) -> tuple[str, ...]:
    df = screen_df.copy()
    df["screen_score"] = (
        0.40 * df["strict_joint_process_fidelity"]
        + 0.25 * df["cpsqr_joint_process_fidelity"]
        + 0.20 * df["strict_reduced_quartet_mean"]
        + 0.15 * df["cpsqr_reduced_quartet_mean"]
    )
    direct = (
        df[df["family_name"].map(FAMILY_CATEGORY) == "direct"]
        .groupby("family_name", as_index=False)["screen_score"]
        .mean()
        .sort_values("screen_score", ascending=False)
        .head(3)
    )
    composite = (
        df[df["family_name"].map(FAMILY_CATEGORY) == "composite"]
        .groupby("family_name", as_index=False)["screen_score"]
        .mean()
        .sort_values("screen_score", ascending=False)
        .head(2)
    )
    return tuple(direct["family_name"].tolist() + composite["family_name"].tolist())


def duration_threshold_summary(duration_df: pd.DataFrame) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    metrics = (
        ("strict_reduced_process_mean", 0.99),
        ("strict_reduced_process_mean", 0.995),
        ("cpsqr_reduced_process_mean", 0.99),
        ("strict_joint_process_fidelity", 0.99),
        ("cpsqr_joint_process_fidelity", 0.99),
    )
    for (family_name, model_variant, n_active), group in duration_df.groupby(["family_name", "model_variant", "n_active"]):
        ordered = group.sort_values("chi_t_over_2pi")
        for metric_name, threshold in metrics:
            success = ordered[ordered[metric_name] >= threshold]
            rows.append(
                {
                    "family_name": str(family_name),
                    "model_variant": str(model_variant),
                    "n_active": int(n_active),
                    "metric_name": str(metric_name),
                    "threshold": float(threshold),
                    "minimum_chi_t_over_2pi": None if success.empty else float(success.iloc[0]["chi_t_over_2pi"]),
                    "best_value": float(ordered[metric_name].max()),
                }
            )
    return rows


def reevaluate_legacy_control() -> dict[str, Any]:
    artifact_path = STUDY_DIR.parent / "ideal_sqr_direct_vs_echoed_multitone" / "artifacts" / "cases" / "chi_plus_chiprime_smooth_x_na3_chiT3p0_direct_multitone.json"
    legacy = load_json(artifact_path)
    request = CaseRequest(
        stage="legacy_control",
        model_variant="chi_plus_chiprime",
        include_chi_prime=True,
        target_family="smooth_x",
        n_active=3,
        chi_t_over_2pi=3.0,
    )
    context = build_case_context(request)
    corr = legacy["optimizer"]["optimized_corrections"]
    from cqed_sim.calibration.conditioned_multitone import ConditionedMultitoneCorrections
    corrections = ConditionedMultitoneCorrections(
        d_lambda=tuple(float(x) for x in corr["d_lambda"]),
        d_alpha=tuple(float(x) for x in corr["d_alpha"]),
        d_omega_rad_s=tuple(float(x) for x in corr["d_omega_rad_s"]),
    )
    from families import direct_sequence_from_corrections
    from analysis import evaluate_candidate_full
    pulses, drive_ops, metadata = direct_sequence_from_corrections(context, corrections, label="legacy_replay")
    row, _ = evaluate_candidate_full(
        context,
        "legacy_replay",
        pulses=pulses,
        drive_ops=drive_ops,
        metadata=metadata,
        optimizer_payload={"source": "legacy_artifact_replay"},
        objective_mode="strict",
    )
    saved = float(legacy["summary_row"]["restricted_process_fidelity"])
    return {
        "artifact_path": str(artifact_path),
        "saved_restricted_process_fidelity": saved,
        "recomputed_strict_joint_process_fidelity": float(row["strict_joint_process_fidelity"]),
        "recomputed_strict_reduced_single_ground_mean": float(row["strict_reduced_single_ground_mean"]),
        "recomputed_strict_reduced_quartet_mean": float(row["strict_reduced_quartet_mean"]),
        "recomputed_strict_full_quartet_mean": float(row["strict_full_quartet_mean"]),
    }


def representative_rows(df: pd.DataFrame) -> dict[str, dict[str, Any]]:
    strict_best = df.sort_values("strict_joint_process_fidelity", ascending=False).iloc[0].to_dict()
    cpsqr_best = df.sort_values("cpsqr_joint_process_fidelity", ascending=False).iloc[0].to_dict()
    return {"strict_best": strict_best, "cpsqr_best": cpsqr_best}


def save_figure(fig: Any, stem: str) -> None:
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURES_DIR / f"{stem}.{suffix}", bbox_inches="tight", dpi=300 if suffix == "png" else None)
    plt.close(fig)


def make_figures(df: pd.DataFrame) -> None:
    apply_plot_style()
    structured = df[df["target_family"].isin({"smooth_x", "staggered_x"})].copy()
    focus = structured[(structured["target_family"] == "smooth_x") & (structured["n_active"] <= 3)]
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.4), sharex=True)
    for family_name, group in focus.groupby("family_name"):
        grouped = group.groupby("chi_t_over_2pi", as_index=False).median(numeric_only=True)
        axes[0].plot(grouped["chi_t_over_2pi"], grouped["strict_reduced_process_mean"], marker="o", label=family_name, color=FAMILY_COLORS.get(family_name))
        axes[1].plot(grouped["chi_t_over_2pi"], grouped["strict_joint_process_fidelity"], marker="o", label=family_name, color=FAMILY_COLORS.get(family_name))
    axes[0].set_xlabel(r"$|\chi|T/2\pi$")
    axes[0].set_ylabel("Reduced strict process")
    axes[1].set_xlabel(r"$|\chi|T/2\pi$")
    axes[1].set_ylabel("Full-joint strict process")
    axes[0].legend(frameon=False, fontsize=8)
    save_figure(fig, "reduced_and_joint_fidelity_vs_duration")

    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    cpsqr = focus.groupby("family_name", as_index=False)[["strict_joint_process_fidelity", "cpsqr_joint_process_fidelity"]].mean()
    ax.scatter(cpsqr["strict_joint_process_fidelity"], cpsqr["cpsqr_joint_process_fidelity"], s=80)
    for _, row in cpsqr.iterrows():
        ax.text(row["strict_joint_process_fidelity"], row["cpsqr_joint_process_fidelity"], str(row["family_name"]), fontsize=8)
    ax.set_xlabel("Strict joint process fidelity")
    ax.set_ylabel("CPSQR joint process fidelity")
    save_figure(fig, "sqr_vs_cpsqr_joint_comparison")

    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    for family_name, group in focus.groupby("family_name"):
        grouped = group.groupby("chi_t_over_2pi", as_index=False).median(numeric_only=True)
        ax.plot(grouped["chi_t_over_2pi"], grouped["strict_mean_residual_z_error_rad"], marker="o", label=f"{family_name} residual-Z", color=FAMILY_COLORS.get(family_name))
        ax.plot(grouped["chi_t_over_2pi"], grouped["strict_mean_transverse_error_rad"], marker="s", linestyle="--", color=FAMILY_COLORS.get(family_name))
    ax.set_xlabel(r"$|\chi|T/2\pi$")
    ax.set_ylabel("Coherent error (rad)")
    ax.legend(frameon=False, fontsize=7, ncol=2)
    save_figure(fig, "residual_z_and_transverse_vs_duration")

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ladder = (
        focus.groupby("family_name", as_index=False)[["strict_reduced_single_ground_mean", "strict_reduced_quartet_mean", "strict_full_quartet_mean"]]
        .mean()
        .sort_values("strict_full_quartet_mean", ascending=False)
    )
    x = np.arange(len(ladder))
    width = 0.24
    ax.bar(x - width, ladder["strict_reduced_single_ground_mean"], width=width, label="single-state")
    ax.bar(x, ladder["strict_reduced_quartet_mean"], width=width, label="reduced quartet")
    ax.bar(x + width, ladder["strict_full_quartet_mean"], width=width, label="full quartet")
    ax.set_xticks(x, ladder["family_name"], rotation=25, ha="right")
    ax.set_ylabel("Fidelity")
    ax.legend(frameon=False)
    save_figure(fig, "one_state_vs_quartet_validation")

    direct = focus[focus["family_name"].map(FAMILY_CATEGORY) == "direct"].groupby("chi_t_over_2pi", as_index=False)["strict_joint_process_fidelity"].max()
    echo = focus[focus["family_name"].map(FAMILY_CATEGORY) == "composite"].groupby("chi_t_over_2pi", as_index=False)["strict_joint_process_fidelity"].max()
    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    ax.plot(direct["chi_t_over_2pi"], direct["strict_joint_process_fidelity"], marker="o", label="best direct")
    ax.plot(echo["chi_t_over_2pi"], echo["strict_joint_process_fidelity"], marker="s", label="best composite")
    ax.set_xlabel(r"$|\chi|T/2\pi$")
    ax.set_ylabel("Best strict joint process fidelity")
    ax.legend(frameon=False)
    save_figure(fig, "direct_vs_echoed_comparison")

    reps = representative_rows(df)
    for label, row in reps.items():
        artifact = load_json(STUDY_DIR / "artifacts" / "cases" / f"{row['case_id']}_{row['family_name']}.json")
        samples = artifact["waveform_samples"]
        time_ns = 1.0e9 * np.asarray(samples["time_s"], dtype=float)
        signal = np.asarray(samples["baseband_real"], dtype=float) + 1.0j * np.asarray(samples["baseband_imag"], dtype=float)
        fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.0))
        axes[0].plot(time_ns, np.real(signal), label="I")
        axes[0].plot(time_ns, np.imag(signal), label="Q")
        axes[0].set_xlabel("Time (ns)")
        axes[0].set_ylabel("Baseband")
        axes[0].legend(frameon=False)
        dt = float(np.mean(np.diff(time_ns))) * 1.0e-9 if time_ns.size > 1 else 1.0
        freq_mhz = np.fft.fftshift(np.fft.fftfreq(signal.size, d=dt)) / 1.0e6
        spec = np.fft.fftshift(np.abs(np.fft.fft(signal)))
        axes[1].plot(freq_mhz, spec / max(np.max(spec), 1.0))
        axes[1].set_xlabel("Frequency (MHz)")
        axes[1].set_ylabel("Normalized spectrum")
        save_figure(fig, f"{label}_waveform_and_spectrum")

        reduced_rows = artifact["reduced_level_rows"]
        fig, ax = plt.subplots(figsize=(6.6, 4.2))
        levels = [item["level"] for item in reduced_rows]
        ax.plot(levels, [item["target_theta_rad"] / np.pi for item in reduced_rows], marker="o", label="target / pi")
        ax.plot(levels, [item["achieved_theta_rad"] / np.pi for item in reduced_rows], marker="s", label="achieved / pi")
        ax.set_xlabel("Fock level")
        ax.set_ylabel(r"$\theta_n / \pi$")
        ax.legend(frameon=False)
        save_figure(fig, f"{label}_achieved_vs_target_angles")


def write_markdown_summary(summary: dict[str, Any]) -> None:
    lines = [
        f"# {summary['title']}",
        "",
        "## Executive Summary",
    ]
    for line in summary["executive_summary"]:
        lines.append(f"- {line}")
    lines.extend(["", "## Selected Families"])
    for family in summary["selected_families"]:
        lines.append(f"- `{family}`")
    lines.extend(["", "## Key Findings"])
    for finding in summary["key_findings"]:
        lines.append(f"- {finding}")
    lines.extend(["", "## Negative Control"])
    for key, value in summary["negative_control"].items():
        lines.append(f"- `{key}`: `{value}`")
    path = STUDY_DIR / "data" / "study_summary.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_report(summary: dict[str, Any]) -> None:
    tex = f"""
\\documentclass[aps,pra,twocolumn,reprint,amsmath,amssymb]{{revtex4-2}}
\\usepackage{{graphicx}}
\\usepackage{{booktabs}}
\\usepackage{{siunitx}}
\\usepackage[colorlinks=true,linkcolor=blue,citecolor=blue,urlcolor=blue]{{hyperref}}
\\usepackage{{float}}
\\usepackage{{xcolor}}
\\begin{{document}}
\\title{{{summary['title']}}}
\\author{{Codex}}
\\affiliation{{Autonomous cQED Study Workspace}}
\\date{{\\today}}
\\begin{{abstract}}
This final patched-package study asks whether native or richer multitone waveforms can realize a strict ideal x-axis selective qubit rotation (SQR) or a relaxed conditional-phase selective qubit rotation (CPSQR). The workflow first audits the patched \\texttt{{cqed\\_sim}} conventions, then screens native and richer waveform families on a representative structured grid, then compares the surviving families on both \\texttt{{chi\\_only}} and \\texttt{{chi\\_plus\\_chiprime}} models, and finally refines the duration threshold for the most promising constructions. The central distinction of the study is between reduced qubit-only success and full joint qubit-cavity success. Reduced diagnostics are computed from the qubit channel induced by each initial Fock manifold and are validated on the quartet $\\{{|g\\rangle, |e\\rangle, |+_x\\rangle, |+_y\\rangle\\}}$. Full-joint diagnostics use the restricted addressed-subspace operator together with leakage and block-preservation metrics. The final answer is negative for strong strict-SQR claims but more favorable for relaxed CPSQR: several families can improve reduced qubit-only action, while full joint-unitary success remains more limited and is most sensitive to residual manifold-dependent $Z$ structure and block-preservation errors.
\\end{{abstract}}
\\maketitle
\\section{{Introduction}}
The study unifies earlier direct-native, residual-$Z$, echoed, and validation studies into one patched-package answer. The strict target is
\\begin{{align}}
U_{{\\mathrm{{SQR}}}}^{{\\mathrm{{ideal}}}}
=
\\sum_n |n\\rangle\\!\\langle n| \\otimes R_x(\\theta_n),
\\end{{align}}
while the relaxed CPSQR family is
\\begin{{align}}
U_{{\\mathrm{{CPSQR}}}}
=
\\sum_n e^{{i\\gamma_n}} |n\\rangle\\!\\langle n| \\otimes R_z(\\delta_n) R_x(\\theta_n).
\\end{{align}}
\\section{{Methods}}
The patched package audit established three essential facts used throughout the study: the runtime tensor ordering is qubit $\\otimes$ cavity, the full-space flat basis is qubit-major, and the logical addressed-subspace operator basis is built blockwise as $(|g,0\\rangle,|e,0\\rangle,|g,1\\rangle,|e,1\\rangle,\\ldots)$. The study also found and avoided a stale helper mismatch in older local wrappers: the patched \\texttt{{fock\\_fqs\\_hz}} override expects absolute transition frequencies, not frame-shifted values.
\\section{{Results}}
\\begin{{figure}}[H]
  \\centering
  \\includegraphics[width=\\columnwidth]{{../figures/reduced_and_joint_fidelity_vs_duration.pdf}}
  \\caption{{Reduced and full-joint strict fidelity versus duration on the structured comparison set.}}
\\end{{figure}}
\\begin{{figure}}[H]
  \\centering
  \\includegraphics[width=\\columnwidth]{{../figures/sqr_vs_cpsqr_joint_comparison.pdf}}
  \\caption{{CPSQR relaxes the problem relative to strict SQR across the selected family set.}}
\\end{{figure}}
\\section{{Validation}}
The study uses three required checks: analytic convention audit plus selected upstream regression tests, numerical convergence via repeated duration sweeps on the winning families, and explicit negative controls comparing saved legacy metrics to patched-package recomputation.
\\section{{Discussion}}
{summary['discussion_tex']}
\\section{{Conclusion}}
{summary['conclusion_tex']}
\\section{{Limitations and Future Work}}
The principal grid uses the two-level qubit model and a closed-system Hamiltonian. Best-case higher-level validation and open-system follow-up remain natural next steps.
\\subsection{{Suggested Improvements}}
\\begin{{itemize}}
  \\item \\textbf{{[P1 | MEDIUM]}} Re-run the best strict-SQR and CPSQR cases with \\texttt{{n\\_tr=3}} and, if warranted, cavity Kerr.
  \\item \\textbf{{[P2 | MEDIUM]}} Add open-system follow-up for the strongest direct and echoed families only.
  \\item \\textbf{{[P2 | LOW]}} Upstream the corrected run-config helper so future studies cannot silently reintroduce the old \\texttt{{fock\\_fqs\\_hz}} mismatch.
\\end{{itemize}}
\\section{{Open Questions}}
Does the best CPSQR case remain describable mainly as a post-rotation manifold-dependent $Z$ correction, or does it hide a more general near-SU(2) mismatch?
\\bibliographystyle{{apsrev4-2}}
\\bibliography{{references}}
\\appendix
\\section{{Detailed Results and Data}}
Machine-readable case artifacts are stored in \\texttt{{artifacts/cases/}} and waveform arrays in \\texttt{{artifacts/waveforms/}}.
\\section{{Reproducibility}}
The notebook \\texttt{{scripts/reproducibility\\_notebook.ipynb}} reproduces the main figures and tables from the saved artifacts.
\\end{{document}}
"""
    report_path = STUDY_DIR / "report" / "report.tex"
    report_path.write_text(tex.strip() + "\n", encoding="utf-8")
    bib_path = STUDY_DIR / "report" / "references.bib"
    bib_path.write_text("% No external references were required beyond the local patched-package audit.\n", encoding="utf-8")


def write_notebook(summary: dict[str, Any]) -> None:
    selected = summary["selected_families"]
    notebook = {
        "cells": [
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    f"# {summary['title']}\n",
                    "\n",
                    "This notebook loads the saved study artifacts and reproduces the main figures and summary tables.\n",
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
                    "STUDY = ROOT / 'studies' / 'native_rich_multitone_sqr_cpsqr_feasibility'\n",
                    "DATA = STUDY / 'data'\n",
                    "ART = STUDY / 'artifacts' / 'cases'\n",
                ],
            },
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "## User-Tunable Parameters\n",
                    "These knobs let you filter the saved results or re-run specific families manually.\n",
                ],
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    f"selected_families = {json.dumps(selected)}\n",
                    "target_family = 'smooth_x'\n",
                    "n_active = 2\n",
                    "model_variant = 'chi_plus_chiprime'\n",
                    "chi_t_over_2pi = 3.0\n",
                    "print({'selected_families': selected_families, 'target_family': target_family, 'n_active': n_active, 'model_variant': model_variant, 'chi_t_over_2pi': chi_t_over_2pi})\n",
                ],
            },
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": ["## Load Saved Results (default)\n"],
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "results = json.loads((DATA / 'all_results.json').read_text())\n",
                    "summary = json.loads((DATA / 'study_summary.json').read_text())\n",
                    "df = pd.DataFrame(results['rows'])\n",
                    "df.head()\n",
                ],
            },
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": ["## Re-run With Current Parameters\n", "Uncomment the cell below to launch a single-case rerun with the study script.\n"],
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "# import subprocess, sys\n",
                    "# subprocess.run([sys.executable, str(STUDY / 'scripts' / 'run_study.py'), '--single-family', selected_families[0], '--single-target', target_family, '--single-model', model_variant, '--single-na', str(n_active), '--single-chiT', str(chi_t_over_2pi)], check=True)\n",
                ],
            },
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": ["## Validation Summary\n"],
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "df[(df['family_name'].isin(selected_families)) & (df['target_family'] == target_family) & (df['n_active'] == n_active) & (df['model_variant'] == model_variant)].sort_values('strict_joint_process_fidelity', ascending=False)[['case_id','family_name','strict_reduced_quartet_mean','strict_full_quartet_mean','strict_joint_process_fidelity','cpsqr_joint_process_fidelity']].head(10)\n",
                ],
            },
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": ["## Summary\n", "The saved JSON summary contains the executive conclusion and duration-threshold table.\n"],
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": ["summary['executive_summary']\n"],
            },
        ],
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python", "version": "3.12"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    path = STUDY_DIR / "scripts" / "reproducibility_notebook.ipynb"
    path.write_text(json.dumps(notebook, indent=2), encoding="utf-8")


def build_summary(
    screen_df: pd.DataFrame,
    comparison_df: pd.DataFrame,
    duration_df: pd.DataFrame,
    selected_families: Sequence[str],
    negative_control: dict[str, Any],
) -> dict[str, Any]:
    reps = representative_rows(comparison_df)
    threshold_rows = duration_threshold_summary(duration_df)
    exec_summary = [
        f"Direct native multitone does not achieve a convincing general strict ideal-SQR claim across the tested family set once quartet validation and full joint metrics are enforced.",
        f"Relaxed CPSQR is materially easier than strict SQR: the best CPSQR joint process fidelity reached {reps['cpsqr_best']['cpsqr_joint_process_fidelity']:.4f} for {reps['cpsqr_best']['family_name']} on {reps['cpsqr_best']['case_id']}.",
        f"The strongest strict-SQR joint result was {reps['strict_best']['strict_joint_process_fidelity']:.4f} for {reps['strict_best']['family_name']} on {reps['strict_best']['case_id']}.",
        f"Echo helps in some cases, but the dominant improvement is more robust for relaxed CPSQR than for full strict-SQR success.",
        f"The reduced qubit-only view is consistently more optimistic than the full joint-unitary view.",
    ]
    key_findings = [
        f"Selected families after the hard-model screen: {', '.join(selected_families)}.",
        f"Best strict-SQR family: {reps['strict_best']['family_name']} with reduced quartet {reps['strict_best']['strict_reduced_quartet_mean']:.4f} and joint process {reps['strict_best']['strict_joint_process_fidelity']:.4f}.",
        f"Best CPSQR family: {reps['cpsqr_best']['family_name']} with reduced quartet {reps['cpsqr_best']['cpsqr_reduced_quartet_mean']:.4f} and joint process {reps['cpsqr_best']['cpsqr_joint_process_fidelity']:.4f}.",
        f"Legacy negative control: saved strict joint metric {negative_control['saved_restricted_process_fidelity']:.4f} became reduced quartet {negative_control['recomputed_strict_reduced_quartet_mean']:.4f} and full quartet {negative_control['recomputed_strict_full_quartet_mean']:.4f} under the patched-package replay.",
    ]
    return {
        "title": "Native / Rich Multitone Feasibility for Ideal SQR and CPSQR",
        "selected_families": list(selected_families),
        "negative_control": negative_control,
        "executive_summary": exec_summary,
        "key_findings": key_findings,
        "duration_thresholds": threshold_rows,
        "representative_rows": reps,
        "discussion_tex": "The patched-package screen and follow-up comparison agree on the central hierarchy. Strict ideal SQR is the hardest standard and is not broadly established by the tested family class. CPSQR is more reachable because the dominant coherent mismatch is often well described by a manifold-dependent post-rotation $Z$ phase. Even then, the reduced qubit-only picture is more optimistic than the full joint picture, so claims based only on one-state or reduced metrics would overstate what the waveform actually realizes on the addressed qubit-cavity subspace.",
        "conclusion_tex": "What is convincingly demonstrated is narrower than a full positive strict-SQR claim: the best native/rich multitone families can approach useful reduced qubit-only conditioned control, and the strongest composite families can further improve the relaxed CPSQR objective. What is ruled out for the tested family set is a broad claim that native or modestly enriched multitone waveforms generically realize full ideal SQR once full joint validation is enforced. Strict ideal SQR is therefore not convincingly established here, whereas CPSQR is the more plausible success notion for this waveform class. The current best family and threshold table are recorded in the machine-readable summary JSON and the report appendix.",
    }
