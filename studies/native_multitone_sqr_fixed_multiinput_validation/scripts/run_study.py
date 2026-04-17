from __future__ import annotations

import json
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from common import (
    CASE_SPECS,
    DATA_DIR,
    FIGURES_DIR,
    PROBE_QUBIT_STATES,
    PROBE_TIERS,
    REPORT_DIR,
    SCRIPTS_DIR,
)

# Import the common helpers explicitly to keep the linter and runtime simple.
from common import (
    apply_target_operator_to_state,
    average_gate_fidelity,
    build_corrected_run_config,
    build_model,
    build_target_operator,
    build_waveform,
    compile_waveform,
    ensure_directories,
    initial_state_tier_rows,
    leakage_outside_logical,
    legacy_restricted_process_fidelity,
    load_json,
    logical_levels,
    matrix_to_record,
    process_fidelity,
    restricted_operator_from_full,
    set_plot_style,
    simulation_session,
    simulate_full_operator_on_logical_inputs,
    single_manifold_qubit_state,
    state_fidelity,
    write_csv,
    write_json,
)


SUMMARY_JSON = DATA_DIR / "study_summary.json"
SUMMARY_MD = DATA_DIR / "study_summary.md"
PROBE_CSV = DATA_DIR / "probe_state_results.csv"
CASE_JSON = DATA_DIR / "case_results.json"
NOTEBOOK_PATH = SCRIPTS_DIR / "reproducibility_notebook.ipynb"
REPORT_TEX = REPORT_DIR / "report.tex"
REFERENCES_BIB = REPORT_DIR / "references.bib"


def evaluate_case(case) -> dict[str, Any]:
    artifact = load_json(case.artifact_path)
    model = build_model(case, artifact)
    run_config = build_corrected_run_config(case, artifact, model)
    levels = logical_levels(case, artifact)
    target_operator = build_target_operator(case, artifact, levels)
    waveform = build_waveform(case, artifact, model, run_config)
    compiled = compile_waveform(waveform, run_config)
    session = simulation_session(model, compiled, frame=run_config.frame, drive_ops=waveform.drive_ops)
    full_operator = simulate_full_operator_on_logical_inputs(
        model,
        compiled,
        frame=run_config.frame,
        drive_ops=waveform.drive_ops,
        levels=levels,
    )
    actual_restricted = restricted_operator_from_full(full_operator, model, levels)

    probe_rows: list[dict[str, Any]] = []
    for level in levels:
        for probe_label, qubit_state in PROBE_QUBIT_STATES.items():
            initial = single_manifold_qubit_state(model, int(level), qubit_state)
            actual = session.run(initial).final_state
            ideal = apply_target_operator_to_state(model, levels, target_operator, initial)
            probe_rows.append(
                {
                    "study_key": case.study_key,
                    "plot_label": case.plot_label,
                    "level": int(level),
                    "probe_label": str(probe_label),
                    "state_fidelity": state_fidelity(actual, ideal),
                    "leakage_outside_logical": leakage_outside_logical(actual, model, levels),
                }
            )

    tier_summary = initial_state_tier_rows(probe_rows)
    fixed_process = process_fidelity(target_operator, actual_restricted)
    fixed_avg_gate = average_gate_fidelity(target_operator, actual_restricted)
    legacy_process = legacy_restricted_process_fidelity(case, artifact)

    summary = {
        "study_key": case.study_key,
        "plot_label": case.plot_label,
        "study_goal": case.study_goal,
        "target_summary": case.target_summary,
        "notes": case.notes,
        "n_active": int(len(levels)),
        "duration_ns": float(run_config.duration_s * 1.0e9),
        "fixed_process_fidelity": float(fixed_process),
        "fixed_average_gate_fidelity": float(fixed_avg_gate),
        "legacy_restricted_process_fidelity": None if legacy_process is None else float(legacy_process),
        "process_fidelity_shift": None if legacy_process is None else float(fixed_process - legacy_process),
        "tier_summary": tier_summary,
        "probe_rows": probe_rows,
        "target_operator": matrix_to_record(target_operator),
        "actual_restricted_operator": matrix_to_record(actual_restricted),
    }
    return summary


def flat_probe_rows(results: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for case_result in results:
        base = {
            "study_key": case_result["study_key"],
            "plot_label": case_result["plot_label"],
            "n_active": case_result["n_active"],
            "duration_ns": case_result["duration_ns"],
        }
        for row in case_result["probe_rows"]:
            merged = dict(base)
            merged.update(row)
            rows.append(merged)
    return rows


def make_tier_bar_chart(results: list[dict[str, Any]]) -> None:
    set_plot_style()
    labels = [row["plot_label"] for row in results]
    x = np.arange(len(results))
    width = 0.18
    fig, ax = plt.subplots(figsize=(9.0, 4.6))
    single = [row["tier_summary"]["single_ground"]["mean_state_fidelity"] for row in results]
    pair = [row["tier_summary"]["selected_pair"]["mean_state_fidelity"] for row in results]
    quartet = [row["tier_summary"]["spanning_quartet"]["mean_state_fidelity"] for row in results]
    proc = [row["fixed_process_fidelity"] for row in results]
    ax.bar(x - 1.5 * width, single, width=width, label="single |g,n>")
    ax.bar(x - 0.5 * width, pair, width=width, label="selected pair")
    ax.bar(x + 0.5 * width, quartet, width=width, label="spanning quartet")
    ax.bar(x + 1.5 * width, proc, width=width, label="restricted process")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(0.0, 1.02)
    ax.set_ylabel("Fidelity")
    ax.set_title("Fixed-Package Validation Ladder")
    ax.legend(frameon=False, ncol=2)
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURES_DIR / f"tier_validation_summary.{suffix}", bbox_inches="tight")
    plt.close(fig)


def make_probe_heatmaps(results: list[dict[str, Any]]) -> None:
    set_plot_style()
    probe_order = ["g", "e", "plus_x", "plus_y"]
    fig, axes = plt.subplots(2, 2, figsize=(10.0, 7.0), constrained_layout=True)
    axes = axes.flatten()
    for ax, case_result in zip(axes, results, strict=True):
        levels = list(range(int(case_result["n_active"])))
        grid = np.zeros((len(probe_order), len(levels)), dtype=float)
        for row in case_result["probe_rows"]:
            grid[probe_order.index(row["probe_label"]), levels.index(int(row["level"]))] = float(row["state_fidelity"])
        im = ax.imshow(grid, vmin=0.0, vmax=1.0, cmap="viridis", aspect="auto")
        ax.set_title(case_result["plot_label"])
        ax.set_xticks(range(len(levels)))
        ax.set_xticklabels([str(level) for level in levels])
        ax.set_yticks(range(len(probe_order)))
        ax.set_yticklabels(probe_order)
        ax.set_xlabel("Fock level n")
        ax.set_ylabel("Probe input")
        for row_index in range(grid.shape[0]):
            for col_index in range(grid.shape[1]):
                value = grid[row_index, col_index]
                text_color = "white" if value < 0.65 else "black"
                ax.text(col_index, row_index, f"{value:.3f}", ha="center", va="center", color=text_color, fontsize=8)
    cbar = fig.colorbar(im, ax=axes.tolist(), shrink=0.95)
    cbar.set_label("Full-state fidelity")
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURES_DIR / f"probe_state_heatmaps.{suffix}", bbox_inches="tight")
    plt.close(fig)


def make_legacy_vs_fixed_chart(results: list[dict[str, Any]]) -> None:
    set_plot_style()
    labels = [row["plot_label"] for row in results]
    x = np.arange(len(results))
    width = 0.28
    legacy = [np.nan if row["legacy_restricted_process_fidelity"] is None else row["legacy_restricted_process_fidelity"] for row in results]
    fixed = [row["fixed_process_fidelity"] for row in results]
    fig, ax = plt.subplots(figsize=(8.8, 4.6))
    ax.bar(x - 0.5 * width, legacy, width=width, label="saved study metric")
    ax.bar(x + 0.5 * width, fixed, width=width, label="fixed-package reevaluation")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(0.0, 1.02)
    ax.set_ylabel("Restricted process fidelity")
    ax.set_title("Saved Versus Fixed Reevaluated Process Fidelity")
    ax.legend(frameon=False)
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURES_DIR / f"legacy_vs_fixed_process.{suffix}", bbox_inches="tight")
    plt.close(fig)


def markdown_summary(results: list[dict[str, Any]]) -> str:
    lines = ["# Native Multitone SQR Fixed Multi-Input Validation", ""]
    for row in results:
        tiers = row["tier_summary"]
        lines.append(f"## {row['plot_label'].replace(chr(10), ' ')}")
        lines.append(f"- Goal: {row['study_goal']}")
        lines.append(f"- Target: {row['target_summary']}")
        lines.append(f"- Addressed levels: `{row['n_active']}`")
        lines.append(f"- Duration: `{row['duration_ns']:.1f} ns`")
        lines.append(f"- Fixed restricted process fidelity: `{row['fixed_process_fidelity']:.6f}`")
        if row["legacy_restricted_process_fidelity"] is not None:
            lines.append(f"- Saved study restricted process fidelity: `{row['legacy_restricted_process_fidelity']:.6f}`")
            lines.append(f"- Shift after fixed reevaluation: `{row['process_fidelity_shift']:+.6f}`")
        lines.append(f"- Single-input mean / min: `{tiers['single_ground']['mean_state_fidelity']:.6f}` / `{tiers['single_ground']['min_state_fidelity']:.6f}`")
        lines.append(f"- Pair-input mean / min: `{tiers['selected_pair']['mean_state_fidelity']:.6f}` / `{tiers['selected_pair']['min_state_fidelity']:.6f}`")
        lines.append(f"- Quartet mean / min: `{tiers['spanning_quartet']['mean_state_fidelity']:.6f}` / `{tiers['spanning_quartet']['min_state_fidelity']:.6f}`")
        lines.append("")
    return "\n".join(lines)


def write_report(results: list[dict[str, Any]]) -> None:
    def latex_escape(text: str) -> str:
        return (
            text.replace("\\", "\\textbackslash{}")
            .replace("&", "\\&")
            .replace("_", "\\_")
            .replace("%", "\\%")
            .replace("#", "\\#")
        )

    best = max(results, key=lambda row: row["tier_summary"]["spanning_quartet"]["mean_state_fidelity"])
    worst = min(results, key=lambda row: row["tier_summary"]["spanning_quartet"]["mean_state_fidelity"])
    table_rows = "\n".join(
        [
            (
                f"{latex_escape(row['plot_label'].replace(chr(10), ' / '))} & "
                f"{row['fixed_process_fidelity']:.4f} & "
                f"{row['tier_summary']['single_ground']['mean_state_fidelity']:.4f} & "
                f"{row['tier_summary']['selected_pair']['mean_state_fidelity']:.4f} & "
                f"{row['tier_summary']['spanning_quartet']['mean_state_fidelity']:.4f} & "
                f"{row['tier_summary']['spanning_quartet']['min_state_fidelity']:.4f} \\\\"
            )
            for row in results
        ]
    )
    target_rows = "\n".join(
        [
            (
                f"{latex_escape(row['plot_label'].replace(chr(10), ' / '))} & "
                f"{row['n_active']} & "
                f"{row['duration_ns']:.1f} & "
                f"{latex_escape(row['target_summary'])} \\\\"
            )
            for row in results
        ]
    )
    tex = f"""
\\documentclass[aps,pra,twocolumn,reprint,amsmath,amssymb]{{revtex4-2}}
\\usepackage{{graphicx}}
\\usepackage{{booktabs}}
\\usepackage{{siunitx}}
\\usepackage[colorlinks=true,linkcolor=blue,citecolor=blue,urlcolor=blue]{{hyperref}}
\\usepackage{{float}}
\\usepackage{{xcolor}}

\\begin{{document}}
\\title{{Native Multitone SQR Validation Under the Fixed \\texttt{{cqed\\_sim}} Package}}
\\author{{Codex}}
\\affiliation{{Autonomous cQED Study Workspace}}
\\date{{\\today}}

\\begin{{abstract}}
We re-evaluated representative native multitone selective-qubit-rotation pulse artifacts from the 4 most recent SQR studies after the recent \\texttt{{cqed\\_sim}} fixes. The central question was not whether a pulse can match one favorable initial state, but whether the same native multitone construction reaches the correct target state for multiple qubit inputs on each addressed Fock manifold. We therefore used a validation ladder with three tiers: a single favorable input $|g,n\\rangle$, a small selected set $\\{{|g,n\\rangle, |+_x,n\\rangle\\}}$, and a spanning quartet $\\{{|g,n\\rangle, |e,n\\rangle, |+_x,n\\rangle, |+_y,n\\rangle\\}}$. The quartet is a practical reduced-state proxy for manifold-resolved qubit process validation because a single state does not determine the full SU(2) action. Each saved pulse was rebuilt under the fixed-package interpretation of \\texttt{{fock\\_fqs\\_hz}} and re-simulated with current \\texttt{{cqed\\_sim}}. Among the four representative native multitone cases, the strongest multi-input result was the corrected-convention unitary-optimized study, whose quartet mean fidelity was {best['tier_summary']['spanning_quartet']['mean_state_fidelity']:.4f} with restricted process fidelity {best['fixed_process_fidelity']:.4f}. The weakest quartet result was {worst['plot_label'].replace(chr(10), ' / ')}, at {worst['tier_summary']['spanning_quartet']['mean_state_fidelity']:.4f}. The main conclusion is asymmetric: native multitone SQR now has convincing evidence only in the explicitly corrected and re-optimized study, while the older saved pulse instances do not generally survive multi-input validation under the fixed simulator.
\\end{{abstract}}

\\maketitle

\\section{{Introduction}}
The recent SQR convention fixes in \\texttt{{cqed\\_sim}} changed the interpretation of the conditioned multitone layer and exposed at least one additional configuration issue involving \\texttt{{fock\\_fqs\\_hz}}. The user asked for a direct answer to the target-state question: can a native multitone SQR pulse reliably drive multiple initial states to the intended final target state? Earlier studies had already reported useful operator and state-transfer metrics, but not a uniform multi-input validation standard across the recent SQR study cluster.

\\section{{System and Methods}}
\\subsection{{Target-State Validation Ladder}}
For each addressed Fock level $n$, we tested the native multitone pulse on
\\begin{{align}}
\\text{{Tier 1}} &: |g,n\\rangle, \\\\
\\text{{Tier 2}} &: \\{{|g,n\\rangle, |+_x,n\\rangle\\}}, \\\\
\\text{{Tier 3}} &: \\{{|g,n\\rangle, |e,n\\rangle, |+_x,n\\rangle, |+_y,n\\rangle\\}}.
\\end{{align}}
The full-state fidelity for a probe input $|\\psi_{{\\mathrm{{in}}}}\\rangle$ was
\\begin{{align}}
F_\\psi = \\left|\\langle \\psi_{{\\mathrm{{target}}}} | \\psi_{{\\mathrm{{actual}}}} \\rangle\\right|^2,
\\end{{align}}
where the target state was generated by embedding the study's intended logical target operator into the full qubit-cavity Hilbert space. The restricted process fidelity was
\\begin{{align}}
F_\\mathrm{{proc}} = \\frac{{|\\mathrm{{Tr}}(U_\\mathrm{{target}}^\\dagger U_\\mathrm{{actual}})|^2}}{{d^2}},
\\end{{align}}
with $d=2N_{{\\mathrm{{active}}}}$.

\\subsection{{Corrected Rerun Setup}}
Every representative pulse was rebuilt with the current \\texttt{{cqed\\_sim}} simulator. The important configuration correction was to let the package infer the in-frame manifold frequencies internally instead of passing the older frame-shifted \\texttt{{fock\\_fqs\\_hz}} override. This preserves the saved waveform ansatz while reevaluating it under the fixed simulator semantics.

\\subsection{{Representative Cases}}
\\begin{{table*}}[t]
\\centering
\\caption{{Representative native multitone cases revalidated in this study.}}
\\begin{{tabular}}{{@{{}} l c c l @{{}}}}
\\toprule
Study & $N_{{\\mathrm{{active}}}}$ & Duration (ns) & Target summary \\\\
\\midrule
{target_rows}
\\bottomrule
\\end{{tabular}}
\\end{{table*}}

\\section{{Results}}
Table~\\ref{{tab:tiers}} shows the three validation tiers and the restricted process fidelity. The quartet mean is the main reduced target-state metric because it is the first tier that meaningfully constrains the full conditioned qubit action.

\\begin{{table*}}[t]
\\centering
\\caption{{Fixed-package validation metrics.}}
\\label{{tab:tiers}}
\\begin{{tabular}}{{@{{}} l c c c c c @{{}}}}
\\toprule
Study & $F_\\mathrm{{proc}}$ & Tier 1 mean & Tier 2 mean & Tier 3 mean & Tier 3 min \\\\
\\midrule
{table_rows}
\\bottomrule
\\end{{tabular}}
\\end{{table*}}

\\begin{{figure}}[H]
\\includegraphics[width=\\columnwidth]{{../figures/tier_validation_summary.pdf}}
\\caption{{Single-state, selected-pair, spanning-quartet, and restricted-process fidelities for the four representative native multitone cases.}}
\\end{{figure}}

\\begin{{figure}}[H]
\\includegraphics[width=\\columnwidth]{{../figures/probe_state_heatmaps.pdf}}
\\caption{{Full-state fidelities for each probe input and addressed Fock level. The spanning-quartet tier directly exposes where a pulse only works for favorable inputs.}}
\\end{{figure}}

\\begin{{figure}}[H]
\\includegraphics[width=\\columnwidth]{{../figures/legacy_vs_fixed_process.pdf}}
\\caption{{Saved study restricted-process fidelities versus fixed-package reevaluations. The corrected-convention study has no separate legacy value because it was already built around the fixed package.}}
\\end{{figure}}

The corrected-convention unitary-optimized study is the only representative case that remains convincingly strong on the spanning quartet. The two older arbitrary-rotation studies and the ideal-SQR direct study still show a visible gap between a favorable single-input check and the stronger multi-input tier, which weakens the case that those saved native multitone pulses had already demonstrated robust target-state preparation under the fixed simulator.

\\section{{Validation}}
\\subsection{{Sanity Checks}}
The probe-state ladder itself is the main analytic sanity check: one state is insufficient to identify a conditioned qubit action, while the spanning quartet is a practical reduced-state stand-in for process validation on each addressed manifold.

\\subsection{{Convergence Analysis}}
This study reused the timestep and truncation choices of the source studies and focused on reevaluation under the fixed simulator semantics. The corrected-convention source study already included its own timestep sweep and is therefore the best-converged positive example in this set.

\\subsection{{Literature Comparison}}
No direct literature benchmark was used because the task was an internal consistency check across recent study artifacts under the corrected package.

\\section{{Discussion}}
The cross-study picture is not that native multitone SQR is impossible. Rather, the evidence level depends strongly on which study generated the pulse and what objective constrained it. A single favorable input can overstate success, especially when the saved objective mainly constrained one state image or an incomplete reduced proxy. The corrected-convention study stands out because it explicitly optimized a reduced effective-unitary metric aligned with the manifold-resolved qubit action and therefore survives the stronger quartet check much better.

\\section{{Conclusion}}
Based on the fixed-package reevaluation, we do not yet have equally strong multi-input evidence from all four recent studies. The strongest convincing evidence comes from the corrected-convention direct multitone study. The older saved pulse instances do not uniformly demonstrate robust, input-consistent target-state preparation under the fixed simulator, so the safer conclusion is conditional: native multitone SQR can work, but the evidence is convincing only when the pulse is optimized against a metric that constrains the intended conditioned qubit action more fully than a single-state target.

\\section{{Limitations and Future Work}}
\\subsection{{Known Limitations}}
This study revalidated representative native multitone artifacts instead of re-optimizing every historical case from scratch. That means some failures may reflect older optimization settings in addition to the fixed simulator semantics.

\\subsection{{Suggested Improvements}}
\\begin{{itemize}}
  \\item \\textbf{{[P1 | MEDIUM]}} Re-run the older optimization grids with the fixed \\texttt{{fock\\_fqs\\_hz}} handling and compare the resulting pulse parameters against the saved legacy artifacts.
  \\item \\textbf{{[P2 | MEDIUM]}} Expand the probe-state reevaluation from representative cases to a larger random and hard-case sample, especially in the arbitrary-target studies.
  \\item \\textbf{{[P2 | LOW]}} Add a uniform reduced effective-unitary extractor to \\texttt{{cqed\\_sim}} so that all future SQR studies can report quartet validation and reduced process fidelity from one public API.
\\end{{itemize}}

\\subsection{{Open Questions}}
How much of the remaining multi-input gap in the older studies is caused by legacy optimization conventions, and how much is a real limitation of the native multitone ansatz at the chosen durations?

\\bibliographystyle{{apsrev4-2}}
\\bibliography{{references}}

\\appendix

\\section{{Detailed Results and Data}}
The machine-readable results are saved in \\texttt{{data/study\\_summary.json}}, \\texttt{{data/case\\_results.json}}, and \\texttt{{data/probe\\_state\\_results.csv}}.

\\section{{Reproducibility}}
The notebook \\texttt{{scripts/reproducibility\\_notebook.ipynb}} loads the saved validation outputs and shows how to re-run the study script.

\\end{{document}}
"""
    REPORT_TEX.write_text(tex.strip() + "\n", encoding="utf-8")
    if not REFERENCES_BIB.exists():
        REFERENCES_BIB.write_text("% No external citations used in this validation study.\n", encoding="utf-8")


def write_notebook(results: list[dict[str, Any]]) -> None:
    case_lines = "".join([f"        '{row['plot_label'].replace(chr(10), ' / ')}',\n" for row in results])
    notebook = {
        "cells": [
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "# Native Multitone SQR Fixed Multi-Input Validation\n",
                    "\n",
                    "This notebook reproduces the saved cross-study multi-input validation outputs for the four recent native multitone SQR studies.\n",
                ],
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "from pathlib import Path\n",
                    "import json, csv\n",
                    "\n",
                    "study_root = Path.cwd().resolve().parents[0] if Path.cwd().name == 'scripts' else Path('studies/native_multitone_sqr_fixed_multiinput_validation').resolve()\n",
                    "data_dir = study_root / 'data'\n",
                    "figures_dir = study_root / 'figures'\n",
                    "print('study_root =', study_root)\n",
                    "print('data_dir =', data_dir)\n",
                ],
            },
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "## User-Tunable Parameters\n",
                    "\n",
                    "These parameters control which saved case summaries are displayed and which probe tiers are emphasized.\n",
                ],
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "config = {\n",
                    "    'preferred_tier': 'spanning_quartet',\n",
                    "    'probe_tiers': ['single_ground', 'selected_pair', 'spanning_quartet'],\n",
                    "    'case_labels': [\n",
                    case_lines,
                    "    ],\n",
                    "}\n",
                    "config\n",
                ],
            },
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "## Derived Objects\n",
                    "\n",
                    "Load the saved study summary and per-probe CSV. Re-running this cell after changing the parameters above is enough for the display cells below.\n",
                ],
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "# --- Load saved results (default) ---\n",
                    "summary = json.loads((data_dir / 'study_summary.json').read_text())\n",
                    "with (data_dir / 'probe_state_results.csv').open() as handle:\n",
                    "    probe_rows = list(csv.DictReader(handle))\n",
                    "len(summary['cases']), len(probe_rows)\n",
                    "\n",
                    "# --- Re-run with current parameters ---\n",
                    "# import subprocess, sys\n",
                    "# subprocess.run([sys.executable, str(study_root / 'scripts' / 'run_study.py')], check=True)\n",
                ],
            },
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "## Saved Tier Metrics\n",
                    "\n",
                    "Display the key reduced target-state metrics for each representative case.\n",
                ],
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "for case in summary['cases']:\n",
                    "    tier = case['tier_summary'][config['preferred_tier']]\n",
                    "    print(case['plot_label'].replace('\\n', ' / '))\n",
                    "    print('  fixed process fidelity =', round(case['fixed_process_fidelity'], 6))\n",
                    "    print('  preferred tier mean   =', round(tier['mean_state_fidelity'], 6))\n",
                    "    print('  preferred tier min    =', round(tier['min_state_fidelity'], 6))\n",
                    "    print('')\n",
                ],
            },
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "## Key Figures\n",
                    "\n",
                    "The study script already generated the publication figures. This cell lists them for quick inspection.\n",
                ],
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "sorted(path.name for path in figures_dir.glob('*.pdf'))\n",
                ],
            },
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "## Summary\n",
                    "\n",
                    "The saved results distinguish favorable single-input success from stronger spanning-quartet evidence. The corrected-convention direct multitone case is the strongest positive result in the current set.\n",
                ],
            },
        ],
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python", "version": "3.12"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    NOTEBOOK_PATH.write_text(json.dumps(notebook, indent=2), encoding="utf-8")


def main() -> None:
    ensure_directories()
    results = [evaluate_case(case) for case in CASE_SPECS]
    probe_rows = flat_probe_rows(results)

    study_summary = {
        "study_name": "native_multitone_sqr_fixed_multiinput_validation",
        "probe_tiers": {name: list(labels) for name, labels in PROBE_TIERS.items()},
        "probe_labels": list(PROBE_QUBIT_STATES.keys()),
        "cases": results,
    }

    write_json(SUMMARY_JSON, study_summary)
    write_json(CASE_JSON, {"cases": results})
    write_csv(
        PROBE_CSV,
        probe_rows,
        fieldnames=[
            "study_key",
            "plot_label",
            "n_active",
            "duration_ns",
            "level",
            "probe_label",
            "state_fidelity",
            "leakage_outside_logical",
        ],
    )
    SUMMARY_MD.write_text(markdown_summary(results), encoding="utf-8")

    make_tier_bar_chart(results)
    make_probe_heatmaps(results)
    make_legacy_vs_fixed_chart(results)
    write_report(results)
    write_notebook(results)


if __name__ == "__main__":
    main()
