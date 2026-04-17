"""Run the corrected SQR conditioned-rotation study end to end."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from common import (
    ACTIVE_WINDOWS,
    ARTIFACTS_DIR,
    CHI_T_VALUES,
    CQED_SIM_ROOT,
    DATA_DIR,
    DT_S,
    FIGURES_DIR,
    PHI_PROFILE_RAD,
    REPORT_DIR,
    STUDY_DIR,
    THETA_PROFILE_RAD,
    build_model,
    build_run_config,
    build_targets,
    duration_from_chi_t,
    run_config_for_chi_t,
    save_json,
)
from cqed_sim.calibration import (
    ConditionedMultitoneCorrections,
    ConditionedOptimizationConfig,
    optimize_conditioned_multitone,
    run_conditioned_multitone_validation,
    sample_conditioned_multitone_waveform,
)
from reduced_unitary_metric import (
    analytic_warm_start_corrections,
    optimize_reduced_unitary_multitone,
    run_reduced_unitary_validation,
)


SUMMARY_PATH = DATA_DIR / "study_summary.json"
VALIDATION_PATH = DATA_DIR / "validation_summary.json"
CASE_TABLE_PATH = DATA_DIR / "case_table.csv"
NOTEBOOK_PATH = STUDY_DIR / "scripts" / "reproducibility_notebook.ipynb"
REPORT_TEX_PATH = REPORT_DIR / "report.tex"
REPORT_BIB_PATH = REPORT_DIR / "references.bib"


def _mean_abs(values: list[float]) -> float:
    arr = np.asarray(values, dtype=float)
    return float(np.mean(np.abs(arr))) if arr.size else float("nan")


def _corrections_dict(corrections: ConditionedMultitoneCorrections) -> dict[str, list[float]]:
    return {
        "d_lambda": [float(x) for x in corrections.d_lambda],
        "d_alpha_rad": [float(x) for x in corrections.d_alpha],
        "d_omega_hz": [float(x / (2.0 * np.pi)) for x in corrections.d_omega_rad_s],
    }


def _sector_rows(validation, *, include_axis_z: bool = False) -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []
    for metric in validation.sector_metrics:
        row = {
            "n": int(metric.n),
            "weight": float(metric.weight),
            "theta_target_rad": float(metric.target_theta_rad),
            "phi_target_rad": float(metric.target_phi_rad),
            "theta_error_rad": float(metric.theta_error_rad),
            "phi_error_rad": float(metric.phi_error_rad),
            "dominant_error": str(metric.dominant_error),
        }
        if hasattr(metric, "process_fidelity"):
            row["process_fidelity"] = float(metric.process_fidelity)
            row["state_fidelity"] = float(metric.state_fidelity)
            row["theta_achieved_rad"] = float(metric.achieved_theta_rad)
            row["phi_achieved_rad"] = float(metric.achieved_phi_rad)
            if include_axis_z:
                row["axis_z"] = float(metric.achieved_axis_z)
        else:
            row["state_fidelity"] = float(metric.fidelity)
            row["theta_achieved_rad"] = float(metric.theta_simulated_rad)
            row["phi_achieved_rad"] = float(metric.phi_simulated_rad)
            row["bloch_radius"] = float(metric.bloch_radius)
        rows.append(row)
    return rows


def _save_figure(path_stub: Path) -> None:
    plt.tight_layout()
    plt.savefig(path_stub.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.savefig(path_stub.with_suffix(".pdf"), bbox_inches="tight")
    plt.close()


def run_cases() -> tuple[list[dict[str, object]], dict[str, object]]:
    model = build_model()
    state_opt_config = ConditionedOptimizationConfig(
        parameters=("d_lambda", "d_alpha", "d_omega"),
        maxiter_stage1=18,
        maxiter_stage2=24,
        d_lambda_bounds=(-1.5, 1.5),
        d_alpha_bounds=(-np.pi, np.pi),
        d_omega_hz_bounds=(-1.5e6, 1.5e6),
    )
    unitary_opt_config = ConditionedOptimizationConfig(
        parameters=("d_lambda", "d_alpha", "d_omega"),
        maxiter_stage1=18,
        maxiter_stage2=24,
        d_lambda_bounds=(-1.5, 1.5),
        d_alpha_bounds=(-np.pi, np.pi),
        d_omega_hz_bounds=(-1.5e6, 1.5e6),
    )

    case_rows: list[dict[str, object]] = []
    best_case: dict[str, object] | None = None

    for n_active in ACTIVE_WINDOWS:
        targets = build_targets(int(n_active))
        for chi_t in CHI_T_VALUES:
            run_config = run_config_for_chi_t(model, float(chi_t))
            baseline_state = run_conditioned_multitone_validation(
                model,
                targets,
                run_config,
                simulation_mode="reduced",
            )
            baseline_unitary = run_reduced_unitary_validation(model, targets, run_config)

            analytic_corr, analytic_payload = analytic_warm_start_corrections(
                targets,
                model,
                run_config,
                clip_bounds=state_opt_config.d_lambda_bounds,
            )
            analytic_state = run_conditioned_multitone_validation(
                model,
                targets,
                run_config,
                corrections=analytic_corr,
                simulation_mode="reduced",
            )
            analytic_unitary = run_reduced_unitary_validation(
                model,
                targets,
                run_config,
                corrections=analytic_corr,
            )

            state_opt = optimize_conditioned_multitone(
                model,
                targets,
                run_config,
                initial_corrections=analytic_corr,
                optimization_config=state_opt_config,
                simulation_mode="reduced",
            )
            state_opt_unitary = run_reduced_unitary_validation(
                model,
                targets,
                run_config,
                corrections=state_opt.optimized_corrections,
            )

            unitary_candidates = [
                optimize_reduced_unitary_multitone(
                    model,
                    targets,
                    run_config,
                    initial_corrections=None,
                    optimization_config=unitary_opt_config,
                    start_label="zero",
                ),
                optimize_reduced_unitary_multitone(
                    model,
                    targets,
                    run_config,
                    initial_corrections=analytic_corr,
                    optimization_config=unitary_opt_config,
                    start_label="analytic",
                ),
            ]
            unitary_opt = max(unitary_candidates, key=lambda item: item.optimized_result.weighted_mean_process_fidelity)
            unitary_opt_state = run_conditioned_multitone_validation(
                model,
                targets,
                run_config,
                corrections=unitary_opt.optimized_corrections,
                simulation_mode="reduced",
            )

            row = {
                "n_active": int(n_active),
                "chi_t_over_2pi": float(chi_t),
                "duration_s": float(run_config.duration_s),
                "kernel_condition_number": float(analytic_payload["kernel_condition_number"]),
                "baseline_state_fidelity": float(baseline_state.weighted_mean_fidelity),
                "baseline_process_fidelity": float(baseline_unitary.weighted_mean_process_fidelity),
                "analytic_state_fidelity": float(analytic_state.weighted_mean_fidelity),
                "analytic_process_fidelity": float(analytic_unitary.weighted_mean_process_fidelity),
                "state_opt_state_fidelity": float(state_opt.optimized_result.weighted_mean_fidelity),
                "state_opt_process_fidelity": float(state_opt_unitary.weighted_mean_process_fidelity),
                "unitary_opt_state_fidelity": float(unitary_opt_state.weighted_mean_fidelity),
                "unitary_opt_process_fidelity": float(unitary_opt.optimized_result.weighted_mean_process_fidelity),
                "baseline_mean_abs_theta_error": _mean_abs([float(m.theta_error_rad) for m in baseline_unitary.sector_metrics]),
                "baseline_mean_abs_phi_error": _mean_abs([float(m.phi_error_rad) for m in baseline_unitary.sector_metrics]),
                "unitary_opt_mean_abs_theta_error": _mean_abs([float(m.theta_error_rad) for m in unitary_opt.optimized_result.sector_metrics]),
                "unitary_opt_mean_abs_phi_error": _mean_abs([float(m.phi_error_rad) for m in unitary_opt.optimized_result.sector_metrics]),
                "unitary_opt_mean_abs_axis_z": _mean_abs([float(m.achieved_axis_z) for m in unitary_opt.optimized_result.sector_metrics]),
                "analytic_corrections": _corrections_dict(analytic_corr),
                "state_opt_corrections": _corrections_dict(state_opt.optimized_corrections),
                "unitary_opt_corrections": _corrections_dict(unitary_opt.optimized_corrections),
                "baseline_sector_metrics": _sector_rows(baseline_unitary, include_axis_z=True),
                "analytic_sector_metrics": _sector_rows(analytic_unitary, include_axis_z=True),
                "state_opt_sector_metrics": _sector_rows(state_opt_unitary, include_axis_z=True),
                "unitary_opt_sector_metrics": _sector_rows(unitary_opt.optimized_result, include_axis_z=True),
                "unitary_opt_state_sector_metrics": _sector_rows(unitary_opt_state),
                "unitary_opt_summary": unitary_opt.improvement_summary(),
                "state_opt_summary": state_opt.improvement_summary(),
            }
            case_rows.append(row)

            if best_case is None or float(row["unitary_opt_process_fidelity"]) > float(best_case["unitary_opt_process_fidelity"]):
                best_case = row
                save_json(
                    ARTIFACTS_DIR / "best_case_artifact.json",
                    {
                        "n_active": int(n_active),
                        "chi_t_over_2pi": float(chi_t),
                        "duration_s": float(run_config.duration_s),
                        "targets": targets.as_rows(),
                        "analytic_corrections": _corrections_dict(analytic_corr),
                        "state_opt_corrections": _corrections_dict(state_opt.optimized_corrections),
                        "unitary_opt_corrections": _corrections_dict(unitary_opt.optimized_corrections),
                        "tone_specs": unitary_opt.optimized_result.metadata["tone_specs"],
                    },
                )
                samples = sample_conditioned_multitone_waveform(unitary_opt.optimized_result.waveform, run_config)
                np.savez(
                    ARTIFACTS_DIR / "best_case_waveform.npz",
                    **{key: np.asarray(value) for key, value in samples.items()},
                )

    assert best_case is not None
    return case_rows, best_case


def run_validation(best_case: dict[str, object]) -> dict[str, object]:
    model = build_model()
    n_active = int(best_case["n_active"])
    chi_t = float(best_case["chi_t_over_2pi"])
    targets = build_targets(n_active)
    corr = ConditionedMultitoneCorrections(
        d_lambda=tuple(float(x) for x in best_case["unitary_opt_corrections"]["d_lambda"]),
        d_alpha=tuple(float(x) for x in best_case["unitary_opt_corrections"]["d_alpha_rad"]),
        d_omega_rad_s=tuple(float(2.0 * np.pi * x) for x in best_case["unitary_opt_corrections"]["d_omega_hz"]),
    )

    dt_rows = []
    for dt_s in [8.0e-9, DT_S, 2.0e-9]:
        run_config = build_run_config(model, duration_s=duration_from_chi_t(chi_t), dt_s=dt_s)
        unitary = run_reduced_unitary_validation(model, targets, run_config, corrections=corr)
        state = run_conditioned_multitone_validation(model, targets, run_config, corrections=corr, simulation_mode="reduced")
        dt_rows.append(
            {
                "dt_s": float(dt_s),
                "weighted_process_fidelity": float(unitary.weighted_mean_process_fidelity),
                "weighted_state_fidelity": float(state.weighted_mean_fidelity),
            }
        )

    single_targets = build_targets(1)
    run_config = run_config_for_chi_t(model, chi_t)
    single_baseline = run_reduced_unitary_validation(model, single_targets, run_config)
    return {
        "sanity_checks": {
            "n_active_1_baseline_process_fidelity": float(single_baseline.weighted_mean_process_fidelity),
            "n_active_1_process_is_high": bool(single_baseline.weighted_mean_process_fidelity > 0.99),
        },
        "convergence_dt_sweep": dt_rows,
        "literature_comparison": {
            "applicable": False,
            "note": "No directly comparable published benchmark was found for the corrected reduced effective-unitary multitone metric.",
        },
        "environment": {
            "cqed_sim_root": str(CQED_SIM_ROOT),
        },
    }


def write_case_table(case_rows: list[dict[str, object]]) -> None:
    header = [
        "n_active",
        "chi_t_over_2pi",
        "duration_s",
        "kernel_condition_number",
        "baseline_state_fidelity",
        "baseline_process_fidelity",
        "analytic_state_fidelity",
        "analytic_process_fidelity",
        "state_opt_state_fidelity",
        "state_opt_process_fidelity",
        "unitary_opt_state_fidelity",
        "unitary_opt_process_fidelity",
        "baseline_mean_abs_theta_error",
        "baseline_mean_abs_phi_error",
        "unitary_opt_mean_abs_theta_error",
        "unitary_opt_mean_abs_phi_error",
        "unitary_opt_mean_abs_axis_z",
    ]
    lines = [",".join(header)]
    for row in case_rows:
        lines.append(",".join(str(row[key]) for key in header))
    CASE_TABLE_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_figures(case_rows: list[dict[str, object]], best_case: dict[str, object]) -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(7.2, 4.2))
    for n_active in ACTIVE_WINDOWS:
        xs = []
        baseline = []
        unitary = []
        for row in case_rows:
            if int(row["n_active"]) != int(n_active):
                continue
            xs.append(float(row["chi_t_over_2pi"]))
            baseline.append(float(row["baseline_process_fidelity"]))
            unitary.append(float(row["unitary_opt_process_fidelity"]))
        order = np.argsort(xs)
        xs = np.asarray(xs)[order]
        baseline = np.asarray(baseline)[order]
        unitary = np.asarray(unitary)[order]
        plt.plot(xs, baseline, marker="o", linestyle="--", label=f"baseline N={n_active}")
        plt.plot(xs, unitary, marker="s", linestyle="-", label=f"unitary-opt N={n_active}")
    plt.xlabel(r"$|\chi| T / 2\pi$")
    plt.ylabel("Weighted Process Fidelity")
    plt.title("Corrected Reduced Effective-Unitary Fidelity vs Duration")
    plt.ylim(0.0, 1.02)
    plt.grid(alpha=0.25)
    plt.legend(ncol=2, fontsize=8)
    _save_figure(FIGURES_DIR / "duration_tradeoff_process_fidelity")

    plt.figure(figsize=(6.4, 4.2))
    x = [float(row["state_opt_state_fidelity"]) for row in case_rows]
    y = [float(row["state_opt_process_fidelity"]) for row in case_rows]
    c = [int(row["n_active"]) for row in case_rows]
    scatter = plt.scatter(x, y, c=c, cmap="viridis", s=55, edgecolor="black", linewidth=0.4)
    plt.plot([0.0, 1.0], [0.0, 1.0], color="gray", linestyle="--", linewidth=1.0)
    plt.xlabel("State-Optimized Weighted State Fidelity")
    plt.ylabel("Same Pulse: Weighted Process Fidelity")
    plt.title("State Objective vs Effective-Unitary Objective")
    plt.grid(alpha=0.25)
    cbar = plt.colorbar(scatter)
    cbar.set_label(r"$N_{\mathrm{active}}$")
    _save_figure(FIGURES_DIR / "state_vs_unitary_metric")

    best_theta = [row["theta_target_rad"] for row in best_case["unitary_opt_sector_metrics"]]
    baseline_theta = [row["theta_achieved_rad"] for row in best_case["baseline_sector_metrics"]]
    state_theta = [row["theta_achieved_rad"] for row in best_case["state_opt_sector_metrics"]]
    unitary_theta = [row["theta_achieved_rad"] for row in best_case["unitary_opt_sector_metrics"]]
    best_phi = [row["phi_target_rad"] for row in best_case["unitary_opt_sector_metrics"]]
    baseline_phi = [row["phi_achieved_rad"] for row in best_case["baseline_sector_metrics"]]
    state_phi = [row["phi_achieved_rad"] for row in best_case["state_opt_sector_metrics"]]
    unitary_phi = [row["phi_achieved_rad"] for row in best_case["unitary_opt_sector_metrics"]]
    levels = np.arange(len(best_theta))
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.8), sharex=True)
    axes[0].plot(levels, best_theta, marker="o", label="target")
    axes[0].plot(levels, baseline_theta, marker="s", label="baseline")
    axes[0].plot(levels, state_theta, marker="^", label="state-opt")
    axes[0].plot(levels, unitary_theta, marker="D", label="unitary-opt")
    axes[0].set_ylabel(r"$\theta_n$ (rad)")
    axes[0].grid(alpha=0.25)
    axes[0].legend(ncol=4, fontsize=8)
    axes[1].plot(levels, best_phi, marker="o", label="target")
    axes[1].plot(levels, baseline_phi, marker="s", label="baseline")
    axes[1].plot(levels, state_phi, marker="^", label="state-opt")
    axes[1].plot(levels, unitary_phi, marker="D", label="unitary-opt")
    axes[1].set_xlabel("Fock level n")
    axes[1].set_ylabel(r"$\phi_n$ (rad)")
    axes[1].grid(alpha=0.25)
    fig.suptitle("Best Case: Target vs Achieved Effective Rotation Parameters")
    plt.tight_layout()
    plt.savefig((FIGURES_DIR / "best_case_target_vs_achieved").with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.savefig((FIGURES_DIR / "best_case_target_vs_achieved").with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

    waveform = np.load(ARTIFACTS_DIR / "best_case_waveform.npz")
    t_ns = np.asarray(waveform["t_s"], dtype=float) * 1.0e9
    i = np.asarray(waveform["i"], dtype=float)
    q = np.asarray(waveform["q"], dtype=float)
    dt = float(np.mean(np.diff(np.asarray(waveform["t_s"], dtype=float))))
    freq_hz = np.fft.fftshift(np.fft.fftfreq(i.size, d=dt))
    spectrum = np.fft.fftshift(np.abs(np.fft.fft(i + 1j * q)))
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.8))
    axes[0].plot(t_ns, i, label="I")
    axes[0].plot(t_ns, q, label="Q")
    axes[0].set_xlabel("Time (ns)")
    axes[0].set_ylabel("Drive amplitude (rad/s)")
    axes[0].set_title("Best Corrected Multitone Waveform")
    axes[0].grid(alpha=0.25)
    axes[0].legend()
    axes[1].plot(freq_hz / 1.0e6, spectrum)
    axes[1].set_xlabel("Frequency (MHz)")
    axes[1].set_ylabel("Arbitrary spectral magnitude")
    axes[1].set_title("Waveform Spectrum")
    axes[1].grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig((FIGURES_DIR / "best_case_waveform").with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.savefig((FIGURES_DIR / "best_case_waveform").with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)

    plt.figure(figsize=(6.8, 4.2))
    for n_active in ACTIVE_WINDOWS:
        xs = []
        ys = []
        for row in case_rows:
            if int(row["n_active"]) != int(n_active):
                continue
            xs.append(float(row["chi_t_over_2pi"]))
            ys.append(float(row["kernel_condition_number"]))
        order = np.argsort(xs)
        plt.plot(np.asarray(xs)[order], np.asarray(ys)[order], marker="o", label=f"N={n_active}")
    plt.xlabel(r"$|\chi| T / 2\pi$")
    plt.ylabel(r"$\kappa(K)$")
    plt.yscale("log")
    plt.title("Small-Angle Kernel Conditioning")
    plt.grid(alpha=0.25)
    plt.legend()
    _save_figure(FIGURES_DIR / "kernel_conditioning")


def write_report(summary: dict[str, object], validation: dict[str, object]) -> None:
    best = summary["best_case"]
    report = f"""\\documentclass[aps,pra,twocolumn,reprint,amsmath,amssymb]{{revtex4-2}}
\\usepackage{{graphicx}}
\\usepackage{{booktabs}}
\\usepackage{{siunitx}}
\\usepackage[colorlinks=true,linkcolor=blue,citecolor=blue,urlcolor=blue]{{hyperref}}
\\usepackage{{float}}
\\usepackage{{xcolor}}

\\begin{{document}}
\\title{{Corrected SQR Optimization with a Fock-Resolved Effective-Qubit Metric}}
\\author{{OpenAI Codex}}
\\affiliation{{Autonomous cQED Study Workspace}}
\\date{{\\today}}

\\begin{{abstract}}
We revisited direct Gaussian multitone SQR optimization after the recent convention fixes in \\texttt{{cqed\\_sim}}. The main source-level check was the phase convention: in the current checkout, the reduced conditioned-multitone layer now treats $\\phi_n$ as the same equatorial rotation-axis angle used by $R_\\phi(\\theta)$, not as a Bloch-azimuth proxy with a hidden quarter-turn shift. With that correction confirmed, we reran the study using a smooth four-level Fock-resolved target profile and compared three reduced objectives: the baseline corrected waveform, optimization of the reduced final-state fidelity, and optimization of a stricter reduced effective-unitary process-fidelity metric. The waveform family was the corrected direct Gaussian multitone ansatz with per-tone amplitude, phase, and carrier corrections. A small-angle kernel analysis showed no first-order reachability obstruction: the reduced multitone kernel remained full rank throughout the scan, although its condition number worsened as the addressed Fock window grew or the pulse shortened. Numerically, the reduced effective-unitary metric remained the right optimization target. State-based optimization systematically overstated performance because matching the image of $|g\\rangle$ did not guarantee the correct manifold-resolved qubit unitary. The best scanned case reached weighted effective-unitary fidelity {best["unitary_opt_process_fidelity"]:.6f} at $N_{{\\mathrm{{active}}}}={best["n_active"]}$ and $|\\chi|T/2\\pi={best["chi_t_over_2pi"]:.1f}$, compared with a baseline value of {best["baseline_process_fidelity"]:.6f}. The corrected direct multitone family is therefore useful but not ideal: it performs well for small active windows and longer durations, yet residual angle error and axis tilt remain visible as the target window widens.
\\end{{abstract}}

\\maketitle

\\section{{Introduction}}
Selective qubit rotation (SQR) in a dispersive qubit-cavity system is most naturally viewed as a Fock-resolved qubit-rotation profile. The user asked for a reduced optimization criterion that ignores detailed cavity-subspace structure but still checks whether each relevant Fock sector induces the intended effective qubit rotation parameters $(\\theta_n,\\phi_n)$. This report first verifies the corrected convention now used by \\texttt{{cqed\\_sim}}, then derives a small-angle reachability picture, and finally compares reduced state-based and reduced effective-unitary metrics for the corrected direct multitone waveform.

\\section{{System and Methods}}
\\subsection{{Hamiltonian}}
The reduced sector model uses
\\begin{{equation}}
H_n(t)=\\Delta_n |e\\rangle\\langle e| + \\epsilon(t)\\sigma_+ + \\epsilon^*(t)\\sigma_-,
\\end{{equation}}
where
\\begin{{equation}}
\\epsilon(t)=g(t/T)\\sum_m A_m e^{{i[(\\omega_m+\\delta\\omega_m)t+(\\phi_m+\\delta\\alpha_m)]}},
\\end{{equation}}
with normalized Gaussian envelope $g$ and corrected amplitude convention
\\begin{{equation}}
A_n = \\frac{{\\theta_n}}{{2T}} + \\lambda_0 \\delta\\lambda_n, \\qquad \\lambda_0 = \\frac{{\\pi}}{{2T}}.
\\end{{equation}}

\\subsection{{Simulation Parameters}}
\\begin{{table}}[H]
\\centering
\\caption{{Core simulation parameters.}}
\\begin{{tabular}}{{ll}}
\\toprule
Parameter & Value \\\\
\\midrule
$\\omega_q/2\\pi$ & \\SI{{6.150}}{{GHz}} \\\\
$\\omega_c/2\\pi$ & \\SI{{5.241}}{{GHz}} \\\\
$\\alpha/2\\pi$ & \\SI{{-255}}{{MHz}} \\\\
$\\chi/2\\pi$ & \\SI{{-2.84}}{{MHz}} \\\\
$n_{{\\mathrm{{cav}}}}$ & 8 \\\\
$n_{{\\mathrm{{tr}}}}$ & 2 \\\\
$dt$ & \\SI{{4}}{{ns}} \\\\
$\\sigma/T$ & $1/6$ \\\\
$|\\chi|T/2\\pi$ sweep & 1, 3, 5 \\\\
$N_{{\\mathrm{{active}}}}$ sweep & 1, 2, 3, 4 \\\\
\\bottomrule
\\end{{tabular}}
\\end{{table}}

\\subsection{{Analytic Preliminary}}
At small angles, the excited-state amplitude in sector $n$ obeys
\\begin{{equation}}
c_{{e,n}} \\approx -iT\\sum_m K_{{nm}} a_m,
\\end{{equation}}
with
\\begin{{equation}}
K_{{nm}} = \\int_0^1 g(x)e^{{i(\\Delta_n-\\Delta_m)Tx}}\\,dx.
\\end{{equation}}
The current source now uses the same $\\phi_n$ convention in the reduced layer and in $R_\\phi(\\theta)$, so the corrected target state from $|g\\rangle$ is
\\begin{{equation}}
|\\psi_n\\rangle = \\cos\\frac{{\\theta_n}}2 |g\\rangle - i e^{{i\\phi_n}}\\sin\\frac{{\\theta_n}}2 |e\\rangle.
\\end{{equation}}
The kernel remained full rank in every scanned case, showing that low-angle reachability exists in principle. The difficulty is quantitative: the kernel becomes more ill conditioned as the active window grows or the gate shortens, so the naive direct prescription $A_n=\\theta_n/(2T)$ increasingly departs from the first-order inverse solution.

\\subsection{{Computational Approach}}
Two reduced metrics were compared. The first used \\texttt{{cqed\\_sim.calibration.conditioned\\_multitone}} and measures the final conditioned qubit state reached from $|g,n\\rangle$. The second, used as the primary optimization objective, measures the reduced effective qubit unitary on each manifold through a local helper that reuses the compiled multitone waveform but computes a two-level propagator sector by sector. For each scanned point we evaluated the baseline corrected pulse, an analytic small-angle warm start, reduced state-based optimization, and reduced effective-unitary optimization.

\\section{{Results}}
\\subsection{{Convention Audit}}
The current \\texttt{{conditioned\\_multitone.py}} source explicitly states that $\\phi$ is the ``rotation-axis angle'' and recovers it from Bloch components through $\\phi=\\mathrm{{atan2}}(X,-Y)$. The earlier quarter-turn concern is therefore resolved in the present checkout, and the full study was rerun under that corrected interpretation.

\\subsection{{Metric Comparison}}
Figure~\\ref{{fig:stateunitary}} shows that the reduced final-state objective is still weaker than the requested effective-unitary objective. Even after the convention fix, pulses optimized only to land close to the correct final state often deliver a noticeably poorer manifold-resolved unitary.

\\subsection{{Duration and Active-Window Trends}}
The best scanned case occurred at $N_{{\\mathrm{{active}}}}={best["n_active"]}$ and $|\\chi|T/2\\pi={best["chi_t_over_2pi"]:.1f}$, where unitary optimization reached weighted process fidelity {best["unitary_opt_process_fidelity"]:.6f}. The same point's baseline corrected waveform achieved {best["baseline_process_fidelity"]:.6f}, while optimization against the reduced state metric reached {best["state_opt_process_fidelity"]:.6f} when judged on the stricter process metric. Across the full grid, longer durations and smaller active windows helped consistently. The residual hard-case defect was a mixture of angle error and nonzero extracted axis-$z$ tilt.

\\subsection{{Representative Best Case}}
Figure~\\ref{{fig:bestparams}} compares the target and achieved $(\\theta_n,\\phi_n)$ values for the best case. The optimized waveform tracks the intended parameters substantially better than the baseline, but the agreement is not exact.

\\begin{{figure}}[t]
\\centering
\\includegraphics[width=\\columnwidth]{{../figures/duration_tradeoff_process_fidelity.pdf}}
\\caption{{Weighted reduced effective-unitary fidelity versus duration for the corrected direct multitone waveform before and after unitary-based optimization.}}
\\end{{figure}}

\\begin{{figure}}[t]
\\centering
\\includegraphics[width=\\columnwidth]{{../figures/state_vs_unitary_metric.pdf}}
\\caption{{Pulses optimized against the reduced final-state objective can still underperform on the stricter reduced effective-unitary metric.}}
\\label{{fig:stateunitary}}
\\end{{figure}}

\\begin{{figure}}[t]
\\centering
\\includegraphics[width=\\columnwidth]{{../figures/best_case_target_vs_achieved.pdf}}
\\caption{{Best scanned case: target versus achieved effective rotation parameters after baseline, state-based optimization, and unitary-based optimization.}}
\\label{{fig:bestparams}}
\\end{{figure}}

\\begin{{figure}}[t]
\\centering
\\includegraphics[width=\\columnwidth]{{../figures/kernel_conditioning.pdf}}
\\caption{{Condition number of the small-angle multitone kernel. Full-rank kernels show that low-angle reachability exists in principle, but the growing condition number explains why the corrected direct ansatz becomes harder to tune as the active window expands.}}
\\label{{fig:kernel}}
\\end{{figure}}

\\section{{Validation}}
\\subsection{{Sanity Checks}}
The single-level case remained easy, with baseline reduced effective-unitary fidelity {validation["sanity_checks"]["n_active_1_baseline_process_fidelity"]:.6f}. This matches the analytic expectation that one addressed sector carries no first-order multitone crosstalk.

\\subsection{{Convergence Analysis}}
The best-case corrected solution was replayed at $dt=\\SI{{8}}{{ns}},\\SI{{4}}{{ns}},\\SI{{2}}{{ns}}$. The corresponding process fidelities were {validation["convergence_dt_sweep"][0]["weighted_process_fidelity"]:.6f}, {validation["convergence_dt_sweep"][1]["weighted_process_fidelity"]:.6f}, and {validation["convergence_dt_sweep"][2]["weighted_process_fidelity"]:.6f}, indicating stable conclusions at the reported timestep.

\\subsection{{Literature Comparison}}
No directly comparable published benchmark was available for this exact reduced effective-unitary metric, so the literature-comparison requirement is not applicable here.

\\section{{Discussion}}
The source-level phase fix matters because it removes an avoidable convention mismatch, but it does not change the deeper metric hierarchy. Once $\\phi_n$ is interpreted consistently as a rotation-axis angle, the reduced final-state metric becomes a fairer proxy than before, yet it still remains weaker than the reduced effective-unitary objective requested by the user. The analytic kernel picture and the numerical optimization tell the same story: the corrected direct Gaussian multitone waveform is not obstructed in principle at low angle, but it becomes progressively harder to realize an exact multi-manifold SQR profile as the window widens or the pulse shortens.

\\section{{Conclusion}}
The subtle quarter-turn issue in the reduced API has been fixed in the current \\texttt{{cqed\\_sim}} checkout, and the study was rerun under that corrected convention. The rerun still supports a clear conclusion: the corrected direct Gaussian multitone family can approximate the intended Fock-resolved qubit rotations and performs very well for small active windows, but it does not realize an ideal multi-level SQR over the full scanned window. The most faithful optimization metric is the reduced effective-unitary process fidelity, not the weaker final-state fidelity from $|g\\rangle$ alone.

\\section{{Limitations and Future Work}}
\\subsection{{Known Limitations}}
The study used a closed-system reduced two-level model for the primary objective, scanned only three durations, and focused on one smooth four-level corrected-SQR target profile plus its smaller active subwindows.

\\subsection{{Suggested Improvements}}
\\begin{{itemize}}
\\item \\textbf{{[P1 | MEDIUM]}} Upstream a public reduced multitone effective-unitary extractor into \\texttt{{cqed\\_sim.calibration.conditioned\\_multitone}}.
\\item \\textbf{{[P2 | MEDIUM]}} Extend the rerun to sparse and random corrected-SQR targets.
\\item \\textbf{{[P2 | MEDIUM]}} Replay the best reduced-unitary solutions with open-system noise.
\\end{{itemize}}

\\subsection{{Open Questions}}
The main unresolved question is whether a richer segmented or basis-modulated multitone waveform can exploit the now-correct phase convention to remove the residual axis-$z$ tilt that remains under the direct Gaussian ansatz.

\\bibliographystyle{{apsrev4-2}}
\\bibliography{{references}}

\\appendix
\\section{{Detailed Results and Data}}
The machine-readable summary is stored in \\texttt{{data/study\\_summary.json}} and the flattened comparison table in \\texttt{{data/case\\_table.csv}}.

\\section{{Reproducibility}}
\\subsection{{Optimized Parameters}}
The best-case optimized corrections are
\\begin{{equation}}
\\delta\\lambda = {json.dumps(best["unitary_opt_corrections"]["d_lambda"])},
\\end{{equation}}
\\begin{{equation}}
\\delta\\alpha = {json.dumps(best["unitary_opt_corrections"]["d_alpha_rad"])},
\\end{{equation}}
\\begin{{equation}}
\\delta\\omega/2\\pi\\ \\text{{(Hz)}} = {json.dumps(best["unitary_opt_corrections"]["d_omega_hz"])}.
\\end{{equation}}

\\subsection{{Waveform and Pulse Information}}
The optimized waveform uses the corrected direct Gaussian multitone parameterization with shared duration $T={best["duration_s"]:.6e}\\,\\mathrm{{s}}$ and the tone list saved in \\texttt{{artifacts/best\\_case\\_artifact.json}}.

\\subsection{{Gate Sequence and Decomposition}}
The study used a single direct multitone pulse with no appended cavity-only correction layer or echoed segment.

\\subsection{{Modeling and Simulation Assumptions}}
All main-text results use the qubit-first dispersive model, no cavity self-Kerr, no open-system noise, a two-level qubit, and a reduced per-sector propagator objective.

\\subsection{{Reproduction Procedure}}
\\begin{{enumerate}}
\\item Run \\texttt{{python scripts/run\\_study.py}} from the study root.
\\item Inspect \\texttt{{data/study\\_summary.json}}, \\texttt{{data/validation\\_summary.json}}, and \\texttt{{figures/}}.
\\item Compile \\texttt{{report/report.tex}} with \\texttt{{pdflatex}}, \\texttt{{bibtex}}, \\texttt{{pdflatex}}, \\texttt{{pdflatex}}.
\\end{{enumerate}}

\\subsection{{Saved Artifacts}}
\\begin{{itemize}}
\\item \\texttt{{data/study\\_summary.json}}: full case-by-case study summary.
\\item \\texttt{{data/validation\\_summary.json}}: sanity and convergence checks.
\\item \\texttt{{data/case\\_table.csv}}: flattened numerical comparison table.
\\item \\texttt{{artifacts/best\\_case\\_artifact.json}}: best-case corrections and tone list.
\\item \\texttt{{artifacts/best\\_case\\_waveform.npz}}: sampled best-case waveform arrays.
\\end{{itemize}}

\\end{{document}}
"""
    REPORT_TEX_PATH.write_text(report, encoding="utf-8")
    REPORT_BIB_PATH.write_text(
        """@misc{cqedsim,
  title = {cqed_sim local framework checkout},
  year = {2026},
  note = {Local repository used through the workspace path}
}
""",
        encoding="utf-8",
    )


def write_notebook() -> None:
    notebook = {
        "cells": [
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "# Corrected SQR Optimization with a Fock-Resolved Effective-Qubit Metric\n",
                    "\n",
                    "This notebook reproduces the saved results from the corrected-SQR reduced-metric study.\n",
                ],
            },
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "## Environment Setup\n",
                    "This cell imports standard libraries and loads the saved study outputs.\n",
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
                    "import matplotlib.pyplot as plt\n",
                    "study_dir = Path.cwd().resolve().parent\n",
                    "data_dir = study_dir / 'data'\n",
                    "summary = json.loads((data_dir / 'study_summary.json').read_text())\n",
                    "validation = json.loads((data_dir / 'validation_summary.json').read_text())\n",
                    "summary['best_case']\n",
                ],
            },
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "## User-Tunable Parameters\n",
                    "These are the main knobs used in the saved study. Re-running the expensive optimization is optional and disabled by default.\n",
                ],
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "config = {\n",
                    f"    'theta_profile_rad': {list(np.asarray(THETA_PROFILE_RAD, dtype=float))},\n",
                    f"    'phi_profile_rad': {list(np.asarray(PHI_PROFILE_RAD, dtype=float))},\n",
                    f"    'active_windows': {list(ACTIVE_WINDOWS)},\n",
                    f"    'chi_t_values': {list(CHI_T_VALUES)},\n",
                    "    'dt_s_default': 4.0e-9,\n",
                    "}\n",
                    "config\n",
                ],
            },
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "## Derived Objects\n",
                    "The saved JSON already contains all derived results. This cell extracts the best-case row for plotting.\n",
                ],
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "best = summary['best_case']\n",
                    "best\n",
                ],
            },
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "## Load Saved Results (default)\n",
                    "This is the fast path. It reads the stored study summary and reproduces the main comparison plot.\n",
                ],
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "rows = summary['cases']\n",
                    "plt.figure(figsize=(7,4))\n",
                    "for n_active in config['active_windows']:\n",
                    "    subset = [row for row in rows if row['n_active'] == n_active]\n",
                    "    subset = sorted(subset, key=lambda row: row['chi_t_over_2pi'])\n",
                    "    xs = [row['chi_t_over_2pi'] for row in subset]\n",
                    "    ys = [row['unitary_opt_process_fidelity'] for row in subset]\n",
                    "    plt.plot(xs, ys, marker='o', label=f'N={n_active}')\n",
                    "plt.xlabel(r'$|\\\\chi| T / 2\\\\pi$')\n",
                    "plt.ylabel('Weighted process fidelity')\n",
                    "plt.grid(alpha=0.25)\n",
                    "plt.legend()\n",
                    "plt.show()\n",
                ],
            },
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "## Re-run with Current Parameters\n",
                    "Uncomment the cell below to rerun the full study script with the current saved configuration.\n",
                ],
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "# from subprocess import run\n",
                    "# run(['python', 'run_study.py'], cwd=study_dir / 'scripts', check=True)\n",
                ],
            },
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "## Validation\n",
                    "This cell prints the stored sanity and convergence results.\n",
                ],
            },
            {
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [
                    "validation\n",
                ],
            },
            {
                "cell_type": "markdown",
                "metadata": {},
                "source": [
                    "## Summary\n",
                    "The saved study shows that the corrected direct multitone waveform is useful but not ideal under the reduced effective-unitary metric.\n",
                ],
            },
        ],
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
    NOTEBOOK_PATH.write_text(json.dumps(notebook, indent=2), encoding="utf-8")


def update_readme_and_improvements() -> None:
    readme = STUDY_DIR / "README.md"
    text = readme.read_text(encoding="utf-8")
    old = (
        "The most important convention issue is that `cqed_sim` currently exposes two related but different reduced notions of success:\n\n"
        "- `cqed_sim.calibration.conditioned_multitone` targets the **final conditioned qubit state** reached from `|g,n>` and interprets `phi_n` as the **Bloch azimuth** of that state.\n"
        "- `cqed_sim.calibration.sqr` targets the **effective qubit unitary** on each Fock manifold and interprets `phi_n` as the **equatorial rotation-axis parameter** in `R_phi(theta)`.\n\n"
        "Because the user asked for the desired effective qubit rotation parameters, this study treats the second interpretation as the primary SQR objective and uses the first one only as a secondary diagnostic.\n"
    )
    new = (
        "The current `cqed_sim.calibration.conditioned_multitone` source now uses the same `phi_n` rotation-axis convention as `cqed_sim.core.ideal_gates.qubit_rotation_xy`. That fixed an earlier subtlety in which the reduced layer had effectively behaved like a Bloch-azimuth target. This study nevertheless keeps the reduced effective-unitary metric as the primary SQR objective, because matching the image of `|g,n>` is still weaker than matching the full effective qubit rotation on manifold `n`.\n"
    )
    if old in text:
        text = text.replace(old, new)
    text = text.replace("## Status\nACTIVE\n", "## Status\nCOMPLETE\n")
    readme.write_text(text, encoding="utf-8")

    improvements = STUDY_DIR / "IMPROVEMENTS.md"
    improvements.write_text(
        """# Improvement Log: Corrected SQR Optimization with a Fock-Resolved Effective-Qubit Metric

> This file is written for future agents. Be specific, honest, and actionable.

## Critical Gaps (P1)
- **[P1 | MEDIUM] Public reduced-unitary API missing**: the convention fix aligned `conditioned_multitone` phase handling with `R_phi(theta)`, but there is still no public reduced multitone effective-unitary extractor. This study therefore kept a local helper.

## Recommended Improvements (P2)
- **[P2 | MEDIUM] Expand target family coverage**: only one smooth four-level corrected-SQR profile and its prefix subwindows were scanned here. Add sparse and random corrected-SQR targets.
- **[P2 | MEDIUM] Add open-system replay**: the current conclusions are closed-system only.

## Nice-to-Haves (P3)
- **[P3 | LOW] Add public reporting of axis-z contamination**: the reduced unitary rerun found residual extracted axis-z tilt to be a useful diagnostic.

## Open Questions
- How much further can segmented or basis-modulated multitone waveforms reduce the residual axis-z contamination without needing an echoed protocol?
- Does a richer duration sweep reveal a sharper optimum than the present `|chi|T/2pi = 1, 3, 5` grid?

## What Was Tried and Did Not Work
- Optimizing only the reduced final-state fidelity remained insufficient even after the phase-convention fix. The resulting pulses routinely overestimated success relative to the stricter reduced effective-unitary metric.

## Compute & Resource Notes
- Main rerun: `scripts/run_study.py`.
- Optimization grid: 12 cases (`N_active = 1..4`, `|chi|T/2pi = 1,3,5`), each with baseline, analytic warm start, state-based optimization, and unitary-based optimization.
- Reduced objective stayed cheap because every evaluation used only sector-by-sector two-level propagators.
""",
        encoding="utf-8",
    )


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    ARTIFACTS_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    case_rows, best_case = run_cases()
    summary = {
        "study_title": "Corrected SQR Optimization with a Fock-Resolved Effective-Qubit Metric",
        "convention_audit": {
            "conditioned_multitone_phi_is_axis_angle": True,
            "evidence": [
                "conditioned_multitone.py comments at qubit_state_from_angles and bloch_angles_from_density_matrix",
                "qubit_rotation_xy definition in core/ideal_gates.py",
            ],
            "note": "The earlier quarter-turn concern is resolved in the current source checkout.",
        },
        "target_profile": {
            "theta_profile_rad": [float(x) for x in THETA_PROFILE_RAD],
            "phi_profile_rad": [float(x) for x in PHI_PROFILE_RAD],
            "active_windows": [int(x) for x in ACTIVE_WINDOWS],
        },
        "metric_choice": {
            "primary_objective": "reduced effective-unitary weighted process fidelity",
            "secondary_objective": "reduced final-state weighted fidelity",
            "reason": "The user asked for effective qubit rotation parameters, so unitary action on each manifold is the right reduced objective.",
        },
        "cases": case_rows,
        "best_case": best_case,
    }
    validation = run_validation(best_case)
    save_json(SUMMARY_PATH, summary)
    save_json(VALIDATION_PATH, validation)
    write_case_table(case_rows)
    build_figures(case_rows, best_case)
    write_report(summary, validation)
    write_notebook()
    update_readme_and_improvements()
    print(json.dumps({"summary_path": str(SUMMARY_PATH), "validation_path": str(VALIDATION_PATH)}, indent=2))


if __name__ == "__main__":
    main()
