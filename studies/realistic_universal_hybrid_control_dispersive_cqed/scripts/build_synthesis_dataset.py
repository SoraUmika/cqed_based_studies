from __future__ import annotations

import math
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from common import (
    ARTIFACTS_DIR,
    DATA_DIR,
    FIGURES_DIR,
    PROJECT_ROOT,
    chi_scaled_duration,
    load_json,
    save_json,
)


CHI_OVER_2PI_MHZ = -2.84
CHI_PRIME_OVER_2PI_MHZ = -0.021
KERR_OVER_2PI_MHZ = -0.028


@dataclass(frozen=True)
class PrimitiveVerdict:
    primitive: str
    ideal_claim: str
    best_structured_family: str
    best_metric_label: str
    best_metric_value: float
    best_duration_ns: float
    chi_t_over_2pi: float
    verdict: str
    main_obstruction: str
    source_path: str


def _load_sources() -> dict:
    waveform_state = load_json(
        PROJECT_ROOT
        / "studies"
        / "waveform_level_gate_realization_dispersive_cqed"
        / "study_state.json"
    )
    waveform_cross_regime = load_json(
        PROJECT_ROOT
        / "studies"
        / "waveform_level_gate_realization_dispersive_cqed"
        / "artifacts"
        / "cross_regime_summary.json"
    )
    literature = load_json(
        PROJECT_ROOT
        / "studies"
        / "literature_informed_selective_primitives"
        / "data"
        / "study_results.json"
    )
    native_rich = load_json(
        PROJECT_ROOT
        / "studies"
        / "native_rich_multitone_sqr_cpsqr_feasibility"
        / "data"
        / "study_summary.json"
    )
    strong_validation = load_json(
        PROJECT_ROOT
        / "studies"
        / "strong_validation_arbitrary_fock_conditional_rotations"
        / "data"
        / "study_summary.json"
    )
    hybrid_runtime_summary_path = (
        PROJECT_ROOT
        / "task_runs"
        / "hybrid_unitary_native_entangling_evolution"
        / "EXECUTION_SUMMARY.md"
    )
    hybrid_runtime_summary = hybrid_runtime_summary_path.read_text(encoding="utf-8")
    return {
        "waveform_state": waveform_state,
        "waveform_cross_regime": waveform_cross_regime,
        "literature": literature,
        "native_rich": native_rich,
        "strong_validation": strong_validation,
        "hybrid_runtime_summary_markdown": hybrid_runtime_summary,
        "hybrid_runtime_summary_path": str(hybrid_runtime_summary_path.relative_to(PROJECT_ROOT)),
    }


def _hybrid_runtime_metrics() -> dict:
    # These values are copied from the accessible execution summary because the
    # source is narrative markdown rather than structured JSON.
    return {
        "candidate": "R2_exact_runtime_to_exact_runtime",
        "process_fidelity": 0.90269,
        "average_probe_fidelity": 0.83495,
        "average_leakage": 0.04318,
        "nominal_noise_average_probe_fidelity": 0.61448,
        "wigner_overlap": 0.96384,
        "duration_ns": 4432.0,
        "source_path": "task_runs/hybrid_unitary_native_entangling_evolution/EXECUTION_SUMMARY.md",
        "interpretation": (
            "Even the best replay-backed hybrid unitary still needs GRAPE-derived "
            "local surrogates and is not hardware-ready once runtime and noise are included."
        ),
    }


def _build_phase_budget() -> list[dict]:
    durations_ns = {
        "fast_displacement": 20.0,
        "vacuum_pi_rotation": 40.0,
        "best_noisy_sqr": 352.11267605633795,
        "best_snap_pi": 774.6478873239435,
        "strict_sqr_threshold": 1056.3380281690138,
        "best_strict_sqr_case": 1760.5633802816897,
        "runtime_hybrid_sequence": 4432.0,
    }
    rows: list[dict] = []
    for label, duration_ns in durations_ns.items():
        for n in range(6):
            qubit_cycles = abs(n * CHI_OVER_2PI_MHZ + n * (n - 1) * CHI_PRIME_OVER_2PI_MHZ) * duration_ns * 1e-3
            kerr_cycles = abs(KERR_OVER_2PI_MHZ) * n * max(n - 1, 0) * duration_ns * 1e-3
            rows.append(
                {
                    "label": label,
                    "duration_ns": duration_ns,
                    "fock_n": n,
                    "qubit_branch_phase_cycles": qubit_cycles,
                    "cavity_kerr_phase_cycles": kerr_cycles,
                }
            )
    return rows


def _primitive_verdicts(sources: dict) -> list[PrimitiveVerdict]:
    waveform_state = sources["waveform_state"]
    literature = sources["literature"]
    native_rich = sources["native_rich"]
    strong_validation = sources["strong_validation"]

    displacement = PrimitiveVerdict(
        primitive="Unconditional cavity displacement",
        ideal_claim=r"$I_q \otimes D(\alpha)$ across both qubit branches",
        best_structured_family="Two-tone branch-compensated displacement",
        best_metric_label="Broad-state mean fidelity",
        best_metric_value=0.9857,
        best_duration_ns=20.0,
        chi_t_over_2pi=chi_scaled_duration(20.0, CHI_OVER_2PI_MHZ),
        verdict="Approximate success on short calibrated pulses",
        main_obstruction=(
            "Branch detuning from chi and finite-pulse phase accumulation; naive echo fails once the "
            "inserted pi pulse becomes manifold dependent."
        ),
        source_path="studies/waveform_level_gate_realization_dispersive_cqed/study_state.json",
    )
    uqr = PrimitiveVerdict(
        primitive="Unconditional qubit rotation",
        ideal_claim=r"$R_q(\theta,\phi)$ independent of cavity occupancy",
        best_structured_family="Short Gaussian/DRAG vacuum-calibrated pulse",
        best_metric_label="Vacuum X_pi fidelity",
        best_metric_value=0.999841,
        best_duration_ns=40.0,
        chi_t_over_2pi=chi_scaled_duration(40.0, CHI_OVER_2PI_MHZ),
        verdict="Only spectator-limited, not truly unconditional",
        main_obstruction=(
            "Photon-number-dependent qubit detuning through chi and chi-prime turns the pulse into a "
            "Fock-selective operation as cavity occupation grows."
        ),
        source_path="studies/waveform_level_gate_realization_dispersive_cqed/study_state.json",
    )
    relaxed_sqr = PrimitiveVerdict(
        primitive="Selective qubit rotation (relaxed CPSQR / cphase-SQR)",
        ideal_claim=r"Blockwise qubit rotation up to removable conditional phases",
        best_structured_family="Short flat-top Gaussian cphase-SQR",
        best_metric_label="Noisy relaxed fidelity",
        best_metric_value=literature["sqr"]["best_noisy_overall"]["noisy_relaxed_avg_state_fidelity"],
        best_duration_ns=literature["sqr"]["best_noisy_overall"]["duration_us"] * 1e3,
        chi_t_over_2pi=literature["sqr"]["best_noisy_overall"]["chi_t"],
        verdict="Robust practical primitive",
        main_obstruction=(
            "Strict logical gauge is not preserved; a later SNAP or virtual-Z cleanup layer is still needed."
        ),
        source_path="studies/literature_informed_selective_primitives/data/study_results.json",
    )
    strict_sqr = PrimitiveVerdict(
        primitive="Strict full-joint ideal SQR",
        ideal_claim=r"$\sum_n |n\rangle\langle n| \otimes R_x(\theta_n)$ with correct inter-manifold phase",
        best_structured_family="Direct reduced-unitary multitone",
        best_metric_label="Joint process fidelity",
        best_metric_value=native_rich["representative_rows"]["strict_best"]["strict_joint_process_fidelity"],
        best_duration_ns=native_rich["representative_rows"]["strict_best"]["duration_ns"],
        chi_t_over_2pi=native_rich["representative_rows"]["strict_best"]["chi_t_over_2pi"],
        verdict="Limited to easy low-window cases",
        main_obstruction=(
            "The same waveform must fit all target angles while cancelling manifold-dependent coherent Z "
            "structure; chi-prime collapses the strict success window as N_active grows."
        ),
        source_path="studies/native_rich_multitone_sqr_cpsqr_feasibility/data/study_summary.json",
    )
    arbitrary_strict = PrimitiveVerdict(
        primitive="Arbitrary Fock-conditional blockwise SU(2)",
        ideal_claim="General strict blockwise SU(2) control across addressed manifolds",
        best_structured_family="Single-pulse Gaussian on structured in-plane subclass",
        best_metric_label="Strict joint process fidelity",
        best_metric_value=strong_validation["best_strict_case"]["strict_joint_process_fidelity"],
        best_duration_ns=strong_validation["best_strict_case"]["duration_ns"],
        chi_t_over_2pi=strong_validation["best_strict_case"]["chi_t_over_2pi"],
        verdict="No general success on hard classes",
        main_obstruction=(
            "Stress and random targets fail even with richer direct families; the obstruction is coherent and "
            "not explained away by Gaussian ansatz rigidity."
        ),
        source_path="studies/strong_validation_arbitrary_fock_conditional_rotations/data/study_summary.json",
    )
    arbitrary_relaxed = PrimitiveVerdict(
        primitive="Arbitrary Fock-conditional control under left-Z gauge",
        ideal_claim="General blockwise SU(2) up to explicit per-block Z gauge",
        best_structured_family="Segmented relaxed CPSQR",
        best_metric_label="Relaxed joint process fidelity",
        best_metric_value=strong_validation["best_relaxed_case"]["relaxed_joint_process_fidelity"],
        best_duration_ns=strong_validation["best_relaxed_case"]["duration_ns"],
        chi_t_over_2pi=strong_validation["best_relaxed_case"]["chi_t_over_2pi"],
        verdict="Strong positive result under explicit gauge relaxation",
        main_obstruction=(
            "The same pulse does not preserve the strict operator; success depends on acknowledging and then "
            "cleaning up the blockwise Z structure."
        ),
        source_path="studies/strong_validation_arbitrary_fock_conditional_rotations/data/study_summary.json",
    )
    snap = PrimitiveVerdict(
        primitive="Selective number-dependent arbitrary phase (SNAP)",
        ideal_claim="Exact cavity phase cleanup layer",
        best_structured_family="Flat-top Gaussian two-pi geometric SNAP",
        best_metric_label="Noisy average state fidelity",
        best_metric_value=literature["snap"]["best_noisy_overall"]["noisy_avg_state_fidelity"],
        best_duration_ns=literature["snap"]["best_noisy_overall"]["total_duration_us"] * 1e3,
        chi_t_over_2pi=literature["snap"]["best_noisy_overall"]["chi_t"],
        verdict="Useful but slow cleanup primitive",
        main_obstruction=(
            "Selectivity requires long total sequence time, so coherence cost is substantially higher than for "
            "practical relaxed SQR."
        ),
        source_path="studies/literature_informed_selective_primitives/data/study_results.json",
    )
    return [displacement, uqr, relaxed_sqr, strict_sqr, arbitrary_strict, arbitrary_relaxed, snap]


def _make_timescale_figure(verdicts: list[PrimitiveVerdict], hybrid_runtime: dict) -> None:
    fig, ax = plt.subplots(figsize=(10.5, 5.7))
    ax.axvspan(0.0, 0.12, color="#d9f0d3", alpha=0.55, label="Broadband / near-unconditional regime")
    ax.axvspan(0.8, 6.0, color="#fee8c8", alpha=0.55, label="Selective / number-resolved regime")
    ax.axvspan(6.0, 15.0, color="#fdd0a2", alpha=0.55, label="Long-sequence regime")

    y_positions = np.arange(len(verdicts), 0, -1)
    colors = {
        "Approximate success on short calibrated pulses": "#1b9e77",
        "Only spectator-limited, not truly unconditional": "#d95f02",
        "Robust practical primitive": "#1b9e77",
        "Limited to easy low-window cases": "#7570b3",
        "No general success on hard classes": "#e7298a",
        "Strong positive result under explicit gauge relaxation": "#66a61e",
        "Useful but slow cleanup primitive": "#a6761d",
    }
    for y, verdict in zip(y_positions, verdicts):
        ax.scatter(
            verdict.chi_t_over_2pi,
            y,
            s=90,
            color=colors[verdict.verdict],
            zorder=3,
        )
        ax.text(
            verdict.chi_t_over_2pi * 1.06,
            y,
            f"{verdict.primitive}\n{verdict.best_structured_family}",
            va="center",
            fontsize=8.5,
        )

    runtime_x = chi_scaled_duration(hybrid_runtime["duration_ns"], CHI_OVER_2PI_MHZ)
    ax.scatter(runtime_x, 0.35, s=110, color="#b2182b", marker="D", zorder=3)
    ax.text(
        runtime_x * 1.02,
        0.35,
        "Replay-backed sequence candidate\nwith GRAPE local surrogates",
        va="center",
        fontsize=8.5,
    )

    ax.set_xscale("log")
    ax.set_xlim(0.01, 20.0)
    ax.set_ylim(-0.2, len(verdicts) + 1.2)
    ax.set_yticks([])
    ax.set_xlabel(r"Normalized duration $|\chi| T / 2\pi$")
    ax.set_title("The same dispersive shift creates incompatible timing demands")
    ax.grid(True, axis="x", alpha=0.25, which="both")
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURES_DIR / f"timescale_hierarchy.{suffix}", dpi=300 if suffix == "png" else None)
    plt.close(fig)


def _make_phase_budget_figure(phase_budget: list[dict]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.8), sharex=True)
    palette = {
        "fast_displacement": "#1b9e77",
        "vacuum_pi_rotation": "#66a61e",
        "best_noisy_sqr": "#7570b3",
        "best_snap_pi": "#e6ab02",
        "strict_sqr_threshold": "#d95f02",
        "best_strict_sqr_case": "#e7298a",
        "runtime_hybrid_sequence": "#b2182b",
    }
    labels = {
        "fast_displacement": "20 ns two-tone displacement",
        "vacuum_pi_rotation": "40 ns vacuum qubit pi",
        "best_noisy_sqr": "352 ns relaxed SQR optimum",
        "best_snap_pi": "775 ns SNAP pi pulse",
        "strict_sqr_threshold": "1056 ns strict-SQR threshold",
        "best_strict_sqr_case": "1761 ns best strict-SQR case",
        "runtime_hybrid_sequence": "4432 ns replay-backed sequence",
    }
    for key in labels:
        rows = [row for row in phase_budget if row["label"] == key]
        n_vals = [row["fock_n"] for row in rows]
        q_vals = [row["qubit_branch_phase_cycles"] for row in rows]
        k_vals = [row["cavity_kerr_phase_cycles"] for row in rows]
        axes[0].plot(n_vals, q_vals, marker="o", color=palette[key], label=labels[key], linewidth=2)
        axes[1].plot(n_vals, k_vals, marker="o", color=palette[key], linewidth=2)
    axes[0].set_title("Qubit branch phase from dispersive detuning")
    axes[1].set_title("Cavity Kerr phase budget")
    axes[0].set_ylabel("Accumulated phase [cycles]")
    for ax in axes:
        ax.set_xlabel("Cavity photon number n")
        ax.grid(True, alpha=0.25)
        ax.set_xticks(range(6))
    axes[0].legend(loc="upper left", fontsize=7.5, frameon=False)
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURES_DIR / f"phase_budget.{suffix}", dpi=300 if suffix == "png" else None)
    plt.close(fig)


def _primitive_table_tex(verdicts: list[PrimitiveVerdict]) -> str:
    def tex_escape(text: str) -> str:
        return (
            text.replace("\\", r"\textbackslash{}")
            .replace("_", r"\_")
            .replace("%", r"\%")
        )

    lines = [
        r"\begin{tabular}{p{2.5cm}p{2.8cm}p{1.6cm}p{1.35cm}p{2.7cm}}",
        r"\toprule",
        r"Primitive & Best structured family & Metric & $|\chi|T/2\pi$ & Verdict \\",
        r"\midrule",
    ]
    for verdict in verdicts:
        metric = f"{tex_escape(verdict.best_metric_label)} {verdict.best_metric_value:.4f}"
        lines.append(
            f"{tex_escape(verdict.primitive)} & {tex_escape(verdict.best_structured_family)} & "
            f"{metric} & {verdict.chi_t_over_2pi:.3f} & {tex_escape(verdict.verdict)} \\\\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}"])
    return "\n".join(lines)


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    ARTIFACTS_DIR.mkdir(parents=True, exist_ok=True)

    sources = _load_sources()
    hybrid_runtime = _hybrid_runtime_metrics()
    phase_budget = _build_phase_budget()
    verdicts = _primitive_verdicts(sources)

    max_safe_fock = (
        sources["waveform_cross_regime"]["artifact_payload"]["qubit_rotation"]["max_fock_levels"]
    )

    summary = {
        "device_point": {
            "omega_q_over_2pi_ghz": 6.150,
            "omega_c_over_2pi_ghz": 5.241,
            "alpha_over_2pi_mhz": -255.0,
            "chi_over_2pi_mhz": CHI_OVER_2PI_MHZ,
            "chi_prime_over_2pi_mhz": CHI_PRIME_OVER_2PI_MHZ,
            "kerr_over_2pi_mhz": KERR_OVER_2PI_MHZ,
        },
        "source_catalog": {
            "waveform_level_gate_realization_dispersive_cqed": {
                "path": "studies/waveform_level_gate_realization_dispersive_cqed/study_state.json",
                "role": "Unconditional displacement, spectator qubit rotations, echoed displacement negative result",
            },
            "literature_informed_selective_primitives": {
                "path": "studies/literature_informed_selective_primitives/data/study_results.json",
                "role": "Best practical noisy SQR and geometric SNAP references",
            },
            "native_rich_multitone_sqr_cpsqr_feasibility": {
                "path": "studies/native_rich_multitone_sqr_cpsqr_feasibility/data/study_summary.json",
                "role": "Strict full-joint ideal SQR versus relaxed CPSQR",
            },
            "strong_validation_arbitrary_fock_conditional_rotations": {
                "path": "studies/strong_validation_arbitrary_fock_conditional_rotations/data/study_summary.json",
                "role": "Arbitrary blockwise SU(2) versus left-Z-gauge-relaxed CPSQR",
            },
            "hybrid_unitary_native_entangling_evolution": {
                "path": hybrid_runtime["source_path"],
                "role": "Sequence-level replay benchmark for hybrid unitary synthesis",
            },
        },
        "derived_timescales": {
            "inverse_abs_chi_over_2pi_ns": 1e3 / abs(CHI_OVER_2PI_MHZ),
            "fast_two_tone_displacement_chi_t_over_2pi": chi_scaled_duration(20.0, CHI_OVER_2PI_MHZ),
            "best_noisy_sqr_chi_t_over_2pi": sources["literature"]["sqr"]["best_noisy_overall"]["chi_t"],
            "best_snap_pi_chi_t_over_2pi": sources["literature"]["snap"]["best_noisy_overall"]["chi_t"],
            "strict_sqr_threshold_chi_t_over_2pi": 3.0,
        },
        "qubit_rotation_spectator_limits": max_safe_fock,
        "phase_budget_rows": phase_budget,
        "primitive_verdicts": [asdict(verdict) for verdict in verdicts],
        "sequence_level_runtime_reference": hybrid_runtime,
        "top_level_verdict": {
            "strict_ideal_gate_set_survives": False,
            "non_grape_structured_universal_control_demonstrated": False,
            "practical_phase_aware_route_supported": True,
            "practical_phase_aware_route_text": (
                "A realistic constructive route is supported only after replacing strict ideal SQR by "
                "phase-aware selective control: short branch-compensated displacement, short spectator-limited "
                "qubit rotations, relaxed CPSQR / cphase-SQR, and a cleanup layer such as SNAP or virtual-Z "
                "correction. Even then, a fully pulse-backed universal architecture is not yet demonstrated."
            ),
            "most_defensible_scientific_statement": (
                "The abstract universal-control claim does not survive literally once chi, chi-prime, Kerr, and "
                "finite-pulse phase accumulation are enforced. What survives is a weaker, phase-aware control "
                "library whose strongest successes are low-dimensional, gauge-relaxed, and still short of a fully "
                "validated non-GRAPE universal stack."
            ),
        },
        "sequence_level_implication": (
            "The best current replay-backed hybrid unitary candidate in the repository still relies on "
            "GRAPE-derived local surrogates and falls to 0.614 noisy average probe fidelity at 4.432 us, so "
            "even the assisted route is architectural rather than deployment-ready."
        ),
    }

    save_json(DATA_DIR / "synthesis_summary.json", summary)
    save_json(
        ARTIFACTS_DIR / "primitive_verdicts.json",
        {"primitive_verdicts": [asdict(verdict) for verdict in verdicts]},
    )
    save_json(
        ARTIFACTS_DIR / "analytic_phase_budget.json",
        {"phase_budget_rows": phase_budget},
    )
    (ARTIFACTS_DIR / "primitive_matrix_table.tex").write_text(
        _primitive_table_tex(verdicts),
        encoding="utf-8",
    )

    _make_timescale_figure(verdicts, hybrid_runtime)
    _make_phase_budget_figure(phase_budget)


if __name__ == "__main__":
    main()
