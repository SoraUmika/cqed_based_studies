"""Run the strict no-detuning multitone-SQR study."""

from __future__ import annotations

import argparse
import time
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import (
    ARTIFACTS_DIR,
    DATA_DIR,
    DEFAULT_DT,
    FIGURES_DIR,
    IDEAL_X_PI,
    PI_PULSE_DURATION_S,
    STUDY_DIR,
    TWO_PI,
    CaseRequest,
    apply_plot_style,
    average_gate_fidelity,
    block_rotation_metrics,
    build_frame,
    build_model,
    build_square_multitone_waveform,
    build_target_operator,
    channel_waveform_samples,
    compile_pulse_sequence,
    corrections_to_vector,
    decoupled_block_operator,
    duration_from_chi_t,
    embed_qubit_operator,
    frobenius_error,
    logical_levels,
    magnus_effective_blocks,
    make_gaussian_qubit_rotation_pulse,
    make_run_config,
    optimize_square_multitone,
    process_fidelity,
    reduced_blockwise_operator,
    restricted_blocks,
    restricted_operator_from_full,
    save_json,
    save_waveform_npz,
    shift_pulse,
    simulate_full_operator_on_logical_inputs,
    target_spec,
)


RESULTS_PATH = DATA_DIR / "study_results.json"
SUMMARY_PATH = DATA_DIR / "study_summary.json"
CSV_PATH = DATA_DIR / "study_results.csv"
MARKDOWN_PATH = DATA_DIR / "study_summary.md"
ANALYTIC_PATH = DATA_DIR / "analytic_summary.json"
VALIDATION_PATH = DATA_DIR / "validation_summary.json"
AUDIT_PATH = DATA_DIR / "prior_audit.json"

CASE_DIR = ARTIFACTS_DIR / "cases"
WAVEFORM_DIR = ARTIFACTS_DIR / "waveforms"

MODEL_VARIANTS = (
    ("chi_only", False),
    ("chi_plus_chiprime", True),
)
BASE_FAMILIES = ("aligned_x", "structured_xy")
ACTIVE_GRID = (2, 3, 4)
DURATION_GRID = (1.0, 3.0, 5.0)
RANDOM_CASES = (
    (2, 3.0, 6101),
    (2, 5.0, 6102),
    (3, 3.0, 6201),
    (3, 5.0, 6202),
)


def prior_audit_payload() -> dict[str, Any]:
    return {
        "strict_scope_mismatches": [
            {
                "study": "studies/ideal_sqr_direct_vs_echoed_multitone",
                "issue": "The main optimizer allowed d_omega, so it does not answer the no-artificial-detuning question directly.",
                "evidence": "scripts/run_study.py uses parameters=(\"d_lambda\", \"d_alpha\", \"d_omega\").",
            },
            {
                "study": "studies/multitone_sqr_arbitrary_fock_conditional_rotations",
                "issue": "Both the direct and echoed follow-up optimizers allowed d_omega.",
                "evidence": "scripts/run_study.py and scripts/run_echo_comparison.py both use parameters=(\"d_lambda\", \"d_alpha\", \"d_omega\").",
            },
            {
                "study": "studies/parameterized_waveform_residual_z_cancellation",
                "issue": "The study explored richer waveform families outside the strict simultaneous shared-line amplitude-plus-azimuth ansatz and also allowed d_omega.",
                "evidence": "README and scripts/run_study.py describe echoed/complex/basis-expanded families with d_omega optimization.",
            },
            {
                "study": "task_runs/the_definitive_ideal_sqr_gate_study",
                "issue": "Repository-level positive strict-SQR claims refer to broader native-rich or echoed constructions rather than the strict no-detuning shared-line multitone model.",
                "evidence": "Execution summary cites native-rich extensions and CPSQR echoed cases rather than the strict simultaneous ansatz studied here.",
            },
        ]
    }


def case_requests() -> list[CaseRequest]:
    rows: list[CaseRequest] = []
    for model_variant, include_chi_prime in MODEL_VARIANTS:
        for family in BASE_FAMILIES:
            for n_active in ACTIVE_GRID:
                for chi_t in DURATION_GRID:
                    rows.append(
                        CaseRequest(
                            model_variant=model_variant,
                            include_chi_prime=include_chi_prime,
                            family=family,
                            n_active=n_active,
                            chi_t_over_2pi=chi_t,
                        )
                    )
        for n_active, chi_t, seed in RANDOM_CASES:
            rows.append(
                CaseRequest(
                    model_variant=model_variant,
                    include_chi_prime=include_chi_prime,
                    family="random_xy",
                    n_active=n_active,
                    chi_t_over_2pi=chi_t,
                    seed=seed,
                )
            )
    return rows


def half_target(spec: Any) -> Any:
    return type(spec)(
        family=f"{spec.family}_half",
        theta_values=tuple(float(value / 2.0) for value in spec.theta_values),
        phi_values=tuple(float(value) for value in spec.phi_values),
        metadata={
            **dict(spec.metadata),
            "description": f"Half-angle target derived from {spec.family}",
            "parent_family": str(spec.family),
            "theta_values_rad": [float(value / 2.0) for value in spec.theta_values],
        },
    )


def row_from_operator(
    request: CaseRequest,
    *,
    construction: str,
    target_operator: np.ndarray,
    actual_operator: np.ndarray,
    extra: dict[str, Any] | None = None,
    validation: Any | None = None,
) -> dict[str, Any]:
    target_blocks = restricted_blocks(target_operator)
    actual_blocks = restricted_blocks(actual_operator)
    block_rows = [block_rotation_metrics(target, actual) for target, actual in zip(target_blocks, actual_blocks, strict=True)]
    row = {
        "case_id": str(request.case_id),
        "model_variant": str(request.model_variant),
        "include_chi_prime": bool(request.include_chi_prime),
        "family": str(request.family),
        "seed": None if request.seed is None else int(request.seed),
        "n_active": int(request.n_active),
        "chi_t_over_2pi": float(request.chi_t_over_2pi),
        "duration_s": float(duration_from_chi_t(request.chi_t_over_2pi)),
        "duration_ns": float(duration_from_chi_t(request.chi_t_over_2pi) * 1.0e9),
        "construction": str(construction),
        "restricted_process_fidelity": float(process_fidelity(target_operator, actual_operator)),
        "restricted_average_gate_fidelity": float(average_gate_fidelity(target_operator, actual_operator)),
        "restricted_fro_error": float(frobenius_error(target_operator, actual_operator)),
        "per_block_process_fidelities": [float(item["process_fidelity"]) for item in block_rows],
        "per_block_average_gate_fidelities": [float(item["average_gate_fidelity"]) for item in block_rows],
        "per_block_rotation_angle_errors_rad": [float(item["rotation_angle_error_rad"]) for item in block_rows],
        "per_block_rotation_axis_errors_rad": [float(item["rotation_axis_error_rad"]) for item in block_rows],
        "per_block_residual_z_error_rad": [float(item["residual_z_error_rad"]) for item in block_rows],
        "per_block_transverse_error_rad": [float(item["transverse_error_rad"]) for item in block_rows],
        "per_block_actual_rotvec_z_rad": [float(item["actual_rotvec_z_rad"]) for item in block_rows],
        "mean_block_average_gate_fidelity": float(np.mean([item["average_gate_fidelity"] for item in block_rows])),
        "worst_block_average_gate_fidelity": float(np.min([item["average_gate_fidelity"] for item in block_rows])),
        "mean_residual_z_error_rad": float(np.mean([item["residual_z_error_rad"] for item in block_rows])),
        "max_residual_z_error_rad": float(np.max([item["residual_z_error_rad"] for item in block_rows])),
        "mean_transverse_error_rad": float(np.mean([item["transverse_error_rad"] for item in block_rows])),
        "max_transverse_error_rad": float(np.max([item["transverse_error_rad"] for item in block_rows])),
    }
    if validation is not None:
        row.update(
            {
                "best_fit_restricted_process_fidelity": float(validation.best_fit_restricted_process_fidelity),
                "uncorrected_restricted_process_fidelity": float(validation.uncorrected_restricted_process_fidelity),
                "same_block_population_mean": float(validation.same_block_population_mean),
                "same_block_population_min": float(validation.same_block_population_min),
                "other_target_population_mean": float(validation.other_target_population_mean),
                "other_target_population_max": float(validation.other_target_population_max),
                "leakage_outside_target_mean": float(validation.leakage_outside_target_mean),
                "leakage_outside_target_max": float(validation.leakage_outside_target_max),
                "weighted_loss": float(validation.weighted_loss),
                "best_fit_block_phase_rms_rad": float(
                    float("nan")
                    if validation.block_phase_diagnostics is None
                    else validation.block_phase_diagnostics.rms_block_phase_error_rad
                ),
            }
        )
    if extra:
        row.update(extra)
    return row


def save_case_artifact(
    path: Path,
    *,
    request: CaseRequest,
    spec: Any,
    construction: str,
    target_operator: np.ndarray,
    actual_operator: np.ndarray,
    tone_specs: Sequence[Any] | None = None,
    corrections_vector: np.ndarray | None = None,
    waveform_samples: dict[str, Any] | None = None,
    metadata: dict[str, Any] | None = None,
) -> None:
    payload = {
        "case_id": str(request.case_id),
        "construction": str(construction),
        "model_variant": str(request.model_variant),
        "family": str(request.family),
        "seed": None if request.seed is None else int(request.seed),
        "n_active": int(request.n_active),
        "chi_t_over_2pi": float(request.chi_t_over_2pi),
        "target_spec": dict(spec.metadata),
        "target_operator": target_operator,
        "restricted_operator": actual_operator,
        "tone_specs": [] if tone_specs is None else [tone.as_dict() for tone in tone_specs],
        "corrections_vector": None if corrections_vector is None else np.asarray(corrections_vector, dtype=float),
        "waveform_samples": waveform_samples,
        "metadata": {} if metadata is None else dict(metadata),
    }
    save_json(path, payload, description=f"Detailed artifact for {construction} on {request.case_id}.")


def figure_duration_tradeoff(df: pd.DataFrame) -> None:
    apply_plot_style()
    subset = df[df["construction"].isin(["full_shared_line", "blockwise_exact_reduced", "decoupled_block"])].copy()
    if subset.empty:
        return
    grouped = (
        subset.groupby(["family", "construction", "chi_t_over_2pi"], as_index=False)["restricted_average_gate_fidelity"]
        .mean()
        .sort_values(["family", "construction", "chi_t_over_2pi"])
    )
    palette = {
        "full_shared_line": "#4477AA",
        "blockwise_exact_reduced": "#CCBB44",
        "decoupled_block": "#228833",
    }
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), sharey=True)
    for axis, family in zip(axes, ["aligned_x", "structured_xy"], strict=True):
        family_df = grouped[grouped["family"] == family]
        for construction in ("full_shared_line", "blockwise_exact_reduced", "decoupled_block"):
            rows = family_df[family_df["construction"] == construction]
            if rows.empty:
                continue
            axis.plot(
                rows["chi_t_over_2pi"],
                rows["restricted_average_gate_fidelity"],
                marker="o",
                linewidth=2.0,
                color=palette[construction],
                label=construction.replace("_", " "),
            )
        axis.set_title(family.replace("_", " "))
        axis.set_xlabel(r"$|\chi| T / 2\pi$")
        axis.set_ylim(0.0, 1.02)
        axis.grid(True, alpha=0.25)
    axes[0].set_ylabel("Mean restricted average gate fidelity")
    axes[1].legend(frameon=False, loc="lower left")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "duration_fidelity_tradeoff.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "duration_fidelity_tradeoff.pdf", bbox_inches="tight")
    plt.close(fig)


def figure_residual_z(df: pd.DataFrame) -> None:
    apply_plot_style()
    subset = df[df["construction"] == "full_shared_line"].copy()
    if subset.empty:
        return
    grouped = (
        subset.groupby(["family", "n_active", "chi_t_over_2pi"], as_index=False)["max_residual_z_error_rad"]
        .mean()
        .sort_values(["family", "n_active", "chi_t_over_2pi"])
    )
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), sharey=True)
    palette = {2: "#4477AA", 3: "#EE6677", 4: "#228833"}
    for axis, family in zip(axes, ["aligned_x", "structured_xy"], strict=True):
        family_df = grouped[grouped["family"] == family]
        for n_active in sorted(family_df["n_active"].unique()):
            rows = family_df[family_df["n_active"] == n_active]
            axis.plot(
                rows["chi_t_over_2pi"],
                rows["max_residual_z_error_rad"],
                marker="o",
                linewidth=2.0,
                color=palette[int(n_active)],
                label=f"N={int(n_active)}",
            )
        axis.set_title(family.replace("_", " "))
        axis.set_xlabel(r"$|\chi| T / 2\pi$")
        axis.grid(True, alpha=0.25)
    axes[0].set_ylabel("Mean max block residual-Z error [rad]")
    axes[1].legend(frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "blockwise_residual_z_vs_duration.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "blockwise_residual_z_vs_duration.pdf", bbox_inches="tight")
    plt.close(fig)


def figure_subspace_scaling(df: pd.DataFrame) -> None:
    apply_plot_style()
    subset = df[df["construction"].isin(["full_shared_line", "decoupled_block"])].copy()
    if subset.empty:
        return
    grouped = (
        subset.groupby(["construction", "n_active"], as_index=False)["restricted_average_gate_fidelity"]
        .mean()
        .sort_values(["construction", "n_active"])
    )
    palette = {"full_shared_line": "#4477AA", "decoupled_block": "#228833"}
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    for construction in ("full_shared_line", "decoupled_block"):
        rows = grouped[grouped["construction"] == construction]
        ax.plot(
            rows["n_active"],
            rows["restricted_average_gate_fidelity"],
            marker="o",
            linewidth=2.0,
            color=palette[construction],
            label=construction.replace("_", " "),
        )
    ax.set_xlabel("Active addressed blocks")
    ax.set_ylabel("Mean restricted average gate fidelity")
    ax.set_ylim(0.0, 1.02)
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, loc="lower left")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "addressed_subspace_scaling.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "addressed_subspace_scaling.pdf", bbox_inches="tight")
    plt.close(fig)


def figure_echo_comparison(df: pd.DataFrame) -> None:
    apply_plot_style()
    subset = df[df["construction"].isin(["full_shared_line", "echo_ideal_instantaneous", "echo_finite_gaussian"])].copy()
    if subset.empty:
        return
    grouped = (
        subset.groupby(["family", "construction"], as_index=False)["restricted_average_gate_fidelity"]
        .mean()
        .sort_values(["family", "construction"])
    )
    order = ["full_shared_line", "echo_ideal_instantaneous", "echo_finite_gaussian"]
    labels = ["plain", "echo ideal", "echo finite"]
    palette = ["#4477AA", "#CCBB44", "#EE6677"]
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.0), sharey=True)
    for axis, family in zip(axes, ["aligned_x", "structured_xy"], strict=True):
        family_df = grouped[grouped["family"] == family].set_index("construction")
        values = [float(family_df.loc[key, "restricted_average_gate_fidelity"]) if key in family_df.index else np.nan for key in order]
        axis.bar(labels, values, color=palette)
        axis.set_title(family.replace("_", " "))
        axis.grid(True, axis="y", alpha=0.25)
        axis.set_ylim(0.0, 1.02)
    axes[0].set_ylabel("Mean restricted average gate fidelity")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "plain_vs_echo_comparison.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "plain_vs_echo_comparison.pdf", bbox_inches="tight")
    plt.close(fig)


def representative_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    df = pd.DataFrame(rows)
    full_df = df[df["construction"] == "full_shared_line"].copy()
    dec_df = df[df["construction"] == "decoupled_block"].copy()
    echo_df = df[df["construction"] == "echo_finite_gaussian"].copy()
    ideal_echo_df = df[df["construction"] == "echo_ideal_instantaneous"].copy()
    random_df = full_df[full_df["family"] == "random_xy"].copy()
    return {
        "strict_shared_line_mean_fidelity": float(full_df["restricted_average_gate_fidelity"].mean()),
        "strict_shared_line_best_fidelity": float(full_df["restricted_average_gate_fidelity"].max()),
        "strict_shared_line_worst_fidelity": float(full_df["restricted_average_gate_fidelity"].min()),
        "strict_shared_line_mean_best_fit_block_gauge": float(full_df["best_fit_restricted_process_fidelity"].mean()),
        "strict_shared_line_mean_max_residual_z_error_rad": float(full_df["max_residual_z_error_rad"].mean()),
        "decoupled_block_mean_fidelity": float(dec_df["restricted_average_gate_fidelity"].mean()),
        "decoupled_block_min_fidelity": float(dec_df["restricted_average_gate_fidelity"].min()),
        "ideal_echo_mean_fidelity": float(ideal_echo_df["restricted_average_gate_fidelity"].mean()) if not ideal_echo_df.empty else float("nan"),
        "finite_echo_mean_fidelity": float(echo_df["restricted_average_gate_fidelity"].mean()) if not echo_df.empty else float("nan"),
        "random_family_mean_fidelity": float(random_df["restricted_average_gate_fidelity"].mean()) if not random_df.empty else float("nan"),
        "best_full_case": (
            None
            if full_df.empty
            else full_df.sort_values("restricted_average_gate_fidelity", ascending=False).iloc[0].to_dict()
        ),
        "worst_full_case": (
            None
            if full_df.empty
            else full_df.sort_values("restricted_average_gate_fidelity", ascending=True).iloc[0].to_dict()
        ),
    }


def validation_payload(case_results: list[dict[str, Any]], echoed_rows: list[dict[str, Any]]) -> dict[str, Any]:
    df = pd.DataFrame(case_results)
    all_df = pd.DataFrame(case_results + echoed_rows)
    rep_case = df[
        (df["construction"] == "full_shared_line")
        & (df["family"] == "structured_xy")
        & (df["n_active"] == 3)
        & (df["chi_t_over_2pi"] == 3.0)
        & (df["model_variant"] == "chi_plus_chiprime")
    ]
    rep = None if rep_case.empty else rep_case.iloc[0].to_dict()
    aligned_plain = all_df[(all_df["family"] == "aligned_x") & (all_df["construction"] == "full_shared_line")]
    aligned_ideal_echo = all_df[(all_df["family"] == "aligned_x") & (all_df["construction"] == "echo_ideal_instantaneous")]
    return {
        "sanity_checks": {
            "decoupled_block_exact_match": float(
                df[df["construction"] == "decoupled_block"]["restricted_average_gate_fidelity"].min()
            ),
            "full_model_vs_magnus_case_count": int(len(df[df["construction"] == "magnus_effective"])),
            "aligned_x_ideal_echo_minus_plain_mean": float(
                aligned_ideal_echo["restricted_average_gate_fidelity"].mean()
                - aligned_plain["restricted_average_gate_fidelity"].mean()
            )
            if (not aligned_plain.empty and not aligned_ideal_echo.empty)
            else float("nan"),
        },
        "convergence_representative_case": rep,
        "notes": {
            "dt_default_s": float(DEFAULT_DT),
            "pi_pulse_duration_s": float(PI_PULSE_DURATION_S),
            "representative_case_selected": rep is not None,
        },
    }


def run_case(request: CaseRequest, *, n_starts: int, maxiter: int) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    spec = target_spec(request.family, request.n_active, seed=request.seed)
    model = build_model(include_chi_prime=request.include_chi_prime, n_active=request.n_active)
    frame = build_frame(model)
    duration_s = duration_from_chi_t(request.chi_t_over_2pi)
    run_config = make_run_config(model, n_active=request.n_active, duration_s=duration_s)
    levels = logical_levels(request.n_active)
    target_operator = build_target_operator(spec, levels)

    optimization = optimize_square_multitone(
        model,
        spec,
        run_config,
        n_starts=n_starts,
        maxiter=maxiter,
        random_seed=1000 + abs(hash(request.case_id)) % 100000,
        label_prefix=request.case_id,
    )
    full_validation = optimization.validation
    full_row = row_from_operator(
        request,
        construction="full_shared_line",
        target_operator=target_operator,
        actual_operator=np.asarray(full_validation.restricted_operator, dtype=np.complex128),
        validation=full_validation,
        extra={
            "optimization_success": bool(optimization.success),
            "optimization_message": str(optimization.message),
            "n_optimizer_history": int(len(optimization.history)),
            "corrections_vector": corrections_to_vector(optimization.corrections, n_active=request.n_active).tolist(),
        },
    )
    waveform_samples = channel_waveform_samples(full_validation.compiled)
    save_case_artifact(
        CASE_DIR / f"{request.case_id}_full_shared_line.json",
        request=request,
        spec=spec,
        construction="full_shared_line",
        target_operator=target_operator,
        actual_operator=np.asarray(full_validation.restricted_operator, dtype=np.complex128),
        tone_specs=optimization.tone_specs,
        corrections_vector=corrections_to_vector(optimization.corrections, n_active=request.n_active),
        waveform_samples=waveform_samples,
        metadata={"optimizer_history_length": int(len(optimization.history))},
    )
    save_waveform_npz(WAVEFORM_DIR / f"{request.case_id}_full_shared_line.npz", waveform_samples)

    full_restricted = np.asarray(full_validation.restricted_operator, dtype=np.complex128)
    reduced_operator = reduced_blockwise_operator(
        model,
        full_validation.compiled,
        optimization.waveform,
        run_config,
        levels=levels,
    )
    reduced_row = row_from_operator(
        request,
        construction="blockwise_exact_reduced",
        target_operator=target_operator,
        actual_operator=reduced_operator,
        extra={
            "source_case": str(request.case_id),
            "reduced_vs_full_restricted_process_fidelity": float(process_fidelity(full_restricted, reduced_operator)),
            "reduced_vs_full_restricted_average_gate_fidelity": float(average_gate_fidelity(full_restricted, reduced_operator)),
        },
    )
    save_case_artifact(
        CASE_DIR / f"{request.case_id}_blockwise_exact_reduced.json",
        request=request,
        spec=spec,
        construction="blockwise_exact_reduced",
        target_operator=target_operator,
        actual_operator=reduced_operator,
        tone_specs=optimization.tone_specs,
        corrections_vector=corrections_to_vector(optimization.corrections, n_active=request.n_active),
        metadata={"source_case": str(request.case_id)},
    )

    magnus_operator = np.zeros_like(target_operator)
    for block_index, block in enumerate(
        magnus_effective_blocks(optimization.tone_specs, model=model, run_config=run_config, levels=levels)
    ):
        magnus_operator[2 * block_index : 2 * block_index + 2, 2 * block_index : 2 * block_index + 2] = block
    magnus_row = row_from_operator(
        request,
        construction="magnus_effective",
        target_operator=target_operator,
        actual_operator=magnus_operator,
        extra={
            "source_case": str(request.case_id),
            "magnus_vs_full_restricted_process_fidelity": float(process_fidelity(full_restricted, magnus_operator)),
            "magnus_vs_full_restricted_average_gate_fidelity": float(average_gate_fidelity(full_restricted, magnus_operator)),
        },
    )
    save_case_artifact(
        CASE_DIR / f"{request.case_id}_magnus_effective.json",
        request=request,
        spec=spec,
        construction="magnus_effective",
        target_operator=target_operator,
        actual_operator=magnus_operator,
        tone_specs=optimization.tone_specs,
        corrections_vector=corrections_to_vector(optimization.corrections, n_active=request.n_active),
        metadata={"source_case": str(request.case_id)},
    )

    _ideal_waveform, ideal_tone_specs = build_square_multitone_waveform(
        model,
        spec,
        run_config,
        corrections=None,
        label=f"{request.case_id}_ideal_target",
    )
    decoupled_operator = decoupled_block_operator(ideal_tone_specs, levels=levels, duration_s=duration_s)
    decoupled_row = row_from_operator(
        request,
        construction="decoupled_block",
        target_operator=target_operator,
        actual_operator=decoupled_operator,
        extra={"same_block_population_mean": 1.0, "leakage_outside_target_mean": 0.0},
    )
    save_case_artifact(
        CASE_DIR / f"{request.case_id}_decoupled_block.json",
        request=request,
        spec=spec,
        construction="decoupled_block",
        target_operator=target_operator,
        actual_operator=decoupled_operator,
        tone_specs=ideal_tone_specs,
        corrections_vector=np.zeros(2 * request.n_active, dtype=float),
    )

    echo_rows: list[dict[str, Any]] = []
    if request.family in BASE_FAMILIES and request.n_active <= 3:
        half_spec = half_target(spec)
        half_run_config = make_run_config(model, n_active=request.n_active, duration_s=0.5 * duration_s)
        half_optimization = optimize_square_multitone(
            model,
            half_spec,
            half_run_config,
            n_starts=max(2, n_starts - 1),
            maxiter=maxiter,
            random_seed=7000 + abs(hash(request.case_id)) % 100000,
            label_prefix=f"{request.case_id}_half",
        )
        half_full_operator = np.asarray(half_optimization.validation.full_operator, dtype=np.complex128)
        x_full = embed_qubit_operator(IDEAL_X_PI, n_cav=int(model.n_cav))
        ideal_echo_full = x_full @ half_full_operator @ x_full @ half_full_operator
        ideal_echo_restricted = restricted_operator_from_full(ideal_echo_full, model, levels)
        ideal_echo_row = row_from_operator(
            request,
            construction="echo_ideal_instantaneous",
            target_operator=target_operator,
            actual_operator=ideal_echo_restricted,
            extra={
                "echo_half_duration_ns": float(0.5 * duration_s * 1.0e9),
                "echo_total_duration_ns": float(duration_s * 1.0e9),
                "echo_pi_duration_ns": 0.0,
            },
        )
        save_case_artifact(
            CASE_DIR / f"{request.case_id}_echo_ideal_instantaneous.json",
            request=request,
            spec=spec,
            construction="echo_ideal_instantaneous",
            target_operator=target_operator,
            actual_operator=ideal_echo_restricted,
            tone_specs=half_optimization.tone_specs,
            corrections_vector=corrections_to_vector(half_optimization.corrections, n_active=request.n_active),
            metadata={"echo_kind": "ideal_instantaneous"},
        )

        pi_pulse = make_gaussian_qubit_rotation_pulse(
            model,
            frame,
            theta=np.pi,
            phase=0.0,
            duration_s=PI_PULSE_DURATION_S,
            manifold_level=0,
            label=f"{request.case_id}_xpi",
        )
        half_pulse = half_optimization.waveform.pulse
        finite_pulses = [
            shift_pulse(half_pulse, t0=0.0, label=f"{request.case_id}_half_1"),
            shift_pulse(pi_pulse, t0=0.5 * duration_s, label=f"{request.case_id}_xpi_1"),
            shift_pulse(half_pulse, t0=0.5 * duration_s + PI_PULSE_DURATION_S, label=f"{request.case_id}_half_2"),
            shift_pulse(pi_pulse, t0=duration_s + PI_PULSE_DURATION_S, label=f"{request.case_id}_xpi_2"),
        ]
        finite_total_duration = duration_s + 2.0 * PI_PULSE_DURATION_S
        finite_compiled = compile_pulse_sequence(
            finite_pulses,
            dt_s=DEFAULT_DT,
            total_duration_s=finite_total_duration,
        )
        finite_full = simulate_full_operator_on_logical_inputs(
            model,
            finite_compiled,
            frame=frame,
            drive_ops={"qubit": "qubit"},
            levels=levels,
        )
        finite_restricted = restricted_operator_from_full(finite_full, model, levels)
        finite_echo_row = row_from_operator(
            request,
            construction="echo_finite_gaussian",
            target_operator=target_operator,
            actual_operator=finite_restricted,
            extra={
                "echo_half_duration_ns": float(0.5 * duration_s * 1.0e9),
                "echo_total_duration_ns": float(finite_total_duration * 1.0e9),
                "echo_pi_duration_ns": float(PI_PULSE_DURATION_S * 1.0e9),
            },
        )
        finite_waveform_samples = channel_waveform_samples(finite_compiled)
        save_case_artifact(
            CASE_DIR / f"{request.case_id}_echo_finite_gaussian.json",
            request=request,
            spec=spec,
            construction="echo_finite_gaussian",
            target_operator=target_operator,
            actual_operator=finite_restricted,
            tone_specs=half_optimization.tone_specs,
            corrections_vector=corrections_to_vector(half_optimization.corrections, n_active=request.n_active),
            waveform_samples=finite_waveform_samples,
            metadata={
                "echo_kind": "finite_gaussian",
                "echo_total_duration_s": float(finite_total_duration),
            },
        )
        save_waveform_npz(WAVEFORM_DIR / f"{request.case_id}_echo_finite_gaussian.npz", finite_waveform_samples)
        echo_rows.extend([ideal_echo_row, finite_echo_row])

    return [full_row, reduced_row, magnus_row, decoupled_row], echo_rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-starts", type=int, default=3, help="Multi-start count for the strict no-detuning optimizer.")
    parser.add_argument("--maxiter", type=int, default=75, help="Powell iteration budget per start.")
    args = parser.parse_args()

    start_time = time.perf_counter()
    CASE_DIR.mkdir(parents=True, exist_ok=True)
    WAVEFORM_DIR.mkdir(parents=True, exist_ok=True)

    save_json(AUDIT_PATH, prior_audit_payload(), description="Repository-scope audit for the strict no-detuning multitone-SQR study.")

    rows: list[dict[str, Any]] = []
    echoed_rows: list[dict[str, Any]] = []
    for request in case_requests():
        case_rows, echo_case_rows = run_case(request, n_starts=args.n_starts, maxiter=args.maxiter)
        rows.extend(case_rows)
        echoed_rows.extend(echo_case_rows)

    all_rows = rows + echoed_rows
    df = pd.DataFrame(all_rows)
    df.to_csv(CSV_PATH, index=False)

    summary = representative_summary(all_rows)
    validation = validation_payload(rows, echoed_rows)
    analytic_summary = {
        "two_block_no_go": {
            "statement": "For two simultaneously addressed blocks with nonzero amplitudes lambda_0 and lambda_1, the second-order spectator-induced Z coefficients are proportional to -lambda_1^2 K(Delta,T) and +lambda_0^2 K(Delta,T). For generic K(Delta,T) != 0, both cannot vanish unless one drive amplitude is zero.",
            "kernel_definition": "K(Delta,T) = 1/Delta - sin(Delta T)/(Delta^2 T) for the square-envelope study pulse.",
        },
        "generic_many_block_claim": {
            "statement": "For more than two blocks, the spectator-induced Z coefficients form a quadratic system zeta_n = -sum_{m != n} lambda_m^2 K(Delta_nm,T). These constraints are independent of azimuth and are not satisfied on an open set of nontrivial target rotations.",
            "generic_meaning": "The exact cancellation set is a measure-zero subset of durations and amplitudes where the kernel matrix has a special null relation or some addressed amplitudes vanish.",
        },
    }

    figure_duration_tradeoff(df)
    figure_residual_z(df)
    figure_subspace_scaling(df)
    figure_echo_comparison(df)

    save_json(
        RESULTS_PATH,
        {
            "case_rows": all_rows,
            "runtime_s": float(time.perf_counter() - start_time),
            "n_base_rows": int(len(rows)),
            "n_echo_rows": int(len(echoed_rows)),
        },
        description="Full machine-readable result table for the strict no-detuning multitone-SQR study.",
    )
    save_json(SUMMARY_PATH, summary, description="Headline summary for the strict no-detuning multitone-SQR study.")
    save_json(ANALYTIC_PATH, analytic_summary, description="Analytic summary for the no-go, decoupled-block, and echo claims.")
    save_json(VALIDATION_PATH, validation, description="Sanity and convergence notes for the strict no-detuning multitone-SQR study.")

    summary_lines = [
        f"# Summary: {STUDY_DIR.name}",
        "",
        f"- Strict shared-line no-detuning mean restricted average gate fidelity: {summary['strict_shared_line_mean_fidelity']:.6f}",
        f"- Strict shared-line best restricted average gate fidelity: {summary['strict_shared_line_best_fidelity']:.6f}",
        f"- Mean strict shared-line max residual-Z error: {summary['strict_shared_line_mean_max_residual_z_error_rad']:.6f} rad",
        f"- Decoupled-block minimum restricted average gate fidelity: {summary['decoupled_block_min_fidelity']:.6f}",
        f"- Ideal instantaneous echo mean fidelity: {summary['ideal_echo_mean_fidelity']:.6f}",
        f"- Finite Gaussian echo mean fidelity: {summary['finite_echo_mean_fidelity']:.6f}",
    ]
    MARKDOWN_PATH.write_text("\n".join(summary_lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
