"""Run the multitone SQR arbitrary Fock-conditional rotation study."""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cqed_sim.calibration.conditioned_multitone import ConditionedOptimizationConfig
from cqed_sim.calibration.targeted_subspace_multitone import (
    TargetedSubspaceObjectiveWeights,
    TargetedSubspaceOptimizationConfig,
    build_spanning_state_transfer_set,
    optimize_targeted_subspace_multitone,
)

from common import (
    ARTIFACTS_DIR,
    CHI,
    CHI_PRIME,
    DATA_DIR,
    FIGURES_DIR,
    STUDY_DIR,
    SUCCESS_TIERS,
    active_subspace_metrics,
    apply_plot_style,
    block_diag_target,
    block_rotation_metrics,
    build_frame,
    build_model,
    classify_failure_mode,
    classify_success,
    conditioned_targets_from_blocks,
    crosstalk_matrix,
    crosstalk_summary,
    duration_from_chi_t,
    json_ready,
    load_json,
    logical_levels,
    make_family_blocks,
    make_run_config,
    manifold_transition_frequencies_hz,
    min_transition_spacing_hz,
    restricted_blocks,
    save_json,
    state_validation_summary,
)


RESULTS_PATH = DATA_DIR / "study_results.json"
MACHINE_SUMMARY_PATH = DATA_DIR / "machine_summary.json"
CHECKPOINT_PATH = DATA_DIR / "study_checkpoint.json"

ACTIVE_GRID = (2, 3, 4, 5)
STRUCTURED_DURATION_GRID = (1.0, 2.0, 3.0, 5.0)
RANDOM_DURATION_GRID = (1.0, 3.0, 5.0)
RANDOM_ENSEMBLE_SIZE = 4
BASE_RANDOM_SEED = 314159

MODEL_VARIANTS = (
    ("chi_only", False),
    ("chi_plus_chiprime", True),
)

STRUCTURED_FAMILIES = ("A", "B", "C")

OBJECTIVE_WEIGHTS = TargetedSubspaceObjectiveWeights(
    qubit_weight=0.15,
    subspace_weight=1.0,
    preservation_weight=0.35,
    leakage_weight=0.35,
)


@dataclass(frozen=True)
class CaseRequest:
    family: str
    model_variant: str
    include_chi_prime: bool
    n_active: int
    duration_chi_t: float
    random_seed: int | None = None

    @property
    def case_id(self) -> str:
        seed_part = "" if self.random_seed is None else f"_seed{int(self.random_seed):03d}"
        duration_label = str(self.duration_chi_t).replace(".", "p")
        return (
            f"{self.model_variant}_na{int(self.n_active)}_chiT{duration_label}_family{self.family}{seed_part}"
        )


def optimization_config(n_active: int) -> TargetedSubspaceOptimizationConfig:
    conditioned = ConditionedOptimizationConfig(
        active_levels=tuple(range(int(n_active))),
        parameters=("d_lambda", "d_alpha", "d_omega"),
        maxiter_stage1=18,
        maxiter_stage2=30,
        d_lambda_bounds=(-0.75, 0.75),
        d_alpha_bounds=(-np.pi, np.pi),
        d_omega_hz_bounds=(-3.0e6, 3.0e6),
        regularization_lambda=5.0e-4,
        regularization_alpha=5.0e-4,
        regularization_omega=5.0e-4,
    )
    return TargetedSubspaceOptimizationConfig(conditioned=conditioned, include_block_phase=False)


def case_requests() -> list[CaseRequest]:
    requests: list[CaseRequest] = []
    for model_variant, include_chi_prime in MODEL_VARIANTS:
        for n_active in ACTIVE_GRID:
            for duration_chi_t in STRUCTURED_DURATION_GRID:
                for family in STRUCTURED_FAMILIES:
                    requests.append(
                        CaseRequest(
                            family=family,
                            model_variant=model_variant,
                            include_chi_prime=include_chi_prime,
                            n_active=n_active,
                            duration_chi_t=duration_chi_t,
                        )
                    )
            for duration_chi_t in RANDOM_DURATION_GRID:
                for offset in range(RANDOM_ENSEMBLE_SIZE):
                    requests.append(
                        CaseRequest(
                            family="D",
                            model_variant=model_variant,
                            include_chi_prime=include_chi_prime,
                            n_active=n_active,
                            duration_chi_t=duration_chi_t,
                            random_seed=BASE_RANDOM_SEED + 1000 * n_active + 100 * int(include_chi_prime) + offset,
                        )
                    )
    return requests


def sample_waveform_payload(waveform, duration_s: float) -> dict[str, Any]:
    tlist = np.linspace(0.0, float(duration_s), 2000)
    samples = np.asarray(waveform.sample(tlist), dtype=np.complex128)
    return {
        "time_s": tlist.tolist(),
        "real": np.real(samples).tolist(),
        "imag": np.imag(samples).tolist(),
    }


def run_case(request: CaseRequest) -> dict[str, Any]:
    levels = logical_levels(request.n_active)
    rng = None if request.random_seed is None else np.random.default_rng(int(request.random_seed))
    blocks, target_meta = make_family_blocks(request.family, request.n_active, rng=rng)
    target_operator = block_diag_target(blocks)
    targets = conditioned_targets_from_blocks(blocks)

    model = build_model(include_chi_prime=request.include_chi_prime, n_active=request.n_active)
    frame = build_frame(model)
    duration_s = duration_from_chi_t(request.duration_chi_t)
    run_config = make_run_config(model, n_active=request.n_active, duration_s=duration_s)
    transfer_set = build_spanning_state_transfer_set(target_operator)

    start = time.perf_counter()
    optimization = optimize_targeted_subspace_multitone(
        model,
        targets,
        run_config,
        logical_levels=levels,
        optimization_config=optimization_config(request.n_active),
        objective_weights=OBJECTIVE_WEIGHTS,
        target_operator=target_operator,
        transfer_set=transfer_set,
        label=request.case_id,
    )
    runtime_s = time.perf_counter() - start

    validation = optimization.optimized_result
    metrics = active_subspace_metrics(target_operator, validation.restricted_operator)
    actual_blocks = restricted_blocks(validation.restricted_operator)
    per_block = [
        {"level": int(level), **block_rotation_metrics(target, actual)}
        for level, target, actual in zip(levels, blocks, actual_blocks, strict=True)
    ]
    residual_z = np.asarray([float(item["residual_z_error_rad"]) for item in per_block], dtype=float)
    transverse = np.asarray([float(item["transverse_error_rad"]) for item in per_block], dtype=float)

    cross = crosstalk_matrix(model, run_config, validation.waveform.tone_specs, levels)
    cross_stats = crosstalk_summary(cross)
    state_summary = state_validation_summary(model, validation.waveform, run_config, levels, target_operator)

    summary = {
        "case_id": request.case_id,
        "family": request.family,
        "model_variant": request.model_variant,
        "include_chi_prime": bool(request.include_chi_prime),
        "random_seed": request.random_seed,
        "n_active": int(request.n_active),
        "logical_levels": [int(level) for level in levels],
        "n_cav": int(model.n_cav),
        "n_tr": int(model.n_tr),
        "pulse_duration_s": float(duration_s),
        "pulse_duration_ns": float(duration_s * 1.0e9),
        "chi_t_over_2pi": float(request.duration_chi_t),
        "parameter_count": int(len(levels) * len(optimization.parameters) + 1),
        "optimized_parameter_count": int(len(optimization.active_levels) * len(optimization.parameters)),
        "runtime_s": float(runtime_s),
        "min_transition_spacing_hz": float(min_transition_spacing_hz(model, levels, frame)),
        "transition_frequencies_hz": manifold_transition_frequencies_hz(model, levels, frame).tolist(),
        "restricted_process_fidelity": float(metrics["process_fidelity"]),
        "average_gate_fidelity": float(metrics["average_gate_fidelity"]),
        "frobenius_error": float(metrics["frobenius_error"]),
        "operator_norm_error": float(metrics["operator_norm_error"]),
        "restricted_unitarity_error": float(metrics["restricted_unitarity_error"]),
        "state_transfer_fidelity_mean": float(validation.state_transfer_fidelity_mean),
        "state_transfer_fidelity_min": float(validation.state_transfer_fidelity_min),
        "same_block_population_mean": float(validation.same_block_population_mean),
        "same_block_population_min": float(validation.same_block_population_min),
        "other_target_population_mean": float(validation.other_target_population_mean),
        "other_target_population_max": float(validation.other_target_population_max),
        "leakage_outside_target_mean": float(validation.leakage_outside_target_mean),
        "leakage_outside_target_max": float(validation.leakage_outside_target_max),
        "weighted_loss": float(validation.weighted_loss),
        "qubit_loss": float(validation.qubit_loss),
        "subspace_loss": float(validation.subspace_loss),
        "preservation_loss": float(validation.preservation_loss),
        "leakage_loss": float(validation.leakage_loss),
        "crosstalk_diagonal_mean": float(cross_stats["diagonal_mean"]),
        "crosstalk_diagonal_min": float(cross_stats["diagonal_min"]),
        "crosstalk_offdiag_mean": float(cross_stats["offdiag_mean"]),
        "crosstalk_offdiag_max": float(cross_stats["offdiag_max"]),
        "per_block_process_fidelities": [float(item["process_fidelity"]) for item in per_block],
        "per_block_average_gate_fidelities": [float(item["average_gate_fidelity"]) for item in per_block],
        "per_block_rotation_angle_errors_rad": [float(item["rotation_angle_error_rad"]) for item in per_block],
        "per_block_rotation_axis_errors_rad": [float(item["rotation_axis_error_rad"]) for item in per_block],
        "per_block_residual_z_errors_rad": [float(item["residual_z_error_rad"]) for item in per_block],
        "per_block_transverse_errors_rad": [float(item["transverse_error_rad"]) for item in per_block],
        "mean_residual_z_error_rad": float(np.mean(residual_z)),
        "max_residual_z_error_rad": float(np.max(residual_z)),
        "mean_transverse_error_rad": float(np.mean(transverse)),
        "max_transverse_error_rad": float(np.max(transverse)),
        "state_validation_ground_fidelity": float(state_summary["states"][0]["state_fidelity"]),
        "state_validation_plus_fidelity": float(state_summary["states"][1]["state_fidelity"]),
        "success_tier": "",
        "notable_failure_mode": "",
    }
    summary["success_tier"] = classify_success(summary["average_gate_fidelity"])
    summary["notable_failure_mode"] = classify_failure_mode(summary)

    waveform_payload = sample_waveform_payload(validation.waveform, duration_s)
    artifact_payload = {
        "study_name": STUDY_DIR.name,
        "date_created": time.strftime("%Y-%m-%d"),
        "description": "Optimized multitone SQR waveform for an arbitrary block-diagonal active-subspace target.",
        "case_request": json_ready(request.__dict__),
        "target_metadata": target_meta,
        "target_operator": target_operator,
        "restricted_operator": validation.restricted_operator,
        "block_metrics": per_block,
        "active_subspace_metrics": metrics,
        "crosstalk_matrix": cross,
        "crosstalk_summary": cross_stats,
        "state_validation": state_summary,
        "optimization_summary": optimization.improvement_summary(),
        "optimization_history": optimization.history,
        "objective_weights": OBJECTIVE_WEIGHTS.as_dict(),
        "optimized_corrections": {
            "d_lambda": list(optimization.optimized_corrections.d_lambda),
            "d_alpha": list(optimization.optimized_corrections.d_alpha),
            "d_omega_rad_s": list(optimization.optimized_corrections.d_omega_rad_s),
            "d_omega_hz": [float(value / (2.0 * np.pi)) for value in optimization.optimized_corrections.d_omega_rad_s],
        },
        "tone_specs": validation.waveform.tone_rows(),
        "waveform_samples": waveform_payload,
        "summary_row": summary,
        "load_instructions": "Load this JSON and inspect `tone_specs`, `optimized_corrections`, `restricted_operator`, and `crosstalk_matrix` to reproduce the reported case.",
    }
    save_json(ARTIFACTS_DIR / "cases" / f"{request.case_id}.json", artifact_payload)
    return summary


def _group_best(df: pd.DataFrame, by: list[str]) -> pd.DataFrame:
    ordered = df.sort_values("average_gate_fidelity", ascending=False)
    return ordered.groupby(by, as_index=False).first()


def build_machine_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    df = pd.DataFrame(rows)
    best_overall = df.sort_values("average_gate_fidelity", ascending=False).iloc[0].to_dict()
    worst_overall = df.sort_values("average_gate_fidelity", ascending=True).iloc[0].to_dict()
    family_summary = _group_best(df, ["model_variant", "family", "n_active"])
    coherent_error_summary = (
        df.groupby(["model_variant", "family", "n_active"], as_index=False)
        .agg(
            avg_fidelity_mean=("average_gate_fidelity", "mean"),
            residual_z_mean=("mean_residual_z_error_rad", "mean"),
            residual_z_best=("mean_residual_z_error_rad", "min"),
            transverse_mean=("mean_transverse_error_rad", "mean"),
            transverse_best=("mean_transverse_error_rad", "min"),
        )
        .sort_values(["model_variant", "family", "n_active"])
    )
    random_df = df[df["family"] == "D"].copy()
    random_stats: list[dict[str, Any]] = []
    for (model_variant, n_active, chi_t), group in random_df.groupby(["model_variant", "n_active", "chi_t_over_2pi"]):
        values = group["average_gate_fidelity"].astype(float).to_numpy()
        random_stats.append(
            {
                "model_variant": model_variant,
                "n_active": int(n_active),
                "chi_t_over_2pi": float(chi_t),
                "count": int(values.size),
                "median": float(np.median(values)),
                "best": float(np.max(values)),
                "worst": float(np.min(values)),
                "mean": float(np.mean(values)),
                "std": float(np.std(values)),
            }
        )
    return {
        "study": STUDY_DIR.name,
        "success_thresholds": SUCCESS_TIERS,
        "n_cases": int(len(rows)),
        "best_overall": best_overall,
        "worst_overall": worst_overall,
        "lowest_residual_z_case": df.sort_values("mean_residual_z_error_rad", ascending=True).iloc[0].to_dict(),
        "lowest_transverse_case": df.sort_values("mean_transverse_error_rad", ascending=True).iloc[0].to_dict(),
        "case_rows": rows,
        "best_by_family_model_nactive": family_summary.to_dict(orient="records"),
        "coherent_error_summary": coherent_error_summary.to_dict(orient="records"),
        "random_ensemble_statistics": random_stats,
    }


def save_results(rows: list[dict[str, Any]]) -> None:
    payload = build_machine_summary(rows)
    save_json(RESULTS_PATH, payload)
    save_json(MACHINE_SUMMARY_PATH, payload)
    save_json(CHECKPOINT_PATH, payload)


def _complex_array_from_json(value: Any) -> np.ndarray:
    if isinstance(value, dict) and {"real", "imag", "shape"}.issubset(value):
        real = np.asarray(value["real"], dtype=float)
        imag = np.asarray(value["imag"], dtype=float)
        shape = tuple(int(item) for item in value["shape"])
        return (real + 1.0j * imag).reshape(shape)
    return np.asarray(value, dtype=np.complex128)


def _enrich_summary_row(summary_row: dict[str, Any], artifact_payload: dict[str, Any]) -> dict[str, Any]:
    row = dict(summary_row)
    block_metrics = artifact_payload.get("block_metrics", [])
    residual_z: np.ndarray | None = None
    transverse: np.ndarray | None = None
    if isinstance(block_metrics, list) and block_metrics and "residual_z_error_rad" in block_metrics[0]:
        residual_z = np.asarray([float(item["residual_z_error_rad"]) for item in block_metrics], dtype=float)
        transverse = np.asarray([float(item["transverse_error_rad"]) for item in block_metrics], dtype=float)
    else:
        target_operator = artifact_payload.get("target_operator")
        restricted_operator = artifact_payload.get("restricted_operator")
        if target_operator is None or restricted_operator is None:
            return row
        target_blocks = restricted_blocks(_complex_array_from_json(target_operator))
        actual_blocks = restricted_blocks(_complex_array_from_json(restricted_operator))
        rebuilt = [
            block_rotation_metrics(target_block, actual_block)
            for target_block, actual_block in zip(target_blocks, actual_blocks, strict=True)
        ]
        residual_z = np.asarray([float(item["residual_z_error_rad"]) for item in rebuilt], dtype=float)
        transverse = np.asarray([float(item["transverse_error_rad"]) for item in rebuilt], dtype=float)
    row.setdefault("per_block_residual_z_errors_rad", residual_z.tolist())
    row.setdefault("per_block_transverse_errors_rad", transverse.tolist())
    row.setdefault("mean_residual_z_error_rad", float(np.mean(residual_z)))
    row.setdefault("max_residual_z_error_rad", float(np.max(residual_z)))
    row.setdefault("mean_transverse_error_rad", float(np.mean(transverse)))
    row.setdefault("max_transverse_error_rad", float(np.max(transverse)))
    return row


def existing_rows() -> list[dict[str, Any]]:
    rows_by_id: dict[str, dict[str, Any]] = {}
    for path in (CHECKPOINT_PATH, RESULTS_PATH, MACHINE_SUMMARY_PATH):
        if not path.exists():
            continue
        try:
            payload = load_json(path)
        except Exception:
            continue
        for row in payload.get("case_rows", []):
            rows_by_id[str(row["case_id"])] = row
    cases_dir = ARTIFACTS_DIR / "cases"
    if cases_dir.exists():
        for artifact_path in sorted(cases_dir.glob("*.json")):
            try:
                payload = load_json(artifact_path)
            except Exception:
                continue
            summary_row = payload.get("summary_row")
            if isinstance(summary_row, dict) and "case_id" in summary_row:
                rows_by_id[str(summary_row["case_id"])] = _enrich_summary_row(summary_row, payload)
    return list(rows_by_id.values())


def load_case_artifact(case_id: str) -> dict[str, Any]:
    return load_json(ARTIFACTS_DIR / "cases" / f"{case_id}.json")


def save_figure(fig, stem: str) -> None:
    fig.savefig(FIGURES_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_transition_frequency_diagram() -> None:
    apply_plot_style()
    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    levels = tuple(range(8))
    for label, include_chi_prime in (("chi only", False), ("chi + chi'", True)):
        model = build_model(include_chi_prime=include_chi_prime, n_active=8)
        frame = build_frame(model)
        freqs_mhz = manifold_transition_frequencies_hz(model, levels, frame) / 1.0e6
        ax.plot(levels, freqs_mhz, marker="o", label=label)
    ax.set_xlabel("Fock level n")
    ax.set_ylabel(r"$\omega_q^{(n)} / 2\pi$ (MHz)")
    ax.set_title("Fock-Conditioned Qubit Transition Frequencies")
    ax.legend()
    save_figure(fig, "transition_frequency_diagram")


def _select_waveform_cases(df: pd.DataFrame) -> tuple[str, str]:
    structured = df[df["family"] != "D"].sort_values("average_gate_fidelity", ascending=False).iloc[0]["case_id"]
    difficult_random = (
        df[(df["family"] == "D") & (df["model_variant"] == "chi_plus_chiprime")]
        .sort_values("average_gate_fidelity", ascending=True)
        .iloc[0]["case_id"]
    )
    return str(structured), str(difficult_random)


def plot_representative_waveforms(df: pd.DataFrame) -> None:
    apply_plot_style()
    structured_id, difficult_id = _select_waveform_cases(df)
    artifacts = [("Best structured", load_case_artifact(structured_id)), ("Hard random", load_case_artifact(difficult_id))]
    fig, axes = plt.subplots(2, 2, figsize=(10.0, 5.5), sharex="col")
    for column, (title, artifact) in enumerate(artifacts):
        waveform = artifact["waveform_samples"]
        time_ns = 1.0e9 * np.asarray(waveform["time_s"], dtype=float)
        real = np.asarray(waveform["real"], dtype=float)
        imag = np.asarray(waveform["imag"], dtype=float)
        axes[0, column].plot(time_ns, real, label="I")
        axes[0, column].plot(time_ns, imag, label="Q")
        axes[0, column].set_title(title)
        axes[0, column].set_ylabel("Drive amplitude (rad/s)")
        axes[0, column].legend()
        sample_dt = float(np.mean(np.diff(np.asarray(waveform["time_s"], dtype=float))))
        spectrum = np.fft.fftshift(np.fft.fft(real + 1.0j * imag))
        freqs = np.fft.fftshift(np.fft.fftfreq(real.size, d=sample_dt)) / 1.0e6
        axes[1, column].plot(freqs, np.abs(spectrum))
        axes[1, column].set_xlabel("Frequency offset (MHz)")
        axes[1, column].set_ylabel("|FFT| (a.u.)")
    fig.suptitle("Representative Multitone Waveforms and Spectra")
    save_figure(fig, "representative_multitone_waveforms")


def plot_blockwise_heatmap(df: pd.DataFrame) -> None:
    apply_plot_style()
    subset = df[(df["model_variant"] == "chi_plus_chiprime") & (df["n_active"] == df["n_active"].max()) & (df["chi_t_over_2pi"] == 5.0)]
    rows = []
    labels = []
    for family in ("A", "B", "C"):
        case = subset[subset["family"] == family].sort_values("average_gate_fidelity", ascending=False).iloc[0]
        rows.append(case["per_block_average_gate_fidelities"])
        labels.append(f"Family {family}")
    random_case = subset[subset["family"] == "D"].sort_values("average_gate_fidelity", ascending=False).iloc[0]
    rows.append(random_case["per_block_average_gate_fidelities"])
    labels.append("Family D")
    matrix = np.asarray(rows, dtype=float)
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    image = ax.imshow(matrix, aspect="auto", vmin=0.0, vmax=1.0, cmap="viridis")
    ax.set_xticks(range(matrix.shape[1]), labels=[str(index) for index in range(matrix.shape[1])])
    ax.set_yticks(range(matrix.shape[0]), labels=labels)
    ax.set_xlabel("Fock manifold n")
    ax.set_ylabel("Target family")
    ax.set_title("Per-Block Average Gate Fidelity")
    fig.colorbar(image, ax=ax, label="Fidelity")
    save_figure(fig, "blockwise_fidelity_heatmap")


def plot_crosstalk_heatmap(df: pd.DataFrame) -> None:
    apply_plot_style()
    case = (
        df[(df["family"] == "D") & (df["model_variant"] == "chi_plus_chiprime") & (df["n_active"] == df["n_active"].max())]
        .sort_values("average_gate_fidelity", ascending=True)
        .iloc[0]
    )
    artifact = load_case_artifact(str(case["case_id"]))
    matrix = np.asarray(artifact["crosstalk_matrix"], dtype=float)
    fig, ax = plt.subplots(figsize=(5.2, 4.2))
    image = ax.imshow(matrix, aspect="auto", cmap="magma")
    ax.set_xlabel("Intended tone manifold")
    ax.set_ylabel("Responding manifold")
    ax.set_title("Single-Tone Crosstalk Matrix")
    fig.colorbar(image, ax=ax, label=r"$P(|g,m\rangle \to |e,m\rangle)$")
    save_figure(fig, "crosstalk_heatmap")


def plot_fidelity_vs_duration(df: pd.DataFrame) -> None:
    apply_plot_style()
    target_n = 4 if 4 in set(df["n_active"]) else int(df["n_active"].max())
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0), sharey=True)
    for axis, (model_variant, title) in zip(axes, (("chi_only", "chi only"), ("chi_plus_chiprime", "chi + chi'")), strict=True):
        subset = df[df["model_variant"] == model_variant]
        for family in ("A", "B", "C"):
            family_df = subset[(subset["family"] == family) & (subset["n_active"] == target_n)]
            grouped = family_df.groupby("chi_t_over_2pi", as_index=False)["average_gate_fidelity"].max()
            axis.plot(grouped["chi_t_over_2pi"], grouped["average_gate_fidelity"], marker="o", label=f"Family {family}")
        random_df = subset[(subset["family"] == "D") & (subset["n_active"] == target_n)]
        grouped = random_df.groupby("chi_t_over_2pi", as_index=False)["average_gate_fidelity"].median()
        axis.plot(grouped["chi_t_over_2pi"], grouped["average_gate_fidelity"], marker="s", linestyle="--", label="Family D median")
        axis.set_title(title)
        axis.set_xlabel(r"$|\chi| T / 2\pi$")
    axes[0].set_ylabel("Average gate fidelity")
    axes[0].legend()
    save_figure(fig, "fidelity_vs_pulse_duration")


def plot_fidelity_vs_nactive(df: pd.DataFrame) -> None:
    apply_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0), sharey=True)
    duration_value = 5.0
    for axis, model_variant in zip(axes, ("chi_only", "chi_plus_chiprime"), strict=True):
        subset = df[(df["model_variant"] == model_variant) & (df["chi_t_over_2pi"] == duration_value)]
        for family in ("A", "B", "C"):
            family_df = subset[subset["family"] == family].groupby("n_active", as_index=False)["average_gate_fidelity"].max()
            axis.plot(family_df["n_active"], family_df["average_gate_fidelity"], marker="o", label=f"Family {family}")
        random_df = subset[subset["family"] == "D"].groupby("n_active", as_index=False)["average_gate_fidelity"].median()
        axis.plot(random_df["n_active"], random_df["average_gate_fidelity"], marker="s", linestyle="--", label="Family D median")
        axis.set_title(model_variant.replace("_", " "))
        axis.set_xlabel(r"$N_{\mathrm{active}}$")
    axes[0].set_ylabel("Average gate fidelity")
    axes[0].legend()
    save_figure(fig, "fidelity_vs_active_subspace_size")


def plot_random_histogram(df: pd.DataFrame) -> None:
    apply_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0), sharey=True)
    n_active = 4 if 4 in set(df["n_active"]) else int(df["n_active"].max())
    duration_value = 3.0
    for axis, model_variant in zip(axes, ("chi_only", "chi_plus_chiprime"), strict=True):
        values = df[
            (df["family"] == "D")
            & (df["model_variant"] == model_variant)
            & (df["n_active"] == n_active)
            & (df["chi_t_over_2pi"] == duration_value)
        ]["average_gate_fidelity"].astype(float)
        axis.hist(values, bins=min(8, max(4, len(values))), alpha=0.8)
        axis.axvline(float(np.median(values)), color="black", linestyle="--", linewidth=1.0, label="median")
        axis.set_title(model_variant.replace("_", " "))
        axis.set_xlabel("Average gate fidelity")
        axis.legend()
    axes[0].set_ylabel("Count")
    save_figure(fig, "random_target_fidelity_histogram")


def plot_state_validation(df: pd.DataFrame) -> None:
    apply_plot_style()
    best_case_id = str(df.sort_values("average_gate_fidelity", ascending=False).iloc[0]["case_id"])
    worst_random_id = str(df[df["family"] == "D"].sort_values("average_gate_fidelity", ascending=True).iloc[0]["case_id"])
    artifacts = [("Best case", load_case_artifact(best_case_id)), ("Worst random", load_case_artifact(worst_random_id))]
    fig, axes = plt.subplots(2, 2, figsize=(10.0, 6.0))
    for row, (title, artifact) in enumerate(artifacts):
        states = artifact["state_validation"]["states"]
        labels = [state["label"] for state in states]
        fidelities = [state["state_fidelity"] for state in states]
        axes[row, 0].bar(labels, fidelities)
        axes[row, 0].set_ylim(0.0, 1.0)
        axes[row, 0].set_ylabel("State fidelity")
        axes[row, 0].set_title(title)
        populations = np.asarray(states[1]["cavity_level_populations"], dtype=float)
        axes[row, 1].bar(np.arange(populations.size), populations)
        axes[row, 1].set_xlabel("Fock level n")
        axes[row, 1].set_ylabel("Population")
        axes[row, 1].set_title(f"{title} conditioned populations")
    save_figure(fig, "state_level_validation")


def plot_best_worst_tables(df: pd.DataFrame) -> None:
    apply_plot_style()
    best = df.sort_values("average_gate_fidelity", ascending=False).head(3)
    worst = df.sort_values("average_gate_fidelity", ascending=True).head(3)
    fig, axes = plt.subplots(2, 1, figsize=(10.0, 4.8))
    for axis, table_df, title in zip(axes, (best, worst), ("Best-performing cases", "Worst-performing cases"), strict=True):
        axis.axis("off")
        table = table_df[["case_id", "family", "model_variant", "n_active", "chi_t_over_2pi", "average_gate_fidelity", "crosstalk_offdiag_max", "notable_failure_mode"]]
        rendered = axis.table(cellText=table.values, colLabels=table.columns, loc="center")
        rendered.auto_set_font_size(False)
        rendered.set_fontsize(7.0)
        rendered.scale(1.0, 1.25)
        axis.set_title(title)
    save_figure(fig, "best_worst_parameter_tables")


def plot_chi_comparison(df: pd.DataFrame) -> None:
    apply_plot_style()
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    compare = (
        df[df["family"] == "D"]
        .groupby(["model_variant", "n_active"], as_index=False)["average_gate_fidelity"]
        .median()
    )
    for model_variant, group in compare.groupby("model_variant"):
        ax.plot(group["n_active"], group["average_gate_fidelity"], marker="o", label=model_variant.replace("_", " "))
    ax.set_xlabel(r"$N_{\mathrm{active}}$")
    ax.set_ylabel("Median random-target fidelity")
    ax.set_title("chi-only versus chi + chi' Comparison")
    ax.legend()
    save_figure(fig, "chi_only_vs_chi_plus_chiprime")


def plot_coherent_error_decomposition(df: pd.DataFrame) -> None:
    apply_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.2), sharey=True)
    palette = {"A": "#4477AA", "B": "#228833", "C": "#CCBB44", "D": "#EE6677"}
    for axis, family_kind, title in zip(axes, ("structured", "random"), ("Structured targets", "Random targets"), strict=True):
        subset = df[df["family"] != "D"] if family_kind == "structured" else df[df["family"] == "D"]
        for family, group in subset.groupby("family"):
            axis.scatter(
                group["mean_residual_z_error_rad"],
                group["mean_transverse_error_rad"],
                s=30,
                alpha=0.7,
                color=palette.get(str(family), "#333333"),
                label=f"Family {family}",
            )
        axis.set_title(title)
        axis.set_xlabel("Mean residual Z error (rad)")
    axes[0].set_ylabel("Mean transverse error (rad)")
    axes[0].legend(frameon=False)
    save_figure(fig, "coherent_error_decomposition")


def generate_figures(rows: list[dict[str, Any]]) -> None:
    df = pd.DataFrame(rows)
    plot_transition_frequency_diagram()
    plot_representative_waveforms(df)
    plot_blockwise_heatmap(df)
    plot_crosstalk_heatmap(df)
    plot_fidelity_vs_duration(df)
    plot_fidelity_vs_nactive(df)
    plot_random_histogram(df)
    plot_state_validation(df)
    plot_best_worst_tables(df)
    plot_chi_comparison(df)
    plot_coherent_error_decomposition(df)


def main() -> None:
    requests = case_requests()
    row_by_id = {str(row["case_id"]): row for row in existing_rows()}
    total = len(requests)
    for index, request in enumerate(requests, start=1):
        if request.case_id in row_by_id:
            print(f"[{index:03d}/{total:03d}] {request.case_id} already present; skipping", flush=True)
            continue
        print(
            f"[{index:03d}/{total:03d}] {request.case_id} "
            f"family={request.family} model={request.model_variant} "
            f"N_active={request.n_active} chiT={request.duration_chi_t}",
            flush=True,
        )
        row = run_case(request)
        row_by_id[row["case_id"]] = row
        ordered_rows = [row_by_id[request_obj.case_id] for request_obj in requests if request_obj.case_id in row_by_id]
        save_results(ordered_rows)
    rows = [row_by_id[request.case_id] for request in requests if request.case_id in row_by_id]
    generate_figures(rows)
    save_results(rows)
    print(f"Saved {RESULTS_PATH}")


if __name__ == "__main__":
    main()