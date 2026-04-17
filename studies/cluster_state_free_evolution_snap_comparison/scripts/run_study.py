from __future__ import annotations

import copy
import time
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from .common import (
        ARTIFACT_DIR,
        DATA_DIR,
        FAMILY_SPECS,
        FIG_DIR,
        FRONTIER_WEIGHTS,
        REFINE_MAXITER,
        REFINE_SEEDS,
        SCREEN_MAXITER,
        SCREEN_SEED,
        SNAP,
        TRUNCATION_LEVELS,
        FreeEvolveCondPhase,
        ablation_sequences,
        base,
        build_family_sequence,
        case_id,
        evaluate_sequence_with_diagnostics,
        sequence_complexity_metrics,
        sequence_for_n_cav,
        sequence_from_payload,
        sequence_phase_budget,
    )
except ImportError:
    from common import (
        ARTIFACT_DIR,
        DATA_DIR,
        FAMILY_SPECS,
        FIG_DIR,
        FRONTIER_WEIGHTS,
        REFINE_MAXITER,
        REFINE_SEEDS,
        SCREEN_MAXITER,
        SCREEN_SEED,
        SNAP,
        TRUNCATION_LEVELS,
        FreeEvolveCondPhase,
        ablation_sequences,
        base,
        build_family_sequence,
        case_id,
        evaluate_sequence_with_diagnostics,
        sequence_complexity_metrics,
        sequence_for_n_cav,
        sequence_from_payload,
        sequence_phase_budget,
    )

STYLE_PATH = Path(__file__).resolve().parents[3] / ".github" / "skills" / "publication-figures" / "assets" / "cqed_style.mplstyle"
if STYLE_PATH.exists():
    plt.style.use(str(STYLE_PATH))


def _json_dump(path: Path, payload: Any) -> None:
    base.save_json(path, payload)


def _family_title(family: str) -> str:
    return str(FAMILY_SPECS[str(family)]["title"])


def _case_summary_dict(
    *,
    family: str,
    blocks: int,
    seed: int,
    maxiter: int,
    duration_weight: float,
    fit_payload: dict[str, Any],
) -> dict[str, Any]:
    sequence = fit_payload["result"].sequence
    evaluation = evaluate_sequence_with_diagnostics(sequence, n_cav=base.DECOMP_N_CAV)
    complexity = sequence_complexity_metrics(sequence)
    phase_budget = sequence_phase_budget(sequence)
    seq_summary = fit_payload["summary"]

    return {
        "case_id": case_id(family, blocks),
        "family": family,
        "family_title": _family_title(family),
        "blocks": int(blocks),
        "seed": int(seed),
        "maxiter": int(maxiter),
        "duration_weight": float(duration_weight),
        "fidelity": float(evaluation["fidelity"]),
        "block_gauge_fidelity": float(evaluation["block_gauge_fidelity"]),
        "best_fit_block_gauge_fidelity": float(evaluation["best_fit_block_gauge_fidelity"]),
        "leakage_average": float(evaluation["leakage_average"]),
        "leakage_worst": float(evaluation["leakage_worst"]),
        "unitarity_error": float(evaluation["unitarity_error"]),
        "block_phases_rad": evaluation["block_phases_rad"],
        "best_fit_correction_phases_rad": evaluation["best_fit_correction_phases_rad"],
        "rms_block_phase_error_rad": float(evaluation["rms_block_phase_error_rad"]),
        "objective": float(fit_payload["objective"]),
        "success": bool(fit_payload["success"]),
        "message": str(fit_payload["message"]),
        "metrics": fit_payload["metrics"],
        "gate_depth": int(seq_summary["gate_depth"]),
        "total_duration_ns": float(seq_summary["total_duration_ns"]),
        "total_wait_time_ns": float(phase_budget["total_wait_time_ns"]),
        "total_fe_logical_delta_phi_rad": float(phase_budget["total_fe_logical_delta_phi_rad"]),
        "total_snap_logical_relative_phase_rad": float(phase_budget["total_snap_logical_relative_phase_rad"]),
        "parameter_count": int(complexity["parameter_count"]),
        "snap_gate_count": int(complexity["snap_gate_count"]),
        "snap_phase_count": int(complexity["snap_phase_count"]),
        "entangling_gate_count": int(complexity["entangling_gate_count"]),
        "wait_gate_count": int(complexity["wait_gate_count"]),
        "sequence_summary": seq_summary,
        "sequence_payload": fit_payload["sequence_payload"],
        "phase_rows": phase_budget["phase_rows"],
        "phase_gate_rows": {
            "free_evolution": phase_budget["fe_rows"],
            "snap": phase_budget["snap_rows"],
        },
    }


def _run_fit(
    *,
    family: str,
    blocks: int,
    seed: int,
    maxiter: int,
    duration_weight: float = 0.0,
    warm_start: Any | None = None,
) -> tuple[dict[str, Any], Any]:
    sequence = build_family_sequence(family=family, blocks=blocks)
    t0 = time.perf_counter()
    fit_payload = base.fit_sequence(
        sequence,
        n_cav=base.DECOMP_N_CAV,
        seed=int(seed),
        init_guess="heuristic",
        multistart=1,
        maxiter=int(maxiter),
        duration_weight=float(duration_weight),
        warm_start=warm_start,
    )
    summary = _case_summary_dict(
        family=family,
        blocks=blocks,
        seed=seed,
        maxiter=maxiter,
        duration_weight=duration_weight,
        fit_payload=fit_payload,
    )
    summary["elapsed_s"] = float(time.perf_counter() - t0)
    return summary, fit_payload["result"]


def _screen_cases() -> dict[str, dict[str, Any]]:
    screen_results: dict[str, dict[str, Any]] = {}
    print("[screen] coarse search", flush=True)
    for family, spec in FAMILY_SPECS.items():
        for blocks in spec["blocks"]:
            print(f"  {case_id(family, blocks)}", flush=True)
            summary, _ = _run_fit(
                family=family,
                blocks=int(blocks),
                seed=SCREEN_SEED,
                maxiter=SCREEN_MAXITER,
            )
            screen_results[summary["case_id"]] = summary
            print(
                f"    F={summary['fidelity']:.4f}  depth={summary['gate_depth']}  "
                f"wait={summary['total_wait_time_ns']:.1f} ns",
                flush=True,
            )
    return screen_results


def _best_case_per_family(results: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    out: dict[str, dict[str, Any]] = {}
    for family in FAMILY_SPECS:
        family_rows = [row for row in results.values() if row["family"] == family]
        out[family] = max(
            family_rows,
            key=lambda row: (float(row["fidelity"]), -float(row["total_wait_time_ns"]), -int(row["gate_depth"])),
        )
    return out


def _refine_cases(screen_results: dict[str, dict[str, Any]]) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
    family_best = _best_case_per_family(screen_results)
    refined_results: dict[str, dict[str, Any]] = {}
    raw_best_results: dict[str, Any] = {}
    print("[refine] targeted search", flush=True)
    for family, coarse in family_best.items():
        best_summary: dict[str, Any] | None = None
        best_raw: Any | None = None
        for seed in REFINE_SEEDS:
            print(f"  {coarse['case_id']} seed={seed}", flush=True)
            summary, raw_result = _run_fit(
                family=family,
                blocks=int(coarse["blocks"]),
                seed=int(seed),
                maxiter=REFINE_MAXITER,
            )
            if best_summary is None or float(summary["fidelity"]) > float(best_summary["fidelity"]):
                best_summary = summary
                best_raw = raw_result
        assert best_summary is not None
        refined_results[family] = best_summary
        raw_best_results[family] = best_raw
        print(
            f"    best F={best_summary['fidelity']:.4f}  depth={best_summary['gate_depth']}  "
            f"wait={best_summary['total_wait_time_ns']:.1f} ns",
            flush=True,
        )
    return refined_results, raw_best_results


def _duration_frontier(
    *,
    base_case: dict[str, Any],
    warm_start: Any,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    print(f"[frontier] {base_case['case_id']}", flush=True)
    warm_payload = warm_start
    for weight in FRONTIER_WEIGHTS:
        summary, raw_result = _run_fit(
            family=str(base_case["family"]),
            blocks=int(base_case["blocks"]),
            seed=int(base_case["seed"]),
            maxiter=SCREEN_MAXITER,
            duration_weight=float(weight),
            warm_start=warm_payload,
        )
        summary["frontier_weight"] = float(weight)
        rows.append(summary)
        warm_payload = raw_result
        print(
            f"  w={weight:.4f}: F={summary['fidelity']:.4f}, "
            f"wait={summary['total_wait_time_ns']:.1f} ns, dur={summary['total_duration_ns']:.1f} ns",
            flush=True,
        )
    return rows


def _threshold_summary(rows: list[dict[str, Any]], *, family: str, thresholds: tuple[float, ...]) -> dict[str, Any]:
    family_rows = [row for row in rows if row["family"] == family]
    out: dict[str, Any] = {}
    for threshold in thresholds:
        feasible = [row for row in family_rows if float(row["fidelity"]) >= float(threshold)]
        if not feasible:
            out[f"{threshold:.3f}"] = None
            continue
        best_depth = min(feasible, key=lambda row: (int(row["gate_depth"]), float(row["total_wait_time_ns"])))
        best_wait = min(feasible, key=lambda row: (float(row["total_wait_time_ns"]), int(row["gate_depth"])))
        out[f"{threshold:.3f}"] = {
            "best_depth": {
                "case_id": best_depth["case_id"],
                "gate_depth": int(best_depth["gate_depth"]),
                "total_wait_time_ns": float(best_depth["total_wait_time_ns"]),
                "fidelity": float(best_depth["fidelity"]),
            },
            "best_wait": {
                "case_id": best_wait["case_id"],
                "gate_depth": int(best_wait["gate_depth"]),
                "total_wait_time_ns": float(best_wait["total_wait_time_ns"]),
                "fidelity": float(best_wait["fidelity"]),
            },
        }
    return out


def _run_ablations(best_case: dict[str, Any]) -> dict[str, Any]:
    sequence = sequence_from_payload(best_case["sequence_payload"], n_cav=base.DECOMP_N_CAV)
    variants = ablation_sequences(sequence)
    out: dict[str, Any] = {}
    for label, variant in variants.items():
        evaluation = evaluate_sequence_with_diagnostics(variant, n_cav=base.DECOMP_N_CAV)
        phase_budget = sequence_phase_budget(variant)
        out[label] = {
            "fidelity": float(evaluation["fidelity"]),
            "block_gauge_fidelity": float(evaluation["block_gauge_fidelity"]),
            "gate_depth": int(len(variant.gates)),
            "total_wait_time_ns": float(phase_budget["total_wait_time_ns"]),
            "snap_gate_count": int(sum(isinstance(g, SNAP) for g in variant.gates)),
            "fe_gate_count": int(sum(isinstance(g, FreeEvolveCondPhase) for g in variant.gates)),
        }
    return out


def _run_truncation_validation(best_case: dict[str, Any]) -> list[dict[str, Any]]:
    base_sequence = sequence_from_payload(best_case["sequence_payload"], n_cav=base.DECOMP_N_CAV)
    rows: list[dict[str, Any]] = []
    for n_cav in TRUNCATION_LEVELS:
        seq = sequence_for_n_cav(base_sequence, n_cav=int(n_cav))
        evaluation = evaluate_sequence_with_diagnostics(seq, n_cav=int(n_cav))
        rows.append(
            {
                "n_cav": int(n_cav),
                "fidelity": float(evaluation["fidelity"]),
                "block_gauge_fidelity": float(evaluation["block_gauge_fidelity"]),
                "leakage_average": float(evaluation["leakage_average"]),
                "leakage_worst": float(evaluation["leakage_worst"]),
            }
        )
    return rows


def _save_best_sequence_artifact(best_case: dict[str, Any], filename: str) -> None:
    _json_dump(ARTIFACT_DIR / filename, best_case)


def _plot_fidelity_vs_depth(screen_results: dict[str, dict[str, Any]], refined_results: dict[str, dict[str, Any]]) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    for family, spec in FAMILY_SPECS.items():
        rows = sorted(
            [row for row in screen_results.values() if row["family"] == family],
            key=lambda row: int(row["gate_depth"]),
        )
        x = [int(row["gate_depth"]) for row in rows]
        y = [float(row["fidelity"]) for row in rows]
        ax.plot(x, y, "o--", color=spec["color"], alpha=0.7, label=f"{spec['title']} (screen)")
        refined = refined_results.get(family)
        if refined is not None:
            ax.plot(
                [int(refined["gate_depth"])],
                [float(refined["fidelity"])],
                marker="*",
                markersize=14,
                color=spec["color"],
            )
    ax.set_xlabel("Gate Depth")
    ax.set_ylabel("Logical Fidelity")
    ax.set_title("Cluster-Unitary Fidelity vs Gate Depth")
    ax.set_ylim(0.55, 1.01)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, ncol=1, loc="lower right")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"fidelity_vs_depth.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def _plot_wait_vs_fidelity(screen_results: dict[str, dict[str, Any]], refined_results: dict[str, dict[str, Any]]) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 5.2))
    for family, spec in FAMILY_SPECS.items():
        rows = [row for row in screen_results.values() if row["family"] == family]
        ax.scatter(
            [float(row["total_wait_time_ns"]) for row in rows],
            [float(row["fidelity"]) for row in rows],
            s=[18 + 4 * int(row["snap_gate_count"]) for row in rows],
            color=spec["color"],
            alpha=0.55,
            label=f"{spec['title']} (screen)",
        )
        refined = refined_results.get(family)
        if refined is not None:
            ax.scatter(
                [float(refined["total_wait_time_ns"])],
                [float(refined["fidelity"])],
                s=180,
                marker="*",
                color=spec["color"],
                edgecolor="black",
                linewidth=0.6,
            )
    ax.set_xlabel("Total Free-Evolution Wait Time (ns)")
    ax.set_ylabel("Logical Fidelity")
    ax.set_title("Fidelity vs Entangling Wait Budget")
    ax.set_ylim(0.55, 1.01)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"wait_vs_fidelity.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def _plot_ablations(ablations: dict[str, dict[str, Any]]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.3), sharey=True)
    for ax, (title, payload) in zip(axes, ablations.items(), strict=True):
        labels = list(payload.keys())
        fids = [float(payload[label]["fidelity"]) for label in labels]
        ax.bar(labels, fids, color=["#4477AA", "#EE6677", "#CCBB44"], edgecolor="black", linewidth=0.5)
        ax.set_title(title)
        ax.set_ylim(0.0, 1.01)
        ax.grid(True, alpha=0.25, axis="y")
        ax.tick_params(axis="x", rotation=18)
    axes[0].set_ylabel("Logical Fidelity")
    fig.suptitle("Ablation of Best Native and Best SNAP-Extended Sequences")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"sequence_ablations.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def _plot_phase_budgets(best_native: dict[str, Any], best_snap: dict[str, Any]) -> None:
    labels = ["Best native FE", "Best FE+SNAP"]
    fe_vals = [
        float(best_native["total_fe_logical_delta_phi_rad"]) / np.pi,
        float(best_snap["total_fe_logical_delta_phi_rad"]) / np.pi,
    ]
    snap_vals = [
        float(best_native["total_snap_logical_relative_phase_rad"]) / np.pi,
        float(best_snap["total_snap_logical_relative_phase_rad"]) / np.pi,
    ]
    x = np.arange(len(labels))
    width = 0.36
    fig, ax = plt.subplots(figsize=(7.4, 4.3))
    ax.bar(x - width / 2, fe_vals, width, label="Free-evolution conditional phase budget", color="#4477AA")
    ax.bar(x + width / 2, snap_vals, width, label="SNAP logical phase budget", color="#228833")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Accumulated logical phase budget (units of pi)")
    ax.set_title("Where the Phase Freedom Comes From")
    ax.grid(True, alpha=0.25, axis="y")
    ax.legend(fontsize=8)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"phase_budget_comparison.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def _plot_truncation_validation(validation: dict[str, list[dict[str, Any]]]) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    colors = {"best_native": "#4477AA", "best_snap": "#228833"}
    labels = {"best_native": "Best native FE", "best_snap": "Best FE+SNAP"}
    for key, rows in validation.items():
        ax.plot(
            [int(row["n_cav"]) for row in rows],
            [float(row["fidelity"]) for row in rows],
            "o-",
            color=colors[key],
            label=labels[key],
        )
    ax.set_xlabel("Cavity Truncation")
    ax.set_ylabel("Logical Fidelity")
    ax.set_title("Truncation Validation")
    ax.set_ylim(0.0, 1.01)
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"truncation_validation.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def _make_summary_table(
    refined_results: dict[str, dict[str, Any]],
    threshold_results: dict[str, Any],
    frontiers: dict[str, list[dict[str, Any]]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for family, row in refined_results.items():
        rows.append(
            {
                "family": family,
                "title": row["family_title"],
                "blocks": int(row["blocks"]),
                "fidelity": float(row["fidelity"]),
                "gate_depth": int(row["gate_depth"]),
                "total_duration_ns": float(row["total_duration_ns"]),
                "total_wait_time_ns": float(row["total_wait_time_ns"]),
                "snap_gate_count": int(row["snap_gate_count"]),
                "snap_phase_count": int(row["snap_phase_count"]),
                "parameter_count": int(row["parameter_count"]),
                "thresholds": threshold_results.get(family, {}),
                "frontier_rows": frontiers.get(family, []),
            }
        )
    return rows


def main() -> None:
    screen_results = _screen_cases()
    _json_dump(DATA_DIR / "screen_results.json", screen_results)

    refined_results, raw_best = _refine_cases(screen_results)
    _json_dump(DATA_DIR / "refined_results.json", refined_results)

    best_native = refined_results["native_fe"]
    best_snap_family = max(
        [family for family in refined_results if family != "native_fe"],
        key=lambda family: float(refined_results[family]["fidelity"]),
    )
    best_snap = refined_results[best_snap_family]

    frontiers = {
        "native_fe": _duration_frontier(base_case=best_native, warm_start=raw_best["native_fe"]),
        best_snap_family: _duration_frontier(base_case=best_snap, warm_start=raw_best[best_snap_family]),
    }
    _json_dump(DATA_DIR / "duration_frontiers.json", frontiers)

    thresholds = (0.90, 0.95, 0.97)
    all_rows = list(screen_results.values()) + list(refined_results.values())
    threshold_results = {
        family: _threshold_summary(all_rows, family=family, thresholds=thresholds)
        for family in refined_results
    }
    _json_dump(DATA_DIR / "threshold_summary.json", threshold_results)

    ablations = {
        "best_native": _run_ablations(best_native),
        "best_snap": _run_ablations(best_snap),
    }
    _json_dump(DATA_DIR / "ablation_results.json", ablations)

    validation = {
        "best_native": _run_truncation_validation(best_native),
        "best_snap": _run_truncation_validation(best_snap),
    }
    _json_dump(DATA_DIR / "validation_summary.json", validation)

    _save_best_sequence_artifact(best_native, "best_native_sequence.json")
    _save_best_sequence_artifact(best_snap, "best_snap_sequence.json")

    comparison_summary = {
        "best_native_family": "native_fe",
        "best_snap_family": best_snap_family,
        "best_native": best_native,
        "best_snap": best_snap,
        "summary_rows": _make_summary_table(refined_results, threshold_results, frontiers),
    }
    _json_dump(DATA_DIR / "study_summary.json", comparison_summary)

    _plot_fidelity_vs_depth(screen_results, refined_results)
    _plot_wait_vs_fidelity(screen_results, refined_results)
    _plot_ablations(ablations)
    _plot_phase_budgets(best_native, best_snap)
    _plot_truncation_validation(validation)

    print("\n[summary]", flush=True)
    print(
        f"  best native: {best_native['case_id']} F={best_native['fidelity']:.4f} "
        f"depth={best_native['gate_depth']} wait={best_native['total_wait_time_ns']:.1f} ns",
        flush=True,
    )
    print(
        f"  best SNAP:   {best_snap['case_id']} F={best_snap['fidelity']:.4f} "
        f"depth={best_snap['gate_depth']} wait={best_snap['total_wait_time_ns']:.1f} ns",
        flush=True,
    )
    print(f"  outputs: {DATA_DIR}", flush=True)


if __name__ == "__main__":
    main()
