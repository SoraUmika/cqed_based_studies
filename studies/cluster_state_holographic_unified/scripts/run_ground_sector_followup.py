from __future__ import annotations

import argparse
import csv
import sys
import time
from datetime import date
from pathlib import Path
from typing import Any, Iterable

import matplotlib.pyplot as plt
import numpy as np
import runtime_compat  # noqa: F401


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import common as c
import run_design_space_study as ds


STYLE_PATH = SCRIPT_DIR.parents[2] / ".github" / "skills" / "publication-figures" / "assets" / "cqed_style.mplstyle"
if STYLE_PATH.exists():
    plt.style.use(str(STYLE_PATH))


SCOPE_SUMMARY_PATH = c.DATA_DIR / "corrected_scope_summary.json"
TRUNCATION_VALUES = (10, 12, 14)
GROUND_THRESHOLD = 0.99
GROUND_REFINE_MAXITER = 16
GROUND_REFINE_SEED = 17
GROUND_OPTIMIZATION_N_CAV = 12
PHYSICAL_ANCILLA_LIMIT = 0.02
PHYSICAL_OUTSIDE_TARGET_LIMIT = 0.08
PHYSICAL_CONVERGENCE_LIMIT = 0.01

FAMILY_LABELS = {
    "drsqr": "D + R + SQR",
    "drcpsqr": "D + R + CPSQR",
}

FAMILY_ARTIFACTS = {
    "drsqr": c.ARTIFACT_DIR / "ground_sector_best_sqr.json",
    "drcpsqr": c.ARTIFACT_DIR / "ground_sector_best_cpsqr.json",
}


def _iter_record_dicts(node: Any) -> Iterable[dict[str, Any]]:
    if isinstance(node, dict):
        required = {"case_id", "family_key", "sequence", "parameter_vector", "time_vector"}
        if required.issubset(node.keys()):
            yield node
        for value in node.values():
            yield from _iter_record_dicts(value)
    elif isinstance(node, list):
        for item in node:
            yield from _iter_record_dicts(item)


def _existing_full_n12_fidelity(record: dict[str, Any]) -> float:
    replay = record.get("replay", {})
    if isinstance(replay, dict):
        n12 = replay.get("12")
        if isinstance(n12, dict) and "fidelity" in n12:
            return float(n12["fidelity"])
    return float(record.get("fidelity", float("-inf")))


def _record_sort_key(record: dict[str, Any]) -> tuple[float, float, float]:
    return (
        _existing_full_n12_fidelity(record),
        float(record.get("fidelity", float("-inf"))),
        -float(record.get("objective", float("inf"))),
    )


def load_scope_records() -> dict[str, dict[str, Any]]:
    payload = c.load_json(SCOPE_SUMMARY_PATH)
    best_by_case: dict[str, dict[str, Any]] = {}
    for record in _iter_record_dicts(payload):
        case_id = str(record["case_id"])
        current = best_by_case.get(case_id)
        candidate = dict(record)
        if current is None or _record_sort_key(candidate) > _record_sort_key(current):
            best_by_case[case_id] = candidate
    return best_by_case


def select_seed_records(records: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    best_by_family_block: dict[tuple[str, int], dict[str, Any]] = {}
    for record in records.values():
        family_key = str(record.get("family_key", ""))
        if family_key not in FAMILY_LABELS:
            continue
        key = (family_key, int(record["blocks"]))
        current = best_by_family_block.get(key)
        if current is None or _record_sort_key(record) > _record_sort_key(current):
            best_by_family_block[key] = record
    ordered = sorted(best_by_family_block.values(), key=lambda row: (str(row["family_key"]), int(row["blocks"])))
    return ordered


def _warm_start_from_record(record: dict[str, Any], *, n_cav: int) -> dict[str, Any]:
    sequence = ds.apply_solution_to_case(record, n_cav=int(n_cav))
    return {
        "parameter_vector": sequence.get_parameter_vector().tolist(),
        "time_raw_vector": sequence.get_time_raw_vector(active_only=True).tolist(),
    }


def evaluate_candidate_record(record: dict[str, Any], *, source: str, objective_variant: str) -> dict[str, Any]:
    evaluations: dict[str, Any] = {}
    for n_cav in TRUNCATION_VALUES:
        sequence = ds.apply_solution_to_case(record, n_cav=int(n_cav))
        full_eval = c.evaluate_sequence(sequence, n_cav=int(n_cav))
        ground_eval = c.evaluate_ground_sector_transfer(sequence, n_cav=int(n_cav))
        evaluations[str(int(n_cav))] = {
            "full_target_fidelity": float(full_eval["fidelity"]),
            "full_target_leakage_worst": float(full_eval["leakage_worst"]),
            "full_target_leakage_average": float(full_eval["leakage_average"]),
            "restricted_ground_fidelity": float(ground_eval["restricted_fidelity"]),
            "ground_subspace_leakage_worst": float(ground_eval["subspace_leakage_worst"]),
            "ground_subspace_leakage_average": float(ground_eval["subspace_leakage_average"]),
            "ancilla_excitation_worst": float(ground_eval["ancilla_excitation_worst"]),
            "ancilla_excitation_average": float(ground_eval["ancilla_excitation_average"]),
            "support_leakage_worst": float(ground_eval["support_leakage_worst"]),
            "support_leakage_average": float(ground_eval["support_leakage_average"]),
            "outside_target_leakage_worst": float(ground_eval["outside_target_leakage_worst"]),
            "outside_target_leakage_average": float(ground_eval["outside_target_leakage_average"]),
            "unitarity_error": float(ground_eval["unitarity_error"]),
        }

    n12 = evaluations["12"]
    ground_values = [float(evaluations[str(value)]["restricted_ground_fidelity"]) for value in TRUNCATION_VALUES]
    full_values = [float(evaluations[str(value)]["full_target_fidelity"]) for value in TRUNCATION_VALUES]
    ancilla_worst = max(float(evaluations[str(value)]["ancilla_excitation_worst"]) for value in TRUNCATION_VALUES)
    outside_worst = max(float(evaluations[str(value)]["outside_target_leakage_worst"]) for value in TRUNCATION_VALUES)
    support_worst = max(float(evaluations[str(value)]["support_leakage_worst"]) for value in TRUNCATION_VALUES)
    ground_convergence_delta = float(max(abs(ground_values[1] - ground_values[0]), abs(ground_values[2] - ground_values[1])))
    full_convergence_delta = float(max(abs(full_values[1] - full_values[0]), abs(full_values[2] - full_values[1])))
    physically_plausible = bool(
        ancilla_worst <= PHYSICAL_ANCILLA_LIMIT
        and outside_worst <= PHYSICAL_OUTSIDE_TARGET_LIMIT
        and ground_convergence_delta <= PHYSICAL_CONVERGENCE_LIMIT
    )
    plausibility_penalty = float(ancilla_worst + outside_worst + 5.0 * ground_convergence_delta)

    return {
        "candidate_id": f"{source}:{record['case_id']}",
        "source": str(source),
        "objective_variant": str(objective_variant),
        "case_id": str(record["case_id"]),
        "family_key": str(record["family_key"]),
        "family_label": str(record["family_label"]),
        "variant_key": str(record["variant_key"]),
        "variant_label": str(record["variant_label"]),
        "order_label": str(record.get("order_label", "")),
        "blocks": int(record["blocks"]),
        "max_tones": int(record["max_tones"]),
        "levels": [int(level) for level in record.get("levels", [])],
        "sequence": c.json_ready(record["sequence"]),
        "summary": {
            "n12_full_target_fidelity": float(n12["full_target_fidelity"]),
            "n12_ground_fidelity": float(n12["restricted_ground_fidelity"]),
            "n12_ground_gain": float(n12["restricted_ground_fidelity"] - n12["full_target_fidelity"]),
            "ground_convergence_delta": float(ground_convergence_delta),
            "full_convergence_delta": float(full_convergence_delta),
            "ancilla_excitation_worst": float(ancilla_worst),
            "support_leakage_worst": float(support_worst),
            "outside_target_leakage_worst": float(outside_worst),
            "physically_plausible": physically_plausible,
            "plausibility_penalty": plausibility_penalty,
        },
        "evaluations": evaluations,
    }


def refine_record_for_ground_sector(seed_record: dict[str, Any], *, maxiter: int, seed: int) -> dict[str, Any]:
    started = time.perf_counter()
    sequence = ds.apply_solution_to_case(seed_record, n_cav=GROUND_OPTIMIZATION_N_CAV)
    fit = c.fit_sequence(
        sequence,
        n_cav=GROUND_OPTIMIZATION_N_CAV,
        seed=int(seed),
        init_guess="heuristic",
        multistart=1,
        maxiter=int(maxiter),
        warm_start=_warm_start_from_record(seed_record, n_cav=GROUND_OPTIMIZATION_N_CAV),
        target_unitary=c.GROUND_SECTOR_TARGET_UNITARY,
        subspace=c.ground_sector_subspace(GROUND_OPTIMIZATION_N_CAV),
    )
    return {
        "case_id": str(seed_record["case_id"]),
        "family_key": str(seed_record["family_key"]),
        "family_label": str(seed_record["family_label"]),
        "variant_key": str(seed_record["variant_key"]),
        "variant_label": str(seed_record["variant_label"]),
        "builder_name": str(seed_record["builder_name"]),
        "builder_kwargs": dict(seed_record.get("builder_kwargs", {})),
        "order_tokens": [str(token) for token in seed_record.get("order_tokens", ())],
        "order_label": str(seed_record.get("order_label", "")),
        "levels": [int(level) for level in seed_record.get("levels", [])],
        "max_tones": int(seed_record["max_tones"]),
        "blocks": int(seed_record["blocks"]),
        "search_phase": "ground_sector_refine",
        "optimization_n_cav": int(GROUND_OPTIMIZATION_N_CAV),
        "seed": int(seed),
        "init_guess": "heuristic",
        "maxiter": int(maxiter),
        "multistart": 1,
        "fidelity": float(fit["fidelity"]),
        "objective": float(fit["objective"]),
        "success": bool(fit["success"]),
        "message": str(fit["message"]),
        "summary": dict(fit["summary"]),
        "metrics": dict(fit["metrics"]),
        "sequence": fit["sequence_payload"],
        "parameter_vector": fit["result"].sequence.get_parameter_vector().tolist(),
        "time_vector": fit["result"].sequence.get_time_vector(active_only=False).tolist(),
        "warm_start_payload": fit["result"].to_payload(include_history=False),
        "elapsed_s": float(time.perf_counter() - started),
        "base_case_id": str(seed_record["case_id"]),
        "base_full_target_n12_fidelity": float(_existing_full_n12_fidelity(seed_record)),
    }


def select_family_outcomes(candidates: list[dict[str, Any]]) -> dict[str, Any]:
    summary: dict[str, Any] = {}
    for family_key, family_label in FAMILY_LABELS.items():
        family_candidates = [candidate for candidate in candidates if candidate["family_key"] == family_key]
        family_candidates.sort(
            key=lambda row: (
                float(row["summary"]["n12_ground_fidelity"]),
                -float(row["summary"]["outside_target_leakage_worst"]),
                -float(row["summary"]["ground_convergence_delta"]),
            ),
            reverse=True,
        )
        best_absolute = family_candidates[0]
        threshold_candidates = [candidate for candidate in family_candidates if float(candidate["summary"]["n12_ground_fidelity"]) >= GROUND_THRESHOLD]
        threshold_candidates.sort(
            key=lambda row: (
                -int(row["blocks"]),
                float(row["summary"]["n12_ground_fidelity"]),
            ),
            reverse=True,
        )
        minimum_depth = None
        if threshold_candidates:
            minimum_depth = min(
                threshold_candidates,
                key=lambda row: (int(row["blocks"]), -float(row["summary"]["n12_ground_fidelity"])),
            )
        plausible_candidates = [candidate for candidate in family_candidates if bool(candidate["summary"]["physically_plausible"])]
        if plausible_candidates:
            most_plausible = max(
                plausible_candidates,
                key=lambda row: (
                    float(row["summary"]["n12_ground_fidelity"]),
                    -float(row["summary"]["outside_target_leakage_worst"]),
                ),
            )
        else:
            most_plausible = min(
                family_candidates,
                key=lambda row: (
                    float(row["summary"]["plausibility_penalty"]),
                    -float(row["summary"]["n12_ground_fidelity"]),
                ),
            )
        summary[family_key] = {
            "family_label": family_label,
            "best_absolute": compact_candidate(best_absolute),
            "minimum_depth_above_threshold": None if minimum_depth is None else compact_candidate(minimum_depth),
            "most_physically_plausible": compact_candidate(most_plausible),
        }
    return summary


def compact_candidate(candidate: dict[str, Any]) -> dict[str, Any]:
    return {
        "candidate_id": str(candidate["candidate_id"]),
        "source": str(candidate["source"]),
        "objective_variant": str(candidate["objective_variant"]),
        "case_id": str(candidate["case_id"]),
        "family_key": str(candidate["family_key"]),
        "family_label": str(candidate["family_label"]),
        "blocks": int(candidate["blocks"]),
        "max_tones": int(candidate["max_tones"]),
        "order_label": str(candidate["order_label"]),
        "levels": [int(level) for level in candidate.get("levels", [])],
        "n12_full_target_fidelity": float(candidate["summary"]["n12_full_target_fidelity"]),
        "n12_ground_fidelity": float(candidate["summary"]["n12_ground_fidelity"]),
        "n12_ground_gain": float(candidate["summary"]["n12_ground_gain"]),
        "ground_convergence_delta": float(candidate["summary"]["ground_convergence_delta"]),
        "ancilla_excitation_worst": float(candidate["summary"]["ancilla_excitation_worst"]),
        "support_leakage_worst": float(candidate["summary"]["support_leakage_worst"]),
        "outside_target_leakage_worst": float(candidate["summary"]["outside_target_leakage_worst"]),
        "physically_plausible": bool(candidate["summary"]["physically_plausible"]),
    }


def write_candidate_csv(path: Path, candidates: list[dict[str, Any]]) -> None:
    fieldnames = [
        "candidate_id",
        "source",
        "objective_variant",
        "case_id",
        "family_key",
        "family_label",
        "blocks",
        "max_tones",
        "order_label",
        "levels",
        "n12_full_target_fidelity",
        "n12_ground_fidelity",
        "n12_ground_gain",
        "ground_convergence_delta",
        "ancilla_excitation_worst",
        "support_leakage_worst",
        "outside_target_leakage_worst",
        "physically_plausible",
    ]
    rows = [compact_candidate(candidate) for candidate in candidates]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            row_copy = dict(row)
            row_copy["levels"] = "-".join(str(level) for level in row_copy["levels"])
            writer.writerow(row_copy)


def plot_block_summary(baseline_candidates: list[dict[str, Any]], refined_candidates: list[dict[str, Any]]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2), sharey=True)
    family_order = ["drsqr", "drcpsqr"]
    for axis, family_key in zip(axes, family_order, strict=True):
        baseline = sorted(
            [candidate for candidate in baseline_candidates if candidate["family_key"] == family_key],
            key=lambda row: int(row["blocks"]),
        )
        refined = sorted(
            [candidate for candidate in refined_candidates if candidate["family_key"] == family_key],
            key=lambda row: int(row["blocks"]),
        )
        axis.axhline(GROUND_THRESHOLD, color="0.35", linestyle="--", linewidth=1.0, label="0.99 threshold")
        axis.plot(
            [int(candidate["blocks"]) for candidate in baseline],
            [float(candidate["summary"]["n12_ground_fidelity"]) for candidate in baseline],
            marker="o",
            linestyle="--",
            linewidth=1.5,
            label="Baseline sequences rescored",
        )
        axis.plot(
            [int(candidate["blocks"]) for candidate in refined],
            [float(candidate["summary"]["n12_ground_fidelity"]) for candidate in refined],
            marker="s",
            linestyle="-",
            linewidth=1.8,
            label="Ground-sector local refine",
        )
        axis.plot(
            [int(candidate["blocks"]) for candidate in refined],
            [float(candidate["summary"]["n12_full_target_fidelity"]) for candidate in refined],
            marker="^",
            linestyle=":",
            linewidth=1.2,
            label="Full-target fidelity after refine",
        )
        axis.set_title(FAMILY_LABELS[family_key])
        axis.set_xlabel("Block count")
        axis.set_xticks(sorted({int(candidate["blocks"]) for candidate in baseline + refined}))
        axis.grid(True, alpha=0.25)
    axes[0].set_ylabel("N=12 fidelity")
    handles, labels = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.92))
    fig.savefig(c.FIG_DIR / "ground_sector_block_summary.png", dpi=300, bbox_inches="tight")
    fig.savefig(c.FIG_DIR / "ground_sector_block_summary.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_objective_scatter(candidates: list[dict[str, Any]]) -> None:
    fig, ax = plt.subplots(figsize=(5.8, 4.6))
    markers = {"drsqr": "o", "drcpsqr": "s"}
    for family_key, family_label in FAMILY_LABELS.items():
        family_candidates = [candidate for candidate in candidates if candidate["family_key"] == family_key]
        ax.scatter(
            [float(candidate["summary"]["n12_full_target_fidelity"]) for candidate in family_candidates],
            [float(candidate["summary"]["n12_ground_fidelity"]) for candidate in family_candidates],
            marker=markers[family_key],
            s=55,
            alpha=0.85,
            label=family_label,
        )
        for candidate in family_candidates:
            ax.annotate(
                f"b{int(candidate['blocks'])}",
                (
                    float(candidate["summary"]["n12_full_target_fidelity"]),
                    float(candidate["summary"]["n12_ground_fidelity"]),
                ),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=8,
            )
    ax.plot([0.0, 1.0], [0.0, 1.0], linestyle="--", color="0.35", linewidth=1.0)
    ax.set_xlabel("N=12 full-target fidelity")
    ax.set_ylabel("N=12 ground-sector fidelity")
    ax.set_xlim(0.0, 1.01)
    ax.set_ylim(0.0, 1.01)
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(c.FIG_DIR / "ground_sector_objective_scatter.png", dpi=300, bbox_inches="tight")
    fig.savefig(c.FIG_DIR / "ground_sector_objective_scatter.pdf", bbox_inches="tight")
    plt.close(fig)


def save_family_artifacts(family_outcomes: dict[str, Any], candidate_lookup: dict[str, dict[str, Any]]) -> None:
    for family_key, path in FAMILY_ARTIFACTS.items():
        best = family_outcomes[family_key]["best_absolute"]
        candidate = candidate_lookup[str(best["candidate_id"])]
        payload = {
            "study_name": "cluster_state_holographic_unified",
            "date_created": date.today().isoformat(),
            "description": f"Ground-sector best candidate for {FAMILY_LABELS[family_key]}",
            "parameters": {
                "family_key": family_key,
                "blocks": int(candidate["blocks"]),
                "max_tones": int(candidate["max_tones"]),
                "order_label": str(candidate["order_label"]),
                "levels": [int(level) for level in candidate.get("levels", [])],
                "objective_variant": str(candidate["objective_variant"]),
            },
            "load_instructions": "import json; from pathlib import Path; payload = json.loads(Path(filename).read_text(encoding='utf-8-sig'))",
            "candidate": compact_candidate(candidate),
            "sequence": candidate["sequence"],
        }
        c.save_json(path, payload)


def build_summary_payload(
    *,
    baseline_candidates: list[dict[str, Any]],
    refined_candidates: list[dict[str, Any]],
    family_outcomes: dict[str, Any],
) -> dict[str, Any]:
    return {
        "study_name": "cluster_state_holographic_unified",
        "date_created": date.today().isoformat(),
        "objective": {
            "description": "Restricted ground-sector transfer objective U_joint(|g>⊗|psi>) ≈ |g>⊗H_c|psi> on {|g,0>, |g,1>}.",
            "target_matrix": c.GROUND_SECTOR_TARGET_UNITARY,
            "threshold": float(GROUND_THRESHOLD),
            "truncation_checks": list(TRUNCATION_VALUES),
            "physical_plausibility_limits": {
                "ancilla_excitation_worst": float(PHYSICAL_ANCILLA_LIMIT),
                "outside_target_leakage_worst": float(PHYSICAL_OUTSIDE_TARGET_LIMIT),
                "ground_convergence_delta": float(PHYSICAL_CONVERGENCE_LIMIT),
            },
        },
        "baseline_rescored_candidates": [compact_candidate(candidate) for candidate in baseline_candidates],
        "ground_refined_candidates": [compact_candidate(candidate) for candidate in refined_candidates],
        "family_outcomes": family_outcomes,
        "cross_family_summary": {
            "best_absolute": max(
                (family_outcomes[family_key]["best_absolute"] for family_key in FAMILY_LABELS),
                key=lambda row: float(row["n12_ground_fidelity"]),
            ),
            "best_plausible": max(
                (family_outcomes[family_key]["most_physically_plausible"] for family_key in FAMILY_LABELS),
                key=lambda row: float(row["n12_ground_fidelity"]),
            ),
            "comparison_note": (
                "Baseline full-target winners were rescored under the ground-sector objective, then the best seed per family/block was locally refined with the restricted target at N_cav = 12."
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skip-refine", action="store_true", help="Only rescore the stored block representatives without local ground-sector refinement.")
    parser.add_argument("--maxiter", type=int, default=GROUND_REFINE_MAXITER, help="Maximum iterations for each local ground-sector refinement.")
    args = parser.parse_args()

    records = load_scope_records()
    seed_records = select_seed_records(records)

    baseline_candidates = [
        evaluate_candidate_record(record, source="baseline_rescore", objective_variant="full_target_seed")
        for record in seed_records
    ]

    refined_candidates: list[dict[str, Any]] = []
    if args.skip_refine:
        refined_candidates = list(baseline_candidates)
    else:
        for index, seed_record in enumerate(seed_records, start=1):
            print(
                f"[ground-refine {index}/{len(seed_records)}] {seed_record['family_key']} blocks={seed_record['blocks']} case={seed_record['case_id']}",
                flush=True,
            )
            refined_record = refine_record_for_ground_sector(seed_record, maxiter=int(args.maxiter), seed=GROUND_REFINE_SEED)
            refined_candidates.append(
                evaluate_candidate_record(refined_record, source="ground_refine", objective_variant="ground_sector")
            )

    family_outcomes = select_family_outcomes(refined_candidates)
    candidate_lookup = {candidate["candidate_id"]: candidate for candidate in refined_candidates}

    write_candidate_csv(c.DATA_DIR / "ground_sector_followup_candidates.csv", refined_candidates)
    plot_block_summary(baseline_candidates, refined_candidates)
    plot_objective_scatter(refined_candidates)
    save_family_artifacts(family_outcomes, candidate_lookup)

    summary_payload = build_summary_payload(
        baseline_candidates=baseline_candidates,
        refined_candidates=refined_candidates,
        family_outcomes=family_outcomes,
    )
    c.save_json(c.DATA_DIR / "ground_sector_followup_summary.json", summary_payload)


if __name__ == "__main__":
    main()