"""Validation checks for the hybrid gate-set comparison study.

This script performs the Step 4 checks required by the study workflow:

1. sanity checks on the nominal benchmark outputs,
2. convergence checks versus cavity truncation and GRAPE time slicing,
3. optimization-stability spot checks for the leading decomposition ansaetze.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from common import (
    DATA_DIR,
    GrapeCase,
    SequenceCase,
    grape_cases,
    json_dump,
    library_a_entangler_sequence,
    library_a_local_sequence,
    run_grape_case,
    run_sequence_case,
)


SUMMARY_PATH = DATA_DIR / "study_summary.json"


def _load_summary() -> dict[str, Any]:
    if not SUMMARY_PATH.exists():
        raise SystemExit("Run scripts/run_study.py first.")
    return json.loads(SUMMARY_PATH.read_text(encoding="utf-8"))


def _results_by_label(summary: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {row["label"]: row for row in summary["results"]}


def _sanity_checks(summary: dict[str, Any]) -> dict[str, Any]:
    rows = _results_by_label(summary)
    checks = {
        "A_local_high_fidelity": float(rows["A_local"]["strict_fidelity"]) > 0.98,
        "A_ent_exact_entangler": float(rows["A_ent"]["strict_fidelity"]) > 0.999,
        "D_ent_exact_entangler": float(rows["D_ent"]["strict_fidelity"]) > 0.999,
        "B_ent_local_equivalence": float(rows["B_ent"]["block_fidelity"]) > 0.99,
        "ECD_not_best_in_fock_local": float(rows["C_local"]["strict_fidelity"]) < float(rows["A_local"]["strict_fidelity"]),
        "GRAPE_local_beats_selective_B_local": float(rows["E_local_320"]["strict_fidelity"]) > float(rows["B_local"]["strict_fidelity"]),
        "SWAP_native_local_beats_selective_B_local": float(rows["F_local"]["strict_fidelity"]) > float(rows["B_local"]["strict_fidelity"]),
        "SWAP_native_entangler_beats_selective_B_ent_strict": float(rows["F_ent"]["strict_fidelity"]) > float(rows["B_ent"]["strict_fidelity"]),
    }
    return {"checks": checks, "passed": all(checks.values())}


def _truncation_checks() -> dict[str, Any]:
    results: dict[str, list[dict[str, float]]] = {}
    passed = True
    for label in ("A_local", "A_ent"):
        rows: list[dict[str, float]] = []
        for n_cav in (6, 8, 10):
            case = (
                SequenceCase(
                    label="A_local",
                    sequence=library_a_local_sequence(n_cav=n_cav),
                    description="Validation rerun for cavity-local benchmark.",
                    target_key="local_h",
                    multistart=3,
                    maxiter=180,
                )
                if label == "A_local"
                else SequenceCase(
                    label="A_ent",
                    sequence=library_a_entangler_sequence(),
                    description="Validation rerun for entangler benchmark.",
                    target_key="cx_c_to_q",
                    multistart=2,
                    maxiter=120,
                )
            )
            record = run_sequence_case(case, n_cav=n_cav)
            rows.append(
                {
                    "n_cav": float(n_cav),
                    "strict_fidelity": float(record["strict_fidelity"]),
                    "block_fidelity": float(record["block_fidelity"]),
                    "leakage_average": float(record["leakage_average"]),
                }
            )
        strict_values = [row["strict_fidelity"] for row in rows]
        spread = max(strict_values) - min(strict_values)
        results[label] = rows
        passed = passed and (spread < 0.02)
    return {"checks": results, "passed": passed}


def _time_grid_checks() -> dict[str, Any]:
    selected = {
        "local_h": {"duration_s": 480.0e-9, "steps": (20, 24, 32), "maxiter": 220},
        "cx_c_to_q": {"duration_s": 400.0e-9, "steps": (20, 24, 28), "maxiter": 220},
    }
    grouped: dict[str, list[dict[str, float]]] = {"local_h": [], "cx_c_to_q": []}
    for target_key, spec in selected.items():
        for steps in spec["steps"]:
            best_record: dict[str, float] | None = None
            for seed in (17, 23, 31):
                case = GrapeCase(
                    label=f"val_{target_key}_{steps}_{seed}",
                    target_key=target_key,
                    duration_s=spec["duration_s"],
                    steps=int(steps),
                    maxiter=int(spec["maxiter"]),
                    seed=int(seed),
                )
                record = run_grape_case(case)
                candidate = {
                    "steps": float(steps),
                    "seed": float(seed),
                    "strict_fidelity": float(record["strict_fidelity"]),
                    "block_fidelity": float(record["block_fidelity"]),
                    "leakage_average": float(record["leakage_average"]),
                }
                if best_record is None or candidate["strict_fidelity"] > best_record["strict_fidelity"]:
                    best_record = candidate
            assert best_record is not None
            grouped[target_key].append(best_record)

    local_values = [row["strict_fidelity"] for row in grouped["local_h"]]
    ent_values = [row["strict_fidelity"] for row in grouped["cx_c_to_q"]]
    passed = (max(local_values) - min(local_values) < 0.03) and (max(ent_values) - min(ent_values) < 0.05)
    return {"checks": grouped, "passed": passed}


def _optimization_stability_checks(summary: dict[str, Any]) -> dict[str, Any]:
    rows = _results_by_label(summary)
    reruns = {
        "A_local_deeper": SequenceCase(
            label="A_local_deeper",
            sequence=library_a_local_sequence(n_cav=8),
            description="Higher-budget rerun for stability.",
            target_key="local_h",
            multistart=5,
            maxiter=260,
        ),
        "A_ent_deeper": SequenceCase(
            label="A_ent_deeper",
            sequence=library_a_entangler_sequence(),
            description="Higher-budget rerun for stability.",
            target_key="cx_c_to_q",
            multistart=4,
            maxiter=180,
        ),
    }
    checks: dict[str, dict[str, float]] = {}
    passed = True
    for label, case in reruns.items():
        record = run_sequence_case(case)
        baseline_label = "A_local" if "local" in label else "A_ent"
        delta = abs(float(record["strict_fidelity"]) - float(rows[baseline_label]["strict_fidelity"]))
        checks[label] = {
            "strict_fidelity": float(record["strict_fidelity"]),
            "baseline_strict_fidelity": float(rows[baseline_label]["strict_fidelity"]),
            "delta": float(delta),
        }
        passed = passed and (delta < 0.01)
    return {"checks": checks, "passed": passed}


def main() -> None:
    summary = _load_summary()
    validation = {
        "sanity": _sanity_checks(summary),
        "convergence_truncation": _truncation_checks(),
        "convergence_time_grid": _time_grid_checks(),
        "optimization_stability": _optimization_stability_checks(summary),
    }
    validation["all_passed"] = all(section["passed"] for section in validation.values())
    json_dump(DATA_DIR / "validation_summary.json", validation)

    print(json.dumps(validation, indent=2, sort_keys=True))
    if not validation["all_passed"]:
        raise SystemExit("Validation checks failed.")


if __name__ == "__main__":
    main()
