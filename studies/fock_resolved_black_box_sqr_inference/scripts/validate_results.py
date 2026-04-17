"""Validate outputs for the Fock-resolved black-box SQR inference study."""

from __future__ import annotations

import json
from pathlib import Path

from common import DATA_DIR, save_json


RESULTS_JSON = DATA_DIR / "study_results.json"
SUMMARY_JSON = DATA_DIR / "study_summary.json"
VALIDATION_JSON = DATA_DIR / "validation_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    results = load_json(RESULTS_JSON)
    summary = load_json(SUMMARY_JSON)
    convergence = results["convergence"]

    checks = [
        {
            "name": "single_qubit_baseline_high_shot",
            "passed": bool(float(summary["best_single_qubit_mle_mean_fidelity"]) >= 0.995),
            "details": {"best_single_qubit_mle_mean_fidelity": float(summary["best_single_qubit_mle_mean_fidelity"])},
        },
        {
            "name": "wait_only_identifiable_displacement_only_not",
            "passed": bool(
                int(summary["identifiability"]["wait_only"]["transverse_rank"]) == 8
                and int(summary["identifiability"]["displacement_only"]["transverse_rank"]) < 8
            ),
            "details": summary["identifiability"],
        },
        {
            "name": "coherence_witness_separates_wait_from_combined",
            "passed": bool(
                float(summary["coherence_wait_residual"]) < 1.0e-10
                and float(summary["coherence_combined_residual"]) > 1.0e-2
            ),
            "details": {
                "coherence_wait_residual": float(summary["coherence_wait_residual"]),
                "coherence_combined_residual": float(summary["coherence_combined_residual"]),
            },
        },
        {
            "name": "explicit_gauge_family_exists",
            "passed": bool(int(summary["exact_gauge_family_count"]) >= 1),
            "details": {"exact_gauge_family_count": int(summary["exact_gauge_family_count"])},
        },
        {
            "name": "exact_recoverable_fit_converges",
            "passed": bool(
                float(convergence["wait_only_exact_rmse"]) < 1.0e-8
                and float(convergence["combined_dense_grid_rmse"]) < 1.0e-8
            ),
            "details": convergence,
        },
        {
            "name": "literature_comparison",
            "passed": True,
            "details": {
                "status": "not_applicable",
                "reason": "This is an original identifiability and protocol-design study rather than a literature reproduction benchmark.",
            },
        },
    ]

    payload = {
        "all_passed": bool(all(check["passed"] for check in checks)),
        "checks": checks,
    }
    save_json(VALIDATION_JSON, payload)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
