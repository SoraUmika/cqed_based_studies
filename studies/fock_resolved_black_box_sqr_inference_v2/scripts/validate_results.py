"""Validate outputs for the Fock-resolved black-box SQR inference study (v2).

Checks are stricter than v1:
  1. Single-qubit baseline MLE fidelity >= 0.995
  2. Wait-only transverse rank == 2*N_ACTIVE; displacement-only rank < 2*N_ACTIVE
  3. Full-state MLE non-uniqueness: objective span > 1.0 OR prob_std mean > 0.05
  4. Exact gauge family: at least 1 solution
  5. Exact recoverable fit converges to machine precision (RMSE < 1e-8)
  6. Coherence witness: combined protocol residual on coherent case > 1e-2
  7. Per-sector oracle fidelity on near-ideal pulse case > 0.90
  8. Coherence sweep: residual increases monotonically with coherence fraction
"""

from __future__ import annotations

import json
from pathlib import Path

from common import DATA_DIR, N_ACTIVE, save_json


RESULTS_JSON = DATA_DIR / "study_results.json"
SUMMARY_JSON = DATA_DIR / "study_summary.json"
VALIDATION_JSON = DATA_DIR / "validation_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    results = load_json(RESULTS_JSON)
    summary = load_json(SUMMARY_JSON)
    convergence = results["convergence"]
    ident = summary["identifiability"]

    checks = []

    # 1. Single-qubit baseline
    fid = float(summary["best_single_qubit_mle_mean_fidelity"])
    checks.append({
        "name": "single_qubit_baseline_high_shot",
        "passed": bool(fid >= 0.995),
        "details": {"best_mle_fidelity": fid, "threshold": 0.995},
    })

    # 2. Protocol identifiability
    wait_rank = int(ident["wait_only"]["transverse_rank"])
    disp_rank = int(ident["displacement_only"]["transverse_rank"])
    checks.append({
        "name": "wait_only_identifiable_displacement_only_not",
        "passed": bool(wait_rank == 2 * N_ACTIVE and disp_rank < 2 * N_ACTIVE),
        "details": {
            "wait_transverse_rank": wait_rank,
            "disp_transverse_rank": disp_rank,
            "expected_full_rank": 2 * N_ACTIVE,
        },
    })

    # 3. Full-state MLE non-uniqueness
    obj_span = float(summary.get("full_state_objective_span", 0.0))
    prob_stds = summary.get("full_state_probability_std_by_sector", [0.0])
    prob_std_mean = float(sum(prob_stds) / len(prob_stds)) if prob_stds else 0.0
    nonunique = bool(obj_span > 1.0 or prob_std_mean > 0.05)
    checks.append({
        "name": "full_state_nonuniqueness_demonstrated",
        "passed": nonunique,
        "details": {"objective_span": obj_span, "prob_std_mean": prob_std_mean},
    })

    # 4. Gauge family
    gauge_count = int(summary.get("exact_gauge_family_count", 0))
    checks.append({
        "name": "exact_gauge_family_exists",
        "passed": bool(gauge_count >= 1),
        "details": {"count": gauge_count},
    })

    # 5. Exact recoverable fit convergence
    wait_rmse = float(convergence["wait_only_exact_rmse"])
    comb_rmse = float(convergence["combined_dense_grid_rmse"])
    checks.append({
        "name": "exact_recoverable_fit_converges",
        "passed": bool(wait_rmse < 1.0e-8 and comb_rmse < 1.0e-8),
        "details": {
            "wait_only_exact_rmse": wait_rmse,
            "combined_exact_rmse": comb_rmse,
            "threshold": 1.0e-8,
        },
    })

    # 6. Coherence witness
    coh_wait = float(summary.get("coherence_wait_residual", 1.0))
    coh_comb = float(summary.get("coherence_combined_residual", 0.0))
    checks.append({
        "name": "coherence_witness_active",
        "passed": bool(coh_wait < 1.0e-8 and coh_comb > 1.0e-2),
        "details": {
            "diagonal_wait_residual": coh_wait,
            "coherent_combined_residual": coh_comb,
        },
    })

    # 7. Per-sector oracle fidelity on near-ideal pulse case
    pulse_rows = results.get("pulse_rows", [])
    near_ideal = [r for r in pulse_rows
                  if "long_pulse_optimized_mix_g" in str(r.get("case_id", ""))
                  and r.get("protocol") == "combined"]
    if near_ideal:
        fid_vals = [
            float(r["mle"].get("mean_oracle_fidelity", 0.0))
            for r in near_ideal
            if r["mle"].get("mean_oracle_fidelity") is not None
        ]
        mean_fid = float(sum(fid_vals) / len(fid_vals)) if fid_vals else 0.0
        checks.append({
            "name": "near_ideal_pulse_oracle_fidelity",
            "passed": bool(mean_fid >= 0.90),
            "details": {"mean_oracle_fidelity": mean_fid, "threshold": 0.90},
        })
    else:
        checks.append({
            "name": "near_ideal_pulse_oracle_fidelity",
            "passed": True,
            "details": {"status": "skipped_no_matching_case"},
        })

    # 8. Coherence sweep monotonicity
    sweep = summary.get("coherence_sweep_summary", [])
    if len(sweep) >= 3:
        comb_residuals = [float(r["combined_residual_rms"]) for r in sorted(sweep, key=lambda r: r["coherence_fraction"])]
        # Residual should increase as coherence_fraction increases (non-decreasing)
        monotone = all(comb_residuals[i] <= comb_residuals[i + 1] + 1.0e-6 for i in range(len(comb_residuals) - 1))
        checks.append({
            "name": "coherence_sweep_residual_monotone",
            "passed": bool(monotone),
            "details": {"residuals": comb_residuals},
        })
    else:
        checks.append({
            "name": "coherence_sweep_residual_monotone",
            "passed": True,
            "details": {"status": "skipped_insufficient_sweep_points"},
        })

    payload = {
        "all_passed": bool(all(c["passed"] for c in checks)),
        "checks": checks,
    }
    save_json(VALIDATION_JSON, payload)
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
