"""Focused validation and regression checks for the extension pass."""

from __future__ import annotations

import numpy as np
import pandas as pd

from common import (
    ARTIFACTS_DIR,
    DATA_DIR,
    build_model,
    build_target_operator,
    corrections_from_vector,
    duration_from_chi_t,
    evaluate_square_multitone,
    load_json,
    logical_levels,
    make_run_config,
    reduced_blockwise_operator,
    target_spec,
    average_gate_fidelity,
    process_fidelity,
)
from run_extension_study import extension_requests, target_for_request
from run_study import CaseRequest, row_from_operator


CASE_DIR = ARTIFACTS_DIR / "cases"
OUTPUT_PATH = DATA_DIR / "extension_validation_summary.json"


def archived_baseline_regression() -> dict[str, float]:
    request = CaseRequest(
        model_variant="chi_only",
        include_chi_prime=False,
        family="aligned_x",
        n_active=2,
        chi_t_over_2pi=1.0,
    )
    archived_results = load_json(DATA_DIR / "study_results.json")["case_rows"]
    archived_row = next(
        row
        for row in archived_results
        if row["case_id"] == request.case_id and row["construction"] == "full_shared_line"
    )
    artifact = load_json(CASE_DIR / "chi_only_aligned_x_na2_chiT1p0_full_shared_line.json")
    spec = target_spec("aligned_x", request.n_active)
    model = build_model(include_chi_prime=request.include_chi_prime, n_active=request.n_active)
    run_config = make_run_config(model, n_active=request.n_active, duration_s=duration_from_chi_t(request.chi_t_over_2pi))
    corrections = corrections_from_vector(np.asarray(artifact["corrections_vector"], dtype=float), n_active=request.n_active)
    validation, _tone_specs, _waveform = evaluate_square_multitone(
        model,
        spec,
        run_config,
        corrections=corrections,
        label="extension_regression_anchor",
    )
    row = row_from_operator(
        request,
        construction="full_shared_line",
        target_operator=build_target_operator(spec, logical_levels(request.n_active)),
        actual_operator=np.asarray(validation.restricted_operator, dtype=np.complex128),
        validation=validation,
    )
    return {
        "restricted_average_gate_fidelity_archived": float(archived_row["restricted_average_gate_fidelity"]),
        "restricted_average_gate_fidelity_rerun": float(row["restricted_average_gate_fidelity"]),
        "restricted_average_gate_fidelity_diff": float(row["restricted_average_gate_fidelity"] - archived_row["restricted_average_gate_fidelity"]),
        "max_residual_z_error_rad_archived": float(archived_row["max_residual_z_error_rad"]),
        "max_residual_z_error_rad_rerun": float(row["max_residual_z_error_rad"]),
        "max_residual_z_error_rad_diff": float(row["max_residual_z_error_rad"] - archived_row["max_residual_z_error_rad"]),
        "best_fit_restricted_process_fidelity_archived": float(archived_row["best_fit_restricted_process_fidelity"]),
        "best_fit_restricted_process_fidelity_rerun": float(row["best_fit_restricted_process_fidelity"]),
        "best_fit_restricted_process_fidelity_diff": float(
            row["best_fit_restricted_process_fidelity"] - archived_row["best_fit_restricted_process_fidelity"]
        ),
    }


def tuned_case_sensitivity() -> dict[str, object]:
    root_chi_t = float(load_json(DATA_DIR / "extension_tuned_set_map.json")["equal_amplitude_roots"][0]["chi_t_over_2pi"])
    request = extension_requests(root_chi_t)[0]
    spec = target_for_request(request)
    model = build_model(include_chi_prime=request.include_chi_prime, n_active=request.n_active)
    levels = logical_levels(request.n_active)
    target_operator = build_target_operator(spec, levels)
    artifact = load_json(CASE_DIR / f"extension_{request.case_id}_full_shared_line.json")
    corrections = corrections_from_vector(np.asarray(artifact["corrections_vector"], dtype=float), n_active=request.n_active)
    duration_s = duration_from_chi_t(request.chi_t_over_2pi)

    rows = []
    for dt_s in (1.0e-9, 2.0e-9, 4.0e-9):
        run_config = make_run_config(model, n_active=request.n_active, duration_s=duration_s, dt_s=dt_s)
        validation, _tone_specs, waveform = evaluate_square_multitone(
            model,
            spec,
            run_config,
            corrections=corrections,
            label=f"extension_validation_dt_{dt_s:.1e}",
        )
        full_restricted = np.asarray(validation.restricted_operator, dtype=np.complex128)
        reduced_operator = reduced_blockwise_operator(
            model,
            validation.compiled,
            waveform,
            run_config,
            levels=levels,
        )
        row = row_from_operator(
            request,
            construction=f"dt_{dt_s:.1e}",
            target_operator=target_operator,
            actual_operator=full_restricted,
            validation=validation,
            extra={
                "reduced_vs_full_restricted_process_fidelity": float(process_fidelity(full_restricted, reduced_operator)),
                "reduced_vs_full_restricted_average_gate_fidelity": float(average_gate_fidelity(full_restricted, reduced_operator)),
            },
        )
        rows.append(
            {
                "dt_s": float(dt_s),
                "restricted_average_gate_fidelity": float(row["restricted_average_gate_fidelity"]),
                "best_fit_restricted_process_fidelity": float(row["best_fit_restricted_process_fidelity"]),
                "max_residual_z_error_rad": float(row["max_residual_z_error_rad"]),
                "reduced_vs_full_restricted_average_gate_fidelity": float(row["reduced_vs_full_restricted_average_gate_fidelity"]),
                "reduced_vs_full_restricted_process_fidelity": float(row["reduced_vs_full_restricted_process_fidelity"]),
            }
        )
    return {"case_id": request.case_id, "dt_rows": rows}


def echo_dominance_audit() -> dict[str, object]:
    rows = load_json(DATA_DIR / "extension_echo_summary.json")["rows"]
    df = pd.DataFrame(rows)
    finite_df = df[df["construction"].isin(["echo_finite_gaussian", "echo_finite_manifold_aware_multitone"])]
    findings = []
    for family in sorted(finite_df["family"].unique()):
        family_df = df[df["family"] == family].copy()
        plain = family_df[family_df["construction"] == "full_shared_line"]
        if plain.empty:
            continue
        plain_row = plain.iloc[0]
        for construction in ("echo_finite_gaussian", "echo_finite_manifold_aware_multitone"):
            candidate = family_df[family_df["construction"] == construction]
            if candidate.empty:
                continue
            candidate_row = candidate.iloc[0]
            findings.append(
                {
                    "family": str(family),
                    "construction": str(construction),
                    "fidelity_delta_vs_plain": float(
                        candidate_row["restricted_average_gate_fidelity"] - plain_row["restricted_average_gate_fidelity"]
                    ),
                    "max_residual_z_delta_vs_plain": float(
                        candidate_row["max_residual_z_error_rad"] - plain_row["max_residual_z_error_rad"]
                    ),
                    "beats_plain_on_both_metrics": bool(
                        candidate_row["restricted_average_gate_fidelity"] > plain_row["restricted_average_gate_fidelity"]
                        and candidate_row["max_residual_z_error_rad"] < plain_row["max_residual_z_error_rad"]
                    ),
                }
            )
    return {
        "rows_checked": int(len(findings)),
        "findings": findings,
        "any_finite_echo_beats_plain_on_both_metrics": bool(
            any(item["beats_plain_on_both_metrics"] for item in findings)
        ),
    }


def main() -> None:
    payload = {
        "archived_baseline_regression": archived_baseline_regression(),
        "tuned_case_sensitivity": tuned_case_sensitivity(),
        "echo_dominance_audit": echo_dominance_audit(),
    }
    from common import save_json

    save_json(
        OUTPUT_PATH,
        payload,
        description="Focused validation and regression checks for the extension-pass outputs.",
    )


if __name__ == "__main__":
    main()