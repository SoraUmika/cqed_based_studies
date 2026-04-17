"""Validation checks for the strict no-detuning multitone-SQR study."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from common import (
    DATA_DIR,
    build_model,
    build_square_multitone_waveform,
    build_target_operator,
    decoupled_block_operator,
    duration_from_chi_t,
    evaluate_square_multitone,
    logical_levels,
    make_run_config,
    optimize_square_multitone,
    save_json,
    target_spec,
)
from run_study import CaseRequest, row_from_operator


OUTPUT_PATH = DATA_DIR / "validation_details.json"


def representative_case() -> CaseRequest:
    return CaseRequest(
        model_variant="chi_plus_chiprime",
        include_chi_prime=True,
        family="structured_xy",
        n_active=3,
        chi_t_over_2pi=3.0,
    )


def aligned_echo_case() -> CaseRequest:
    return CaseRequest(
        model_variant="chi_only",
        include_chi_prime=False,
        family="aligned_x",
        n_active=2,
        chi_t_over_2pi=1.0,
    )


def main() -> None:
    rep = representative_case()
    rep_spec = target_spec(rep.family, rep.n_active)
    rep_model = build_model(include_chi_prime=rep.include_chi_prime, n_active=rep.n_active)
    rep_duration_s = duration_from_chi_t(rep.chi_t_over_2pi)
    rep_levels = logical_levels(rep.n_active)
    rep_target = build_target_operator(rep_spec, rep_levels)

    low_budget = optimize_square_multitone(rep_model, rep_spec, make_run_config(rep_model, n_active=rep.n_active, duration_s=rep_duration_s), n_starts=1, maxiter=5, random_seed=9123, label_prefix="validation_low")
    high_budget = optimize_square_multitone(rep_model, rep_spec, make_run_config(rep_model, n_active=rep.n_active, duration_s=rep_duration_s), n_starts=3, maxiter=60, random_seed=9123, label_prefix="validation_high")
    low_row = row_from_operator(
        rep,
        construction="validation_low_budget",
        target_operator=rep_target,
        actual_operator=np.asarray(low_budget.validation.restricted_operator, dtype=np.complex128),
        validation=low_budget.validation,
    )
    high_row = row_from_operator(
        rep,
        construction="validation_high_budget",
        target_operator=rep_target,
        actual_operator=np.asarray(high_budget.validation.restricted_operator, dtype=np.complex128),
        validation=high_budget.validation,
    )

    dt_rows = {}
    for dt_s in (1.0e-9, 2.0e-9, 4.0e-9):
        validation, _tone_specs, _waveform = evaluate_square_multitone(
            rep_model,
            rep_spec,
            make_run_config(rep_model, n_active=rep.n_active, duration_s=rep_duration_s, dt_s=dt_s),
            corrections=high_budget.corrections,
            label=f"validation_dt_{dt_s:.1e}",
        )
        row = row_from_operator(
            rep,
            construction=f"dt_{dt_s:.1e}",
            target_operator=rep_target,
            actual_operator=np.asarray(validation.restricted_operator, dtype=np.complex128),
            validation=validation,
        )
        dt_rows[f"{dt_s:.1e}"] = {
            "restricted_average_gate_fidelity": row["restricted_average_gate_fidelity"],
            "best_fit_restricted_process_fidelity": row["best_fit_restricted_process_fidelity"],
            "max_residual_z_error_rad": row["max_residual_z_error_rad"],
        }

    ntr_rows = {}
    for n_tr in (2, 3):
        model = build_model(include_chi_prime=rep.include_chi_prime, n_active=rep.n_active, n_tr=n_tr)
        validation, _tone_specs, _waveform = evaluate_square_multitone(
            model,
            rep_spec,
            make_run_config(model, n_active=rep.n_active, duration_s=rep_duration_s),
            corrections=high_budget.corrections,
            label=f"validation_ntr_{n_tr}",
        )
        row = row_from_operator(
            rep,
            construction=f"ntr_{n_tr}",
            target_operator=rep_target,
            actual_operator=np.asarray(validation.restricted_operator, dtype=np.complex128),
            validation=validation,
        )
        ntr_rows[str(n_tr)] = {
            "restricted_average_gate_fidelity": row["restricted_average_gate_fidelity"],
            "best_fit_restricted_process_fidelity": row["best_fit_restricted_process_fidelity"],
            "max_residual_z_error_rad": row["max_residual_z_error_rad"],
        }

    ideal_waveform, ideal_tone_specs = build_square_multitone_waveform(
        rep_model,
        rep_spec,
        make_run_config(rep_model, n_active=rep.n_active, duration_s=rep_duration_s),
        corrections=None,
        label="validation_decoupled",
    )
    decoupled = decoupled_block_operator(ideal_tone_specs, levels=rep_levels, duration_s=rep_duration_s)
    decoupled_row = row_from_operator(
        rep,
        construction="decoupled_validation",
        target_operator=rep_target,
        actual_operator=decoupled,
    )

    aligned = aligned_echo_case()
    aligned_spec = target_spec(aligned.family, aligned.n_active)
    aligned_model = build_model(include_chi_prime=aligned.include_chi_prime, n_active=aligned.n_active)
    aligned_duration_s = duration_from_chi_t(aligned.chi_t_over_2pi)
    aligned_run = make_run_config(aligned_model, n_active=aligned.n_active, duration_s=aligned_duration_s)
    aligned_opt = optimize_square_multitone(aligned_model, aligned_spec, aligned_run, n_starts=3, maxiter=60, random_seed=9555, label_prefix="validation_aligned")

    payload = {
        "representative_case": {
            "case_id": rep.case_id,
            "low_budget_restricted_average_gate_fidelity": float(low_row["restricted_average_gate_fidelity"]),
            "high_budget_restricted_average_gate_fidelity": float(high_row["restricted_average_gate_fidelity"]),
            "low_budget_best_fit_process": float(low_row["best_fit_restricted_process_fidelity"]),
            "high_budget_best_fit_process": float(high_row["best_fit_restricted_process_fidelity"]),
            "high_budget_max_residual_z_error_rad": float(max(high_row["per_block_residual_z_error_rad"])),
        },
        "dt_convergence": dt_rows,
        "ntr_convergence": ntr_rows,
        "decoupled_exact_match": {
            "restricted_average_gate_fidelity": float(decoupled_row["restricted_average_gate_fidelity"]),
            "restricted_process_fidelity": float(decoupled_row["restricted_process_fidelity"]),
        },
        "aligned_x_reference": {
            "case_id": aligned.case_id,
            "plain_restricted_average_gate_fidelity": float(
                row_from_operator(
                    aligned,
                    construction="validation_aligned",
                    target_operator=build_target_operator(aligned_spec, logical_levels(aligned.n_active)),
                    actual_operator=np.asarray(aligned_opt.validation.restricted_operator, dtype=np.complex128),
                    validation=aligned_opt.validation,
                )["restricted_average_gate_fidelity"]
            ),
            "plain_best_fit_restricted_process_fidelity": float(aligned_opt.validation.best_fit_restricted_process_fidelity),
        },
    }
    save_json(OUTPUT_PATH, payload, description="Validation-detail checks for the strict no-detuning multitone-SQR study.")


if __name__ == "__main__":
    main()
