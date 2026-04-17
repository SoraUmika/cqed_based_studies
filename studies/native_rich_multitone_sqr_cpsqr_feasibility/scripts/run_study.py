"""Run the final native/rich multitone SQR/CPSQR feasibility study."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

import runtime_compat  # noqa: F401

from analysis import CaseRequest, build_case_context, save_case_artifact
from families import family_runner
from reporting import (
    build_summary,
    make_figures,
    reevaluate_legacy_control,
    select_families,
    write_markdown_summary,
    write_notebook,
    write_report,
)


STAGE1_MODELS = (("chi_plus_chiprime", True),)
STAGE2_MODELS = (("chi_only", False), ("chi_plus_chiprime", True))
STAGE_STRUCTURED_TARGETS = ("smooth_x", "staggered_x")
STAGE_DURATIONS = (1.0, 3.0, 5.0)
DURATION_REFINEMENT = (0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0)
RANDOM_X_SEEDS = (4101, 4102)

SCREEN_FAMILIES = (
    "gaussian_seed",
    "native_direct_strict",
    "reduced_unitary_direct",
    "symmetric_two_segment",
    "complex_envelope",
    "basis_expanded",
    "echoed_symmetric",
    "echoed_independent",
    "echoed_asymmetric",
    "echoed_cpsqr",
)

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
RESULTS_PATH = DATA_DIR / "all_results.json"
SCREEN_PATH = DATA_DIR / "screen_results.json"
COMPARISON_PATH = DATA_DIR / "comparison_results.json"
DURATION_PATH = DATA_DIR / "duration_results.json"
SUMMARY_PATH = DATA_DIR / "study_summary.json"
CSV_PATH = DATA_DIR / "study_results.csv"
NEGATIVE_CONTROL_PATH = DATA_DIR / "negative_controls.json"


def screen_case_requests() -> list[CaseRequest]:
    rows: list[CaseRequest] = []
    for model_variant, include_chi_prime in STAGE1_MODELS:
        for n_active in (2, 3):
            for chi_t in STAGE_DURATIONS:
                for target_family in STAGE_STRUCTURED_TARGETS:
                    rows.append(
                        CaseRequest(
                            stage="screen",
                            model_variant=model_variant,
                            include_chi_prime=include_chi_prime,
                            target_family=target_family,
                            n_active=n_active,
                            chi_t_over_2pi=chi_t,
                        )
                    )
    return rows


def comparison_case_requests() -> list[CaseRequest]:
    rows: list[CaseRequest] = []
    for model_variant, include_chi_prime in STAGE2_MODELS:
        for n_active in (2, 3, 4):
            for chi_t in STAGE_DURATIONS:
                for target_family in STAGE_STRUCTURED_TARGETS:
                    rows.append(
                        CaseRequest(
                            stage="comparison",
                            model_variant=model_variant,
                            include_chi_prime=include_chi_prime,
                            target_family=target_family,
                            n_active=n_active,
                            chi_t_over_2pi=chi_t,
                        )
                    )
                if n_active == 4:
                    for seed in RANDOM_X_SEEDS:
                        rows.append(
                            CaseRequest(
                                stage="comparison",
                                model_variant=model_variant,
                                include_chi_prime=include_chi_prime,
                                target_family="random_x",
                                n_active=n_active,
                                chi_t_over_2pi=chi_t,
                                random_seed=int(seed),
                            )
                        )
    return rows


def duration_case_requests() -> list[CaseRequest]:
    rows: list[CaseRequest] = []
    for model_variant, include_chi_prime in STAGE2_MODELS:
        for n_active in (2, 3):
            for chi_t in DURATION_REFINEMENT:
                rows.append(
                    CaseRequest(
                        stage="duration",
                        model_variant=model_variant,
                        include_chi_prime=include_chi_prime,
                        target_family="smooth_x",
                        n_active=n_active,
                        chi_t_over_2pi=chi_t,
                    )
                )
    return rows


def execute_stage(requests: list[CaseRequest], families: tuple[str, ...]) -> list[dict]:
    rows: list[dict] = []
    for request in requests:
        print(f"[case] {request.stage}: {request.case_id}")
        context = build_case_context(request)
        for family_name in families:
            print(f"  [family] {family_name}")
            row, artifact = family_runner(family_name)(context)
            rows.append(row)
            save_case_artifact(request.case_id, family_name, artifact, artifact.get("waveform_samples"))
            print(
                "    strict_joint={:.4f} cpsqr_joint={:.4f} strict_reduced_quartet={:.4f} full_quartet={:.4f}".format(
                    row["strict_joint_process_fidelity"],
                    row["cpsqr_joint_process_fidelity"],
                    row["strict_reduced_quartet_mean"],
                    row["strict_full_quartet_mean"],
                )
            )
    return rows


def duration_families_from_comparison(comparison_df: pd.DataFrame) -> tuple[str, ...]:
    scores = (
        comparison_df.groupby("family_name", as_index=False)[["strict_joint_process_fidelity", "cpsqr_joint_process_fidelity"]]
        .mean()
        .assign(score=lambda frame: frame[["strict_joint_process_fidelity", "cpsqr_joint_process_fidelity"]].max(axis=1))
        .sort_values("score", ascending=False)
    )
    return tuple(scores.head(2)["family_name"].tolist())


def run_single(args: argparse.Namespace) -> None:
    request = CaseRequest(
        stage="single",
        model_variant=str(args.single_model),
        include_chi_prime=str(args.single_model) == "chi_plus_chiprime",
        target_family=str(args.single_target),
        n_active=int(args.single_na),
        chi_t_over_2pi=float(args.single_chiT),
    )
    context = build_case_context(request)
    row, artifact = family_runner(str(args.single_family))(context)
    save_case_artifact(request.case_id, str(args.single_family), artifact, artifact.get("waveform_samples"))
    print(json.dumps(row, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--single-family", type=str, default="")
    parser.add_argument("--single-target", type=str, default="smooth_x")
    parser.add_argument("--single-model", type=str, default="chi_plus_chiprime")
    parser.add_argument("--single-na", type=int, default=2)
    parser.add_argument("--single-chiT", type=float, default=3.0)
    args = parser.parse_args()

    if args.single_family:
        run_single(args)
        return

    screen_rows = execute_stage(screen_case_requests(), SCREEN_FAMILIES)
    screen_df = pd.DataFrame(screen_rows)
    selected_families = select_families(screen_df)
    print(f"[selected] {selected_families}")

    comparison_rows = execute_stage(comparison_case_requests(), selected_families)
    comparison_df = pd.DataFrame(comparison_rows)
    duration_families = duration_families_from_comparison(comparison_df)
    print(f"[duration families] {duration_families}")

    duration_rows = execute_stage(duration_case_requests(), duration_families)
    duration_df = pd.DataFrame(duration_rows)

    negative_control = reevaluate_legacy_control()
    NEGATIVE_CONTROL_PATH.write_text(json.dumps(negative_control, indent=2), encoding="utf-8")

    summary = build_summary(screen_df, comparison_df, duration_df, selected_families, negative_control)
    RESULTS_PATH.write_text(json.dumps({"rows": screen_rows + comparison_rows + duration_rows}, indent=2), encoding="utf-8")
    SCREEN_PATH.write_text(json.dumps({"rows": screen_rows}, indent=2), encoding="utf-8")
    COMPARISON_PATH.write_text(json.dumps({"rows": comparison_rows}, indent=2), encoding="utf-8")
    DURATION_PATH.write_text(json.dumps({"rows": duration_rows}, indent=2), encoding="utf-8")
    SUMMARY_PATH.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    pd.DataFrame(screen_rows + comparison_rows + duration_rows).to_csv(CSV_PATH, index=False)

    make_figures(comparison_df)
    write_markdown_summary(summary)
    write_report(summary)
    write_notebook(summary)


if __name__ == "__main__":
    main()
