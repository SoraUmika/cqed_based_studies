"""Resume the staged multitone SQR/CPSQR study from saved case artifacts."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from analysis import CaseRequest, build_case_context, save_case_artifact
from common import load_json
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
from run_study import (
    COMPARISON_PATH,
    CSV_PATH,
    DATA_DIR,
    DURATION_PATH,
    NEGATIVE_CONTROL_PATH,
    RESULTS_PATH,
    SCREEN_FAMILIES,
    SCREEN_PATH,
    SUMMARY_PATH,
    comparison_case_requests,
    duration_case_requests,
    duration_families_from_comparison,
    screen_case_requests,
)


ARTIFACT_CASE_DIR = Path(__file__).resolve().parents[1] / "artifacts" / "cases"


def artifact_path(case_id: str, family_name: str) -> Path:
    return ARTIFACT_CASE_DIR / f"{case_id}_{family_name}.json"


def load_or_run_case(request: CaseRequest, family_name: str, context_cache: dict[str, object]) -> dict:
    path = artifact_path(request.case_id, family_name)
    if path.exists():
        payload = load_json(path)
        print(f"  [resume-hit] {family_name}")
        return dict(payload["summary_row"])

    print(f"  [resume-run] {family_name}")
    if request.case_id not in context_cache:
        context_cache[request.case_id] = build_case_context(request)
    context = context_cache[request.case_id]
    row, artifact = family_runner(family_name)(context)  # type: ignore[arg-type]
    save_case_artifact(request.case_id, family_name, artifact, artifact.get("waveform_samples"))
    return row


def execute_stage_with_resume(requests: list[CaseRequest], families: tuple[str, ...]) -> list[dict]:
    rows: list[dict] = []
    context_cache: dict[str, object] = {}
    for request in requests:
        print(f"[case] {request.stage}: {request.case_id}")
        for family_name in families:
            row = load_or_run_case(request, family_name, context_cache)
            rows.append(row)
            print(
                "    strict_joint={:.4f} cpsqr_joint={:.4f} strict_reduced_quartet={:.4f} full_quartet={:.4f}".format(
                    row["strict_joint_process_fidelity"],
                    row["cpsqr_joint_process_fidelity"],
                    row["strict_reduced_quartet_mean"],
                    row["strict_full_quartet_mean"],
                )
            )
    return rows


def write_outputs(screen_rows: list[dict], comparison_rows: list[dict], duration_rows: list[dict]) -> None:
    screen_df = pd.DataFrame(screen_rows)
    comparison_df = pd.DataFrame(comparison_rows)
    duration_df = pd.DataFrame(duration_rows)

    negative_control = reevaluate_legacy_control()
    NEGATIVE_CONTROL_PATH.write_text(json.dumps(negative_control, indent=2), encoding="utf-8")

    selected_families = select_families(screen_df)
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


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    screen_rows = execute_stage_with_resume(screen_case_requests(), SCREEN_FAMILIES)
    screen_df = pd.DataFrame(screen_rows)
    selected_families = select_families(screen_df)
    print(f"[selected] {selected_families}")

    comparison_rows = execute_stage_with_resume(comparison_case_requests(), selected_families)
    comparison_df = pd.DataFrame(comparison_rows)
    duration_families = duration_families_from_comparison(comparison_df)
    print(f"[duration families] {duration_families}")

    duration_rows = execute_stage_with_resume(duration_case_requests(), duration_families)
    write_outputs(screen_rows, comparison_rows, duration_rows)
    print("Resume pass complete.")


if __name__ == "__main__":
    main()
