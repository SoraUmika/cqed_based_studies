from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFINITIVE_STUDY = Path(__file__).resolve().parents[1]
SOURCE_ROOT = DEFINITIVE_STUDY.parent


STUDY_PATHS = {
    "study1": SOURCE_ROOT / "multitone_sqr_arbitrary_fock_conditional_rotations",
    "study2": SOURCE_ROOT / "parameterized_waveform_residual_z_cancellation",
    "study3": SOURCE_ROOT / "ideal_sqr_direct_vs_echoed_multitone",
    "study4": SOURCE_ROOT / "native_rich_multitone_sqr_cpsqr_feasibility",
}


def _approx_parameter_count(construction: str, n_active: int) -> int:
    name = str(construction)
    n = int(n_active)
    if name in {"direct_multitone", "baseline_multitone", "gaussian_seed", "native_direct_strict", "reduced_unitary_direct"}:
        return 3 * n
    if name in {"complex_envelope", "basis_expanded"}:
        return 3 * n + 6
    if name in {"symmetric_two_segment"}:
        return 6 * n
    if name in {"echoed_symmetric"}:
        return 6 * n + 1
    if name in {"echoed_independent", "echoed_asymmetric", "echoed_multitone", "echoed_cpsqr"}:
        return 6 * n + 3
    if "echo" in name:
        return 6 * n + 2
    return 3 * n


def construction_display(name: str) -> str:
    mapping = {
        "direct_multitone": "direct",
        "baseline_multitone": "direct",
        "native_direct_strict": "direct",
        "reduced_unitary_direct": "direct",
        "gaussian_seed": "direct seed",
        "symmetric_two_segment": "two-segment direct",
        "complex_envelope": "rich envelope",
        "basis_expanded": "basis-expanded",
        "echoed_symmetric": "symmetric echo",
        "echoed_independent": "independent echo",
        "echoed_asymmetric": "asymmetric echo",
        "echoed_multitone": "echoed",
        "echoed_cpsqr": "CPSQR echo",
        "single_pulse": "single pulse",
        "echoed_fixed_total": "symmetric echo",
        "echoed_fixed_active": "asymmetric echo",
    }
    return mapping.get(str(name), str(name).replace("_", " "))


def construction_family(name: str) -> str:
    label = construction_display(name).lower()
    if "independent echo" in label:
        return "independent echo"
    if "symmetric echo" in label or label == "echoed":
        return "symmetric echo"
    if "asymmetric echo" in label:
        return "asymmetric echo"
    if "hybrid" in label:
        return "hybrid"
    return "direct"


def _save_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _normalize_study1() -> pd.DataFrame:
    source = STUDY_PATHS["study1"] / "data" / "echo_comparison_results.csv"
    df = pd.read_csv(source)
    rows = []
    for item in df.to_dict(orient="records"):
        artifact_rel = item.get("baseline_artifact")
        artifact_path = None if pd.isna(artifact_rel) else STUDY_PATHS["study1"] / str(artifact_rel)
        construction = item.get("sequence_family", item.get("family", "single_pulse"))
        rows.append(
            {
                "study": "multitone_sqr_arbitrary_fock_conditional_rotations",
                "study_role": "prior",
                "study_short": "arbitrary SU(2) prior",
                "case_id": item["case_id"],
                "construction": construction,
                "construction_display": construction_display(str(construction)),
                "construction_family": construction_family(str(construction)),
                "ansatz": str(item.get("sequence_family", "single_pulse")),
                "envelope": "gaussian",
                "duration_ns": item.get("pulse_duration_ns"),
                "total_gate_duration_ns": item.get("total_gate_duration_ns", item.get("pulse_duration_ns")),
                "model_variant": item.get("model_variant"),
                "include_chi_prime": item.get("include_chi_prime"),
                "n_active": item.get("n_active"),
                "target_family": item.get("family"),
                "target_kind": "arbitrary_su2",
                "strict_process_fidelity": item.get("restricted_process_fidelity"),
                "avg_gate_fidelity": item.get("average_gate_fidelity"),
                "cpsqr_process_fidelity": None,
                "reduced_quartet_fidelity": None,
                "full_quartet_fidelity": None,
                "mean_residual_z_error_rad": item.get("mean_residual_z_error_rad"),
                "mean_transverse_error_rad": item.get("mean_transverse_error_rad"),
                "parameter_count": _approx_parameter_count(str(construction), int(item.get("n_active", 1))),
                "target_parameter_count": 2 * int(item.get("n_active", 1)),
                "artifact_path": None if artifact_path is None else str(artifact_path.resolve()),
                "notes": item.get("target_operator_rule"),
                "source_csv": str(source.resolve()),
            }
        )
    return pd.DataFrame(rows)


def _normalize_study2() -> pd.DataFrame:
    source = STUDY_PATHS["study2"] / "data" / "study_results.csv"
    df = pd.read_csv(source)
    rows = []
    for item in df.to_dict(orient="records"):
        artifact_path = STUDY_PATHS["study2"] / "artifacts" / "cases" / f"{item['case_id']}_{item['waveform_family']}.json"
        rows.append(
            {
                "study": "parameterized_waveform_residual_z_cancellation",
                "study_role": "prior",
                "study_short": "residual-Z prior",
                "case_id": item["case_id"],
                "construction": item["waveform_family"],
                "construction_display": construction_display(str(item["waveform_family"])),
                "construction_family": construction_family(str(item["waveform_family"])),
                "ansatz": str(item["waveform_family"]).replace("_", " "),
                "envelope": "gaussian" if "complex" not in str(item["waveform_family"]) else "rich envelope",
                "duration_ns": item.get("pulse_duration_ns"),
                "total_gate_duration_ns": item.get("pulse_duration_ns"),
                "model_variant": item.get("model_variant"),
                "include_chi_prime": item.get("include_chi_prime"),
                "n_active": item.get("n_active"),
                "target_family": item.get("target_family"),
                "target_kind": "arbitrary_su2",
                "strict_process_fidelity": item.get("restricted_process_fidelity"),
                "avg_gate_fidelity": item.get("average_gate_fidelity"),
                "cpsqr_process_fidelity": None,
                "reduced_quartet_fidelity": None,
                "full_quartet_fidelity": None,
                "mean_residual_z_error_rad": item.get("mean_residual_z_error_rad"),
                "mean_transverse_error_rad": item.get("mean_transverse_error_rad"),
                "parameter_count": _approx_parameter_count(str(item["waveform_family"]), int(item.get("n_active", 1))),
                "target_parameter_count": 2 * int(item.get("n_active", 1)),
                "artifact_path": str(artifact_path.resolve()) if artifact_path.exists() else None,
                "notes": "Residual-Z cancellation prior study",
                "source_csv": str(source.resolve()),
            }
        )
    return pd.DataFrame(rows)


def _normalize_study3() -> pd.DataFrame:
    source = STUDY_PATHS["study3"] / "data" / "study_results.csv"
    df = pd.read_csv(source)
    rows = []
    for item in df.to_dict(orient="records"):
        artifact_path = STUDY_PATHS["study3"] / "artifacts" / "cases" / f"{item['case_id']}_{item['construction']}.json"
        rows.append(
            {
                "study": "ideal_sqr_direct_vs_echoed_multitone",
                "study_role": "prior",
                "study_short": "ideal-SQR baseline",
                "case_id": item["case_id"],
                "construction": item["construction"],
                "construction_display": construction_display(str(item["construction"])),
                "construction_family": construction_family(str(item["construction"])),
                "ansatz": str(item["construction"]).replace("_", " "),
                "envelope": "gaussian",
                "duration_ns": item.get("pulse_duration_ns"),
                "total_gate_duration_ns": item.get("total_gate_duration_ns"),
                "model_variant": item.get("model_variant"),
                "include_chi_prime": item.get("include_chi_prime"),
                "n_active": item.get("n_active"),
                "target_family": item.get("target_family"),
                "target_kind": "ideal_sqr",
                "strict_process_fidelity": item.get("restricted_process_fidelity"),
                "avg_gate_fidelity": item.get("average_gate_fidelity"),
                "cpsqr_process_fidelity": None,
                "reduced_quartet_fidelity": None,
                "full_quartet_fidelity": None,
                "mean_residual_z_error_rad": item.get("mean_residual_z_error_rad"),
                "mean_transverse_error_rad": item.get("mean_transverse_error_rad"),
                "parameter_count": _approx_parameter_count(str(item["construction"]), int(item.get("n_active", 1))),
                "target_parameter_count": 2 * int(item.get("n_active", 1)),
                "artifact_path": str(artifact_path.resolve()) if artifact_path.exists() else None,
                "notes": item.get("target_operator_rule"),
                "source_csv": str(source.resolve()),
            }
        )
    return pd.DataFrame(rows)


def _normalize_study4() -> pd.DataFrame:
    source = STUDY_PATHS["study4"] / "data" / "study_results.csv"
    df = pd.read_csv(source)
    rows = []
    for item in df.to_dict(orient="records"):
        artifact_path = STUDY_PATHS["study4"] / "artifacts" / "cases" / f"{item['case_id']}_{item['family_name']}.json"
        rows.append(
            {
                "study": "native_rich_multitone_sqr_cpsqr_feasibility",
                "study_role": "new",
                "study_short": "native-rich extension",
                "case_id": item["case_id"],
                "construction": item["family_name"],
                "construction_display": construction_display(str(item["family_name"])),
                "construction_family": construction_family(str(item["family_name"])),
                "ansatz": str(item["family_name"]).replace("_", " "),
                "envelope": "gaussian" if item["family_name"] in {"gaussian_seed", "native_direct_strict", "reduced_unitary_direct"} else "rich envelope",
                "duration_ns": item.get("duration_ns"),
                "total_gate_duration_ns": item.get("total_gate_duration_ns"),
                "model_variant": item.get("model_variant"),
                "include_chi_prime": item.get("include_chi_prime"),
                "n_active": item.get("n_active"),
                "target_family": item.get("target_family"),
                "target_kind": "ideal_sqr",
                "strict_process_fidelity": item.get("strict_joint_process_fidelity"),
                "avg_gate_fidelity": item.get("strict_joint_average_gate_fidelity"),
                "cpsqr_process_fidelity": item.get("cpsqr_joint_process_fidelity"),
                "reduced_quartet_fidelity": item.get("strict_reduced_quartet_mean"),
                "full_quartet_fidelity": item.get("strict_full_quartet_mean"),
                "mean_residual_z_error_rad": item.get("strict_mean_residual_z_error_rad"),
                "mean_transverse_error_rad": item.get("strict_mean_transverse_error_rad"),
                "parameter_count": _approx_parameter_count(str(item["family_name"]), int(item.get("n_active", 1))),
                "target_parameter_count": 2 * int(item.get("n_active", 1)),
                "artifact_path": str(artifact_path.resolve()) if artifact_path.exists() else None,
                "notes": item.get("classification_label"),
                "source_csv": str(source.resolve()),
            }
        )
    return pd.DataFrame(rows)


def load_normalized_results() -> dict[str, pd.DataFrame]:
    return {
        "study1": _normalize_study1(),
        "study2": _normalize_study2(),
        "study3": _normalize_study3(),
        "study4": _normalize_study4(),
    }


def save_snapshots(normalized: dict[str, pd.DataFrame]) -> None:
    prior_dir = DEFINITIVE_STUDY / "data" / "prior_studies"
    _save_json(prior_dir / "study1_results.json", {"rows": normalized["study1"].to_dict(orient="records")})
    _save_json(prior_dir / "study2_results.json", {"rows": normalized["study2"].to_dict(orient="records")})
    _save_json(prior_dir / "study3_baseline_results.json", {"rows": normalized["study3"].to_dict(orient="records")})


def combined_dataframe(normalized: dict[str, pd.DataFrame]) -> pd.DataFrame:
    return pd.concat([normalized["study1"], normalized["study2"], normalized["study3"], normalized["study4"]], ignore_index=True)


def main() -> None:
    normalized = load_normalized_results()
    save_snapshots(normalized)
    combined = combined_dataframe(normalized)
    output = DEFINITIVE_STUDY / "data" / "prior_studies" / "combined_prior_results.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output, index=False)


if __name__ == "__main__":
    main()
