from __future__ import annotations

import csv
import dataclasses
import importlib.util
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Any, Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np
import qutip as qt

from cqed_sim import SimulationConfig, prepare_simulation
from cqed_sim.calibration import ConditionedQubitTargets, build_block_rotation_target_operator
from cqed_sim.calibration.conditioned_multitone import (
    ConditionedMultitoneCorrections,
    ConditionedMultitoneRunConfig,
    MultitoneTone,
    build_conditioned_multitone_waveform,
    compile_conditioned_multitone_waveform,
)
from cqed_sim.core.conventions import qubit_cavity_block_indices


ROOT = Path(__file__).resolve().parents[3]
STUDIES_ROOT = ROOT / "studies"
STUDY_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = STUDY_ROOT / "data"
FIGURES_DIR = STUDY_ROOT / "figures"
ARTIFACTS_DIR = STUDY_ROOT / "artifacts"
REPORT_DIR = STUDY_ROOT / "report"
SCRIPTS_DIR = STUDY_ROOT / "scripts"

PROBE_QUBIT_STATES: dict[str, np.ndarray] = {
    "g": np.asarray([1.0, 0.0], dtype=np.complex128),
    "e": np.asarray([0.0, 1.0], dtype=np.complex128),
    "plus_x": np.asarray([1.0, 1.0], dtype=np.complex128) / np.sqrt(2.0),
    "plus_y": np.asarray([1.0, 1.0j], dtype=np.complex128) / np.sqrt(2.0),
}

PROBE_TIERS: dict[str, tuple[str, ...]] = {
    "single_ground": ("g",),
    "selected_pair": ("g", "plus_x"),
    "spanning_quartet": ("g", "e", "plus_x", "plus_y"),
}


@dataclass(frozen=True)
class CaseSpec:
    study_key: str
    study_dir: str
    case_artifact: str
    plot_label: str
    study_goal: str
    target_summary: str
    notes: str

    @property
    def artifact_path(self) -> Path:
        return ROOT / self.case_artifact

    @property
    def study_root(self) -> Path:
        return ROOT / self.study_dir

    @property
    def common_py(self) -> Path:
        return self.study_root / "scripts" / "common.py"


CASE_SPECS: tuple[CaseSpec, ...] = (
    CaseSpec(
        study_key="parameterized",
        study_dir="studies/parameterized_waveform_residual_z_cancellation",
        case_artifact="studies/parameterized_waveform_residual_z_cancellation/artifacts/cases/chi_plus_chiprime_na4_chiT3p0_targetD_seed94340_baseline_multitone.json",
        plot_label="Residual-Z study\nbaseline multitone",
        study_goal="Compare richer waveform families for residual-Z suppression while still attempting the target conditional rotations.",
        target_summary="Random target-D block-diagonal conditional qubit rotations on 4 addressed Fock levels.",
        notes="Representative hard native multitone baseline from the richer-waveform comparison.",
    ),
    CaseSpec(
        study_key="arbitrary",
        study_dir="studies/multitone_sqr_arbitrary_fock_conditional_rotations",
        case_artifact="studies/multitone_sqr_arbitrary_fock_conditional_rotations/artifacts/cases/chi_plus_chiprime_na4_chiT3p0_familyC.json",
        plot_label="Arbitrary-rotation study\ndirect multitone",
        study_goal="Fit arbitrary block-diagonal conditional qubit rotations with a native multitone SQR waveform.",
        target_summary="Structured family-C arbitrary conditional SU(2) target on 4 addressed Fock levels.",
        notes="Representative direct multitone case from the arbitrary-target study.",
    ),
    CaseSpec(
        study_key="ideal",
        study_dir="studies/ideal_sqr_direct_vs_echoed_multitone",
        case_artifact="studies/ideal_sqr_direct_vs_echoed_multitone/artifacts/cases/chi_plus_chiprime_smooth_x_na3_chiT3p0_direct_multitone.json",
        plot_label="Ideal-SQR study\ndirect multitone",
        study_goal="Test whether direct or echoed multitone constructions can realize an ideal x-axis SQR profile.",
        target_summary="Smooth x-axis ideal-SQR target on 3 addressed Fock levels.",
        notes="Representative direct native multitone case from the ideal-SQR feasibility study.",
    ),
    CaseSpec(
        study_key="corrected",
        study_dir="studies/corrected_sqr_conditioned_rotation_metric",
        case_artifact="studies/corrected_sqr_conditioned_rotation_metric/artifacts/representative_case_artifact.json",
        plot_label="Corrected-metric study\nunitary-optimized",
        study_goal="Re-optimize a corrected direct multitone SQR under the fixed phase convention using a reduced effective-unitary metric.",
        target_summary="Smooth corrected-SQR profile with explicit per-level (theta_n, phi_n) on 4 addressed Fock levels.",
        notes="Representative 4-level corrected-convention case, already re-optimized under the fixed package.",
    ),
)


_MODULE_CACHE: dict[str, ModuleType] = {}


def ensure_directories() -> None:
    for path in (DATA_DIR, FIGURES_DIR, ARTIFACTS_DIR, REPORT_DIR, SCRIPTS_DIR):
        path.mkdir(parents=True, exist_ok=True)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, payload: Any) -> None:
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def write_csv(path: Path, rows: Sequence[dict[str, Any]], fieldnames: Sequence[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames))
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def matrix_from_record(record: dict[str, Any]) -> np.ndarray:
    real = np.asarray(record["real"], dtype=float)
    imag = np.asarray(record["imag"], dtype=float)
    return np.asarray(real + 1.0j * imag, dtype=np.complex128)


def matrix_to_record(matrix: np.ndarray) -> dict[str, Any]:
    arr = np.asarray(matrix, dtype=np.complex128)
    return {
        "shape": list(arr.shape),
        "real": np.real(arr).tolist(),
        "imag": np.imag(arr).tolist(),
    }


def load_study_module(case: CaseSpec) -> ModuleType:
    cached = _MODULE_CACHE.get(case.study_key)
    if cached is not None:
        return cached
    module_path = case.common_py
    module_name = f"native_validation_{case.study_key}_common"
    sys.path.insert(0, str(module_path.parent))
    try:
        spec = importlib.util.spec_from_file_location(module_name, module_path)
        if spec is None or spec.loader is None:
            raise RuntimeError(f"Could not load module spec from {module_path}")
        module = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = module
        spec.loader.exec_module(module)
    finally:
        sys.path.pop(0)
    _MODULE_CACHE[case.study_key] = module
    return module


def n_active_from_artifact(case: CaseSpec, artifact: dict[str, Any]) -> int:
    if case.study_key == "corrected":
        return int(artifact["n_active"])
    case_request = artifact.get("case_request", {})
    return int(case_request["n_active"])


def include_chi_prime_from_artifact(case: CaseSpec, artifact: dict[str, Any]) -> bool:
    if case.study_key == "corrected":
        return False
    case_request = artifact.get("case_request", {})
    if "include_chi_prime" in case_request:
        return bool(case_request["include_chi_prime"])
    return str(case_request.get("model_variant", "")).lower() == "chi_plus_chiprime"


def duration_s_from_artifact(case: CaseSpec, artifact: dict[str, Any]) -> float:
    if case.study_key == "corrected":
        return float(artifact["duration_s"])
    summary = artifact.get("summary_row", {})
    if "pulse_duration_s" in summary:
        return float(summary["pulse_duration_s"])
    if "pulse_duration_ns" in summary:
        return 1.0e-9 * float(summary["pulse_duration_ns"])
    raise KeyError(f"No duration field found for {case.study_key}")


def build_model(case: CaseSpec, artifact: dict[str, Any]) -> Any:
    module = load_study_module(case)
    n_active = n_active_from_artifact(case, artifact)
    if case.study_key == "corrected":
        return module.build_model()
    return module.build_model(include_chi_prime=include_chi_prime_from_artifact(case, artifact), n_active=n_active)


def build_frame(model: Any) -> Any:
    from cqed_sim import FrameSpec

    return FrameSpec(omega_q_frame=float(model.omega_q), omega_c_frame=float(model.omega_c))


def build_corrected_run_config(case: CaseSpec, artifact: dict[str, Any], model: Any) -> ConditionedMultitoneRunConfig:
    duration_s = duration_s_from_artifact(case, artifact)
    frame = build_frame(model)
    dt_s = 4.0e-9
    sigma_fraction = 1.0 / 6.0
    base = ConditionedMultitoneRunConfig(
        frame=frame,
        duration_s=float(duration_s),
        dt_s=float(dt_s),
        sigma_fraction=float(sigma_fraction),
        tone_cutoff=1.0e-12,
        include_all_levels=False,
        max_step_s=float(dt_s),
        fock_fqs_hz=None,
    )
    module = load_study_module(case)
    if case.study_key == "corrected":
        legacy = module.build_run_config(model, duration_s=float(duration_s), dt_s=float(dt_s), sigma_fraction=float(sigma_fraction))
    elif hasattr(module, "make_run_config"):
        legacy = module.make_run_config(
            model,
            n_active=n_active_from_artifact(case, artifact),
            duration_s=float(duration_s),
            dt_s=float(dt_s),
            sigma_fraction=float(sigma_fraction),
        )
    else:
        legacy = base
    return dataclasses.replace(legacy, fock_fqs_hz=None)


def logical_levels(case: CaseSpec, artifact: dict[str, Any]) -> tuple[int, ...]:
    return tuple(range(n_active_from_artifact(case, artifact)))


def logical_indices(model: Any, levels: Sequence[int]) -> np.ndarray:
    indices: list[int] = []
    for level in levels:
        indices.extend(int(idx) for idx in qubit_cavity_block_indices(int(model.n_cav), int(level)))
    return np.asarray(indices, dtype=int)


def multitone_tones_from_rows(rows: Sequence[dict[str, Any]]) -> tuple[MultitoneTone, ...]:
    tone_specs: list[MultitoneTone] = []
    for row in rows:
        tone_specs.append(
            MultitoneTone(
                manifold=int(row.get("manifold", row.get("n"))),
                omega_rad_s=float(row["omega_rad_s"]),
                amp_rad_s=float(row["amp_rad_s"]),
                phase_rad=float(row["phase_rad"]),
            )
        )
    return tuple(tone_specs)


def corrections_from_record(record: dict[str, Any]) -> ConditionedMultitoneCorrections:
    d_omega_rad_s = record.get("d_omega_rad_s")
    if d_omega_rad_s is None:
        d_omega_rad_s = [2.0 * math.pi * float(x) for x in record.get("d_omega_hz", [])]
    return ConditionedMultitoneCorrections(
        d_lambda=tuple(float(x) for x in record.get("d_lambda", [])),
        d_alpha=tuple(float(x) for x in record.get("d_alpha", [])),
        d_omega_rad_s=tuple(float(x) for x in d_omega_rad_s),
    )


def build_target_operator(case: CaseSpec, artifact: dict[str, Any], levels: Sequence[int]) -> np.ndarray:
    if "target_operator" in artifact:
        return matrix_from_record(artifact["target_operator"])
    if case.study_key != "corrected":
        raise KeyError(f"No target operator found for {case.study_key}")
    rows = [
        (float(item["theta_target_rad"]), float(item["phi_target_rad"]))
        for item in artifact["targets"]
    ]
    weights = [float(item.get("weight", 1.0 / len(rows))) for item in artifact["targets"]]
    targets = ConditionedQubitTargets.from_spec(rows, n_levels=len(rows), weights=weights)
    return np.asarray(build_block_rotation_target_operator(targets, logical_levels=levels), dtype=np.complex128)


def build_waveform(case: CaseSpec, artifact: dict[str, Any], model: Any, run_config: ConditionedMultitoneRunConfig) -> Any:
    if case.study_key == "parameterized":
        rows = artifact["validation"]["metadata"]["tone_specs"]
        return build_conditioned_multitone_waveform(multitone_tones_from_rows(rows), run_config, label=f"{case.study_key}_fixed_validation")
    if case.study_key == "arbitrary":
        rows = artifact["tone_specs"]
        return build_conditioned_multitone_waveform(multitone_tones_from_rows(rows), run_config, label=f"{case.study_key}_fixed_validation")
    if case.study_key == "ideal":
        module = load_study_module(case)
        family = str(artifact["target_spec"]["family"])
        spec = module.target_spec(family, n_active_from_artifact(case, artifact))
        corrections = corrections_from_record(artifact["optimizer"]["optimized_corrections"])
        waveform, _ = module.build_multitone_waveform_from_corrections(
            model,
            spec,
            run_config,
            corrections=corrections,
            label=f"{case.study_key}_fixed_validation",
        )
        return waveform
    if case.study_key == "corrected":
        rows = artifact["tone_specs"]
        return build_conditioned_multitone_waveform(multitone_tones_from_rows(rows), run_config, label=f"{case.study_key}_fixed_validation")
    raise ValueError(f"Unsupported case type {case.study_key}")


def compile_waveform(waveform: Any, run_config: ConditionedMultitoneRunConfig) -> Any:
    return compile_conditioned_multitone_waveform(waveform, run_config)


def simulation_session(model: Any, compiled: Any, *, frame: Any, drive_ops: dict[str, str]) -> Any:
    config = SimulationConfig(frame=frame, store_states=False)
    return prepare_simulation(model, compiled, drive_ops, config=config, e_ops={})


def simulate_full_operator_on_logical_inputs(
    model: Any,
    compiled: Any,
    *,
    frame: Any,
    drive_ops: dict[str, str],
    levels: Sequence[int],
) -> np.ndarray:
    session = simulation_session(model, compiled, frame=frame, drive_ops=drive_ops)
    full_dim = int(model.n_tr) * int(model.n_cav)
    operator = np.zeros((full_dim, full_dim), dtype=np.complex128)
    for level in levels:
        for qubit_level, input_index in enumerate(qubit_cavity_block_indices(int(model.n_cav), int(level))):
            result = session.run(model.basis_state(int(qubit_level), int(level)))
            operator[:, int(input_index)] = np.asarray(result.final_state.full(), dtype=np.complex128).reshape(-1)
    return operator


def restricted_operator_from_full(full_operator: np.ndarray, model: Any, levels: Sequence[int]) -> np.ndarray:
    indices = logical_indices(model, levels)
    return np.asarray(full_operator[np.ix_(indices, indices)], dtype=np.complex128)


def single_manifold_qubit_state(model: Any, level: int, qubit_state: Sequence[complex]) -> qt.Qobj:
    qvec = np.asarray(qubit_state, dtype=np.complex128)
    full = np.zeros(int(model.n_tr) * int(model.n_cav), dtype=np.complex128)
    idx_g, idx_e = qubit_cavity_block_indices(int(model.n_cav), int(level))
    full[int(idx_g)] = qvec[0]
    full[int(idx_e)] = qvec[1]
    norm_val = float(np.linalg.norm(full))
    if norm_val > 0.0:
        full = full / norm_val
    return qt.Qobj(full.reshape((-1, 1)), dims=[[int(model.n_tr), int(model.n_cav)], [1, 1]])


def apply_target_operator_to_state(model: Any, levels: Sequence[int], target_operator: np.ndarray, state: qt.Qobj) -> qt.Qobj:
    indices = logical_indices(model, levels)
    vector = np.asarray(state.full(), dtype=np.complex128).reshape(-1)
    logical_input = vector[indices]
    logical_output = np.asarray(target_operator, dtype=np.complex128) @ logical_input
    full_output = np.zeros_like(vector)
    full_output[indices] = logical_output
    return qt.Qobj(full_output.reshape((-1, 1)), dims=state.dims)


def state_fidelity(actual: qt.Qobj, target: qt.Qobj) -> float:
    return float(qt.fidelity(actual, target) ** 2)


def leakage_outside_logical(actual: qt.Qobj, model: Any, levels: Sequence[int]) -> float:
    indices = logical_indices(model, levels)
    vector = np.asarray(actual.full(), dtype=np.complex128).reshape(-1)
    logical_pop = float(np.sum(np.abs(vector[indices]) ** 2))
    return max(0.0, 1.0 - logical_pop)


def process_fidelity(target_operator: np.ndarray, actual_operator: np.ndarray) -> float:
    target = np.asarray(target_operator, dtype=np.complex128)
    actual = np.asarray(actual_operator, dtype=np.complex128)
    dim = float(target.shape[0])
    overlap = np.trace(target.conj().T @ actual)
    return float(np.clip(abs(overlap) ** 2 / (dim * dim), 0.0, 1.0))


def average_gate_fidelity(target_operator: np.ndarray, actual_operator: np.ndarray) -> float:
    dim = float(np.asarray(target_operator).shape[0])
    proc = process_fidelity(target_operator, actual_operator)
    return float((dim * proc + 1.0) / (dim + 1.0))


def mean_value(values: Iterable[float]) -> float:
    rows = list(float(v) for v in values)
    return float(np.mean(rows)) if rows else float("nan")


def min_value(values: Iterable[float]) -> float:
    rows = list(float(v) for v in values)
    return float(np.min(rows)) if rows else float("nan")


def max_value(values: Iterable[float]) -> float:
    rows = list(float(v) for v in values)
    return float(np.max(rows)) if rows else float("nan")


def legacy_restricted_process_fidelity(case: CaseSpec, artifact: dict[str, Any]) -> float | None:
    if case.study_key == "corrected":
        return None
    summary = artifact.get("summary_row", {})
    value = summary.get("restricted_process_fidelity")
    return None if value is None else float(value)


def initial_state_tier_rows(probe_rows: Sequence[dict[str, Any]]) -> dict[str, dict[str, float]]:
    tier_summary: dict[str, dict[str, float]] = {}
    for tier_name, labels in PROBE_TIERS.items():
        rows = [row for row in probe_rows if row["probe_label"] in labels]
        tier_summary[tier_name] = {
            "mean_state_fidelity": mean_value(row["state_fidelity"] for row in rows),
            "min_state_fidelity": min_value(row["state_fidelity"] for row in rows),
            "mean_leakage": mean_value(row["leakage_outside_logical"] for row in rows),
            "max_leakage": max_value(row["leakage_outside_logical"] for row in rows),
        }
    return tier_summary


def set_plot_style() -> None:
    plt.rcParams.update(
        {
            "figure.figsize": (7.0, 4.2),
            "font.size": 10,
            "axes.grid": False,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "savefig.dpi": 300,
        }
    )
