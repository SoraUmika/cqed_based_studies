"""Follow-up optimization and implementation study for hybrid qubit-cavity control.

This script extends the consolidated hybrid study with four follow-up tasks:

1. Restart-based diagnostics for the key decomposition-level candidates.
2. A deeper ansatz-size check for the structured U_target synthesis.
3. Pulse-level realization data for SNAP-, SQR-, and GRAPE-backed controls.
4. A numerical test of the constructive "selective pi + broadband pi" shortcut.

The script only uses cqed_sim-backed simulation and optimization. Study-local
helpers are imported only where cqed_sim still lacks first-class SNAP /
ConditionalPhaseSQR pulse builders.
"""

from __future__ import annotations

import copy
import importlib.util
import json
import statistics
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


STUDY_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = STUDY_ROOT / "data" / "followup_optimization"
FIG_DIR = STUDY_ROOT / "figures" / "followup_optimization"
DATA_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

GATE_SET_SCRIPTS = STUDY_ROOT / "scripts" / "gate_set_comparison"
LIT_SCRIPTS = STUDY_ROOT.parent / "literature_informed_selective_primitives" / "scripts"

for path in (GATE_SET_SCRIPTS, LIT_SCRIPTS):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

from cqed_sim.core import FrameSpec
from cqed_sim.io.gates import DisplacementGate, RotationGate, SQRGate
from cqed_sim.optimal_control import GrapeConfig, GrapeSolver, LeakagePenalty as OCLeakagePenalty
from cqed_sim.optimal_control import UnitaryObjective as OCUnitaryObjective
from cqed_sim.optimal_control import build_control_problem_from_model
from cqed_sim.pulses.builders import build_displacement_pulse, build_rotation_pulse, build_sqr_multitone_pulse
from cqed_sim.sequence import SequenceCompiler
from cqed_sim.sim import SimulationConfig, prepare_simulation
from cqed_sim.unitary_synthesis import (
    ConditionalPhaseSQR,
    Displacement,
    DriftPhaseModel,
    ExecutionOptions,
    GateSequence,
    LeakagePenalty,
    MultiObjective,
    QubitRotation,
    SQR,
    Subspace,
    TargetUnitary,
    UnitarySynthesizer,
    leakage_metrics,
    simulate_sequence,
    subspace_unitary_fidelity,
)

_GATE_COMMON_SPEC = importlib.util.spec_from_file_location("gate_common", GATE_SET_SCRIPTS / "common.py")
if _GATE_COMMON_SPEC is None or _GATE_COMMON_SPEC.loader is None:
    raise RuntimeError("Unable to load gate-set helper module.")
gate_common = importlib.util.module_from_spec(_GATE_COMMON_SPEC)
sys.modules["gate_common"] = gate_common
_GATE_COMMON_SPEC.loader.exec_module(gate_common)

CHI = gate_common.CHI
CHI_PRIME = gate_common.CHI_PRIME
KERR = gate_common.KERR
N_CAV_DEFAULT = gate_common.N_CAV_DEFAULT
_base_synth_kwargs = gate_common._base_synth_kwargs
_synthesis_target = gate_common._synthesis_target
build_frame = gate_common.build_frame
build_model = gate_common.build_model
grape_cases = gate_common.grape_cases
library_a_entangler_sequence = gate_common.library_a_entangler_sequence
library_a_local_sequence = gate_common.library_a_local_sequence
logical_subspace = gate_common.logical_subspace
per_fock_block_slices = gate_common.per_fock_block_slices
replay_grape_operator = gate_common.replay_grape_operator
target_matrix = gate_common.target_matrix

_LIT_COMMON_SPEC = importlib.util.spec_from_file_location("lit_common", LIT_SCRIPTS / "common.py")
if _LIT_COMMON_SPEC is None or _LIT_COMMON_SPEC.loader is None:
    raise RuntimeError("Unable to load literature selective-primitives helper module.")
lit_common = importlib.util.module_from_spec(_LIT_COMMON_SPEC)
sys.modules["lit_common"] = lit_common
_LIT_COMMON_SPEC.loader.exec_module(lit_common)

LIT_DT = lit_common.DT
average_target_state_fidelity = lit_common.average_target_state_fidelity
build_noise_spec = lit_common.build_noise_spec
build_selective_qubit_pulse = lit_common.build_selective_qubit_pulse
build_session = lit_common.build_session
cavity_ground_indices = lit_common.cavity_ground_indices
duration_from_chi_t = lit_common.duration_from_chi_t
embed_model_logical_state = lit_common.embed_model_logical_state
extract_restricted_operator = lit_common.extract_restricted_operator
lit_logical_indices = lit_common.logical_indices
qobj_probability_in_indices = lit_common.qobj_probability_in_indices
sample_total_waveform = lit_common.sample_total_waveform
snap_probe_state_vectors = lit_common.snap_probe_state_vectors
snap_target_operator = lit_common.snap_target_operator


TWO_PI = 2.0 * np.pi
FOLLOWUP_SEEDS = {
    "A_local": [1101, 1102, 1103, 1104],
    "D_ent": [1201, 1202, 1203, 1204],
    "L1c": [1301, 1302, 1303, 1304],
    "L1d": [1311, 1312, 1313],
    "L2c": [1401, 1402, 1403, 1404],
    "L2d": [1411, 1412, 1413],
}


@dataclass
class RestartResult:
    label: str
    seed: int
    result: Any
    process_fidelity: float
    block_fidelity: float
    operator_error: float
    leakage_average: float
    leakage_worst: float
    duration_ns: float
    parameter_count: int
    history: list[dict[str, Any]]


def _json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, complex):
        return {"real": float(value.real), "imag": float(value.imag)}
    if isinstance(value, (float, int, str, bool)) or value is None:
        return value
    return repr(value)


def _best_global_phase(target: np.ndarray, actual: np.ndarray) -> complex:
    overlap = np.trace(target.conj().T @ actual)
    if abs(overlap) < 1.0e-15:
        return 1.0 + 0.0j
    return overlap / abs(overlap)


def operator_error_norm(target: np.ndarray, actual: np.ndarray) -> float:
    phase = _best_global_phase(target, actual)
    delta = actual - phase * target
    return float(np.linalg.norm(delta, ord=2))


def average_gate_fidelity_from_process(f_process: float, dim: int) -> float:
    return float((dim * float(f_process) + 1.0) / (dim + 1.0))


def summarize_restart_batch(batch: list[RestartResult]) -> dict[str, Any]:
    process_values = [row.process_fidelity for row in batch]
    block_values = [row.block_fidelity for row in batch]
    leak_values = [row.leakage_average for row in batch]
    duration_values = [row.duration_ns for row in batch]
    sorted_batch = sorted(batch, key=lambda row: (row.process_fidelity, -row.leakage_average), reverse=True)
    best = sorted_batch[0]
    return {
        "best_seed": int(best.seed),
        "best_process_fidelity": float(best.process_fidelity),
        "median_process_fidelity": float(statistics.median(process_values)),
        "worst_process_fidelity": float(min(process_values)),
        "best_block_fidelity": float(best.block_fidelity),
        "median_block_fidelity": float(statistics.median(block_values)),
        "worst_block_fidelity": float(min(block_values)),
        "best_leakage_average": float(best.leakage_average),
        "median_leakage_average": float(statistics.median(leak_values)),
        "worst_leakage_average": float(max(leak_values)),
        "best_duration_ns": float(best.duration_ns),
        "median_duration_ns": float(statistics.median(duration_values)),
        "operator_error_best": float(best.operator_error),
        "parameter_count": int(best.parameter_count),
    }


def build_utarget_target() -> TargetUnitary:
    s = 1.0 / np.sqrt(2.0)
    u_target = np.array(
        [
            [s, 0.0, s, 0.0],
            [s, 0.0, -s, 0.0],
            [0.0, s, 0.0, s],
            [0.0, -s, 0.0, s],
        ],
        dtype=np.complex128,
    )
    return TargetUnitary(u_target, ignore_global_phase=True)


def build_utarget_subspace(n_cav: int = 8) -> Subspace:
    return Subspace.custom(
        2 * int(n_cav),
        [0, 1, int(n_cav), int(n_cav) + 1],
        ["|g,0>", "|g,1>", "|e,0>", "|e,1>"],
    )


def build_utarget_sequences() -> dict[str, GateSequence]:
    drift = DriftPhaseModel(chi=CHI, chi2=CHI_PRIME, kerr=KERR)
    n_cav = 8
    t_rot = 100e-9
    t_sqr = 2.0e-6
    t_disp = 200e-9
    t_cp = 100e-9

    l1c = GateSequence(
        gates=[
            Displacement(name="D1", alpha=0.3 + 0.0j, duration=t_disp),
            QubitRotation(name="R1", theta=np.pi / 2.0, phi=0.0, duration=t_rot),
            SQR(name="S1", theta_n=[0.0] * n_cav, phi_n=[0.0] * n_cav, drift_model=drift, duration=t_sqr),
            QubitRotation(name="R2", theta=np.pi / 2.0, phi=np.pi / 2.0, duration=t_rot),
            Displacement(name="D2", alpha=0.3 + 0.0j, duration=t_disp),
            QubitRotation(name="R3", theta=np.pi / 2.0, phi=0.0, duration=t_rot),
            SQR(name="S2", theta_n=[np.pi / 2.0] * n_cav, phi_n=[0.0] * n_cav, drift_model=drift, duration=t_sqr),
            QubitRotation(name="R4", theta=np.pi / 2.0, phi=np.pi / 2.0, duration=t_rot),
            Displacement(name="D3", alpha=0.3 + 0.0j, duration=t_disp),
            QubitRotation(name="R5", theta=np.pi / 4.0, phi=0.0, duration=t_rot),
            SQR(name="S3", theta_n=[np.pi / 4.0] * n_cav, phi_n=[0.0] * n_cav, drift_model=drift, duration=t_sqr),
        ],
        n_cav=n_cav,
    )
    l1d = GateSequence(gates=copy.deepcopy(l1c.gates) + [Displacement(name="D4", alpha=0.1 + 0.0j, duration=t_disp)], n_cav=n_cav)

    l2c = GateSequence(
        gates=[
            Displacement(name="D1", alpha=0.3 + 0.0j, duration=t_disp),
            QubitRotation(name="R1", theta=np.pi / 2.0, phi=0.0, duration=t_rot),
            ConditionalPhaseSQR(
                name="CP1",
                phases_n=[0.0, np.pi] + [0.0] * (n_cav - 2),
                drift_model=drift,
                duration=t_cp,
            ),
            QubitRotation(name="R2", theta=np.pi / 2.0, phi=np.pi / 2.0, duration=t_rot),
            Displacement(name="D2", alpha=0.3 + 0.0j, duration=t_disp),
            QubitRotation(name="R3", theta=np.pi / 2.0, phi=0.0, duration=t_rot),
            ConditionalPhaseSQR(
                name="CP2",
                phases_n=[0.0, np.pi / 2.0] + [0.0] * (n_cav - 2),
                drift_model=drift,
                duration=t_cp,
            ),
            QubitRotation(name="R4", theta=np.pi / 2.0, phi=np.pi / 2.0, duration=t_rot),
            Displacement(name="D3", alpha=0.3 + 0.0j, duration=t_disp),
            QubitRotation(name="R5", theta=np.pi / 4.0, phi=0.0, duration=t_rot),
            ConditionalPhaseSQR(
                name="CP3",
                phases_n=[0.0, np.pi / 4.0] + [0.0] * (n_cav - 2),
                drift_model=drift,
                duration=t_cp,
            ),
        ],
        n_cav=n_cav,
    )
    l2d = GateSequence(gates=copy.deepcopy(l2c.gates) + [Displacement(name="D4", alpha=0.1 + 0.0j, duration=t_disp)], n_cav=n_cav)
    return {"L1c": l1c, "L1d": l1d, "L2c": l2c, "L2d": l2d}


def synthesizer_for_sequence(
    *,
    sequence: GateSequence,
    target: TargetUnitary,
    subspace: Subspace,
    model,
    seed: int,
    optimizer: str = "powell",
) -> UnitarySynthesizer:
    kwargs = dict(_base_synth_kwargs())
    kwargs["optimizer"] = optimizer
    kwargs["optimizer_options"] = {"xatol": 1.0e-4, "fatol": 1.0e-5}
    kwargs["progress"] = {"enabled": True, "every": 1, "live": False}
    kwargs["seed"] = int(seed)
    kwargs["execution"] = ExecutionOptions(engine="auto", use_fast_path=True)
    return UnitarySynthesizer(
        subspace=subspace,
        primitives=copy.deepcopy(sequence.gates),
        target=target,
        model=model,
        **kwargs,
    )


def run_single_restart(
    *,
    label: str,
    sequence: GateSequence,
    target_key: str | None,
    target: TargetUnitary | None,
    subspace: Subspace,
    seed: int,
    maxiter: int,
    model,
    optimizer: str = "powell",
) -> RestartResult:
    resolved_target = target if target is not None else _synthesis_target(target_key)
    synth = synthesizer_for_sequence(
        sequence=sequence,
        target=resolved_target,
        subspace=subspace,
        model=model,
        seed=seed,
        optimizer=optimizer,
    )
    result = synth.fit(init_guess="random", multistart=1, maxiter=maxiter)
    actual = np.asarray(result.simulation.subspace_operator, dtype=np.complex128)
    target_matrix_local = np.asarray(resolved_target.matrix, dtype=np.complex128)
    block = float(
        subspace_unitary_fidelity(
            actual,
            target_matrix_local,
            gauge="block",
            block_slices=per_fock_block_slices(),
        )
    )
    leak = leakage_metrics(np.asarray(result.simulation.full_operator, dtype=np.complex128), subspace)
    return RestartResult(
        label=label,
        seed=int(seed),
        result=result,
        process_fidelity=float(subspace_unitary_fidelity(actual, target_matrix_local, gauge="global")),
        block_fidelity=block,
        operator_error=operator_error_norm(target_matrix_local, actual),
        leakage_average=float(leak.average),
        leakage_worst=float(leak.worst),
        duration_ns=float(result.sequence.total_duration() * 1.0e9),
        parameter_count=int(result.sequence.get_parameter_vector().size),
        history=list(result.history),
    )


def local_refine(
    *,
    label: str,
    seed: int,
    best_restart: RestartResult,
    target: TargetUnitary,
    subspace: Subspace,
    model,
    maxiter: int = 120,
    optimizer: str = "nelder_mead",
) -> RestartResult:
    synth = synthesizer_for_sequence(
        sequence=best_restart.result.sequence,
        target=target,
        subspace=subspace,
        model=model,
        seed=seed,
        optimizer=optimizer,
    )
    refined = synth.fit(init_guess="heuristic", multistart=1, maxiter=maxiter)
    actual = np.asarray(refined.simulation.subspace_operator, dtype=np.complex128)
    target_matrix_local = np.asarray(target.matrix, dtype=np.complex128)
    block = float(
        subspace_unitary_fidelity(
            actual,
            target_matrix_local,
            gauge="block",
            block_slices=per_fock_block_slices(),
        )
    )
    leak = leakage_metrics(np.asarray(refined.simulation.full_operator, dtype=np.complex128), subspace)
    return RestartResult(
        label=label,
        seed=int(seed),
        result=refined,
        process_fidelity=float(subspace_unitary_fidelity(actual, target_matrix_local, gauge="global")),
        block_fidelity=block,
        operator_error=operator_error_norm(target_matrix_local, actual),
        leakage_average=float(leak.average),
        leakage_worst=float(leak.worst),
        duration_ns=float(refined.sequence.total_duration() * 1.0e9),
        parameter_count=int(refined.sequence.get_parameter_vector().size),
        history=list(refined.history),
    )


def batch_to_payload(batch: list[RestartResult]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for row in batch:
        rows.append(
            {
                "label": row.label,
                "seed": int(row.seed),
                "process_fidelity": float(row.process_fidelity),
                "block_fidelity": float(row.block_fidelity),
                "avg_gate_fidelity": average_gate_fidelity_from_process(row.process_fidelity, 4),
                "operator_error": float(row.operator_error),
                "leakage_average": float(row.leakage_average),
                "leakage_worst": float(row.leakage_worst),
                "duration_ns": float(row.duration_ns),
                "parameter_count": int(row.parameter_count),
                "objective": float(row.result.objective),
                "sequence": row.result.sequence.serialize(),
                "history": row.history,
                "report": row.result.report,
            }
        )
    return rows


def make_restart_figure(payload: dict[str, Any]) -> None:
    labels = []
    best_vals = []
    med_vals = []
    worst_vals = []
    for label, record in payload.items():
        labels.append(label)
        best_vals.append(record["summary"]["best_process_fidelity"])
        med_vals.append(record["summary"]["median_process_fidelity"])
        worst_vals.append(record["summary"]["worst_process_fidelity"])

    x = np.arange(len(labels), dtype=float)
    fig, ax = plt.subplots(figsize=(10.5, 4.4))
    ax.bar(x, med_vals, color="#7fb3d5", edgecolor="black", linewidth=0.7, label="median")
    ax.scatter(x, best_vals, color="#1f78b4", s=55, zorder=3, label="best")
    ax.scatter(x, worst_vals, color="#e31a1c", s=55, zorder=3, label="worst")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=15, ha="right")
    ax.set_ylabel("Process fidelity")
    ax.set_ylim(0.45, 1.02)
    ax.set_title("Follow-up multistart diagnostics for key structured candidates")
    ax.legend(frameon=False)
    ax.grid(True, axis="y", alpha=0.25, linestyle=":")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "multistart_statistics.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / "multistart_statistics.pdf", bbox_inches="tight")
    plt.close(fig)


def make_convergence_figure(payload: dict[str, Any]) -> None:
    focus_labels = [label for label in ("A_local", "L1c", "L1d", "L2c", "L2d") if label in payload]
    ncols = 3
    nrows = int(np.ceil(len(focus_labels) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.1 * ncols, 3.8 * nrows))
    axes = np.atleast_1d(axes).ravel()
    for ax, label in zip(axes, focus_labels):
        runs = payload[label]["restarts"]
        for row in runs:
            history = row["history"]
            xs = [event["iteration"] for event in history]
            ys = [event["metrics"]["fidelity_subspace"] for event in history]
            ax.plot(xs, ys, alpha=0.45, linewidth=1.0, color="#4c78a8")
        best = max(runs, key=lambda row: row["process_fidelity"])
        xs_best = [event["iteration"] for event in best["history"]]
        ys_best = [event["metrics"]["fidelity_subspace"] for event in best["history"]]
        ax.plot(xs_best, ys_best, linewidth=2.3, color="#d62728", label="best run")
        ax.set_title(label)
        ax.set_xlabel("Iteration")
        ax.set_ylabel("Subspace fidelity")
        ax.set_ylim(0.45, 1.02)
        ax.grid(True, alpha=0.25, linestyle=":")
        ax.legend(frameon=False, loc="lower right")
    for ax in axes[len(focus_labels) :]:
        ax.axis("off")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "convergence_histories.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / "convergence_histories.pdf", bbox_inches="tight")
    plt.close(fig)


def make_ansatz_scaling_figure(existing_phase4: dict[str, Any], followup: dict[str, Any]) -> None:
    labels = ["L1a_d7", "L1b_d9", "L1c_d11", "L1d_d12", "L2a_d7", "L2b_d9", "L2c_d11", "L2d_d12"]
    values = [
        existing_phase4["L1a_D_R_SQR_d7"]["F_proj"],
        existing_phase4["L1b_D_R_SQR_d9"]["F_proj"],
        followup["L1c"]["summary"]["best_process_fidelity"],
        followup["L1d"]["summary"]["best_process_fidelity"],
        existing_phase4["L2a_D_R_CP_d7"]["F_proj"],
        existing_phase4["L2b_D_R_CP_d9"]["F_proj"],
        followup["L2c"]["summary"]["best_process_fidelity"],
        followup["L2d"]["summary"]["best_process_fidelity"],
    ]
    colors = ["#1f77b4"] * 4 + ["#ff7f0e"] * 4
    fig, ax = plt.subplots(figsize=(10.4, 4.4))
    ax.bar(labels, values, color=colors, edgecolor="black", linewidth=0.7)
    ax.set_ylabel("Process fidelity")
    ax.set_ylim(0.45, 1.02)
    ax.set_title("Ansatz-size scaling for structured U_target synthesis")
    ax.tick_params(axis="x", rotation=20)
    ax.grid(True, axis="y", alpha=0.25, linestyle=":")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "ansatz_scaling.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / "ansatz_scaling.pdf", bbox_inches="tight")
    plt.close(fig)


def sample_spectrum(time_s: np.ndarray, waveform: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    dt = float(np.mean(np.diff(time_s)))
    freqs_hz = np.fft.fftfreq(time_s.size, d=dt)
    spec = np.fft.fft(waveform)
    keep = freqs_hz >= 0.0
    return freqs_hz[keep], np.abs(spec[keep])


def sequence_from_serialized(
    serialized_sequence: list[dict[str, Any]],
    *,
    chi: float = CHI,
    chi_prime: float = CHI_PRIME,
    kerr: float = KERR,
    amp_scale: float = 1.0,
    dur_scale: float = 1.0,
    phase_offset: float = 0.0,
) -> GateSequence:
    drift = DriftPhaseModel(chi=float(chi), chi2=float(chi_prime), kerr=float(kerr))
    gates = []
    for gate in serialized_sequence:
        gate_type = str(gate["type"])
        name = str(gate["name"])
        duration = float(gate["duration"]) * float(dur_scale)
        params = list(gate.get("parameters", []))
        if gate_type == "Displacement":
            gates.append(
                Displacement(
                    name=name,
                    alpha=complex(float(params[0]) * float(amp_scale), float(params[1]) * float(amp_scale)),
                    duration=duration,
                    optimize_time=False,
                )
            )
        elif gate_type == "QubitRotation":
            gates.append(
                QubitRotation(
                    name=name,
                    theta=float(params[0]) * float(amp_scale),
                    phi=float(params[1]) + float(phase_offset),
                    duration=duration,
                    optimize_time=False,
                )
            )
        elif gate_type == "SQR":
            half = len(params) // 2
            theta_n = [float(x) * float(amp_scale) for x in params[:half]]
            phi_n = [float(x) + float(phase_offset) for x in params[half:]]
            gates.append(
                SQR(
                    name=name,
                    theta_n=theta_n,
                    phi_n=phi_n,
                    drift_model=drift,
                    duration=duration,
                    optimize_time=False,
                )
            )
        elif gate_type == "ConditionalPhaseSQR":
            phases_n = [float(x) + float(phase_offset) for x in params]
            gates.append(
                ConditionalPhaseSQR(
                    name=name,
                    phases_n=phases_n,
                    drift_model=drift,
                    duration=duration,
                    optimize_time=False,
                )
            )
        else:
            continue
    return GateSequence(gates=gates, n_cav=N_CAV_DEFAULT)


def evaluate_serialized_sequence(
    serialized_sequence: list[dict[str, Any]],
    *,
    target: TargetUnitary,
    subspace: Subspace,
    chi: float = CHI,
    chi_prime: float = CHI_PRIME,
    kerr: float = KERR,
    amp_scale: float = 1.0,
    dur_scale: float = 1.0,
    phase_offset: float = 0.0,
) -> dict[str, float]:
    seq = sequence_from_serialized(
        serialized_sequence,
        chi=chi,
        chi_prime=chi_prime,
        kerr=kerr,
        amp_scale=amp_scale,
        dur_scale=dur_scale,
        phase_offset=phase_offset,
    )
    sim = simulate_sequence(seq, subspace=subspace, backend="ideal")
    target_matrix_local = np.asarray(target.matrix, dtype=np.complex128)
    actual = np.asarray(sim.subspace_operator, dtype=np.complex128)
    leak = leakage_metrics(np.asarray(sim.full_operator, dtype=np.complex128), subspace)
    return {
        "process_fidelity": float(subspace_unitary_fidelity(actual, target_matrix_local, gauge="global")),
        "block_fidelity": float(
            subspace_unitary_fidelity(actual, target_matrix_local, gauge="block", block_slices=per_fock_block_slices())
        ),
        "leakage_average": float(leak.average),
    }


def _rms_drop(sweep: dict[float, dict[str, float]], nominal_key: str = "process_fidelity", band: float = 0.05) -> float:
    nominal = float(sweep[0.0][nominal_key])
    vals = [float(row[nominal_key]) - nominal for eps, row in sweep.items() if abs(float(eps)) <= float(band)]
    if not vals:
        return float("nan")
    arr = np.asarray(vals, dtype=float)
    return float(np.sqrt(np.mean(arr**2)))


def compute_structured_robustness(
    restart_payload: dict[str, Any],
    *,
    target: TargetUnitary,
    subspace: Subspace,
) -> dict[str, Any]:
    payload: dict[str, Any] = {}
    for label in ("L1d", "L2d"):
        if label not in restart_payload:
            continue
        best = max(restart_payload[label]["restarts"], key=lambda row: row["process_fidelity"])
        serialized = best["sequence"]
        chi_eps = np.linspace(-0.10, 0.10, 21)
        amp_eps = np.linspace(-0.10, 0.10, 21)
        dur_eps = np.linspace(-0.15, 0.15, 21)
        phi_eps = np.linspace(-0.20, 0.20, 21)
        chi_sweep = {
            float(eps): evaluate_serialized_sequence(
                serialized,
                target=target,
                subspace=subspace,
                chi=CHI * (1.0 + float(eps)),
            )
            for eps in chi_eps
        }
        amp_sweep = {
            float(eps): evaluate_serialized_sequence(
                serialized,
                target=target,
                subspace=subspace,
                amp_scale=1.0 + float(eps),
            )
            for eps in amp_eps
        }
        dur_sweep = {
            float(eps): evaluate_serialized_sequence(
                serialized,
                target=target,
                subspace=subspace,
                dur_scale=1.0 + float(eps),
            )
            for eps in dur_eps
        }
        phi_sweep = {
            float(eps): evaluate_serialized_sequence(
                serialized,
                target=target,
                subspace=subspace,
                phase_offset=float(eps),
            )
            for eps in phi_eps
        }
        payload[label] = {
            "nominal": chi_sweep[0.0],
            "chi_sweep": chi_sweep,
            "amplitude_sweep": amp_sweep,
            "duration_sweep": dur_sweep,
            "phase_sweep": phi_sweep,
            "sensitivity_summary": {
                "chi_rms_drop_pm5pct": _rms_drop(chi_sweep, band=0.05),
                "amplitude_rms_drop_pm5pct": _rms_drop(amp_sweep, band=0.05),
                "duration_rms_drop_pm5pct": _rms_drop(dur_sweep, band=0.05),
                "phase_rms_drop_pm0p1rad": _rms_drop(phi_sweep, band=0.10),
            },
        }
    return payload


def make_structured_robustness_figure(payload: dict[str, Any]) -> None:
    labels = [label for label in ("L1d", "L2d") if label in payload]
    fig, axes = plt.subplots(1, len(labels), figsize=(5.5 * max(1, len(labels)), 4.2))
    axes = np.atleast_1d(axes).ravel()
    colors = {
        "chi_sweep": "#1f77b4",
        "amplitude_sweep": "#d62728",
        "duration_sweep": "#2ca02c",
        "phase_sweep": "#9467bd",
    }
    for ax, label in zip(axes, labels, strict=True):
        record = payload[label]
        for key, legend in (
            ("chi_sweep", r"$\chi$ scale"),
            ("amplitude_sweep", "amp scale"),
            ("duration_sweep", "duration scale"),
            ("phase_sweep", "phase offset"),
        ):
            xs = []
            ys = []
            for raw_x, row in sorted(record[key].items(), key=lambda item: float(item[0])):
                x = float(raw_x)
                if key == "phase_sweep":
                    xs.append(x)
                else:
                    xs.append(100.0 * x)
                ys.append(float(row["process_fidelity"]))
            ax.plot(xs, ys, linewidth=1.8, color=colors[key], label=legend)
        ax.set_title(f"{label} local sensitivity")
        ax.set_xlabel("Relative error (%) or phase offset (rad)")
        ax.set_ylabel("Process fidelity")
        ax.set_ylim(0.70, 1.01)
        ax.grid(True, alpha=0.25, linestyle=":")
        ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "structured_robustness.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / "structured_robustness.pdf", bbox_inches="tight")
    plt.close(fig)


def make_pulse_landscape_figure(snap_followup: dict[str, Any], constructive_followup: dict[str, Any]) -> None:
    snap_rows = [
        row
        for row in snap_followup["candidates"]
        if str(row["family"]) == "gaussian" and abs(float(row["amplitude_scale"]) - 1.0) < 1.0e-12
    ]
    snap_chi = sorted({float(row["chi_t"]) for row in snap_rows})
    snap_shape = sorted({float(row["shape_parameter"]) for row in snap_rows})
    snap_grid = np.full((len(snap_shape), len(snap_chi)), np.nan, dtype=float)
    for row in snap_rows:
        i = snap_shape.index(float(row["shape_parameter"]))
        j = snap_chi.index(float(row["chi_t"]))
        snap_grid[i, j] = float(row["noisy_avg_state_fidelity"])

    cons_rows = [
        row
        for row in constructive_followup["candidates"]
        if str(row["family"]) == "gaussian" and abs(float(row["amplitude_scale"]) - 1.0) < 1.0e-12
    ]
    cons_broad = sorted({float(row["broad_duration_ns"]) for row in cons_rows})
    cons_chi = sorted({float(row["chi_t"]) for row in cons_rows})
    cons_grid = np.full((len(cons_broad), len(cons_chi)), np.nan, dtype=float)
    for i, broad in enumerate(cons_broad):
        for j, chi_t in enumerate(cons_chi):
            matches = [
                float(row["noisy_snap_like_avg_state_fidelity"])
                for row in cons_rows
                if abs(float(row["broad_duration_ns"]) - broad) < 1.0e-12 and abs(float(row["chi_t"]) - chi_t) < 1.0e-12
            ]
            if matches:
                cons_grid[i, j] = max(matches)

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4))
    im = axes[0].imshow(
        snap_grid,
        origin="lower",
        aspect="auto",
        extent=[min(snap_chi), max(snap_chi), min(snap_shape), max(snap_shape)],
        cmap="viridis",
        vmin=0.90,
        vmax=1.00,
    )
    axes[0].set_title("Logical SNAP noisy-fidelity landscape")
    axes[0].set_xlabel(r"$\chi T/2\pi$")
    axes[0].set_ylabel("Shape parameter")
    fig.colorbar(im, ax=axes[0], fraction=0.046, pad=0.04, label="Noisy avg state fidelity")

    im = axes[1].imshow(
        cons_grid,
        origin="lower",
        aspect="auto",
        extent=[min(cons_chi), max(cons_chi), min(cons_broad), max(cons_broad)],
        cmap="magma",
        vmin=0.45,
        vmax=0.70,
    )
    axes[1].set_title("Constructive shortcut best-over-shape landscape")
    axes[1].set_xlabel(r"$\chi T/2\pi$")
    axes[1].set_ylabel("Broadband $\\pi$ duration (ns)")
    fig.colorbar(im, ax=axes[1], fraction=0.046, pad=0.04, label="Noisy avg state fidelity")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "pulse_landscapes.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / "pulse_landscapes.pdf", bbox_inches="tight")
    plt.close(fig)


def build_pulse_list_from_sequence(
    serialized_sequence: list[dict[str, Any]],
    *,
    model,
    frame: FrameSpec,
    sqr_sigma_fraction: float = 0.18,
) -> tuple[list[Any], dict[str, str], list[dict[str, Any]]]:
    pulses = []
    drive_ops: dict[str, str] = {}
    pulse_metadata: list[dict[str, Any]] = []
    t_cursor = 0.0
    sqr_config = {
        "duration_sqr_s": None,
        "sqr_sigma_fraction": sqr_sigma_fraction,
        "sqr_theta_cutoff": 1.0e-8,
        "use_rotating_frame": True,
        "fock_fqs_hz": None,
    }
    for index, gate in enumerate(serialized_sequence):
        gate_type = str(gate["type"])
        duration = float(gate["duration"])
        params = gate.get("parameters", [])
        if gate_type == "Displacement":
            io_gate = DisplacementGate(index=index, name=str(gate["name"]), re=float(params[0]), im=float(params[1]))
            built, mapping, meta = build_displacement_pulse(io_gate, {"duration_displacement_s": duration})
        elif gate_type == "QubitRotation":
            io_gate = RotationGate(index=index, name=str(gate["name"]), theta=float(params[0]), phi=float(params[1]))
            built, mapping, meta = build_rotation_pulse(
                io_gate,
                {"duration_rotation_s": duration, "rotation_sigma_fraction": 0.18},
            )
        elif gate_type == "SQR":
            half = len(params) // 2
            io_gate = SQRGate(
                index=index,
                name=str(gate["name"]),
                theta=tuple(float(x) for x in params[:half]),
                phi=tuple(float(x) for x in params[half:]),
            )
            local_cfg = dict(sqr_config)
            local_cfg["duration_sqr_s"] = duration
            built, mapping, meta = build_sqr_multitone_pulse(io_gate, model, local_cfg, frame=frame)
        else:
            continue
        shifted = []
        for pulse in built:
            shifted.append(
                pulse.__class__(
                    channel=pulse.channel,
                    t0=float(t_cursor + pulse.t0),
                    duration=float(pulse.duration),
                    envelope=pulse.envelope,
                    amp=float(pulse.amp),
                    phase=float(pulse.phase),
                    carrier=float(pulse.carrier),
                    drag=float(getattr(pulse, "drag", 0.0)),
                    label=pulse.label,
                )
            )
        pulses.extend(shifted)
        drive_ops.update(mapping)
        meta_row = dict(meta)
        meta_row["gate_name"] = str(gate["name"])
        meta_row["gate_type"] = gate_type
        meta_row["t0_s"] = float(t_cursor)
        pulse_metadata.append(meta_row)
        t_cursor += duration
    return pulses, drive_ops, pulse_metadata


def evaluate_snap_phase_grid(
    *,
    phase_angle: float,
    model,
    frame: FrameSpec,
    noise_spec,
    target_branch: int = 1,
) -> dict[str, Any]:
    families = ("gaussian", "flat_top_gaussian")
    chi_t_values = (1.3, 1.7, 2.2, 3.0)
    shape_values = {"gaussian": (0.16, 0.20, 0.24), "flat_top_gaussian": (0.12, 0.18, 0.24)}
    amp_values = (0.90, 1.00, 1.10)
    logical_n = 2
    target = snap_target_operator(logical_n, target_branch, phase_angle)
    indices = cavity_ground_indices(model, logical_n)
    rows: list[dict[str, Any]] = []
    for family in families:
        for chi_t in chi_t_values:
            duration_s = duration_from_chi_t(chi_t)
            for shape in shape_values[family]:
                for amp_scale in amp_values:
                    first = build_selective_qubit_pulse(
                        model,
                        frame,
                        family=family,
                        branch=target_branch,
                        theta=np.pi,
                        phi=0.0,
                        duration_s=duration_s,
                        shape_parameter=shape,
                        amplitude_scale=amp_scale,
                        t0=0.0,
                        label=f"snap_{family}_a",
                    )
                    second = build_selective_qubit_pulse(
                        model,
                        frame,
                        family=family,
                        branch=target_branch,
                        theta=np.pi,
                        phi=float(np.pi - phase_angle),
                        duration_s=duration_s,
                        shape_parameter=shape,
                        amplitude_scale=amp_scale,
                        t0=duration_s,
                        label=f"snap_{family}_b",
                    )
                    pulses = [first, second]
                    closed = build_session(
                        model,
                        frame,
                        pulses,
                        {"qubit": "qubit"},
                        total_duration_s=2.0 * duration_s,
                        noise=None,
                        dt=LIT_DT,
                    )
                    noisy = build_session(
                        model,
                        frame,
                        pulses,
                        {"qubit": "qubit"},
                        total_duration_s=2.0 * duration_s,
                        noise=noise_spec,
                        dt=LIT_DT,
                    )
                    probes = snap_probe_state_vectors(logical_n)
                    closed_avg, _ = average_target_state_fidelity(
                        closed,
                        probes,
                        target,
                        model=model,
                        indices=indices,
                    )
                    noisy_avg, noisy_details = average_target_state_fidelity(
                        noisy,
                        probes,
                        target,
                        model=model,
                        indices=indices,
                    )
                    rows.append(
                        {
                            "family": family,
                            "chi_t": float(chi_t),
                            "shape_parameter": float(shape),
                            "amplitude_scale": float(amp_scale),
                            "pi_pulse_duration_us": float(duration_s * 1.0e6),
                            "total_duration_us": float(2.0 * duration_s * 1.0e6),
                            "closed_avg_state_fidelity": float(closed_avg),
                            "noisy_avg_state_fidelity": float(noisy_avg),
                            "noisy_details": noisy_details,
                        }
                    )
    best = max(rows, key=lambda row: row["noisy_avg_state_fidelity"])
    return {"candidates": rows, "best": best}


def evaluate_constructive_variant(
    *,
    model,
    frame: FrameSpec,
    noise_spec,
    target_branch: int = 1,
) -> dict[str, Any]:
    families = ("gaussian", "flat_top_gaussian")
    chi_t_values = (1.0, 1.3, 1.7, 2.2)
    shape_values = {"gaussian": (0.16, 0.20, 0.24), "flat_top_gaussian": (0.12, 0.18, 0.24)}
    amp_values = (0.90, 1.00, 1.10)
    broad_durations = (40.0e-9, 60.0e-9, 80.0e-9)
    logical_n = 2
    snap_like_target = snap_target_operator(logical_n, target_branch, np.pi)
    cavity_indices = cavity_ground_indices(model, logical_n)
    rows: list[dict[str, Any]] = []
    for family in families:
        for chi_t in chi_t_values:
            duration_s = duration_from_chi_t(chi_t)
            for shape in shape_values[family]:
                for amp_scale in amp_values:
                    selective = build_selective_qubit_pulse(
                        model,
                        frame,
                        family=family,
                        branch=target_branch,
                        theta=np.pi,
                        phi=0.0,
                        duration_s=duration_s,
                        shape_parameter=shape,
                        amplitude_scale=amp_scale,
                        t0=0.0,
                        label=f"constructive_sel_{family}",
                    )
                    for broad_duration in broad_durations:
                        io_gate = RotationGate(
                            index=1,
                            name=f"constructive_broad_{int(round(broad_duration * 1.0e9))}",
                            theta=np.pi,
                            phi=np.pi / 2.0,
                        )
                        broad_pulses, mapping, _meta = build_rotation_pulse(
                            io_gate,
                            {"duration_rotation_s": broad_duration, "rotation_sigma_fraction": 0.18},
                        )
                        shifted = []
                        for pulse in broad_pulses:
                            shifted.append(
                                pulse.__class__(
                                    channel=pulse.channel,
                                    t0=float(duration_s + pulse.t0),
                                    duration=float(pulse.duration),
                                    envelope=pulse.envelope,
                                    amp=float(pulse.amp),
                                    phase=float(pulse.phase),
                                    carrier=float(pulse.carrier),
                                    drag=float(getattr(pulse, "drag", 0.0)),
                                    label=pulse.label,
                                )
                            )
                        pulses = [selective] + shifted
                        total_duration = duration_s + broad_duration
                        closed = build_session(
                            model,
                            frame,
                            pulses,
                            {"qubit": "qubit", **mapping},
                            total_duration_s=total_duration,
                            noise=None,
                            dt=LIT_DT,
                        )
                        noisy = build_session(
                            model,
                            frame,
                            pulses,
                            {"qubit": "qubit", **mapping},
                            total_duration_s=total_duration,
                            noise=noise_spec,
                            dt=LIT_DT,
                        )
                        probes = snap_probe_state_vectors(logical_n)
                        strict_closed, _ = average_target_state_fidelity(
                            closed,
                            probes,
                            snap_like_target,
                            model=model,
                            indices=cavity_indices,
                        )
                        strict_noisy, noisy_details = average_target_state_fidelity(
                            noisy,
                            probes,
                            snap_like_target,
                            model=model,
                            indices=cavity_indices,
                        )
                        rows.append(
                            {
                                "family": family,
                                "chi_t": float(chi_t),
                                "shape_parameter": float(shape),
                                "amplitude_scale": float(amp_scale),
                                "broad_duration_ns": float(broad_duration * 1.0e9),
                                "total_duration_us": float(total_duration * 1.0e6),
                                "closed_snap_like_avg_state_fidelity": float(strict_closed),
                                "noisy_snap_like_avg_state_fidelity": float(strict_noisy),
                                "noisy_details": noisy_details,
                            }
                        )
    best = max(rows, key=lambda row: row["noisy_snap_like_avg_state_fidelity"])
    return {"candidates": rows, "best": best}


def run_grape_multiseed_followup() -> dict[str, Any]:
    case = grape_cases()["E_local_320"]
    model = build_model()
    frame = build_frame(model)
    subspace = logical_subspace()
    problem = build_control_problem_from_model(
        model,
        frame=frame,
        time_grid=gate_common.PiecewiseConstantTimeGrid.uniform(
            steps=case.steps,
            dt_s=case.duration_s / case.steps,
        ),
        channel_specs=(
            gate_common.ModelControlChannelSpec(
                name="storage",
                target="storage",
                quadratures=("I", "Q"),
                amplitude_bounds=(-case.amp_bound_rad_s, case.amp_bound_rad_s),
                export_channel="storage",
            ),
            gate_common.ModelControlChannelSpec(
                name="qubit",
                target="qubit",
                quadratures=("I", "Q"),
                amplitude_bounds=(-case.amp_bound_rad_s, case.amp_bound_rad_s),
                export_channel="qubit",
            ),
        ),
        objectives=(
            OCUnitaryObjective(
                target_operator=target_matrix(case.target_key),
                subspace=subspace,
                ignore_global_phase=True,
                name=case.label,
            ),
        ),
        penalties=(OCLeakagePenalty(weight=0.02, subspace=subspace),),
    )
    runs = []
    for seed in (17, 23, 29):
        solver = GrapeSolver(GrapeConfig(maxiter=90, seed=seed, random_scale=0.40))
        result = solver.solve(problem)
        sub_op, leakage = replay_grape_operator(result=result, problem=problem, model=model, frame=frame)
        strict = float(subspace_unitary_fidelity(sub_op, target_matrix(case.target_key), gauge="global"))
        pulses, _drive_ops, pulse_meta = result.to_pulses()
        time_axis, waveform = sample_total_waveform(pulses)
        runs.append(
            {
                "seed": int(seed),
                "strict_fidelity": strict,
                "leakage_average": float(np.mean(leakage)),
                "objective": float(result.objective_value),
                "max_abs_command_rad_s": float(result.schedule.max_abs_amplitude()),
                "pulse_meta": pulse_meta,
                "time_s": time_axis.tolist(),
                "waveform_re": np.real(waveform).tolist(),
                "waveform_im": np.imag(waveform).tolist(),
            }
        )
    best = max(runs, key=lambda row: row["strict_fidelity"])
    return {"runs": runs, "best": best}


def make_constructive_comparison_figure(snap_best: dict[str, Any], constructive_best: dict[str, Any]) -> None:
    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    labels = ["Standard SNAP", "Selective pi + broad pi"]
    closed_vals = [
        float(snap_best["closed_avg_state_fidelity"]),
        float(constructive_best["closed_snap_like_avg_state_fidelity"]),
    ]
    noisy_vals = [
        float(snap_best["noisy_avg_state_fidelity"]),
        float(constructive_best["noisy_snap_like_avg_state_fidelity"]),
    ]
    x = np.arange(2, dtype=float)
    width = 0.36
    ax.bar(x - width / 2.0, closed_vals, width=width, color="#4c78a8", label="closed")
    ax.bar(x + width / 2.0, noisy_vals, width=width, color="#f58518", label="noisy")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(0.0, 1.02)
    ax.set_ylabel("Average state fidelity")
    ax.set_title("Constructive shortcut versus tuned logical-window SNAP")
    ax.legend(frameon=False)
    ax.grid(True, axis="y", alpha=0.25, linestyle=":")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "constructive_variant_comparison.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / "constructive_variant_comparison.pdf", bbox_inches="tight")
    plt.close(fig)


def make_waveform_figure(
    *,
    best_l1: dict[str, Any],
    snap_best: dict[str, Any],
    constructive_best: dict[str, Any],
    grape_best: dict[str, Any],
) -> None:
    model = build_model(n_tr=3)
    frame = build_frame(model)
    l1_sequence = best_l1["sequence"]
    _pulses, _drive_ops, metadata = build_pulse_list_from_sequence(l1_sequence, model=model, frame=frame)
    sqr_meta = [row for row in metadata if row["gate_type"] == "SQR"]
    l1_waveforms = []
    for row in sqr_meta[:3]:
        gate_name = row["gate_name"]
        gate = next(item for item in l1_sequence if item["name"] == gate_name)
        sub_pulses, _mapping, _meta = build_pulse_list_from_sequence([gate], model=model, frame=frame)
        t_s, w_s = sample_total_waveform(sub_pulses)
        l1_waveforms.append((gate_name, t_s, w_s, row))

    snap_duration = duration_from_chi_t(float(snap_best["chi_t"]))
    snap_a = build_selective_qubit_pulse(
        model,
        frame,
        family=str(snap_best["family"]),
        branch=1,
        theta=np.pi,
        phi=0.0,
        duration_s=snap_duration,
        shape_parameter=float(snap_best["shape_parameter"]),
        amplitude_scale=float(snap_best["amplitude_scale"]),
        t0=0.0,
        label="snap_a",
    )
    snap_b = build_selective_qubit_pulse(
        model,
        frame,
        family=str(snap_best["family"]),
        branch=1,
        theta=np.pi,
        phi=0.0,
        duration_s=snap_duration,
        shape_parameter=float(snap_best["shape_parameter"]),
        amplitude_scale=float(snap_best["amplitude_scale"]),
        t0=snap_duration,
        label="snap_b",
    )
    t_snap, w_snap = sample_total_waveform([snap_a, snap_b])
    f_snap, s_snap = sample_spectrum(t_snap, w_snap)

    broad_duration = float(constructive_best["broad_duration_ns"]) * 1.0e-9
    cons_selective = build_selective_qubit_pulse(
        model,
        frame,
        family=str(constructive_best["family"]),
        branch=1,
        theta=np.pi,
        phi=0.0,
        duration_s=duration_from_chi_t(float(constructive_best["chi_t"])),
        shape_parameter=float(constructive_best["shape_parameter"]),
        amplitude_scale=float(constructive_best["amplitude_scale"]),
        t0=0.0,
        label="constructive_sel",
    )
    io_gate = RotationGate(index=1, name="constructive_broad", theta=np.pi, phi=np.pi / 2.0)
    broad_pulses, _, _ = build_rotation_pulse(io_gate, {"duration_rotation_s": broad_duration, "rotation_sigma_fraction": 0.18})
    cons_broad = []
    for pulse in broad_pulses:
        cons_broad.append(
            pulse.__class__(
                channel=pulse.channel,
                t0=float(cons_selective.duration + pulse.t0),
                duration=float(pulse.duration),
                envelope=pulse.envelope,
                amp=float(pulse.amp),
                phase=float(pulse.phase),
                carrier=float(pulse.carrier),
                drag=float(getattr(pulse, "drag", 0.0)),
                label=pulse.label,
            )
        )
    t_cons, w_cons = sample_total_waveform([cons_selective] + cons_broad)

    t_grape = np.asarray(grape_best["time_s"], dtype=float)
    w_grape = np.asarray(grape_best["waveform_re"], dtype=float) + 1j * np.asarray(grape_best["waveform_im"], dtype=float)

    fig, axes = plt.subplots(3, 2, figsize=(11.6, 10.8))
    ax = axes[0, 0]
    ax.plot(t_snap * 1.0e9, np.real(w_snap), color="#1f77b4", label="I")
    ax.plot(t_snap * 1.0e9, np.imag(w_snap), color="#d62728", label="Q")
    ax.set_title("Logical-window SNAP realization")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Drive amplitude")
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.25, linestyle=":")

    ax = axes[0, 1]
    mask = f_snap < 30.0e6
    ax.plot(f_snap[mask] / 1.0e6, s_snap[mask], color="#4c78a8")
    ax.set_title("SNAP spectrum")
    ax.set_xlabel("Frequency (MHz)")
    ax.set_ylabel("|FFT|")
    ax.grid(True, alpha=0.25, linestyle=":")

    ax = axes[1, 0]
    for gate_name, t_s, w_s, _row in l1_waveforms:
        ax.plot(t_s * 1.0e6, np.real(w_s), linewidth=1.3, label=f"{gate_name} I")
    ax.set_title("Optimized L1d SQR waveforms")
    ax.set_xlabel("Time (us)")
    ax.set_ylabel("Drive amplitude")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(True, alpha=0.25, linestyle=":")

    ax = axes[1, 1]
    tone_lines = []
    for row in sqr_meta[:1]:
        for tone in row.get("active_tones", []):
            tone_lines.append((int(tone.get("manifold", tone.get("n", -1))), float(tone["amp_rad_s"]) / TWO_PI / 1.0e6))
    tone_lines.sort(key=lambda item: item[0])
    if tone_lines:
        ax.bar([f"n={n}" for n, _ in tone_lines], [amp for _, amp in tone_lines], color="#f58518", alpha=0.85)
        ax.set_ylabel("Tone amplitude / 2pi (MHz)")
        ax.set_title("Representative multitone SQR decomposition (S1)")
        ax.tick_params(axis="x", rotation=20)
    ax.grid(True, axis="y", alpha=0.25, linestyle=":")

    ax = axes[2, 0]
    ax.plot(t_cons * 1.0e9, np.real(w_cons), color="#2ca02c", label="I")
    ax.plot(t_cons * 1.0e9, np.imag(w_cons), color="#9467bd", label="Q")
    ax.set_title("Constructive selective-pi + broadband-pi variant")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Drive amplitude")
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.25, linestyle=":")

    ax = axes[2, 1]
    ax.plot(t_grape * 1.0e9, np.real(w_grape), color="#8c564b", label="I")
    ax.plot(t_grape * 1.0e9, np.imag(w_grape), color="#17becf", label="Q")
    ax.set_title("GRAPE lower-bound waveform (best seed)")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Drive amplitude")
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.25, linestyle=":")

    fig.tight_layout()
    fig.savefig(FIG_DIR / "waveform_realizations.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIG_DIR / "waveform_realizations.pdf", bbox_inches="tight")
    plt.close(fig)


def write_generated_tables(
    *,
    restart_payload: dict[str, Any],
    robustness_payload: dict[str, Any],
    snap_best: dict[str, Any],
    constructive_best: dict[str, Any],
    grape_best: dict[str, Any],
) -> None:
    best_a = max(restart_payload["A_local"]["restarts"], key=lambda row: row["process_fidelity"])
    best_d = max(restart_payload["D_ent"]["restarts"], key=lambda row: row["process_fidelity"])
    best_l1 = max(restart_payload["L1d"]["restarts"], key=lambda row: row["process_fidelity"])
    best_l2 = max(restart_payload["L2d"]["restarts"], key=lambda row: row["process_fidelity"])
    model = build_model(n_tr=3)
    frame = build_frame(model)
    _pulses, _mapping, best_l1_meta = build_pulse_list_from_sequence(best_l1["sequence"], model=model, frame=frame)
    best_l1_s1 = next((row for row in best_l1_meta if row["gate_name"] == "S1"), {})

    def _gate_rows(sequence: list[dict[str, Any]]) -> list[str]:
        rows = []
        for gate in sequence:
            gate_type = str(gate["type"])
            duration_ns = float(gate["duration"]) * 1.0e9
            params = gate.get("parameters", [])
            if gate_type == "Displacement":
                desc = f"$\\alpha={params[0]:+.4f}{params[1]:+.4f}i$"
            elif gate_type == "QubitRotation":
                desc = f"$\\theta={params[0]:.4f}$, $\\phi={params[1]:.4f}$"
            elif gate_type == "SQR":
                half = len(params) // 2
                active = [idx for idx, val in enumerate(params[:half]) if abs(float(val)) > 1.0e-3]
                desc = f"{len(active)} active tones; $n={','.join(str(idx) for idx in active[:6])}$"
            elif gate_type == "ConditionalPhaseSQR":
                active = [idx for idx, val in enumerate(params) if abs(float(val)) > 1.0e-3]
                desc = f"{len(active)} active phases; $n={','.join(str(idx) for idx in active[:6])}$"
            elif gate_type == "FreeEvolveCondPhase":
                desc = "native dispersive wait"
            else:
                desc = gate_type
            rows.append(f"{gate['name']} & {gate_type} & {duration_ns:.1f} & {desc} \\\\")
        return rows

    tex = []
    tex.append("% Auto-generated by run_followup_optimization.py")
    tex.append("\\begin{table}[H]")
    tex.append("\\centering")
    tex.append("\\caption{Actual optimized follow-up candidates. Sequence-level metrics are quoted on the logical target subspace used by each optimization, while the GRAPE row remains a pulse-level lower bound.}")
    tex.append("\\begin{tabular}{@{} l l c c c c c @{} }")
    tex.append("\\toprule")
    tex.append("Candidate & Level & Fidelity & Leakage & Duration & Params & Multistart best/median/worst \\\\")
    tex.append("\\midrule")
    for label, level in [
        ("A_local", "ideal D-SNAP-D"),
        ("D_ent", "ideal R-W-R"),
        ("L1c", "ideal D-R-SQR"),
        ("L1d", "ideal D-R-SQR (+D)"),
        ("L2c", "ideal D-R-CP"),
        ("L2d", "ideal D-R-CP (+D)"),
    ]:
        summary = restart_payload[label]["summary"]
        label_tex = label.replace("_", "\\_")
        tex.append(
            f"{label_tex} & {level} & {summary['best_process_fidelity']:.4f} & "
            f"{summary['best_leakage_average']:.4f} & {summary['best_duration_ns']:.1f} ns & "
            f"{summary['parameter_count']} & {summary['best_process_fidelity']:.4f}/"
            f"{summary['median_process_fidelity']:.4f}/{summary['worst_process_fidelity']:.4f} \\\\"
        )
    tex.append(
        f"logical SNAP & pulse-backed & {snap_best['closed_avg_state_fidelity']:.4f} & -- & "
        f"{snap_best['total_duration_us']:.3f} $\\mu$s & 4 & -- \\\\"
    )
    tex.append(
        f"constructive & pulse-backed & {constructive_best['closed_snap_like_avg_state_fidelity']:.4f} & -- & "
        f"{constructive_best['total_duration_us']:.3f} $\\mu$s & 5 & -- \\\\"
    )
    tex.append(
        f"E\\_local\\_320 & GRAPE lower bound & {grape_best['strict_fidelity']:.4f} & "
        f"{grape_best['leakage_average']:.4f} & 320.0 ns & 64 & -- \\\\"
    )
    tex.append("\\bottomrule")
    tex.append("\\end{tabular}")
    tex.append("\\end{table}")
    tex.append("")
    tex.append("\\begin{table}[H]")
    tex.append("\\centering")
    tex.append("\\caption{Actual optimized gate list for the best local and native-entangler solutions.}")
    tex.append("\\begin{tabular}{@{} l l c l @{} }")
    tex.append("\\toprule")
    tex.append("Gate & Type & Duration (ns) & Optimized parameters \\\\")
    tex.append("\\midrule")
    tex.extend(_gate_rows(best_a["sequence"]))
    tex.append("\\midrule")
    tex.extend(_gate_rows(best_d["sequence"]))
    tex.append("\\bottomrule")
    tex.append("\\end{tabular}")
    tex.append("\\end{table}")
    tex.append("")
    tex.append("\\begin{table}[H]")
    tex.append("\\centering")
    tex.append("\\caption{Actual optimized gate list for the best structured $U_{\\mathrm{target}}$ solutions.}")
    tex.append("\\begin{tabular}{@{} l l c l @{} }")
    tex.append("\\toprule")
    tex.append("Gate & Type & Duration (ns) & Parameter summary \\\\")
    tex.append("\\midrule")
    tex.extend(_gate_rows(best_l1["sequence"]))
    tex.append("\\midrule")
    tex.extend(_gate_rows(best_l2["sequence"]))
    tex.append("\\bottomrule")
    tex.append("\\end{tabular}")
    tex.append("\\end{table}")
    tex.append("")
    tex.append("\\begin{table}[H]")
    tex.append("\\centering")
    tex.append("\\caption{Representative multitone decomposition for the first SQR gate in the best L1d sequence.}")
    tex.append("\\begin{tabular}{@{} c c c c @{} }")
    tex.append("\\toprule")
    tex.append("Fock branch $n$ & Detuning$/2\\pi$ (MHz) & Tone amp$/2\\pi$ (MHz) & Phase (rad) \\\\")
    tex.append("\\midrule")
    for tone in best_l1_s1.get("active_tones", []):
        tex.append(
            f"{int(tone.get('manifold', tone.get('n', -1)))} & "
            f"{float(tone['omega_rad_s']) / TWO_PI / 1.0e6:.4f} & "
            f"{float(tone['amp_rad_s']) / TWO_PI / 1.0e6:.4f} & "
            f"{float(tone['phase_rad']):.4f} \\\\"
        )
    tex.append("\\bottomrule")
    tex.append("\\end{tabular}")
    tex.append("\\end{table}")
    tex.append("")
    tex.append("\\begin{table}[H]")
    tex.append("\\centering")
    tex.append("\\caption{Local perturbation robustness for the best structured $U_{\\mathrm{target}}$ candidates.}")
    tex.append("\\begin{tabular}{@{} l c c c c @{} }")
    tex.append("\\toprule")
    tex.append("Candidate & RMS drop ($\\chi$, $\\pm 5\\%$) & RMS drop (amp, $\\pm 5\\%$) & RMS drop (dur, $\\pm 5\\%$) & RMS drop (phase, $\\pm 0.1$ rad) \\\\")
    tex.append("\\midrule")
    for label in ("L1d", "L2d"):
        summary = robustness_payload[label]["sensitivity_summary"]
        tex.append(
            f"{label} & {summary['chi_rms_drop_pm5pct']:.4f} & {summary['amplitude_rms_drop_pm5pct']:.4f} & "
            f"{summary['duration_rms_drop_pm5pct']:.4f} & {summary['phase_rms_drop_pm0p1rad']:.4f} \\\\"
        )
    tex.append("\\bottomrule")
    tex.append("\\end{tabular}")
    tex.append("\\end{table}")
    tex.append("")
    tex.append("\\begin{table}[H]")
    tex.append("\\centering")
    tex.append("\\caption{Pulse-backed selective control settings used in the follow-up waveform section.}")
    tex.append("\\begin{tabular}{@{} l c c c c c @{} }")
    tex.append("\\toprule")
    tex.append("Primitive & Family & $\\chi T/2\\pi$ & Shape & Amp scale & Duration \\\\")
    tex.append("\\midrule")
    tex.append(
        f"logical SNAP & {snap_best['family'].replace('_', ' ')} & {snap_best['chi_t']:.1f} & "
        f"{snap_best['shape_parameter']:.2f} & {snap_best['amplitude_scale']:.2f} & "
        f"{snap_best['total_duration_us']:.3f} $\\mu$s \\\\"
    )
    tex.append(
        f"constructive selective stage & {constructive_best['family'].replace('_', ' ')} & {constructive_best['chi_t']:.1f} & "
        f"{constructive_best['shape_parameter']:.2f} & {constructive_best['amplitude_scale']:.2f} & "
        f"{constructive_best['total_duration_us']:.3f} $\\mu$s total \\\\"
    )
    tex.append("\\bottomrule")
    tex.append("\\end{tabular}")
    tex.append("\\end{table}")
    (DATA_DIR / "generated_tables.tex").write_text("\n".join(tex), encoding="utf-8")


def main() -> None:
    t0 = time.time()
    subspace_gate = logical_subspace()
    target_utarget = build_utarget_target()
    subspace_utarget = build_utarget_subspace()
    restart_payload: dict[str, Any]
    structured_path = DATA_DIR / "structured_restarts.json"
    snap_path = DATA_DIR / "logical_snap_followup.json"
    constructive_path = DATA_DIR / "constructive_variant.json"
    grape_path = DATA_DIR / "grape_followup.json"
    robustness_path = DATA_DIR / "structured_robustness.json"

    local_variants = {
        "A_local": GateSequence(gates=library_a_local_sequence(), n_cav=N_CAV_DEFAULT),
        "D_ent": GateSequence(gates=library_a_entangler_sequence(), n_cav=N_CAV_DEFAULT),
    }
    utarget_variants = build_utarget_sequences()
    sequence_runs = {
        "A_local": {"sequence": local_variants["A_local"], "target_key": "local_h", "subspace": subspace_gate, "target": None, "maxiter": 90, "model": build_model()},
        "D_ent": {"sequence": local_variants["D_ent"], "target_key": "cx_c_to_q", "subspace": subspace_gate, "target": None, "maxiter": 70, "model": build_model()},
        "L1c": {"sequence": utarget_variants["L1c"], "target_key": None, "subspace": subspace_utarget, "target": target_utarget, "maxiter": 90, "model": None},
        "L1d": {"sequence": utarget_variants["L1d"], "target_key": None, "subspace": subspace_utarget, "target": target_utarget, "maxiter": 80, "model": None},
        "L2c": {"sequence": utarget_variants["L2c"], "target_key": None, "subspace": subspace_utarget, "target": target_utarget, "maxiter": 90, "model": None},
        "L2d": {"sequence": utarget_variants["L2d"], "target_key": None, "subspace": subspace_utarget, "target": target_utarget, "maxiter": 80, "model": None},
    }

    if structured_path.exists():
        restart_payload = json.loads(structured_path.read_text(encoding="utf-8"))
    else:
        restart_payload = {}
        for label, spec in sequence_runs.items():
            batch: list[RestartResult] = []
            print(f"[followup] running {label} restart sweep")
            for seed in FOLLOWUP_SEEDS[label]:
                batch.append(
                    run_single_restart(
                        label=label,
                        sequence=spec["sequence"],
                        target_key=spec["target_key"],
                        target=spec["target"],
                        subspace=spec["subspace"],
                        seed=seed,
                        maxiter=spec["maxiter"],
                        model=spec["model"],
                    )
                )
            if label in {"A_local", "L1c", "L2c"}:
                best = max(batch, key=lambda row: row.process_fidelity)
                batch.append(
                    local_refine(
                        label=label,
                        seed=best.seed + 77,
                        best_restart=best,
                        target=spec["target"] if spec["target"] is not None else _synthesis_target(spec["target_key"]),
                        subspace=spec["subspace"],
                        model=spec["model"],
                        maxiter=120,
                        optimizer="nelder_mead",
                    )
                )
            restart_payload[label] = {
                "summary": summarize_restart_batch(batch),
                "restarts": batch_to_payload(batch),
            }
            structured_path.write_text(json.dumps(_json_ready(restart_payload), indent=2), encoding="utf-8")

    print("[followup] selective pulse follow-up")
    pulse_model = build_model(n_tr=3)
    pulse_frame = build_frame(pulse_model)
    noise_spec = build_noise_spec()
    if snap_path.exists():
        snap_followup = json.loads(snap_path.read_text(encoding="utf-8"))
    else:
        snap_followup = evaluate_snap_phase_grid(
            phase_angle=np.pi,
            model=pulse_model,
            frame=pulse_frame,
            noise_spec=noise_spec,
        )
        snap_path.write_text(json.dumps(_json_ready(snap_followup), indent=2), encoding="utf-8")
    if constructive_path.exists():
        constructive_followup = json.loads(constructive_path.read_text(encoding="utf-8"))
    else:
        constructive_followup = evaluate_constructive_variant(
            model=pulse_model,
            frame=pulse_frame,
            noise_spec=noise_spec,
        )
        constructive_path.write_text(json.dumps(_json_ready(constructive_followup), indent=2), encoding="utf-8")

    print("[followup] GRAPE lower-bound rerun")
    if grape_path.exists():
        grape_followup = json.loads(grape_path.read_text(encoding="utf-8"))
    else:
        grape_followup = run_grape_multiseed_followup()
        grape_path.write_text(json.dumps(_json_ready(grape_followup), indent=2), encoding="utf-8")
    existing_phase4 = json.loads((STUDY_ROOT / "data" / "utarget_decomposition" / "phase4_results.json").read_text(encoding="utf-8"))
    if robustness_path.exists():
        structured_robustness = json.loads(robustness_path.read_text(encoding="utf-8"))
    else:
        structured_robustness = compute_structured_robustness(
            restart_payload,
            target=target_utarget,
            subspace=subspace_utarget,
        )
        robustness_path.write_text(json.dumps(_json_ready(structured_robustness), indent=2), encoding="utf-8")

    make_restart_figure(restart_payload)
    make_convergence_figure(restart_payload)
    make_ansatz_scaling_figure(existing_phase4, restart_payload)
    make_constructive_comparison_figure(snap_followup["best"], constructive_followup["best"])
    make_structured_robustness_figure(structured_robustness)
    make_pulse_landscape_figure(snap_followup, constructive_followup)
    make_waveform_figure(
        best_l1=max(restart_payload["L1d"]["restarts"], key=lambda row: row["process_fidelity"]),
        snap_best=snap_followup["best"],
        constructive_best=constructive_followup["best"],
        grape_best=grape_followup["best"],
    )
    write_generated_tables(
        restart_payload=restart_payload,
        robustness_payload=structured_robustness,
        snap_best=snap_followup["best"],
        constructive_best=constructive_followup["best"],
        grape_best=grape_followup["best"],
    )

    validation = {
        "figures_exist": all(
            (FIG_DIR / name).exists()
            for name in (
                "multistart_statistics.pdf",
                "convergence_histories.pdf",
                "ansatz_scaling.pdf",
                "constructive_variant_comparison.pdf",
                "structured_robustness.pdf",
                "pulse_landscapes.pdf",
                "waveform_realizations.pdf",
            )
        ),
        "generated_tables_exist": (DATA_DIR / "generated_tables.tex").exists(),
        "restart_batches": {
            key: {
                "best_process_fidelity": float(value["summary"]["best_process_fidelity"]),
                "median_process_fidelity": float(value["summary"]["median_process_fidelity"]),
                "runs": len(value["restarts"]),
            }
            for key, value in restart_payload.items()
        },
        "snap_best_noisy_avg_state_fidelity": float(snap_followup["best"]["noisy_avg_state_fidelity"]),
        "constructive_best_noisy_snap_like_avg_state_fidelity": float(
            constructive_followup["best"]["noisy_snap_like_avg_state_fidelity"]
        ),
        "grape_best_strict_fidelity": float(grape_followup["best"]["strict_fidelity"]),
        "structured_robustness": _json_ready(structured_robustness),
    }
    payload = {
        "metadata": {
            "runtime_s": float(time.time() - t0),
            "study_root": str(STUDY_ROOT),
        },
        "structured_restarts": _json_ready(restart_payload),
        "logical_snap_followup": _json_ready(snap_followup),
        "constructive_variant": _json_ready(constructive_followup),
        "grape_followup": _json_ready(grape_followup),
        "validation": _json_ready(validation),
    }
    (DATA_DIR / "followup_results.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    (DATA_DIR / "validation_summary.json").write_text(json.dumps(_json_ready(validation), indent=2), encoding="utf-8")
    print(f"Saved follow-up outputs to {DATA_DIR}")


if __name__ == "__main__":
    main()
