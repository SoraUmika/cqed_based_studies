from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import Any, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import runtime_compat  # noqa: F401
import qutip as qt


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import common as c


STYLE_PATH = SCRIPT_DIR.parents[2] / ".github" / "skills" / "publication-figures" / "assets" / "cqed_style.mplstyle"
if STYLE_PATH.exists():
    plt.style.use(str(STYLE_PATH))


OPTIMIZATION_N_CAV = 8
VALIDATION_N_CAV = 12
DEFAULT_MAXITER = 200
DEFAULT_SEEDS = (17, 42, 73)
WIGNER_XVEC = np.linspace(-3.0, 3.0, 151)

JSON_PATH = c.DATA_DIR / "primitive_grape_diagnostics.json"
ARRAY_PATH = c.ARTIFACT_DIR / "primitive_grape_diagnostics_arrays.npz"
SUMMARY_FIG_PATH = c.FIG_DIR / "primitive_grape_summary"
WIGNER_FIG_PATH = c.FIG_DIR / "primitive_grape_wigner"

SOURCE_ARTIFACTS = {
    "sqr": c.ARTIFACT_DIR / "corrected_best_sqr.json",
    "cpsqr": c.ARTIFACT_DIR / "corrected_best_cpsqr.json",
}


def full_subspace(*, n_cav: int) -> c.Subspace:
    indices = tuple(range(2 * int(n_cav)))
    labels = tuple(
        [f"|g,{level}>" for level in range(int(n_cav))]
        + [f"|e,{level}>" for level in range(int(n_cav))]
    )
    return c.Subspace.custom(2 * int(n_cav), indices, labels)


def selected_level_subspace(levels: Sequence[int], *, n_cav: int) -> c.Subspace:
    ordered_levels = tuple(int(level) for level in levels)
    indices = tuple(list(ordered_levels) + [int(n_cav) + int(level) for level in ordered_levels])
    labels = tuple(
        [f"|g,{level}>" for level in ordered_levels]
        + [f"|e,{level}>" for level in ordered_levels]
    )
    return c.Subspace.custom(2 * int(n_cav), indices, labels)


def full_basis_states(model: c.DispersiveTransmonCavityModel) -> list[tuple[str, qt.Qobj]]:
    rows: list[tuple[str, qt.Qobj]] = []
    for level in range(int(model.n_cav)):
        rows.append((f"|g,{level}>", model.basis_state(0, level)))
    for level in range(int(model.n_cav)):
        rows.append((f"|e,{level}>", model.basis_state(1, level)))
    return rows


def selected_basis_states(
    model: c.DispersiveTransmonCavityModel,
    *,
    levels: Sequence[int],
) -> list[tuple[str, qt.Qobj]]:
    rows: list[tuple[str, qt.Qobj]] = []
    for level in levels:
        rows.append((f"|g,{int(level)}>", model.basis_state(0, int(level))))
    for level in levels:
        rows.append((f"|e,{int(level)}>", model.basis_state(1, int(level))))
    return rows


def reduced_cavity_density(state: qt.Qobj, *, n_cav: int) -> qt.Qobj:
    if state.isket:
        ket = state
        if ket.dims != [[2, int(n_cav)], [1]] and ket.dims != [[2, int(n_cav)], [1, 1]]:
            ket = qt.Qobj(ket.full(), dims=[[2, int(n_cav)], [1, 1]])
        rho = qt.ket2dm(ket)
    else:
        rho = state
        if rho.dims != [[2, int(n_cav)], [2, int(n_cav)]]:
            rho = qt.Qobj(rho.full(), dims=[[2, int(n_cav)], [2, int(n_cav)]])
    rho = rho / rho.tr()
    return rho.ptrace(1)


def normalized_state(state: qt.Qobj, *, n_cav: int) -> qt.Qobj:
    if state.isket:
        vec = np.asarray(state.full(), dtype=np.complex128).reshape(-1)
        norm = np.linalg.norm(vec)
        if norm > 0.0:
            vec = vec / norm
        return qt.Qobj(vec, dims=[[2, int(n_cav)], [1, 1]])
    rho = state / state.tr()
    return qt.Qobj(rho.full(), dims=[[2, int(n_cav)], [2, int(n_cav)]])


def pure_target_fidelity(target_state: qt.Qobj, actual_state: qt.Qobj) -> float:
    target = normalized_state(target_state, n_cav=target_state.dims[0][1] if len(target_state.dims[0]) > 1 else int(target_state.shape[0] // 2))
    if actual_state.isket:
        actual = normalized_state(actual_state, n_cav=target.dims[0][1])
        overlap = np.vdot(
            np.asarray(target.full(), dtype=np.complex128).reshape(-1),
            np.asarray(actual.full(), dtype=np.complex128).reshape(-1),
        )
        return float(abs(overlap) ** 2)
    rho = normalized_state(actual_state, n_cav=target.dims[0][1])
    target_dm = qt.ket2dm(target)
    value = (target.dag() * rho * target).full()[0, 0]
    return float(np.real(value))


def apply_ideal_operator(target_operator: np.ndarray, state: qt.Qobj, *, n_cav: int) -> qt.Qobj:
    vector = np.asarray(state.full(), dtype=np.complex128).reshape(-1)
    final = np.asarray(target_operator, dtype=np.complex128) @ vector
    return qt.Qobj(final, dims=[[2, int(n_cav)], [1, 1]])


def plus_cavity_superposition_state(
    model: c.DispersiveTransmonCavityModel,
    levels: Sequence[int],
    *,
    coeffs: Sequence[complex] | None = None,
) -> qt.Qobj:
    n_cav = int(model.n_cav)
    level_list = [int(level) for level in levels]
    if coeffs is None:
        coeff_list = [1.0] * len(level_list)
    else:
        coeff_list = [complex(value) for value in coeffs]
    vector = np.zeros(2 * n_cav, dtype=np.complex128)
    for level, coeff in zip(level_list, coeff_list):
        vector[int(level)] += coeff / np.sqrt(2.0)
        vector[n_cav + int(level)] += coeff / np.sqrt(2.0)
    norm = np.linalg.norm(vector)
    if norm <= 0.0:
        raise ValueError("Probe-state norm is zero.")
    vector = vector / norm
    return qt.Qobj(vector, dims=[[2, n_cav], [1, 1]])


def short_probe_label(levels: Sequence[int]) -> str:
    joined = ",".join(str(int(level)) for level in levels)
    return f"+ : ({joined})"


def sqr_probe_specs() -> list[dict[str, Any]]:
    return [
        {"levels": (0, 5), "coeffs": None, "plot_label": short_probe_label((0, 5)), "detail_label": "|+>x(|0>+|5>)/sqrt(2)"},
        {"levels": (1, 3), "coeffs": None, "plot_label": short_probe_label((1, 3)), "detail_label": "|+>x(|1>+|3>)/sqrt(2)"},
        {"levels": (0, 1, 2, 3), "coeffs": None, "plot_label": short_probe_label((0, 1, 2, 3)), "detail_label": "|+>x(|0>+|1>+|2>+|3>)/2"},
        {"levels": (2, 5), "coeffs": None, "plot_label": short_probe_label((2, 5)), "detail_label": "|+>x(|2>+|5>)/sqrt(2)"},
    ]


def cpsqr_probe_specs() -> list[dict[str, Any]]:
    return [
        {"levels": (1, 2), "coeffs": None, "plot_label": short_probe_label((1, 2)), "detail_label": "|+>x(|1>+|2>)/sqrt(2)"},
        {"levels": (0, 1), "coeffs": None, "plot_label": short_probe_label((0, 1)), "detail_label": "|+>x(|0>+|1>)/sqrt(2)"},
        {"levels": (2, 3), "coeffs": None, "plot_label": short_probe_label((2, 3)), "detail_label": "|+>x(|2>+|3>)/sqrt(2)"},
        {"levels": (0, 1, 2, 3), "coeffs": None, "plot_label": short_probe_label((0, 1, 2, 3)), "detail_label": "|+>x(|0>+|1>+|2>+|3>)/2"},
    ]


def diagnostic_levels(gate_key: str) -> tuple[int, ...]:
    specs = sqr_probe_specs() if gate_key == "sqr" else cpsqr_probe_specs()
    return tuple(sorted({int(level) for spec in specs for level in spec["levels"]}))


def load_record(path: Path) -> dict[str, Any]:
    payload = c.load_json(path)
    return dict(payload["record"])


def select_primitive_gate(record: dict[str, Any], *, family_key: str) -> tuple[dict[str, Any], str]:
    desired_kind = "MaskedSQR" if family_key == "sqr" else "MaskedCPSQR"
    candidates = [
        dict(gate)
        for gate in record["sequence"]
        if gate.get("type") == "PrimitiveGate"
        and isinstance(gate.get("metadata"), dict)
        and gate["metadata"].get("ideal_kind") == desired_kind
    ]
    if not candidates:
        raise RuntimeError(f"No {desired_kind} primitive found in {record['case_id']}.")
    if family_key == "sqr":
        return candidates[0], "first selective block from retained SQR winner"
    return max(candidates, key=lambda gate: float(gate["duration"])), "longest selective block from retained CPSQR winner"


def target_matrix_from_gate(gate: dict[str, Any], *, n_cav: int) -> np.ndarray:
    metadata = gate.get("metadata", {}) if isinstance(gate.get("metadata"), dict) else {}
    levels = tuple(int(level) for level in metadata.get("levels", []))
    params = [float(value) for value in gate.get("parameters", [])]
    duration = float(gate["duration"])
    kind = str(metadata.get("ideal_kind", "PrimitiveGate"))
    if kind == "MaskedSQR":
        tone_count = len(levels)
        theta = params[:tone_count]
        phi = params[tone_count:]
        return c._masked_sqr_matrix(
            theta=theta,
            phi=phi,
            levels=levels,
            n_cav=int(n_cav),
            duration=duration,
            include_conditional_phase=bool(metadata.get("include_conditional_phase", False)),
            drift_model=c.IDEAL_DRIFT,
        )
    if kind == "MaskedCPSQR":
        return c._masked_cpsqr_matrix(
            phases=params,
            levels=levels,
            n_cav=int(n_cav),
            duration=duration,
            include_drift=bool(metadata.get("include_drift", True)),
            drift_model=c.PHYSICAL_DRIFT,
        )
    raise ValueError(f"Unsupported primitive kind {kind}.")


def restrict_operator_to_levels(operator: np.ndarray, *, levels: Sequence[int], n_cav: int) -> np.ndarray:
    ordered_levels = [int(level) for level in levels]
    indices = ordered_levels + [int(n_cav) + int(level) for level in ordered_levels]
    full = np.asarray(operator, dtype=np.complex128)
    return np.asarray(full[np.ix_(indices, indices)], dtype=np.complex128)


def build_primitive_grape_problem(
    *,
    model: c.DispersiveTransmonCavityModel,
    target_operator: np.ndarray,
    duration_s: float,
    amp_bound_rad_s: float,
    name: str,
    subspace: c.Subspace,
) -> Any:
    frame = c.build_frame(model)
    steps = max(10, round(float(duration_s) / c.GRAPE_DT_S))
    return c.build_control_problem_from_model(
        model,
        frame=frame,
        time_grid=c.PiecewiseConstantTimeGrid.uniform(steps=int(steps), dt_s=float(duration_s) / int(steps)),
        channel_specs=(
            c.ModelControlChannelSpec(
                name="qubit",
                target="qubit",
                quadratures=("I", "Q"),
                amplitude_bounds=(-float(amp_bound_rad_s), float(amp_bound_rad_s)),
                export_channel="qubit",
            ),
        ),
        objectives=(
            c.OCUnitaryObjective(
                target_operator=np.asarray(target_operator, dtype=np.complex128),
                subspace=subspace,
                ignore_global_phase=True,
                name=str(name),
            ),
        ),
    )


def channel_payload(compiled: Any, *, channel_name: str) -> dict[str, Any]:
    channel = compiled.channels.get(str(channel_name)) if hasattr(compiled, "channels") else None
    if channel is None:
        return {
            "time_ns": np.asarray([], dtype=float),
            "baseband": np.asarray([], dtype=np.complex128),
            "distorted": np.asarray([], dtype=np.complex128),
        }
    tlist = np.asarray(compiled.tlist, dtype=float)
    baseband = np.asarray(channel.baseband, dtype=np.complex128)
    distorted = np.asarray(getattr(channel, "distorted", channel.baseband), dtype=np.complex128)
    count = min(int(tlist.size), int(baseband.size), int(distorted.size))
    return {
        "time_ns": 1.0e9 * tlist[:count],
        "baseband": baseband[:count],
        "distorted": distorted[:count],
    }


def analyze_probe_states(
    *,
    gate_key: str,
    compiled: Any,
    drive_ops: dict[str, Any],
    target_operator_full: np.ndarray,
    model: c.DispersiveTransmonCavityModel,
) -> dict[str, Any]:
    specs = sqr_probe_specs() if gate_key == "sqr" else cpsqr_probe_specs()
    basis_states = [
        (
            spec["detail_label"],
            plus_cavity_superposition_state(model, spec["levels"], coeffs=spec.get("coeffs")),
        )
        for spec in specs
    ]
    rows = c.replay_compiled_sequence(
        model=model,
        compiled=compiled,
        drive_ops=drive_ops,
        basis_states=basis_states,
        noise=None,
        store_states=False,
    )
    metrics: list[dict[str, Any]] = []
    worst_entry: dict[str, Any] | None = None
    for spec, row in zip(specs, rows):
        initial_state = basis_states[specs.index(spec)][1]
        ideal_state = normalized_state(
            apply_ideal_operator(target_operator_full, initial_state, n_cav=int(model.n_cav)),
            n_cav=int(model.n_cav),
        )
        actual_state = normalized_state(row["simulation"].final_state, n_cav=int(model.n_cav))
        full_state_fidelity = pure_target_fidelity(ideal_state, actual_state)
        ideal_rho = reduced_cavity_density(ideal_state, n_cav=int(model.n_cav))
        actual_rho = reduced_cavity_density(actual_state, n_cav=int(model.n_cav))
        reduced_state_fidelity = float(qt.fidelity(ideal_rho, actual_rho))
        target_wigner = qt.wigner(ideal_rho, WIGNER_XVEC, WIGNER_XVEC)
        actual_wigner = qt.wigner(actual_rho, WIGNER_XVEC, WIGNER_XVEC)
        delta_wigner = actual_wigner - target_wigner
        wigner_rms = float(np.sqrt(np.mean(np.square(delta_wigner))))
        entry = {
            "plot_label": spec["plot_label"],
            "detail_label": spec["detail_label"],
            "levels": [int(level) for level in spec["levels"]],
            "full_state_fidelity": float(full_state_fidelity),
            "reduced_cavity_fidelity": float(reduced_state_fidelity),
            "wigner_rms": float(wigner_rms),
        }
        if worst_entry is None or (
            entry["reduced_cavity_fidelity"],
            entry["full_state_fidelity"],
        ) < (
            worst_entry["reduced_cavity_fidelity"],
            worst_entry["full_state_fidelity"],
        ):
            worst_entry = {
                **entry,
                "target_wigner": target_wigner,
                "actual_wigner": actual_wigner,
                "delta_wigner": delta_wigner,
            }
        metrics.append(entry)
    if worst_entry is None:
        raise RuntimeError("No primitive probe metrics were generated.")
    return {
        "probe_metrics": metrics,
        "mean_full_state_fidelity": float(np.mean([entry["full_state_fidelity"] for entry in metrics])),
        "min_full_state_fidelity": float(np.min([entry["full_state_fidelity"] for entry in metrics])),
        "mean_reduced_cavity_fidelity": float(np.mean([entry["reduced_cavity_fidelity"] for entry in metrics])),
        "min_reduced_cavity_fidelity": float(np.min([entry["reduced_cavity_fidelity"] for entry in metrics])),
        "worst_probe": worst_entry,
    }


def run_gate_diagnostic(
    *,
    gate_key: str,
    record: dict[str, Any],
    gate_entry: dict[str, Any],
    source_selection: str,
    seeds: Sequence[int],
    maxiter: int,
    engine: str,
    jax_device: str | None,
    optimization_n_cav: int,
    validation_n_cav: int,
    amp_bound_rad_s: float,
) -> dict[str, Any]:
    opt_model = c.build_model(n_cav=int(optimization_n_cav))
    val_model = c.build_model(n_cav=int(validation_n_cav))
    kind = str(gate_entry["metadata"]["ideal_kind"])
    levels_for_diagnostics = diagnostic_levels(gate_key)
    target_opt_full = target_matrix_from_gate(gate_entry, n_cav=int(optimization_n_cav))
    target_val_full = target_matrix_from_gate(gate_entry, n_cav=int(validation_n_cav))
    target_opt = restrict_operator_to_levels(target_opt_full, levels=levels_for_diagnostics, n_cav=int(optimization_n_cav))
    target_val = restrict_operator_to_levels(target_val_full, levels=levels_for_diagnostics, n_cav=int(validation_n_cav))
    opt_subspace = selected_level_subspace(levels_for_diagnostics, n_cav=int(optimization_n_cav))
    val_subspace = selected_level_subspace(levels_for_diagnostics, n_cav=int(validation_n_cav))
    problem = build_primitive_grape_problem(
        model=opt_model,
        target_operator=target_opt,
        duration_s=float(gate_entry["duration"]),
        amp_bound_rad_s=float(amp_bound_rad_s),
        name=f"{gate_key}_{gate_entry['name']}",
        subspace=opt_subspace,
    )
    val_basis = selected_basis_states(val_model, levels=levels_for_diagnostics)
    seed_rows: list[dict[str, Any]] = []
    best_row: dict[str, Any] | None = None
    for seed in seeds:
        print(f"[primitive-grape] {gate_key}: seed={int(seed)} start", flush=True)
        started = time.perf_counter()
        result = c.run_grape_seed(problem, seed=int(seed), maxiter=int(maxiter), engine=engine, jax_device=jax_device)
        replay = c.replay_grape_subspace(
            result=result,
            problem=problem,
            model=val_model,
            subspace=val_subspace,
            target_unitary=target_val,
            basis_states=val_basis,
            noise=None,
            store_states=False,
        )
        nominal_fidelity = float(getattr(result, "metrics", {}).get("nominal_fidelity", getattr(result, "metrics", {}).get("fidelity", np.nan)))
        validation_subspace_fidelity = float(replay["fidelity"])
        elapsed_s = float(time.perf_counter() - started)
        row = {
            "seed": int(seed),
            "engine": c.resolve_grape_engine(engine),
            "success": bool(getattr(result, "success", True)),
            "message": str(getattr(result, "message", "")),
            "iterations": int(getattr(result, "nit", 0)),
            "nominal_fidelity": nominal_fidelity,
            "validation_subspace_fidelity": validation_subspace_fidelity,
            "elapsed_s": elapsed_s,
        }
        seed_rows.append(row)
        print(
            f"[primitive-grape] {gate_key}: seed={int(seed)} nominal={nominal_fidelity:.6f} validation_subspace={validation_subspace_fidelity:.6f} elapsed={elapsed_s:.1f}s",
            flush=True,
        )
        if best_row is None or row["validation_subspace_fidelity"] > best_row["validation_subspace_fidelity"]:
            best_row = {**row, "result": result, "replay": replay}
    if best_row is None:
        raise RuntimeError(f"No GRAPE results were produced for {gate_key}.")
    replay = best_row["replay"]
    waveform = channel_payload(replay["compiled"], channel_name="qubit")
    peak_mhz = float(np.max(np.abs(waveform["distorted"])) / (2.0 * np.pi * 1.0e6)) if waveform["distorted"].size else 0.0
    probe_payload = analyze_probe_states(
        gate_key=gate_key,
        compiled=replay["compiled"],
        drive_ops=replay["drive_ops"],
        target_operator_full=target_val_full,
        model=val_model,
    )
    metadata = gate_entry.get("metadata", {}) if isinstance(gate_entry.get("metadata"), dict) else {}
    parameter_values = [float(value) for value in gate_entry.get("parameters", [])]
    return {
        "family_key": str(record["family_key"]),
        "record_case_id": str(record["case_id"]),
        "source_artifact": SOURCE_ARTIFACTS[gate_key].name,
        "source_selection": source_selection,
        "gate_name": str(gate_entry["name"]),
        "ideal_kind": kind,
        "levels": [int(level) for level in metadata.get("levels", [])],
        "diagnostic_levels": [int(level) for level in levels_for_diagnostics],
        "duration_ns": float(gate_entry["duration"]) * 1.0e9,
        "parameters": parameter_values,
        "metadata": metadata,
        "optimization_n_cav": int(optimization_n_cav),
        "validation_n_cav": int(validation_n_cav),
        "best_seed": int(best_row["seed"]),
        "best_engine": str(best_row["engine"]),
        "best_nominal_fidelity": float(best_row["nominal_fidelity"]),
        "best_validation_subspace_fidelity": float(best_row["validation_subspace_fidelity"]),
        "best_iterations": int(best_row["iterations"]),
        "peak_qubit_drive_mhz": peak_mhz,
        "seed_runs": seed_rows,
        "waveform": {
            "time_ns": waveform["time_ns"],
            "baseband": waveform["baseband"],
            "distorted": waveform["distorted"],
        },
        **probe_payload,
    }


def plot_summary(results: dict[str, dict[str, Any]]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 7.8), constrained_layout=True)
    order = [("sqr", "Representative SQR primitive"), ("cpsqr", "Representative CPSQR primitive")]
    for row_index, (gate_key, title) in enumerate(order):
        payload = results[gate_key]
        waveform = payload["waveform"]
        wave_axis = axes[row_index, 0]
        amp_scale = 2.0 * np.pi * 1.0e6
        distorted_mhz = np.asarray(waveform["distorted"], dtype=complex) / amp_scale
        wave_axis.plot(waveform["time_ns"], np.real(distorted_mhz), linewidth=1.4, color="#4477AA", label="I")
        wave_axis.plot(waveform["time_ns"], np.imag(distorted_mhz), linewidth=1.4, color="#EE6677", label="Q")
        wave_axis.axhline(c.GRAPE_AMP_BOUND / amp_scale, color="0.4", linewidth=1.0, linestyle="--")
        wave_axis.axhline(-c.GRAPE_AMP_BOUND / amp_scale, color="0.4", linewidth=1.0, linestyle="--")
        wave_axis.set_ylabel("Drive amplitude (MHz)")
        wave_axis.set_xlabel("Time (ns)")
        wave_axis.set_title(
            f"{title}: {payload['gate_name']}\n"
            f"F_sub={payload['best_validation_subspace_fidelity']:.4f}, peak={payload['peak_qubit_drive_mhz']:.1f} MHz"
        )
        if row_index == 0:
            wave_axis.legend(frameon=False, ncol=2, loc="upper right")

        fidelity_axis = axes[row_index, 1]
        probe_rows = payload["probe_metrics"]
        x = np.arange(len(probe_rows), dtype=float)
        full_values = [float(entry["full_state_fidelity"]) for entry in probe_rows]
        reduced_values = [float(entry["reduced_cavity_fidelity"]) for entry in probe_rows]
        fidelity_axis.bar(x - 0.18, full_values, width=0.36, color="#4477AA", label="Full-state fidelity")
        fidelity_axis.bar(x + 0.18, reduced_values, width=0.36, color="#228833", label="Reduced cavity fidelity")
        fidelity_axis.set_xticks(x)
        fidelity_axis.set_xticklabels([entry["plot_label"] for entry in probe_rows], rotation=15)
        fidelity_axis.set_ylim(0.0, 1.005)
        fidelity_axis.set_ylabel("Fidelity")
        fidelity_axis.set_title(
            f"{title}: probe-state diagnostics\n"
            f"min(full)={payload['min_full_state_fidelity']:.4f}, min(cav)={payload['min_reduced_cavity_fidelity']:.4f}"
        )
        if row_index == 0:
            fidelity_axis.legend(frameon=False, loc="lower right")

    fig.savefig(SUMMARY_FIG_PATH.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(SUMMARY_FIG_PATH.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_wigner(results: dict[str, dict[str, Any]]) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(10.2, 6.8), constrained_layout=True)
    row_titles = {
        "sqr": "Representative SQR primitive",
        "cpsqr": "Representative CPSQR primitive",
    }
    mesh = None
    vmax = 0.35
    delta_max = 0.12
    for row_index, gate_key in enumerate(("sqr", "cpsqr")):
        worst = results[gate_key]["worst_probe"]
        panel_specs = [
            ("Ideal cavity state", worst["target_wigner"], "RdBu_r", -vmax, vmax),
            ("GRAPE replay cavity state", worst["actual_wigner"], "RdBu_r", -vmax, vmax),
            ("Difference", worst["delta_wigner"], "PuOr_r", -delta_max, delta_max),
        ]
        for col_index, (title, data, cmap, vmin, vmax_local) in enumerate(panel_specs):
            axis = axes[row_index, col_index]
            mesh = axis.pcolormesh(WIGNER_XVEC, WIGNER_XVEC, np.asarray(data), shading="auto", cmap=cmap, vmin=vmin, vmax=vmax_local, rasterized=True)
            axis.set_aspect("equal")
            axis.set_title(title)
            if col_index == 0:
                axis.set_ylabel(
                    f"{row_titles[gate_key]}\n{worst['plot_label']}\nIm(alpha)"
                )
            if row_index == 1:
                axis.set_xlabel("Re(alpha)")
        axes[row_index, 1].text(
            0.02,
            0.03,
            (
                f"probe={worst['detail_label']}\n"
                f"full fidelity={worst['full_state_fidelity']:.4f}\n"
                f"cavity fidelity={worst['reduced_cavity_fidelity']:.4f}\n"
                f"Wigner RMS={worst['wigner_rms']:.4e}"
            ),
            transform=axes[row_index, 1].transAxes,
            fontsize=8.5,
            bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "alpha": 0.85, "edgecolor": "0.8"},
        )

    fig.savefig(WIGNER_FIG_PATH.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(WIGNER_FIG_PATH.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_artifacts(results: dict[str, dict[str, Any]], *, args: argparse.Namespace) -> None:
    payload = {
        "study_name": "cluster_state_holographic_unified",
        "date_created": time.strftime("%Y-%m-%d"),
        "description": "Representative primitive-level GRAPE diagnostics for the retained SQR and CPSQR winners.",
        "parameters": {
            "optimization_n_cav": int(args.optimization_n_cav),
            "validation_n_cav": int(args.validation_n_cav),
            "maxiter": int(args.maxiter),
            "seeds": [int(seed) for seed in args.seeds],
            "engine_requested": str(args.engine),
            "engine_resolved": c.resolve_grape_engine(args.engine),
            "jax_device": None if args.jax_device is None else str(args.jax_device),
            "amp_bound_mhz": float(args.amp_bound_mhz),
        },
        "load_instructions": (
            "import json, numpy as np; from pathlib import Path; summary = json.loads(Path('data/primitive_grape_diagnostics.json').read_text(encoding='utf-8')); arrays = np.load('artifacts/primitive_grape_diagnostics_arrays.npz')"
        ),
        "results": c.json_ready(results),
    }
    c.save_json(JSON_PATH, payload)

    np.savez(
        ARRAY_PATH,
        wigner_xvec=WIGNER_XVEC,
        sqr_waveform_time_ns=np.asarray(results["sqr"]["waveform"]["time_ns"], dtype=float),
        sqr_waveform_real=np.real(np.asarray(results["sqr"]["waveform"]["distorted"], dtype=np.complex128)),
        sqr_waveform_imag=np.imag(np.asarray(results["sqr"]["waveform"]["distorted"], dtype=np.complex128)),
        cpsqr_waveform_time_ns=np.asarray(results["cpsqr"]["waveform"]["time_ns"], dtype=float),
        cpsqr_waveform_real=np.real(np.asarray(results["cpsqr"]["waveform"]["distorted"], dtype=np.complex128)),
        cpsqr_waveform_imag=np.imag(np.asarray(results["cpsqr"]["waveform"]["distorted"], dtype=np.complex128)),
        sqr_target_wigner=np.asarray(results["sqr"]["worst_probe"]["target_wigner"], dtype=float),
        sqr_actual_wigner=np.asarray(results["sqr"]["worst_probe"]["actual_wigner"], dtype=float),
        sqr_delta_wigner=np.asarray(results["sqr"]["worst_probe"]["delta_wigner"], dtype=float),
        cpsqr_target_wigner=np.asarray(results["cpsqr"]["worst_probe"]["target_wigner"], dtype=float),
        cpsqr_actual_wigner=np.asarray(results["cpsqr"]["worst_probe"]["actual_wigner"], dtype=float),
        cpsqr_delta_wigner=np.asarray(results["cpsqr"]["worst_probe"]["delta_wigner"], dtype=float),
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Optimize representative SQR/CPSQR primitives with GRAPE and generate cavity-state diagnostics.")
    parser.add_argument("--optimization-n-cav", type=int, default=OPTIMIZATION_N_CAV)
    parser.add_argument("--validation-n-cav", type=int, default=VALIDATION_N_CAV)
    parser.add_argument("--maxiter", type=int, default=DEFAULT_MAXITER)
    parser.add_argument("--engine", type=str, default="auto")
    parser.add_argument("--jax-device", type=str, default=None)
    parser.add_argument("--amp-bound-mhz", type=float, default=float(c.GRAPE_AMP_BOUND / (2.0 * np.pi * 1.0e6)))
    parser.add_argument("--seeds", type=int, nargs="+", default=list(DEFAULT_SEEDS))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    resolved_engine = c.resolve_grape_engine(args.engine)
    print(
        (
            f"[primitive-grape] optimization_n_cav={int(args.optimization_n_cav)} "
            f"validation_n_cav={int(args.validation_n_cav)} maxiter={int(args.maxiter)} "
            f"engine={resolved_engine} seeds={list(int(seed) for seed in args.seeds)}"
        ),
        flush=True,
    )

    sqr_record = load_record(SOURCE_ARTIFACTS["sqr"])
    cpsqr_record = load_record(SOURCE_ARTIFACTS["cpsqr"])
    sqr_gate, sqr_selection = select_primitive_gate(sqr_record, family_key="sqr")
    cpsqr_gate, cpsqr_selection = select_primitive_gate(cpsqr_record, family_key="cpsqr")

    amp_bound_rad_s = 2.0 * np.pi * float(args.amp_bound_mhz) * 1.0e6
    results = {
        "sqr": run_gate_diagnostic(
            gate_key="sqr",
            record=sqr_record,
            gate_entry=sqr_gate,
            source_selection=sqr_selection,
            seeds=args.seeds,
            maxiter=args.maxiter,
            engine=args.engine,
            jax_device=args.jax_device,
            optimization_n_cav=args.optimization_n_cav,
            validation_n_cav=args.validation_n_cav,
            amp_bound_rad_s=amp_bound_rad_s,
        ),
        "cpsqr": run_gate_diagnostic(
            gate_key="cpsqr",
            record=cpsqr_record,
            gate_entry=cpsqr_gate,
            source_selection=cpsqr_selection,
            seeds=args.seeds,
            maxiter=args.maxiter,
            engine=args.engine,
            jax_device=args.jax_device,
            optimization_n_cav=args.optimization_n_cav,
            validation_n_cav=args.validation_n_cav,
            amp_bound_rad_s=amp_bound_rad_s,
        ),
    }
    plot_summary(results)
    plot_wigner(results)
    save_artifacts(results, args=args)
    print(f"[primitive-grape] wrote {JSON_PATH.relative_to(c.STUDY_ROOT)}", flush=True)
    print(f"[primitive-grape] wrote {ARRAY_PATH.relative_to(c.STUDY_ROOT)}", flush=True)
    print(f"[primitive-grape] wrote {SUMMARY_FIG_PATH.with_suffix('.pdf').relative_to(c.STUDY_ROOT)}", flush=True)
    print(f"[primitive-grape] wrote {WIGNER_FIG_PATH.with_suffix('.pdf').relative_to(c.STUDY_ROOT)}", flush=True)


if __name__ == "__main__":
    main()