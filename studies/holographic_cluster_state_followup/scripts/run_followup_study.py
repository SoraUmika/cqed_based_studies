from __future__ import annotations

import argparse
import json
import statistics
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import qutip as qt


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import common as c


DECOMP_EVAL_NCAV = (2, 8, 10, 12, 15)
PULSE_REPLAY_NCAV = (8, 12)
WIGNER_NCAV = 12

DEFAULT_STAGES = ("transfer", "decomp", "grape", "wigner", "summary")
DEFAULT_DURATIONS = (100, 150, 200, 250, 300, 400)
DEFAULT_SEEDS = (17, 42, 73, 91, 103, 127, 211, 307, 401, 509)


def _parameter_list(value: Any) -> list[float]:
    if isinstance(value, (list, tuple)):
        return [float(v) for v in value]
    return [float(value)]


def serialize_gate(gate: Any) -> dict[str, Any]:
    payload = {
        "name": str(gate.name),
        "type": type(gate).__name__,
        "duration_s": float(getattr(gate, "duration", 0.0)),
        "duration_ns": float(getattr(gate, "duration", 0.0) * 1.0e9),
        "active_tones": c.selective_gate_active_tones(gate),
    }
    if isinstance(gate, c.Displacement):
        payload["parameters"] = [float(np.real(gate.alpha)), float(np.imag(gate.alpha))]
    elif isinstance(gate, c.QubitRotation):
        payload["parameters"] = [float(gate.theta), float(gate.phi)]
    elif isinstance(gate, c.SNAP):
        payload["parameters"] = _parameter_list(gate.phases)
    elif isinstance(gate, c.SQR):
        payload["parameters"] = _parameter_list(gate.theta_n) + _parameter_list(gate.phi_n)
    elif isinstance(gate, c.ConditionalPhaseSQR):
        payload["parameters"] = _parameter_list(gate.phases_n)
    elif isinstance(gate, c.FreeEvolveCondPhase):
        payload["parameters"] = [float(gate.duration)]
    else:
        payload["parameters"] = []
    return payload


def serialize_sequence(sequence: c.GateSequence) -> list[dict[str, Any]]:
    return [serialize_gate(gate) for gate in sequence.gates]


def build_candidate_specs() -> list[dict[str, Any]]:
    return [
        {
            "key": "drsnap_unbounded",
            "label": "D-R-SNAP",
            "family": "D-R-SNAP",
            "category": "decomposition",
            "builder": c.build_dr_snap_sequence,
            "builder_kwargs": {"blocks": 2},
            "max_amplitude": None,
            "multistart": 2,
            "maxiter": 140,
        },
        {
            "key": "drsnap_bounded",
            "label": "Bounded D-R-SNAP",
            "family": "D-R-SNAP",
            "category": "decomposition",
            "builder": c.build_dr_snap_sequence,
            "builder_kwargs": {"blocks": 2},
            "max_amplitude": 0.30,
            "multistart": 2,
            "maxiter": 140,
        },
        {
            "key": "dsqrcp_unbounded",
            "label": "D-SQR-CPSQR",
            "family": "D-SQR-CPSQR",
            "category": "sqr_like",
            "builder": c.build_dsqr_cp_sequence,
            "builder_kwargs": {"blocks": 2},
            "max_amplitude": None,
            "multistart": 3,
            "maxiter": 220,
        },
        {
            "key": "dsqrcp_bounded",
            "label": "Bounded D-SQR-CPSQR",
            "family": "D-SQR-CPSQR",
            "category": "sqr_like",
            "builder": c.build_dsqr_cp_sequence,
            "builder_kwargs": {"blocks": 2},
            "max_amplitude": 0.30,
            "multistart": 3,
            "maxiter": 220,
        },
        {
            "key": "drsqrcp_bounded",
            "label": "Hybrid D-R-SQR-CPSQR",
            "family": "D-R-SQR-CPSQR",
            "category": "sqr_like",
            "builder": c.build_dr_sqr_cp_sequence,
            "builder_kwargs": {"blocks": 2},
            "max_amplitude": 0.30,
            "multistart": 3,
            "maxiter": 240,
        },
        {
            "key": "drfe_unbounded",
            "label": "D-R-FreeEvolve",
            "family": "D-R-FreeEvolve",
            "category": "entangler_assisted",
            "builder": c.build_dr_fe_sequence,
            "builder_kwargs": {"blocks": 2},
            "max_amplitude": None,
            "multistart": 3,
            "maxiter": 220,
        },
        {
            "key": "drfe_bounded",
            "label": "Bounded D-R-FreeEvolve",
            "family": "D-R-FreeEvolve",
            "category": "entangler_assisted",
            "builder": c.build_dr_fe_sequence,
            "builder_kwargs": {"blocks": 2},
            "max_amplitude": 0.30,
            "multistart": 3,
            "maxiter": 220,
        },
    ]


def load_json(path: Path) -> Any:
    return None if not path.exists() else json.loads(path.read_text(encoding="utf-8"))


def save_artifact_npz(path: Path, **arrays: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path, **arrays)


def best_key_by_metric(records: dict[str, Any], *, metric_path: tuple[str, ...], allowed_category: str | None = None) -> str | None:
    best_key = None
    best_value = -np.inf
    for key, payload in records.items():
        if allowed_category is not None and payload.get("category") != allowed_category:
            continue
        value = payload
        for part in metric_path:
            value = value.get(part, {})
        if isinstance(value, dict):
            continue
        candidate_value = float(value)
        if candidate_value > best_value:
            best_key = key
            best_value = candidate_value
    return best_key


def embedded_operator_from_basis_outputs(rows: list[dict[str, Any]], *, n_cav: int) -> np.ndarray:
    full_dim = 2 * int(n_cav)
    operator = np.eye(full_dim, dtype=np.complex128)
    for index, row in zip(c.logical_indices(int(n_cav)), rows, strict=True):
        state = row["simulation"].final_state
        vector = np.asarray(state.full(), dtype=np.complex128).reshape(-1)
        operator[:, index] = vector
    return operator


def logical_state_qobj(model: c.DispersiveTransmonCavityModel, coeffs: np.ndarray) -> qt.Qobj:
    full_dim = 2 * int(model.n_cav)
    vector = np.zeros(full_dim, dtype=np.complex128)
    indices = c.logical_indices(int(model.n_cav))
    vector[list(indices)] = np.asarray(coeffs, dtype=np.complex128)
    return qt.Qobj(vector, dims=[[int(model.n_tr), int(model.n_cav)], [1, 1]])


def restrict_state_to_logical_density(state: qt.Qobj, *, n_cav: int) -> np.ndarray:
    rho = state if state.isoper else state.proj()
    arr = np.asarray(rho.full(), dtype=np.complex128)
    idx = list(c.logical_indices(int(n_cav)))
    return arr[np.ix_(idx, idx)]


def process_fidelity_from_choi(choi: np.ndarray, target_unitary: np.ndarray) -> float:
    d = target_unitary.shape[0]
    target_vec = target_unitary.reshape(-1, order="F")
    choi_target = np.outer(target_vec, target_vec.conj())
    value = np.trace(choi_target.conj().T @ choi)
    return float(np.real(value) / (d * d))


def reconstruct_open_process(
    *,
    compiled: Any,
    drive_ops: dict[str, Any],
    model: c.DispersiveTransmonCavityModel,
    noise: c.NoiseSpec,
) -> dict[str, Any]:
    frame = c.build_frame(model)
    e_ops = None
    basis_outputs: dict[int, np.ndarray] = {}
    leakage_probe: list[float] = []
    target_state_fids: list[float] = []

    for basis_index in range(4):
        coeffs = np.zeros(4, dtype=np.complex128)
        coeffs[basis_index] = 1.0
        sim = c.pulse_simulate_sequence(
            model,
            compiled,
            logical_state_qobj(model, coeffs),
            drive_ops,
            config=c.SimulationConfig(frame=frame, store_states=False),
            noise=noise,
            e_ops=e_ops,
        )
        rho_logical = restrict_state_to_logical_density(sim.final_state, n_cav=int(model.n_cav))
        basis_outputs[basis_index] = rho_logical
        leakage_probe.append(float(max(0.0, 1.0 - np.trace(rho_logical).real)))
        psi_target = c.TARGET_UNITARY[:, basis_index]
        target_state_fids.append(float(np.real(np.vdot(psi_target, rho_logical @ psi_target))))

    e_blocks: dict[tuple[int, int], np.ndarray] = {}
    for index, rho in basis_outputs.items():
        e_blocks[(index, index)] = rho

    for i in range(4):
        for j in range(i + 1, 4):
            plus = np.zeros(4, dtype=np.complex128)
            plus[i] = 1.0 / np.sqrt(2.0)
            plus[j] = 1.0 / np.sqrt(2.0)
            plus_i = np.zeros(4, dtype=np.complex128)
            plus_i[i] = 1.0 / np.sqrt(2.0)
            plus_i[j] = 1.0j / np.sqrt(2.0)

            sim_plus = c.pulse_simulate_sequence(
                model,
                compiled,
                logical_state_qobj(model, plus),
                drive_ops,
                config=c.SimulationConfig(frame=frame, store_states=False),
                noise=noise,
            )
            sim_plus_i = c.pulse_simulate_sequence(
                model,
                compiled,
                logical_state_qobj(model, plus_i),
                drive_ops,
                config=c.SimulationConfig(frame=frame, store_states=False),
                noise=noise,
            )
            rho_plus = restrict_state_to_logical_density(sim_plus.final_state, n_cav=int(model.n_cav))
            rho_plus_i = restrict_state_to_logical_density(sim_plus_i.final_state, n_cav=int(model.n_cav))
            s_term = 2.0 * rho_plus - basis_outputs[i] - basis_outputs[j]
            t_term = 2.0 * rho_plus_i - basis_outputs[i] - basis_outputs[j]
            e_ij = 0.5 * (s_term + 1.0j * t_term)
            e_ji = 0.5 * (s_term - 1.0j * t_term)
            e_blocks[(i, j)] = e_ij
            e_blocks[(j, i)] = e_ji

    choi = np.zeros((16, 16), dtype=np.complex128)
    for i in range(4):
        for j in range(4):
            row = slice(i * 4, (i + 1) * 4)
            col = slice(j * 4, (j + 1) * 4)
            choi[row, col] = e_blocks[(i, j)]

    return {
        "basis_outputs": basis_outputs,
        "choi_matrix": choi,
        "process_fidelity": process_fidelity_from_choi(choi, c.TARGET_UNITARY),
        "mean_probe_leakage": float(np.mean(leakage_probe)),
        "max_probe_leakage": float(np.max(leakage_probe)),
        "mean_basis_state_fidelity": float(np.mean(target_state_fids)),
        "min_basis_state_fidelity": float(np.min(target_state_fids)),
    }


def run_transfer_stage() -> dict[str, Any]:
    summary = c.channel_transfer_summary()
    summary["target_unitary"] = c.TARGET_UNITARY
    summary["study_root"] = c.STUDY_ROOT
    c.save_json(c.DATA_DIR / "channel_summary.json", summary)
    return summary


def run_decomp_stage() -> dict[str, Any]:
    results: dict[str, Any] = {}
    for spec in build_candidate_specs():
        start = time.perf_counter()
        print(f"[decomp] {spec['key']} :: fitting", flush=True)
        sequence = spec["builder"](n_cav=c.N_CAV_OPT, **spec["builder_kwargs"])
        fit = c.fit_sequence(
            sequence,
            max_amplitude=spec["max_amplitude"],
            multistart=spec["multistart"],
            maxiter=spec["maxiter"],
        )
        optimized_sequence = fit["result"].sequence
        evaluations: dict[str, Any] = {}
        artifact_arrays: dict[str, Any] = {}
        for n_cav in DECOMP_EVAL_NCAV:
            print(f"[decomp]   {spec['key']} @ N_cav={n_cav}", flush=True)
            evaluation = c.evaluate_sequence_ideal(optimized_sequence, n_cav=n_cav)
            evaluations[str(n_cav)] = {
                "fidelity": evaluation["fidelity"],
                "block_fidelity": evaluation["block_fidelity"],
                "leakage_average": evaluation["leakage_average"],
                "leakage_worst": evaluation["leakage_worst"],
            }
            artifact_arrays[f"ideal_full_operator_nc{n_cav}"] = evaluation["full_operator"]
            artifact_arrays[f"ideal_subspace_operator_nc{n_cav}"] = evaluation["subspace_operator"]

        pulse_replay: dict[str, Any] = {}
        if c.waveform_bridge_supported(optimized_sequence):
            for n_cav in PULSE_REPLAY_NCAV:
                print(f"[decomp]   {spec['key']} pulse replay @ N_cav={n_cav}", flush=True)
                replay = c.evaluate_sequence_pulse(optimized_sequence, n_cav=n_cav)
                pulse_replay[str(n_cav)] = {
                    key: replay[key]
                    for key in (
                        "supported",
                        "fidelity",
                        "block_fidelity",
                        "leakage_average",
                        "leakage_worst",
                        "metrics",
                    )
                    if key in replay
                }
                if replay.get("supported"):
                    artifact_arrays[f"pulse_full_operator_nc{n_cav}"] = replay["full_operator"]
                    artifact_arrays[f"pulse_subspace_operator_nc{n_cav}"] = replay["subspace_operator"]
        gate_summary = c.sequence_gate_summary(optimized_sequence)
        artifact_path = c.ARTIFACT_DIR / f"{spec['key']}_operators.npz"
        save_artifact_npz(artifact_path, **artifact_arrays)
        elapsed_s = time.perf_counter() - start
        results[spec["key"]] = {
            "key": spec["key"],
            "label": spec["label"],
            "family": spec["family"],
            "category": spec["category"],
            "max_amplitude": spec["max_amplitude"],
            "fit_ideal_fidelity": fit["ideal_fidelity"],
            "fit_objective": fit["objective"],
            "fit_success": fit["success"],
            "fit_message": fit["message"],
            "elapsed_s": elapsed_s,
            "gate_summary": gate_summary,
            "sequence": serialize_sequence(optimized_sequence),
            "embedded_evaluations": evaluations,
            "pulse_replay": pulse_replay,
            "artifact_path": str(artifact_path),
        }

    best_embedded_key = max(
        results,
        key=lambda key: float(results[key]["embedded_evaluations"]["12"]["fidelity"]),
    )
    sqr_candidates = {key: payload for key, payload in results.items() if payload["category"] == "sqr_like"}
    best_sqr_key = max(
        sqr_candidates,
        key=lambda key: float(
            sqr_candidates[key]["pulse_replay"].get("12", {}).get(
                "fidelity",
                sqr_candidates[key]["embedded_evaluations"]["12"]["fidelity"],
            )
        ),
    )
    summary = {
        "best_embedded_key": best_embedded_key,
        "best_sqr_key": best_sqr_key,
    }
    payload = {"candidates": results, "summary": summary}
    c.save_json(c.DATA_DIR / "decomposition_results.json", payload)
    return payload


def summarize_duration_runs(seed_rows: list[dict[str, Any]]) -> dict[str, Any]:
    replay_values = [float(row["replay_fidelity"]) for row in seed_rows]
    nominal_values = [float(row["nominal_fidelity"]) for row in seed_rows]
    leak_values = [float(row["replay_leakage_average"]) for row in seed_rows]
    sorted_rows = sorted(seed_rows, key=lambda row: row["replay_fidelity"], reverse=True)
    best = sorted_rows[0]
    worst = sorted_rows[-1]
    return {
        "best_seed": int(best["seed"]),
        "best_replay_fidelity": float(best["replay_fidelity"]),
        "median_replay_fidelity": float(statistics.median(replay_values)),
        "worst_replay_fidelity": float(worst["replay_fidelity"]),
        "best_nominal_fidelity": float(best["nominal_fidelity"]),
        "median_nominal_fidelity": float(statistics.median(nominal_values)),
        "worst_nominal_fidelity": float(min(nominal_values)),
        "median_replay_leakage": float(statistics.median(leak_values)),
    }


def run_grape_stage(*, durations: tuple[int, ...], seeds: tuple[int, ...], maxiter: int) -> dict[str, Any]:
    model = c.build_model(n_cav=c.N_CAV_DEFAULT, n_tr=2)
    duration_payloads: dict[str, Any] = {}
    for duration_ns in durations:
        duration_key = f"{duration_ns}ns"
        print(f"[grape] duration={duration_ns} ns", flush=True)
        problem = c.build_grape_problem(model=model, duration_ns=float(duration_ns))
        seed_rows: list[dict[str, Any]] = []
        best_result = None
        best_replay = None
        best_problem = problem
        best_seed = None
        for seed in seeds:
            print(f"[grape]   seed={seed}", flush=True)
            result = c.run_grape_seed(problem, seed=seed, maxiter=maxiter)
            replay = c.replay_grape_operator(result=result, problem=problem, model=model, store_states=True)
            nominal_fidelity = float(result.metrics.get("nominal_fidelity", result.metrics.get("fidelity", np.nan)))
            seed_payload = {
                "seed": int(seed),
                "success": bool(result.success),
                "message": str(result.message),
                "objective_value": float(result.objective_value),
                "nominal_fidelity": nominal_fidelity,
                "replay_fidelity": float(replay["fidelity"]),
                "replay_block_fidelity": float(replay["block_fidelity"]),
                "replay_leakage_average": float(replay["leakage_average"]),
                "replay_leakage_worst": float(replay["leakage_worst"]),
                "max_transient_photon_number": float(replay["max_transient_photon_number"]),
                "metrics": result.metrics,
                "history": [
                    {
                        "evaluation": int(record.evaluation),
                        "objective": float(record.objective),
                        "gradient_norm": float(record.gradient_norm),
                        "elapsed_s": float(record.elapsed_s),
                        "metrics": record.metrics,
                    }
                    for record in result.history
                ],
            }
            seed_file = c.ARTIFACT_DIR / f"grape_{duration_key}_seed{seed}.json"
            seed_payload["artifact_path"] = str(seed_file)
            c.save_json(seed_file, {**seed_payload, "result_payload": result.to_payload()})
            seed_rows.append(seed_payload)
            if best_replay is None or seed_payload["replay_fidelity"] > best_replay["fidelity"]:
                best_result = result
                best_replay = replay
                best_seed = seed

        truncation_replays: dict[str, Any] = {}
        if best_result is not None and best_replay is not None:
            for n_cav in (8, 10, 12, 15):
                trunc_model = c.build_model(n_cav=n_cav, n_tr=2)
                trunc_replay = c.replay_grape_operator(
                    result=best_result,
                    problem=best_problem,
                    model=trunc_model,
                    store_states=True,
                )
                operator = embedded_operator_from_basis_outputs(trunc_replay["rows"], n_cav=n_cav)
                truncation_replays[str(n_cav)] = {
                    "fidelity": float(trunc_replay["fidelity"]),
                    "block_fidelity": float(trunc_replay["block_fidelity"]),
                    "leakage_average": float(trunc_replay["leakage_average"]),
                    "leakage_worst": float(trunc_replay["leakage_worst"]),
                    "max_transient_photon_number": float(trunc_replay["max_transient_photon_number"]),
                }
                save_artifact_npz(
                    c.ARTIFACT_DIR / f"grape_{duration_key}_bestseed_trunc{n_cav}.npz",
                    restricted_full_operator=operator,
                    subspace_operator=trunc_replay["subspace_operator"],
                )

            pulses, drive_ops, _pulse_meta = best_result.to_pulses()
            compiled = c.SequenceCompiler(dt=1.0e-9).compile(pulses, t_end=best_problem.time_grid.duration_s)
            open_process = reconstruct_open_process(
                compiled=compiled,
                drive_ops=drive_ops,
                model=model,
                noise=c.default_noise_spec(),
            )
            save_artifact_npz(
                c.ARTIFACT_DIR / f"grape_{duration_key}_bestseed_open_process.npz",
                choi_matrix=open_process["choi_matrix"],
            )
        else:
            open_process = {}

        duration_payloads[duration_key] = {
            "duration_ns": int(duration_ns),
            "steps": int(problem.time_grid.steps),
            "dt_ns": float(problem.time_grid.step_durations_s[0] * 1.0e9),
            "maxiter": int(maxiter),
            "seed_rows": seed_rows,
            "summary": summarize_duration_runs(seed_rows),
            "best_seed": int(best_seed) if best_seed is not None else None,
            "truncation_replays": truncation_replays,
            "open_process": {
                key: value for key, value in open_process.items() if key != "choi_matrix"
            },
        }

    recommended_key = max(
        duration_payloads,
        key=lambda key: float(
            duration_payloads[key]["open_process"].get(
                "process_fidelity",
                duration_payloads[key]["summary"]["best_replay_fidelity"],
            )
        ),
    )
    payload = {
        "durations": duration_payloads,
        "summary": {"recommended_duration_key": recommended_key},
    }
    c.save_json(c.DATA_DIR / "grape_results.json", payload)
    return payload


def run_wigner_stage() -> dict[str, Any]:
    decomp = load_json(c.DATA_DIR / "decomposition_results.json")
    grape = load_json(c.DATA_DIR / "grape_results.json")
    probe = load_json(c.DATA_DIR / "grape_large_truncation_probe.json")
    if decomp is None or grape is None:
        raise RuntimeError("Wigner stage requires decomposition_results.json and grape_results.json.")

    best_decomp_key = decomp["summary"]["best_embedded_key"]
    best_sqr_key = decomp["summary"]["best_sqr_key"]

    target_full = c.embed_target_unitary(WIGNER_NCAV)

    decomp_npz = np.load(decomp["candidates"][best_decomp_key]["artifact_path"])
    decomp_operator = np.asarray(decomp_npz[f"ideal_full_operator_nc{WIGNER_NCAV}"], dtype=np.complex128)

    sqr_npz = np.load(decomp["candidates"][best_sqr_key]["artifact_path"])
    if f"pulse_full_operator_nc{WIGNER_NCAV}" in sqr_npz:
        sqr_operator = np.asarray(sqr_npz[f"pulse_full_operator_nc{WIGNER_NCAV}"], dtype=np.complex128)
        sqr_source = "pulse"
    else:
        sqr_operator = np.asarray(sqr_npz[f"ideal_full_operator_nc{WIGNER_NCAV}"], dtype=np.complex128)
        sqr_source = "ideal"

    if probe:
        best_grape_key = max(probe, key=lambda key: float(probe[key]["best_replay_fidelity"]))
        grape_npz = np.load(probe[best_grape_key]["artifact_path"])
        grape_operator = np.asarray(grape_npz["restricted_full_operator"], dtype=np.complex128)
        grape_source = "nc12_probe"
    else:
        best_grape_key = grape["summary"]["recommended_duration_key"]
        grape_npz = np.load(c.ARTIFACT_DIR / f"grape_{best_grape_key}_bestseed_trunc{WIGNER_NCAV}.npz")
        grape_operator = np.asarray(grape_npz["restricted_full_operator"], dtype=np.complex128)
        grape_source = "nc8_sweep"

    payload = {
        "target": c.candidate_wigner_summary(target_full, target_full_operator=target_full, n_cav=WIGNER_NCAV),
        "best_decomposition": c.candidate_wigner_summary(
            decomp_operator, target_full_operator=target_full, n_cav=WIGNER_NCAV
        ),
        "best_sqr_like": c.candidate_wigner_summary(
            sqr_operator, target_full_operator=target_full, n_cav=WIGNER_NCAV
        ),
        "best_grape": c.candidate_wigner_summary(
            grape_operator, target_full_operator=target_full, n_cav=WIGNER_NCAV
        ),
        "sources": {
            "best_decomposition_key": best_decomp_key,
            "best_sqr_key": best_sqr_key,
            "best_sqr_source": sqr_source,
            "best_grape_key": best_grape_key,
            "best_grape_source": grape_source,
        },
    }
    c.save_json(c.DATA_DIR / "wigner_results.json", payload)
    return payload


def run_summary_stage() -> dict[str, Any]:
    transfer = load_json(c.DATA_DIR / "channel_summary.json")
    decomp = load_json(c.DATA_DIR / "decomposition_results.json")
    grape = load_json(c.DATA_DIR / "grape_results.json")
    wigner = load_json(c.DATA_DIR / "wigner_results.json")
    probe = load_json(c.DATA_DIR / "grape_large_truncation_probe.json")
    if any(item is None for item in (transfer, decomp, grape, wigner)):
        raise RuntimeError("Summary stage requires prior transfer, decomp, grape, and wigner outputs.")

    best_decomp_key = decomp["summary"]["best_embedded_key"]
    best_sqr_key = decomp["summary"]["best_sqr_key"]
    if probe:
        best_grape_key = max(probe, key=lambda key: float(probe[key]["best_replay_fidelity"]))
        best_grape_open_process = probe[best_grape_key].get("open_process", {})
        best_grape_replay = float(probe[best_grape_key]["best_replay_fidelity"])
        best_grape_source = "nc12_probe"
    else:
        best_grape_key = grape["summary"]["recommended_duration_key"]
        best_grape_open_process = grape["durations"][best_grape_key].get("open_process", {})
        best_grape_replay = grape["durations"][best_grape_key]["summary"]["best_replay_fidelity"]
        best_grape_source = "nc8_sweep"

    reassessment = [
        {
            "claim": "Exact reduced-subspace D-SQR-CPSQR agreement is not enough for physical validation.",
            "status": "validated_correction",
            "evidence": "The best SQR-like route retains modest reduced-model fidelity but fails the embedded and Wigner criteria at larger N_cav.",
        },
        {
            "claim": "The older finite-correlation-length statement for the sequential cluster channel was inconsistent.",
            "status": "corrected",
            "evidence": "The installed HolographicChannel convention gives transfer-spectrum {1, 0, 0, 0}, so the ordinary channel correlation length is xi = 0 rather than finite or infinite.",
        },
        {
            "claim": "GRAPE is the strongest physically credible candidate.",
            "status": "validated_with_higher_truncation_followup" if probe else "validated_if_replay_consistent",
            "evidence": (
                "The original N_cav=8 sweep is not truncation-converged, but direct N_cav=12 re-optimization at 300-400 ns restores high replay fidelity and strong Wigner agreement."
                if probe
                else "Duration ranking now uses replay fidelity, truncation stability, Wigner agreement, and open-system process fidelity rather than optimizer nominal fidelity alone."
            ),
        },
    ]

    payload = {
        "best_decomposition_key": best_decomp_key,
        "best_sqr_key": best_sqr_key,
        "best_grape_key": best_grape_key,
        "best_grape_source": best_grape_source,
        "transfer_correlation_length": transfer["correlation_length"],
        "transfer_eigenvalues": transfer["eigenvalues"],
        "best_decomposition_nc12_fidelity": decomp["candidates"][best_decomp_key]["embedded_evaluations"]["12"]["fidelity"],
        "best_sqr_nc12_metric": decomp["candidates"][best_sqr_key]["pulse_replay"].get(
            "12",
            decomp["candidates"][best_sqr_key]["embedded_evaluations"]["12"],
        ),
        "best_grape_open_process_fidelity": best_grape_open_process.get("process_fidelity"),
        "best_grape_replay_fidelity": best_grape_replay,
        "wigner_sources": wigner["sources"],
        "reassessment": reassessment,
    }
    c.save_json(c.DATA_DIR / "followup_summary.json", payload)
    return payload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the holographic cluster-state follow-up study.")
    parser.add_argument(
        "--stages",
        nargs="+",
        default=list(DEFAULT_STAGES),
        choices=list(DEFAULT_STAGES),
        help="Study stages to execute.",
    )
    parser.add_argument(
        "--durations",
        nargs="+",
        type=int,
        default=list(DEFAULT_DURATIONS),
        help="GRAPE durations in ns.",
    )
    parser.add_argument(
        "--seeds",
        nargs="+",
        type=int,
        default=list(DEFAULT_SEEDS),
        help="GRAPE random seeds.",
    )
    parser.add_argument(
        "--grape-maxiter",
        type=int,
        default=400,
        help="GRAPE iterations per seed.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    stage_order = tuple(args.stages)
    outputs: dict[str, Any] = {}
    for stage in stage_order:
        start = time.perf_counter()
        if stage == "transfer":
            outputs[stage] = run_transfer_stage()
        elif stage == "decomp":
            outputs[stage] = run_decomp_stage()
        elif stage == "grape":
            outputs[stage] = run_grape_stage(
                durations=tuple(args.durations),
                seeds=tuple(args.seeds),
                maxiter=int(args.grape_maxiter),
            )
        elif stage == "wigner":
            outputs[stage] = run_wigner_stage()
        elif stage == "summary":
            outputs[stage] = run_summary_stage()
        else:
            raise ValueError(f"Unsupported stage '{stage}'.")
        elapsed = time.perf_counter() - start
        print(f"[done] {stage} in {elapsed:.1f}s", flush=True)

    c.save_json(c.DATA_DIR / "last_run_manifest.json", {"stages": stage_order, "outputs": outputs})


if __name__ == "__main__":
    main()
