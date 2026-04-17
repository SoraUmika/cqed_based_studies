from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Any

import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import common as c
from run_followup_study import reconstruct_open_process


DURATIONS_NS = (300, 400)
SEEDS = (17, 42, 73)
MAXITER = 200
N_CAV = 12


def embedded_operator_from_rows(rows: list[dict[str, Any]], *, n_cav: int) -> np.ndarray:
    full_dim = 2 * int(n_cav)
    operator = np.eye(full_dim, dtype=np.complex128)
    for index, row in zip(c.logical_indices(int(n_cav)), rows, strict=True):
        state = row["simulation"].final_state
        vector = np.asarray(state.full(), dtype=np.complex128).reshape(-1)
        operator[:, index] = vector
    return operator


def main() -> None:
    model = c.build_model(n_cav=N_CAV, n_tr=2)
    results: dict[str, Any] = {}
    for duration_ns in DURATIONS_NS:
        print(f"[probe] N_cav={N_CAV}, duration={duration_ns} ns", flush=True)
        problem = c.build_grape_problem(model=model, duration_ns=float(duration_ns))
        rows: list[dict[str, Any]] = []
        best_replay = -1.0
        best_result = None
        best_payload = None
        for seed in SEEDS:
            start = time.perf_counter()
            result = c.run_grape_seed(problem, seed=seed, maxiter=MAXITER)
            replay = c.replay_grape_operator(result=result, problem=problem, model=model, store_states=True)
            elapsed_s = time.perf_counter() - start
            seed_artifact = c.ARTIFACT_DIR / f"grape_nc{N_CAV}_{duration_ns}ns_seed{seed}_probe.json"
            payload = {
                "seed": int(seed),
                "elapsed_s": elapsed_s,
                "nominal_fidelity": float(result.metrics.get("nominal_fidelity", result.metrics.get("fidelity", np.nan))),
                "replay_fidelity": float(replay["fidelity"]),
                "replay_block_fidelity": float(replay["block_fidelity"]),
                "replay_leakage_average": float(replay["leakage_average"]),
                "replay_leakage_worst": float(replay["leakage_worst"]),
                "max_transient_photon_number": float(replay["max_transient_photon_number"]),
                "success": bool(result.success),
                "message": str(result.message),
                "artifact_path": str(seed_artifact),
            }
            rows.append(payload)
            c.save_json(seed_artifact, {"summary": payload, "result_payload": result.to_payload()})
            if payload["replay_fidelity"] > best_replay:
                best_replay = payload["replay_fidelity"]
                best_result = result
                best_payload = replay
        if best_result is None or best_payload is None:
            continue

        operator = embedded_operator_from_rows(best_payload["rows"], n_cav=N_CAV)
        artifact_path = c.ARTIFACT_DIR / f"grape_nc{N_CAV}_{duration_ns}ns_probe.npz"
        np.savez(
            artifact_path,
            restricted_full_operator=operator,
            subspace_operator=best_payload["subspace_operator"],
        )
        pulses, drive_ops, _pulse_meta = best_result.to_pulses()
        compiled = c.SequenceCompiler(dt=1.0e-9).compile(pulses, t_end=problem.time_grid.duration_s)
        open_process = reconstruct_open_process(
            compiled=compiled,
            drive_ops=drive_ops,
            model=model,
            noise=c.default_noise_spec(),
        )
        open_artifact = c.ARTIFACT_DIR / f"grape_nc{N_CAV}_{duration_ns}ns_open_process_probe.npz"
        np.savez(open_artifact, choi_matrix=open_process["choi_matrix"])
        wigner = c.candidate_wigner_summary(
            operator,
            target_full_operator=c.embed_target_unitary(N_CAV),
            n_cav=N_CAV,
        )
        results[f"{duration_ns}ns"] = {
            "duration_ns": int(duration_ns),
            "seed_rows": rows,
            "best_replay_fidelity": float(best_replay),
            "best_seed": int(max(rows, key=lambda row: row["replay_fidelity"])["seed"]),
            "artifact_path": str(artifact_path),
            "open_process": {
                key: value for key, value in open_process.items() if key != "choi_matrix"
            },
            "open_artifact_path": str(open_artifact),
            "wigner_summary": {
                key: {
                    "cavity_fidelity": value["cavity_fidelity"],
                    "wigner_l2": value["wigner_l2"],
                }
                for key, value in wigner.items()
            },
        }

    c.save_json(c.DATA_DIR / "grape_large_truncation_probe.json", results)


if __name__ == "__main__":
    main()
