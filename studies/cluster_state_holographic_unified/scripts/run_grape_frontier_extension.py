from __future__ import annotations

import argparse
import time
from pathlib import Path
from typing import Any

import numpy as np
import qutip as qt

import common as c


DEFAULT_DURATIONS_NS = (500, 600)
DEFAULT_SEEDS = (17, 42, 73)
DEFAULT_MAXITER = 250
DEFAULT_N_CAV = 12
DEFAULT_ENGINE = "auto"


def time_grid_dt_ns(problem: Any) -> float:
    return float(problem.time_grid.duration_s) * 1.0e9 / float(problem.time_grid.steps)


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
    dim = target_unitary.shape[0]
    target_vec = target_unitary.reshape(-1, order="F")
    choi_target = np.outer(target_vec, target_vec.conj())
    value = np.trace(choi_target.conj().T @ choi)
    return float(np.real(value) / (dim * dim))


def reconstruct_open_process(*, compiled: Any, drive_ops: dict[str, Any], model: c.DispersiveTransmonCavityModel, noise: c.NoiseSpec) -> dict[str, Any]:
    frame = c.build_frame(model)
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
        "process_fidelity": process_fidelity_from_choi(choi, c.TARGET_UNITARY),
        "mean_probe_leakage": float(np.mean(leakage_probe)),
        "max_probe_leakage": float(np.max(leakage_probe)),
        "mean_basis_state_fidelity": float(np.mean(target_state_fids)),
        "min_basis_state_fidelity": float(np.min(target_state_fids)),
        "choi_matrix": choi,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Extend the validated N_cav=12 GRAPE frontier.")
    parser.add_argument("--durations-ns", type=int, nargs="+", default=list(DEFAULT_DURATIONS_NS))
    parser.add_argument("--seeds", type=int, nargs="+", default=list(DEFAULT_SEEDS))
    parser.add_argument("--maxiter", type=int, default=DEFAULT_MAXITER)
    parser.add_argument("--n-cav", type=int, default=DEFAULT_N_CAV)
    parser.add_argument("--engine", choices=("auto", "numpy", "jax"), default=DEFAULT_ENGINE)
    parser.add_argument("--jax-device", type=str, default=None)
    args = parser.parse_args()

    n_cav = int(args.n_cav)
    requested_engine = str(args.engine).lower()
    resolved_engine = c.resolve_grape_engine(requested_engine)
    jax_device = str(args.jax_device) if args.jax_device else None
    model = c.build_model(n_cav=n_cav, n_tr=2)
    basis_states = list(zip(c.LOGICAL_LABELS, c.logical_basis_states(model), strict=True))
    output_path = c.DATA_DIR / "grape_frontier_extension.json"
    results: dict[str, Any] = c.load_json(output_path) if output_path.exists() else {}

    print(
        f"[grape-frontier] engine={resolved_engine} requested={requested_engine} jax_available={c.jax_available()} jax_device={jax_device}",
        flush=True,
    )

    for duration_ns in [int(value) for value in args.durations_ns]:
        print(f"[grape-frontier] N_cav={n_cav}, duration={duration_ns} ns", flush=True)
        problem = c.build_grape_problem(model=model, duration_ns=float(duration_ns))
        print(
            f"[grape-frontier]   built problem steps={int(problem.time_grid.steps)} dt_ns={time_grid_dt_ns(problem):.3f}",
            flush=True,
        )
        seed_rows: list[dict[str, Any]] = []
        best_row: dict[str, Any] | None = None
        best_result = None
        best_replay = None
        for seed in [int(value) for value in args.seeds]:
            print(f"[grape-frontier]   seed={seed} start", flush=True)
            start = time.perf_counter()
            result = c.run_grape_seed(
                problem,
                seed=seed,
                maxiter=int(args.maxiter),
                engine=resolved_engine,
                jax_device=jax_device,
            )
            replay = c.replay_grape_subspace(
                result=result,
                problem=problem,
                model=model,
                subspace=c.logical_subspace(n_cav),
                target_unitary=c.TARGET_UNITARY,
                basis_states=basis_states,
                noise=None,
                store_states=False,
            )
            elapsed_s = time.perf_counter() - start
            row = {
                "seed": int(seed),
                "elapsed_s": float(elapsed_s),
                "nominal_fidelity": float(c.summarise_grape_result(result, n_cav=n_cav)["nominal_fidelity"]),
                "replay_fidelity": float(replay["fidelity"]),
                "replay_leakage_average": float(replay["leakage_average"]),
                "replay_leakage_worst": float(replay["leakage_worst"]),
                "success": bool(result.success),
                "message": str(result.message),
            }
            seed_rows.append(row)
            results[f"{duration_ns}ns"] = {
                "duration_ns": int(duration_ns),
                "status": "seed-scan",
                "n_cav": int(n_cav),
                "steps": int(problem.time_grid.steps),
                "dt_ns": time_grid_dt_ns(problem),
                "maxiter": int(args.maxiter),
                "engine_requested": requested_engine,
                "engine_resolved": resolved_engine,
                "jax_available": bool(c.jax_available()),
                "jax_device": jax_device if resolved_engine == "jax" else None,
                "seed_rows": seed_rows,
            }
            c.save_json(output_path, results)
            print(
                f"[grape-frontier]   seed={seed} nominal={row['nominal_fidelity']:.6f} replay={row['replay_fidelity']:.6f} "
                f"leak={row['replay_leakage_worst']:.6f} elapsed={elapsed_s:.1f}s",
                flush=True,
            )
            if best_row is None or row["replay_fidelity"] > best_row["replay_fidelity"]:
                best_row = row
                best_result = result
                best_replay = replay

        if best_row is None or best_result is None or best_replay is None:
            continue

        print(f"[grape-frontier]   reconstructing open process for best seed {int(best_row['seed'])}", flush=True)
        open_process = reconstruct_open_process(
            compiled=best_replay["compiled"],
            drive_ops=best_replay["drive_ops"],
            model=model,
            noise=c.default_noise_spec(),
        )
        artifact_path = c.ARTIFACT_DIR / f"grape_frontier_nc{n_cav}_{duration_ns}ns_best.npz"
        np.savez(
            artifact_path,
            subspace_operator=np.asarray(best_replay["subspace_operator"], dtype=np.complex128),
            choi_matrix=np.asarray(open_process["choi_matrix"], dtype=np.complex128),
        )
        results[f"{duration_ns}ns"] = {
            "duration_ns": int(duration_ns),
            "status": "complete",
            "n_cav": int(n_cav),
            "steps": int(problem.time_grid.steps),
            "dt_ns": time_grid_dt_ns(problem),
            "maxiter": int(args.maxiter),
            "engine_requested": requested_engine,
            "engine_resolved": resolved_engine,
            "jax_available": bool(c.jax_available()),
            "jax_device": jax_device if resolved_engine == "jax" else None,
            "seed_rows": seed_rows,
            "best_seed": int(best_row["seed"]),
            "best_nominal_fidelity": float(best_row["nominal_fidelity"]),
            "best_replay_fidelity": float(best_row["replay_fidelity"]),
            "best_replay_leakage_worst": float(best_row["replay_leakage_worst"]),
            "open_process": {
                key: value
                for key, value in open_process.items()
                if key != "choi_matrix"
            },
            "artifact_path": str(artifact_path),
        }
        c.save_json(output_path, results)
        print(f"[grape-frontier] checkpointed {duration_ns} ns", flush=True)

    c.save_json(output_path, results)
    print("[grape-frontier] wrote data/grape_frontier_extension.json", flush=True)


if __name__ == "__main__":
    main()