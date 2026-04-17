from __future__ import annotations

import argparse
import json
import time
from datetime import date
from pathlib import Path
from typing import Any

import matplotlib
import numpy as np
import qutip as qt
from scipy.linalg import logm

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cqed_sim.sim import NoiseSpec
from cqed_sim.unitary_synthesis import CQEDSystemAdapter, subspace_unitary_fidelity
from cqed_sim.unitary_synthesis.backends import simulate_sequence

from common import ARTIFACT_DIR, DATA_DIR, FIG_DIR, FAMILY_SPECS, REFINE_SEEDS, REFINE_MAXITER, base, build_family_sequence, evaluate_sequence_with_diagnostics, sequence_complexity_metrics, sequence_for_n_cav, sequence_from_payload, sequence_phase_budget


STYLE_PATH = Path(__file__).resolve().parents[3] / ".github" / "skills" / "publication-figures" / "assets" / "cqed_style.mplstyle"
if STYLE_PATH.exists():
    plt.style.use(str(STYLE_PATH))


STUDY_NAME = "cluster_state_free_evolution_snap_comparison"
TARGET_CASES = (
    ("native_fe", 6),
    ("snap_interleaved", 6),
)
REOPT_N_CAV = 6
VALIDATION_LEVELS = (6, 8)
PULSE_DT_S = float(base.GRAPE_DT_S)


def _json_dump(path: Path, payload: Any) -> None:
    base.save_json(path, payload)


def _logical_probe_states(n_cav: int) -> list[dict[str, Any]]:
    single_qubit = {
        "g": np.array([1.0, 0.0], dtype=np.complex128),
        "e": np.array([0.0, 1.0], dtype=np.complex128),
        "+x": np.array([1.0, 1.0], dtype=np.complex128) / np.sqrt(2.0),
        "+y": np.array([1.0, 1.0j], dtype=np.complex128) / np.sqrt(2.0),
    }
    cavity = {
        "0": np.array([1.0, 0.0], dtype=np.complex128),
        "1": np.array([0.0, 1.0], dtype=np.complex128),
        "+": np.array([1.0, 1.0], dtype=np.complex128) / np.sqrt(2.0),
        "+i": np.array([1.0, 1.0j], dtype=np.complex128) / np.sqrt(2.0),
    }
    logical_subspace = base.logical_subspace(int(n_cav))
    probes: list[dict[str, Any]] = []
    for qubit_label, qubit_vec in single_qubit.items():
        for cavity_label, cavity_vec in cavity.items():
            logical_vec = np.kron(qubit_vec, cavity_vec)
            full_vec = logical_subspace.embed(logical_vec)
            probes.append(
                {
                    "label": f"{qubit_label}|{cavity_label}",
                    "logical_state": logical_vec,
                    "input_state": qt.Qobj(full_vec.reshape(-1, 1), dims=[[2, int(n_cav)], [1]]),
                }
            )
    return probes


def _logical_subspace_population(state: qt.Qobj, *, n_cav: int) -> float:
    indices = tuple(int(idx) for idx in base.logical_indices(int(n_cav)))
    dense = np.asarray(state.full(), dtype=np.complex128)
    if state.isket:
        vector = dense.reshape(-1)
        return float(np.sum(np.abs(vector[list(indices)]) ** 2))
    return float(np.real(np.trace(dense[np.ix_(indices, indices)])))


def _target_state_from_logical(logical_vec: np.ndarray, *, n_cav: int) -> qt.Qobj:
    logical_subspace = base.logical_subspace(int(n_cav))
    target_vec = base.TARGET_UNITARY @ np.asarray(logical_vec, dtype=np.complex128)
    full_vec = logical_subspace.embed(target_vec)
    return qt.Qobj(full_vec.reshape(-1, 1), dims=[[2, int(n_cav)], [1]])


def _pulse_replay_metrics(sequence, *, n_cav: int, noise: NoiseSpec | None = None) -> dict[str, Any]:
    logical_subspace = base.logical_subspace(int(n_cav))
    model = base.build_model(n_cav=int(n_cav), n_tr=base.N_TR)
    system = CQEDSystemAdapter(model=model)
    probes = _logical_probe_states(int(n_cav))
    result = simulate_sequence(
        sequence,
        logical_subspace,
        backend="pulse",
        system=system,
        dt=PULSE_DT_S,
        frame=base.build_frame(model),
        need_operator=noise is None,
        state_inputs=[row["input_state"] for row in probes],
        noise=noise,
    )
    fidelities: list[float] = []
    leakages: list[float] = []
    for probe, output in zip(probes, result.state_outputs or [], strict=True):
        target_state = _target_state_from_logical(probe["logical_state"], n_cav=int(n_cav))
        fidelities.append(float(qt.metrics.fidelity(target_state, output)))
        leakages.append(1.0 - _logical_subspace_population(output, n_cav=int(n_cav)))
    payload = {
        "n_cav": int(n_cav),
        "probe_count": len(probes),
        "probe_fidelity_mean": float(np.mean(fidelities)),
        "probe_fidelity_min": float(np.min(fidelities)),
        "logical_leakage_mean": float(np.mean(leakages)),
        "logical_leakage_max": float(np.max(leakages)),
        "unitarity_error": float(result.metrics.get("unitarity_error", np.nan)),
    }
    if noise is None and result.subspace_operator is not None:
        payload["operator_fidelity"] = float(
            subspace_unitary_fidelity(
                np.asarray(result.subspace_operator, dtype=np.complex128),
                base.TARGET_UNITARY,
                gauge="global",
            )
        )
    return payload


def _noise_collapse_ops(*, n_cav: int, noise: NoiseSpec) -> list[qt.Qobj]:
    ops: list[qt.Qobj] = []
    qeye_c = qt.qeye(int(n_cav))
    if noise.t1 is not None:
        ops.append(np.sqrt(1.0 / float(noise.t1)) * qt.tensor(qt.destroy(int(base.N_TR)), qeye_c))
    if noise.tphi is not None:
        ops.append(np.sqrt(0.5 / float(noise.tphi)) * qt.tensor(qt.sigmaz(), qeye_c))
    if noise.kappa is not None:
        ops.append(np.sqrt(float(noise.kappa)) * qt.tensor(qt.qeye(int(base.N_TR)), qt.destroy(int(n_cav))))
    return ops


def _lindblad_surrogate_metrics(sequence, *, n_cav: int, noise: NoiseSpec) -> dict[str, Any]:
    probes = _logical_probe_states(int(n_cav))
    collapse_ops = _noise_collapse_ops(n_cav=int(n_cav), noise=noise)
    fidelities: list[float] = []
    leakages: list[float] = []
    for probe in probes:
        state: qt.Qobj = qt.ket2dm(probe["input_state"])
        for gate in sequence.gates:
            duration = max(float(gate.duration), 1.0e-12)
            unitary = np.asarray(gate.pulse_unitary(int(n_cav)).full(), dtype=np.complex128)
            generator = 1.0j * logm(unitary) / duration
            generator = 0.5 * (generator + generator.conj().T)
            hamiltonian = qt.Qobj(generator, dims=[[int(base.N_TR), int(n_cav)], [int(base.N_TR), int(n_cav)]])
            state = qt.mesolve(hamiltonian, state, [0.0, duration], c_ops=collapse_ops).states[-1]
        target_state = _target_state_from_logical(probe["logical_state"], n_cav=int(n_cav))
        fidelities.append(float(qt.metrics.fidelity(target_state, state)))
        leakages.append(1.0 - _logical_subspace_population(state, n_cav=int(n_cav)))
    return {
        "method": "local_lindblad_surrogate",
        "n_cav": int(n_cav),
        "probe_count": len(probes),
        "probe_fidelity_mean": float(np.mean(fidelities)),
        "probe_fidelity_min": float(np.min(fidelities)),
        "logical_leakage_mean": float(np.mean(leakages)),
        "logical_leakage_max": float(np.max(leakages)),
    }


def _sequence_summary_dict(*, family: str, blocks: int, n_cav: int, seed: int, maxiter: int, fit_payload: dict[str, Any]) -> dict[str, Any]:
    sequence = fit_payload["result"].sequence
    evaluation = evaluate_sequence_with_diagnostics(sequence, n_cav=int(n_cav))
    complexity = sequence_complexity_metrics(sequence)
    phase_budget = sequence_phase_budget(sequence)
    seq_summary = fit_payload["summary"]
    return {
        "study_name": STUDY_NAME,
        "family": family,
        "family_title": str(FAMILY_SPECS[str(family)]["title"]),
        "case_id": f"{family}_b{int(blocks)}_n{int(n_cav)}",
        "blocks": int(blocks),
        "optimized_n_cav": int(n_cav),
        "seed": int(seed),
        "maxiter": int(maxiter),
        "fidelity": float(evaluation["fidelity"]),
        "block_gauge_fidelity": float(evaluation["block_gauge_fidelity"]),
        "best_fit_block_gauge_fidelity": float(evaluation["best_fit_block_gauge_fidelity"]),
        "leakage_average": float(evaluation["leakage_average"]),
        "leakage_worst": float(evaluation["leakage_worst"]),
        "unitarity_error": float(evaluation["unitarity_error"]),
        "objective": float(fit_payload["objective"]),
        "success": bool(fit_payload["success"]),
        "message": str(fit_payload["message"]),
        "metrics": fit_payload["metrics"],
        "gate_depth": int(seq_summary["gate_depth"]),
        "total_duration_ns": float(seq_summary["total_duration_ns"]),
        "total_wait_time_ns": float(phase_budget["total_wait_time_ns"]),
        "total_fe_logical_delta_phi_rad": float(phase_budget["total_fe_logical_delta_phi_rad"]),
        "total_snap_logical_relative_phase_rad": float(phase_budget["total_snap_logical_relative_phase_rad"]),
        "parameter_count": int(complexity["parameter_count"]),
        "snap_gate_count": int(complexity["snap_gate_count"]),
        "snap_phase_count": int(complexity["snap_phase_count"]),
        "entangling_gate_count": int(complexity["entangling_gate_count"]),
        "wait_gate_count": int(complexity["wait_gate_count"]),
        "sequence_summary": seq_summary,
        "sequence_payload": fit_payload["sequence_payload"],
    }


def _reoptimize_case(*, family: str, blocks: int, n_cav: int, seeds: tuple[int, ...], maxiter: int) -> dict[str, Any]:
    best: dict[str, Any] | None = None
    for seed in seeds:
        sequence = build_family_sequence(family=family, blocks=blocks, n_cav=int(n_cav))
        t0 = time.perf_counter()
        fit_payload = base.fit_sequence(
            sequence,
            n_cav=int(n_cav),
            seed=int(seed),
            init_guess="heuristic",
            multistart=1,
            maxiter=int(maxiter),
        )
        summary = _sequence_summary_dict(
            family=family,
            blocks=blocks,
            n_cav=n_cav,
            seed=seed,
            maxiter=maxiter,
            fit_payload=fit_payload,
        )
        summary["elapsed_s"] = float(time.perf_counter() - t0)
        if best is None or float(summary["fidelity"]) > float(best["fidelity"]):
            best = summary
    assert best is not None
    return best


def _evaluate_saved_baseline(*, artifact_name: str, n_cav: int) -> dict[str, Any]:
    payload = json.loads((ARTIFACT_DIR / artifact_name).read_text(encoding="utf-8"))
    sequence = sequence_from_payload(payload["sequence_payload"], n_cav=int(n_cav))
    ideal = evaluate_sequence_with_diagnostics(sequence, n_cav=int(n_cav))
    pulse = _pulse_replay_metrics(sequence, n_cav=int(n_cav), noise=None)
    noisy = _lindblad_surrogate_metrics(sequence, n_cav=int(n_cav), noise=_nominal_noise())
    return {
        "artifact_name": artifact_name,
        "n_cav": int(n_cav),
        "ideal": {
            "fidelity": float(ideal["fidelity"]),
            "block_gauge_fidelity": float(ideal["block_gauge_fidelity"]),
            "leakage_average": float(ideal["leakage_average"]),
            "leakage_worst": float(ideal["leakage_worst"]),
        },
        "pulse": pulse,
        "noisy": noisy,
    }


def _evaluate_sequence_family(summary: dict[str, Any], *, validation_levels: tuple[int, ...]) -> dict[str, Any]:
    base_sequence = sequence_from_payload(summary["sequence_payload"], n_cav=int(summary["optimized_n_cav"]))
    out: dict[str, Any] = {"optimized_n_cav": int(summary["optimized_n_cav"]), "by_n_cav": {}}
    for n_cav in validation_levels:
        sequence = sequence_for_n_cav(base_sequence, n_cav=int(n_cav))
        ideal = evaluate_sequence_with_diagnostics(sequence, n_cav=int(n_cav))
        out["by_n_cav"][str(int(n_cav))] = {
            "ideal": {
                "fidelity": float(ideal["fidelity"]),
                "block_gauge_fidelity": float(ideal["block_gauge_fidelity"]),
                "leakage_average": float(ideal["leakage_average"]),
                "leakage_worst": float(ideal["leakage_worst"]),
            },
            "pulse": _pulse_replay_metrics(sequence, n_cav=int(n_cav), noise=None),
            "noisy": _lindblad_surrogate_metrics(sequence, n_cav=int(n_cav), noise=_nominal_noise()),
        }
    return out


def _nominal_noise() -> NoiseSpec:
    return NoiseSpec(t1=30.0e-6, tphi=20.0e-6, kappa=1.0 / 200.0e-6)


def _artifact_payload(summary: dict[str, Any], *, description: str, parameters: dict[str, Any]) -> dict[str, Any]:
    payload = dict(summary)
    payload.update(
        {
            "study_name": STUDY_NAME,
            "date_created": date.today().isoformat(),
            "description": description,
            "parameters": parameters,
            "load_instructions": (
                "import json; from pathlib import Path; payload = json.loads(Path(filename).read_text(encoding='utf-8'))"
            ),
        }
    )
    return payload


def _plot_reoptimization_comparison(summary: dict[str, Any]) -> None:
    families = ["native_fe", "snap_interleaved"]
    labels = [str(FAMILY_SPECS[family]["title"]) for family in families]
    baseline_vals = [float(summary["baseline_replay"][family]["ideal"]["fidelity"]) for family in families]
    reopt_vals = [float(summary["reoptimized"][family]["evaluation"]["by_n_cav"][str(REOPT_N_CAV)]["ideal"]["fidelity"]) for family in families]
    noisy_vals = [float(summary["reoptimized"][family]["evaluation"]["by_n_cav"][str(REOPT_N_CAV)]["noisy"]["probe_fidelity_mean"]) for family in families]

    x = np.arange(len(families))
    width = 0.24
    fig, ax = plt.subplots(figsize=(8.0, 4.6))
    ax.bar(x - width, baseline_vals, width, label="Baseline replay at n_cav=6", color="#4477AA")
    ax.bar(x, reopt_vals, width, label="Re-optimized ideal at n_cav=6", color="#228833")
    ax.bar(x + width, noisy_vals, width, label="Re-optimized noisy replay", color="#CCBB44")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Fidelity")
    ax.set_title("Larger-Truncation Improvement Pass")
    ax.set_ylim(0.0, 1.01)
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"improvement_reoptimization_comparison.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def _plot_truncation_and_noise(summary: dict[str, Any]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4))
    colors = {"native_fe": "#4477AA", "snap_interleaved": "#228833"}
    titles = {"native_fe": "Native FE", "snap_interleaved": "Interleaved SNAP"}

    for family in ("native_fe", "snap_interleaved"):
        xs = [int(level) for level in summary["reoptimized"][family]["evaluation"]["by_n_cav"].keys()]
        ideal = [
            float(summary["reoptimized"][family]["evaluation"]["by_n_cav"][str(level)]["ideal"]["fidelity"])
            for level in xs
        ]
        pulse = [
            float(summary["reoptimized"][family]["evaluation"]["by_n_cav"][str(level)]["pulse"]["operator_fidelity"])
            for level in xs
        ]
        noisy = [
            float(summary["reoptimized"][family]["evaluation"]["by_n_cav"][str(level)]["noisy"]["probe_fidelity_mean"])
            for level in xs
        ]
        axes[0].plot(xs, ideal, "o-", color=colors[family], label=f"{titles[family]} ideal")
        axes[0].plot(xs, pulse, "s--", color=colors[family], alpha=0.8, label=f"{titles[family]} pulse")
        axes[1].plot(xs, noisy, "o-", color=colors[family], label=titles[family])

    axes[0].set_xlabel("Cavity truncation")
    axes[0].set_ylabel("Fidelity")
    axes[0].set_title("Ideal vs Pulse Replay")
    axes[0].set_ylim(0.0, 1.01)
    axes[0].grid(True, alpha=0.25)
    axes[0].legend(fontsize=8)

    axes[1].set_xlabel("Cavity truncation")
    axes[1].set_ylabel("Mean probe-state fidelity")
    axes[1].set_title("Nominal Noise Replay")
    axes[1].set_ylim(0.0, 1.01)
    axes[1].grid(True, alpha=0.25)
    axes[1].legend(fontsize=8)

    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"improvement_truncation_noise_summary.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the cluster-state free-evolution/SNAP improvement pass.")
    parser.add_argument("--maxiter", type=int, default=REFINE_MAXITER)
    args = parser.parse_args()

    baseline_replay = {
        "native_fe": _evaluate_saved_baseline(artifact_name="best_native_sequence.json", n_cav=REOPT_N_CAV),
        "snap_interleaved": _evaluate_saved_baseline(artifact_name="best_snap_sequence.json", n_cav=REOPT_N_CAV),
    }

    reoptimized: dict[str, Any] = {}
    for family, blocks in TARGET_CASES:
        print(f"[reopt] {family} blocks={blocks} n_cav={REOPT_N_CAV}", flush=True)
        artifact_name = f"{family}_b{blocks}_n{REOPT_N_CAV}_reoptimized.json"
        artifact_path = ARTIFACT_DIR / artifact_name
        if artifact_path.exists():
            best = json.loads(artifact_path.read_text(encoding="utf-8"))
        else:
            best = _reoptimize_case(
                family=family,
                blocks=blocks,
                n_cav=REOPT_N_CAV,
                seeds=tuple(int(seed) for seed in REFINE_SEEDS),
                maxiter=int(args.maxiter),
            )
        evaluation = _evaluate_sequence_family(best, validation_levels=VALIDATION_LEVELS)
        reoptimized[family] = {"summary": best, "evaluation": evaluation}
        _json_dump(
            ARTIFACT_DIR / artifact_name,
            _artifact_payload(
                best,
                description=(
                    f"Best {family} sequence re-optimized at n_cav={REOPT_N_CAV} for the improvement pass. "
                    "Noisy follow-up uses a local Lindblad surrogate because cqed_sim gate-sequence replay does not apply NoiseSpec for these gates."
                ),
                parameters={"family": family, "blocks": blocks, "optimized_n_cav": REOPT_N_CAV, "seeds": list(REFINE_SEEDS), "maxiter": int(args.maxiter)},
            ),
        )

    summary = {
        "study_name": STUDY_NAME,
        "date_created": date.today().isoformat(),
        "description": "Improvement pass covering larger-truncation re-optimization, surrogate pulse replay, and Lindblad-noise surrogate replay.",
        "parameters": {
            "reoptimization_n_cav": REOPT_N_CAV,
            "validation_levels": list(VALIDATION_LEVELS),
            "seeds": list(REFINE_SEEDS),
            "maxiter": int(args.maxiter),
            "noise": {"t1_s": 30.0e-6, "tphi_s": 20.0e-6, "kappa_hz": 1.0 / 200.0e-6},
            "noise_method": "local_lindblad_surrogate_from_gate_pulse_unitaries",
            "pulse_dt_s": PULSE_DT_S,
        },
        "baseline_replay": baseline_replay,
        "reoptimized": reoptimized,
    }
    _json_dump(DATA_DIR / "improvement_pass_summary.json", summary)

    _plot_reoptimization_comparison(summary)
    _plot_truncation_and_noise(summary)

    print("[done] wrote improvement_pass_summary.json and updated figures/artifacts", flush=True)


if __name__ == "__main__":
    main()