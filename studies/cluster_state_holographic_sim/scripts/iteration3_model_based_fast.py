"""Iteration 3 — Efficient model-based verification.

Streamlined version: lower multistart/maxiter for manageable compute at N_cav=8.
All output goes to both stdout and a log file for monitoring.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
STUDY_DIR = SCRIPT_DIR.parent
DATA_DIR = STUDY_DIR / "data"
FIG_DIR = STUDY_DIR / "figures"
ARTIFACTS_DIR = STUDY_DIR / "artifacts"
LOG_FILE = DATA_DIR / "iteration3_model_based.log"

CQED_SIM = Path(
    r"C:\Users\dazzl\Box\Shyam Shankar Quantum Circuits Group"
    r"\Users\Users_JianJun\cQED_simulation"
)
if str(CQED_SIM) not in sys.path:
    sys.path.insert(0, str(CQED_SIM))

from cqed_sim.unitary_synthesis import (
    Displacement, QubitRotation, SQR, ConditionalPhaseSQR,
    FreeEvolveCondPhase, Subspace, TargetUnitary, UnitarySynthesizer,
    GateSequence, DriftPhaseModel, LeakagePenalty, MultiObjective,
    ExecutionOptions, subspace_unitary_fidelity, simulate_sequence,
)
from cqed_sim.unitary_synthesis.targets import make_target

TWO_PI = 2 * np.pi
CHI = TWO_PI * (-2.84e6)
CHIP = TWO_PI * (-21e3)
KERR = TWO_PI * (-28e3)

no_drift = DriftPhaseModel(chi=0.0, chi2=0.0, kerr=0.0)
full_drift = DriftPhaseModel(chi=CHI, chi2=CHIP, kerr=KERR)

U_target = make_target("cluster", n_match=1)
TAU_CZ = np.pi / abs(CHI)

DATA_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)
ARTIFACTS_DIR.mkdir(parents=True, exist_ok=True)


def log(msg: str) -> None:
    """Print and append to log file."""
    print(msg, flush=True)
    with open(LOG_FILE, "a", encoding="utf-8") as f:
        f.write(msg + "\n")


def make_subspace(n_cav: int) -> Subspace:
    full_dim = 2 * n_cav
    return Subspace.custom(full_dim, [0, 1, n_cav, n_cav + 1],
                           ["|g,0>", "|g,1>", "|e,0>", "|e,1>"])


def run_synthesis(
    sequence: GateSequence,
    label: str,
    n_cav: int,
    multistart: int = 6,
    maxiter: int = 400,
) -> dict:
    subspace = make_subspace(n_cav)
    target = TargetUnitary(U_target, ignore_global_phase=True)
    log(f"\n  [{label}] n_cav={n_cav} ms={multistart} maxiter={maxiter}")
    t0 = time.perf_counter()
    synth = UnitarySynthesizer(
        primitives=sequence.gates,
        subspace=subspace,
        objectives=MultiObjective(fidelity_weight=1.0, leakage_weight=0.05),
        leakage_penalty=LeakagePenalty(weight=0.05),
        execution=ExecutionOptions(engine="auto", use_fast_path=True),
    )
    result = synth.fit(
        target=target,
        init_guess="heuristic",
        multistart=multistart,
        maxiter=maxiter,
    )
    dt = time.perf_counter() - t0
    F = subspace_unitary_fidelity(
        result.simulation.subspace_operator,
        U_target, gauge="global",
    )
    log(f"    F_proj = {F:.6f}   objective = {result.objective:.6f}  ({dt:.1f}s)")
    return {
        "label": label,
        "fidelity": float(F),
        "objective": float(result.objective),
        "success": bool(result.success),
        "elapsed_s": float(dt),
        "n_cav": n_cav,
        "sequence": result.sequence.serialize(),
    }


def build_strategy_B(n_blocks: int, drift: DriftPhaseModel, n_cav: int) -> GateSequence:
    gates = []
    for i in range(n_blocks):
        gates.append(Displacement(name=f"D{i}", alpha=0.3+0j, duration=200e-9))
        gates.append(SQR(name=f"S{2*i}", theta_n=[0.0]*n_cav,
                         phi_n=[0.0]*n_cav, drift_model=drift, duration=400e-9))
        gates.append(ConditionalPhaseSQR(name=f"CP{i}", phases_n=[0.0]*n_cav,
                         drift_model=drift, duration=200e-9))
        gates.append(SQR(name=f"S{2*i+1}", theta_n=[0.0]*n_cav,
                         phi_n=[0.0]*n_cav, drift_model=drift, duration=400e-9))
    gates.append(Displacement(name=f"D{n_blocks}", alpha=0.3+0j, duration=200e-9))
    return GateSequence(gates=gates, n_cav=n_cav)


def build_strategy_D(n_blocks: int, n_cav: int) -> GateSequence:
    gates = []
    for i in range(n_blocks):
        gates.append(Displacement(name=f"D{i}", alpha=0.3+0j, duration=200e-9))
        gates.append(QubitRotation(name=f"R{i}", theta=np.pi/2, phi=0.0, duration=100e-9))
        gates.append(FreeEvolveCondPhase(
            name=f"FE{i}", duration=200e-9,
            drift_model=DriftPhaseModel(chi=abs(CHI), chi2=0.0, kerr=0.0),
            optimize_time=True,
        ))
    gates.append(Displacement(name=f"D{n_blocks}", alpha=0.3+0j, duration=200e-9))
    gates.append(QubitRotation(name=f"R{n_blocks}", theta=np.pi/2, phi=np.pi/2, duration=100e-9))
    return GateSequence(gates=gates, n_cav=n_cav)


# ═══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    # Clear log
    LOG_FILE.write_text("", encoding="utf-8")

    log("=" * 70)
    log("  Iteration 3: Model-Based + Convergence + FE Analysis (efficient)")
    log("=" * 70)
    results = {}
    t_start = time.perf_counter()

    # ── Expt 1a: Strategy B model-based, N_cav=8 ──
    log("\n== Expt 1a: Strategy B (D+SQR+CP, 2 blocks) model-based, N_cav=8 ==")
    seq_B = build_strategy_B(2, full_drift, 8)
    res_B = run_synthesis(seq_B, "B_model_2blk_ncav8", n_cav=8, multistart=6, maxiter=400)
    res_B["strategy"] = "D+SQR+CP (model-based)"
    res_B["n_blocks"] = 2
    results["B_model_2blk_ncav8"] = res_B

    # ── Expt 1b: Strategy D model-based, N_cav=8 ──
    log("\n== Expt 1b: Strategy D (D+R+FE, 2 blocks) model-based, N_cav=8 ==")
    seq_D = build_strategy_D(2, 8)
    res_D = run_synthesis(seq_D, "D_model_2blk_ncav8", n_cav=8, multistart=6, maxiter=400)
    res_D["strategy"] = "D+R+FE (model-based)"
    res_D["n_blocks"] = 2
    results["D_model_2blk_ncav8"] = res_D

    # ── Expt 4: FE wait-time analysis ──
    log("\n== Expt 4: FE wait-time extraction ==")
    fe_times = []
    for gate_data in res_D["sequence"]:
        if gate_data.get("type") == "FreeEvolveCondPhase":
            dur_s = float(gate_data["duration"])
            dur_ns = dur_s * 1e9
            ratio = dur_s / TAU_CZ
            fe_times.append({
                "name": gate_data["name"],
                "duration_s": dur_s,
                "duration_ns": dur_ns,
                "tau_cz_ratio": ratio,
            })
            log(f"    {gate_data['name']}: {dur_ns:.1f} ns  ({ratio:.3f} x tau_CZ)")
    results["fe_wait_time_analysis"] = {
        "tau_cz_ns": float(TAU_CZ * 1e9),
        "chi_MHz": -2.84,
        "gates": fe_times,
    }

    # ── Expt 3: Hilbert space convergence ──
    log("\n== Expt 3: Hilbert space convergence (N_cav=12) ==")
    convergence = {"8": res_B["fidelity"]}

    seq_B12 = build_strategy_B(2, full_drift, 12)
    res_B12 = run_synthesis(seq_B12, "B_model_2blk_ncav12", n_cav=12, multistart=4, maxiter=300)
    convergence["12"] = res_B12["fidelity"]
    results["B_model_2blk_ncav12"] = res_B12

    results["hilbert_convergence"] = convergence

    # ── Save ──
    t_total = time.perf_counter() - t_start
    log(f"\nTotal wall time: {t_total:.1f}s ({t_total/60:.1f} min)")
    log("\n" + "=" * 70)
    log("  SUMMARY")
    log("=" * 70)
    log(f"  B model-based (N_cav=8):  F = {res_B['fidelity']:.6f}  ({res_B['elapsed_s']:.1f}s)")
    log(f"  D model-based (N_cav=8):  F = {res_D['fidelity']:.6f}  ({res_D['elapsed_s']:.1f}s)")
    log(f"  B convergence (N_cav=12): F = {res_B12['fidelity']:.6f}")
    log(f"  tau_CZ = {TAU_CZ*1e9:.1f} ns")
    if fe_times:
        for ft in fe_times:
            log(f"  {ft['name']}: {ft['duration_ns']:.1f} ns ({ft['tau_cz_ratio']:.3f} x tau_CZ)")
    log("=" * 70)

    out_path = DATA_DIR / "iteration3_model_based.json"
    out_path.write_text(json.dumps(results, indent=2, default=str), encoding="utf-8")
    log(f"Results saved to {out_path}")
    log("Done.")
