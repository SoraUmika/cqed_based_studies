"""Iteration 3 — Experiments 1, 3, 4: Model-based verification, Hilbert space
convergence, and FreeEvolveCondPhase wait-time analysis.

Addresses:
  [P1 | MEDIUM] Ideal-mode results must be verified in full dispersive model
  [P2 | LOW]    Increase Hilbert space truncation (convergence check)
  [P2 | MEDIUM] FreeEvolveCondPhase wait-time analysis
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
WORKSPACE_ROOT = STUDY_DIR.parent.parent

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

# ── Physical constants ────────────────────────────────────────────────────────
TWO_PI = 2 * np.pi
CHI = TWO_PI * (-2.84e6)     # rad/s
CHIP = TWO_PI * (-21e3)      # rad/s
KERR = TWO_PI * (-28e3)      # rad/s

# Subspace builder for a given n_cav
def make_subspace(n_cav: int) -> Subspace:
    full_dim = 2 * n_cav
    labels = ["|g,0>", "|g,1>", "|e,0>", "|e,1>"]
    return Subspace.custom(full_dim, [0, 1, n_cav, n_cav + 1], labels)


# ── Target ────────────────────────────────────────────────────────────────────
U_target = make_target("cluster", n_match=1)

# τ_CZ: the CZ-equivalent interaction time for Fock level n=1
TAU_CZ = np.pi / abs(CHI)  # ≈ 176 ns

# ── Drift models ──────────────────────────────────────────────────────────────
no_drift = DriftPhaseModel(chi=0.0, chi2=0.0, kerr=0.0)
full_drift = DriftPhaseModel(chi=CHI, chi2=CHIP, kerr=KERR)

DATA_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)
ARTIFACTS_DIR.mkdir(parents=True, exist_ok=True)


def run_synthesis(
    sequence: GateSequence,
    label: str,
    n_cav: int,
    multistart: int = 12,
    maxiter: int = 800,
) -> dict:
    """Run UnitarySynthesizer and return result dict."""
    subspace = make_subspace(n_cav)
    target = TargetUnitary(U_target, ignore_global_phase=True)
    print(f"\n  [{label}] n_cav={n_cav} multistart={multistart} maxiter={maxiter}")
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
    print(f"    F_proj = {F:.6f}   objective = {result.objective:.6f}  ({dt:.1f}s)")
    return {
        "label": label,
        "fidelity": float(F),
        "objective": float(result.objective),
        "success": bool(result.success),
        "elapsed_s": float(dt),
        "n_cav": n_cav,
        "sequence": result.sequence.serialize(),
    }


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("  Iteration 3: Model-Based Verification + Convergence + FE Analysis")
print("=" * 70)

results = {}

# ═══════════════════════════════════════════════════════════════════════════════
# Experiment 1a: Strategy B (D+SQR+CP, 2 blocks) — model-based, N_cav=8
# ═══════════════════════════════════════════════════════════════════════════════
print("\n── Expt 1a: Strategy B model-based (D+SQR+CP, 2 blocks, N_cav=8) ──")

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

seq_B_model = build_strategy_B(2, full_drift, 8)
res_B_model = run_synthesis(seq_B_model, "B_model_2blk_ncav8", n_cav=8)
res_B_model["strategy"] = "D+SQR+CP (model-based)"
res_B_model["n_blocks"] = 2
res_B_model["drift"] = "full (chi, chi2, kerr)"
results["B_model_2blk_ncav8"] = res_B_model


# ═══════════════════════════════════════════════════════════════════════════════
# Experiment 1b: Strategy D (D+R+FE, 2 blocks) — model-based, N_cav=8
# ═══════════════════════════════════════════════════════════════════════════════
print("\n── Expt 1b: Strategy D model-based (D+R+FE, 2 blocks, N_cav=8) ──")

def build_strategy_D(n_blocks: int, drift: DriftPhaseModel, n_cav: int) -> GateSequence:
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

seq_D_model = build_strategy_D(2, full_drift, 8)
res_D_model = run_synthesis(seq_D_model, "D_model_2blk_ncav8", n_cav=8)
res_D_model["strategy"] = "D+R+FE (model-based)"
res_D_model["n_blocks"] = 2
res_D_model["drift"] = "full (chi for FE)"
results["D_model_2blk_ncav8"] = res_D_model

# Extract FE wait times from the optimized sequence
print("\n── Expt 4: FreeEvolveCondPhase wait-time analysis ──")
fe_times = []
for gate_data in res_D_model["sequence"]:
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
        print(f"    {gate_data['name']}: {dur_ns:.1f} ns  (= {ratio:.3f} × τ_CZ)")

results["fe_wait_time_analysis"] = {
    "tau_cz_ns": TAU_CZ * 1e9,
    "tau_cz_formula": "pi / |chi|",
    "chi_MHz": -2.84,
    "gates": fe_times,
}


# ═══════════════════════════════════════════════════════════════════════════════
# Experiment 3: Hilbert space convergence
# ═══════════════════════════════════════════════════════════════════════════════
print("\n── Expt 3: Hilbert space convergence — Strategy B ──")

convergence = {"8": res_B_model["fidelity"]}

for n_cav_test in [12, 15]:
    print(f"\n  N_cav={n_cav_test}")
    seq_B_conv = build_strategy_B(2, full_drift, n_cav_test)
    res_conv = run_synthesis(
        seq_B_conv, f"B_model_2blk_ncav{n_cav_test}", n_cav=n_cav_test,
        multistart=8, maxiter=600,
    )
    res_conv["strategy"] = "D+SQR+CP (model-based)"
    res_conv["n_blocks"] = 2
    convergence[str(n_cav_test)] = res_conv["fidelity"]
    results[f"B_model_2blk_ncav{n_cav_test}"] = res_conv

print("\n  Hilbert space convergence summary:")
for nc, fid in sorted(convergence.items(), key=lambda x: int(x[0])):
    print(f"    N_cav={nc}: F={fid:.6f}")

results["hilbert_convergence"] = convergence


# ═══════════════════════════════════════════════════════════════════════════════
# Save results
# ═══════════════════════════════════════════════════════════════════════════════
out_path = DATA_DIR / "iteration3_model_based.json"
out_path.write_text(json.dumps(results, indent=2, default=str), encoding="utf-8")
print(f"\nResults saved to {out_path}")

print("\n" + "=" * 70)
print("  SUMMARY")
print("=" * 70)
print(f"  B model-based (2 blocks, N_cav=8):  F = {res_B_model['fidelity']:.6f}")
print(f"  D model-based (2 blocks, N_cav=8):  F = {res_D_model['fidelity']:.6f}")
print(f"  τ_CZ = {TAU_CZ*1e9:.1f} ns (= π/|χ|)")
if fe_times:
    print(f"  FE0 wait time: {fe_times[0]['duration_ns']:.1f} ns ({fe_times[0]['tau_cz_ratio']:.3f} × τ_CZ)")
    if len(fe_times) > 1:
        print(f"  FE1 wait time: {fe_times[1]['duration_ns']:.1f} ns ({fe_times[1]['tau_cz_ratio']:.3f} × τ_CZ)")
print(f"  Convergence: N_cav=8: {convergence.get('8', 'N/A'):.6f}"
      f"  N_cav=12: {convergence.get('12', 'N/A'):.6f}"
      f"  N_cav=15: {convergence.get('15', 'N/A'):.6f}")
print("=" * 70)
print("Done.")
