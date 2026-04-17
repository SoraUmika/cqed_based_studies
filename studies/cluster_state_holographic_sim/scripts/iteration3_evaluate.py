"""Iteration 3 — Model-based evaluation of ideal-mode decompositions.

Instead of re-optimizing from scratch (very slow at N_cav=8), we:
1. Load the ideal-mode optimized parameters from decomposition_comparison.json
2. Reconstruct the gate sequences in N_cav=8 Hilbert space with full drift
3. Propagate with simulate_sequence (fast — no optimization)
4. Compute fidelity to see the impact of the dispersive model

This directly answers the P1 question: "Do ideal-mode results hold up
in the full dispersive model?"

If fidelity degrades, we then try a small refinement optimization.
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
    leakage_metrics,
)
from cqed_sim.unitary_synthesis.targets import make_target

TWO_PI = 2 * np.pi
CHI = TWO_PI * (-2.84e6)
CHIP = TWO_PI * (-21e3)
KERR = TWO_PI * (-28e3)
TAU_CZ = np.pi / abs(CHI)

no_drift = DriftPhaseModel(chi=0.0, chi2=0.0, kerr=0.0)
full_drift = DriftPhaseModel(chi=CHI, chi2=CHIP, kerr=KERR)

U_target = make_target("cluster", n_match=1)

DATA_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)
ARTIFACTS_DIR.mkdir(parents=True, exist_ok=True)


def make_subspace(n_cav: int) -> Subspace:
    full_dim = 2 * n_cav
    return Subspace.custom(full_dim, [0, 1, n_cav, n_cav + 1],
                           ["|g,0>", "|g,1>", "|e,0>", "|e,1>"])


print("=" * 70, flush=True)
print("  Iteration 3: Model-Based Evaluation of Ideal Decompositions", flush=True)
print("=" * 70, flush=True)

# ── Load ideal-mode results ──────────────────────────────────────────────────
comp_path = DATA_DIR / "decomposition_comparison.json"
comp_data = json.loads(comp_path.read_text(encoding="utf-8"))
print(f"Loaded decomposition_comparison.json", flush=True)

results = {}


# ═══════════════════════════════════════════════════════════════════════════════
# Helper: reconstruct a gate sequence from serialized parameters
# ═══════════════════════════════════════════════════════════════════════════════

def rebuild_strategy_B(
    serialized: list[dict],
    drift: DriftPhaseModel,
    n_cav: int,
) -> GateSequence:
    """Rebuild Strategy B from serialized gate parameters.

    The ideal-mode SQR has 2 theta + 2 phi parameters. For n_cav > 2,
    pad the extra Fock levels with zeros.
    """
    gates = []
    for g in serialized:
        dur = float(g["duration"])
        p = g["parameters"]
        if g["type"] == "Displacement":
            gates.append(Displacement(name=g["name"], alpha=complex(p[0], p[1]),
                                      duration=dur))
        elif g["type"] == "SQR":
            n_ideal = len(p) // 2
            theta_n = list(p[:n_ideal]) + [0.0] * (n_cav - n_ideal)
            phi_n = list(p[n_ideal:]) + [0.0] * (n_cav - n_ideal)
            gates.append(SQR(name=g["name"], theta_n=theta_n, phi_n=phi_n,
                             drift_model=drift, duration=dur))
        elif g["type"] == "ConditionalPhaseSQR":
            n_ideal = len(p)
            phases_n = list(p) + [0.0] * (n_cav - n_ideal)
            gates.append(ConditionalPhaseSQR(name=g["name"], phases_n=phases_n,
                                              drift_model=drift, duration=dur))
    return GateSequence(gates=gates, n_cav=n_cav)


def rebuild_strategy_D(
    serialized: list[dict],
    n_cav: int,
) -> GateSequence:
    """Rebuild Strategy D from serialized gate parameters."""
    gates = []
    for g in serialized:
        dur = float(g["duration"])
        p = g["parameters"]
        if g["type"] == "Displacement":
            gates.append(Displacement(name=g["name"], alpha=complex(p[0], p[1]),
                                      duration=dur))
        elif g["type"] == "QubitRotation":
            gates.append(QubitRotation(name=g["name"], theta=p[0], phi=p[1],
                                       duration=dur))
        elif g["type"] == "FreeEvolveCondPhase":
            gates.append(FreeEvolveCondPhase(
                name=g["name"], duration=dur,
                drift_model=DriftPhaseModel(chi=abs(CHI), chi2=0.0, kerr=0.0),
                optimize_time=False,  # Fixed to ideal-mode optimized value
            ))
    return GateSequence(gates=gates, n_cav=n_cav)


# ═══════════════════════════════════════════════════════════════════════════════
# Part 1: Evaluate ideal-mode Strategy B in different settings
# ═══════════════════════════════════════════════════════════════════════════════
print("\n── Strategy B (D+SQR+CP, 2 blocks): Evaluation ──", flush=True)
b2_serial = comp_data["B_D_SQR_CP_blocks2"]["sequence"]

# 1a. Reproduce ideal mode at N_cav=2 (sanity check)
sub2 = make_subspace(2)
seq_B_ideal = rebuild_strategy_B(b2_serial, no_drift, 2)
sim_B_ideal = simulate_sequence(seq_B_ideal, subspace=sub2, backend="ideal")
F_B_ideal = float(subspace_unitary_fidelity(
    sim_B_ideal.subspace_operator, U_target, gauge="global"))
print(f"  B ideal  (N_cav=2, no drift):   F = {F_B_ideal:.6f}  (sanity check)", flush=True)

# 1b. Embed in N_cav=8 with NO drift (check embedding doesn't break it)
sub8 = make_subspace(8)
seq_B_embed = rebuild_strategy_B(b2_serial, no_drift, 8)
sim_B_embed = simulate_sequence(seq_B_embed, subspace=sub8, backend="ideal")
F_B_embed = float(subspace_unitary_fidelity(
    sim_B_embed.subspace_operator, U_target, gauge="global"))
lm_B_embed = leakage_metrics(sim_B_embed.full_operator, sub8)
print(f"  B embed  (N_cav=8, no drift):   F = {F_B_embed:.6f}  leak = {lm_B_embed.average:.6f}", flush=True)

# 1c. Embed in N_cav=8 with FULL drift (the P1 question!)
seq_B_model = rebuild_strategy_B(b2_serial, full_drift, 8)
sim_B_model = simulate_sequence(seq_B_model, subspace=sub8, backend="ideal")
F_B_model = float(subspace_unitary_fidelity(
    sim_B_model.subspace_operator, U_target, gauge="global"))
lm_B_model = leakage_metrics(sim_B_model.full_operator, sub8)
print(f"  B model  (N_cav=8, full drift): F = {F_B_model:.6f}  leak = {lm_B_model.average:.6f}", flush=True)

results["B_ideal_ncav2"] = {"fidelity": F_B_ideal, "n_cav": 2, "drift": "none"}
results["B_embed_ncav8"] = {"fidelity": F_B_embed, "n_cav": 8, "drift": "none",
                            "leakage": float(lm_B_embed.average)}
results["B_model_ncav8"] = {"fidelity": F_B_model, "n_cav": 8, "drift": "full",
                            "leakage": float(lm_B_model.average)}

# 1d. Hilbert space convergence: evaluate at N_cav=12
sub12 = make_subspace(12)
seq_B_12 = rebuild_strategy_B(b2_serial, full_drift, 12)
sim_B_12 = simulate_sequence(seq_B_12, subspace=sub12, backend="ideal")
F_B_12 = float(subspace_unitary_fidelity(
    sim_B_12.subspace_operator, U_target, gauge="global"))
lm_B_12 = leakage_metrics(sim_B_12.full_operator, sub12)
print(f"  B model  (N_cav=12, full drift):F = {F_B_12:.6f}  leak = {lm_B_12.average:.6f}", flush=True)
results["B_model_ncav12"] = {"fidelity": F_B_12, "n_cav": 12, "drift": "full",
                             "leakage": float(lm_B_12.average)}

# 1e. N_cav=15
sub15 = make_subspace(15)
seq_B_15 = rebuild_strategy_B(b2_serial, full_drift, 15)
sim_B_15 = simulate_sequence(seq_B_15, subspace=sub15, backend="ideal")
F_B_15 = float(subspace_unitary_fidelity(
    sim_B_15.subspace_operator, U_target, gauge="global"))
print(f"  B model  (N_cav=15, full drift):F = {F_B_15:.6f}", flush=True)
results["B_model_ncav15"] = {"fidelity": F_B_15, "n_cav": 15, "drift": "full"}

results["B_convergence"] = {
    "2_nodrift": F_B_ideal,
    "8_nodrift": F_B_embed,
    "8_drift": F_B_model,
    "12_drift": F_B_12,
    "15_drift": F_B_15,
}


# ═══════════════════════════════════════════════════════════════════════════════
# Part 2: Evaluate ideal-mode Strategy D in different settings
# ═══════════════════════════════════════════════════════════════════════════════
print("\n── Strategy D (D+R+FE, 2 blocks): Evaluation ──", flush=True)
d2_serial = comp_data["D_D_R_FE_blocks2"]["sequence"]

# 2a. Reproduce ideal mode at N_cav=2
seq_D_ideal = rebuild_strategy_D(d2_serial, 2)
sim_D_ideal = simulate_sequence(seq_D_ideal, subspace=sub2, backend="ideal")
F_D_ideal = float(subspace_unitary_fidelity(
    sim_D_ideal.subspace_operator, U_target, gauge="global"))
print(f"  D ideal  (N_cav=2):  F = {F_D_ideal:.6f}  (sanity check)", flush=True)

# 2b. Embed in N_cav=8
seq_D_embed = rebuild_strategy_D(d2_serial, 8)
sim_D_embed = simulate_sequence(seq_D_embed, subspace=sub8, backend="ideal")
F_D_embed = float(subspace_unitary_fidelity(
    sim_D_embed.subspace_operator, U_target, gauge="global"))
lm_D_embed = leakage_metrics(sim_D_embed.full_operator, sub8)
print(f"  D embed  (N_cav=8):  F = {F_D_embed:.6f}  leak = {lm_D_embed.average:.6f}", flush=True)

# 2c. N_cav=12
seq_D_12 = rebuild_strategy_D(d2_serial, 12)
sim_D_12 = simulate_sequence(seq_D_12, subspace=sub12, backend="ideal")
F_D_12 = float(subspace_unitary_fidelity(
    sim_D_12.subspace_operator, U_target, gauge="global"))
print(f"  D embed  (N_cav=12): F = {F_D_12:.6f}", flush=True)

results["D_ideal_ncav2"] = {"fidelity": F_D_ideal, "n_cav": 2}
results["D_embed_ncav8"] = {"fidelity": F_D_embed, "n_cav": 8,
                            "leakage": float(lm_D_embed.average)}
results["D_embed_ncav12"] = {"fidelity": F_D_12, "n_cav": 12}


# ═══════════════════════════════════════════════════════════════════════════════
# Part 3: FreeEvolveCondPhase wait-time analysis
# ═══════════════════════════════════════════════════════════════════════════════
print("\n── FE wait-time analysis ──", flush=True)
fe_times = []
for g in d2_serial:
    if g["type"] == "FreeEvolveCondPhase":
        dur_s = float(g["duration"])
        dur_ns = dur_s * 1e9
        ratio = dur_s / TAU_CZ
        fe_times.append({
            "name": g["name"],
            "duration_s": dur_s,
            "duration_ns": dur_ns,
            "tau_cz_ratio": ratio,
        })
        print(f"  {g['name']}: {dur_ns:.1f} ns  ({ratio:.3f} x tau_CZ = {TAU_CZ*1e9:.1f} ns)", flush=True)

results["fe_wait_time_analysis"] = {
    "tau_cz_ns": float(TAU_CZ * 1e9),
    "chi_rad_per_s": float(CHI),
    "chi_MHz": -2.84,
    "gates": fe_times,
}


# ═══════════════════════════════════════════════════════════════════════════════
# Part 4: Model-based re-optimization (if fidelity dropped significantly)
# ═══════════════════════════════════════════════════════════════════════════════
REOPT_THRESHOLD = 0.99

need_reopt_B = F_B_model < REOPT_THRESHOLD
need_reopt_D = F_D_embed < REOPT_THRESHOLD

if need_reopt_B or need_reopt_D:
    print(f"\n── Fidelity dropped below {REOPT_THRESHOLD}; attempting re-optimization ──", flush=True)

    if need_reopt_B:
        print(f"  Need re-opt for B (F={F_B_model:.6f} < {REOPT_THRESHOLD})", flush=True)
        # Re-optimize Strategy B with model-based drift at N_cav=8
        # Use only N_cav=2 parameters (embed higher levels as zero)
        n_cav = 8
        gates_reopt = []
        for i in range(2):
            gates_reopt.append(Displacement(name=f"D{i}", alpha=0.3+0j, duration=200e-9))
            gates_reopt.append(SQR(name=f"S{2*i}", theta_n=[0.0, 0.0]+[0.0]*(n_cav-2),
                                   phi_n=[0.0, 0.0]+[0.0]*(n_cav-2),
                                   drift_model=full_drift, duration=400e-9))
            gates_reopt.append(ConditionalPhaseSQR(name=f"CP{i}",
                                   phases_n=[0.0, 0.0]+[0.0]*(n_cav-2),
                                   drift_model=full_drift, duration=200e-9))
            gates_reopt.append(SQR(name=f"S{2*i+1}", theta_n=[0.0, 0.0]+[0.0]*(n_cav-2),
                                   phi_n=[0.0, 0.0]+[0.0]*(n_cav-2),
                                   drift_model=full_drift, duration=400e-9))
        gates_reopt.append(Displacement(name="D2", alpha=0.3+0j, duration=200e-9))
        seq_reopt = GateSequence(gates=gates_reopt, n_cav=n_cav)

        t0 = time.perf_counter()
        synth = UnitarySynthesizer(
            primitives=seq_reopt.gates,
            subspace=sub8,
            objectives=MultiObjective(fidelity_weight=1.0, leakage_weight=0.05),
            leakage_penalty=LeakagePenalty(weight=0.05),
            execution=ExecutionOptions(engine="auto", use_fast_path=True),
        )
        result = synth.fit(
            target=TargetUnitary(U_target, ignore_global_phase=True),
            init_guess="heuristic",
            multistart=4,
            maxiter=300,
        )
        dt = time.perf_counter() - t0
        F_reopt = float(subspace_unitary_fidelity(
            result.simulation.subspace_operator, U_target, gauge="global"))
        print(f"  B re-optimized (N_cav=8, drift): F = {F_reopt:.6f}  ({dt:.1f}s)", flush=True)
        results["B_reopt_ncav8"] = {
            "fidelity": F_reopt,
            "elapsed_s": dt,
            "n_cav": 8,
            "sequence": result.sequence.serialize(),
        }
else:
    print(f"\n  Both strategies above {REOPT_THRESHOLD} threshold — no re-optimization needed.", flush=True)


# ═══════════════════════════════════════════════════════════════════════════════
# Save results
# ═══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 70, flush=True)
print("  SUMMARY", flush=True)
print("=" * 70, flush=True)
print(f"  Strategy B (D+SQR+CP, 2 blocks):", flush=True)
print(f"    Ideal N_cav=2:               F = {F_B_ideal:.6f}", flush=True)
print(f"    Embedded N_cav=8 (no drift):  F = {F_B_embed:.6f}  leak = {lm_B_embed.average:.6f}", flush=True)
print(f"    Model N_cav=8 (full drift):   F = {F_B_model:.6f}  leak = {lm_B_model.average:.6f}", flush=True)
print(f"    Model N_cav=12 (full drift):  F = {F_B_12:.6f}", flush=True)
print(f"    Model N_cav=15 (full drift):  F = {F_B_15:.6f}", flush=True)
print(f"  Strategy D (D+R+FE, 2 blocks):", flush=True)
print(f"    Ideal N_cav=2:               F = {F_D_ideal:.6f}", flush=True)
print(f"    Embedded N_cav=8:            F = {F_D_embed:.6f}  leak = {lm_D_embed.average:.6f}", flush=True)
print(f"    Embedded N_cav=12:           F = {F_D_12:.6f}", flush=True)
print(f"  FE wait times:", flush=True)
for ft in fe_times:
    print(f"    {ft['name']}: {ft['duration_ns']:.1f} ns ({ft['tau_cz_ratio']:.3f} x tau_CZ)", flush=True)
print("=" * 70, flush=True)

out_path = DATA_DIR / "iteration3_model_based.json"
out_path.write_text(json.dumps(results, indent=2, default=str), encoding="utf-8")
print(f"\nResults saved to {out_path}", flush=True)
print("Done.", flush=True)
