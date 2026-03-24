"""Phase 2 — Idealized logical synthesis in the 4D subspace.

Objectives
----------
Use the cqed_sim UnitarySynthesizer to find short ideal gate decompositions
of U_target in the strict 4D logical space {|g,0>,|g,1>,|e,0>,|e,1>}.

Key finding: U_target requires cavity Fock-level mixing (H_c, CNOT_{q->c}),
which cannot be achieved with SQR + ConditionalPhaseSQR alone.
Displacement gates are required.

Ansatz families tested:
  A. Exact 3-gate factorization (PrimitiveGate): verifies simulator correctness, F=1.
  B. Displacement+SQR+CP (structure-guided): optimizer finds F=1 with physical primitives.
  C. SQR + CP only (controllability probe): max achievable F ~ 0.5 (Fock-mixing limit).
  D. Depth sweep: Library B (D+SQR+CP) at depths 3, 5, 7, 9.

Basis convention: qubit-first {|g,0>,|g,1>,|e,0>,|e,1>} = full-space indices [0,1,N,N+1].
We use Subspace.custom to preserve this ordering (qubit_cavity_block uses a different order).
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── path setup ────────────────────────────────────────────────────────────
STUDY_ROOT = Path(__file__).resolve().parents[2]
COMPONENT_NAME = Path(__file__).resolve().parent.name
SIM_ROOT   = Path("C:/Users/dazzl/Box/Shyam Shankar Quantum Circuits Group"
                  "/Users/Users_JianJun/cQED_simulation")
if str(SIM_ROOT) not in sys.path:
    sys.path.insert(0, str(SIM_ROOT))

from cqed_sim.unitary_synthesis import (
    GateSequence, QubitRotation, SQR, Displacement,
    FreeEvolveCondPhase, ConditionalPhaseSQR, SNAP,
    DriftPhaseModel, Subspace, TargetUnitary,
    LeakagePenalty, MultiObjective, ExecutionOptions,
    UnitarySynthesizer, subspace_unitary_fidelity, simulate_sequence,
    PrimitiveGate,
)

DATA_DIR = STUDY_ROOT / "data" / COMPONENT_NAME
FIG_DIR  = STUDY_ROOT / "figures" / COMPONENT_NAME
DATA_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ── Target ───────────────────────────────────────────────────────────────
s = 1.0 / np.sqrt(2)
U_target_4x4 = np.array([
    [ s, 0,  s, 0],
    [ s, 0, -s, 0],
    [ 0, s,  0, s],
    [ 0,-s,  0, s],
], dtype=np.complex128)

# n_cav = 2 for Phase 2 (minimal logical space); no physical drift
N_CAV = 2
FULL_DIM = 2 * N_CAV

# Logical subspace in qubit-first ordering: {|g,0>=0, |g,1>=1, |e,0>=2, |e,1>=3}
# Using Subspace.custom to preserve this basis order (qubit_cavity_block uses a
# different order: {|g,0>,|e,0>,|g,1>,|e,1>} which would require U_target reordering).
LOGICAL_LABELS = ["|g,0>", "|g,1>", "|e,0>", "|e,1>"]
subspace = Subspace.custom(FULL_DIM, [0, 1, 2, 3], LOGICAL_LABELS)

# target with ignore_global_phase=True since physical gates accumulate phases
target = TargetUnitary(U_target_4x4, ignore_global_phase=True)
target_exact = TargetUnitary(U_target_4x4, ignore_global_phase=False)

no_drift = DriftPhaseModel(chi=0.0, chi2=0.0, kerr=0.0)

results = {}


# ════════════════════════════════════════════════════════════════════════════
# Helper: run UnitarySynthesizer and return best fidelity
# ════════════════════════════════════════════════════════════════════════════
def run_synthesis(
    sequence: GateSequence,
    label: str,
    multistart: int = 8,
    maxiter: int = 500,
    ign_phase: bool = True,
) -> dict:
    print(f"\n  [{label}] multistart={multistart} maxiter={maxiter}")
    t0 = time.perf_counter()
    tgt = target if ign_phase else target_exact
    synth = UnitarySynthesizer(
        primitives=sequence.gates,
        subspace=subspace,
        objectives=MultiObjective(fidelity_weight=1.0, leakage_weight=0.05),
        leakage_penalty=LeakagePenalty(weight=0.05),
        execution=ExecutionOptions(engine="auto", use_fast_path=True),
    )
    result = synth.fit(
        target=tgt,
        init_guess="heuristic",
        multistart=multistart,
        maxiter=maxiter,
    )
    dt = time.perf_counter() - t0

    F = subspace_unitary_fidelity(
        result.simulation.subspace_operator,
        U_target_4x4,
        gauge="global",
    )
    print(f"    F_proj = {F:.6f}   objective = {result.objective:.6f}  ({dt:.1f}s)")
    return {
        "label": label,
        "fidelity": float(F),
        "objective": float(result.objective),
        "success": bool(result.success),
        "elapsed_s": float(dt),
        "sequence": result.sequence.serialize(),
    }


# ════════════════════════════════════════════════════════════════════════════
# A.  Exact 3-gate factorization: verify simulator and subspace ordering
#     U_target = H_c · CNOT_{c→q} · CNOT_{q→c}  (all in qubit-first basis)
# ════════════════════════════════════════════════════════════════════════════
print("\n── A. Exact 3-gate factorization (PrimitiveGate) ────────────────────")

CNOT_q2c_4x4 = np.array([[1,0,0,0],[0,1,0,0],[0,0,0,1],[0,0,1,0]], dtype=np.complex128)
CNOT_c2q_4x4 = np.array([[1,0,0,0],[0,0,0,1],[0,0,1,0],[0,1,0,0]], dtype=np.complex128)
H2 = np.array([[1,1],[1,-1]], dtype=np.complex128) / np.sqrt(2)
Hc_4x4 = np.kron(np.eye(2, dtype=np.complex128), H2)

U_check = Hc_4x4 @ CNOT_c2q_4x4 @ CNOT_q2c_4x4
err_check = np.linalg.norm(U_check - U_target_4x4, ord="fro")
print(f"  Analytic factorization error: {err_check:.2e}")

seq_A_exact = GateSequence(
    gates=[
        PrimitiveGate(name="CNOT_q2c", duration=1e-9, matrix=CNOT_q2c_4x4,
                      optimize_time=False, hilbert_dim=FULL_DIM),
        PrimitiveGate(name="CNOT_c2q", duration=1e-9, matrix=CNOT_c2q_4x4,
                      optimize_time=False, hilbert_dim=FULL_DIM),
        PrimitiveGate(name="IqHc",     duration=1e-9, matrix=Hc_4x4,
                      optimize_time=False, hilbert_dim=FULL_DIM),
    ],
    n_cav=N_CAV,
)

sim_A = simulate_sequence(seq_A_exact, subspace=subspace, backend="ideal")
F_A = subspace_unitary_fidelity(sim_A.subspace_operator, U_target_4x4, gauge="global")
print(f"  Exact 3-gate F = {F_A:.8f}")
results["A_exact_3gate"] = {
    "label": "A_exact_3gate",
    "fidelity": float(F_A),
    "objective": float(1.0 - F_A),
    "success": F_A > 0.9999,
    "notes": "PrimitiveGate analytic matrices; verifies simulator + subspace basis",
}


# ════════════════════════════════════════════════════════════════════════════
# B.  Structure-guided: Displacement + SQR + ConditionalPhaseSQR
#     U_target requires cavity Fock-level mixing (H_c, CNOT_{q->c}) which
#     SQR alone cannot provide; Displacement is the Fock-mixing primitive.
#     Ansatz: D · SQR · CP · SQR · D · SQR · CP · SQR · D (depth-9)
# ════════════════════════════════════════════════════════════════════════════
print("\n── B. D+SQR+CP structure-guided (Library B+D) ───────────────────────")

seq_B = GateSequence(
    gates=[
        Displacement(name="D1", alpha=0.3+0j, duration=200e-9),
        SQR(name="S1", theta_n=[np.pi/2, np.pi], phi_n=[0.0, 0.0],
            drift_model=no_drift, duration=400e-9),
        ConditionalPhaseSQR(name="CP1", phases_n=[0.0, np.pi],
            drift_model=no_drift, duration=200e-9),
        SQR(name="S2", theta_n=[np.pi/2, np.pi], phi_n=[0.0, 0.0],
            drift_model=no_drift, duration=400e-9),
        Displacement(name="D2", alpha=0.3+0j, duration=200e-9),
        SQR(name="S3", theta_n=[np.pi/2, np.pi/2], phi_n=[0.0, 0.0],
            drift_model=no_drift, duration=400e-9),
        ConditionalPhaseSQR(name="CP2", phases_n=[0.0, np.pi],
            drift_model=no_drift, duration=200e-9),
        SQR(name="S4", theta_n=[np.pi/2, np.pi/2], phi_n=[0.0, 0.0],
            drift_model=no_drift, duration=400e-9),
        Displacement(name="D3", alpha=0.3+0j, duration=200e-9),
    ],
    n_cav=N_CAV,
)

res_B = run_synthesis(seq_B, "B_D_SQR_CP_depth9", multistart=16, maxiter=800)
results["B_D_SQR_CP_depth9"] = res_B


# ════════════════════════════════════════════════════════════════════════════
# C.  SQR + CP only (controllability probe)
#     SQR gates preserve Fock-number sectors: they cannot mix |n=0> and |n=1>.
#     U_target requires such mixing; hence max achievable F < 1.
# ════════════════════════════════════════════════════════════════════════════
print("\n── C. SQR + CP only (Fock-mixing limitation) ────────────────────────")

seq_C = GateSequence(
    gates=[
        QubitRotation(name="R1", theta=np.pi/2, phi=0.0, duration=100e-9),
        SQR(name="S1", theta_n=[np.pi/2, np.pi], phi_n=[0.0, 0.0],
            drift_model=no_drift, duration=400e-9),
        ConditionalPhaseSQR(name="CP1", phases_n=[0.0, np.pi],
            drift_model=no_drift, duration=200e-9),
        SQR(name="S2", theta_n=[np.pi/2, np.pi/2], phi_n=[0.0, 0.0],
            drift_model=no_drift, duration=400e-9),
        ConditionalPhaseSQR(name="CP2", phases_n=[0.0, np.pi/2],
            drift_model=no_drift, duration=200e-9),
        QubitRotation(name="R2", theta=np.pi/2, phi=np.pi/2, duration=100e-9),
    ],
    n_cav=N_CAV,
)

res_C = run_synthesis(seq_C, "C_SQR_CP_only", multistart=12, maxiter=600)
results["C_SQR_CP_only"] = res_C
print("  (Expected F < 1.0: SQR+CP cannot mix Fock levels, max ≈ 0.5)")


# ════════════════════════════════════════════════════════════════════════════
# D.  Depth sweep: Library B+D (D+SQR+CP) at varying depths
#     Explore minimum depth needed to achieve F > 0.99
# ════════════════════════════════════════════════════════════════════════════
print("\n── D. Depth sweep (Library B+D) ─────────────────────────────────────")

depth_fidelities = {}
for n_d in [1, 2, 3, 4]:
    n_sqr = n_d + 1         # SQR gates
    n_cp  = n_d             # CP gates
    n_dis = n_d + 1         # Displacement gates (surround SQR)
    print(f"\n  n_displacement={n_dis}, n_sqr={n_sqr}, n_cp={n_cp}:")
    gates = []
    for i in range(n_d):
        gates.append(Displacement(name=f"D{i}", alpha=0.3+0j, duration=200e-9))
        gates.append(SQR(name=f"S{2*i}", theta_n=[0.0, np.pi/2], phi_n=[0.0, 0.0],
                         drift_model=no_drift, duration=400e-9))
        gates.append(ConditionalPhaseSQR(name=f"CP{i}", phases_n=[0.0, 0.0],
                         drift_model=no_drift, duration=200e-9))
        gates.append(SQR(name=f"S{2*i+1}", theta_n=[0.0, np.pi/2], phi_n=[0.0, 0.0],
                         drift_model=no_drift, duration=400e-9))
    gates.append(Displacement(name=f"D{n_d}", alpha=0.3+0j, duration=200e-9))
    seq = GateSequence(gates=gates, n_cav=N_CAV)
    synth = UnitarySynthesizer(
        primitives=seq.gates, subspace=subspace,
        objectives=MultiObjective(fidelity_weight=1.0),
        execution=ExecutionOptions(use_fast_path=True),
    )
    res = synth.fit(target=target, init_guess="heuristic", multistart=8, maxiter=400)
    F_d = subspace_unitary_fidelity(res.simulation.subspace_operator, U_target_4x4, gauge="global")
    depth_fidelities[n_d] = float(F_d)
    print(f"    F = {F_d:.5f}")

results["D_depth_sweep"] = depth_fidelities


# ════════════════════════════════════════════════════════════════════════════
# E.  Convergence plot: fidelity vs depth
# ════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

# Left: depth sweep
ax = axes[0]
depths = sorted(depth_fidelities)
fids   = [depth_fidelities[d] for d in depths]
ax.plot(depths, fids, "o-", color="steelblue", lw=2, markersize=8)
ax.axhline(0.99, ls="--", color="gray", lw=1, label="F=0.99")
ax.axhline(0.999, ls=":", color="gray", lw=1, label="F=0.999")
ax.set_xlabel("Number of D·(SQR·CP·SQR)·D blocks", fontsize=12)
ax.set_ylabel("Projected fidelity $F_{\\rm proj}$", fontsize=12)
ax.set_title("Phase 2: Depth sweep — Library B+D", fontsize=11)
ax.set_xticks(depths)
ax.set_ylim(0, 1.05)
ax.legend(fontsize=10)
ax.grid(True, alpha=0.3)

# Right: comparison of methods
ax = axes[1]
labels_plot = ["A (exact\nfactorization)", "B (D+SQR+CP\ndepth-9)", "C (SQR+CP\nonly)"]
fids_plot  = [results["A_exact_3gate"]["fidelity"],
              results["B_D_SQR_CP_depth9"]["fidelity"],
              results["C_SQR_CP_only"]["fidelity"]]
colors = ["steelblue", "darkorange", "tomato"]
bars = ax.bar(range(3), fids_plot, color=colors, alpha=0.85, edgecolor="k", linewidth=0.8)
ax.axhline(1.0, ls="--", color="gray", lw=1)
ax.set_ylabel("Projected fidelity $F_{\\rm proj}$", fontsize=12)
ax.set_title("Phase 2: Ansatz comparison", fontsize=11)
ax.set_xticks(range(3))
ax.set_xticklabels(labels_plot, fontsize=9)
ax.set_ylim(0, 1.12)
for bar, f in zip(bars, fids_plot):
    ax.text(bar.get_x() + bar.get_width()/2, f + 0.02,
            f"{f:.3f}", ha="center", va="bottom", fontsize=10)
ax.grid(True, alpha=0.3, axis="y")

fig.tight_layout()
fig.savefig(FIG_DIR / "phase2_synthesis.png", dpi=150)
fig.savefig(FIG_DIR / "phase2_synthesis.pdf")
plt.close(fig)


# ════════════════════════════════════════════════════════════════════════════
# Summary
# ════════════════════════════════════════════════════════════════════════════
print("\n══════════ Phase 2 Summary ══════════")
print(f"{'Ansatz':45s} {'F_proj':>10s} {'Success':>10s}")
print("-" * 67)
for key, r in results.items():
    if isinstance(r, dict) and "fidelity" in r:
        label   = r.get("label", key)
        fidelity = r.get("fidelity", 0.0)
        success  = r.get("success", False)
        print(f"{label:45s} {fidelity:10.6f} {str(success):>10s}")

print(f"\nDepth sweep (Library B+D): {depth_fidelities}")

with open(DATA_DIR / "phase2_results.json", "w") as f:
    json.dump(results, f, indent=2, default=str)

print(f"\nPhase 2 results written to {DATA_DIR / 'phase2_results.json'}")
print("Phase 2 complete.")
