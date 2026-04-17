"""Multi-decomposition comparison for the cluster-state per-site unitary.

Compares decomposition strategies for the 4x4 cluster-state transfer
matrix U = SWAP . CZ . (H x I):

  Strategy A: D + R + SNAP (prior result — SNAP is cavity-only, F~0.5)
  Strategy B: D + SQR + ConditionalPhaseSQR (structure-guided, depth sweep)
  Strategy C: R + SQR + CP only (no Displacement — controllability probe)
  Strategy D: D + R + FreeEvolveCondPhase (native chi-wait entangler)
  Strategy E: GRAPE full-target (prior results imported)

Also loads cross-study results from hybrid_qubit_cavity_control.

Output: data/decomposition_comparison.json, figures/decomposition_*.{png,pdf}
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

STUDY_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = STUDY_ROOT.parents[1]
DATA_DIR = STUDY_ROOT / "data"
FIG_DIR = STUDY_ROOT / "figures"
ARTIFACTS_DIR = STUDY_ROOT / "artifacts"
DATA_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)
ARTIFACTS_DIR.mkdir(parents=True, exist_ok=True)

STYLE_PATH = WORKSPACE_ROOT / ".github" / "skills" / "publication-figures" / "assets" / "cqed_style.mplstyle"
if STYLE_PATH.exists():
    plt.style.use(str(STYLE_PATH))

SIM_ROOT = Path(
    "C:/Users/dazzl/Box/Shyam Shankar Quantum Circuits Group"
    "/Users/Users_JianJun/cQED_simulation"
)
if str(SIM_ROOT) not in sys.path:
    sys.path.insert(0, str(SIM_ROOT))

# ── cqed_sim imports ──────────────────────────────────────────────────────
from cqed_sim.unitary_synthesis import (
    Displacement, QubitRotation, SQR, ConditionalPhaseSQR,
    FreeEvolveCondPhase, SNAP, Subspace, TargetUnitary,
    UnitarySynthesizer, GateSequence, DriftPhaseModel,
    LeakagePenalty, MultiObjective, ExecutionOptions,
    subspace_unitary_fidelity, simulate_sequence,
)
from cqed_sim.unitary_synthesis.targets import make_target

# ── constants ─────────────────────────────────────────────────────────────
TWO_PI = 2.0 * np.pi
CHI = TWO_PI * (-2.84e6)
CHIP = TWO_PI * (-21e3)
KERR = TWO_PI * (-28e3)
OMEGA_Q = TWO_PI * 6.150e9
OMEGA_C = TWO_PI * 5.241e9
ALPHA = TWO_PI * (-255e6)

# Ideal mode: n_cav=2 for logical subspace, no physical drift
N_CAV = 2
FULL_DIM = 2 * N_CAV

LOGICAL_LABELS = ["|g,0>", "|g,1>", "|e,0>", "|e,1>"]
subspace = Subspace.custom(FULL_DIM, [0, 1, 2, 3], LOGICAL_LABELS)

no_drift = DriftPhaseModel(chi=0.0, chi2=0.0, kerr=0.0)

COLORS = ['#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB']

# ── target unitary ────────────────────────────────────────────────────────
print("=" * 60)
print("  Multi-Decomposition Comparison: Cluster-State U_target")
print("=" * 60)

U_target = make_target("cluster", n_match=1)
print(f"Target shape: {U_target.shape}")
print(np.array2string(U_target, precision=4, suppress_small=True))
target = TargetUnitary(U_target, ignore_global_phase=True)

results = {}


# ── helper ────────────────────────────────────────────────────────────────
def run_synthesis(
    sequence: GateSequence,
    label: str,
    multistart: int = 8,
    maxiter: int = 500,
) -> dict:
    """Run UnitarySynthesizer and return {label, fidelity, objective, ...}."""
    print(f"\n  [{label}] multistart={multistart} maxiter={maxiter}")
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
        U_target,
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


# ═══════════════════════════════════════════════════════════════════════════
# Strategy A: D + R + SNAP (expected failure — SNAP is cavity-only)
# Uses prior result: SNAP = I_q x diag(e^{i theta_n}), cannot entangle.
# ═══════════════════════════════════════════════════════════════════════════
print("\n── A. D + R + SNAP (SNAP is cavity-only — expected F ~ 0.5) ──")
seq_A = GateSequence(
    gates=[
        Displacement(name="D1", alpha=0.3+0j, duration=200e-9),
        QubitRotation(name="R1", theta=np.pi/2, phi=0.0, duration=100e-9),
        SNAP(name="S1", phases=[0.0, np.pi], duration=200e-9),
        Displacement(name="D2", alpha=0.3+0j, duration=200e-9),
        QubitRotation(name="R2", theta=np.pi/2, phi=np.pi/2, duration=100e-9),
        SNAP(name="S2", phases=[0.0, np.pi/2], duration=200e-9),
        Displacement(name="D3", alpha=0.3+0j, duration=200e-9),
    ],
    n_cav=N_CAV,
)
res_A = run_synthesis(seq_A, "A_D_R_SNAP", multistart=8, maxiter=400)
res_A["strategy"] = "D + R + SNAP"
res_A["note"] = "SNAP is cavity-only; cannot produce qubit-cavity entanglement"
results["A_D_R_SNAP"] = res_A


# ═══════════════════════════════════════════════════════════════════════════
# Strategy B: D + SQR + ConditionalPhaseSQR (structure-guided)
# From hybrid study: SQR provides qubit-conditional rotation, Displacement
# provides Fock-level mixing. Test depth sweep.
# Ansatz: [D · SQR · CP · SQR]^n · D
# ═══════════════════════════════════════════════════════════════════════════
print("\n── B. D + SQR + CP depth sweep (Library B+D from hybrid study) ──")

for n_blocks in [1, 2, 3, 4]:
    gates_B = []
    for i in range(n_blocks):
        gates_B.append(Displacement(name=f"D{i}", alpha=0.3+0j, duration=200e-9))
        gates_B.append(SQR(name=f"S{2*i}", theta_n=[0.0, np.pi/2],
                           phi_n=[0.0, 0.0], drift_model=no_drift, duration=400e-9))
        gates_B.append(ConditionalPhaseSQR(name=f"CP{i}", phases_n=[0.0, 0.0],
                           drift_model=no_drift, duration=200e-9))
        gates_B.append(SQR(name=f"S{2*i+1}", theta_n=[0.0, np.pi/2],
                           phi_n=[0.0, 0.0], drift_model=no_drift, duration=400e-9))
    gates_B.append(Displacement(name=f"D{n_blocks}", alpha=0.3+0j, duration=200e-9))
    seq_B = GateSequence(gates=gates_B, n_cav=N_CAV)
    label = f"B_D_SQR_CP_blocks{n_blocks}"
    res_B = run_synthesis(seq_B, label, multistart=8, maxiter=500)
    res_B["strategy"] = "D + SQR + CP"
    res_B["n_blocks"] = n_blocks
    res_B["n_gates"] = len(gates_B)
    results[label] = res_B


# ═══════════════════════════════════════════════════════════════════════════
# Strategy C: R + SQR + CP only (no Displacement — controllability limit)
# SQR preserves Fock-number sectors and cannot mix |n=0> and |n=1>.
# U_target requires such mixing → max F < 1.
# ═══════════════════════════════════════════════════════════════════════════
print("\n── C. SQR + CP only (no Displacement — Fock-mixing limitation) ──")
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
res_C["strategy"] = "R + SQR + CP (no Displacement)"
res_C["note"] = "No Fock-mixing gate; expected F < 1"
results["C_SQR_CP_only"] = res_C


# ═══════════════════════════════════════════════════════════════════════════
# Strategy D: D + R + FreeEvolveCondPhase (native chi-wait entangler)
# FreeEvolveCondPhase uses the dispersive interaction (chi-wait) as the
# entangling primitive instead of SQR/ConditionalPhaseSQR.
# Ansatz: [D · R · FE]^n · D · R
# ═══════════════════════════════════════════════════════════════════════════
print("\n── D. D + R + FreeEvolveCondPhase (native chi-wait) ──")

for n_blocks in [2, 3, 4]:
    gates_D = []
    for i in range(n_blocks):
        gates_D.append(Displacement(name=f"D{i}", alpha=0.3+0j, duration=200e-9))
        gates_D.append(QubitRotation(name=f"R{i}", theta=np.pi/2, phi=0.0, duration=100e-9))
        gates_D.append(FreeEvolveCondPhase(
            name=f"FE{i}", duration=200e-9,
            drift_model=DriftPhaseModel(chi=abs(CHI), chi2=0.0, kerr=0.0),
            optimize_time=True,
        ))
    gates_D.append(Displacement(name=f"D{n_blocks}", alpha=0.3+0j, duration=200e-9))
    gates_D.append(QubitRotation(name=f"R{n_blocks}", theta=np.pi/2, phi=np.pi/2, duration=100e-9))
    seq_D = GateSequence(gates=gates_D, n_cav=N_CAV)
    label = f"D_D_R_FE_blocks{n_blocks}"
    res_D = run_synthesis(seq_D, label, multistart=8, maxiter=500)
    res_D["strategy"] = "D + R + FreeEvolveCondPhase"
    res_D["n_blocks"] = n_blocks
    res_D["n_gates"] = len(gates_D)
    results[label] = res_D


# ═══════════════════════════════════════════════════════════════════════════
# Strategy E: GRAPE (prior sweep results)
# ═══════════════════════════════════════════════════════════════════════════
print("\n── E. GRAPE (prior sweep results) ──")
prior_grape = {
    50: 0.633663, 100: 0.949427, 150: 0.956114,
    200: 0.996610, 300: 0.995730, 400: 0.998963,
}
for dur_ns, fid in sorted(prior_grape.items()):
    results[f"E_GRAPE_{dur_ns}ns"] = {
        "label": f"E_GRAPE_{dur_ns}ns",
        "strategy": "GRAPE",
        "duration_ns": dur_ns,
        "fidelity": fid,
        "source": "prior_sweep",
    }
    print(f"  {dur_ns}ns: F={fid:.6f}")
results["prior_grape_sweep"] = {str(k): v for k, v in prior_grape.items()}

# ═══════════════════════════════════════════════════════════════════════════
# Cross-study comparison: load hybrid_qubit_cavity_control results
# ═══════════════════════════════════════════════════════════════════════════
print("\n─── Cross-Study: hybrid_qubit_cavity_control ───")
HYBRID_DIR = WORKSPACE_ROOT / "studies" / "hybrid_qubit_cavity_control"
cross_study = {}
for subdir, glob_pat in [
    ("utarget_decomposition", "phase4_results.json"),
    ("speed_limit_feasibility", "strategy_summary_refined.json"),
    ("extension_pass_2", "grape_results.json"),
]:
    try:
        fpath = HYBRID_DIR / "data" / subdir / glob_pat
        if fpath.exists():
            raw_data = json.loads(fpath.read_text(encoding="utf-8"))
            if isinstance(raw_data, list):
                for entry in raw_data:
                    lbl = entry.get("label", "unknown")
                    fid = entry.get("strict_fidelity", entry.get("fidelity", 0))
                    cross_study[f"{subdir}_{lbl}"] = {"fidelity": float(fid), "source": subdir}
            elif isinstance(raw_data, dict):
                for lbl, entry in raw_data.items():
                    if isinstance(entry, dict):
                        fid = entry.get("strict_fidelity", entry.get("fidelity", 0))
                        cross_study[f"{subdir}_{lbl}"] = {"fidelity": float(fid), "source": subdir}
            print(f"  Loaded {subdir}: {len(raw_data) if isinstance(raw_data, (list, dict)) else 0} entries")
    except Exception as e:
        print(f"  {subdir} load failed: {e}")
results["cross_study_hybrid"] = cross_study

# ═══════════════════════════════════════════════════════════════════════════
# Save combined results
# ═══════════════════════════════════════════════════════════════════════════
out_path = DATA_DIR / "decomposition_comparison.json"
out_path.write_text(json.dumps(results, indent=2, default=str), encoding="utf-8")
print(f"\nResults saved to {out_path}")

# ═══════════════════════════════════════════════════════════════════════════
# Save target unitary as artifact
# ═══════════════════════════════════════════════════════════════════════════
np.savez(
    str(ARTIFACTS_DIR / "target_unitary.npz"),
    target=U_target,
    study_name="cluster_state_holographic_sim",
    description="4x4 cluster-state per-site MPS isometry: SWAP.CZ.(HxI)",
    load_instructions="np.load(file)['target'] gives the 4x4 complex128 unitary",
)
print("Target unitary saved to artifacts/target_unitary.npz")

# ═══════════════════════════════════════════════════════════════════════════
# Generate figures
# ═══════════════════════════════════════════════════════════════════════════
print("\n─── Generating figures ───")

# Collect strategy results (exclude metadata-only keys)
strat_items = []
for key, val in sorted(results.items()):
    if isinstance(val, dict) and "fidelity" in val and "strategy" in val:
        strat_items.append((key, val))

# Figure 1: Strategy comparison bar chart
fig, ax = plt.subplots(figsize=(10, 5))
labels_bar = [s[0] for s in strat_items]
fids_bar = [s[1]["fidelity"] for s in strat_items]
cat_map = {
    "D + R + SNAP": COLORS[0],
    "D + SQR + CP": COLORS[1],
    "R + SQR + CP (no Displacement)": COLORS[2],
    "D + R + FreeEvolveCondPhase": COLORS[3],
    "GRAPE": COLORS[4],
}
colors_bar = [cat_map.get(s[1]["strategy"], COLORS[6]) for s in strat_items]
ax.bar(range(len(labels_bar)), fids_bar, color=colors_bar, alpha=0.85,
       edgecolor="black", lw=0.5)
ax.set_xticks(range(len(labels_bar)))
ax.set_xticklabels(labels_bar, fontsize=7, rotation=45, ha="right")
ax.set_ylabel("Fidelity")
ax.set_title("Decomposition Strategy Comparison — Cluster-State Target")
ax.axhline(0.99, ls="--", color="red", alpha=0.5, label="99%")
ax.axhline(0.999, ls=":", color="darkred", alpha=0.5, label="99.9%")
ax.legend(fontsize=8)
ax.set_ylim(0, 1.05)
ax.grid(axis="y", alpha=0.3)
fig.tight_layout()
for fmt in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"decomposition_comparison.{fmt}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("  decomposition_comparison done")

# Figure 2: GRAPE fidelity vs duration
all_grape_dns = sorted(prior_grape.keys())
all_grape_fids = [prior_grape[d] for d in all_grape_dns]
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(all_grape_dns, all_grape_fids, "o-", color=COLORS[4], ms=7, lw=2, label="GRAPE (3 seeds)")
ax.axhline(0.999, ls="--", color="red", alpha=0.5, label="99.9%")
ax.axhline(0.99, ls=":", color="orange", alpha=0.5, label="99%")
ax.set_xlabel("Pulse Duration (ns)")
ax.set_ylabel("Fidelity")
ax.set_title("GRAPE Fidelity vs Duration — Cluster-State Target")
ax.legend(fontsize=9)
ax.set_ylim(0.6, 1.005)
ax.grid(True, alpha=0.3)
for fmt in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"grape_fidelity_comparison.{fmt}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("  grape_fidelity_comparison done")

# Figure 3: Strategy ranking (sorted bar chart)
ranked = sorted(strat_items, key=lambda x: -x[1]["fidelity"])
fig, ax = plt.subplots(figsize=(8, max(4, len(ranked) * 0.35)))
labels_r = [r[0] for r in ranked]
fids_r = [r[1]["fidelity"] for r in ranked]
colors_r = [cat_map.get(r[1]["strategy"], COLORS[6]) for r in ranked]
ax.barh(range(len(ranked)), fids_r, color=colors_r, alpha=0.85,
        edgecolor="black", lw=0.5)
ax.set_yticks(range(len(ranked)))
ax.set_yticklabels(labels_r, fontsize=8)
ax.set_xlabel("Fidelity")
ax.set_title("Strategy Ranking by Fidelity")
ax.axvline(0.99, ls="--", color="red", alpha=0.5)
ax.axvline(0.999, ls=":", color="darkred", alpha=0.5)
ax.set_xlim(0, 1.05)
ax.grid(axis="x", alpha=0.3)
ax.invert_yaxis()
fig.tight_layout()
for fmt in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"strategy_ranking.{fmt}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("  strategy_ranking done")

# Figure 4: Depth scaling for D+SQR+CP
depth_keys = [k for k in sorted(results.keys()) if k.startswith("B_D_SQR_CP_blocks")]
if depth_keys:
    blocks = [results[k]["n_blocks"] for k in depth_keys]
    fids_depth = [results[k]["fidelity"] for k in depth_keys]
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(blocks, fids_depth, "s-", color=COLORS[1], ms=8, lw=2)
    ax.set_xlabel("Number of [D·SQR·CP·SQR] blocks")
    ax.set_ylabel("Fidelity")
    ax.set_title("D + SQR + CP: Fidelity vs Depth")
    ax.axhline(0.99, ls="--", color="red", alpha=0.5, label="99%")
    ax.axhline(0.999, ls=":", color="darkred", alpha=0.5, label="99.9%")
    ax.legend(fontsize=9)
    ax.set_ylim(0, 1.05)
    ax.grid(True, alpha=0.3)
    for fmt in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"depth_scaling_D_SQR_CP.{fmt}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  depth_scaling_D_SQR_CP done")

print("\n" + "=" * 60)
print("  SUMMARY")
print("=" * 60)
for key, val in sorted(results.items()):
    if isinstance(val, dict) and "fidelity" in val and "strategy" in val:
        print(f"  {key:35s}  F={val['fidelity']:.6f}  ({val['strategy']})")
print("=" * 60)
print("Done.")
