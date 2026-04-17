"""Decomposition comparison (streamlined) for cluster-state target.

Incorporates results from the main decomposition_comparison.py run
and tests remaining strategies with reduced multistart/maxiter.

Strategies tested:
  A: D + R + SNAP            (F=0.500, from main run)
  B: D + SQR + CP blocks1    (F=0.707, from main run)
  B: D + SQR + CP blocks2    (F=1.000, from main run)
  C: R + SQR + CP only       (no Displacement)
  D: D + R + FreeEvolveCondPhase (native chi-wait)
  E: GRAPE at multiple durations (prior sweep)
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

from cqed_sim.unitary_synthesis import (
    Displacement, QubitRotation, SQR, ConditionalPhaseSQR,
    FreeEvolveCondPhase, SNAP, Subspace, TargetUnitary,
    UnitarySynthesizer, GateSequence, DriftPhaseModel,
    LeakagePenalty, MultiObjective, ExecutionOptions,
    subspace_unitary_fidelity,
)
from cqed_sim.unitary_synthesis.targets import make_target

TWO_PI = 2.0 * np.pi
CHI = TWO_PI * (-2.84e6)

N_CAV = 2
FULL_DIM = 2 * N_CAV
LOGICAL_LABELS = ["|g,0>", "|g,1>", "|e,0>", "|e,1>"]
subspace = Subspace.custom(FULL_DIM, [0, 1, 2, 3], LOGICAL_LABELS)
no_drift = DriftPhaseModel(chi=0.0, chi2=0.0, kerr=0.0)

U_target = make_target("cluster", n_match=1)
target = TargetUnitary(U_target, ignore_global_phase=True)

COLORS = ['#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB']

results = {}


def run_synthesis(seq, label, multistart=4, maxiter=300):
    print(f"  [{label}] multistart={multistart} maxiter={maxiter} ...", end="", flush=True)
    t0 = time.perf_counter()
    synth = UnitarySynthesizer(
        primitives=seq.gates,
        subspace=subspace,
        objectives=MultiObjective(fidelity_weight=1.0, leakage_weight=0.05),
        leakage_penalty=LeakagePenalty(weight=0.05),
        execution=ExecutionOptions(engine="auto", use_fast_path=True),
    )
    result = synth.fit(target=target, init_guess="heuristic",
                       multistart=multistart, maxiter=maxiter)
    dt = time.perf_counter() - t0
    F = subspace_unitary_fidelity(
        result.simulation.subspace_operator, U_target, gauge="global")
    print(f" F={F:.6f} ({dt:.1f}s)")
    return {
        "label": label, "fidelity": float(F),
        "objective": float(result.objective),
        "success": bool(result.success), "elapsed_s": float(dt),
        "sequence": result.sequence.serialize(),
    }


# Pre-computed results from the main run
print("== Using precomputed results from main run ==")
results["A_D_R_SNAP"] = {
    "label": "A_D_R_SNAP", "strategy": "D + R + SNAP",
    "fidelity": 0.500000, "objective": 0.500000,
    "note": "SNAP is cavity-only; cannot produce qubit-cavity entanglement",
    "elapsed_s": 17.9,
}
results["B_D_SQR_CP_blocks1"] = {
    "label": "B_D_SQR_CP_blocks1", "strategy": "D + SQR + CP",
    "fidelity": 0.707107, "objective": 0.292893,
    "n_blocks": 1, "n_gates": 5, "elapsed_s": 33.5,
}
results["B_D_SQR_CP_blocks2"] = {
    "label": "B_D_SQR_CP_blocks2", "strategy": "D + SQR + CP",
    "fidelity": 1.000000, "objective": 0.000000,
    "n_blocks": 2, "n_gates": 9, "elapsed_s": 70.8,
}
for k, v in results.items():
    print(f"  {k}: F={v['fidelity']:.6f}")


# Strategy C: R + SQR + CP only (no Displacement)
print("\n== C: SQR + CP only (no Displacement) ==")
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
res_C = run_synthesis(seq_C, "C_SQR_CP_only", multistart=6, maxiter=400)
res_C["strategy"] = "R + SQR + CP (no D)"
res_C["note"] = "No Fock-mixing gate; expected F < 1"
results["C_SQR_CP_only"] = res_C


# Strategy D: D + R + FreeEvolveCondPhase
print("\n== D: D + R + FreeEvolveCondPhase (native chi-wait) ==")
for n_blocks in [2, 3]:
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
    res_D = run_synthesis(seq_D, label, multistart=4, maxiter=400)
    res_D["strategy"] = "D + R + FreeEvolve"
    res_D["n_blocks"] = n_blocks
    res_D["n_gates"] = len(gates_D)
    results[label] = res_D


# Strategy E: GRAPE (prior sweep results)
print("\n== E: GRAPE (prior sweep results) ==")
prior_grape = {
    50: 0.633663, 100: 0.949427, 150: 0.956114,
    200: 0.996610, 300: 0.995730, 400: 0.998963,
}
for dur_ns, fid in sorted(prior_grape.items()):
    results[f"E_GRAPE_{dur_ns}ns"] = {
        "label": f"E_GRAPE_{dur_ns}ns", "strategy": "GRAPE",
        "duration_ns": dur_ns, "fidelity": fid, "source": "prior_sweep",
    }
    print(f"  {dur_ns}ns: F={fid:.6f}")


# Cross-study comparison
print("\n== Cross-study: hybrid_qubit_cavity_control ==")
HYBRID_DIR = WORKSPACE_ROOT / "studies" / "hybrid_qubit_cavity_control"
cross_study = {}
for subdir, filename in [
    ("utarget_decomposition", "phase4_results.json"),
    ("speed_limit_feasibility", "strategy_summary_refined.json"),
    ("extension_pass_2", "grape_results.json"),
]:
    try:
        fpath = HYBRID_DIR / "data" / subdir / filename
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
            print(f"  {subdir}: loaded {len(raw_data) if isinstance(raw_data, (list, dict)) else 0} entries")
    except Exception as e:
        print(f"  {subdir}: failed - {e}")
results["cross_study_hybrid"] = cross_study


# Save results
out_path = DATA_DIR / "decomposition_comparison.json"
out_path.write_text(json.dumps(results, indent=2, default=str), encoding="utf-8")
print(f"\nResults saved to {out_path}")

# Save target unitary
np.savez(
    str(ARTIFACTS_DIR / "target_unitary.npz"),
    target=U_target,
    study_name="cluster_state_holographic_sim",
    description="4x4 cluster-state per-site MPS isometry: SWAP.CZ.(HxI)",
)
print("Target unitary saved to artifacts/target_unitary.npz")


# == Figures ==
print("\n== Generating figures ==")

# Collect strategy items
strat_items = []
for key, val in sorted(results.items()):
    if isinstance(val, dict) and "fidelity" in val and "strategy" in val:
        strat_items.append((key, val))

cat_map = {
    "D + R + SNAP": COLORS[0],
    "D + SQR + CP": COLORS[1],
    "R + SQR + CP (no D)": COLORS[2],
    "D + R + FreeEvolve": COLORS[3],
    "GRAPE": COLORS[4],
}

# Fig 1: Bar chart
fig, ax = plt.subplots(figsize=(10, 5))
labels = [s[0] for s in strat_items]
fids = [s[1]["fidelity"] for s in strat_items]
colors = [cat_map.get(s[1]["strategy"], COLORS[6]) for s in strat_items]
ax.bar(range(len(labels)), fids, color=colors, alpha=0.85, edgecolor="black", lw=0.5)
ax.set_xticks(range(len(labels)))
ax.set_xticklabels(labels, fontsize=7, rotation=45, ha="right")
ax.set_ylabel("Fidelity")
ax.set_title("Decomposition Strategy Comparison - Cluster-State Target")
ax.axhline(0.99, ls="--", color="red", alpha=0.5, label="99%")
ax.axhline(0.999, ls=":", color="darkred", alpha=0.5, label="99.9%")
ax.legend(fontsize=8)
ax.set_ylim(0, 1.05)
ax.grid(axis="y", alpha=0.3)
fig.tight_layout()
for fmt in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"decomposition_comparison.{fmt}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("  decomposition_comparison")

# Fig 2: GRAPE vs duration
dns = sorted(prior_grape.keys())
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(dns, [prior_grape[d] for d in dns], "o-", color=COLORS[4], ms=7, lw=2, label="GRAPE (3 seeds)")
ax.axhline(0.999, ls="--", color="red", alpha=0.5, label="99.9%")
ax.axhline(0.99, ls=":", color="orange", alpha=0.5, label="99%")
ax.set_xlabel("Pulse Duration (ns)")
ax.set_ylabel("Fidelity")
ax.set_title("GRAPE Fidelity vs Duration")
ax.legend(fontsize=9); ax.grid(True, alpha=0.3); ax.set_ylim(0.6, 1.005)
for fmt in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"grape_fidelity_comparison.{fmt}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("  grape_fidelity_comparison")

# Fig 3: Strategy ranking
ranked = sorted(strat_items, key=lambda x: -x[1]["fidelity"])
fig, ax = plt.subplots(figsize=(8, max(4, len(ranked) * 0.35)))
lbl_r = [r[0] for r in ranked]
fid_r = [r[1]["fidelity"] for r in ranked]
col_r = [cat_map.get(r[1]["strategy"], COLORS[6]) for r in ranked]
ax.barh(range(len(ranked)), fid_r, color=col_r, alpha=0.85, edgecolor="black", lw=0.5)
ax.set_yticks(range(len(ranked)))
ax.set_yticklabels(lbl_r, fontsize=8)
ax.set_xlabel("Fidelity")
ax.set_title("Strategy Ranking by Fidelity")
ax.axvline(0.99, ls="--", color="red", alpha=0.5)
ax.axvline(0.999, ls=":", color="darkred", alpha=0.5)
ax.set_xlim(0, 1.05); ax.grid(axis="x", alpha=0.3); ax.invert_yaxis()
fig.tight_layout()
for fmt in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"strategy_ranking.{fmt}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("  strategy_ranking")

# Fig 4: Depth scaling D+SQR+CP
depth_keys = [k for k in sorted(results.keys()) if k.startswith("B_D_SQR_CP_blocks")]
if depth_keys:
    blocks = [results[k].get("n_blocks", int(k.split("blocks")[1])) for k in depth_keys]
    fids_d = [results[k]["fidelity"] for k in depth_keys]
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.plot(blocks, fids_d, "s-", color=COLORS[1], ms=8, lw=2)
    ax.set_xlabel("Number of [D.SQR.CP.SQR] blocks")
    ax.set_ylabel("Fidelity")
    ax.set_title("D + SQR + CP: Fidelity vs Depth")
    ax.axhline(0.99, ls="--", color="red", alpha=0.5, label="99%")
    ax.axhline(0.999, ls=":", color="darkred", alpha=0.5, label="99.9%")
    ax.legend(fontsize=9); ax.grid(True, alpha=0.3); ax.set_ylim(0, 1.05)
    for fmt in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"depth_scaling_D_SQR_CP.{fmt}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  depth_scaling_D_SQR_CP")

# Fig 5: Cross-strategy summary (best per family)
print("\n== Saving best-per-family artifact ==")
best_per_family = {}
for key, val in results.items():
    if isinstance(val, dict) and "fidelity" in val and "strategy" in val:
        strat = val["strategy"]
        if strat not in best_per_family or val["fidelity"] > best_per_family[strat]["fidelity"]:
            best_per_family[strat] = {"label": key, "fidelity": val["fidelity"]}

artifact = {
    "study_name": "cluster_state_holographic_sim",
    "description": "Multi-decomposition comparison: best fidelity per strategy family",
    "target": "SWAP.CZ.(HxI) cluster-state 4x4 isometry",
    "best_per_family": best_per_family,
    "full_results_file": "data/decomposition_comparison.json",
}
(ARTIFACTS_DIR / "decomposition_best.json").write_text(
    json.dumps(artifact, indent=2), encoding="utf-8")

# Final summary
print("\n" + "=" * 60)
print("  SUMMARY")
print("=" * 60)
for key, val in sorted(results.items()):
    if isinstance(val, dict) and "fidelity" in val and "strategy" in val:
        print(f"  {key:35s}  F={val['fidelity']:.6f}  ({val['strategy']})")
print("=" * 60)
print("Done.")
