"""Final analysis and figure generation for decomposition comparison.

Uses all results obtained from the UnitarySynthesizer runs:
  A: D + R + SNAP           F=0.500  (SNAP is cavity-only)
  B1: D + SQR + CP (1 blk)  F=0.707  (insufficient depth)
  B2: D + SQR + CP (2 blk)  F=1.000  (perfect ideal)
  C: R + SQR + CP (no D)    F=0.500  (no Fock mixing)
  D2: D + R + FE (2 blk)    F=0.9999 (near-perfect)
  E: GRAPE                   F=0.634-0.999 (prior sweep)
"""
from __future__ import annotations

import json
import sys
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

SIM_ROOT = Path(
    "C:/Users/dazzl/Box/Shyam Shankar Quantum Circuits Group"
    "/Users/Users_JianJun/cQED_simulation"
)
if str(SIM_ROOT) not in sys.path:
    sys.path.insert(0, str(SIM_ROOT))
from cqed_sim.unitary_synthesis.targets import make_target

STYLE_PATH = WORKSPACE_ROOT / ".github" / "skills" / "publication-figures" / "assets" / "cqed_style.mplstyle"
if STYLE_PATH.exists():
    plt.style.use(str(STYLE_PATH))

COLORS = {
    "D + R + SNAP":     '#4477AA',
    "D + SQR + CP":     '#EE6677',
    "R + SQR + CP":     '#228833',
    "D + R + FE":       '#CCBB44',
    "GRAPE":            '#66CCEE',
}

# =====================================================================
# Results from UnitarySynthesizer runs
# =====================================================================
results = {
    "A_D_R_SNAP": {
        "label": "D+R+SNAP (7 gates)", "strategy": "D + R + SNAP",
        "fidelity": 0.500000, "elapsed_s": 17.9,
        "note": "SNAP = I_q x diag(e^{i*theta_n}): cavity-only, no entanglement",
    },
    "B_D_SQR_CP_blocks1": {
        "label": "D+SQR+CP (1 block, 5g)", "strategy": "D + SQR + CP",
        "fidelity": 0.707107, "n_blocks": 1, "n_gates": 5, "elapsed_s": 33.5,
    },
    "B_D_SQR_CP_blocks2": {
        "label": "D+SQR+CP (2 blocks, 9g)", "strategy": "D + SQR + CP",
        "fidelity": 1.000000, "n_blocks": 2, "n_gates": 9, "elapsed_s": 70.8,
    },
    "C_SQR_CP_only": {
        "label": "R+SQR+CP (no D, 6g)", "strategy": "R + SQR + CP",
        "fidelity": 0.500000, "elapsed_s": 3.4,
        "note": "SQR preserves Fock sectors; cannot mix |n=0> and |n=1>",
    },
    "D_D_R_FE_blocks2": {
        "label": "D+R+FE (2 blocks, 8g)", "strategy": "D + R + FE",
        "fidelity": 0.999869, "n_blocks": 2, "n_gates": 8, "elapsed_s": 54.1,
    },
}

prior_grape = {50: 0.6337, 100: 0.9494, 150: 0.9561, 200: 0.9966, 300: 0.9957, 400: 0.9990}
for dur_ns, fid in sorted(prior_grape.items()):
    results[f"E_GRAPE_{dur_ns}ns"] = {
        "label": f"GRAPE {dur_ns}ns", "strategy": "GRAPE",
        "duration_ns": dur_ns, "fidelity": fid,
    }

# =====================================================================
# Cross-study: load hybrid_qubit_cavity_control results
# =====================================================================
HYBRID_DIR = WORKSPACE_ROOT / "studies" / "hybrid_qubit_cavity_control"
cross_study = {}
for subdir, fname in [
    ("utarget_decomposition", "phase4_results.json"),
    ("speed_limit_feasibility", "strategy_summary_refined.json"),
]:
    try:
        fpath = HYBRID_DIR / "data" / subdir / fname
        if fpath.exists():
            raw = json.loads(fpath.read_text(encoding="utf-8"))
            if isinstance(raw, list):
                for e in raw:
                    lbl = e.get("label", "unknown")
                    fid = e.get("strict_fidelity", e.get("fidelity", 0))
                    cross_study[f"{subdir}_{lbl}"] = float(fid)
            elif isinstance(raw, dict):
                for lbl, e in raw.items():
                    if isinstance(e, dict):
                        fid = e.get("strict_fidelity", e.get("fidelity", 0))
                        cross_study[f"{subdir}_{lbl}"] = float(fid)
            print(f"  Cross-study: {subdir} loaded")
    except Exception as e:
        print(f"  Cross-study: {subdir} failed: {e}")

# =====================================================================
# Save results
# =====================================================================
out = {k: v for k, v in results.items()}
out["cross_study_hybrid"] = cross_study
out["prior_grape_sweep"] = {str(k): v for k, v in prior_grape.items()}
(DATA_DIR / "decomposition_comparison.json").write_text(
    json.dumps(out, indent=2, default=str), encoding="utf-8")
print("Saved decomposition_comparison.json")

# Save target unitary
U_target = make_target("cluster", n_match=1)
np.savez(str(ARTIFACTS_DIR / "target_unitary.npz"), target=U_target,
         study_name="cluster_state_holographic_sim",
         description="4x4 cluster-state: SWAP.CZ.(HxI)")
print("Saved artifacts/target_unitary.npz")

# Best-per-family artifact
best = {}
for k, v in results.items():
    s = v["strategy"]
    if s not in best or v["fidelity"] > best[s]["fidelity"]:
        best[s] = {"label": k, "fidelity": v["fidelity"]}
(ARTIFACTS_DIR / "decomposition_best.json").write_text(
    json.dumps({"best_per_family": best, "target": "SWAP.CZ.(HxI)"}, indent=2),
    encoding="utf-8")
print("Saved artifacts/decomposition_best.json")

# =====================================================================
# Figure 1: Strategy comparison bar chart
# =====================================================================
discrete_items = [
    ("A_D_R_SNAP", results["A_D_R_SNAP"]),
    ("B_blocks1", results["B_D_SQR_CP_blocks1"]),
    ("B_blocks2", results["B_D_SQR_CP_blocks2"]),
    ("C_no_Disp", results["C_SQR_CP_only"]),
    ("D_FE_blk2", results["D_D_R_FE_blocks2"]),
    ("E_GRAPE_200ns", results["E_GRAPE_200ns"]),
    ("E_GRAPE_400ns", results["E_GRAPE_400ns"]),
]

fig, ax = plt.subplots(figsize=(9, 4.5))
xlabels = [it[0] for it in discrete_items]
fids = [it[1]["fidelity"] for it in discrete_items]
cols = [COLORS.get(it[1]["strategy"], "#BBBBBB") for it in discrete_items]
bars = ax.bar(range(len(xlabels)), fids, color=cols, alpha=0.85,
              edgecolor="black", lw=0.5)
# Add value labels on bars
for bar, f in zip(bars, fids):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.015,
            f"{f:.3f}", ha="center", va="bottom", fontsize=8)
ax.set_xticks(range(len(xlabels)))
ax.set_xticklabels(xlabels, fontsize=8, rotation=30, ha="right")
ax.set_ylabel("Subspace Fidelity")
ax.set_title("Decomposition Strategy Comparison -- Cluster-State Target")
ax.axhline(0.99, ls="--", color="red", alpha=0.5, label="99%")
ax.axhline(0.999, ls=":", color="darkred", alpha=0.5, label="99.9%")
ax.legend(fontsize=8, loc="lower right")
ax.set_ylim(0, 1.12)
ax.grid(axis="y", alpha=0.3)
fig.tight_layout()
for fmt in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"decomposition_comparison.{fmt}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("Fig 1: decomposition_comparison")

# =====================================================================
# Figure 2: GRAPE fidelity vs duration
# =====================================================================
dns = sorted(prior_grape.keys())
fvals = [prior_grape[d] for d in dns]
fig, ax = plt.subplots(figsize=(5.5, 3.8))
ax.plot(dns, fvals, "o-", color=COLORS["GRAPE"], ms=7, lw=2, label="GRAPE (3 seeds)")
ax.axhline(0.999, ls="--", color="red", alpha=0.5, label="F = 99.9%")
ax.axhline(0.99, ls=":", color="orange", alpha=0.5, label="F = 99%")
ax.set_xlabel("Pulse Duration (ns)")
ax.set_ylabel("Fidelity")
ax.set_title("GRAPE Fidelity vs Duration")
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_ylim(0.6, 1.005)
fig.tight_layout()
for fmt in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"grape_fidelity_vs_duration.{fmt}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("Fig 2: grape_fidelity_vs_duration")

# =====================================================================
# Figure 3: Strategy ranking (horizontal bars, sorted)
# =====================================================================
all_items = sorted(
    [(k, v) for k, v in results.items()
     if isinstance(v, dict) and "fidelity" in v and "strategy" in v],
    key=lambda x: -x[1]["fidelity"])
fig, ax = plt.subplots(figsize=(7.5, max(3.5, len(all_items)*0.32)))
lbl_r = [it[1].get("label", it[0]) for it in all_items]
fid_r = [it[1]["fidelity"] for it in all_items]
col_r = [COLORS.get(it[1]["strategy"], "#BBBBBB") for it in all_items]
ax.barh(range(len(all_items)), fid_r, color=col_r, alpha=0.85,
        edgecolor="black", lw=0.4)
ax.set_yticks(range(len(all_items)))
ax.set_yticklabels(lbl_r, fontsize=7)
ax.set_xlabel("Fidelity")
ax.set_title("Strategy Ranking by Fidelity")
ax.axvline(0.99, ls="--", color="red", alpha=0.5)
ax.axvline(0.999, ls=":", color="darkred", alpha=0.5)
ax.set_xlim(0, 1.08)
ax.grid(axis="x", alpha=0.3)
ax.invert_yaxis()
# Add fidelity text
for i, f in enumerate(fid_r):
    ax.text(f + 0.008, i, f"{f:.4f}", va="center", fontsize=7)
fig.tight_layout()
for fmt in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"strategy_ranking.{fmt}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("Fig 3: strategy_ranking")

# =====================================================================
# Figure 4: Depth scaling for D+SQR+CP
# =====================================================================
depth_data = [(1, 0.707107), (2, 1.000000)]
fig, ax = plt.subplots(figsize=(4.5, 3.5))
bx, by = zip(*depth_data)
ax.plot(bx, by, "s-", color=COLORS["D + SQR + CP"], ms=9, lw=2)
ax.set_xlabel("Number of [D . SQR . CP . SQR] blocks")
ax.set_ylabel("Fidelity")
ax.set_title("D + SQR + CP: Fidelity vs Depth")
ax.axhline(0.99, ls="--", color="red", alpha=0.5, label="99%")
ax.axhline(0.999, ls=":", color="darkred", alpha=0.5, label="99.9%")
ax.legend(fontsize=8)
ax.set_ylim(0.4, 1.05)
ax.set_xticks([1, 2])
ax.grid(True, alpha=0.3)
fig.tight_layout()
for fmt in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"depth_scaling_D_SQR_CP.{fmt}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("Fig 4: depth_scaling_D_SQR_CP")

# =====================================================================
# Figure 5: Infidelity comparison (log scale)
# =====================================================================
selected = [
    ("B_blocks2\n(D+SQR+CP)", 1.0 - 1.000000),
    ("D_FE_blk2\n(D+R+FE)", 1.0 - 0.999869),
    ("GRAPE 400ns", 1.0 - 0.9990),
    ("GRAPE 200ns", 1.0 - 0.9966),
    ("GRAPE 100ns", 1.0 - 0.9494),
    ("A (D+R+SNAP)", 1.0 - 0.500),
    ("C (no Disp)", 1.0 - 0.500),
]
# Filter out zero infidelity for log scale
sel_labels = []
sel_inf = []
for lbl, inf in selected:
    if inf > 0:
        sel_labels.append(lbl)
        sel_inf.append(inf)
    else:
        sel_labels.append(lbl)
        sel_inf.append(1e-8)  # placeholder for perfect fidelity

fig, ax = plt.subplots(figsize=(6.5, 4))
bar_cols = [COLORS["D + SQR + CP"], COLORS["D + R + FE"],
            COLORS["GRAPE"], COLORS["GRAPE"], COLORS["GRAPE"],
            COLORS["D + R + SNAP"], COLORS["R + SQR + CP"]]
ax.barh(range(len(sel_labels)), sel_inf, color=bar_cols, alpha=0.85,
        edgecolor="black", lw=0.4)
ax.set_xscale("log")
ax.set_yticks(range(len(sel_labels)))
ax.set_yticklabels(sel_labels, fontsize=8)
ax.set_xlabel("Infidelity (1 - F)")
ax.set_title("Infidelity Comparison (log scale)")
ax.axvline(1e-2, ls="--", color="red", alpha=0.5, label="1%")
ax.axvline(1e-3, ls=":", color="darkred", alpha=0.5, label="0.1%")
ax.legend(fontsize=8)
ax.set_xlim(1e-9, 1.0)
ax.grid(axis="x", alpha=0.3)
ax.invert_yaxis()
fig.tight_layout()
for fmt in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"infidelity_comparison.{fmt}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("Fig 5: infidelity_comparison")

print("\n" + "=" * 60)
print("  FINAL SUMMARY")
print("=" * 60)
print(f"  {'Strategy':<35s}  {'Fidelity':>10s}")
print(f"  {'-'*35}  {'-'*10}")
for k, v in sorted(results.items(), key=lambda x: -x[1].get("fidelity", 0)):
    if "fidelity" in v:
        print(f"  {v.get('label', k):<35s}  {v['fidelity']:>10.6f}")
print("=" * 60)
print("\nKey finding: D+SQR+CP with 2 blocks achieves F=1.000 (ideal)")
print("Key finding: D+R+FreeEvolveCondPhase with 2 blocks achieves F=0.9999")
print("Key finding: SNAP and no-Displacement limited to F=0.500")
print("Key finding: GRAPE at 400ns achieves F=0.999 (model-based)")
print("\nDone.")
