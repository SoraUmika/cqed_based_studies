"""Iteration 4 v2: Focused bounded-displacement optimization at N_cav=12.

Lean version with reduced sweep for faster results:
  - Strategy B: amp=[0.3, 0.5, 1.0], blocks=[2, 3]
  - Strategy D: amp=[0.5], blocks=[3]
  - multistart=3, maxiter=300 (faster convergence check)
  - Plus Wigner functions for best result

Output:
  data/iteration4_results.json
  figures/bounded_displacement_sweep.{png,pdf}
  figures/wigner_comparison.{png,pdf}
  figures/iteration4_ranking.{png,pdf}
  artifacts/best_strategy_B.json
  artifacts/best_strategy_D.json
"""
from __future__ import annotations

import json
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── paths ─────────────────────────────────────────────────────────────────
STUDY_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = STUDY_ROOT.parents[1]
DATA_DIR = STUDY_ROOT / "data"
FIG_DIR = STUDY_ROOT / "figures"
ARTIFACTS_DIR = STUDY_ROOT / "artifacts"
for d in (DATA_DIR, FIG_DIR, ARTIFACTS_DIR):
    d.mkdir(parents=True, exist_ok=True)

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
print("Importing cqed_sim...", flush=True)
from cqed_sim.unitary_synthesis import (
    Displacement, QubitRotation, SQR, ConditionalPhaseSQR,
    FreeEvolveCondPhase, Subspace, TargetUnitary,
    UnitarySynthesizer, GateSequence, DriftPhaseModel,
    LeakagePenalty, MultiObjective, ExecutionOptions,
    SynthesisConstraints,
    subspace_unitary_fidelity,
)
from cqed_sim.unitary_synthesis.targets import make_target
print("Imports complete.", flush=True)

# ── constants ─────────────────────────────────────────────────────────────
TWO_PI = 2.0 * np.pi
CHI   = TWO_PI * (-2.84e6)
CHIP  = TWO_PI * (-21e3)
KERR  = TWO_PI * (-28e3)

N_CAV = 12
FULL_DIM = 2 * N_CAV  # 24

LOGICAL_LABELS = ["|g,0>", "|g,1>", "|e,0>", "|e,1>"]
subspace = Subspace.custom(FULL_DIM, [0, 1, N_CAV, N_CAV + 1], LOGICAL_LABELS)
no_drift = DriftPhaseModel(chi=0.0, chi2=0.0, kerr=0.0)

COLORS = ['#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB']

# ── target ────────────────────────────────────────────────────────────────
U_target = make_target("cluster", n_match=1)
target = TargetUnitary(U_target, ignore_global_phase=True)
print(f"Target: {U_target.shape}", flush=True)

results = {}

def run_opt(seq, label, max_amp=None, multistart=3, maxiter=300, dur_wt=0.0):
    """Run a single bounded synthesis and return dict."""
    print(f"\n  [{label}] amp<={max_amp} ms={multistart} mi={maxiter}", flush=True)
    cst = SynthesisConstraints(max_amplitude=max_amp) if max_amp else None
    obj = MultiObjective(fidelity_weight=1.0, leakage_weight=0.05, duration_weight=dur_wt)
    t0 = time.perf_counter()
    try:
        synth = UnitarySynthesizer(
            primitives=seq.gates,
            subspace=subspace,
            objectives=obj,
            leakage_penalty=LeakagePenalty(weight=0.05),
            synthesis_constraints=cst,
            execution=ExecutionOptions(engine="auto", use_fast_path=True),
        )
        res = synth.fit(target=target, init_guess="heuristic",
                        multistart=multistart, maxiter=maxiter)
        dt = time.perf_counter() - t0
        F = subspace_unitary_fidelity(res.simulation.subspace_operator,
                                       U_target, gauge="global")
        print(f"    F={F:.6f}  obj={res.objective:.6f}  ({dt:.0f}s)", flush=True)

        disp_amps = []
        for g in res.sequence.gates:
            if hasattr(g, 'alpha'):
                disp_amps.append(float(abs(g.alpha)))

        return {"label": label, "fidelity": float(F), "objective": float(res.objective),
                "success": bool(res.success), "elapsed_s": float(dt),
                "max_amp": max_amp, "disp_amps": disp_amps,
                "sequence": res.sequence.serialize(), "_synth_result": res}
    except Exception as e:
        dt = time.perf_counter() - t0
        print(f"    FAILED: {e} ({dt:.0f}s)", flush=True)
        traceback.print_exc()
        return {"label": label, "fidelity": 0.0, "success": False,
                "elapsed_s": float(dt), "error": str(e), "max_amp": max_amp}


# ═══════════════════════════════════════════════════════════════════════════
# Strategy B (D + SQR + CP): bounded displacement sweep
# ═══════════════════════════════════════════════════════════════════════════
print("\n=== Strategy B: D + SQR + CP ===", flush=True)

def make_B(nb, nc):
    gates = []
    th = [0.0]*nc; th[0] = np.pi/2
    if nc > 1: th[1] = np.pi/4
    ph = [0.0]*nc
    cp_ph = [0.0]*nc
    for i in range(nb):
        gates.append(Displacement(name=f"D{i}", alpha=0.3+0j, duration=200e-9))
        gates.append(SQR(name=f"S{2*i}", theta_n=th, phi_n=ph,
                         drift_model=no_drift, duration=400e-9))
        gates.append(ConditionalPhaseSQR(name=f"CP{i}", phases_n=cp_ph,
                         drift_model=no_drift, duration=200e-9))
        gates.append(SQR(name=f"S{2*i+1}", theta_n=th, phi_n=ph,
                         drift_model=no_drift, duration=400e-9))
    gates.append(Displacement(name=f"D{nb}", alpha=0.3+0j, duration=200e-9))
    return GateSequence(gates=gates, n_cav=nc)

B_CONFIGS = [
    (0.3, 2), (0.3, 3),
    (0.5, 2), (0.5, 3),
    (1.0, 2), (1.0, 3),
]
for amp, nb in B_CONFIGS:
    seq = make_B(nb, N_CAV)
    lbl = f"B_amp{amp}_blk{nb}"
    r = run_opt(seq, lbl, max_amp=amp)
    r["strategy"] = "D+SQR+CP"; r["n_blocks"] = nb; r["n_cav"] = N_CAV
    results[lbl] = r

# ═══════════════════════════════════════════════════════════════════════════
# Strategy D (D + R + FE): bounded displacement sweep
# ═══════════════════════════════════════════════════════════════════════════
print("\n=== Strategy D: D + R + FE ===", flush=True)

def make_D(nb, nc):
    gates = []
    for i in range(nb):
        gates.append(Displacement(name=f"D{i}", alpha=0.3+0j, duration=200e-9))
        gates.append(QubitRotation(name=f"R{i}", theta=np.pi/2, phi=0.0, duration=100e-9))
        gates.append(FreeEvolveCondPhase(
            name=f"FE{i}", duration=200e-9,
            drift_model=DriftPhaseModel(chi=abs(CHI), chi2=0.0, kerr=0.0),
            optimize_time=True))
    gates.append(Displacement(name=f"D{nb}", alpha=0.3+0j, duration=200e-9))
    gates.append(QubitRotation(name=f"R{nb}", theta=np.pi/2, phi=np.pi/2, duration=100e-9))
    return GateSequence(gates=gates, n_cav=nc)

D_CONFIGS = [
    (0.5, 2), (0.5, 3),
    (1.0, 3),
]
for amp, nb in D_CONFIGS:
    seq = make_D(nb, N_CAV)
    lbl = f"D_amp{amp}_blk{nb}"
    r = run_opt(seq, lbl, max_amp=amp)
    r["strategy"] = "D+R+FE"; r["n_blocks"] = nb; r["n_cav"] = N_CAV
    results[lbl] = r

# ═══════════════════════════════════════════════════════════════════════════
# Also do Strategy B unconstrained (no amplitude bound) for comparison
# ═══════════════════════════════════════════════════════════════════════════
print("\n=== Strategy B unconstrained (N_cav=12) ===", flush=True)
for nb in [2, 3]:
    seq = make_B(nb, N_CAV)
    lbl = f"B_unconstrained_blk{nb}"
    r = run_opt(seq, lbl, max_amp=None)
    r["strategy"] = "D+SQR+CP"; r["n_blocks"] = nb; r["n_cav"] = N_CAV
    r["note"] = "No amplitude constraint (baseline)"
    results[lbl] = r


# ═══════════════════════════════════════════════════════════════════════════
# Duration optimization for best result
# ═══════════════════════════════════════════════════════════════════════════
print("\n=== Duration optimization ===", flush=True)
best_key = max((k for k in results if "fidelity" in results[k] and results[k]["fidelity"] > 0.8),
               key=lambda k: results[k]["fidelity"], default=None)
if best_key:
    bv = results[best_key]
    print(f"  Best so far: {best_key} F={bv['fidelity']:.6f}", flush=True)
    nb = bv.get("n_blocks", 3)
    amp = bv.get("max_amp")
    strat = bv.get("strategy", "")
    seq = make_B(nb, N_CAV) if "SQR" in strat else make_D(nb, N_CAV)
    r = run_opt(seq, "best_dur_opt", max_amp=amp, multistart=4, maxiter=400, dur_wt=0.01)
    r["strategy"] = strat; r["n_blocks"] = nb; r["n_cav"] = N_CAV
    r["note"] = "Duration-optimized"
    results["best_dur_opt"] = r


# ═══════════════════════════════════════════════════════════════════════════
# Wigner functions
# ═══════════════════════════════════════════════════════════════════════════
print("\n=== Wigner functions ===", flush=True)
try:
    from cqed_sim.sim.extractors import cavity_wigner, reduced_cavity_state
    import qutip as qt
    HAS_WIGNER = True
except ImportError as e:
    HAS_WIGNER = False
    print(f"  Wigner import failed: {e}", flush=True)

# Find best overall result
ranked = sorted(
    [(k, v["fidelity"]) for k, v in results.items()
     if isinstance(v, dict) and v.get("fidelity", 0) > 0.5],
    key=lambda x: x[1], reverse=True)

print(f"\n  Top results:", flush=True)
for k, f in ranked[:5]:
    print(f"    {k}: F={f:.6f}", flush=True)

wigner_data = {}
if HAS_WIGNER and ranked:
    best_k, best_f = ranked[0]
    sr = results[best_k].get("_synth_result")
    if sr is not None:
        print(f"\n  Computing Wigner for {best_k} (F={best_f:.6f})", flush=True)
        U_full = sr.simulation.full_operator
        U_full_np = U_full.full() if isinstance(U_full, qt.Qobj) else np.asarray(U_full)

        # Target unitary embedded in full space
        U_tgt_full = np.eye(FULL_DIM, dtype=complex)
        lidx = [0, 1, N_CAV, N_CAV + 1]
        for i, li in enumerate(lidx):
            for j, lj in enumerate(lidx):
                U_tgt_full[li, lj] = U_target[i, j]

        basis_labels = ["|g,0>", "|g,1>", "|e,0>", "|e,1>"]
        n_pts = 51; ext = 3.0

        for bi, bl in zip(lidx, basis_labels):
            psi0 = np.zeros(FULL_DIM, dtype=complex)
            psi0[bi] = 1.0
            psi_t = U_tgt_full @ psi0
            psi_a = U_full_np @ psi0

            rho_t = qt.Qobj(np.outer(psi_t, psi_t.conj()), dims=[[2, N_CAV], [2, N_CAV]])
            rho_a = qt.Qobj(np.outer(psi_a, psi_a.conj()), dims=[[2, N_CAV], [2, N_CAV]])

            rho_ct = reduced_cavity_state(rho_t)
            rho_ca = reduced_cavity_state(rho_a)

            xv, yv, Wt = cavity_wigner(rho_ct, n_points=n_pts, extent=ext)
            _,  _,  Wa = cavity_wigner(rho_ca, n_points=n_pts, extent=ext)

            # Cavity state fidelity
            try:
                fcav = float(qt.fidelity(rho_ct, rho_ca)**2)
            except Exception:
                fcav = float(abs(np.trace(rho_ct.full().conj().T @ rho_ca.full())))

            dW = Wt - Wa
            dx = xv[1] - xv[0] if len(xv) > 1 else 1.0
            dy = yv[1] - yv[0] if len(yv) > 1 else 1.0
            l2 = float(np.sqrt(np.sum(dW**2) * dx * dy))
            print(f"    {bl}: F_cav={fcav:.6f}  L2={l2:.6f}", flush=True)

            wigner_data[bl] = {
                "xvec": xv.tolist(), "yvec": yv.tolist(),
                "W_target": Wt.tolist(), "W_achieved": Wa.tolist(),
                "cavity_fidelity": fcav, "l2_distance": l2}

        results["wigner"] = {
            "best_key": best_k, "best_fid": best_f,
            "cavity_fidelities": {bl: wigner_data[bl]["cavity_fidelity"] for bl in basis_labels},
            "l2_distances": {bl: wigner_data[bl]["l2_distance"] for bl in basis_labels}}

        # ── Wigner figure ─────────────────────────────────────────────
        fig, axes = plt.subplots(2, 4, figsize=(14, 7))
        vmax_all = max(
            max(np.max(np.abs(wigner_data[bl]["W_target"])),
                np.max(np.abs(wigner_data[bl]["W_achieved"])))
            for bl in basis_labels)
        for col, bl in enumerate(basis_labels):
            wd = wigner_data[bl]
            X, Y = np.meshgrid(wd["xvec"], wd["yvec"])
            Wt = np.array(wd["W_target"])
            Wa = np.array(wd["W_achieved"])
            axes[0, col].pcolormesh(X, Y, Wt, cmap="RdBu_r",
                                     vmin=-vmax_all, vmax=vmax_all, shading="auto")
            axes[0, col].set_title(f"Target: {bl}", fontsize=9)
            axes[0, col].set_aspect("equal")
            im = axes[1, col].pcolormesh(X, Y, Wa, cmap="RdBu_r",
                                          vmin=-vmax_all, vmax=vmax_all, shading="auto")
            axes[1, col].set_title(f"Achieved: $F_{{cav}}$={wd['cavity_fidelity']:.4f}", fontsize=9)
            axes[1, col].set_aspect("equal")
            axes[1, col].set_xlabel("Re($\\alpha$)")
            if col == 0:
                axes[0, col].set_ylabel("Im($\\alpha$)")
                axes[1, col].set_ylabel("Im($\\alpha$)")
        fig.suptitle(f"Wigner Comparison — {best_k} (F={best_f:.4f})", fontsize=12, fontweight="bold")
        fig.colorbar(im, ax=axes.ravel().tolist(), label="$W(\\alpha)$", shrink=0.8, pad=0.02)
        fig.tight_layout(rect=[0, 0, 0.92, 0.95])
        for fmt in ("png", "pdf"):
            fig.savefig(FIG_DIR / f"wigner_comparison.{fmt}", dpi=300, bbox_inches="tight")
        plt.close(fig)
        print("  Wigner figure saved.", flush=True)


# ═══════════════════════════════════════════════════════════════════════════
# Save
# ═══════════════════════════════════════════════════════════════════════════
print("\n=== Saving results ===", flush=True)

# Clean non-serializable objects
results_json = {}
for k, v in results.items():
    if isinstance(v, dict):
        results_json[k] = {kk: vv for kk, vv in v.items() if kk != "_synth_result"}
    else:
        results_json[k] = v

out = DATA_DIR / "iteration4_results.json"
out.write_text(json.dumps(results_json, indent=2, default=str), encoding="utf-8")
print(f"  Results -> {out}", flush=True)

# Save artifacts
for strat_tag, prefix in [("D+SQR+CP", "best_strategy_B"), ("D+R+FE", "best_strategy_D")]:
    best_e = max(
        ((k, v) for k, v in results.items()
         if isinstance(v, dict) and strat_tag in v.get("strategy", "") and v.get("fidelity", 0) > 0.3),
        key=lambda x: x[1]["fidelity"], default=(None, None))
    if best_e[0]:
        ek, ev = best_e
        sr = ev.get("_synth_result")
        ap = ARTIFACTS_DIR / f"{prefix}.json"
        if sr is not None:
            try:
                sr.save(str(ap), include_history=True)
                print(f"  Artifact: {ap.name} (F={ev['fidelity']:.6f})", flush=True)
            except Exception as e:
                print(f"  Artifact save warning: {e}", flush=True)
                ap.write_text(json.dumps({"label": ek, "fidelity": ev["fidelity"],
                                           "sequence": ev.get("sequence")},
                                          indent=2, default=str), encoding="utf-8")
        else:
            ap.write_text(json.dumps({"label": ek, "fidelity": ev["fidelity"],
                                       "sequence": ev.get("sequence")},
                                      indent=2, default=str), encoding="utf-8")
            print(f"  Fallback artifact: {ap.name}", flush=True)


# ═══════════════════════════════════════════════════════════════════════════
# Summary figures
# ═══════════════════════════════════════════════════════════════════════════
print("\n=== Figures ===", flush=True)

# 1. Bounded sweep: fidelity vs max_amplitude
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
for ax, strat_tag, configs, title in [
    (ax1, "D+SQR+CP", B_CONFIGS, "Strategy B: D + SQR + CP"),
    (ax2, "D+R+FE", D_CONFIGS, "Strategy D: D + R + FE"),
]:
    # Group by n_blocks
    blocks_set = sorted(set(c[1] for c in configs))
    for i, nb in enumerate(blocks_set):
        amps, fids = [], []
        for amp, nb2 in configs:
            if nb2 != nb:
                continue
            k = f"{'B' if 'SQR' in strat_tag else 'D'}_amp{amp}_blk{nb}"
            if k in results and results[k].get("fidelity", 0) > 0:
                amps.append(amp)
                fids.append(results[k]["fidelity"])
        if amps:
            ax.plot(amps, fids, "o-", color=COLORS[i], ms=7, lw=2, label=f"{nb} blocks")
    # Add unconstrained baselines
    for i, nb in enumerate(blocks_set):
        uk = f"B_unconstrained_blk{nb}" if "SQR" in strat_tag else None
        if uk and uk in results:
            ax.axhline(results[uk]["fidelity"], ls=":", color=COLORS[i], alpha=0.5,
                       label=f"{nb} blk (no bound)")
    ax.set_xlabel("Max displacement $|\\alpha|_{\\max}$")
    ax.set_ylabel("Subspace fidelity")
    ax.set_title(title)
    ax.axhline(0.99, ls="--", color="red", alpha=0.3, lw=0.5)
    ax.axhline(0.999, ls="--", color="darkred", alpha=0.3, lw=0.5)
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.3)
fig.suptitle(f"Bounded-Displacement Optimization — $N_{{cav}}={N_CAV}$",
             fontsize=13, fontweight="bold")
fig.tight_layout()
for fmt in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"bounded_displacement_sweep.{fmt}", dpi=300, bbox_inches="tight")
plt.close(fig)

# 2. Ranking bar chart
ranked_all = sorted(
    [(k, v) for k, v in results.items()
     if isinstance(v, dict) and "fidelity" in v and "strategy" in v and v["fidelity"] > 0],
    key=lambda x: x[1]["fidelity"], reverse=True)
if ranked_all:
    fig, ax = plt.subplots(figsize=(10, max(4, len(ranked_all)*0.4)))
    lbs = [s[0] for s in ranked_all]
    fds = [s[1]["fidelity"] for s in ranked_all]
    cs = [COLORS[1] if "SQR" in s[1].get("strategy","") else COLORS[3] for s in ranked_all]
    ax.barh(range(len(lbs)), fds, color=cs, alpha=0.85, edgecolor="black", lw=0.5)
    ax.set_yticks(range(len(lbs)))
    ax.set_yticklabels(lbs, fontsize=8)
    ax.set_xlabel("Subspace Fidelity")
    ax.set_title(f"Iteration 4 Results at $N_{{cav}}={N_CAV}$")
    ax.axvline(0.99, ls="--", color="red", alpha=0.4)
    ax.axvline(0.999, ls=":", color="darkred", alpha=0.4)
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()
    for fmt in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"iteration4_ranking.{fmt}", dpi=300, bbox_inches="tight")
    plt.close(fig)

print("\n=== Summary ===", flush=True)
for k, f in ranked[:8] if ranked else []:
    print(f"  {k}: F={f:.6f}", flush=True)
print("\nDone.", flush=True)
