"""Iteration 4 v4: Minimal bounded-displacement study.

Fastest possible path to results:
  1. Sweep max_amplitude at N_cav=4 (dim=8, ultra-fast) for trend.
  2. Validate best at N_cav=6 (dim=12, fast).
  3. Compare to existing GRAPE data from Iter 3.
  4. Wigner comparison using best achievable unitary.
  5. Artifacts.

Physics: With |alpha|<=0.5, P(n>=4)<5e-5, so N_cav=4 is essentially exact.
We confirm convergence by running the best config also at N_cav=6,12.
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

STYLE_PATH = (WORKSPACE_ROOT / ".github" / "skills" / "publication-figures"
              / "assets" / "cqed_style.mplstyle")
if STYLE_PATH.exists():
    plt.style.use(str(STYLE_PATH))

SIM_ROOT = Path(
    "C:/Users/dazzl/Box/Shyam Shankar Quantum Circuits Group"
    "/Users/Users_JianJun/cQED_simulation"
)
if str(SIM_ROOT) not in sys.path:
    sys.path.insert(0, str(SIM_ROOT))

print("Importing cqed_sim...", flush=True)
t_import = time.perf_counter()
from cqed_sim.unitary_synthesis import (
    Displacement, QubitRotation, SQR, ConditionalPhaseSQR,
    FreeEvolveCondPhase, Subspace, TargetUnitary,
    UnitarySynthesizer, GateSequence, DriftPhaseModel,
    LeakagePenalty, MultiObjective, ExecutionOptions,
    SynthesisConstraints,
    subspace_unitary_fidelity,
)
from cqed_sim.unitary_synthesis.targets import make_target
print(f"Import done ({time.perf_counter()-t_import:.1f}s).", flush=True)

# ── constants ─────────────────────────────────────────────────────────────
TWO_PI = 2.0 * np.pi
CHI   = TWO_PI * (-2.84e6)
U_target = make_target("cluster", n_match=1)
target = TargetUnitary(U_target, ignore_global_phase=True)
print(f"Target: {U_target.shape}", flush=True)

no_drift = DriftPhaseModel(chi=0.0, chi2=0.0, kerr=0.0)
results = {}
COLORS = ['#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB']

def make_subspace(n_cav):
    full_dim = 2 * n_cav
    return Subspace.custom(full_dim, [0, 1, n_cav, n_cav + 1],
                           ["|g,0>", "|g,1>", "|e,0>", "|e,1>"])

def make_B(nb, nc):
    """Strategy B: [D · SQR · CP · SQR]^nb · D"""
    gates = []
    th = [0.0]*nc; th[0] = np.pi/2
    if nc > 1: th[1] = np.pi/4
    ph = [0.0]*nc
    for i in range(nb):
        gates.append(Displacement(name=f"D{i}", alpha=0.3+0j, duration=200e-9))
        gates.append(SQR(name=f"S{2*i}", theta_n=th[:], phi_n=ph[:],
                         drift_model=no_drift, duration=400e-9))
        gates.append(ConditionalPhaseSQR(name=f"CP{i}", phases_n=[0.0]*nc,
                         drift_model=no_drift, duration=200e-9))
        gates.append(SQR(name=f"S{2*i+1}", theta_n=th[:], phi_n=ph[:],
                         drift_model=no_drift, duration=400e-9))
    gates.append(Displacement(name=f"D{nb}", alpha=0.3+0j, duration=200e-9))
    return GateSequence(gates=gates, n_cav=nc)

def make_D(nb, nc):
    """Strategy D: [D · R · FE]^nb · D · R"""
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

def run_opt(seq, label, subsp, max_amp=None, multistart=2, maxiter=200, dur_wt=0.0):
    print(f"  [{label}] amp<={max_amp} ms={multistart} mi={maxiter} nc={seq.n_cav}", flush=True)
    cst = SynthesisConstraints(max_amplitude=max_amp) if max_amp else None
    obj = MultiObjective(fidelity_weight=1.0, leakage_weight=0.05, duration_weight=dur_wt)
    t0 = time.perf_counter()
    try:
        synth = UnitarySynthesizer(
            primitives=seq.gates, subspace=subsp, objectives=obj,
            leakage_penalty=LeakagePenalty(weight=0.05),
            synthesis_constraints=cst,
            execution=ExecutionOptions(engine="auto", use_fast_path=True))
        res = synth.fit(target=target, init_guess="heuristic",
                        multistart=multistart, maxiter=maxiter)
        dt = time.perf_counter() - t0
        F = subspace_unitary_fidelity(res.simulation.subspace_operator,
                                       U_target, gauge="global")
        disp_amps = [float(abs(g.alpha)) for g in res.sequence.gates if hasattr(g, 'alpha')]
        print(f"    F={F:.6f}  obj={res.objective:.6f}  ({dt:.0f}s)", flush=True)
        return {"label": label, "fidelity": float(F), "objective": float(res.objective),
                "success": bool(res.success), "elapsed_s": float(dt),
                "max_amp": max_amp, "disp_amps": disp_amps,
                "n_cav": seq.n_cav, "_synth": res}
    except Exception as e:
        dt = time.perf_counter() - t0
        print(f"    FAILED: {e} ({dt:.0f}s)", flush=True)
        traceback.print_exc()
        return {"label": label, "fidelity": 0.0, "success": False,
                "elapsed_s": float(dt), "error": str(e), "n_cav": seq.n_cav}


# ═══════════════════════════════════════════════════════════════════════════
# PHASE 1: N_cav=4 bounded sweep (ultra-fast)
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70, flush=True)
print("  PHASE 1: Bounded Displacement Sweep at N_cav=4 (dim=8)", flush=True)
print("="*70, flush=True)

NC4 = 4
sub4 = make_subspace(NC4)

AMPS = [0.2, 0.3, 0.5, 0.7, 1.0, None]
BLOCKS = [2, 3]

print("\n--- Strategy B: D + SQR + CP ---", flush=True)
for nb in BLOCKS:
    for amp in AMPS:
        seq = make_B(nb, NC4)
        amp_s = f"{amp}" if amp else "none"
        lbl = f"B_nc4_amp{amp_s}_blk{nb}"
        r = run_opt(seq, lbl, sub4, max_amp=amp, multistart=3, maxiter=200)
        r["strategy"] = "D+SQR+CP"; r["n_blocks"] = nb
        results[lbl] = r

print("\n--- Strategy D: D + R + FE ---", flush=True)
for nb in [2, 3]:
    for amp in [0.3, 0.5, 1.0, None]:
        seq = make_D(nb, NC4)
        amp_s = f"{amp}" if amp else "none"
        lbl = f"D_nc4_amp{amp_s}_blk{nb}"
        r = run_opt(seq, lbl, sub4, max_amp=amp, multistart=3, maxiter=200)
        r["strategy"] = "D+R+FE"; r["n_blocks"] = nb
        results[lbl] = r

# Save Phase 1 intermediate
p1_json = {k: {kk: vv for kk, vv in v.items() if kk != "_synth"}
            for k, v in results.items() if isinstance(v, dict)}
(DATA_DIR / "iteration4_phase1_nc4.json").write_text(
    json.dumps(p1_json, indent=2, default=str), encoding="utf-8")
print(f"\n  Phase 1: {len(results)} configs done. Saved.", flush=True)


# ═══════════════════════════════════════════════════════════════════════════
# PHASE 2: Validate Top 3 at N_cav=6 and N_cav=12
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70, flush=True)
print("  PHASE 2: Validate Best at N_cav=6 and N_cav=12", flush=True)
print("="*70, flush=True)

# Rank Phase 1
ranked_p1 = sorted(
    [(k, v["fidelity"]) for k, v in results.items()
     if isinstance(v, dict) and v.get("fidelity", 0) > 0.3],
    key=lambda x: x[1], reverse=True)
print(f"\n  Phase 1 ranking:", flush=True)
for k, f in ranked_p1[:6]:
    print(f"    {k}: F={f:.6f}", flush=True)

# Top 3 to validate
top3 = ranked_p1[:3]
for NC_V, ms, mi in [(6, 3, 250), (12, 2, 200)]:
    sub_v = make_subspace(NC_V)
    print(f"\n  === Validating at N_cav={NC_V} ===", flush=True)
    for k, f_nc4 in top3:
        bv = results[k]
        nb = bv["n_blocks"]
        amp = bv.get("max_amp")
        strat = bv.get("strategy", "")
        seq = make_B(nb, NC_V) if "SQR" in strat else make_D(nb, NC_V)
        lbl = f"{k.replace('nc4', f'nc{NC_V}')}_val"
        r = run_opt(seq, lbl, sub_v, max_amp=amp, multistart=ms, maxiter=mi)
        r["strategy"] = strat; r["n_blocks"] = nb
        r["note"] = f"Validation of {k} (F@nc4={f_nc4:.6f})"
        results[lbl] = r


# ═══════════════════════════════════════════════════════════════════════════
# PHASE 3: Duration Optimization (best at N_cav=6)
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70, flush=True)
print("  PHASE 3: Duration Optimization", flush=True)
print("="*70, flush=True)

overall_best = max(
    (k for k in results if results[k].get("fidelity", 0) > 0.5),
    key=lambda k: results[k]["fidelity"], default=None)
if overall_best:
    bv = results[overall_best]
    nb = bv["n_blocks"]; amp = bv.get("max_amp")
    strat = bv.get("strategy", "")
    nc = bv.get("n_cav", NC4)
    sub = make_subspace(nc)
    seq = make_B(nb, nc) if "SQR" in strat else make_D(nb, nc)
    print(f"  Duration opt: {overall_best} (F={bv['fidelity']:.6f})", flush=True)
    r = run_opt(seq, "duration_opt", sub, max_amp=amp, multistart=3, maxiter=300, dur_wt=0.01)
    r["strategy"] = strat; r["n_blocks"] = nb
    r["note"] = "Duration-optimised version"
    results["duration_opt"] = r


# ═══════════════════════════════════════════════════════════════════════════
# PHASE 4: Wigner Function Comparison
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70, flush=True)
print("  PHASE 4: Wigner Function Comparison", flush=True)
print("="*70, flush=True)

wigner_data = {}
try:
    from cqed_sim.sim.extractors import cavity_wigner, reduced_cavity_state
    import qutip as qt
    HAS_WIGNER = True
except ImportError as e:
    HAS_WIGNER = False
    print(f"  Wigner import failed: {e}", flush=True)

if HAS_WIGNER:
    # Find best result with stored SynthesisResult
    best_with_synth = max(
        (k for k in results
         if isinstance(results[k], dict) and results[k].get("_synth") is not None
         and results[k].get("fidelity", 0) > 0.3),
        key=lambda k: results[k]["fidelity"], default=None)

    if best_with_synth:
        bv = results[best_with_synth]
        sr = bv["_synth"]
        nc_w = bv["n_cav"]
        full_dim = 2 * nc_w
        F_best = bv["fidelity"]

        # Try to get full unitary
        U_full_np = None
        try:
            U_full = sr.simulation.full_operator
            U_full_np = U_full.full() if hasattr(U_full, 'full') else np.asarray(U_full)
        except Exception as e:
            print(f"  Could not extract full unitary: {e}", flush=True)

        if U_full_np is not None:
            print(f"  Wigner source: {best_with_synth} (F={F_best:.6f}, nc={nc_w})", flush=True)

            # Target in full space
            U_tgt_full = np.eye(full_dim, dtype=complex)
            lidx = [0, 1, nc_w, nc_w + 1]
            for i, li in enumerate(lidx):
                for j, lj in enumerate(lidx):
                    U_tgt_full[li, lj] = U_target[i, j]

            basis_labels = ["|g,0>", "|g,1>", "|e,0>", "|e,1>"]
            n_pts = 51; ext = 3.0

            for bi, bl in zip(lidx, basis_labels):
                psi0 = np.zeros(full_dim, dtype=complex)
                psi0[bi] = 1.0
                psi_t = U_tgt_full @ psi0
                psi_a = U_full_np @ psi0

                rho_t = qt.Qobj(np.outer(psi_t, psi_t.conj()),
                                dims=[[2, nc_w], [2, nc_w]])
                rho_a = qt.Qobj(np.outer(psi_a, psi_a.conj()),
                                dims=[[2, nc_w], [2, nc_w]])
                rho_ct = reduced_cavity_state(rho_t)
                rho_ca = reduced_cavity_state(rho_a)

                xv, yv, Wt = cavity_wigner(rho_ct, n_points=n_pts, extent=ext)
                _,  _,  Wa = cavity_wigner(rho_ca, n_points=n_pts, extent=ext)

                try:
                    fcav = float(qt.fidelity(rho_ct, rho_ca)**2)
                except Exception:
                    fcav = float(np.abs(np.trace(rho_ct.full().conj().T @ rho_ca.full())))

                dx = xv[1] - xv[0] if len(xv) > 1 else 1.0
                dy = yv[1] - yv[0] if len(yv) > 1 else 1.0
                l2 = float(np.sqrt(np.sum((Wt - Wa)**2) * dx * dy))
                print(f"    {bl}: F_cav={fcav:.6f}  L2={l2:.6f}", flush=True)

                wigner_data[bl] = {
                    "xvec": xv.tolist(), "yvec": yv.tolist(),
                    "W_target": Wt.tolist(), "W_achieved": Wa.tolist(),
                    "cavity_fidelity": fcav, "l2_distance": l2}

            results["wigner"] = {
                "source": best_with_synth, "fidelity": F_best, "n_cav": nc_w,
                "cavity_fidelities": {bl: wigner_data[bl]["cavity_fidelity"]
                                       for bl in basis_labels},
                "l2_distances": {bl: wigner_data[bl]["l2_distance"]
                                  for bl in basis_labels}}

            # ── Wigner figure ─────────────────────────────────────────────
            fig, axes = plt.subplots(3, 4, figsize=(14, 10))
            vmax = max(max(np.max(np.abs(np.array(wigner_data[bl]["W_target"]))),
                           np.max(np.abs(np.array(wigner_data[bl]["W_achieved"]))))
                       for bl in basis_labels)
            for col, bl in enumerate(basis_labels):
                wd = wigner_data[bl]
                X, Y = np.meshgrid(wd["xvec"], wd["yvec"])
                Wt = np.array(wd["W_target"])
                Wa = np.array(wd["W_achieved"])
                dW = Wt - Wa

                axes[0, col].pcolormesh(X, Y, Wt, cmap="RdBu_r",
                                         vmin=-vmax, vmax=vmax, shading="auto")
                axes[0, col].set_title(f"Target: {bl}", fontsize=9)
                axes[0, col].set_aspect("equal")

                im = axes[1, col].pcolormesh(X, Y, Wa, cmap="RdBu_r",
                                              vmin=-vmax, vmax=vmax, shading="auto")
                axes[1, col].set_title(f"Achieved: $F_{{cav}}$={wd['cavity_fidelity']:.4f}",
                                       fontsize=9)
                axes[1, col].set_aspect("equal")

                dmax = max(np.max(np.abs(dW)), 1e-10)
                axes[2, col].pcolormesh(X, Y, dW, cmap="RdBu_r",
                                         vmin=-dmax, vmax=dmax, shading="auto")
                axes[2, col].set_title(f"$\\Delta W$ (L2={wd['l2_distance']:.4f})",
                                       fontsize=9)
                axes[2, col].set_aspect("equal")
                axes[2, col].set_xlabel("Re($\\alpha$)")
                for row in range(3):
                    if col == 0:
                        axes[row, col].set_ylabel("Im($\\alpha$)")

            fig.suptitle(
                f"Wigner Comparison — {best_with_synth}\n"
                f"$\\mathcal{{F}}$={F_best:.4f}, $N_{{cav}}$={nc_w}",
                fontsize=12, fontweight="bold")
            fig.colorbar(im, ax=axes[:2].ravel().tolist(), label="$W(\\alpha)$",
                         shrink=0.7, pad=0.02)
            fig.tight_layout(rect=[0, 0, 0.93, 0.93])
            for fmt in ("png", "pdf"):
                fig.savefig(FIG_DIR / f"wigner_comparison.{fmt}",
                            dpi=300, bbox_inches="tight")
            plt.close(fig)
            print("  Wigner figure saved.", flush=True)
        else:
            print("  No full unitary available for Wigner.", flush=True)
    else:
        print("  No synthesis result available for Wigner.", flush=True)


# ═══════════════════════════════════════════════════════════════════════════
# PHASE 5: Save ALL Results, Artifacts, Figures
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70, flush=True)
print("  PHASE 5: Save Results & Artifacts", flush=True)
print("="*70, flush=True)

# Final results JSON (without internal objects)
results_json = {}
for k, v in results.items():
    if isinstance(v, dict):
        results_json[k] = {kk: vv for kk, vv in v.items() if kk != "_synth"}
    else:
        results_json[k] = v
(DATA_DIR / "iteration4_results.json").write_text(
    json.dumps(results_json, indent=2, default=str), encoding="utf-8")
print("  iteration4_results.json saved.", flush=True)

# Wigner NPZ
if wigner_data:
    npz = {}
    for bl, wd in wigner_data.items():
        safe = bl.replace("|", "").replace(">", "").replace(",", "_")
        for key in ("xvec", "yvec", "W_target", "W_achieved"):
            npz[f"{safe}_{key}"] = np.array(wd[key])
        npz[f"{safe}_cavity_fidelity"] = np.array([wd["cavity_fidelity"]])
        npz[f"{safe}_l2_distance"] = np.array([wd["l2_distance"]])
    np.savez(str(ARTIFACTS_DIR / "wigner_comparison.npz"), **npz)
    print("  wigner_comparison.npz saved.", flush=True)

# Synthesis artifacts
for strat_tag, prefix in [("D+SQR+CP", "best_strategy_B"), ("D+R+FE", "best_strategy_D")]:
    best_e = max(
        ((k, v) for k, v in results.items()
         if isinstance(v, dict) and strat_tag in v.get("strategy", "")
         and v.get("fidelity", 0) > 0.3),
        key=lambda x: x[1]["fidelity"], default=(None, None))
    if best_e[0]:
        ek, ev = best_e
        sr = ev.get("_synth")
        ap = ARTIFACTS_DIR / f"{prefix}.json"
        if sr is not None:
            try:
                sr.save(str(ap), include_history=True)
                print(f"  {ap.name}: F={ev['fidelity']:.6f}", flush=True)
            except Exception as e:
                # fallback: save metadata
                ap.write_text(json.dumps(
                    {"label": ek, "fidelity": ev["fidelity"],
                     "max_amp": ev.get("max_amp"), "n_cav": ev.get("n_cav"),
                     "n_blocks": ev.get("n_blocks"), "strategy": ev.get("strategy")},
                    indent=2, default=str), encoding="utf-8")
                print(f"  {ap.name} (fallback): {e}", flush=True)


# ── Bounded displacement sweep figure ─────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

for i, nb in enumerate(BLOCKS):
    amps_b, fids_b = [], []
    for amp_val in [0.2, 0.3, 0.5, 0.7, 1.0]:
        k = f"B_nc4_amp{amp_val}_blk{nb}"
        if k in results and results[k].get("fidelity", 0) > 0:
            amps_b.append(amp_val)
            fids_b.append(results[k]["fidelity"])
    if amps_b:
        ax1.plot(amps_b, fids_b, "o-", color=COLORS[i], ms=7, lw=2, label=f"{nb} blocks")
    uk = f"B_nc4_ampnone_blk{nb}"
    if uk in results and results[uk].get("fidelity", 0) > 0:
        ax1.axhline(results[uk]["fidelity"], ls=":", color=COLORS[i], alpha=0.6,
                    label=f"{nb}blk unconstrained")

ax1.set_xlabel("Max displacement $|\\alpha|_{\\max}$")
ax1.set_ylabel("Subspace fidelity")
ax1.set_title("Strategy B: D + SQR + CP ($N_{cav}$=4)")
ax1.axhline(0.99, ls="--", color="red", alpha=0.3, lw=0.5)
ax1.legend(fontsize=8)
ax1.grid(alpha=0.3)

for i, nb in enumerate([2, 3]):
    amps_d, fids_d = [], []
    for amp_val in [0.3, 0.5, 1.0]:
        k = f"D_nc4_amp{amp_val}_blk{nb}"
        if k in results and results[k].get("fidelity", 0) > 0:
            amps_d.append(amp_val)
            fids_d.append(results[k]["fidelity"])
    if amps_d:
        ax2.plot(amps_d, fids_d, "o-", color=COLORS[i+2], ms=7, lw=2, label=f"{nb} blocks")
    uk = f"D_nc4_ampnone_blk{nb}"
    if uk in results and results[uk].get("fidelity", 0) > 0:
        ax2.axhline(results[uk]["fidelity"], ls=":", color=COLORS[i+2], alpha=0.6,
                    label=f"{nb}blk unconstrained")

ax2.set_xlabel("Max displacement $|\\alpha|_{\\max}$")
ax2.set_ylabel("Subspace fidelity")
ax2.set_title("Strategy D: D + R + FE ($N_{cav}$=4)")
ax2.axhline(0.99, ls="--", color="red", alpha=0.3, lw=0.5)
ax2.legend(fontsize=8)
ax2.grid(alpha=0.3)

fig.suptitle("Bounded-Displacement Optimisation: Fidelity vs Max $|\\alpha|$",
             fontsize=13, fontweight="bold")
fig.tight_layout()
for fmt in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"bounded_displacement_sweep.{fmt}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("  bounded_displacement_sweep figure saved.", flush=True)

# ── Truncation convergence figure (best at NC=4,6,12) ─────────────────────
conv_data = {}
for tag in ["B", "D"]:
    for nc_val in [4, 6, 12]:
        keys = [k for k in results if k.startswith(f"{tag}_nc{nc_val}")
                and results[k].get("fidelity", 0) > 0]
        if keys:
            bk = max(keys, key=lambda k: results[k]["fidelity"])
            conv_data.setdefault(tag, {})[nc_val] = {
                "key": bk, "fidelity": results[bk]["fidelity"],
                "max_amp": results[bk].get("max_amp")}

if conv_data:
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for i, (tag, ncdata) in enumerate(sorted(conv_data.items())):
        ncs_sorted = sorted(ncdata.keys())
        fids = [ncdata[nc]["fidelity"] for nc in ncs_sorted]
        ax.plot(ncs_sorted, fids, "o-", color=COLORS[i], ms=8, lw=2,
                label=f"Strategy {tag}")
    ax.set_xlabel("$N_{cav}$ (cavity truncation)")
    ax.set_ylabel("Best subspace fidelity")
    ax.set_title("Truncation Convergence: Best Bounded-Displacement Results")
    ax.axhline(0.99, ls="--", color="red", alpha=0.3, lw=0.5)
    ax.legend(fontsize=9); ax.grid(alpha=0.3)
    fig.tight_layout()
    for fmt in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"truncation_convergence_bounded.{fmt}",
                    dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  truncation_convergence_bounded figure saved.", flush=True)

# ── Ranking bar chart ─────────────────────────────────────────────────────
rank_all = sorted(
    [(k, v) for k, v in results.items()
     if isinstance(v, dict) and "fidelity" in v and "strategy" in v
     and v["fidelity"] > 0],
    key=lambda x: x[1]["fidelity"], reverse=True)
if rank_all:
    top_n = rank_all[:15]
    fig, ax = plt.subplots(figsize=(10, max(4, len(top_n)*0.35)))
    lbs = [s[0] for s in top_n]
    fds = [s[1]["fidelity"] for s in top_n]
    cs = [COLORS[1] if "SQR" in s[1].get("strategy", "") else COLORS[3] for s in top_n]
    ax.barh(range(len(lbs)), fds, color=cs, alpha=0.85, edgecolor="black", lw=0.5)
    ax.set_yticks(range(len(lbs)))
    ax.set_yticklabels(lbs, fontsize=7)
    ax.set_xlabel("Subspace Fidelity")
    ax.set_title("Iteration 4: All Optimisation Results")
    ax.axvline(0.99, ls="--", color="red", alpha=0.4, label="99%")
    ax.axvline(0.999, ls=":", color="darkred", alpha=0.4, label="99.9%")
    ax.legend(fontsize=8); ax.invert_yaxis(); ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()
    for fmt in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"iteration4_ranking.{fmt}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  iteration4_ranking figure saved.", flush=True)


# ═══════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70, flush=True)
print("  ITERATION 4: FINAL SUMMARY", flush=True)
print("="*70, flush=True)

ranked_final = sorted(
    [(k, v["fidelity"]) for k, v in results.items()
     if isinstance(v, dict) and v.get("fidelity", 0) > 0
     and "strategy" in v],
    key=lambda x: x[1], reverse=True)
for i, (k, f) in enumerate(ranked_final[:15]):
    nc = results[k].get("n_cav", "?")
    strat = results[k].get("strategy", "?")
    amp = results[k].get("max_amp", "none")
    nb = results[k].get("n_blocks", "?")
    print(f"  {i+1:2d}. {k}: F={f:.6f} ({strat}, nc={nc}, amp={amp}, blk={nb})", flush=True)

if wigner_data:
    print(f"\nWigner comparison source: {results.get('wigner', {}).get('source', 'N/A')}", flush=True)
    for bl in ["|g,0>", "|g,1>", "|e,0>", "|e,1>"]:
        if bl in wigner_data:
            print(f"  {bl}: F_cav={wigner_data[bl]['cavity_fidelity']:.6f}  L2={wigner_data[bl]['l2_distance']:.6f}", flush=True)

print(f"\nFiles generated:", flush=True)
for d in (DATA_DIR, FIG_DIR, ARTIFACTS_DIR):
    for f in sorted(d.glob("iteration4*")) if "data" in str(d) else sorted(d.glob("*")):
        if f.is_file() and ("iteration4" in f.name or "wigner" in f.name
                            or "bounded" in f.name or "truncation" in f.name
                            or "best_strategy" in f.name or "grape" in f.name):
            print(f"  {f.relative_to(STUDY_ROOT)}", flush=True)

print("\n  DONE.", flush=True)
