"""Iteration 4 v6: Lean bounded-displacement study.

Optimise 6 core configs at N_cav=2 (fast), evaluate at N_cav=4,6,8,12.
Wigner comparison for best result.  Reduced multistart/maxiter for speed.
"""
from __future__ import annotations
import json, sys, time, traceback
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

STUDY_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR   = STUDY_ROOT / "data"
FIG_DIR    = STUDY_ROOT / "figures"
ART_DIR    = STUDY_ROOT / "artifacts"
for d in (DATA_DIR, FIG_DIR, ART_DIR):
    d.mkdir(parents=True, exist_ok=True)

STYLE = STUDY_ROOT.parents[1] / ".github/skills/publication-figures/assets/cqed_style.mplstyle"
if STYLE.exists(): plt.style.use(str(STYLE))

SIM = Path("C:/Users/dazzl/Box/Shyam Shankar Quantum Circuits Group/Users/Users_JianJun/cQED_simulation")
if str(SIM) not in sys.path: sys.path.insert(0, str(SIM))

print("Importing cqed_sim...", flush=True)
t0 = time.perf_counter()
from cqed_sim.unitary_synthesis import (
    Displacement, QubitRotation, SQR, ConditionalPhaseSQR,
    FreeEvolveCondPhase, Subspace, TargetUnitary,
    UnitarySynthesizer, GateSequence, DriftPhaseModel,
    LeakagePenalty, MultiObjective, ExecutionOptions,
    SynthesisConstraints, subspace_unitary_fidelity, simulate_sequence,
)
from cqed_sim.unitary_synthesis.targets import make_target
print(f"Import: {time.perf_counter()-t0:.1f}s", flush=True)

TWO_PI = 2*np.pi
CHI = TWO_PI * (-2.84e6)
COLORS = ['#4477AA','#EE6677','#228833','#CCBB44','#66CCEE','#AA3377','#BBBBBB']

U_target = make_target("cluster", n_match=1)
target = TargetUnitary(U_target, ignore_global_phase=True)
no_drift = DriftPhaseModel(chi=0.0, chi2=0.0, kerr=0.0)
print(f"Target shape: {U_target.shape}", flush=True)

results = {}

def make_subspace(nc):
    return Subspace.custom(2*nc, [0,1,nc,nc+1], ["|g,0>","|g,1>","|e,0>","|e,1>"])

def make_B(nb, nc):
    gates = []
    th = [0.0]*nc; th[0]=np.pi/2
    if nc>1: th[1]=np.pi/4
    ph = [0.0]*nc
    for i in range(nb):
        gates.append(Displacement(f"D{i}", alpha=0.3+0j, duration=200e-9))
        gates.append(SQR(f"S{2*i}", theta_n=th[:], phi_n=ph[:], drift_model=no_drift, duration=400e-9))
        gates.append(ConditionalPhaseSQR(f"CP{i}", phases_n=[0.0]*nc, drift_model=no_drift, duration=200e-9))
        gates.append(SQR(f"S{2*i+1}", theta_n=th[:], phi_n=ph[:], drift_model=no_drift, duration=400e-9))
    gates.append(Displacement(f"D{nb}", alpha=0.3+0j, duration=200e-9))
    return GateSequence(gates=gates, n_cav=nc)

def make_D(nb, nc):
    gates = []
    for i in range(nb):
        gates.append(Displacement(f"D{i}", alpha=0.3+0j, duration=200e-9))
        gates.append(QubitRotation(f"R{i}", theta=np.pi/2, phi=0.0, duration=100e-9))
        gates.append(FreeEvolveCondPhase(
            f"FE{i}", duration=200e-9,
            drift_model=DriftPhaseModel(chi=abs(CHI), chi2=0.0, kerr=0.0),
            optimize_time=True))
    gates.append(Displacement(f"D{nb}", alpha=0.3+0j, duration=200e-9))
    gates.append(QubitRotation(f"R{nb}", theta=np.pi/2, phi=np.pi/2, duration=100e-9))
    return GateSequence(gates=gates, n_cav=nc)

def run_opt(seq, label, subsp, max_amp=None, multistart=3, maxiter=200, dur_wt=0.0):
    print(f"\n  [{label}] amp<={max_amp} ms={multistart} mi={maxiter}", flush=True)
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
        F = subspace_unitary_fidelity(res.simulation.subspace_operator, U_target, gauge="global")
        disp_amps = [float(abs(g.alpha)) for g in res.sequence.gates if hasattr(g,'alpha')]
        print(f"    F={F:.6f} obj={res.objective:.6f} ({dt:.0f}s)", flush=True)
        print(f"    |alpha|: {[round(a,4) for a in disp_amps]}", flush=True)
        return {"label": label, "fidelity": float(F), "objective": float(res.objective),
                "success": bool(res.success), "elapsed_s": dt, "max_amp": max_amp,
                "disp_amps": disp_amps, "_result": res}
    except Exception as e:
        dt = time.perf_counter() - t0
        print(f"    FAILED: {e} ({dt:.0f}s)", flush=True)
        traceback.print_exc()
        return {"label": label, "fidelity": 0.0, "success": False, "elapsed_s": dt, "error": str(e)}

def rescale_sequence(seq, target_nc):
    new_gates = []
    for g in seq.gates:
        if isinstance(g, SQR):
            th = list(g.theta_n)+[0.0]*(target_nc-len(g.theta_n))
            ph = list(g.phi_n)+[0.0]*(target_nc-len(g.phi_n))
            new_gates.append(SQR(g.name, theta_n=th[:target_nc], phi_n=ph[:target_nc],
                                  drift_model=g.drift_model, duration=g.duration))
        elif isinstance(g, ConditionalPhaseSQR):
            ps = list(g.phases_n)+[0.0]*(target_nc-len(g.phases_n))
            new_gates.append(ConditionalPhaseSQR(g.name, phases_n=ps[:target_nc],
                                  drift_model=g.drift_model, duration=g.duration))
        elif isinstance(g, Displacement):
            new_gates.append(Displacement(g.name, alpha=g.alpha, duration=g.duration))
        elif isinstance(g, QubitRotation):
            new_gates.append(QubitRotation(g.name, theta=g.theta, phi=g.phi, duration=g.duration))
        elif isinstance(g, FreeEvolveCondPhase):
            new_gates.append(FreeEvolveCondPhase(g.name, duration=g.duration,
                                  drift_model=g.drift_model, optimize_time=False))
        else:
            new_gates.append(g)
    return GateSequence(gates=new_gates, n_cav=target_nc)

def evaluate_at_ncav(result, target_nc):
    sr = result.get("_result")
    if sr is None: return None
    seq_ext = rescale_sequence(sr.sequence, target_nc)
    sub = make_subspace(target_nc)
    t0 = time.perf_counter()
    try:
        sim_res = simulate_sequence(seq_ext, sub)
        dt = time.perf_counter() - t0
        U_sub = sim_res.subspace_operator
        U_sub_np = U_sub.full() if hasattr(U_sub,'full') else np.asarray(U_sub)
        F = subspace_unitary_fidelity(U_sub_np, U_target, gauge="global")
        U_full = sim_res.full_operator
        U_full_np = U_full.full() if hasattr(U_full,'full') else np.asarray(U_full)
        full_dim = 2*target_nc
        lidx = [0,1,target_nc,target_nc+1]
        leak_sum = 0.0
        for bi in lidx:
            psi0 = np.zeros(full_dim, dtype=complex); psi0[bi]=1.0
            psi_out = U_full_np @ psi0
            sub_pop = sum(abs(psi_out[li])**2 for li in lidx)
            leak_sum += 1.0 - sub_pop
        return {"n_cav": target_nc, "fidelity": float(F),
                "avg_leakage": float(leak_sum/len(lidx)), "eval_time": dt,
                "full_unitary": U_full_np}
    except Exception as e:
        print(f"      eval@nc{target_nc} FAILED: {e}", flush=True)
        return {"n_cav": target_nc, "fidelity": 0.0, "error": str(e)}

# ═══════════════════════════════════════════════════════════════
# PHASE 1: 6 core configs at N_cav=2
# ═══════════════════════════════════════════════════════════════
print("\n"+"="*60+"\n  PHASE 1: Synthesis at N_cav=2\n"+"="*60, flush=True)
NC_OPT = 2
sub_opt = make_subspace(NC_OPT)

CONFIGS = [
    ("B", 2, 0.3, "B2_amp0.3"),
    ("B", 2, 0.5, "B2_amp0.5"),
    ("B", 2, None, "B2_unbounded"),
    ("D", 3, 0.3, "D3_amp0.3"),
    ("D", 3, 0.5, "D3_amp0.5"),
    ("D", 3, None, "D3_unbounded"),
]

for strat, nb, amp, lbl in CONFIGS:
    seq = make_B(nb, NC_OPT) if strat=="B" else make_D(nb, NC_OPT)
    r = run_opt(seq, lbl, sub_opt, max_amp=amp, multistart=3, maxiter=200)
    r["strategy"] = "D+SQR+CP" if strat=="B" else "D+R+FE"
    r["n_blocks"] = nb; r["n_cav_opt"] = NC_OPT
    results[lbl] = r
    # Save incremental results
    inc = {k: {kk:vv for kk,vv in v.items() if kk!="_result"} for k,v in results.items()}
    (DATA_DIR/"iteration4_incremental.json").write_text(json.dumps(inc, indent=2, default=str))

# ═══════════════════════════════════════════════════════════════
# PHASE 2: Evaluate at N_cav=4,6,8,12
# ═══════════════════════════════════════════════════════════════
print("\n"+"="*60+"\n  PHASE 2: Cross-N_cav Evaluation\n"+"="*60, flush=True)
EVAL_NCS = [4, 6, 8, 12]
cross_nc = {}

for lbl, r in results.items():
    if r.get("fidelity",0) < 0.3 or r.get("_result") is None:
        continue
    cross_nc[lbl] = {"nc2": r["fidelity"]}
    for nc in EVAL_NCS:
        ev = evaluate_at_ncav(r, nc)
        if ev:
            cross_nc[lbl][f"nc{nc}"] = ev.get("fidelity",0.0)
            cross_nc[lbl][f"leak_nc{nc}"] = ev.get("avg_leakage",0.0)
            if nc == 12: cross_nc[lbl]["U_full_nc12"] = ev.get("full_unitary")
            print(f"    [{lbl}@nc{nc}] F={ev.get('fidelity',0):.6f} leak={ev.get('avg_leakage',0):.6f}", flush=True)

print("\n  Cross-N_cav Summary:", flush=True)
hdr = f"  {'Config':<20s}"+"".join(f" {'nc'+str(n):>8s}" for n in [2]+EVAL_NCS)+" {'leak12':>8s}"
print(hdr, flush=True)
for lbl in sorted(cross_nc.keys()):
    d = cross_nc[lbl]
    row = f"  {lbl:<20s}"
    for k in ["nc2"]+[f"nc{n}" for n in EVAL_NCS]:
        v = d.get(k)
        row += f" {v:8.4f}" if v is not None else "      ---"
    row += f" {d.get('leak_nc12',0):8.6f}"
    print(row, flush=True)

# ═══════════════════════════════════════════════════════════════
# PHASE 3: Duration-optimised variant of best bounded
# ═══════════════════════════════════════════════════════════════
print("\n"+"="*60+"\n  PHASE 3: Duration Optimisation\n"+"="*60, flush=True)
best_key = max(
    ((l,d) for l,d in cross_nc.items() if results[l].get("max_amp") is not None and d.get("nc12",0)>0.3),
    key=lambda x: x[1].get("nc12",0), default=(None,None))

if best_key[0]:
    bk = best_key[0]; bv = results[bk]
    nb, amp, strat = bv["n_blocks"], bv.get("max_amp"), bv["strategy"]
    print(f"  Best bounded@nc12: {bk} (F@nc12={best_key[1].get('nc12',0):.6f})", flush=True)
    seq = make_B(nb, NC_OPT) if "SQR" in strat else make_D(nb, NC_OPT)
    r = run_opt(seq, "dur_opt", sub_opt, max_amp=amp, multistart=3, maxiter=200, dur_wt=0.01)
    r["strategy"] = strat; r["n_blocks"] = nb; r["n_cav_opt"] = NC_OPT
    results["dur_opt"] = r
    ev12 = evaluate_at_ncav(r, 12)
    if ev12: print(f"    dur_opt@nc12: F={ev12.get('fidelity',0):.6f}", flush=True)
else:
    print("  No bounded config with F>0.3 at nc12.", flush=True)

# ═══════════════════════════════════════════════════════════════
# PHASE 4: Wigner Comparison
# ═══════════════════════════════════════════════════════════════
print("\n"+"="*60+"\n  PHASE 4: Wigner Comparison\n"+"="*60, flush=True)
wigner_data = {}
try:
    from cqed_sim.sim.extractors import cavity_wigner, reduced_cavity_state
    import qutip as qt
    HAS_WIG = True; print("  Wigner imports OK.", flush=True)
except ImportError as e:
    HAS_WIG = False; print(f"  Wigner import failed: {e}", flush=True)

if HAS_WIG:
    best_wig_key, best_wig_F, best_wig_U = None, 0, None
    nc_wig = 12
    for lbl, d in cross_nc.items():
        U = d.get("U_full_nc12")
        f12 = d.get("nc12",0)
        if U is not None and f12 > best_wig_F:
            best_wig_key, best_wig_F, best_wig_U = lbl, f12, U

    if best_wig_key and best_wig_U is not None:
        print(f"  Source: {best_wig_key} (F@nc12={best_wig_F:.6f})", flush=True)
        full_dim = 2*nc_wig
        U_tgt_full = np.eye(full_dim, dtype=complex)
        lidx = [0,1,nc_wig,nc_wig+1]
        for i,li in enumerate(lidx):
            for j,lj in enumerate(lidx):
                U_tgt_full[li,lj] = U_target[i,j]

        basis_labels = ["|g,0>","|g,1>","|e,0>","|e,1>"]
        n_pts, ext = 51, 3.0

        for bi, bl in zip(lidx, basis_labels):
            psi0 = np.zeros(full_dim, dtype=complex); psi0[bi]=1.0
            psi_t = U_tgt_full @ psi0
            psi_a = best_wig_U @ psi0
            rho_t = qt.Qobj(np.outer(psi_t, psi_t.conj()), dims=[[2,nc_wig],[2,nc_wig]])
            rho_a = qt.Qobj(np.outer(psi_a, psi_a.conj()), dims=[[2,nc_wig],[2,nc_wig]])
            rho_ct = reduced_cavity_state(rho_t)
            rho_ca = reduced_cavity_state(rho_a)
            xv, yv, Wt = cavity_wigner(rho_ct, n_points=n_pts, extent=ext)
            _,  _,  Wa = cavity_wigner(rho_ca, n_points=n_pts, extent=ext)
            try:
                fcav = float(qt.fidelity(rho_ct, rho_ca)**2)
            except Exception:
                ov = rho_ct.full().conj().T @ rho_ca.full()
                fcav = float(np.abs(np.trace(ov)))
            dx = xv[1]-xv[0] if len(xv)>1 else 1.0
            dy = yv[1]-yv[0] if len(yv)>1 else 1.0
            l2 = float(np.sqrt(np.sum((Wt-Wa)**2)*dx*dy))
            print(f"    {bl}: F_cav={fcav:.6f}  L2={l2:.6f}", flush=True)
            wigner_data[bl] = {"xvec":xv.tolist(),"yvec":yv.tolist(),
                               "W_target":Wt.tolist(),"W_achieved":Wa.tolist(),
                               "cavity_fidelity":fcav,"l2_distance":l2}

        # Wigner figure
        fig, axes = plt.subplots(3, 4, figsize=(14, 10))
        vmax = max(max(np.max(np.abs(np.array(wigner_data[bl]["W_target"]))),
                       np.max(np.abs(np.array(wigner_data[bl]["W_achieved"]))))
                   for bl in basis_labels if bl in wigner_data)
        for col, bl in enumerate(basis_labels):
            if bl not in wigner_data: continue
            wd = wigner_data[bl]
            X,Y = np.meshgrid(wd["xvec"], wd["yvec"])
            Wt_ = np.array(wd["W_target"]); Wa_ = np.array(wd["W_achieved"])
            dW = Wt_ - Wa_
            axes[0,col].pcolormesh(X,Y,Wt_,cmap="RdBu_r",vmin=-vmax,vmax=vmax,shading="auto")
            axes[0,col].set_title(f"Target: {bl}",fontsize=9); axes[0,col].set_aspect("equal")
            im = axes[1,col].pcolormesh(X,Y,Wa_,cmap="RdBu_r",vmin=-vmax,vmax=vmax,shading="auto")
            axes[1,col].set_title(f"Achieved ($F_c$={wd['cavity_fidelity']:.4f})",fontsize=9)
            axes[1,col].set_aspect("equal")
            dmax = max(np.max(np.abs(dW)),1e-10)
            axes[2,col].pcolormesh(X,Y,dW,cmap="RdBu_r",vmin=-dmax,vmax=dmax,shading="auto")
            axes[2,col].set_title(f"$\\Delta W$ (L2={wd['l2_distance']:.4f})",fontsize=9)
            axes[2,col].set_aspect("equal")
            axes[2,col].set_xlabel("Re($\\alpha$)")
            for row in range(3):
                if col==0: axes[row,col].set_ylabel("Im($\\alpha$)")
        fig.suptitle(f"Wigner — {best_wig_key}  $\\mathcal{{F}}_{{sub}}$={best_wig_F:.4f}  ($N_{{cav}}$={nc_wig})",
                     fontsize=12, fontweight="bold")
        fig.colorbar(im, ax=axes[:2].ravel().tolist(), label="$W(\\alpha)$", shrink=0.7, pad=0.02)
        fig.tight_layout(rect=[0,0,0.93,0.93])
        for fmt in ("png","pdf"):
            fig.savefig(FIG_DIR/f"wigner_comparison.{fmt}", dpi=300, bbox_inches="tight")
        plt.close(fig)
        print("  Wigner figure saved.", flush=True)
        npz_d = {}
        for bl, wd in wigner_data.items():
            s = bl.replace("|","").replace(">","").replace(",","_")
            for key in ("xvec","yvec","W_target","W_achieved"):
                npz_d[f"{s}_{key}"] = np.array(wd[key])
            npz_d[f"{s}_cavity_fidelity"] = np.array([wd["cavity_fidelity"]])
        np.savez(str(ART_DIR/"wigner_comparison.npz"), **npz_d)
        print("  Wigner NPZ saved.", flush=True)
    else:
        print("  No unitary for Wigner.", flush=True)

# ═══════════════════════════════════════════════════════════════
# PHASE 5: Summary Figures & Artifacts
# ═══════════════════════════════════════════════════════════════
print("\n"+"="*60+"\n  PHASE 5: Figures & Artifacts\n"+"="*60, flush=True)

# Figure 1: Bounded vs Unbounded at nc2 and nc12
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
for tag, col, mk in [("B2",COLORS[0],"o"),("D3",COLORS[2],"^")]:
    amps_f2, amps_f12, amps_lk = [],[],[]
    for lbl in sorted(cross_nc.keys()):
        if not lbl.startswith(tag+"_amp"): continue
        a = results[lbl].get("max_amp")
        if a is None: continue
        amps_f2.append((a, results[lbl]["fidelity"]))
        amps_f12.append((a, cross_nc[lbl].get("nc12",0)))
        amps_lk.append((a, cross_nc[lbl].get("leak_nc12",0)))
    if amps_f2:
        amps_f2.sort(); amps_f12.sort(); amps_lk.sort()
        axes[0].plot([x[0] for x in amps_f2],[x[1] for x in amps_f2],f"{mk}-",color=col,ms=7,lw=2,label=tag)
        axes[1].plot([x[0] for x in amps_f12],[x[1] for x in amps_f12],f"{mk}-",color=col,ms=7,lw=2,label=tag)
        axes[2].plot([x[0] for x in amps_lk],[x[1] for x in amps_lk],f"{mk}-",color=col,ms=7,lw=2,label=tag)
    ub = f"{tag}_unbounded"
    if ub in results and results[ub].get("fidelity",0)>0:
        axes[0].axhline(results[ub]["fidelity"],ls=":",color=col,alpha=0.5)
    if ub in cross_nc:
        f12 = cross_nc[ub].get("nc12",0)
        if f12>0: axes[1].axhline(f12,ls=":",color=col,alpha=0.5)
        l12 = cross_nc[ub].get("leak_nc12",0)
        if l12>0: axes[2].axhline(l12,ls=":",color=col,alpha=0.5)

axes[0].set(xlabel="Max $|\\alpha|$",ylabel="$\\mathcal{F}$ (N$_{cav}$=2)",title="(a) Ideal-mode fidelity")
axes[1].set(xlabel="Max $|\\alpha|$",ylabel="$\\mathcal{F}$ (N$_{cav}$=12)",title="(b) Physical fidelity (N$_{cav}$=12)")
axes[2].set(xlabel="Max $|\\alpha|$",ylabel="Avg leakage",title="(c) Leakage at N$_{cav}$=12")
for ax in axes: ax.legend(fontsize=8); ax.grid(alpha=0.3)
axes[2].set_yscale("log")
fig.suptitle("Bounded-Displacement: Ideal-Mode Optimisation $\\rightarrow$ N$_{cav}$=12 Evaluation",fontsize=13,fontweight="bold")
fig.tight_layout()
for fmt in ("png","pdf"):
    fig.savefig(FIG_DIR/f"bounded_displacement_sweep.{fmt}",dpi=300,bbox_inches="tight")
plt.close(fig)
print("  bounded_displacement_sweep saved.", flush=True)

# Figure 2: Truncation convergence
fig, ax = plt.subplots(figsize=(8, 5))
for i, lbl in enumerate(sorted(cross_nc.keys())):
    d = cross_nc[lbl]; r = results[lbl]
    amp = r.get("max_amp")
    ncs = [2]+EVAL_NCS
    fids = [d.get(f"nc{n}",0) for n in ncs]
    if all(f>0 for f in fids):
        label_str = f"{lbl}" + (f" (|α|≤{amp})" if amp else " (unbound)")
        ls = "-" if amp else "--"
        ax.plot(ncs, fids, f"o{ls}", ms=6, lw=1.5, label=label_str, color=COLORS[i%len(COLORS)])
ax.set(xlabel="$N_{cav}$",ylabel="Subspace fidelity",title="Truncation Convergence: Bounded vs Unbounded")
ax.axhline(0.99,ls="--",color="gray",alpha=0.3); ax.legend(fontsize=7,ncol=2); ax.grid(alpha=0.3)
fig.tight_layout()
for fmt in ("png","pdf"):
    fig.savefig(FIG_DIR/f"truncation_convergence_bounded.{fmt}",dpi=300,bbox_inches="tight")
plt.close(fig)
print("  truncation_convergence_bounded saved.", flush=True)

# Save JSON results
output = {"configs":{}, "cross_ncav":{}, "wigner_summary":{}}
for lbl, r in results.items():
    output["configs"][lbl] = {k:v for k,v in r.items() if k!="_result"}
for lbl, d in cross_nc.items():
    output["cross_ncav"][lbl] = {k:v for k,v in d.items() if not isinstance(v,np.ndarray)}
for bl, wd in wigner_data.items():
    output["wigner_summary"][bl] = {"cavity_fidelity": wd["cavity_fidelity"], "l2_distance": wd["l2_distance"]}
(DATA_DIR/"iteration4_results.json").write_text(json.dumps(output, indent=2, default=str))
print("  iteration4_results.json saved.", flush=True)

# Save best artifacts
for tag, prefix in [("B","best_strategy_B"),("D","best_strategy_D")]:
    best12 = max(
        ((l, d.get("nc12",0)) for l,d in cross_nc.items()
         if l.startswith(tag) and results[l].get("max_amp") is not None),
        key=lambda x: x[1], default=(None,0))
    if best12[0] and best12[1]>0.1:
        sr = results[best12[0]].get("_result")
        ap = ART_DIR/f"{prefix}.json"
        if sr:
            try:
                sr.save(str(ap), include_history=True)
                print(f"  {ap.name}: {best12[0]} F@nc12={best12[1]:.6f}", flush=True)
            except Exception as e:
                ap.write_text(json.dumps({
                    "label":best12[0],"fidelity_nc2":results[best12[0]]["fidelity"],
                    "fidelity_nc12":best12[1],"max_amp":results[best12[0]].get("max_amp"),
                },indent=2))
                print(f"  {ap.name} (fallback): {e}", flush=True)

# Final summary
print("\n"+"="*60+"\n  FINAL SUMMARY\n"+"="*60, flush=True)
print(f"  {'Config':<20s} {'nc2':>8s} {'nc4':>8s} {'nc6':>8s} {'nc8':>8s} {'nc12':>8s} {'leak12':>8s}", flush=True)
for lbl in sorted(cross_nc.keys()):
    d = cross_nc[lbl]
    row = f"  {lbl:<20s}"
    for k in ["nc2"]+[f"nc{n}" for n in EVAL_NCS]:
        v = d.get(k); row += f" {v:8.4f}" if v is not None else "      ---"
    row += f" {d.get('leak_nc12',0):8.6f}"
    print(row, flush=True)

if wigner_data:
    print(f"\n  Wigner (source={best_wig_key}):", flush=True)
    for bl in basis_labels:
        if bl in wigner_data:
            print(f"    {bl}: F_cav={wigner_data[bl]['cavity_fidelity']:.6f}", flush=True)

print("\n  Output files:", flush=True)
for d in (DATA_DIR, FIG_DIR, ART_DIR):
    for f in sorted(d.iterdir()):
        if f.is_file() and any(x in f.name for x in ("iteration4","wigner","bounded","truncation","best_strategy")):
            print(f"    {f.relative_to(STUDY_ROOT)}", flush=True)

print("\n  DONE.", flush=True)
