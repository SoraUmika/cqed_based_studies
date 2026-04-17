"""Cluster-state holographic simulation — corrected execution script.

Phases:
  1. Verify cqed_sim cluster target vs canonical definition.
  2. Compute ideal cluster-state observables.
  3. Decompose per-site unitary into QubitRotation + Displacement + SNAP.
  4. GRAPE-optimise each SNAP gate for minimum pulse duration.
  5. Timing analysis.
  6. GRAPE full-unitary reference.
"""
from __future__ import annotations

import json, os, sys, time, traceback
from pathlib import Path
from typing import Any

os.environ["PYTHONIOENCODING"] = "utf-8"
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import qutip as qt

# ── paths ────────────────────────────────────────────────────────────────────
STUDY   = Path(__file__).resolve().parents[1]
WS      = STUDY.parents[1]
DATA    = STUDY / "data";   DATA.mkdir(parents=True, exist_ok=True)
FIG     = STUDY / "figures"; FIG.mkdir(parents=True, exist_ok=True)
SIM     = Path("C:/Users/dazzl/Box/Shyam Shankar Quantum Circuits Group"
               "/Users/Users_JianJun/cQED_simulation")

if str(SIM) not in sys.path:
    sys.path.insert(0, str(SIM))
STYLE = WS / ".github/skills/publication-figures/assets/cqed_style.mplstyle"
if STYLE.exists():
    plt.style.use(str(STYLE))

TWO_PI = 2.0 * np.pi
COLORS = ['#4477AA','#EE6677','#228833','#CCBB44','#66CCEE','#AA3377','#BBBBBB']

# ── cqed_sim ─────────────────────────────────────────────────────────────────
from cqed_sim import (
    DispersiveTransmonCavityModel, FrameSpec,
    ModelControlChannelSpec, PiecewiseConstantTimeGrid,
)
from cqed_sim.gates.bosonic import snap as snap_gate
from cqed_sim.optimal_control import (
    GrapeConfig, GrapeSolver,
    LeakagePenalty   as OCLeak,
    UnitaryObjective as OCObj,
    build_control_problem_from_model,
)
from cqed_sim.unitary_synthesis import (
    Displacement  as SDisp,
    QubitRotation as SRot,
    SNAP          as SSNAP,
    Subspace, TargetUnitary, UnitarySynthesizer,
    subspace_unitary_fidelity,
)
from cqed_sim.unitary_synthesis.targets import make_target

# ── device parameters (AGENTS.md) ───────────────────────────────────────────
OQ = TWO_PI*6.150e9;  OC = TWO_PI*5.241e9
AL = TWO_PI*(-255e6);  CH = TWO_PI*(-2.84e6)
CP = TWO_PI*(-21e3);   KR = TWO_PI*(-28e3)
NCAV = 8;  NTR = 2;  NSITES = 6

SNAP_THR  = 0.999
AMP_BND   = TWO_PI*50e6
GR_ITER   = 200
GR_SEEDS  = [17, 42]
ROT_NS    = 16.0
DISP_NS   = 48.0

# ── helpers ──────────────────────────────────────────────────────────────────
def _jdump(path, obj):
    def _d(o):
        if isinstance(o, np.ndarray): return o.tolist()
        if isinstance(o, np.floating): return float(o)
        if isinstance(o, np.integer):  return int(o)
        if isinstance(o, (complex, np.complexfloating)):
            return {"re": float(np.real(o)), "im": float(np.imag(o))}
        return str(o)
    path.write_text(json.dumps(obj, indent=2, default=_d), encoding="utf-8")

def _model(nc=NCAV, nt=NTR):
    return DispersiveTransmonCavityModel(
        omega_c=OC, omega_q=OQ, alpha=AL,
        chi=CH, chi_higher=(CP,), kerr=KR, n_cav=nc, n_tr=nt)

def _frame(m):
    return FrameSpec(omega_c_frame=m.omega_c, omega_q_frame=m.omega_q)

def _sub(nc=NCAV):
    return Subspace.custom(NTR*nc, (0,1,nc,nc+1),
                           ("|g,0>","|g,1>","|e,0>","|e,1>"))

def _sfid(a, b): return float(abs(np.vdot(a,b))**2)
def _ufid(U, V): d=U.shape[0]; return float(abs(np.trace(U.conj().T@V))**2/d**2)

# ═════════════════════════════════════════════════════════════════════════════
# PHASE 1 — Target verification
# ═════════════════════════════════════════════════════════════════════════════
def phase1():
    print("\n"+"="*70+"\n  PHASE 1: Target Verification\n"+"="*70)
    R = {}
    U = make_target("cluster", n_match=1)
    R["shape"] = list(U.shape)
    R["unitarity_err"] = float(np.linalg.norm(U.conj().T@U - np.eye(4)))
    print(f"  ||U†U-I|| = {R['unitarity_err']:.2e}")

    H  = np.array([[1,1],[1,-1]], dtype=complex)/np.sqrt(2)
    CZ = np.diag([1.,1.,1.,-1.]).astype(complex)
    SW = np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]], dtype=complex)
    Uns = CZ @ np.kron(H, np.eye(2))
    Usw = SW @ CZ @ np.kron(H, np.eye(2))

    R["match_swap"]   = bool(np.allclose(U, Usw, atol=1e-12))
    R["match_noswap"] = bool(np.allclose(U, Uns, atol=1e-12))
    R["fid_swap"]     = _ufid(U, Usw)
    R["fid_noswap"]   = _ufid(U, Uns)
    print(f"  Matches SWAP·CZ·(H⊗I): {R['match_swap']}")
    print(f"  Matches CZ·(H⊗I):      {R['match_noswap']}")

    # MPS tensors
    def _mps_ns(M):
        A={}
        for s in range(2):
            A[s]=np.zeros((2,2),dtype=complex)
            for a in range(2):
                for b in range(2): A[s][a,b]=M[s*2+b,a]
        return A
    def _mps_sw(M):
        B={}
        for s in range(2):
            B[s]=np.zeros((2,2),dtype=complex)
            for a in range(2):
                for b in range(2): B[s][a,b]=M[b*2+s,a]
        return B

    A0=np.eye(2)/np.sqrt(2); A1=np.diag([1,-1]).astype(complex)/np.sqrt(2)
    Ans = _mps_ns(Uns); Bsw = _mps_sw(U)
    R["mps_ns_ok"] = bool(np.allclose(Ans[0],A0) and np.allclose(Ans[1],A1))
    R["mps_sw_ok"] = bool(np.allclose(Bsw[0],A0) and np.allclose(Bsw[1],A1))
    print(f"  MPS(noswap) A^0=I/√2,A^1=Z/√2: {R['mps_ns_ok']}")
    print(f"  MPS(swap)   B^0=I/√2,B^1=Z/√2: {R['mps_sw_ok']}")

    # Holographic construction
    def _holo(tens, n):
        dim=2**n; psi=np.zeros(dim, dtype=complex)
        for idx in range(dim):
            bits=[(idx>>(n-1-k))&1 for k in range(n)]
            v=np.array([1,0],dtype=complex)
            for s in bits: v=tens[s].T@v
            psi[idx]=v.sum()
        nm=np.linalg.norm(psi)
        return psi/nm if nm>1e-15 else psi

    def _canon(n):
        dim=2**n; psi=np.zeros(dim,dtype=complex); psi[0]=1.
        Hf=np.array([1.],dtype=complex)
        for _ in range(n): Hf=np.kron(Hf,H)
        psi=Hf@psi
        for i in range(n-1):
            for idx in range(dim):
                bits=[(idx>>(n-1-k))&1 for k in range(n)]
                if bits[i]==1 and bits[i+1]==1: psi[idx]*=-1
        return psi

    psi_ns = _holo(_mps_ns(Uns), NSITES)
    psi_sw = _holo(_mps_sw(U),   NSITES)
    psi_cn = _canon(NSITES)
    R["fid_holo_ns"] = _sfid(psi_ns, psi_cn)
    R["fid_holo_sw"] = _sfid(psi_sw, psi_cn)
    R["nsites"]      = NSITES
    ok = R["fid_holo_ns"]>0.999 or R["fid_holo_sw"]>0.999
    R["passed"] = bool(ok)
    print(f"  F(holo-ns,canon)={R['fid_holo_ns']:.12f}")
    print(f"  F(holo-sw,canon)={R['fid_holo_sw']:.12f}")
    print(f"  PHASE 1: {'PASS' if ok else 'FAIL'}")
    return R, U, psi_cn

# ═════════════════════════════════════════════════════════════════════════════
# PHASE 2 — Ideal observables
# ═════════════════════════════════════════════════════════════════════════════
def phase2(psi, n):
    print("\n"+"="*70+"\n  PHASE 2: Ideal Observables\n"+"="*70)
    R = {"n": n}
    I2=np.eye(2,dtype=complex)
    X=np.array([[0,1],[1,0]],dtype=complex)
    Y=np.array([[0,-1j],[1j,0]],dtype=complex)
    Z=np.array([[1,0],[0,-1]],dtype=complex)

    def _nq(ops):
        m=np.array([1.],dtype=complex)
        for o in ops: m=np.kron(m,o)
        return m
    def _ev(op): return float(np.real(psi.conj()@op@psi))

    # 2a single-site Pauli
    print("  2a. Single-site Pauli:")
    pauli=[]
    for i in range(n):
        ox=[I2]*n; ox[i]=X; oy=[I2]*n; oy[i]=Y; oz=[I2]*n; oz[i]=Z
        ex,ey,ez = _ev(_nq(ox)),_ev(_nq(oy)),_ev(_nq(oz))
        pauli.append({"site":i,"X":ex,"Y":ey,"Z":ez})
        print(f"    site {i}: <X>={ex:+.8f} <Y>={ey:+.8f} <Z>={ez:+.8f}")
    R["pauli"]=pauli
    R["max_pauli"]=max(max(abs(r["X"]),abs(r["Y"]),abs(r["Z"])) for r in pauli)

    # 2b stabilisers
    print("  2b. Stabilisers:")
    stab=[]
    for i in range(n):
        ops=[I2]*n; ops[i]=X
        if i>0: ops[i-1]=Z
        if i<n-1: ops[i+1]=Z
        ki=_ev(_nq(ops))
        stab.append({"site":i,"K":ki})
        print(f"    K_{i} = {ki:+.10f}")
    R["stab"]=stab
    R["max_stab_dev"]=max(abs(s["K"]-1.) for s in stab)

    # 2c bulk ZXZ
    print("  2c. ZXZ:")
    zxz=[]
    for i in range(1,n-1):
        ops=[I2]*n; ops[i-1]=Z; ops[i]=X; ops[i+1]=Z
        v=_ev(_nq(ops)); zxz.append({"i":i,"v":v})
        print(f"    Z_{i-1}X_{i}Z_{i+1} = {v:+.10f}")
    R["zxz"]=zxz

    # 2d string-order
    print("  2d. String-order:")
    strs=[]
    for i in range(n):
        for j in range(i+2,n):
            ops=[I2]*n; ops[i]=Z; ops[j]=Z
            for k in range(i+1,j): ops[k]=X
            v=_ev(_nq(ops)); strs.append({"i":i,"j":j,"v":v})
            mid=" ".join(f"X_{k}" for k in range(i+1,j))
            print(f"    Z_{i} {mid} Z_{j} = {v:+.10f}")
    R["strings"]=strs
    print("  PHASE 2 COMPLETE")
    return R

# ═════════════════════════════════════════════════════════════════════════════
# PHASE 3 — Decomposition
# ═════════════════════════════════════════════════════════════════════════════
def phase3(Ut):
    print("\n"+"="*70+"\n  PHASE 3: Decomposition\n"+"="*70)
    sub = _sub()
    Ds = DISP_NS*1e-9; Rs = ROT_NS*1e-9; Ss = 200e-9

    best=None; bfid=0.; bns=0
    for ns in [2,3,4]:
        print(f"  {ns}-SNAP …")
        gl = []
        for k in range(ns):
            gl.append(SDisp(name=f"D{k+1}",alpha=.1+0j,duration=Ds,optimize_time=False))
            gl.append(SRot(name=f"R{k+1}",theta=.1,phi=0.,duration=Rs,optimize_time=False))
            gl.append(SSNAP(name=f"SNAP{k+1}",phases=[0.]*NCAV,duration=Ss,optimize_time=False))
        gl.append(SDisp(name=f"D{ns+1}",alpha=.1+0j,duration=Ds,optimize_time=False))
        gl.append(SRot(name=f"R{ns+1}",theta=.1,phi=0.,duration=Rs,optimize_time=False))

        tgt = TargetUnitary(Ut, ignore_global_phase=True)
        try:
            syn = UnitarySynthesizer(
                primitives=gl, subspace=sub, target=tgt,
                seed=42, optimize_times=False, optimizer="powell")
            sr = syn.fit(multistart=5, maxiter=500)
            fid = float(sr.report.get("fidelity", 1.0 - sr.objective))
            print(f"    F={fid:.8f}  obj={sr.objective:.6e}")
            if fid > bfid: bfid=fid; best=sr; bns=ns
            if fid > 0.9999:
                print(f"    ✓ sufficient"); break
        except Exception: traceback.print_exc()

    if best is None: return {"error":"no decomposition"}

    seq = best.sequence
    gp=[]; snaps=[]; nr=nd=0
    for g in seq.gates:
        gt=type(g).__name__; e={"name":g.name,"type":gt}
        if gt=="SNAP":
            ph=list(g.phases) if hasattr(g,"phases") else []
            e["phases"]=ph; snaps.append({"name":g.name,"phases":ph})
        elif gt=="QubitRotation":
            e["theta"]=float(g.theta); e["phi"]=float(g.phi); nr+=1
        elif gt=="Displacement":
            e["alpha_re"]=float(np.real(g.alpha))
            e["alpha_im"]=float(np.imag(g.alpha)); nd+=1
        gp.append(e)

    sub_U = best.simulation.subspace_operator if best.simulation else None

    R = {"n_snap":bns, "fid":bfid, "obj":float(best.objective),
         "gates":gp, "snaps":snaps, "nr":nr, "nd":nd}

    print(f"  Best: {bns} SNAP, F={bfid:.8f}")
    for e in gp:
        if e["type"]=="SNAP":
            print(f"    {e['name']}: {[f'{p:.3f}' for p in e['phases'][:4]]}")
        elif e["type"]=="QubitRotation":
            print(f"    {e['name']}: θ={e['theta']:.4f} φ={e['phi']:.4f}")
        elif e["type"]=="Displacement":
            print(f"    {e['name']}: α={e['alpha_re']:.4f}+{e['alpha_im']:.4f}j")
    return R

# ═════════════════════════════════════════════════════════════════════════════
# PHASE 4 — SNAP GRAPE
# ═════════════════════════════════════════════════════════════════════════════
def phase4(snaps):
    print("\n"+"="*70+"\n  PHASE 4: SNAP GRAPE Optimisation\n"+"="*70)
    mod=_model(); fr=_frame(mod); sub=_sub()
    R={}
    for si in snaps:
        sn=si["name"]; ph=si["phases"]
        print(f"\n  {sn} …")
        snap_q = qt.tensor(qt.qeye(NTR), snap_gate(ph[:NCAV], dim=NCAV))
        snap_s = np.asarray(sub.restrict_operator(snap_q.full()), dtype=complex)

        durs=[400,300,200,150,100]; sweep=[]
        for dns in durs:
            ns=max(10,dns//4); dt=dns*1e-9/ns
            bf=0.; br={}
            for sd in GR_SEEDS:
                try:
                    p=build_control_problem_from_model(
                        mod, frame=fr,
                        time_grid=PiecewiseConstantTimeGrid.uniform(steps=ns,dt_s=dt),
                        channel_specs=(ModelControlChannelSpec(
                            name="qubit",target="qubit",
                            quadratures=("I","Q"),
                            amplitude_bounds=(-AMP_BND,AMP_BND)),),
                        objectives=(OCObj(
                            target_operator=snap_s, subspace=sub,
                            ignore_global_phase=True, name=f"{sn}_{dns}"),),
                        penalties=(OCLeak(subspace=sub, weight=0.05),))
                    g=GrapeSolver(GrapeConfig(maxiter=GR_ITER,seed=sd,random_scale=.3)).solve(p)
                    f=float(g.metrics.get("nominal_fidelity",g.metrics.get("fidelity",0.)))
                    if f>bf: bf=f; br={"seed":sd,"fid":f}
                except Exception as e:
                    print(f"      sd={sd}@{dns}ns: {e}")
            ok=bf>=SNAP_THR
            sweep.append({"dns":dns,"fid":bf,"ok":ok}); sweep[-1].update(br)
            print(f"    {dns:4d}ns: F={bf:.6f} {'✓' if ok else '✗'}")
            if not ok and len(sweep)>=2 and sweep[-2]["ok"]:
                print(f"    → min≈{sweep[-2]['dns']}ns"); break

        passing=[s for s in sweep if s["ok"]]
        md=min((s["dns"] for s in passing),default=sweep[0]["dns"])
        R[sn]={"phases":ph,"sweep":sweep,"min_ns":md,
               "best_fid":max((s["fid"] for s in sweep),default=0.)}
    R["threshold"]=SNAP_THR
    print("  PHASE 4 COMPLETE"); return R

# ═════════════════════════════════════════════════════════════════════════════
# PHASE 5 — Timing
# ═════════════════════════════════════════════════════════════════════════════
def phase5(dc, sn):
    print("\n"+"="*70+"\n  PHASE 5: Timing\n"+"="*70)
    sd=[]
    for si in dc["snaps"]:
        nm=si["name"]; d=sn.get(nm,{}).get("min_ns",200.)
        sd.append({"name":nm,"ns":d})
    ts=sum(d["ns"] for d in sd)
    tr=dc["nr"]*ROT_NS; td=dc["nd"]*DISP_NS
    tt=ts+tr+td
    R={"snap":{"total":ts,"each":sd},
       "full":{"nr":dc["nr"],"nd":dc["nd"],"ns":dc["n_snap"],
               "tr":tr,"td":td,"ts":ts,"tt":tt}}
    print(f"  {dc['nr']}×Rq@{ROT_NS}ns = {tr:.0f}ns")
    print(f"  {dc['nd']}×D@{DISP_NS}ns  = {td:.0f}ns")
    for d in sd: print(f"  {d['name']}: {d['ns']:.0f}ns")
    print(f"  Total: {tt:.0f}ns"); return R

# ═════════════════════════════════════════════════════════════════════════════
# PHASE 6 — GRAPE full unitary
# ═════════════════════════════════════════════════════════════════════════════
def phase6(Ut):
    print("\n"+"="*70+"\n  PHASE 6: GRAPE Full-Unitary\n"+"="*70)
    mod=_model(); fr=_frame(mod); sub=_sub(); R={}
    for tns,lab in [(200,"200ns"),(400,"400ns")]:
        ns=max(20,tns//4); dt=tns*1e-9/ns; bf=0.
        for sd in GR_SEEDS:
            try:
                p=build_control_problem_from_model(
                    mod, frame=fr,
                    time_grid=PiecewiseConstantTimeGrid.uniform(steps=ns,dt_s=dt),
                    channel_specs=(
                        ModelControlChannelSpec(name="storage",target="storage",
                            quadratures=("I","Q"),amplitude_bounds=(-AMP_BND,AMP_BND)),
                        ModelControlChannelSpec(name="qubit",target="qubit",
                            quadratures=("I","Q"),amplitude_bounds=(-AMP_BND,AMP_BND)),),
                    objectives=(OCObj(target_operator=Ut, subspace=sub,
                        ignore_global_phase=True, name=f"full_{lab}"),),
                    penalties=(OCLeak(subspace=sub, weight=0.02),))
                g=GrapeSolver(GrapeConfig(maxiter=GR_ITER,seed=sd,random_scale=.3)).solve(p)
                f=float(g.metrics.get("nominal_fidelity",g.metrics.get("fidelity",0.)))
                if f>bf: bf=f
            except Exception as e: print(f"    sd={sd}@{lab}: {e}")
        R[lab]={"ns":tns,"fid":bf}
        print(f"  GRAPE {lab}: F={bf:.8f}")
    return R

# ═════════════════════════════════════════════════════════════════════════════
# Figures
# ═════════════════════════════════════════════════════════════════════════════
def figures(A):
    print("\n  Generating figures …")
    p2=A.get("p2",{}); p4=A.get("p4",{})

    # stabilisers
    st=p2.get("stab",[])
    if st:
        fig,ax=plt.subplots(figsize=(7,3.5))
        ax.plot([s["site"] for s in st],[s["K"] for s in st],"o-",color=COLORS[0],ms=7)
        ax.axhline(1,ls=":",color="gray",alpha=.5)
        ax.set(xlabel="Site $i$",ylabel=r"$\langle K_i\rangle$",
               title="Stabiliser Expectations",ylim=(0.9,1.05))
        for e in ("png","pdf"):
            fig.savefig(FIG/f"stabilisers.{e}",dpi=300,bbox_inches="tight")
        plt.close(fig)

    # SNAP sweep
    sk=[k for k in p4 if k.startswith("SNAP")]
    if sk:
        fig,ax=plt.subplots(figsize=(7,3.5))
        for i,k in enumerate(sk):
            sw=p4[k]["sweep"]
            ax.plot([s["dns"] for s in sw],[s["fid"] for s in sw],
                    "o-",color=COLORS[i%len(COLORS)],label=k)
        ax.axhline(SNAP_THR,ls="--",color="red",label=f"Thr={SNAP_THR}")
        ax.set(xlabel="Duration (ns)",ylabel="Fidelity",title="SNAP GRAPE Sweep")
        ax.legend()
        for e in ("png","pdf"):
            fig.savefig(FIG/f"snap_sweep.{e}",dpi=300,bbox_inches="tight")
        plt.close(fig)

    # Pauli
    pa=p2.get("pauli",[])
    if pa:
        fig,ax=plt.subplots(figsize=(7,3.5))
        s=[d["site"] for d in pa]
        ax.plot(s,[d["X"] for d in pa],"o-",color=COLORS[0],label=r"$\langle X\rangle$")
        ax.plot(s,[d["Y"] for d in pa],"s-",color=COLORS[1],label=r"$\langle Y\rangle$")
        ax.plot(s,[d["Z"] for d in pa],"^-",color=COLORS[2],label=r"$\langle Z\rangle$")
        ax.axhline(0,ls=":",color="gray",alpha=.5)
        ax.set(xlabel="Site",ylabel="Value",title="Single-Site Pauli")
        ax.legend()
        for e in ("png","pdf"):
            fig.savefig(FIG/f"pauli.{e}",dpi=300,bbox_inches="tight")
        plt.close(fig)

    # String-order heatmap
    ss=p2.get("strings",[]); nn=p2.get("n",6)
    if ss:
        M=np.full((nn,nn),np.nan)
        for e in ss: M[e["i"],e["j"]]=e["v"]
        fig,ax=plt.subplots(figsize=(5.5,4.5))
        im=ax.imshow(M,cmap="RdBu_r",vmin=-1,vmax=1,origin="upper")
        ax.set(xlabel="$j$",ylabel="$i$",title=r"String-Order $\langle Z_i\prod X_k\,Z_j\rangle$")
        plt.colorbar(im,ax=ax)
        for e in ("png","pdf"):
            fig.savefig(FIG/f"string_order.{e}",dpi=300,bbox_inches="tight")
        plt.close(fig)
    print("  Done.")

# ═════════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════════
def main():
    t0=time.perf_counter(); A={}

    r1, Ut, psi = phase1()
    A["p1"]=r1
    if not r1["passed"]:
        print("  ABORT"); _jdump(DATA/"results.json",A); return

    A["p2"] = phase2(psi, NSITES)
    r3 = phase3(Ut); A["p3"]={k:v for k,v in r3.items() if k!="_sub_U"}
    if "error" in r3:
        print("  ABORT"); _jdump(DATA/"results.json",A); return

    A["p4"] = phase4(r3["snaps"])
    A["p5"] = phase5(r3, A["p4"])
    A["p6"] = phase6(Ut)

    figures(A)
    _jdump(DATA/"results.json", A)

    dt=time.perf_counter()-t0
    print(f"\n{'='*70}\n  DONE in {dt:.0f}s\n{'='*70}")
    print(f"  Verification: {'PASS' if r1['passed'] else 'FAIL'}")
    print(f"  Decomposition: F={r3['fid']:.8f} ({r3['n_snap']} SNAP)")
    fs=A["p5"]["full"]
    print(f"  Total: {fs['tt']:.0f}ns (SNAP={fs['ts']:.0f}ns)")
    for k,v in A["p6"].items():
        if isinstance(v,dict): print(f"  GRAPE {k}: F={v.get('fid','?')}")

if __name__=="__main__":
    main()
