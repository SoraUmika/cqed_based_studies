"""
Cluster-state holographic simulation — final study script (v4).

Executes all phases:
  1. Target verification (matrix-level)
  2. Ideal observables on the canonical OBC cluster state
  3. Analytical decomposition analysis + SQR-based synthesis attempt
  4. GRAPE optimization for the full target at various durations
  5. Holographic correlator cross-check
  6. Figure generation

Key finding: the cqed_sim ideal SNAP is cavity-only (not qubit-conditional),
so D-R-SNAP cannot implement entangling gates. We use SQR-based decomposition
and direct GRAPE instead.
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
    Displacement as SDisp, QubitRotation as SRot, SNAP as SSNAP,
    SQR, Subspace, TargetUnitary, UnitarySynthesizer,
    subspace_unitary_fidelity,
)
from cqed_sim.unitary_synthesis.targets import make_target
from cqed_sim.quantum_algorithms import HolographicChannel, HolographicSampler

# ── device parameters ────────────────────────────────────────────────────────
OQ = TWO_PI*6.150e9;  OC = TWO_PI*5.241e9
AL = TWO_PI*(-255e6);  CH = TWO_PI*(-2.84e6)
CP = TWO_PI*(-21e3);   KR = TWO_PI*(-28e3)
NCAV = 8;  NTR = 2;  NSITES = 6

AMP_BND   = TWO_PI*50e6
GR_ITER   = 300
GR_SEEDS  = [17, 42, 73]
ROT_NS    = 16.0
DISP_NS   = 48.0
SEL_NS    = 200.0   # SQR / selective gate duration

# ── JSON helper ──────────────────────────────────────────────────────────────
def _jdump(path, obj):
    def _d(o):
        if isinstance(o, np.ndarray): return o.tolist()
        if isinstance(o, np.floating): return float(o)
        if isinstance(o, np.integer):  return int(o)
        if isinstance(o, (complex, np.complexfloating)):
            return {"re": float(np.real(o)), "im": float(np.imag(o))}
        return str(o)
    path.write_text(json.dumps(obj, indent=2, default=_d), encoding="utf-8")

# ── builders ─────────────────────────────────────────────────────────────────
def _model(nc=NCAV, nt=NTR):
    return DispersiveTransmonCavityModel(
        omega_c=OC, omega_q=OQ, alpha=AL,
        chi=CH, chi_higher=(CP,), kerr=KR, n_cav=nc, n_tr=nt)

def _frame(m):
    return FrameSpec(omega_c_frame=m.omega_c, omega_q_frame=m.omega_q)

def _sub(nc=NCAV):
    return Subspace.custom(NTR*nc, (0,1,nc,nc+1),
                           ("|g,0>","|g,1>","|e,0>","|e,1>"))

def _ufid(U, V):
    d = U.shape[0]
    return float(abs(np.trace(U.conj().T @ V))**2 / d**2)


# ═════════════════════════════════════════════════════════════════════════════
# PHASE 1 — Target verification
# ═════════════════════════════════════════════════════════════════════════════
def phase1():
    print("\n" + "="*70 + "\n  PHASE 1: Target Verification\n" + "="*70)
    R = {}
    U = make_target("cluster", n_match=1)
    R["shape"] = list(U.shape)
    R["unitarity_err"] = float(np.linalg.norm(U.conj().T @ U - np.eye(4)))
    print(f"  ||U†U - I|| = {R['unitarity_err']:.2e}")

    # Build reference: SWAP · CZ · (H ⊗ I)
    H = np.array([[1,1],[1,-1]], dtype=complex) / np.sqrt(2)
    CZ = np.diag([1.,1.,1.,-1.]).astype(complex)
    SW = np.array([[1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]], dtype=complex)
    U_ref = SW @ CZ @ np.kron(H, np.eye(2))

    R["match_swap"] = bool(np.allclose(U, U_ref, atol=1e-12))
    R["fid_swap"] = _ufid(U, U_ref)
    print(f"  Matches SWAP·CZ·(H⊗I): {R['match_swap']} (F={R['fid_swap']:.10f})")

    # MPS tensors (CZ·(H⊗I) convention)
    V = CZ @ np.kron(H, np.eye(2))
    A0 = np.eye(2, dtype=complex) / np.sqrt(2)
    A1 = np.diag([1., -1.]).astype(complex) / np.sqrt(2)
    A0_ex = np.zeros((2,2), dtype=complex)
    A1_ex = np.zeros((2,2), dtype=complex)
    for a in range(2):
        for b in range(2):
            A0_ex[b,a] = V[0*2+b, a]
            A1_ex[b,a] = V[1*2+b, a]
    R["mps_A0_eq_I"] = bool(np.allclose(A0_ex, A0))
    R["mps_A1_eq_Z"] = bool(np.allclose(A1_ex, A1))
    print(f"  MPS: A^0=I/√2: {R['mps_A0_eq_I']}, A^1=Z/√2: {R['mps_A1_eq_Z']}")

    # HolographicChannel
    ch = HolographicChannel.from_unitary(U, physical_dim=2)
    R["holo_bond_dim"] = 2**ch.num_bond_qubits
    R["holo_phys_dim"] = 2**ch.num_physical_qubits
    R["holo_ref_state"] = ch.reference_state.tolist()
    print(f"  HolographicChannel: bond={R['holo_bond_dim']}, phys={R['holo_phys_dim']}")

    R["passed"] = bool(R["match_swap"] and R["mps_A0_eq_I"] and R["mps_A1_eq_Z"])
    print(f"  PHASE 1: {'PASS' if R['passed'] else 'FAIL'}")
    return R, U


# ═════════════════════════════════════════════════════════════════════════════
# PHASE 2 — Ideal observables
# ═════════════════════════════════════════════════════════════════════════════
def phase2(n):
    print("\n" + "="*70 + "\n  PHASE 2: Ideal Observables (N=%d)\n" % n + "="*70)

    # Build canonical OBC cluster state
    H = np.array([[1,1],[1,-1]], dtype=complex) / np.sqrt(2)
    dim = 2**n
    Hfull = np.array([1.], dtype=complex)
    for _ in range(n):
        Hfull = np.kron(Hfull, H)
    psi = Hfull @ np.eye(dim, 1, dtype=complex).ravel()  # H^N |0...0>
    # Apply CZ gates
    for i in range(n-1):
        for idx in range(dim):
            bits = [(idx >> (n-1-k)) & 1 for k in range(n)]
            if bits[i] == 1 and bits[i+1] == 1:
                psi[idx] *= -1

    I2 = np.eye(2, dtype=complex)
    X = np.array([[0,1],[1,0]], dtype=complex)
    Y = np.array([[0,-1j],[1j,0]], dtype=complex)
    Z = np.array([[1,0],[0,-1]], dtype=complex)

    def _nq(ops):
        m = np.array([1.], dtype=complex)
        for o in ops: m = np.kron(m, o)
        return m

    def _ev(op): return float(np.real(psi.conj() @ op @ psi))

    R = {"n": n}

    # Single-site Pauli
    pauli = []
    for i in range(n):
        ox=[I2]*n; ox[i]=X; oy=[I2]*n; oy[i]=Y; oz=[I2]*n; oz[i]=Z
        pauli.append({"site":i, "X":_ev(_nq(ox)), "Y":_ev(_nq(oy)), "Z":_ev(_nq(oz))})
    R["pauli"] = pauli
    R["max_pauli"] = max(max(abs(r["X"]),abs(r["Y"]),abs(r["Z"])) for r in pauli)
    print(f"  max |<sigma>| = {R['max_pauli']:.2e}")

    # Stabilisers K_i
    stab = []
    for i in range(n):
        ops = [I2]*n; ops[i] = X
        if i > 0:   ops[i-1] = Z
        if i < n-1: ops[i+1] = Z
        ki = _ev(_nq(ops))
        stab.append({"site":i, "K":ki})
    R["stab"] = stab
    R["max_stab_dev"] = max(abs(s["K"]-1.) for s in stab)
    print(f"  max |<Ki>-1| = {R['max_stab_dev']:.2e}")

    # ZXZ correlators (bulk)
    zxz = []
    for i in range(1, n-1):
        ops = [I2]*n; ops[i-1]=Z; ops[i]=X; ops[i+1]=Z
        zxz.append({"i":i, "v":_ev(_nq(ops))})
    R["zxz"] = zxz

    # String-order correlators
    strings = []
    for i in range(n):
        for j in range(i+2, n):
            ops = [I2]*n; ops[i]=Z; ops[j]=Z
            for k in range(i+1, j): ops[k] = X
            strings.append({"i":i, "j":j, "v":_ev(_nq(ops))})
    R["strings"] = strings
    # Print summary
    nn_str = [s for s in strings if s["j"]-s["i"]==2]
    nnn_str = [s for s in strings if s["j"]-s["i"]==3]
    print(f"  NN string-order: {[f'{s['v']:+.6f}' for s in nn_str[:4]]}")
    print(f"  NNN string-order: {[f'{s['v']:+.6f}' for s in nnn_str[:4]]}")

    # Two-point correlators
    twopt = []
    for i in range(n):
        for j in range(i+1, n):
            ops_zz = [I2]*n; ops_zz[i]=Z; ops_zz[j]=Z
            ops_xx = [I2]*n; ops_xx[i]=X; ops_xx[j]=X
            twopt.append({"i":i, "j":j, "ZZ":_ev(_nq(ops_zz)), "XX":_ev(_nq(ops_xx))})
    R["twopt"] = twopt
    R["max_twopt"] = max(max(abs(t["ZZ"]),abs(t["XX"])) for t in twopt)
    print(f"  max |two-point| = {R['max_twopt']:.2e}")

    print("  PHASE 2 COMPLETE")
    return R


# ═════════════════════════════════════════════════════════════════════════════
# PHASE 3 — Decomposition analysis
# ═════════════════════════════════════════════════════════════════════════════
def phase3(Ut):
    print("\n" + "="*70 + "\n  PHASE 3: Decomposition Analysis\n" + "="*70)
    R = {}
    sub = _sub()

    # 3a. Document ideal SNAP = cavity-only
    s = SSNAP(name='test', phases=[0.,np.pi]+[0.]*(NCAV-2),
              duration=SEL_NS*1e-9, optimize_time=False)
    U_snap = s.ideal_unitary(n_cav=NCAV, n_tr=NTR).full()
    diag_snap = np.diag(U_snap)
    # Check: g-manifold and e-manifold get SAME phases
    g_phases = [float(np.angle(diag_snap[n])) for n in range(NCAV)]
    e_phases = [float(np.angle(diag_snap[NCAV+n])) for n in range(NCAV)]
    R["snap_cavity_only"] = bool(np.allclose(g_phases, e_phases, atol=1e-10))
    print(f"  SNAP is cavity-only (g==e phases): {R['snap_cavity_only']}")
    R["snap_g_phases"] = g_phases[:4]
    R["snap_e_phases"] = e_phases[:4]

    # 3b. Analytical decomposition: U = SWAP · CZ · (H ⊗ I)
    # In the qubit-conditional SNAP convention:
    #   CZ = SNAP_QC(0, pi) where SNAP_QC(θ0,θ1) = diag(1,1,e^{iθ0},e^{iθ1})
    #   H_q = R_q(π/2, 0) (Hadamard on qubit)
    #   SWAP = 3 × CNOT = 3 × (R ⊗ D_H · CZ · R ⊗ D_H)
    # Minimum: 4 selective layers + rotations + displacements
    R["analytical"] = {
        "target": "SWAP · CZ · (H ⊗ I)",
        "CZ_as_SNAP_QC": "diag(1,1,1,-1) = SNAP_QC(0, π)",
        "H_as_Rq": "R_q(π/2, 0)",
        "SWAP_min_CNOT": 3,
        "min_selective_layers": 4,
        "note": "cqed_sim ideal SNAP is I_q ⊗ S_c (cavity-only), not qubit-conditional"
    }

    # 3c. SQR-based numerical synthesis (SQR IS qubit-conditional)
    Rs = ROT_NS*1e-9; Ds = DISP_NS*1e-9; Ss = SEL_NS*1e-9
    tgt = TargetUnitary(Ut, ignore_global_phase=True)

    best = None; bfid = 0; blab = ""
    for nlayers, label in [(2, "2-SQR"), (3, "3-SQR")]:
        print(f"  SQR synthesis: {label} ...")
        gl = []
        for k in range(nlayers):
            gl.append(SRot(name=f'R{2*k+1}', theta=np.pi/2, phi=0.,
                           duration=Rs, optimize_time=False))
            gl.append(SQR(name=f'SQR{k+1}',
                          theta_n=[0.1]*NCAV, phi_n=[0.]*NCAV,
                          duration=Ss, optimize_time=False))
            gl.append(SRot(name=f'R{2*k+2}', theta=np.pi/2, phi=np.pi/2,
                           duration=Rs, optimize_time=False))
            if k < nlayers - 1:
                gl.append(SDisp(name=f'D{k+1}', alpha=.2+0j,
                                duration=Ds, optimize_time=False))
        t0 = time.time()
        try:
            syn = UnitarySynthesizer(
                primitives=gl, subspace=sub, target=tgt,
                seed=42, optimize_times=False, optimizer='powell')
            sr = syn.fit(multistart=3, maxiter=300)
            fid = float(sr.report.get('metrics',{}).get('fidelity', 1-sr.objective))
            dt = time.time() - t0
            print(f"    F={fid:.6f} obj={sr.objective:.6e} time={dt:.1f}s")
            R[label] = {"fidelity": fid, "objective": float(sr.objective), "time_s": dt}
            if fid > bfid:
                bfid = fid; best = sr; blab = label
        except Exception as e:
            print(f"    FAILED: {e}")
            R[label] = {"error": str(e)}

    R["best_label"] = blab
    R["best_fidelity"] = bfid
    if best is not None:
        gates = []
        for g in best.sequence.gates:
            gt = type(g).__name__
            e = {"name": g.name, "type": gt}
            if gt == "SQR":
                e["theta_n"] = [float(t) for t in g.theta_n[:4]]
                e["phi_n"] = [float(p) for p in g.phi_n[:4]]
            elif gt == "QubitRotation":
                e["theta"] = float(g.theta); e["phi"] = float(g.phi)
            elif gt == "Displacement":
                e["alpha_re"] = float(np.real(g.alpha))
                e["alpha_im"] = float(np.imag(g.alpha))
            gates.append(e)
        R["best_gates"] = gates

    print(f"  Best: {blab} F={bfid:.6f}")
    print("  PHASE 3 COMPLETE")
    return R


# ═════════════════════════════════════════════════════════════════════════════
# PHASE 4 — GRAPE full-target optimization
# ═════════════════════════════════════════════════════════════════════════════
def phase4(Ut):
    print("\n" + "="*70 + "\n  PHASE 4: GRAPE Full-Target Optimization\n" + "="*70)
    mod = _model(); fr = _frame(mod); sub = _sub()
    R = {}

    durations_ns = [100, 200, 300, 400, 600]
    for tns in durations_ns:
        ns = max(20, tns // 4)
        dt = tns * 1e-9 / ns
        bf = 0.; br = None
        t0 = time.time()
        for sd in GR_SEEDS:
            try:
                p = build_control_problem_from_model(
                    mod, frame=fr,
                    time_grid=PiecewiseConstantTimeGrid.uniform(steps=ns, dt_s=dt),
                    channel_specs=(
                        ModelControlChannelSpec(name="storage", target="storage",
                            quadratures=("I","Q"),
                            amplitude_bounds=(-AMP_BND, AMP_BND)),
                        ModelControlChannelSpec(name="qubit", target="qubit",
                            quadratures=("I","Q"),
                            amplitude_bounds=(-AMP_BND, AMP_BND)),),
                    objectives=(OCObj(target_operator=Ut, subspace=sub,
                        ignore_global_phase=True, name=f"cluster_{tns}ns"),),
                    penalties=(OCLeak(subspace=sub, weight=0.02),))
                g = GrapeSolver(GrapeConfig(
                    maxiter=GR_ITER, seed=sd, random_scale=0.3)).solve(p)
                m = g.metrics
                f = float(m.get("nominal_fidelity", m.get("fidelity", 0.)))
                leak = float(m.get("leakage_average", 0.))
                if f > bf:
                    bf = f
                    br = {"seed": sd, "fid": f, "leak": leak,
                          "obj": float(g.objective_value)}
            except Exception as e:
                print(f"    sd={sd}@{tns}ns: {e}")
        dt_wall = time.time() - t0
        entry = {"dns": tns, "fid": bf, "wall_s": dt_wall}
        if br:
            entry.update(br)
        R[f"{tns}ns"] = entry
        ok = bf >= 0.999
        print(f"  {tns:4d}ns: F={bf:.6f} leak={br.get('leak',0):.2e} "
              f"{'PASS' if ok else ''} ({dt_wall:.0f}s)")

    print("  PHASE 4 COMPLETE")
    return R


# ═════════════════════════════════════════════════════════════════════════════
# PHASE 5 — Holographic channel analysis
# ═════════════════════════════════════════════════════════════════════════════
def phase5(Ut):
    print("\n" + "="*70 + "\n  PHASE 5: Holographic Channel Analysis\n" + "="*70)
    R = {}

    ch = HolographicChannel.from_unitary(Ut, physical_dim=2)
    mps = ch.mps_matrices
    ref = ch.reference_state

    # Report MPS matrices
    R["mps_A0"] = mps[0].tolist()
    R["mps_A1"] = mps[1].tolist()
    R["ref_state"] = ref.tolist()
    R["bond_dim"] = 2**ch.num_bond_qubits

    # Transfer matrix eigenvalues (characterise correlations)
    # T = sum_s A^s_conj ⊗ A^s
    d = mps[0].shape[0]
    T = np.zeros((d**2, d**2), dtype=complex)
    for s in range(len(mps)):
        T += np.kron(mps[s].conj(), mps[s])
    evals = np.linalg.eigvals(T)
    idx = np.argsort(-np.abs(evals))
    evals_sorted = evals[idx]
    R["transfer_eigenvalues"] = [{"re": float(np.real(e)), "im": float(np.imag(e)),
                                  "abs": float(np.abs(e))} for e in evals_sorted]
    print(f"  Transfer matrix eigenvalues (|λ|):")
    for i, e in enumerate(evals_sorted):
        print(f"    λ_{i} = {np.real(e):+.6f}{np.imag(e):+.6f}i  (|λ|={np.abs(e):.6f})")

    # Correlation length
    if len(evals_sorted) >= 2 and np.abs(evals_sorted[1]) > 1e-12:
        xi = -1.0 / np.log(np.abs(evals_sorted[1]) / np.abs(evals_sorted[0]))
        R["correlation_length"] = float(xi)
        print(f"  Correlation length ξ = {xi:.4f} sites")
    else:
        R["correlation_length"] = float('inf')
        print(f"  Correlation length ξ = ∞ (degenerate leading eigenvalue)")

    # Quantum channel check: E(ρ) = Σ A^s ρ A^{s†}, completeness: Σ A^{s†}A^s = I
    comp_err = np.linalg.norm(sum(m.conj().T @ m for m in mps) - np.eye(d))
    R["completeness_error"] = float(comp_err)
    print(f"  Completeness error ||Σ A†A - I|| = {comp_err:.2e}")

    # Fixed-point analysis
    rho_fixed = np.outer(ref, ref.conj())
    rho_out = sum(mps[s] @ rho_fixed @ mps[s].conj().T for s in range(len(mps)))
    fix_err = np.linalg.norm(rho_out - rho_fixed)
    R["fixed_point_err"] = float(fix_err)
    is_fixed = fix_err < 1e-10
    R["ref_is_fixed_point"] = is_fixed
    print(f"  |E(ρ_ref) - ρ_ref| = {fix_err:.2e} (fixed point: {is_fixed})")

    print("  PHASE 5 COMPLETE")
    return R


# ═════════════════════════════════════════════════════════════════════════════
# FIGURES
# ═════════════════════════════════════════════════════════════════════════════
def generate_figures(A):
    print("\n  Generating figures ...")
    p2 = A.get("p2", {}); p4 = A.get("p4", {})

    # Fig 1: Stabiliser expectations
    stab = p2.get("stab", [])
    if stab:
        fig, ax = plt.subplots(figsize=(6, 3))
        sites = [s["site"] for s in stab]
        vals = [s["K"] for s in stab]
        ax.bar(sites, vals, color=COLORS[0], alpha=0.8, edgecolor='black', linewidth=0.5)
        ax.axhline(1, ls=":", color="gray", alpha=0.5)
        ax.set(xlabel="Site $i$", ylabel=r"$\langle K_i \rangle$",
               title=f"Stabiliser Expectations ($N={p2.get('n',6)}$)",
               ylim=(0.95, 1.05))
        for fmt in ("png","pdf"):
            fig.savefig(FIG / f"stabilisers.{fmt}", dpi=300, bbox_inches="tight")
        plt.close(fig)

    # Fig 2: Single-site Pauli
    pauli = p2.get("pauli", [])
    if pauli:
        fig, ax = plt.subplots(figsize=(6, 3))
        sites = [d["site"] for d in pauli]
        for op, mk, c, lab in [("X","o",COLORS[0],r"$\langle X\rangle$"),
                                 ("Y","s",COLORS[1],r"$\langle Y\rangle$"),
                                 ("Z","^",COLORS[2],r"$\langle Z\rangle$")]:
            ax.plot(sites, [d[op] for d in pauli], mk+"-", color=c, ms=6, label=lab)
        ax.axhline(0, ls=":", color="gray", alpha=0.5)
        ax.set(xlabel="Site $i$", ylabel="Expectation value",
               title="Single-Site Pauli Expectations")
        ax.legend(); ax.set_ylim(-0.3, 0.3)
        for fmt in ("png","pdf"):
            fig.savefig(FIG / f"pauli.{fmt}", dpi=300, bbox_inches="tight")
        plt.close(fig)

    # Fig 3: String-order heatmap
    strings = p2.get("strings", []); nn = p2.get("n", 6)
    if strings:
        M = np.full((nn, nn), np.nan)
        for e in strings:
            M[e["i"], e["j"]] = e["v"]
        fig, ax = plt.subplots(figsize=(5, 4))
        im = ax.imshow(M, cmap="RdBu_r", vmin=-1, vmax=1, origin="upper")
        ax.set(xlabel="$j$", ylabel="$i$",
               title=r"String-Order $\langle Z_i \prod_{k} X_k\, Z_j \rangle$")
        plt.colorbar(im, ax=ax, label="Correlator value")
        for fmt in ("png","pdf"):
            fig.savefig(FIG / f"string_order.{fmt}", dpi=300, bbox_inches="tight")
        plt.close(fig)

    # Fig 4: GRAPE fidelity vs duration
    if p4:
        entries = sorted([v for v in p4.values() if isinstance(v, dict) and "dns" in v],
                         key=lambda x: x["dns"])
        if entries:
            fig, ax = plt.subplots(figsize=(6, 3.5))
            dns = [e["dns"] for e in entries]
            fids = [e["fid"] for e in entries]
            ax.plot(dns, fids, "o-", color=COLORS[0], ms=7, linewidth=2)
            ax.axhline(0.999, ls="--", color="red", alpha=0.7, label="99.9% threshold")
            ax.set(xlabel="Pulse Duration (ns)", ylabel="Fidelity",
                   title="GRAPE Fidelity vs. Duration")
            ax.legend()
            ax.set_ylim(min(fids)-0.02, 1.005)
            for fmt in ("png","pdf"):
                fig.savefig(FIG / f"grape_fidelity.{fmt}", dpi=300, bbox_inches="tight")
            plt.close(fig)

    # Fig 5: Transfer matrix eigenvalue spectrum
    p5 = A.get("p5", {})
    te = p5.get("transfer_eigenvalues", [])
    if te:
        fig, ax = plt.subplots(figsize=(5, 5))
        theta = np.linspace(0, 2*np.pi, 100)
        ax.plot(np.cos(theta), np.sin(theta), 'k-', alpha=0.2, linewidth=0.5)
        for i, e in enumerate(te):
            ax.plot(e["re"], e["im"], 'o', color=COLORS[i%len(COLORS)],
                    ms=10, label=f"|λ_{i}|={e['abs']:.4f}")
        ax.set(xlabel="Re(λ)", ylabel="Im(λ)",
               title="Transfer Matrix Eigenvalues", aspect='equal')
        ax.legend(fontsize=8)
        ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.3, 1.3)
        for fmt in ("png","pdf"):
            fig.savefig(FIG / f"transfer_eigenvalues.{fmt}", dpi=300, bbox_inches="tight")
        plt.close(fig)

    print("  Done.")


# ═════════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════════
def main():
    t0 = time.perf_counter()
    A = {}

    # Phase 1
    r1, Ut = phase1()
    A["p1"] = r1
    if not r1["passed"]:
        print("  ABORT: Phase 1 failed")
        _jdump(DATA / "results.json", A)
        return

    # Phase 2
    A["p2"] = phase2(NSITES)

    # Phase 3
    A["p3"] = phase3(Ut)

    # Phase 4
    A["p4"] = phase4(Ut)

    # Phase 5
    A["p5"] = phase5(Ut)

    # Figures
    generate_figures(A)

    # Save
    _jdump(DATA / "results.json", A)

    dt = time.perf_counter() - t0
    print(f"\n{'='*70}\n  ALL PHASES COMPLETE in {dt:.0f}s\n{'='*70}")
    print(f"  Phase 1: {'PASS' if r1['passed'] else 'FAIL'}")
    p3 = A.get("p3", {})
    print(f"  Phase 3: SQR decomp best={p3.get('best_label','-')} F={p3.get('best_fidelity',0):.6f}")
    p4 = A.get("p4", {})
    for k, v in p4.items():
        if isinstance(v, dict) and "fid" in v:
            print(f"  Phase 4 GRAPE {k}: F={v['fid']:.6f}")
    p5 = A.get("p5", {})
    print(f"  Phase 5: ξ={p5.get('correlation_length','?')} sites")


if __name__ == "__main__":
    main()
