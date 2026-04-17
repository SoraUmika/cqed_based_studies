"""
Completion script: GRAPE fidelity vs duration sweep + holographic analysis +
figures. Loads existing results from v3 and augments with new data.
"""
from __future__ import annotations
import json, os, sys, time, traceback
from pathlib import Path
import numpy as np
os.environ["PYTHONIOENCODING"] = "utf-8"
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

STUDY = Path(__file__).resolve().parents[1]
DATA  = STUDY / "data"; DATA.mkdir(exist_ok=True)
FIG   = STUDY / "figures"; FIG.mkdir(exist_ok=True)
SIM   = Path("C:/Users/dazzl/Box/Shyam Shankar Quantum Circuits Group"
             "/Users/Users_JianJun/cQED_simulation")
if str(SIM) not in sys.path:
    sys.path.insert(0, str(SIM))
TWO_PI = 2*np.pi
COLORS = ['#4477AA','#EE6677','#228833','#CCBB44','#66CCEE','#AA3377','#BBBBBB']

from cqed_sim import (
    DispersiveTransmonCavityModel, FrameSpec,
    ModelControlChannelSpec, PiecewiseConstantTimeGrid,
)
from cqed_sim.optimal_control import (
    GrapeConfig, GrapeSolver,
    LeakagePenalty as OCLeak,
    UnitaryObjective as OCObj,
    build_control_problem_from_model,
)
from cqed_sim.unitary_synthesis import Subspace
from cqed_sim.unitary_synthesis.targets import make_target
from cqed_sim.quantum_algorithms import HolographicChannel

OQ=TWO_PI*6.150e9; OC=TWO_PI*5.241e9
AL=TWO_PI*(-255e6); CH=TWO_PI*(-2.84e6)
CP=TWO_PI*(-21e3);  KR=TWO_PI*(-28e3)
NCAV=8; NTR=2; AMP=TWO_PI*50e6
GR_ITER=300; GR_SEEDS=[17,42,73]

def _model():
    return DispersiveTransmonCavityModel(
        omega_c=OC, omega_q=OQ, alpha=AL,
        chi=CH, chi_higher=(CP,), kerr=KR, n_cav=NCAV, n_tr=NTR)

def _frame(m):
    return FrameSpec(omega_c_frame=m.omega_c, omega_q_frame=m.omega_q)

def _sub():
    return Subspace.custom(NTR*NCAV, (0,1,NCAV,NCAV+1),
                           ("|g,0>","|g,1>","|e,0>","|e,1>"))

def _jdump(path, obj):
    def _d(o):
        if isinstance(o, np.ndarray): return o.tolist()
        if isinstance(o, np.floating): return float(o)
        if isinstance(o, np.integer): return int(o)
        if isinstance(o, (complex, np.complexfloating)):
            return {"re": float(np.real(o)), "im": float(np.imag(o))}
        return str(o)
    path.write_text(json.dumps(obj, indent=2, default=_d), encoding="utf-8")


# ═══════════════════════════════════════════════════════════════════════════
# GRAPE sweep
# ═══════════════════════════════════════════════════════════════════════════
def grape_sweep(Ut):
    print("="*60 + "\n  GRAPE Fidelity vs Duration Sweep\n" + "="*60)
    mod = _model(); fr = _frame(mod); sub = _sub()
    results = []
    durations = [50, 100, 150, 200, 300, 400, 500, 600, 800]
    for tns in durations:
        ns = max(20, tns // 5)
        dt = tns * 1e-9 / ns
        bf = 0.; bleak = 0.; bobj = 0.
        t0 = time.time()
        for sd in GR_SEEDS:
            try:
                p = build_control_problem_from_model(
                    mod, frame=fr,
                    time_grid=PiecewiseConstantTimeGrid.uniform(steps=ns, dt_s=dt),
                    channel_specs=(
                        ModelControlChannelSpec(name="storage", target="storage",
                            quadratures=("I","Q"),
                            amplitude_bounds=(-AMP, AMP)),
                        ModelControlChannelSpec(name="qubit", target="qubit",
                            quadratures=("I","Q"),
                            amplitude_bounds=(-AMP, AMP)),),
                    objectives=(OCObj(target_operator=Ut, subspace=sub,
                        ignore_global_phase=True, name=f"cluster_{tns}ns"),),
                    penalties=(OCLeak(subspace=sub, weight=0.02),))
                g = GrapeSolver(GrapeConfig(
                    maxiter=GR_ITER, seed=sd, random_scale=0.3)).solve(p)
                m = g.metrics
                f = float(m.get("nominal_fidelity", m.get("fidelity", 0.)))
                lk = float(m.get("leakage_average", 0.))
                if f > bf:
                    bf = f; bleak = lk; bobj = float(g.objective_value)
            except Exception as e:
                print(f"    (sd={sd}@{tns}ns: {e})")
        wall = time.time() - t0
        entry = {"dns": tns, "fid": bf, "leak": bleak, "obj": bobj, "wall_s": wall}
        results.append(entry)
        ok = "PASS" if bf >= 0.999 else ""
        print(f"  {tns:4d}ns: F={bf:.6f}  leak={bleak:.2e}  {ok}  ({wall:.0f}s)")
    return results


# ═══════════════════════════════════════════════════════════════════════════
# Holographic channel analysis
# ═══════════════════════════════════════════════════════════════════════════
def holographic_analysis(Ut):
    print("\n" + "="*60 + "\n  Holographic Channel Analysis\n" + "="*60)
    R = {}
    ch = HolographicChannel.from_unitary(Ut, physical_dim=2)
    mps = ch.mps_matrices
    ref = ch.reference_state
    R["mps_A0"] = mps[0].tolist()
    R["mps_A1"] = mps[1].tolist()
    R["ref_state"] = ref.tolist()
    R["bond_dim"] = 2**ch.num_bond_qubits
    d = mps[0].shape[0]

    # Transfer matrix
    T = np.zeros((d**2, d**2), dtype=complex)
    for s in range(len(mps)):
        T += np.kron(mps[s].conj(), mps[s])
    evals = np.linalg.eigvals(T)
    idx = np.argsort(-np.abs(evals))
    evals_s = evals[idx]
    R["transfer_evals"] = [{"re": float(np.real(e)), "im": float(np.imag(e)),
                             "abs": float(np.abs(e))} for e in evals_s]
    print("  Transfer matrix eigenvalues:")
    for i, e in enumerate(evals_s):
        print(f"    lambda_{i} = {np.real(e):+.6f}{np.imag(e):+.6f}i  "
              f"(|lambda|={np.abs(e):.6f})")

    # Correlation length
    if len(evals_s) >= 2 and np.abs(evals_s[1]) > 1e-12:
        xi = -1.0 / np.log(np.abs(evals_s[1]) / np.abs(evals_s[0]))
        R["correlation_length"] = float(xi)
    else:
        R["correlation_length"] = float('inf')
    print(f"  Correlation length: {R['correlation_length']:.4f} sites")

    # Completeness
    comp = np.linalg.norm(sum(m.conj().T @ m for m in mps) - np.eye(d))
    R["completeness_err"] = float(comp)
    print(f"  Completeness ||sum A^dag A - I|| = {comp:.2e}")

    # Fixed point
    rho = np.outer(ref, ref.conj())
    rho_out = sum(mps[s] @ rho @ mps[s].conj().T for s in range(len(mps)))
    fe = np.linalg.norm(rho_out - rho)
    R["fixed_point_err"] = float(fe)
    R["ref_is_fixed_point"] = bool(fe < 1e-10)
    print(f"  |E(rho_ref) - rho_ref| = {fe:.2e} (fixed: {R['ref_is_fixed_point']})")

    # MPS in alternative (no-SWAP) convention: check A^0=I/sqrt(2), A^1=Z/sqrt(2)
    H = np.array([[1,1],[1,-1]], dtype=complex)/np.sqrt(2)
    CZ = np.diag([1.,1.,1.,-1.]).astype(complex)
    V = CZ @ np.kron(H, np.eye(2))
    A0_ns = np.array([[V[b, a] for a in range(2)] for b in range(2)])
    A1_ns = np.array([[V[2+b, a] for a in range(2)] for b in range(2)])
    A0_ref = np.eye(2, dtype=complex)/np.sqrt(2)
    A1_ref = np.diag([1.,-1.]).astype(complex)/np.sqrt(2)
    R["noswap_A0_eq_I"] = bool(np.allclose(A0_ns, A0_ref))
    R["noswap_A1_eq_Z"] = bool(np.allclose(A1_ns, A1_ref))
    print(f"  No-SWAP MPS: A^0=I/sqrt2: {R['noswap_A0_eq_I']}, "
          f"A^1=Z/sqrt2: {R['noswap_A1_eq_Z']}")

    return R


# ═══════════════════════════════════════════════════════════════════════════
# Figures
# ═══════════════════════════════════════════════════════════════════════════
def make_figures(grape_data, holo_data, old_results):
    print("\n  Generating figures ...")
    p2 = old_results.get("p2", {})
    nn = p2.get("n", 6)

    # 1. Stabiliser bar chart
    stab = p2.get("stab", [])
    if stab:
        fig, ax = plt.subplots(figsize=(5.5, 3))
        sites = [s["site"] for s in stab]
        vals = [s["K"] for s in stab]
        ax.bar(sites, vals, color=COLORS[0], alpha=0.85, edgecolor='black', lw=0.5)
        ax.axhline(1, ls=":", color="gray", alpha=0.5)
        ax.set(xlabel="Site $i$", ylabel=r"$\langle K_i \rangle$",
               title=f"Cluster-State Stabiliser Expectations ($N={nn}$)")
        ax.set_ylim(0.95, 1.02)
        for f in ("png","pdf"):
            fig.savefig(FIG/f"stabilisers.{f}", dpi=300, bbox_inches="tight")
        plt.close(fig)

    # 2. String-order heatmap
    strings = p2.get("strings", [])
    if strings:
        M = np.full((nn, nn), np.nan)
        for e in strings: M[e["i"], e["j"]] = e["v"]
        fig, ax = plt.subplots(figsize=(4.5, 3.8))
        im = ax.imshow(M, cmap="RdBu_r", vmin=-1, vmax=1, origin="upper")
        ax.set(xlabel="$j$", ylabel="$i$",
               title=r"String-Order $\langle Z_i \prod_k X_k Z_j \rangle$")
        plt.colorbar(im, ax=ax, shrink=0.85)
        for f in ("png","pdf"):
            fig.savefig(FIG/f"string_order.{f}", dpi=300, bbox_inches="tight")
        plt.close(fig)

    # 3. GRAPE fidelity vs duration
    if grape_data:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
        dns = [e["dns"] for e in grape_data]
        fids = [e["fid"] for e in grape_data]
        leaks = [e["leak"] for e in grape_data]

        ax1.plot(dns, fids, "o-", color=COLORS[0], ms=7, lw=2)
        ax1.axhline(0.999, ls="--", color="red", alpha=0.6, label="$99.9\\%$")
        ax1.axhline(0.99, ls=":", color="orange", alpha=0.6, label="$99\\%$")
        ax1.set(xlabel="Pulse Duration (ns)", ylabel="Fidelity",
                title="GRAPE Fidelity vs Duration")
        ax1.legend(fontsize=9)
        ax1.set_ylim(min(fids) - 0.02, 1.005)

        ax2.semilogy(dns, leaks, "s-", color=COLORS[1], ms=6, lw=2)
        ax2.set(xlabel="Pulse Duration (ns)", ylabel="Average Leakage",
                title="Leakage vs Duration")

        fig.tight_layout()
        for f in ("png","pdf"):
            fig.savefig(FIG/f"grape_fidelity.{f}", dpi=300, bbox_inches="tight")
        plt.close(fig)

    # 4. Transfer matrix eigenvalue plot
    te = holo_data.get("transfer_evals", [])
    if te:
        fig, ax = plt.subplots(figsize=(4.5, 4.5))
        theta = np.linspace(0, 2*np.pi, 100)
        ax.plot(np.cos(theta), np.sin(theta), 'k-', alpha=0.15, lw=0.5)
        for i, e in enumerate(te):
            ax.plot(e["re"], e["im"], 'o', color=COLORS[i % len(COLORS)],
                    ms=10, label=f"$|\\lambda_{i}|={e['abs']:.3f}$")
        ax.set(xlabel="Re($\\lambda$)", ylabel="Im($\\lambda$)",
               title="Transfer Matrix Eigenvalues", aspect='equal')
        ax.legend(fontsize=9); ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.3, 1.3)
        for f in ("png","pdf"):
            fig.savefig(FIG/f"transfer_eigenvalues.{f}", dpi=300, bbox_inches="tight")
        plt.close(fig)

    # 5. Pauli expectations
    pauli = p2.get("pauli", [])
    if pauli:
        fig, ax = plt.subplots(figsize=(5.5, 3))
        sites = [d["site"] for d in pauli]
        for op, mk, c, lb in [("X","o",COLORS[0],r"$\langle X\rangle$"),
                                ("Y","s",COLORS[1],r"$\langle Y\rangle$"),
                                ("Z","^",COLORS[2],r"$\langle Z\rangle$")]:
            ax.plot(sites, [d[op] for d in pauli], mk+"-", color=c, ms=5, label=lb)
        ax.axhline(0, ls=":", color="gray", alpha=0.5)
        ax.set(xlabel="Site", ylabel="Expectation value",
               title="Single-Site Pauli (all vanish)")
        ax.legend(); ax.set_ylim(-0.1, 0.1)
        for f in ("png","pdf"):
            fig.savefig(FIG/f"pauli.{f}", dpi=300, bbox_inches="tight")
        plt.close(fig)

    # 6. String-order decay curve
    if strings:
        # Group by separation d = j - i
        from collections import defaultdict
        by_sep = defaultdict(list)
        for s in strings:
            by_sep[s["j"] - s["i"]].append(abs(s["v"]))
        seps = sorted(by_sep.keys())
        means = [np.mean(by_sep[s]) for s in seps]
        fig, ax = plt.subplots(figsize=(5.5, 3.5))
        ax.plot(seps, means, "o-", color=COLORS[2], ms=7, lw=2)
        ax.set(xlabel="Separation $|j-i|$", ylabel=r"$|\langle Z_i \prod X_k Z_j \rangle|$",
               title="String-Order Correlator vs Separation")
        for f in ("png","pdf"):
            fig.savefig(FIG/f"string_order_decay.{f}", dpi=300, bbox_inches="tight")
        plt.close(fig)

    print("  Figures saved.")


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════
def main():
    t0 = time.perf_counter()

    # Load existing results
    rpath = DATA / "results.json"
    if rpath.exists():
        old = json.loads(rpath.read_text(encoding="utf-8"))
        print(f"  Loaded existing results: {list(old.keys())}")
    else:
        old = {}

    Ut = make_target("cluster", n_match=1)

    # GRAPE sweep
    grape_data = grape_sweep(Ut)

    # Holographic analysis
    holo_data = holographic_analysis(Ut)

    # Save combined results
    combined = dict(old)
    combined["grape_sweep"] = grape_data
    combined["holographic"] = holo_data
    combined["snap_cavity_only"] = True
    combined["analytical_decomposition"] = {
        "target": "SWAP * CZ * (H x I)",
        "CZ_equivalent": "SNAP_QC(0, pi) = diag(1,1,1,-1) in 4D subspace",
        "H_equivalent": "R_q(pi/2, 0)",
        "SWAP_decomposition": "3 CNOT = 6 Hadamard + 3 CZ",
        "min_selective_layers": 4,
        "framework_limitation": "cqed_sim ideal SNAP is I_q x diag(phases) "
            "(cavity-only), not |g><g| x I + |e><e| x diag(phases) "
            "(qubit-conditional). SQR IS qubit-conditional."
    }
    _jdump(DATA / "results_combined.json", combined)

    # Figures
    make_figures(grape_data, holo_data, old)

    dt = time.perf_counter() - t0
    print(f"\n{'='*60}\n  COMPLETE in {dt:.0f}s\n{'='*60}")
    # Summary
    for e in grape_data:
        ok = "PASS" if e["fid"] >= 0.999 else ""
        print(f"  GRAPE {e['dns']:4d}ns: F={e['fid']:.6f} {ok}")
    print(f"  Correlation length: {holo_data.get('correlation_length','?')} sites")


if __name__ == "__main__":
    main()
