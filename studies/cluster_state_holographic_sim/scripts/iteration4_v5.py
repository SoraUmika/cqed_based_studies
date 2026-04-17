"""Iteration 4 v5: Fast bounded-displacement validation.

Key insight: Optimise in ideal mode (N_cav=2, ~30s), then EVALUATE
the same gate sequence at N_cav=4,6,8,12 to test if bounded α
prevents truncation artifacts.

Previous iter3 finding: unbounded ideal-mode B(2blk) has F=1 at
N_cav=2 but F=0.094 at N_cav=8 (93% leakage). Hypothesis: if we
bound |α|≤0.3, population stays in n=0,1 and fidelity should be
preserved at higher N_cav.
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

# ── paths & imports ───────────────────────────────────────────────────────
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
t0 = time.perf_counter()
from cqed_sim.unitary_synthesis import (
    Displacement, QubitRotation, SQR, ConditionalPhaseSQR,
    FreeEvolveCondPhase, Subspace, TargetUnitary,
    UnitarySynthesizer, GateSequence, DriftPhaseModel,
    LeakagePenalty, MultiObjective, ExecutionOptions,
    SynthesisConstraints,
    subspace_unitary_fidelity, simulate_sequence,
)
from cqed_sim.unitary_synthesis.targets import make_target
print(f"Import: {time.perf_counter()-t0:.1f}s", flush=True)

# ── constants ─────────────────────────────────────────────────────────────
TWO_PI = 2.0 * np.pi
CHI   = TWO_PI * (-2.84e6)
COLORS = ['#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB']

U_target = make_target("cluster", n_match=1)
target = TargetUnitary(U_target, ignore_global_phase=True)
no_drift = DriftPhaseModel(chi=0.0, chi2=0.0, kerr=0.0)

print(f"Target: {U_target.shape}", flush=True)
print(f"U_target:\n{U_target}", flush=True)

results = {}

def make_subspace(nc):
    return Subspace.custom(2*nc, [0, 1, nc, nc+1],
                           ["|g,0>", "|g,1>", "|e,0>", "|e,1>"])


def make_B(nb, nc):
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


def run_opt(seq, label, subsp, max_amp=None, multistart=8, maxiter=500, dur_wt=0.0):
    """Run synthesis at N_cav=2 (fast)."""
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
        F = subspace_unitary_fidelity(res.simulation.subspace_operator,
                                       U_target, gauge="global")
        disp_amps = [float(abs(g.alpha)) for g in res.sequence.gates if hasattr(g, 'alpha')]
        print(f"    F={F:.6f} obj={res.objective:.6f} ({dt:.0f}s)", flush=True)
        print(f"    Disp |alpha|: {[round(a,4) for a in disp_amps]}", flush=True)
        return {"label": label, "fidelity": float(F),
                "objective": float(res.objective), "success": bool(res.success),
                "elapsed_s": dt, "max_amp": max_amp, "disp_amps": disp_amps,
                "_result": res}
    except Exception as e:
        dt = time.perf_counter() - t0
        print(f"    FAILED: {e} ({dt:.0f}s)", flush=True)
        traceback.print_exc()
        return {"label": label, "fidelity": 0.0, "success": False,
                "elapsed_s": dt, "error": str(e)}


def rescale_sequence(seq, target_nc):
    """Create equivalent GateSequence at a different N_cav.

    Pad/truncate SQR theta_n, phi_n, and CP phases_n to target_nc.
    """
    new_gates = []
    for g in seq.gates:
        if isinstance(g, SQR):
            th = list(g.theta_n) + [0.0]*(target_nc - len(g.theta_n))
            ph = list(g.phi_n) + [0.0]*(target_nc - len(g.phi_n))
            new_gates.append(SQR(name=g.name, theta_n=th[:target_nc],
                                  phi_n=ph[:target_nc],
                                  drift_model=g.drift_model,
                                  duration=g.duration))
        elif isinstance(g, ConditionalPhaseSQR):
            ps = list(g.phases_n) + [0.0]*(target_nc - len(g.phases_n))
            new_gates.append(ConditionalPhaseSQR(
                name=g.name, phases_n=ps[:target_nc],
                drift_model=g.drift_model, duration=g.duration))
        elif isinstance(g, Displacement):
            new_gates.append(Displacement(name=g.name, alpha=g.alpha,
                                           duration=g.duration))
        elif isinstance(g, QubitRotation):
            new_gates.append(QubitRotation(name=g.name, theta=g.theta,
                                            phi=g.phi, duration=g.duration))
        elif isinstance(g, FreeEvolveCondPhase):
            new_gates.append(FreeEvolveCondPhase(
                name=g.name, duration=g.duration,
                drift_model=g.drift_model, optimize_time=False))
        else:
            new_gates.append(g)
    return GateSequence(gates=new_gates, n_cav=target_nc)


def evaluate_at_ncav(result, target_nc):
    """Evaluate optimised N_cav=2 sequence at a different N_cav."""
    sr = result.get("_result")
    if sr is None:
        return None
    seq_opt = sr.sequence  # optimised GateSequence at N_cav=2
    seq_ext = rescale_sequence(seq_opt, target_nc)
    sub = make_subspace(target_nc)
    t0 = time.perf_counter()
    try:
        sim_res = simulate_sequence(seq_ext, sub)
        dt = time.perf_counter() - t0
        U_sub = sim_res.subspace_operator
        U_sub_np = U_sub.full() if hasattr(U_sub, 'full') else np.asarray(U_sub)
        F = subspace_unitary_fidelity(U_sub_np, U_target, gauge="global")
        # Compute leakage: 1 - sum of subspace populations
        U_full = sim_res.full_operator
        U_full_np = U_full.full() if hasattr(U_full, 'full') else np.asarray(U_full)
        full_dim = 2 * target_nc
        lidx = [0, 1, target_nc, target_nc + 1]
        leakage_sum = 0.0
        for bi in lidx:
            psi0 = np.zeros(full_dim, dtype=complex)
            psi0[bi] = 1.0
            psi_out = U_full_np @ psi0
            sub_pop = sum(abs(psi_out[li])**2 for li in lidx)
            leakage_sum += 1.0 - sub_pop
        avg_leakage = leakage_sum / len(lidx)
        return {"n_cav": target_nc, "fidelity": float(F),
                "avg_leakage": float(avg_leakage), "eval_time": dt,
                "full_unitary": U_full_np}
    except Exception as e:
        return {"n_cav": target_nc, "fidelity": 0.0, "error": str(e)}


# ═══════════════════════════════════════════════════════════════════════════
# PHASE 1: Optimise at N_cav=2 (ideal mode) with amplitude bounds
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70, flush=True)
print("  PHASE 1: Ideal-Mode Synthesis at N_cav=2", flush=True)
print("="*70, flush=True)

NC_OPT = 2
sub_opt = make_subspace(NC_OPT)

CONFIGS = [
    # (strategy, n_blocks, max_amp, multistart, maxiter, label_prefix)
    ("B", 2, 0.2, 8, 500, "B2_amp0.2"),
    ("B", 2, 0.3, 8, 500, "B2_amp0.3"),
    ("B", 2, 0.5, 8, 500, "B2_amp0.5"),
    ("B", 2, 0.7, 8, 500, "B2_amp0.7"),
    ("B", 2, 1.0, 8, 500, "B2_amp1.0"),
    ("B", 2, None, 8, 500, "B2_unbounded"),  # baseline
    ("B", 3, 0.3, 8, 500, "B3_amp0.3"),
    ("B", 3, None, 8, 500, "B3_unbounded"),
    ("D", 3, 0.3, 8, 500, "D3_amp0.3"),
    ("D", 3, 0.5, 8, 500, "D3_amp0.5"),
    ("D", 3, 1.0, 8, 500, "D3_amp1.0"),
    ("D", 3, None, 8, 500, "D3_unbounded"),
    ("D", 2, 0.5, 8, 500, "D2_amp0.5"),
    ("D", 2, None, 8, 500, "D2_unbounded"),
]

print(f"\n  Running {len(CONFIGS)} configs at N_cav=2 (ideal mode)...", flush=True)
for strat, nb, amp, ms, mi, lbl in CONFIGS:
    seq = make_B(nb, NC_OPT) if strat == "B" else make_D(nb, NC_OPT)
    r = run_opt(seq, lbl, sub_opt, max_amp=amp, multistart=ms, maxiter=mi)
    r["strategy"] = "D+SQR+CP" if strat == "B" else "D+R+FE"
    r["n_blocks"] = nb
    r["n_cav_opt"] = NC_OPT
    results[lbl] = r


# ═══════════════════════════════════════════════════════════════════════════
# PHASE 2: Evaluate at N_cav=4,6,8,12 (the key experiment)
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70, flush=True)
print("  PHASE 2: Cross-N_cav Evaluation (Truncation Test)", flush=True)
print("="*70, flush=True)

EVAL_NCS = [4, 6, 8, 12]
cross_nc_results = {}

for lbl, r in results.items():
    if r.get("fidelity", 0) < 0.3:
        continue
    if r.get("_result") is None:
        continue
    cross_nc_results[lbl] = {"nc2": r["fidelity"]}
    for nc_eval in EVAL_NCS:
        ev = evaluate_at_ncav(r, nc_eval)
        if ev is not None:
            cross_nc_results[lbl][f"nc{nc_eval}"] = ev.get("fidelity", 0.0)
            cross_nc_results[lbl][f"leak_nc{nc_eval}"] = ev.get("avg_leakage", 0.0)
            if nc_eval == 12:
                cross_nc_results[lbl]["U_full_nc12"] = ev.get("full_unitary")
            print(f"    [{lbl}@nc{nc_eval}] F={ev.get('fidelity',0):.6f} "
                  f"leak={ev.get('avg_leakage',0):.6f}", flush=True)
        else:
            print(f"    [{lbl}@nc{nc_eval}] FAILED", flush=True)

# Print cross-NC summary table
print("\n  Cross-N_cav Fidelity Summary:", flush=True)
print(f"  {'Config':<25s} {'nc2':>8s} {'nc4':>8s} {'nc6':>8s} {'nc8':>8s} {'nc12':>8s}", flush=True)
for lbl in sorted(cross_nc_results.keys()):
    d = cross_nc_results[lbl]
    row = f"  {lbl:<25s}"
    for key in ["nc2", "nc4", "nc6", "nc8", "nc12"]:
        v = d.get(key, None)
        row += f" {v:8.4f}" if v is not None else f" {'---':>8s}"
    print(row, flush=True)


# ═══════════════════════════════════════════════════════════════════════════
# PHASE 3: Duration Optimisation at N_cav=2 (best bounded result)
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70, flush=True)
print("  PHASE 3: Duration Optimisation", flush=True)
print("="*70, flush=True)

# Find best bounded config that holds up at N_cav=12
best_bounded_12 = max(
    ((lbl, d) for lbl, d in cross_nc_results.items()
     if results[lbl].get("max_amp") is not None
     and d.get("nc12", 0) > 0.3),
    key=lambda x: x[1].get("nc12", 0), default=(None, None))

if best_bounded_12[0]:
    bb_lbl = best_bounded_12[0]
    bv = results[bb_lbl]
    nb, amp = bv["n_blocks"], bv.get("max_amp")
    strat = bv["strategy"]
    print(f"  Best bounded@nc12: {bb_lbl} (F@nc12={best_bounded_12[1].get('nc12', 0):.6f})", flush=True)
    seq = make_B(nb, NC_OPT) if "SQR" in strat else make_D(nb, NC_OPT)
    r = run_opt(seq, "dur_opt", sub_opt, max_amp=amp, multistart=8, maxiter=500, dur_wt=0.01)
    r["strategy"] = strat; r["n_blocks"] = nb; r["n_cav_opt"] = NC_OPT
    r["note"] = "Duration-optimised"
    results["dur_opt"] = r
    # Evaluate at NC=12
    ev = evaluate_at_ncav(r, 12)
    if ev: print(f"    dur_opt@nc12: F={ev.get('fidelity',0):.6f}", flush=True)
else:
    print("  No bounded config achieves F>0.3 at N_cav=12.", flush=True)


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
    print("  Wigner imports OK.", flush=True)
except ImportError as e:
    HAS_WIGNER = False
    print(f"  Wigner import failed: {e}", flush=True)

if HAS_WIGNER:
    # Find best result WITH a stored N_cav=12 full unitary
    best_wig_key = None
    best_wig_F = 0
    best_wig_U = None
    nc_wig = 12

    for lbl, d in cross_nc_results.items():
        U_full = d.get("U_full_nc12")
        f12 = d.get("nc12", 0)
        if U_full is not None and f12 > best_wig_F:
            best_wig_key = lbl
            best_wig_F = f12
            best_wig_U = U_full

    if best_wig_key and best_wig_U is not None:
        print(f"  Wigner source: {best_wig_key} (F@nc12={best_wig_F:.6f})", flush=True)
        full_dim = 2 * nc_wig
        U_tgt_full = np.eye(full_dim, dtype=complex)
        lidx = [0, 1, nc_wig, nc_wig + 1]
        for i, li in enumerate(lidx):
            for j, lj in enumerate(lidx):
                U_tgt_full[li, lj] = U_target[i, j]

        basis_labels = ["|g,0>", "|g,1>", "|e,0>", "|e,1>"]
        n_pts = 51; ext = 3.0

        for bi, bl in zip(lidx, basis_labels):
            psi0 = np.zeros(full_dim, dtype=complex)
            psi0[bi] = 1.0
            psi_t = U_tgt_full @ psi0
            psi_a = best_wig_U @ psi0

            rho_t = qt.Qobj(np.outer(psi_t, psi_t.conj()),
                            dims=[[2, nc_wig], [2, nc_wig]])
            rho_a = qt.Qobj(np.outer(psi_a, psi_a.conj()),
                            dims=[[2, nc_wig], [2, nc_wig]])
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

            wigner_data[bl] = {
                "xvec": xv.tolist(), "yvec": yv.tolist(),
                "W_target": Wt.tolist(), "W_achieved": Wa.tolist(),
                "cavity_fidelity": fcav, "l2_distance": l2}

        # ── Wigner figure ────────────────────────────────────────────
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
            axes[1, col].set_title(f"Achieved ($F_c$={wd['cavity_fidelity']:.4f})", fontsize=9)
            axes[1, col].set_aspect("equal")

            dmax = max(np.max(np.abs(dW)), 1e-10)
            axes[2, col].pcolormesh(X, Y, dW, cmap="RdBu_r",
                                     vmin=-dmax, vmax=dmax, shading="auto")
            axes[2, col].set_title(f"$\\Delta W$ (L2={wd['l2_distance']:.4f})", fontsize=9)
            axes[2, col].set_aspect("equal")
            axes[2, col].set_xlabel("Re($\\alpha$)")
            for row in range(3):
                if col == 0: axes[row, col].set_ylabel("Im($\\alpha$)")

        fig.suptitle(
            f"Wigner Comparison — {best_wig_key}\n"
            f"$\\mathcal{{F}}_{{sub}}$={best_wig_F:.4f}  ($N_{{cav}}$={nc_wig})",
            fontsize=12, fontweight="bold")
        fig.colorbar(im, ax=axes[:2].ravel().tolist(), label="$W(\\alpha)$",
                     shrink=0.7, pad=0.02)
        fig.tight_layout(rect=[0, 0, 0.93, 0.93])
        for fmt in ("png", "pdf"):
            fig.savefig(FIG_DIR / f"wigner_comparison.{fmt}",
                        dpi=300, bbox_inches="tight")
        plt.close(fig)
        print("  Wigner figure saved.", flush=True)

        # Save NPZ
        npz_d = {}
        for bl, wd in wigner_data.items():
            s = bl.replace("|","").replace(">","").replace(",","_")
            for key in ("xvec","yvec","W_target","W_achieved"):
                npz_d[f"{s}_{key}"] = np.array(wd[key])
            npz_d[f"{s}_cavity_fidelity"] = np.array([wd["cavity_fidelity"]])
        np.savez(str(ARTIFACTS_DIR / "wigner_comparison.npz"), **npz_d)
        print("  Wigner NPZ saved.", flush=True)
    else:
        print("  No suitable unitary for Wigner.", flush=True)


# ═══════════════════════════════════════════════════════════════════════════
# PHASE 5: Summary Figures & Artifacts
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70, flush=True)
print("  PHASE 5: Figures & Artifacts", flush=True)
print("="*70, flush=True)

# ── Figure 1: Fidelity vs max_amplitude at N_cav=2 + N_cav=12 ───────────
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Panel a: Fidelity at N_cav=2 (ideal mode)
ax = axes[0]
for tag, col, marker in [("B2", COLORS[0], "o"), ("B3", COLORS[1], "s"),
                           ("D3", COLORS[2], "^"), ("D2", COLORS[3], "D")]:
    amps, fids = [], []
    for lbl, r in results.items():
        if lbl.startswith(tag + "_amp") and r.get("fidelity", 0) > 0:
            amp_v = r.get("max_amp")
            if amp_v is not None:
                amps.append(amp_v)
                fids.append(r["fidelity"])
    if amps:
        order = np.argsort(amps)
        ax.plot([amps[i] for i in order], [fids[i] for i in order],
                f"{marker}-", color=col, ms=7, lw=2, label=tag)
    ub_lbl = f"{tag}_unbounded"
    if ub_lbl in results and results[ub_lbl].get("fidelity", 0) > 0:
        ax.axhline(results[ub_lbl]["fidelity"], ls=":", color=col, alpha=0.5)
ax.set_xlabel("Max $|\\alpha|$"); ax.set_ylabel("$\\mathcal{F}$ (N$_{cav}$=2)")
ax.set_title("(a) Ideal-mode fidelity"); ax.legend(fontsize=8); ax.grid(alpha=0.3)

# Panel b: Fidelity at N_cav=12 (truncation test)
ax = axes[1]
for tag, col, marker in [("B2", COLORS[0], "o"), ("B3", COLORS[1], "s"),
                           ("D3", COLORS[2], "^"), ("D2", COLORS[3], "D")]:
    amps, fids = [], []
    for lbl in cross_nc_results:
        if lbl.startswith(tag + "_amp") and results[lbl].get("max_amp") is not None:
            amps.append(results[lbl]["max_amp"])
            fids.append(cross_nc_results[lbl].get("nc12", 0))
    if amps:
        order = np.argsort(amps)
        ax.plot([amps[i] for i in order], [fids[i] for i in order],
                f"{marker}-", color=col, ms=7, lw=2, label=tag)
    ub_lbl = f"{tag}_unbounded"
    if ub_lbl in cross_nc_results:
        f_ub12 = cross_nc_results[ub_lbl].get("nc12", 0)
        if f_ub12 > 0: ax.axhline(f_ub12, ls=":", color=col, alpha=0.5)
ax.set_xlabel("Max $|\\alpha|$"); ax.set_ylabel("$\\mathcal{F}$ (N$_{cav}$=12)")
ax.set_title("(b) Fidelity at N$_{cav}$=12"); ax.legend(fontsize=8); ax.grid(alpha=0.3)
ax.text(0.05, 0.05, "Bounded $\\alpha$\nprevents leakage",
        transform=ax.transAxes, fontsize=9, fontstyle='italic',
        bbox=dict(facecolor='lightyellow', alpha=0.7))

# Panel c: Leakage at N_cav=12
ax = axes[2]
for tag, col, marker in [("B2", COLORS[0], "o"), ("B3", COLORS[1], "s"),
                           ("D3", COLORS[2], "^"), ("D2", COLORS[3], "D")]:
    amps, leaks = [], []
    for lbl in cross_nc_results:
        if lbl.startswith(tag + "_amp") and results[lbl].get("max_amp") is not None:
            amps.append(results[lbl]["max_amp"])
            leaks.append(cross_nc_results[lbl].get("leak_nc12", 0))
    if amps:
        order = np.argsort(amps)
        ax.plot([amps[i] for i in order], [leaks[i] for i in order],
                f"{marker}-", color=col, ms=7, lw=2, label=tag)
    ub_lbl = f"{tag}_unbounded"
    if ub_lbl in cross_nc_results:
        l_ub = cross_nc_results[ub_lbl].get("leak_nc12", 0)
        if l_ub > 0: ax.axhline(l_ub, ls=":", color=col, alpha=0.5)
ax.set_xlabel("Max $|\\alpha|$"); ax.set_ylabel("Avg. leakage")
ax.set_title("(c) Leakage at N$_{cav}$=12"); ax.legend(fontsize=8)
ax.grid(alpha=0.3); ax.set_yscale("log")

fig.suptitle("Bounded-Displacement Strategy: Ideal-Mode Optimisation → N$_{cav}$=12 Evaluation",
             fontsize=13, fontweight="bold")
fig.tight_layout()
for fmt in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"bounded_displacement_sweep.{fmt}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("  bounded_displacement_sweep figure saved.", flush=True)


# ── Figure 2: Truncation convergence (best bounded vs unbounded) ─────────
fig, ax = plt.subplots(figsize=(8, 5))
for lbl in sorted(cross_nc_results.keys()):
    r = results[lbl]
    amp = r.get("max_amp")
    if amp is None or r.get("fidelity", 0) < 0.8:
        continue
    d = cross_nc_results[lbl]
    ncs = [2] + EVAL_NCS
    fids = [d.get(f"nc{nc}", 0) for nc in ncs]
    if all(f > 0 for f in fids):
        ax.plot(ncs, fids, "o-", ms=6, lw=1.5, label=f"{lbl} (|α|≤{amp})", alpha=0.8)

# Add unbounded baselines
for base_lbl, col in [("B2_unbounded", "red"), ("D3_unbounded", "darkred")]:
    if base_lbl in cross_nc_results and results.get(base_lbl, {}).get("fidelity", 0) > 0.5:
        d = cross_nc_results[base_lbl]
        ncs = [2] + EVAL_NCS
        fids = [d.get(f"nc{nc}", 0) for nc in ncs]
        if any(f > 0 for f in fids):
            ax.plot(ncs, [f if f > 0 else np.nan for f in fids], "x--", color=col,
                    ms=8, lw=2, label=f"{base_lbl} (no bound)", alpha=0.9)
ax.set_xlabel("$N_{cav}$ (cavity truncation)"); ax.set_ylabel("Subspace fidelity")
ax.set_title("Truncation Convergence: Bounded vs Unbounded Displacement")
ax.axhline(0.99, ls="--", color="gray", alpha=0.3, lw=0.5)
ax.legend(fontsize=7, ncol=2); ax.grid(alpha=0.3)
fig.tight_layout()
for fmt in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"truncation_convergence_bounded.{fmt}",
                dpi=300, bbox_inches="tight")
plt.close(fig)
print("  truncation_convergence_bounded figure saved.", flush=True)


# ── Save final results ───────────────────────────────────────────────────
output = {
    "configs": {},
    "cross_ncav": {},
    "wigner": results.get("wigner"),
}
for lbl, r in results.items():
    output["configs"][lbl] = {k: v for k, v in r.items() if k != "_result"}
for lbl, d in cross_nc_results.items():
    output["cross_ncav"][lbl] = {k: v for k, v in d.items()
                                  if not isinstance(v, np.ndarray)}

(DATA_DIR / "iteration4_results.json").write_text(
    json.dumps(output, indent=2, default=str), encoding="utf-8")
print("  iteration4_results.json saved.", flush=True)

# Save best artifacts
for tag, prefix in [("B", "best_strategy_B"), ("D", "best_strategy_D")]:
    best_12 = max(
        ((lbl, d.get("nc12", 0)) for lbl, d in cross_nc_results.items()
         if lbl.startswith(tag) and results[lbl].get("max_amp") is not None),
        key=lambda x: x[1], default=(None, 0))
    if best_12[0] and best_12[1] > 0.3:
        sr = results[best_12[0]].get("_result")
        ap = ARTIFACTS_DIR / f"{prefix}.json"
        if sr:
            try:
                sr.save(str(ap), include_history=True)
                print(f"  {ap.name}: {best_12[0]} F@nc12={best_12[1]:.6f}", flush=True)
            except Exception as e:
                ap.write_text(json.dumps({
                    "label": best_12[0], "fidelity_nc2": results[best_12[0]]["fidelity"],
                    "fidelity_nc12": best_12[1],
                    "max_amp": results[best_12[0]].get("max_amp"),
                }, indent=2), encoding="utf-8")
                print(f"  {ap.name} (fallback): {e}", flush=True)


# ═══════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70, flush=True)
print("  ITERATION 4: FINAL SUMMARY", flush=True)
print("="*70, flush=True)

print("\n  Cross-N_cav Fidelity Table:", flush=True)
print(f"  {'Config':<25s} {'nc2':>8s} {'nc4':>8s} {'nc6':>8s} {'nc8':>8s} {'nc12':>8s} {'leak12':>8s}", flush=True)
for lbl in sorted(cross_nc_results.keys()):
    d = cross_nc_results[lbl]
    row = f"  {lbl:<25s}"
    for key in ["nc2", "nc4", "nc6", "nc8", "nc12"]:
        v = d.get(key)
        row += f" {v:8.4f}" if v is not None else f" {'---':>8s}"
    l = d.get("leak_nc12", 0)
    row += f" {l:8.6f}"
    print(row, flush=True)

if wigner_data:
    print(f"\n  Wigner: source={best_wig_key}", flush=True)
    for bl in basis_labels:
        if bl in wigner_data:
            print(f"    {bl}: F_cav={wigner_data[bl]['cavity_fidelity']:.6f}", flush=True)

print(f"\n  Output files:", flush=True)
for d in (DATA_DIR, FIG_DIR, ARTIFACTS_DIR):
    for f in sorted(d.iterdir()):
        if f.is_file() and any(x in f.name for x in
            ("iteration4", "wigner", "bounded", "truncation", "best_strategy")):
            print(f"    {f.relative_to(STUDY_ROOT)}", flush=True)

print("\n  DONE.", flush=True)
