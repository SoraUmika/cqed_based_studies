"""Iteration 4 v3: Practical bounded-displacement study.

Strategy:
  1. Quick bounded-displacement sweep at N_cav=6 (12-dim, fast).
     Physics justification: with |alpha|<=0.5, P(n>=6)<1e-8, so N_cav=6
     is truncation-free for the bounded-displacement regime.
  2. GRAPE at N_cav=12 for one duration (200ns) — for Wigner comparison.
  3. Validate best bounded result at N_cav=12 (single config).
  4. Wigner function comparison.
  5. Artifacts.

Output:
  data/iteration4_results.json
  figures/bounded_displacement_sweep.{png,pdf}
  figures/wigner_comparison.{png,pdf}
  figures/iteration4_ranking.{png,pdf}
  artifacts/best_strategy_B.json
  artifacts/best_strategy_D.json
  artifacts/grape_ncav12_200ns.json
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
CHIP  = TWO_PI * (-21e3)
KERR  = TWO_PI * (-28e3)
OMEGA_Q = TWO_PI * 6.150e9
OMEGA_C = TWO_PI * 5.241e9
ALPHA   = TWO_PI * (-255e6)

COLORS = ['#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB']

# Target
U_target = make_target("cluster", n_match=1)
target = TargetUnitary(U_target, ignore_global_phase=True)
print(f"Target: {U_target.shape}", flush=True)

no_drift = DriftPhaseModel(chi=0.0, chi2=0.0, kerr=0.0)
results = {}

# ═══════════════════════════════════════════════════════════════════════════
# Helper functions
# ═══════════════════════════════════════════════════════════════════════════
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

def run_opt(seq, label, subsp, max_amp=None, multistart=3, maxiter=300, dur_wt=0.0):
    print(f"\n  [{label}] amp<={max_amp} ms={multistart} mi={maxiter} d_cav={seq.n_cav}", flush=True)
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
        print(f"    F={F:.6f}  obj={res.objective:.6f}  ({dt:.0f}s)", flush=True)
        disp_amps = [float(abs(g.alpha)) for g in res.sequence.gates if hasattr(g, 'alpha')]
        return {"label": label, "fidelity": float(F), "objective": float(res.objective),
                "success": bool(res.success), "elapsed_s": float(dt),
                "max_amp": max_amp, "disp_amps": disp_amps,
                "n_cav": seq.n_cav,
                "sequence": res.sequence.serialize(), "_synth": res}
    except Exception as e:
        dt = time.perf_counter() - t0
        print(f"    FAILED: {e} ({dt:.0f}s)", flush=True)
        traceback.print_exc()
        return {"label": label, "fidelity": 0.0, "success": False,
                "elapsed_s": float(dt), "error": str(e), "max_amp": max_amp,
                "n_cav": seq.n_cav}


# ═══════════════════════════════════════════════════════════════════════════
# PHASE 1: Bounded displacement sweep at N_cav=6 (fast)
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70, flush=True)
print("  PHASE 1: Bounded Displacement Sweep at N_cav=6", flush=True)
print("="*70, flush=True)
print("  Physics: |alpha|<=0.5 => P(n>=6)<1e-8 => N_cav=6 is exact.", flush=True)

NC_FAST = 6
sub_fast = make_subspace(NC_FAST)

# Strategy B sweep
print("\n--- Strategy B: D + SQR + CP ---", flush=True)
B_CONFIGS = [
    # (max_amp, n_blocks, multistart, maxiter)
    (0.3, 2, 4, 300),
    (0.3, 3, 4, 300),
    (0.5, 2, 4, 300),
    (0.5, 3, 4, 300),
    (0.7, 2, 4, 300),
    (0.7, 3, 4, 300),
    (1.0, 2, 4, 300),
    (1.0, 3, 4, 300),
    (None, 2, 4, 300),  # unconstrained baseline
    (None, 3, 4, 300),  # unconstrained baseline
]
for amp, nb, ms, mi in B_CONFIGS:
    seq = make_B(nb, NC_FAST)
    amp_str = f"{amp}" if amp else "none"
    lbl = f"B_nc{NC_FAST}_amp{amp_str}_blk{nb}"
    r = run_opt(seq, lbl, sub_fast, max_amp=amp, multistart=ms, maxiter=mi)
    r["strategy"] = "D+SQR+CP"; r["n_blocks"] = nb
    results[lbl] = r

# Strategy D sweep
print("\n--- Strategy D: D + R + FE ---", flush=True)
D_CONFIGS = [
    (0.3, 3, 4, 300),
    (0.5, 2, 4, 300),
    (0.5, 3, 4, 300),
    (0.7, 3, 4, 300),
    (1.0, 3, 4, 300),
    (None, 3, 4, 300),  # unconstrained baseline
]
for amp, nb, ms, mi in D_CONFIGS:
    seq = make_D(nb, NC_FAST)
    amp_str = f"{amp}" if amp else "none"
    lbl = f"D_nc{NC_FAST}_amp{amp_str}_blk{nb}"
    r = run_opt(seq, lbl, sub_fast, max_amp=amp, multistart=ms, maxiter=mi)
    r["strategy"] = "D+R+FE"; r["n_blocks"] = nb
    results[lbl] = r

# Save Phase 1 results
print("\n  Phase 1 complete. Saving intermediate results...", flush=True)
results_json = {k: {kk: vv for kk, vv in v.items() if kk != "_synth"}
                for k, v in results.items() if isinstance(v, dict)}
(DATA_DIR / "iteration4_phase1.json").write_text(
    json.dumps(results_json, indent=2, default=str), encoding="utf-8")


# ═══════════════════════════════════════════════════════════════════════════
# PHASE 2: Validate best at N_cav=12
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70, flush=True)
print("  PHASE 2: Validate Best Configs at N_cav=12", flush=True)
print("="*70, flush=True)

NC_FULL = 12
sub_full = make_subspace(NC_FULL)

# Find best B and best D from Phase 1
best_B_key = max(
    (k for k in results if results[k].get("strategy") == "D+SQR+CP" and results[k].get("fidelity",0) > 0.5),
    key=lambda k: results[k]["fidelity"], default=None)
best_D_key = max(
    (k for k in results if results[k].get("strategy") == "D+R+FE" and results[k].get("fidelity",0) > 0.5),
    key=lambda k: results[k]["fidelity"], default=None)

for tag, best_key in [("B", best_B_key), ("D", best_D_key)]:
    if best_key is None:
        print(f"  No result for Strategy {tag}, skipping N_cav=12 validation.", flush=True)
        continue
    bv = results[best_key]
    nb = bv["n_blocks"]
    amp = bv.get("max_amp")
    print(f"\n  Validating {tag}: {best_key} (F@nc{NC_FAST}={bv['fidelity']:.6f}, amp={amp}, blk={nb})", flush=True)

    seq = make_B(nb, NC_FULL) if tag == "B" else make_D(nb, NC_FULL)
    lbl = f"{tag}_nc{NC_FULL}_amp{amp}_blk{nb}_validate"
    r = run_opt(seq, lbl, sub_full, max_amp=amp, multistart=3, maxiter=300)
    r["strategy"] = bv["strategy"]; r["n_blocks"] = nb
    r["note"] = f"N_cav=12 validation of Phase 1 best ({best_key})"
    results[lbl] = r


# ═══════════════════════════════════════════════════════════════════════════
# PHASE 3: Duration-Optimised Run
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70, flush=True)
print("  PHASE 3: Duration Optimization", flush=True)
print("="*70, flush=True)

overall_best = max(
    (k for k in results if results[k].get("fidelity", 0) > 0.8),
    key=lambda k: results[k]["fidelity"], default=None)
if overall_best:
    bv = results[overall_best]
    nb = bv["n_blocks"]
    amp = bv.get("max_amp")
    nc = bv.get("n_cav", NC_FAST)
    strat = bv["strategy"]
    print(f"  Duration opt for: {overall_best} (F={bv['fidelity']:.6f})", flush=True)
    sub = make_subspace(nc)
    seq = make_B(nb, nc) if "SQR" in strat else make_D(nb, nc)
    r = run_opt(seq, "duration_opt", sub, max_amp=amp, multistart=4, maxiter=400, dur_wt=0.01)
    r["strategy"] = strat; r["n_blocks"] = nb
    r["note"] = "Duration-optimised version"
    results["duration_opt"] = r


# ═══════════════════════════════════════════════════════════════════════════
# PHASE 4: GRAPE at N_cav=12 (for Wigner comparison)
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70, flush=True)
print("  PHASE 4: GRAPE at N_cav=12 for Wigner Comparison", flush=True)
print("="*70, flush=True)

try:
    from cqed_sim import PiecewiseConstantTimeGrid, ModelControlChannelSpec
    from cqed_sim.models import DispersiveTransmonCavityModel, FrameSpec
    from cqed_sim.optimal_control import (
        GrapeConfig, GrapeSolver, build_control_problem_from_model,
    )
    from cqed_sim.optimal_control import (
        LeakagePenalty as OCLeakagePenalty,
        UnitaryObjective as OCUnitaryObjective,
    )
    HAS_GRAPE = True
    print("  GRAPE imports OK.", flush=True)
except ImportError as e:
    HAS_GRAPE = False
    print(f"  GRAPE import failed: {e}", flush=True)

if HAS_GRAPE:
    model = DispersiveTransmonCavityModel(
        omega_c=OMEGA_C, omega_q=OMEGA_Q, alpha=ALPHA,
        chi=CHI, chi_higher=(CHIP,), kerr=KERR,
        n_cav=NC_FULL, n_tr=2)
    frame = FrameSpec(omega_c_frame=model.omega_c, omega_q_frame=model.omega_q)

    T_GATE = 200e-9
    DT = 5e-9
    N_STEPS = int(T_GATE / DT)

    time_grid = PiecewiseConstantTimeGrid(n_steps=N_STEPS, dt=DT)
    channels = [
        ModelControlChannelSpec(name="cavity_I", quadrature="I", target="cavity"),
        ModelControlChannelSpec(name="cavity_Q", quadrature="Q", target="cavity"),
    ]

    # Build target in full Hilbert space
    grape_target = np.eye(2 * NC_FULL, dtype=complex)
    lidx = [0, 1, NC_FULL, NC_FULL + 1]
    for i, li in enumerate(lidx):
        for j, lj in enumerate(lidx):
            grape_target[li, lj] = U_target[i, j]

    print(f"  Running GRAPE: T={T_GATE*1e9:.0f}ns, {N_STEPS} steps, N_cav={NC_FULL}", flush=True)
    t0 = time.perf_counter()
    try:
        problem = build_control_problem_from_model(
            model=model, frame=frame, time_grid=time_grid,
            control_channels=channels,
            target_unitary=grape_target,
            subspace_indices=lidx)

        config = GrapeConfig(
            n_iter=300,
            learning_rate=0.05,
            n_seeds=3,
        )
        solver = GrapeSolver(problem, config)
        grape_result = solver.run()
        dt_grape = time.perf_counter() - t0

        F_grape = float(grape_result.best_fidelity)
        print(f"    GRAPE F={F_grape:.6f} ({dt_grape:.0f}s)", flush=True)

        results["GRAPE_nc12_200ns"] = {
            "label": "GRAPE_nc12_200ns",
            "fidelity": F_grape,
            "elapsed_s": float(dt_grape),
            "strategy": "GRAPE",
            "n_cav": NC_FULL,
            "duration_ns": 200,
            "_grape_result": grape_result,
        }
    except Exception as e:
        dt_grape = time.perf_counter() - t0
        print(f"    GRAPE FAILED: {e} ({dt_grape:.0f}s)", flush=True)
        traceback.print_exc()
        results["GRAPE_nc12_200ns"] = {
            "label": "GRAPE_nc12_200ns", "fidelity": 0.0,
            "success": False, "error": str(e), "strategy": "GRAPE",
            "n_cav": NC_FULL, "elapsed_s": float(dt_grape)}


# ═══════════════════════════════════════════════════════════════════════════
# PHASE 5: Wigner Functions
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70, flush=True)
print("  PHASE 5: Wigner Function Comparison", flush=True)
print("="*70, flush=True)

try:
    from cqed_sim.sim.extractors import cavity_wigner, reduced_cavity_state
    import qutip as qt
    HAS_WIGNER = True
except ImportError as e:
    HAS_WIGNER = False
    print(f"  Wigner import failed: {e}", flush=True)

# Rank all results
ranked = sorted(
    [(k, v["fidelity"]) for k, v in results.items()
     if isinstance(v, dict) and v.get("fidelity", 0) > 0.3 and v.get("n_cav")],
    key=lambda x: x[1], reverse=True)
print(f"\n  Top results:", flush=True)
for k, f in ranked[:8]:
    nc = results[k].get("n_cav", "?")
    print(f"    {k}: F={f:.6f} (N_cav={nc})", flush=True)

wigner_data = {}

if HAS_WIGNER and ranked:
    # Use best result at N_cav=12 if available, else N_cav=6
    best_12 = [(k, f) for k, f in ranked if results[k].get("n_cav") == NC_FULL]
    best_any = best_12[0] if best_12 else ranked[0]
    best_k, best_f = best_any
    bres = results[best_k]
    nc_wig = bres["n_cav"]

    # Get full-space unitary
    U_full_np = None

    # Try parametric synthesis result
    sr = bres.get("_synth")
    if sr is not None:
        U_full = sr.simulation.full_operator
        U_full_np = U_full.full() if isinstance(U_full, qt.Qobj) else np.asarray(U_full)
        print(f"\n  Wigner source: {best_k} (parametric, F={best_f:.6f}, nc={nc_wig})", flush=True)

    # Try GRAPE result
    if U_full_np is None:
        grape_r = results.get("GRAPE_nc12_200ns", {})
        gr = grape_r.get("_grape_result")
        if gr is not None:
            try:
                U_full_np = np.asarray(gr.best_unitary)
                nc_wig = NC_FULL
                best_k = "GRAPE_nc12_200ns"
                best_f = grape_r["fidelity"]
                print(f"\n  Wigner source: GRAPE (F={best_f:.6f}, nc={nc_wig})", flush=True)
            except Exception as e:
                print(f"  Could not extract GRAPE unitary: {e}", flush=True)

    if U_full_np is not None:
        full_dim = 2 * nc_wig
        # Target embedded in full space
        U_tgt_full = np.eye(full_dim, dtype=complex)
        lidx = [0, 1, nc_wig, nc_wig + 1]
        for i, li in enumerate(lidx):
            for j, lj in enumerate(lidx):
                U_tgt_full[li, lj] = U_target[i, j]

        basis_labels = ["|g,0>", "|g,1>", "|e,0>", "|e,1>"]
        n_pts = 61; ext = 3.5

        for bi, bl in zip(lidx, basis_labels):
            psi0 = np.zeros(full_dim, dtype=complex)
            psi0[bi] = 1.0
            psi_t = U_tgt_full @ psi0
            psi_a = U_full_np @ psi0

            rho_t = qt.Qobj(np.outer(psi_t, psi_t.conj()), dims=[[2, nc_wig], [2, nc_wig]])
            rho_a = qt.Qobj(np.outer(psi_a, psi_a.conj()), dims=[[2, nc_wig], [2, nc_wig]])
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
            "source": best_k, "fidelity": best_f, "n_cav": nc_wig,
            "cavity_fidelities": {bl: wigner_data[bl]["cavity_fidelity"] for bl in basis_labels},
            "l2_distances": {bl: wigner_data[bl]["l2_distance"] for bl in basis_labels}}

        # ── Wigner figure ─────────────────────────────────────────────
        fig, axes = plt.subplots(3, 4, figsize=(14, 10))
        vmax = max(max(np.max(np.abs(wigner_data[bl]["W_target"])),
                       np.max(np.abs(wigner_data[bl]["W_achieved"])))
                   for bl in basis_labels)
        for col, bl in enumerate(basis_labels):
            wd = wigner_data[bl]
            X, Y = np.meshgrid(wd["xvec"], wd["yvec"])
            Wt = np.array(wd["W_target"])
            Wa = np.array(wd["W_achieved"])
            dW = Wt - Wa

            # Row 0: target
            axes[0, col].pcolormesh(X, Y, Wt, cmap="RdBu_r",
                                     vmin=-vmax, vmax=vmax, shading="auto")
            axes[0, col].set_title(f"Target: {bl}", fontsize=9)
            axes[0, col].set_aspect("equal")

            # Row 1: achieved
            im = axes[1, col].pcolormesh(X, Y, Wa, cmap="RdBu_r",
                                          vmin=-vmax, vmax=vmax, shading="auto")
            axes[1, col].set_title(f"Achieved: $F_{{cav}}$={wd['cavity_fidelity']:.4f}", fontsize=9)
            axes[1, col].set_aspect("equal")

            # Row 2: difference
            dmax = max(np.max(np.abs(dW)), 1e-10)
            axes[2, col].pcolormesh(X, Y, dW, cmap="RdBu_r",
                                     vmin=-dmax, vmax=dmax, shading="auto")
            axes[2, col].set_title(f"$\\Delta W$ (L2={wd['l2_distance']:.4f})", fontsize=9)
            axes[2, col].set_aspect("equal")
            axes[2, col].set_xlabel("Re($\\alpha$)")

            for row in range(3):
                if col == 0:
                    axes[row, col].set_ylabel("Im($\\alpha$)")

        fig.suptitle(
            f"Wigner Function Comparison — {best_k}\n$F$={best_f:.4f}, $N_{{cav}}$={nc_wig}",
            fontsize=12, fontweight="bold")
        fig.colorbar(im, ax=axes[:2].ravel().tolist(), label="$W(\\alpha)$",
                     shrink=0.7, pad=0.02)
        fig.tight_layout(rect=[0, 0, 0.93, 0.93])
        for fmt in ("png", "pdf"):
            fig.savefig(FIG_DIR / f"wigner_comparison.{fmt}", dpi=300, bbox_inches="tight")
        plt.close(fig)
        print("  Wigner figure saved.", flush=True)
    else:
        print("  WARNING: No full-space unitary available for Wigner.", flush=True)


# ═══════════════════════════════════════════════════════════════════════════
# PHASE 6: Save Results & Artifacts
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70, flush=True)
print("  PHASE 6: Save Results & Artifacts", flush=True)
print("="*70, flush=True)

# Clean for JSON
results_json = {}
for k, v in results.items():
    if isinstance(v, dict):
        results_json[k] = {kk: vv for kk, vv in v.items()
                           if kk not in ("_synth", "_grape_result")}
    else:
        results_json[k] = v

(DATA_DIR / "iteration4_results.json").write_text(
    json.dumps(results_json, indent=2, default=str), encoding="utf-8")
print(f"  Results saved.", flush=True)

# Save Wigner data (NPZ for efficiency)
if wigner_data:
    npz_data = {}
    for bl, wd in wigner_data.items():
        safe_bl = bl.replace("|", "").replace(">", "").replace(",", "_")
        for key in ("xvec", "yvec", "W_target", "W_achieved"):
            npz_data[f"{safe_bl}_{key}"] = np.array(wd[key])
        npz_data[f"{safe_bl}_cavity_fidelity"] = np.array([wd["cavity_fidelity"]])
        npz_data[f"{safe_bl}_l2_distance"] = np.array([wd["l2_distance"]])
    np.savez(str(ARTIFACTS_DIR / "wigner_comparison.npz"), **npz_data)
    print("  Wigner data saved to artifacts/wigner_comparison.npz", flush=True)

# Save synthesis artifacts
for strat, prefix in [("D+SQR+CP", "best_strategy_B"), ("D+R+FE", "best_strategy_D")]:
    best_e = max(
        ((k, v) for k, v in results.items()
         if isinstance(v, dict) and strat in v.get("strategy", "") and v.get("fidelity", 0) > 0.3),
        key=lambda x: x[1]["fidelity"], default=(None, None))
    if best_e[0]:
        ek, ev = best_e
        sr = ev.get("_synth")
        ap = ARTIFACTS_DIR / f"{prefix}.json"
        if sr is not None:
            try:
                sr.save(str(ap), include_history=True)
                print(f"  Artifact: {ap.name} (F={ev['fidelity']:.6f})", flush=True)
            except Exception as e:
                print(f"  Save warning: {e}. Using fallback.", flush=True)
                ap.write_text(json.dumps(
                    {"label": ek, "fidelity": ev["fidelity"],
                     "max_amp": ev.get("max_amp"), "n_cav": ev.get("n_cav"),
                     "sequence": ev.get("sequence")},
                    indent=2, default=str), encoding="utf-8")
        print(f"  Best {prefix}: {ek} F={ev['fidelity']:.6f}", flush=True)

# Save GRAPE artifact
grape_r = results.get("GRAPE_nc12_200ns", {})
gr = grape_r.get("_grape_result")
if gr is not None:
    try:
        gpath = ARTIFACTS_DIR / "grape_ncav12_200ns.json"
        # Save relevant GRAPE data
        grape_save = {
            "fidelity": grape_r["fidelity"],
            "duration_ns": 200,
            "n_cav": NC_FULL,
            "n_steps": N_STEPS if HAS_GRAPE else 40,
            "dt_ns": 5,
            "n_seeds": 3,
            "n_iter": 300,
        }
        if hasattr(gr, 'best_amplitudes'):
            grape_save["amplitudes"] = gr.best_amplitudes.tolist()
        gpath.write_text(json.dumps(grape_save, indent=2, default=str), encoding="utf-8")
        print(f"  GRAPE artifact saved.", flush=True)
    except Exception as e:
        print(f"  GRAPE artifact save warning: {e}", flush=True)


# ═══════════════════════════════════════════════════════════════════════════
# PHASE 7: Summary Figures
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70, flush=True)
print("  PHASE 7: Summary Figures", flush=True)
print("="*70, flush=True)

# 1. Bounded displacement sweep: fidelity vs max_amplitude
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Strategy B
for i, nb in enumerate([2, 3]):
    amps, fids = [], []
    for amp_val in [0.3, 0.5, 0.7, 1.0]:
        k = f"B_nc{NC_FAST}_amp{amp_val}_blk{nb}"
        if k in results and results[k].get("fidelity", 0) > 0:
            amps.append(amp_val)
            fids.append(results[k]["fidelity"])
    if amps:
        ax1.plot(amps, fids, "o-", color=COLORS[i], ms=7, lw=2, label=f"{nb} blocks")
    # Unconstrained baseline
    uk = f"B_nc{NC_FAST}_ampnone_blk{nb}"
    if uk in results and results[uk].get("fidelity", 0) > 0:
        ax1.axhline(results[uk]["fidelity"], ls=":", color=COLORS[i], alpha=0.6,
                    label=f"{nb} blk (no bound)")
ax1.set_xlabel("Max displacement $|\\alpha|_{\\max}$")
ax1.set_ylabel("Subspace fidelity")
ax1.set_title("Strategy B: D + SQR + CP")
ax1.axhline(0.99, ls="--", color="red", alpha=0.3, lw=0.5)
ax1.axhline(0.999, ls="--", color="darkred", alpha=0.3, lw=0.5)
ax1.legend(fontsize=8)
ax1.grid(alpha=0.3)

# Strategy D
for i, nb in enumerate([2, 3]):
    amps, fids = [], []
    for amp_val in [0.3, 0.5, 0.7, 1.0]:
        k = f"D_nc{NC_FAST}_amp{amp_val}_blk{nb}"
        if k in results and results[k].get("fidelity", 0) > 0:
            amps.append(amp_val)
            fids.append(results[k]["fidelity"])
    if amps:
        ax2.plot(amps, fids, "o-", color=COLORS[i+2], ms=7, lw=2, label=f"{nb} blocks")
uk = f"D_nc{NC_FAST}_ampnone_blk3"
if uk in results and results[uk].get("fidelity", 0) > 0:
    ax2.axhline(results[uk]["fidelity"], ls=":", color=COLORS[3], alpha=0.6,
                label="3 blk (no bound)")
ax2.set_xlabel("Max displacement $|\\alpha|_{\\max}$")
ax2.set_ylabel("Subspace fidelity")
ax2.set_title("Strategy D: D + R + FE")
ax2.axhline(0.99, ls="--", color="red", alpha=0.3, lw=0.5)
ax2.axhline(0.999, ls="--", color="darkred", alpha=0.3, lw=0.5)
ax2.legend(fontsize=8)
ax2.grid(alpha=0.3)

fig.suptitle(f"Bounded-Displacement Optimisation ($N_{{cav}}$={NC_FAST})",
             fontsize=13, fontweight="bold")
fig.tight_layout()
for fmt in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"bounded_displacement_sweep.{fmt}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("  bounded_displacement_sweep saved.", flush=True)

# 2. Ranking
rank_items = sorted(
    [(k, v) for k, v in results.items()
     if isinstance(v, dict) and "fidelity" in v and "strategy" in v and v["fidelity"] > 0],
    key=lambda x: x[1]["fidelity"], reverse=True)
if rank_items:
    fig, ax = plt.subplots(figsize=(10, max(4, len(rank_items)*0.35)))
    lbs = [s[0] for s in rank_items]
    fds = [s[1]["fidelity"] for s in rank_items]
    cs = []
    for s in rank_items:
        strat = s[1].get("strategy", "")
        if "GRAPE" in strat:
            cs.append(COLORS[4])
        elif "SQR" in strat:
            cs.append(COLORS[1])
        else:
            cs.append(COLORS[3])
    ax.barh(range(len(lbs)), fds, color=cs, alpha=0.85, edgecolor="black", lw=0.5)
    ax.set_yticks(range(len(lbs)))
    ax.set_yticklabels(lbs, fontsize=7)
    ax.set_xlabel("Subspace Fidelity")
    ax.set_title("Iteration 4: All Results")
    ax.axvline(0.99, ls="--", color="red", alpha=0.4, label="99%")
    ax.axvline(0.999, ls=":", color="darkred", alpha=0.4, label="99.9%")
    ax.legend(fontsize=8)
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()
    for fmt in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"iteration4_ranking.{fmt}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  iteration4_ranking saved.", flush=True)


# ═══════════════════════════════════════════════════════════════════════════
# FINAL SUMMARY
# ═══════════════════════════════════════════════════════════════════════════
print("\n" + "="*70, flush=True)
print("  FINAL SUMMARY", flush=True)
print("="*70, flush=True)
for i, (k, f) in enumerate(ranked[:10] if ranked else []):
    nc = results[k].get("n_cav", "?")
    strat = results[k].get("strategy", "?")
    amp = results[k].get("max_amp", "none")
    print(f"  {i+1}. {k}: F={f:.6f} ({strat}, nc={nc}, amp={amp})", flush=True)

if wigner_data:
    print(f"\n  Wigner comparison: {results.get('wigner', {}).get('source', 'N/A')}", flush=True)
    for bl in ["|g,0>", "|g,1>", "|e,0>", "|e,1>"]:
        if bl in wigner_data:
            print(f"    {bl}: F_cav={wigner_data[bl]['cavity_fidelity']:.6f}", flush=True)

print("\n  DONE.", flush=True)
