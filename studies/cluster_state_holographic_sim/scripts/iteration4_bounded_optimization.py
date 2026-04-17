"""Iteration 4: Bounded-displacement re-optimization at N_cav=12.

Key ideas from user:
  1. Bound displacement amplitude to keep cavity population low (n=0,1 only).
  2. Work at N_cav=12 for physical truncation.
  3. Try to minimize SQR/SNAP gate duration via MultiObjective.
  4. Compute Wigner functions for best result vs target.
  5. Save all optimized artifacts as JSON.

Physics reasoning: With |alpha| small, D(alpha)|0> stays mostly in n=0,1.
SQR gates can then effectively control the 4-level logical subspace
without leaking into n>=2:
  |alpha|=0.3 -> P(n>=2) ~ 0.4%
  |alpha|=0.5 -> P(n>=2) ~ 2.6%
  |alpha|=0.7 -> P(n>=2) ~ 9%
  |alpha|=1.0 -> P(n>=2) ~ 26%

Output:
  data/iteration4_bounded_optimization.json
  figures/bounded_displacement_sweep.{png,pdf}
  figures/wigner_comparison.{png,pdf}
  figures/duration_optimization.{png,pdf}
  artifacts/best_bounded_strategy_B.json
  artifacts/best_bounded_strategy_D.json
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

# ── cqed_sim imports ──────────────────────────────────────────────────────
from cqed_sim.unitary_synthesis import (
    Displacement, QubitRotation, SQR, ConditionalPhaseSQR,
    FreeEvolveCondPhase, SNAP, Subspace, TargetUnitary,
    UnitarySynthesizer, GateSequence, DriftPhaseModel,
    LeakagePenalty, MultiObjective, ExecutionOptions,
    SynthesisConstraints,
    subspace_unitary_fidelity, simulate_sequence,
)
from cqed_sim.unitary_synthesis.targets import make_target

# ── constants ─────────────────────────────────────────────────────────────
TWO_PI = 2.0 * np.pi
CHI   = TWO_PI * (-2.84e6)
CHIP  = TWO_PI * (-21e3)
KERR  = TWO_PI * (-28e3)
OMEGA_Q = TWO_PI * 6.150e9
OMEGA_C = TWO_PI * 5.241e9
ALPHA   = TWO_PI * (-255e6)

N_CAV = 12             # physical truncation requested by user
FULL_DIM = 2 * N_CAV   # qubit (2) x cavity (12) = 24
N_TR = 2               # transmon levels

LOGICAL_LABELS = ["|g,0>", "|g,1>", "|e,0>", "|e,1>"]
subspace = Subspace.custom(FULL_DIM, [0, 1, N_CAV, N_CAV + 1], LOGICAL_LABELS)

# Drift models
no_drift = DriftPhaseModel(chi=0.0, chi2=0.0, kerr=0.0)
physical_drift = DriftPhaseModel(chi=abs(CHI), chi2=abs(CHIP), kerr=abs(KERR))

COLORS = ['#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE', '#AA3377', '#BBBBBB']

# ── target unitary ────────────────────────────────────────────────────────
print("=" * 70)
print("  Iteration 4: Bounded-Displacement Optimization at N_cav=12")
print("=" * 70)

U_target = make_target("cluster", n_match=1)
print(f"Target shape: {U_target.shape}")
target = TargetUnitary(U_target, ignore_global_phase=True)

results = {}
log_lines = []

def log(msg: str):
    """Print and record log line."""
    print(msg)
    log_lines.append(msg)


# ═══════════════════════════════════════════════════════════════════════════
# Part 1: Bounded-Displacement Sweep for Strategy B (D + SQR + CP)
# ═══════════════════════════════════════════════════════════════════════════
log("\n" + "=" * 70)
log("  Part 1: Strategy B (D + SQR + CP) — Bounded Displacement Sweep")
log("=" * 70)

MAX_AMP_VALUES = [0.3, 0.5, 0.7, 1.0]
BLOCK_COUNTS_B = [2, 3, 4]

def build_strategy_B(n_blocks: int, n_cav: int) -> GateSequence:
    """Build D+SQR+CP sequence with n_blocks at given n_cav."""
    gates = []
    # SQR parameters: theta_n and phi_n of length n_cav
    theta_init = [0.0] * n_cav
    theta_init[0] = np.pi / 2
    if n_cav > 1:
        theta_init[1] = np.pi / 4
    phi_init = [0.0] * n_cav
    phases_init = [0.0] * n_cav

    for i in range(n_blocks):
        gates.append(Displacement(name=f"D{i}", alpha=0.3 + 0j, duration=200e-9))
        gates.append(SQR(name=f"S{2*i}", theta_n=theta_init,
                         phi_n=phi_init, drift_model=no_drift, duration=400e-9))
        gates.append(ConditionalPhaseSQR(name=f"CP{i}", phases_n=phases_init,
                         drift_model=no_drift, duration=200e-9))
        gates.append(SQR(name=f"S{2*i+1}", theta_n=theta_init,
                         phi_n=phi_init, drift_model=no_drift, duration=400e-9))
    gates.append(Displacement(name=f"D{n_blocks}", alpha=0.3 + 0j, duration=200e-9))
    return GateSequence(gates=gates, n_cav=n_cav)


def run_bounded_synthesis(
    sequence: GateSequence,
    label: str,
    max_amplitude: float | None = None,
    multistart: int = 6,
    maxiter: int = 500,
    duration_weight: float = 0.0,
) -> dict:
    """Run synthesis with optional amplitude bound and duration objective."""
    log(f"\n  [{label}] max_amp={max_amplitude}, multistart={multistart}, maxiter={maxiter}")
    constraints = None
    if max_amplitude is not None:
        constraints = SynthesisConstraints(max_amplitude=max_amplitude)

    objectives = MultiObjective(
        fidelity_weight=1.0,
        leakage_weight=0.05,
        duration_weight=duration_weight,
    )

    t0 = time.perf_counter()
    try:
        synth = UnitarySynthesizer(
            primitives=sequence.gates,
            subspace=subspace,
            objectives=objectives,
            leakage_penalty=LeakagePenalty(weight=0.05),
            synthesis_constraints=constraints,
            execution=ExecutionOptions(engine="auto", use_fast_path=True),
        )
        result = synth.fit(
            target=target,
            init_guess="heuristic",
            multistart=multistart,
            maxiter=maxiter,
        )
        dt = time.perf_counter() - t0
        F = subspace_unitary_fidelity(
            result.simulation.subspace_operator,
            U_target,
            gauge="global",
        )
        log(f"    F_proj = {F:.6f}   objective = {result.objective:.6f}  ({dt:.1f}s)")

        # Extract displacement amplitudes from optimized sequence
        disp_amplitudes = []
        for g in result.sequence.gates:
            if hasattr(g, 'alpha'):
                disp_amplitudes.append(float(abs(g.alpha)))

        return {
            "label": label,
            "fidelity": float(F),
            "objective": float(result.objective),
            "success": bool(result.success),
            "elapsed_s": float(dt),
            "max_amplitude_constraint": max_amplitude,
            "disp_amplitudes": disp_amplitudes,
            "max_disp_achieved": float(max(disp_amplitudes)) if disp_amplitudes else 0.0,
            "sequence": result.sequence.serialize(),
            "_result": result,  # keep for artifact saving
        }
    except Exception as e:
        dt = time.perf_counter() - t0
        log(f"    FAILED: {e}")
        traceback.print_exc()
        return {
            "label": label,
            "fidelity": 0.0,
            "objective": float('inf'),
            "success": False,
            "elapsed_s": float(dt),
            "max_amplitude_constraint": max_amplitude,
            "error": str(e),
        }


# Run Strategy B sweep
for max_amp in MAX_AMP_VALUES:
    for n_blocks in BLOCK_COUNTS_B:
        seq = build_strategy_B(n_blocks, N_CAV)
        label = f"B_ncav{N_CAV}_amp{max_amp}_blocks{n_blocks}"
        res = run_bounded_synthesis(seq, label, max_amplitude=max_amp)
        res["strategy"] = "D + SQR + CP"
        res["n_blocks"] = n_blocks
        res["n_cav"] = N_CAV
        results[label] = res


# ═══════════════════════════════════════════════════════════════════════════
# Part 2: Bounded-Displacement Sweep for Strategy D (D + R + FE)
# ═══════════════════════════════════════════════════════════════════════════
log("\n" + "=" * 70)
log("  Part 2: Strategy D (D + R + FE) — Bounded Displacement Sweep")
log("=" * 70)

BLOCK_COUNTS_D = [2, 3, 4]

def build_strategy_D(n_blocks: int, n_cav: int) -> GateSequence:
    """Build D+R+FE sequence with n_blocks at given n_cav."""
    gates = []
    for i in range(n_blocks):
        gates.append(Displacement(name=f"D{i}", alpha=0.3 + 0j, duration=200e-9))
        gates.append(QubitRotation(name=f"R{i}", theta=np.pi / 2, phi=0.0, duration=100e-9))
        gates.append(FreeEvolveCondPhase(
            name=f"FE{i}", duration=200e-9,
            drift_model=DriftPhaseModel(chi=abs(CHI), chi2=0.0, kerr=0.0),
            optimize_time=True,
        ))
    gates.append(Displacement(name=f"D{n_blocks}", alpha=0.3 + 0j, duration=200e-9))
    gates.append(QubitRotation(name=f"R{n_blocks}", theta=np.pi / 2, phi=np.pi / 2, duration=100e-9))
    return GateSequence(gates=gates, n_cav=n_cav)


for max_amp in MAX_AMP_VALUES:
    for n_blocks in BLOCK_COUNTS_D:
        seq = build_strategy_D(n_blocks, N_CAV)
        label = f"D_ncav{N_CAV}_amp{max_amp}_blocks{n_blocks}"
        res = run_bounded_synthesis(seq, label, max_amplitude=max_amp)
        res["strategy"] = "D + R + FreeEvolveCondPhase"
        res["n_blocks"] = n_blocks
        res["n_cav"] = N_CAV
        results[label] = res


# ═══════════════════════════════════════════════════════════════════════════
# Part 3: Duration Optimization for Best Configs
# ═══════════════════════════════════════════════════════════════════════════
log("\n" + "=" * 70)
log("  Part 3: Duration Optimization (minimize gate time for best configs)")
log("=" * 70)

# Find best Strategy B and D results
best_B = max(
    ((k, v) for k, v in results.items() if v.get("strategy") == "D + SQR + CP" and v.get("fidelity", 0) > 0.9),
    key=lambda x: x[1]["fidelity"],
    default=(None, None),
)
best_D = max(
    ((k, v) for k, v in results.items() if v.get("strategy") == "D + R + FreeEvolveCondPhase" and v.get("fidelity", 0) > 0.9),
    key=lambda x: x[1]["fidelity"],
    default=(None, None),
)

# Re-optimize best configs with duration penalty
for tag, (best_key, best_val) in [("B", best_B), ("D", best_D)]:
    if best_key is None:
        log(f"  No high-fidelity result for Strategy {tag}; skipping duration opt.")
        continue
    n_blocks = best_val["n_blocks"]
    max_amp = best_val.get("max_amplitude_constraint")
    log(f"\n  Duration opt for Strategy {tag}: {best_key} (F={best_val['fidelity']:.6f})")

    if tag == "B":
        seq = build_strategy_B(n_blocks, N_CAV)
    else:
        seq = build_strategy_D(n_blocks, N_CAV)

    label = f"{tag}_duration_opt_ncav{N_CAV}"
    res = run_bounded_synthesis(
        seq, label,
        max_amplitude=max_amp,
        multistart=8,
        maxiter=600,
        duration_weight=0.01,  # small weight to gently push for shorter gates
    )
    res["strategy"] = best_val["strategy"]
    res["n_blocks"] = n_blocks
    res["n_cav"] = N_CAV
    res["note"] = "Duration-optimized version of best bounded config"
    results[label] = res


# ═══════════════════════════════════════════════════════════════════════════
# Part 4: Wigner Function Comparison
# ═══════════════════════════════════════════════════════════════════════════
log("\n" + "=" * 70)
log("  Part 4: Wigner Function Comparison (target vs achieved)")
log("=" * 70)

try:
    from cqed_sim.sim.extractors import cavity_wigner, reduced_cavity_state
    HAS_WIGNER = True
    log("  Wigner function API imported successfully.")
except ImportError as e:
    HAS_WIGNER = False
    log(f"  WARNING: Could not import Wigner API: {e}")

# Find overall best result
all_fid = [(k, v["fidelity"]) for k, v in results.items()
           if isinstance(v, dict) and "fidelity" in v and v["fidelity"] > 0.5]
all_fid.sort(key=lambda x: x[1], reverse=True)
log(f"\n  Top 5 results:")
for k, f in all_fid[:5]:
    log(f"    {k}: F={f:.6f}")

if HAS_WIGNER and all_fid:
    import qutip as qt

    best_key, best_fid = all_fid[0]
    best_res = results[best_key]
    log(f"\n  Computing Wigner functions for: {best_key} (F={best_fid:.6f})")

    # Get the full-space unitary from the synthesis result
    synth_result = best_res.get("_result")
    if synth_result is not None:
        U_full = synth_result.simulation.full_operator  # 24x24 numpy array
        if isinstance(U_full, qt.Qobj):
            U_full_np = U_full.full()
        else:
            U_full_np = np.asarray(U_full)
        log(f"  Full operator shape: {U_full_np.shape}")

        # Build target unitary embedded in full space
        U_target_full = np.eye(FULL_DIM, dtype=complex)
        logical_indices = [0, 1, N_CAV, N_CAV + 1]
        for i, li in enumerate(logical_indices):
            for j, lj in enumerate(logical_indices):
                U_target_full[li, lj] = U_target[i, j]

        # Computational basis states in full 24-dim Hilbert space
        basis_labels = ["|g,0>", "|g,1>", "|e,0>", "|e,1>"]
        basis_indices = logical_indices

        wigner_data = {}
        n_points = 51
        extent = 3.0

        for idx, (bi, bl) in enumerate(zip(basis_indices, basis_labels)):
            # Initial state
            psi0 = np.zeros(FULL_DIM, dtype=complex)
            psi0[bi] = 1.0

            # Target output state
            psi_target = U_target_full @ psi0

            # Achieved output state
            psi_achieved = U_full_np @ psi0

            # Convert to density matrices in full space
            rho_target = np.outer(psi_target, psi_target.conj())
            rho_achieved = np.outer(psi_achieved, psi_achieved.conj())

            # Convert to qutip Qobj with tensor structure [qubit, cavity]
            rho_target_qt = qt.Qobj(rho_target, dims=[[2, N_CAV], [2, N_CAV]])
            rho_achieved_qt = qt.Qobj(rho_achieved, dims=[[2, N_CAV], [2, N_CAV]])

            # Reduced cavity states (trace out qubit)
            rho_cav_target = reduced_cavity_state(rho_target_qt)
            rho_cav_achieved = reduced_cavity_state(rho_achieved_qt)

            # Compute Wigner functions
            xvec_t, yvec_t, W_target = cavity_wigner(
                rho_cav_target, n_points=n_points, extent=extent
            )
            xvec_a, yvec_a, W_achieved = cavity_wigner(
                rho_cav_achieved, n_points=n_points, extent=extent
            )

            # Wigner overlap (L2 distance)
            dW = W_target - W_achieved
            dx = xvec_t[1] - xvec_t[0] if len(xvec_t) > 1 else 1.0
            dy = yvec_t[1] - yvec_t[0] if len(yvec_t) > 1 else 1.0
            l2_dist = np.sqrt(np.sum(dW**2) * dx * dy)

            # State fidelity (cavity)
            if rho_cav_target.isket:
                fid_cav = float(abs(rho_cav_achieved.overlap(rho_cav_target))**2)
            else:
                fid_cav = float(qt.fidelity(rho_cav_target, rho_cav_achieved)**2)

            log(f"    {bl}: cavity fidelity={fid_cav:.6f}, Wigner L2 dist={l2_dist:.6f}")

            wigner_data[bl] = {
                "xvec": xvec_t.tolist(),
                "yvec": yvec_t.tolist(),
                "W_target": W_target.tolist(),
                "W_achieved": W_achieved.tolist(),
                "cavity_fidelity": fid_cav,
                "wigner_l2_distance": l2_dist,
            }

        results["wigner_comparison"] = {
            "best_strategy": best_key,
            "best_fidelity": best_fid,
            "basis_states": basis_labels,
            "cavity_fidelities": {bl: wigner_data[bl]["cavity_fidelity"] for bl in basis_labels},
            "wigner_l2_distances": {bl: wigner_data[bl]["wigner_l2_distance"] for bl in basis_labels},
        }

        # ── Plot Wigner comparison ────────────────────────────────────────
        fig, axes = plt.subplots(2, 4, figsize=(14, 7))
        vmax = max(
            max(abs(wigner_data[bl]["W_target"]).max(), abs(wigner_data[bl]["W_achieved"]).max())
            for bl in basis_labels
        ) if wigner_data else 1.0

        for col, bl in enumerate(basis_labels):
            wd = wigner_data[bl]
            X, Y = np.meshgrid(wd["xvec"], wd["yvec"])
            W_t = np.array(wd["W_target"])
            W_a = np.array(wd["W_achieved"])

            # Target (top row)
            im_t = axes[0, col].pcolormesh(X, Y, W_t, cmap="RdBu_r",
                                            vmin=-vmax, vmax=vmax, shading="auto")
            axes[0, col].set_title(f"Target: {bl}", fontsize=9)
            axes[0, col].set_aspect("equal")
            if col == 0:
                axes[0, col].set_ylabel("Im($\\alpha$)")

            # Achieved (bottom row)
            im_a = axes[1, col].pcolormesh(X, Y, W_a, cmap="RdBu_r",
                                            vmin=-vmax, vmax=vmax, shading="auto")
            axes[1, col].set_title(f"Achieved: F_cav={wd['cavity_fidelity']:.4f}", fontsize=9)
            axes[1, col].set_aspect("equal")
            axes[1, col].set_xlabel("Re($\\alpha$)")
            if col == 0:
                axes[1, col].set_ylabel("Im($\\alpha$)")

        fig.suptitle(
            f"Wigner Function Comparison — {best_key} (F={best_fid:.4f})",
            fontsize=12, fontweight="bold",
        )
        fig.colorbar(im_t, ax=axes.ravel().tolist(), label="$W(\\alpha)$",
                     shrink=0.8, pad=0.02)
        fig.tight_layout(rect=[0, 0, 0.92, 0.95])
        for fmt in ("png", "pdf"):
            fig.savefig(FIG_DIR / f"wigner_comparison.{fmt}", dpi=300, bbox_inches="tight")
        plt.close(fig)
        log("  Wigner comparison figure saved.")
    else:
        log("  WARNING: No SynthesisResult object available for Wigner computation.")


# ═══════════════════════════════════════════════════════════════════════════
# Part 5: Save Results and Artifacts
# ═══════════════════════════════════════════════════════════════════════════
log("\n" + "=" * 70)
log("  Part 5: Save Results and Artifacts")
log("=" * 70)

# Remove non-serializable _result objects before saving JSON
results_clean = {}
for k, v in results.items():
    if isinstance(v, dict):
        results_clean[k] = {kk: vv for kk, vv in v.items() if kk != "_result"}
    else:
        results_clean[k] = v

out_path = DATA_DIR / "iteration4_bounded_optimization.json"
out_path.write_text(json.dumps(results_clean, indent=2, default=str), encoding="utf-8")
log(f"  Results saved to {out_path}")

# Save SynthesisResult artifacts for top results
for tag, prefix in [("D + SQR + CP", "best_bounded_strategy_B"),
                     ("D + R + FreeEvolveCondPhase", "best_bounded_strategy_D")]:
    best_entry = max(
        ((k, v) for k, v in results.items()
         if isinstance(v, dict) and v.get("strategy") == tag and v.get("fidelity", 0) > 0.5),
        key=lambda x: x[1]["fidelity"],
        default=(None, None),
    )
    if best_entry[0] is not None:
        entry_key, entry_val = best_entry
        synth_res = entry_val.get("_result")
        if synth_res is not None:
            artifact_path = ARTIFACTS_DIR / f"{prefix}.json"
            try:
                synth_res.save(str(artifact_path), include_history=True)
                log(f"  Saved artifact: {artifact_path.name} (F={entry_val['fidelity']:.6f})")
            except Exception as e:
                log(f"  WARNING: Could not save artifact {prefix}: {e}")
                # Fallback: save sequence as JSON
                fallback_data = {
                    "label": entry_key,
                    "fidelity": entry_val["fidelity"],
                    "strategy": tag,
                    "sequence": entry_val.get("sequence"),
                    "n_cav": N_CAV,
                    "max_amplitude_constraint": entry_val.get("max_amplitude_constraint"),
                }
                artifact_path.write_text(
                    json.dumps(fallback_data, indent=2, default=str), encoding="utf-8"
                )
                log(f"  Saved fallback artifact: {artifact_path.name}")


# ═══════════════════════════════════════════════════════════════════════════
# Part 6: Generate Summary Figures
# ═══════════════════════════════════════════════════════════════════════════
log("\n" + "=" * 70)
log("  Part 6: Summary Figures")
log("=" * 70)

# Figure: Bounded displacement sweep — fidelity vs max_amplitude for each block count
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

for ax, strategy_tag, block_counts, title in [
    (ax1, "D + SQR + CP", BLOCK_COUNTS_B, "Strategy B: D + SQR + CP"),
    (ax2, "D + R + FreeEvolveCondPhase", BLOCK_COUNTS_D, "Strategy D: D + R + FE"),
]:
    for i, nb in enumerate(block_counts):
        amps = []
        fids = []
        for max_amp in MAX_AMP_VALUES:
            key = None
            for k, v in results.items():
                if (isinstance(v, dict)
                    and v.get("strategy") == strategy_tag
                    and v.get("n_blocks") == nb
                    and v.get("max_amplitude_constraint") == max_amp
                    and v.get("n_cav") == N_CAV
                    and "duration_opt" not in k):
                    key = k
                    break
            if key:
                amps.append(max_amp)
                fids.append(results[key]["fidelity"])
        if amps:
            ax.plot(amps, fids, "o-", color=COLORS[i], ms=7, lw=2,
                    label=f"{nb} blocks")

    ax.set_xlabel("Max displacement amplitude $|\\alpha|_{\\max}$")
    ax.set_ylabel("Subspace fidelity")
    ax.set_title(title)
    ax.axhline(0.99, ls="--", color="red", alpha=0.4, label="99%")
    ax.axhline(0.999, ls=":", color="darkred", alpha=0.4, label="99.9%")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    ax.set_ylim(0, 1.05)

fig.suptitle(f"Bounded-Displacement Optimization at $N_{{cav}}={N_CAV}$",
             fontsize=13, fontweight="bold")
fig.tight_layout()
for fmt in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"bounded_displacement_sweep.{fmt}", dpi=300, bbox_inches="tight")
plt.close(fig)
log("  bounded_displacement_sweep figure saved.")

# Figure: Combined ranking of all Iter4 results
strat_items = sorted(
    [(k, v) for k, v in results.items()
     if isinstance(v, dict) and "fidelity" in v and "strategy" in v],
    key=lambda x: x[1]["fidelity"],
    reverse=True,
)
if strat_items:
    fig, ax = plt.subplots(figsize=(12, 6))
    labels = [s[0] for s in strat_items[:20]]  # top 20
    fids = [s[1]["fidelity"] for s in strat_items[:20]]
    cat_colors = {
        "D + SQR + CP": COLORS[1],
        "D + R + FreeEvolveCondPhase": COLORS[3],
    }
    colors = [cat_colors.get(s[1].get("strategy", ""), COLORS[6]) for s in strat_items[:20]]
    ax.barh(range(len(labels)), fids, color=colors, alpha=0.85,
            edgecolor="black", lw=0.5)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Subspace Fidelity")
    ax.set_title(f"Iteration 4 Strategy Ranking — N_cav={N_CAV}")
    ax.axvline(0.99, ls="--", color="red", alpha=0.5, label="99%")
    ax.axvline(0.999, ls=":", color="darkred", alpha=0.5, label="99.9%")
    ax.legend(fontsize=8)
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.3)
    fig.tight_layout()
    for fmt in ("png", "pdf"):
        fig.savefig(FIG_DIR / f"iteration4_ranking.{fmt}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    log("  iteration4_ranking figure saved.")


# ═══════════════════════════════════════════════════════════════════════════
# Part 7: Summary
# ═══════════════════════════════════════════════════════════════════════════
log("\n" + "=" * 70)
log("  SUMMARY")
log("=" * 70)

for strategy_tag, strat_name in [("D + SQR + CP", "Strategy B"),
                                  ("D + R + FreeEvolveCondPhase", "Strategy D")]:
    log(f"\n  {strat_name}:")
    strat_results = sorted(
        [(k, v) for k, v in results.items()
         if isinstance(v, dict) and v.get("strategy") == strategy_tag],
        key=lambda x: x[1].get("fidelity", 0),
        reverse=True,
    )
    for k, v in strat_results[:5]:
        ma = v.get("max_amplitude_constraint", "?")
        nb = v.get("n_blocks", "?")
        f = v.get("fidelity", 0)
        t = v.get("elapsed_s", 0)
        log(f"    {k}: F={f:.6f} (|α|_max={ma}, {nb} blocks, {t:.0f}s)")

log(f"\nOverall best: {all_fid[0][0]} = {all_fid[0][1]:.6f}" if all_fid else "\nNo results.")

# Save log
log_path = DATA_DIR / "iteration4_log.txt"
log_path.write_text("\n".join(log_lines), encoding="utf-8")
log(f"\nLog saved to {log_path}")

print("\n  Iteration 4 complete.")
