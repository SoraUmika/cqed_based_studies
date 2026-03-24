"""Phase 4 — Variational synthesis with physical drift.

Two primitive libraries are tested, each at multiple depths:

  Library 1: {D(α), QubitRotation, SQR_n}
  Library 2: {D(α), QubitRotation, ConditionalPhaseSQR}

Key finding from Phase 2: U_target requires Fock-level mixing (H_c and
CNOT_{q→c}); Displacement is the essential primitive in both libraries.

Optimization metric
-------------------
Minimize: (1 - F_proj) + λ_leak · L
where F_proj = |Tr(U†_target · U_log)| / 4

Hyperparameters
---------------
  multistart:        8
  maxiter:         500
  cavity truncation: N_CAV = 8

Basis convention: qubit-first {|g,0⟩,|g,1⟩,|e,0⟩,|e,1⟩}.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--only", default=None,
                    help="Run only this ansatz label (e.g. L1a_D_R_SQR_d7)")
_args = parser.parse_args()

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── path setup ─────────────────────────────────────────────────────────────
STUDY_ROOT = Path(__file__).resolve().parents[2]
COMPONENT_NAME = Path(__file__).resolve().parent.name
SIM_ROOT   = Path("C:/Users/dazzl/Box/Shyam Shankar Quantum Circuits Group"
                  "/Users/Users_JianJun/cQED_simulation")
if str(SIM_ROOT) not in sys.path:
    sys.path.insert(0, str(SIM_ROOT))

from cqed_sim.unitary_synthesis import (
    ConditionalPhaseSQR, Displacement, DriftPhaseModel, GateSequence,
    LeakagePenalty, MultiObjective, ExecutionOptions,
    QubitRotation, SQR, Subspace,
    TargetUnitary, UnitarySynthesizer,
    leakage_metrics, subspace_unitary_fidelity,
)

DATA_DIR = STUDY_ROOT / "data" / COMPONENT_NAME
FIG_DIR  = STUDY_ROOT / "figures" / COMPONENT_NAME
DATA_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ── Physical parameters ──────────────────────────────────────────────────────
CHI  = 2 * np.pi * (-2.84e6)   # rad/s
CHIP = 2 * np.pi * (-21e3)     # rad/s  (χ′)
KERR = 2 * np.pi * (-28e3)     # rad/s  (K)

N_CAV    = 8
FULL_DIM = 2 * N_CAV

# Logical subspace in qubit-first ordering: |g,n⟩=n, |e,n⟩=N_CAV+n
LOGICAL_LABELS = ["|g,0>", "|g,1>", "|e,0>", "|e,1>"]
subspace = Subspace.custom(FULL_DIM, [0, 1, N_CAV, N_CAV + 1], LOGICAL_LABELS)

# Drift model
drift    = DriftPhaseModel(chi=CHI, chi2=CHIP, kerr=KERR)

# Target
s = 1.0 / np.sqrt(2)
U_target = np.array([
    [ s, 0,  s, 0],
    [ s, 0, -s, 0],
    [ 0, s,  0, s],
    [ 0,-s,  0, s],
], dtype=np.complex128)
target = TargetUnitary(U_target, ignore_global_phase=True)

# Hyperparameters
MULTISTART   = 2
MAXITER      = 150
LAMBDA_LEAK  = 0.1

# Gate durations
T_ROT  = 100e-9    # 100 ns qubit rotation
T_SQR  = 2000e-9   # 2 μs SQR
T_DISP = 200e-9    # 200 ns displacement
T_CP   = 100e-9    # 100 ns ConditionalPhaseSQR

# Load any previous checkpoint
_ckpt = DATA_DIR / "phase4_results.json"
if _ckpt.exists():
    with open(_ckpt) as _f:
        all_results: dict = json.load(_f)
    print(f"[resume] loaded {len(all_results)} existing results from checkpoint")
else:
    all_results: dict = {}


# ════════════════════════════════════════════════════════════════════════════
# Helper: run synthesis and extract metrics
# ════════════════════════════════════════════════════════════════════════════
def run_synth(sequence: GateSequence, label: str) -> dict:
    print(f"\n  Synthesizing: {label}")
    print(f"    Gates: {[g.name for g in sequence.gates]}")
    t0 = time.perf_counter()

    synth = UnitarySynthesizer(
        primitives=sequence.gates,
        subspace=subspace,
        objectives=MultiObjective(
            fidelity_weight=1.0,
            leakage_weight=LAMBDA_LEAK,
        ),
        leakage_penalty=LeakagePenalty(weight=LAMBDA_LEAK),
        execution=ExecutionOptions(engine="auto", use_fast_path=True),
    )

    result = synth.fit(
        target=target,
        init_guess="heuristic",
        multistart=MULTISTART,
        maxiter=MAXITER,
    )
    dt = time.perf_counter() - t0

    U_sub  = result.simulation.subspace_operator
    F_proj = subspace_unitary_fidelity(U_sub, U_target, gauge="global")

    U_full = result.simulation.full_operator
    if U_full is not None:
        lm      = leakage_metrics(U_full, subspace)
        L_avg   = lm.average
        L_worst = lm.worst
    else:
        L_avg = L_worst = float("nan")

    gate_list = [{"type": g.__class__.__name__, "name": g.name}
                 for g in result.sequence.gates]

    rec = {
        "label":          label,
        "F_proj":         float(F_proj),
        "leakage_avg":    float(L_avg),
        "leakage_worst":  float(L_worst),
        "objective":      float(result.objective),
        "success":        bool(result.success),
        "depth":          len(sequence.gates),
        "n_cav":          N_CAV,
        "elapsed_s":      float(dt),
        "gate_sequence":  gate_list,
        "sequence_params": result.sequence.serialize(),
        "hyperparams": {
            "multistart":   MULTISTART,
            "maxiter":      MAXITER,
            "lambda_leak":  LAMBDA_LEAK,
        },
    }
    print(f"    F_proj={F_proj:.5f}  L_avg={L_avg:.4f}  "
          f"L_worst={L_worst:.4f}  t={dt:.1f}s")
    sys.stdout.flush()
    return rec


# ════════════════════════════════════════════════════════════════════════════
# LIBRARY 1: {D(α), QubitRotation, SQR_n}
#
# D provides Fock-level mixing; R_q gives unconditional qubit rotations;
# SQR_n gives Fock-conditioned qubit rotations.
# ════════════════════════════════════════════════════════════════════════════
print("\n══════════ Library 1: D + QubitRotation + SQR ══════════")

# L1a: depth-7   D R SQR D R SQR D
seq_L1a = GateSequence(gates=[
    Displacement(name="D1",  alpha=0.3+0j, duration=T_DISP),
    QubitRotation(name="R1", theta=np.pi/2, phi=0.0, duration=T_ROT),
    SQR(name="S1", theta_n=[0.0]*N_CAV, phi_n=[0.0]*N_CAV,
        drift_model=drift, duration=T_SQR),
    Displacement(name="D2",  alpha=0.3+0j, duration=T_DISP),
    QubitRotation(name="R2", theta=np.pi/2, phi=0.0, duration=T_ROT),
    SQR(name="S2", theta_n=[np.pi/2]*N_CAV, phi_n=[0.0]*N_CAV,
        drift_model=drift, duration=T_SQR),
    Displacement(name="D3",  alpha=0.3+0j, duration=T_DISP),
], n_cav=N_CAV)

if "L1a_D_R_SQR_d7" not in all_results and (_args.only is None or _args.only == "L1a_D_R_SQR_d7"):
    all_results["L1a_D_R_SQR_d7"] = run_synth(seq_L1a, "L1a_D_R_SQR_d7")
    with open(DATA_DIR / "phase4_results.json", "w") as _f:
        json.dump(all_results, _f, indent=2, default=str)

# L1b: depth-9   D R SQR R D R SQR R D
seq_L1b = GateSequence(gates=[
    Displacement(name="D1",  alpha=0.3+0j, duration=T_DISP),
    QubitRotation(name="R1", theta=np.pi/2, phi=0.0,    duration=T_ROT),
    SQR(name="S1", theta_n=[0.0]*N_CAV, phi_n=[0.0]*N_CAV,
        drift_model=drift, duration=T_SQR),
    QubitRotation(name="R2", theta=np.pi/2, phi=np.pi/2, duration=T_ROT),
    Displacement(name="D2",  alpha=0.3+0j, duration=T_DISP),
    QubitRotation(name="R3", theta=np.pi/2, phi=0.0,    duration=T_ROT),
    SQR(name="S2", theta_n=[np.pi/2]*N_CAV, phi_n=[0.0]*N_CAV,
        drift_model=drift, duration=T_SQR),
    QubitRotation(name="R4", theta=np.pi/2, phi=np.pi/2, duration=T_ROT),
    Displacement(name="D3",  alpha=0.3+0j, duration=T_DISP),
], n_cav=N_CAV)

if "L1b_D_R_SQR_d9" not in all_results and (_args.only is None or _args.only == "L1b_D_R_SQR_d9"):
    all_results["L1b_D_R_SQR_d9"] = run_synth(seq_L1b, "L1b_D_R_SQR_d9")
    with open(DATA_DIR / "phase4_results.json", "w") as _f:
        json.dump(all_results, _f, indent=2, default=str)

# L1c: depth-11  D R SQR R D R SQR R D R SQR
seq_L1c = GateSequence(gates=[
    Displacement(name="D1",  alpha=0.3+0j, duration=T_DISP),
    QubitRotation(name="R1", theta=np.pi/2, phi=0.0,    duration=T_ROT),
    SQR(name="S1", theta_n=[0.0]*N_CAV, phi_n=[0.0]*N_CAV,
        drift_model=drift, duration=T_SQR),
    QubitRotation(name="R2", theta=np.pi/2, phi=np.pi/2, duration=T_ROT),
    Displacement(name="D2",  alpha=0.3+0j, duration=T_DISP),
    QubitRotation(name="R3", theta=np.pi/2, phi=0.0,    duration=T_ROT),
    SQR(name="S2", theta_n=[np.pi/2]*N_CAV, phi_n=[0.0]*N_CAV,
        drift_model=drift, duration=T_SQR),
    QubitRotation(name="R4", theta=np.pi/2, phi=np.pi/2, duration=T_ROT),
    Displacement(name="D3",  alpha=0.3+0j, duration=T_DISP),
    QubitRotation(name="R5", theta=np.pi/4, phi=0.0,    duration=T_ROT),
    SQR(name="S3", theta_n=[np.pi/4]*N_CAV, phi_n=[0.0]*N_CAV,
        drift_model=drift, duration=T_SQR),
], n_cav=N_CAV)

if "L1c_D_R_SQR_d11" not in all_results and (_args.only is None or _args.only == "L1c_D_R_SQR_d11"):
    all_results["L1c_D_R_SQR_d11"] = run_synth(seq_L1c, "L1c_D_R_SQR_d11")
    with open(DATA_DIR / "phase4_results.json", "w") as _f:
        json.dump(all_results, _f, indent=2, default=str)


# ════════════════════════════════════════════════════════════════════════════
# LIBRARY 2: {D(α), QubitRotation, ConditionalPhaseSQR}
#
# D provides Fock mixing; R_q gives qubit rotations;
# CP applies cavity-number-selective phases (no qubit rotation).
# ════════════════════════════════════════════════════════════════════════════
print("\n══════════ Library 2: D + QubitRotation + ConditionalPhaseSQR ══════════")

# L2a: depth-7   D R CP D R CP D
seq_L2a = GateSequence(gates=[
    Displacement(name="D1",  alpha=0.3+0j, duration=T_DISP),
    QubitRotation(name="R1", theta=np.pi/2, phi=0.0, duration=T_ROT),
    ConditionalPhaseSQR(name="CP1", phases_n=[0.0, np.pi] + [0.0]*(N_CAV-2),
                        drift_model=drift, duration=T_CP),
    Displacement(name="D2",  alpha=0.3+0j, duration=T_DISP),
    QubitRotation(name="R2", theta=np.pi/2, phi=0.0, duration=T_ROT),
    ConditionalPhaseSQR(name="CP2", phases_n=[0.0, np.pi/2] + [0.0]*(N_CAV-2),
                        drift_model=drift, duration=T_CP),
    Displacement(name="D3",  alpha=0.3+0j, duration=T_DISP),
], n_cav=N_CAV)

if "L2a_D_R_CP_d7" not in all_results and (_args.only is None or _args.only == "L2a_D_R_CP_d7"):
    all_results["L2a_D_R_CP_d7"] = run_synth(seq_L2a, "L2a_D_R_CP_d7")
    with open(DATA_DIR / "phase4_results.json", "w") as _f:
        json.dump(all_results, _f, indent=2, default=str)

# L2b: depth-9   D R CP R D R CP R D
seq_L2b = GateSequence(gates=[
    Displacement(name="D1",  alpha=0.3+0j, duration=T_DISP),
    QubitRotation(name="R1", theta=np.pi/2, phi=0.0,    duration=T_ROT),
    ConditionalPhaseSQR(name="CP1", phases_n=[0.0, np.pi] + [0.0]*(N_CAV-2),
                        drift_model=drift, duration=T_CP),
    QubitRotation(name="R2", theta=np.pi/2, phi=np.pi/2, duration=T_ROT),
    Displacement(name="D2",  alpha=0.3+0j, duration=T_DISP),
    QubitRotation(name="R3", theta=np.pi/2, phi=0.0,    duration=T_ROT),
    ConditionalPhaseSQR(name="CP2", phases_n=[0.0, np.pi/2] + [0.0]*(N_CAV-2),
                        drift_model=drift, duration=T_CP),
    QubitRotation(name="R4", theta=np.pi/2, phi=np.pi/2, duration=T_ROT),
    Displacement(name="D3",  alpha=0.3+0j, duration=T_DISP),
], n_cav=N_CAV)

if "L2b_D_R_CP_d9" not in all_results and (_args.only is None or _args.only == "L2b_D_R_CP_d9"):
    all_results["L2b_D_R_CP_d9"] = run_synth(seq_L2b, "L2b_D_R_CP_d9")
    with open(DATA_DIR / "phase4_results.json", "w") as _f:
        json.dump(all_results, _f, indent=2, default=str)

# L2c: depth-11  D R CP R D R CP R D R CP
seq_L2c = GateSequence(gates=[
    Displacement(name="D1",  alpha=0.3+0j, duration=T_DISP),
    QubitRotation(name="R1", theta=np.pi/2, phi=0.0,    duration=T_ROT),
    ConditionalPhaseSQR(name="CP1", phases_n=[0.0, np.pi] + [0.0]*(N_CAV-2),
                        drift_model=drift, duration=T_CP),
    QubitRotation(name="R2", theta=np.pi/2, phi=np.pi/2, duration=T_ROT),
    Displacement(name="D2",  alpha=0.3+0j, duration=T_DISP),
    QubitRotation(name="R3", theta=np.pi/2, phi=0.0,    duration=T_ROT),
    ConditionalPhaseSQR(name="CP2", phases_n=[0.0, np.pi/2] + [0.0]*(N_CAV-2),
                        drift_model=drift, duration=T_CP),
    QubitRotation(name="R4", theta=np.pi/2, phi=np.pi/2, duration=T_ROT),
    Displacement(name="D3",  alpha=0.3+0j, duration=T_DISP),
    QubitRotation(name="R5", theta=np.pi/4, phi=0.0,    duration=T_ROT),
    ConditionalPhaseSQR(name="CP3", phases_n=[0.0, np.pi/4] + [0.0]*(N_CAV-2),
                        drift_model=drift, duration=T_CP),
], n_cav=N_CAV)

if "L2c_D_R_CP_d11" not in all_results and (_args.only is None or _args.only == "L2c_D_R_CP_d11"):
    all_results["L2c_D_R_CP_d11"] = run_synth(seq_L2c, "L2c_D_R_CP_d11")
    with open(DATA_DIR / "phase4_results.json", "w") as _f:
        json.dump(all_results, _f, indent=2, default=str)


# ════════════════════════════════════════════════════════════════════════════
# Summary (only printed when all ansatzes are available)
# ════════════════════════════════════════════════════════════════════════════
_ALL_LABELS = ["L1a_D_R_SQR_d7", "L1b_D_R_SQR_d9", "L1c_D_R_SQR_d11",
               "L2a_D_R_CP_d7",  "L2b_D_R_CP_d9",  "L2c_D_R_CP_d11"]
_done = [k for k in _ALL_LABELS if k in all_results]
print(f"\nProgress: {len(_done)}/{len(_ALL_LABELS)} ansatzes complete: {_done}")
if len(_done) < len(_ALL_LABELS):
    print("Run again (or with --only LABEL) to complete remaining ansatzes.")
    sys.exit(0)

print(f"\n\n══════════ Phase 4 Summary ══════════")
header = f"{'Label':<30} {'F_proj':>8} {'L_avg':>7} {'L_worst':>8} {'Depth':>6} {'t(s)':>8}"
print(header)
print("-" * len(header))
for label, rec in all_results.items():
    print(f"{rec['label']:<30} {rec['F_proj']:>8.5f} {rec['leakage_avg']:>7.4f} "
          f"{rec['leakage_worst']:>8.4f} {rec['depth']:>6} {rec['elapsed_s']:>8.1f}")

print("\nBest per library:")
for lib_prefix, lib_name in [("L1", "Library 1 (D+R+SQR)"),
                              ("L2", "Library 2 (D+R+CP)")]:
    lib_recs = {k: v for k, v in all_results.items() if k.startswith(lib_prefix)}
    if lib_recs:
        best = max(lib_recs.items(), key=lambda kv: kv[1]["F_proj"])
        print(f"  {lib_name}: {best[0]}  F={best[1]['F_proj']:.5f}")

with open(DATA_DIR / "phase4_results.json", "w") as f:
    json.dump(all_results, f, indent=2, default=str)
print(f"\nResults written to {DATA_DIR / 'phase4_results.json'}")


# ════════════════════════════════════════════════════════════════════════════
# Figures
# ════════════════════════════════════════════════════════════════════════════
labels  = list(all_results.keys())
F_vals  = [all_results[k]["F_proj"] for k in labels]
L_vals  = [all_results[k]["leakage_avg"] for k in labels]
depths  = [all_results[k]["depth"] for k in labels]
colors  = ["steelblue"  if k.startswith("L1") else "darkorange" for k in labels]

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

ax = axes[0]
bars = ax.bar(range(len(labels)), F_vals, color=colors, alpha=0.8, edgecolor="k")
ax.axhline(0.95, ls=":", color="red", label="F=0.95")
ax.set_xticks(range(len(labels)))
ax.set_xticklabels([all_results[k]["label"] for k in labels],
                   rotation=20, ha="right", fontsize=9)
ax.set_ylabel("$F_{\\rm proj}$", fontsize=12)
ax.set_title("Phase 4: Subspace fidelity by ansatz\n"
             "(blue=Library 1 D+R+SQR, orange=Library 2 D+R+CP)", fontsize=10)
ax.set_ylim(0, 1.05)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3, axis="y")
for bar, v in zip(bars, F_vals):
    ax.text(bar.get_x() + bar.get_width()/2, v + 0.005, f"{v:.3f}",
            ha="center", va="bottom", fontsize=8)

# Mark depth on x-axis
for i, (bar, d) in enumerate(zip(bars, depths)):
    ax.text(bar.get_x() + bar.get_width()/2, 0.01, f"d={d}",
            ha="center", va="bottom", fontsize=7, color="white")

ax = axes[1]
ax.bar(range(len(labels)), L_vals, color=colors, alpha=0.8, edgecolor="k")
ax.set_xticks(range(len(labels)))
ax.set_xticklabels([all_results[k]["label"] for k in labels],
                   rotation=20, ha="right", fontsize=9)
ax.set_ylabel("$L_{\\rm avg}$", fontsize=12)
ax.set_title("Average leakage by ansatz", fontsize=10)
ax.set_ylim(0, max(L_vals) * 1.2 if L_vals else 1)
ax.grid(True, alpha=0.3, axis="y")
for bar, v in zip(ax.patches, L_vals):
    ax.text(bar.get_x() + bar.get_width()/2, v + 0.002, f"{v:.3f}",
            ha="center", va="bottom", fontsize=8)

fig.tight_layout()
fig.savefig(FIG_DIR / "phase4_synthesis_results.png", dpi=150)
fig.savefig(FIG_DIR / "phase4_synthesis_results.pdf")
plt.close(fig)

print(f"Figures saved to   {FIG_DIR}")
print("Phase 4 complete.")
