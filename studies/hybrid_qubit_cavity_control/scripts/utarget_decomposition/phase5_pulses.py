"""Phase 5 — Finite-duration pulse effects.

Compare the best Phase 4 result (L1c: D+QubitRotation+SQR, depth-11, F=0.919)
with and without dispersive drift accumulation during SQR gates.  Quantifies
how much the finite 2-μs SQR pulses contribute to the residual infidelity,
and analyses the coherence budget of the sequence.

Sections
--------
1. Ideal vs Physical:  L1c params with full drift vs no drift during gates
2. SQR duration sweep: dispersive phase accumulated as a function of t_SQR
3. Coherence budget:   total sequence time vs T1, T2
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

STUDY_ROOT = Path(__file__).resolve().parents[2]
COMPONENT_NAME = Path(__file__).resolve().parent.name
SIM_ROOT   = Path("C:/Users/dazzl/Box/Shyam Shankar Quantum Circuits Group"
                  "/Users/Users_JianJun/cQED_simulation")
if str(SIM_ROOT) not in sys.path:
    sys.path.insert(0, str(SIM_ROOT))

from cqed_sim.unitary_synthesis import (
    DriftPhaseModel, Displacement, GateSequence, QubitRotation, SQR, Subspace,
    simulate_sequence, subspace_unitary_fidelity, leakage_metrics,
)

DATA_DIR = STUDY_ROOT / "data" / COMPONENT_NAME
FIG_DIR  = STUDY_ROOT / "figures" / COMPONENT_NAME

# ── Physical parameters ──────────────────────────────────────────────────────
CHI  = 2 * np.pi * (-2.84e6)
CHIP = 2 * np.pi * (-21e3)
KERR = 2 * np.pi * (-28e3)

N_CAV    = 8
FULL_DIM = 2 * N_CAV

# Logical subspace: qubit-first {|g,0>, |g,1>, |e,0>, |e,1>}
LOGICAL_LABELS = ["|g,0>", "|g,1>", "|e,0>", "|e,1>"]
subspace = Subspace.custom(FULL_DIM, [0, 1, N_CAV, N_CAV + 1], LOGICAL_LABELS)

s = 1.0 / np.sqrt(2)
U_target = np.array([
    [ s, 0, s, 0], [ s, 0, -s, 0], [0, s, 0, s], [0, -s, 0, s],
], dtype=np.complex128)

drift_full = DriftPhaseModel(chi=CHI,  chi2=CHIP, kerr=KERR)
no_drift   = DriftPhaseModel(chi=0.0,  chi2=0.0,  kerr=0.0)

results: dict = {}


# ── Load L1c optimised parameters from Phase 4 ──────────────────────────────
with open(DATA_DIR / "phase4_results.json") as f:
    p4 = json.load(f)

l1c_raw = p4["L1c_D_R_SQR_d11"]["sequence_params"]

def _p(name: str) -> list:
    for g in l1c_raw:
        if g["name"] == name:
            return g["parameters"]
    raise KeyError(name)

def _dur(name: str) -> float:
    for g in l1c_raw:
        if g["name"] == name:
            return float(g["duration"])
    raise KeyError(name)


def make_l1c(
    drift_sqr: DriftPhaseModel,
    amp_scale: float = 1.0,
    dur_scale: float = 1.0,
    phi_offset: float = 0.0,
) -> GateSequence:
    """Reconstruct L1c (D+QubitRotation+SQR, depth-11) with perturbations.

    Parameters
    ----------
    drift_sqr  : DriftPhaseModel applied to each SQR gate
    amp_scale  : multiplicative factor on SQR theta_n (amplitude calibration)
    dur_scale  : multiplicative factor on all gate durations
    phi_offset : additive offset (rad) on all SQR phi_n
    """
    def sqr(name: str) -> SQR:
        p = _p(name)
        theta = [v * amp_scale for v in p[:N_CAV]]
        phi   = [v + phi_offset for v in p[N_CAV:]]
        return SQR(
            name=name, theta_n=theta, phi_n=phi,
            drift_model=drift_sqr,
            duration=_dur(name) * dur_scale,
        )

    def rot(name: str) -> QubitRotation:
        p = _p(name)
        return QubitRotation(
            name=name, theta=p[0] * amp_scale, phi=p[1] + phi_offset,
            duration=_dur(name) * dur_scale,
        )

    def disp(name: str) -> Displacement:
        p = _p(name)
        return Displacement(name=name, alpha=complex(p[0], p[1]),
                            duration=_dur(name) * dur_scale)

    # L1c gate order: D1 R1 S1 R2 D2 R3 S2 R4 D3 R5 S3
    return GateSequence(gates=[
        disp("D1"), rot("R1"), sqr("S1"), rot("R2"),
        disp("D2"), rot("R3"), sqr("S2"), rot("R4"),
        disp("D3"), rot("R5"), sqr("S3"),
    ], n_cav=N_CAV)


# ════════════════════════════════════════════════════════════════════════════
# 1. Ideal (no drift during gates) vs Physical (full drift during gates)
# ════════════════════════════════════════════════════════════════════════════
print("\n── 1. Ideal vs Physical drift during gates ──────────────────────────")

seq_physical = make_l1c(drift_sqr=drift_full)
seq_ideal    = make_l1c(drift_sqr=no_drift)

sim_physical = simulate_sequence(seq_physical, subspace=subspace, backend="ideal")
sim_ideal    = simulate_sequence(seq_ideal,    subspace=subspace, backend="ideal")

F_physical = float(subspace_unitary_fidelity(
    sim_physical.subspace_operator, U_target, gauge="global"))
F_ideal    = float(subspace_unitary_fidelity(
    sim_ideal.subspace_operator,    U_target, gauge="global"))

lm_physical = leakage_metrics(sim_physical.full_operator, subspace)
lm_ideal    = leakage_metrics(sim_ideal.full_operator,    subspace)

delta_infid = abs((1 - F_physical) - (1 - F_ideal))

print(f"  F_physical (drift during SQR, confirms Phase 4): {F_physical:.5f}")
print(f"  F_ideal    (drift switched off):                  {F_ideal:.5f}")
print(f"  Infidelity penalty from finite-duration drift:    {delta_infid:.4f}")
print(f"  Leakage physical / ideal: {lm_physical.average:.4f} / {lm_ideal.average:.4f}")

results["ideal_vs_physical"] = {
    "F_physical":   F_physical,
    "F_ideal":      F_ideal,
    "delta_infid":  delta_infid,
    "L_physical":   float(lm_physical.average),
    "L_ideal":      float(lm_ideal.average),
}


# ════════════════════════════════════════════════════════════════════════════
# 2. Dispersive phase accumulated during SQR as a function of pulse duration
# ════════════════════════════════════════════════════════════════════════════
print("\n── 2. SQR duration sweep (analytical) ──────────────────────────────")

t_sqr_vals = np.logspace(np.log10(100e-9), np.log10(20e-6), 50)

phi_n1 = np.abs(CHI * 1 * t_sqr_vals)
phi_n2 = np.abs(CHI * 2 * t_sqr_vals + CHIP * 2 * t_sqr_vals + KERR * t_sqr_vals)
delta_phi_n2 = np.abs((CHIP * 2 + KERR) * t_sqr_vals)  # chi', K excess at n=2
F_approx = np.maximum(0.0, 1.0 - delta_phi_n2**2 / 4)

t_ref = _dur("S1")  # nominal SQR duration from optimizer
delta_phi_ref = float(np.abs((CHIP * 2 + KERR) * t_ref))
F_ref = float(max(0.0, 1.0 - delta_phi_ref**2 / 4))
print(f"  At t_SQR = {t_ref*1e6:.2f} μs: δφ(n=2) = {delta_phi_ref:.3f} rad, F_approx = {F_ref:.4f}")

results["sqr_duration_sweep"] = {
    "t_sqr_us": [float(t * 1e6) for t in t_sqr_vals],
    "phi_n1_rad": phi_n1.tolist(),
    "phi_n2_rad": phi_n2.tolist(),
    "delta_phi_n2_rad": delta_phi_n2.tolist(),
    "F_approx": F_approx.tolist(),
    "t_ref_us": float(t_ref * 1e6),
    "delta_phi_ref_rad": delta_phi_ref,
}


# ════════════════════════════════════════════════════════════════════════════
# 3. Coherence budget
# ════════════════════════════════════════════════════════════════════════════
print("\n── 3. Coherence budget ─────────────────────────────────────────────")

gate_dur = {g["name"]: float(g["duration"]) for g in l1c_raw}
total_time = sum(gate_dur.values())

T1_qubit  = 9.8e-6
T2_qubit  = 6.3e-6
T1_cavity = 200e-6

budget = {
    "total_time_us": float(total_time * 1e6),
    "gate_durations_us": {k: float(v * 1e6) for k, v in gate_dur.items()},
    "qubit_T1_loss":   float(1.0 - np.exp(-total_time / T1_qubit)),
    "qubit_T2_loss":   float(1.0 - np.exp(-total_time / T2_qubit)),
    "cavity_T1_loss":  float(1.0 - np.exp(-total_time / T1_cavity)),
    "T1_qubit_us":     float(T1_qubit  * 1e6),
    "T2_qubit_us":     float(T2_qubit  * 1e6),
    "T1_cavity_us":    float(T1_cavity * 1e6),
}

print(f"  Total sequence time: {budget['total_time_us']:.2f} μs")
for g, t in gate_dur.items():
    print(f"    {g:10s}: {t*1e6:.2f} μs")
print(f"  Qubit T1 loss:   {budget['qubit_T1_loss']*100:.1f}%")
print(f"  Qubit T2 loss:   {budget['qubit_T2_loss']*100:.1f}%")
print(f"  Cavity T1 loss:  {budget['cavity_T1_loss']*100:.2f}%")
results["coherence_budget"] = budget


# ════════════════════════════════════════════════════════════════════════════
# 4. Plots
# ════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

# Panel A: bar chart – ideal vs physical fidelity
ax = axes[0]
labels_bar = ["Physical\n(drift during SQR)", "Ideal\n(no drift during SQR)"]
vals   = [F_physical, F_ideal]
colors = ["steelblue", "forestgreen"]
bars = ax.bar(labels_bar, vals, color=colors, alpha=0.8, edgecolor="k")
ax.axhline(F_physical, ls=":", color="gray", label=f"F₀={F_physical:.4f}")
for bar, v in zip(bars, vals):
    ax.text(bar.get_x() + bar.get_width()/2, v + 0.005, f"{v:.4f}",
            ha="center", va="bottom", fontsize=10)
ax.set_ylabel("$F_{\\rm proj}$", fontsize=11)
ax.set_title("Drift during SQR gates\n(L1c optimised params)", fontsize=10)
ax.set_ylim(0, 1.05)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3, axis="y")

# Panel B: dispersive phase vs SQR duration
ax = axes[1]
ax.semilogx([t * 1e6 for t in t_sqr_vals], phi_n1, "b-",
            lw=2, label=r"$|\phi_{n=1}| = |\chi| t$")
ax.semilogx([t * 1e6 for t in t_sqr_vals], phi_n2, "r--",
            lw=2, label=r"$|\phi_{n=2}|$")
ax.semilogx([t * 1e6 for t in t_sqr_vals], delta_phi_n2, "g:",
            lw=2, label=r"$|\delta\phi_{n=2}|$ ($\chi'$+K excess)")
ax.axvline(t_ref * 1e6, ls="--", color="gray", label=f"Ref: {t_ref*1e6:.2f} μs")
ax.set_xlabel("SQR duration (μs)", fontsize=11)
ax.set_ylabel("Phase (rad)", fontsize=11)
ax.set_title("Dispersive phase during SQR pulse", fontsize=10)
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)

# Panel C: gate duration bar chart with T1 reference
ax = axes[2]
g_names = list(gate_dur.keys())
g_times = [gate_dur[g] * 1e6 for g in g_names]
gate_colors = ["steelblue"  if "D" in g else
               ("darkorange" if "S" in g else "mediumpurple")
               for g in g_names]
ax.barh(g_names, g_times, color=gate_colors, alpha=0.8, edgecolor="k")
ax.axvline(T1_qubit * 1e6, ls="--", color="red",
           label=f"T₁ = {T1_qubit*1e6:.1f} μs")
ax.axvline(total_time * 1e6, ls=":", color="black",
           label=f"Total = {total_time*1e6:.1f} μs")
ax.set_xlabel("Duration (μs)", fontsize=11)
ax.set_title("L1c gate duration breakdown", fontsize=10)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3, axis="x")

fig.tight_layout()
fig.savefig(FIG_DIR / "phase5_pulse_effects.png", dpi=150)
fig.savefig(FIG_DIR / "phase5_pulse_effects.pdf")
plt.close(fig)

with open(DATA_DIR / "phase5_results.json", "w") as f:
    json.dump(results, f, indent=2, default=str)

print(f"\nResults written to {DATA_DIR / 'phase5_results.json'}")
print("Phase 5 complete.")
