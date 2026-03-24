"""Phase 6 — Robustness analysis.

Evaluate sensitivity of the best Phase 4 sequence (L1c: D+QubitRotation+SQR,
depth-11, F=0.919) to systematic calibration errors:

  1. χ calibration error   (±10%)
  2. Pulse amplitude error (±10%)
  3. Pulse duration error  (±15%)
  4. SQR phase offset      (±0.2 rad)

For each perturbation, the optimised L1c parameters from Phase 4 are held
fixed and fidelity is re-evaluated without re-optimising.
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

# ── Nominal physical parameters ──────────────────────────────────────────────
CHI_NOM  = 2 * np.pi * (-2.84e6)
CHIP_NOM = 2 * np.pi * (-21e3)
KERR_NOM = 2 * np.pi * (-28e3)

N_CAV_NOM = 8

s = 1.0 / np.sqrt(2)
U_target = np.array([
    [ s, 0, s, 0], [ s, 0, -s, 0], [0, s, 0, s], [0, -s, 0, s],
], dtype=np.complex128)

LOGICAL_LABELS = ["|g,0>", "|g,1>", "|e,0>", "|e,1>"]

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
    chi: float,
    chip: float,
    kerr: float,
    amp_scale: float = 1.0,
    dur_scale: float = 1.0,
    phi_offset: float = 0.0,
) -> GateSequence:
    """Build L1c with perturbed parameters."""
    drift = DriftPhaseModel(chi=chi, chi2=chip, kerr=kerr)

    def sqr(name: str) -> SQR:
        p = _p(name)
        theta = [v * amp_scale for v in p[:N_CAV_NOM]]
        phi   = [v + phi_offset for v in p[N_CAV_NOM:]]
        return SQR(name=name, theta_n=theta, phi_n=phi,
                   drift_model=drift, duration=_dur(name) * dur_scale)

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
    ], n_cav=N_CAV_NOM)


def eval_seq(seq: GateSequence) -> dict:
    sub = Subspace.custom(2 * N_CAV_NOM, [0, 1, N_CAV_NOM, N_CAV_NOM + 1], LOGICAL_LABELS)
    sim = simulate_sequence(seq, subspace=sub, backend="ideal")
    F = float(subspace_unitary_fidelity(sim.subspace_operator, U_target, gauge="global"))
    lm = leakage_metrics(sim.full_operator, sub)
    return {"F_proj": F, "leakage": float(lm.average)}


# ── Nominal (confirms Phase 4) ───────────────────────────────────────────────
nom_seq = make_l1c(CHI_NOM, CHIP_NOM, KERR_NOM)
nom_res = eval_seq(nom_seq)
F_nom   = nom_res["F_proj"]
print(f"\nNominal L1c (should match Phase 4 ≈ 0.9193): "
      f"F={F_nom:.5f}  L={nom_res['leakage']:.4f}")
results["nominal"] = nom_res


# ════════════════════════════════════════════════════════════════════════════
# 1. χ calibration error  (±10%)
# ════════════════════════════════════════════════════════════════════════════
print("\n── 1. χ calibration error ───────────────────────────────────────────")
chi_eps = np.linspace(-0.10, 0.10, 21)
chi_sweep: dict[float, dict] = {}
for eps in chi_eps:
    seq = make_l1c(CHI_NOM * (1 + eps), CHIP_NOM, KERR_NOM)
    chi_sweep[float(eps)] = eval_seq(seq)

results["chi_error_sweep"] = chi_sweep
F_chi = [chi_sweep[e]["F_proj"] for e in sorted(chi_sweep)]
print(f"  F at ε=-5%: {chi_sweep[min(chi_sweep, key=lambda e: abs(e+0.05))]['F_proj']:.5f}")
print(f"  F at ε= 0%: {chi_sweep[0.0]['F_proj']:.5f}")
print(f"  F at ε=+5%: {chi_sweep[min(chi_sweep, key=lambda e: abs(e-0.05))]['F_proj']:.5f}")


# ════════════════════════════════════════════════════════════════════════════
# 2. Pulse amplitude error (±10%)
# ════════════════════════════════════════════════════════════════════════════
print("\n── 2. Pulse amplitude error ─────────────────────────────────────────")
amp_eps = np.linspace(-0.10, 0.10, 21)
amp_sweep: dict[float, dict] = {}
for eps in amp_eps:
    seq = make_l1c(CHI_NOM, CHIP_NOM, KERR_NOM, amp_scale=1 + eps)
    amp_sweep[float(eps)] = eval_seq(seq)

results["amp_error_sweep"] = amp_sweep
F_amp = [amp_sweep[e]["F_proj"] for e in sorted(amp_sweep)]
print(f"  F at ε=-5%: {amp_sweep[min(amp_sweep, key=lambda e: abs(e+0.05))]['F_proj']:.5f}")
print(f"  F at ε= 0%: {amp_sweep[0.0]['F_proj']:.5f}")
print(f"  F at ε=+5%: {amp_sweep[min(amp_sweep, key=lambda e: abs(e-0.05))]['F_proj']:.5f}")


# ════════════════════════════════════════════════════════════════════════════
# 3. Pulse duration error (±15%)
# ════════════════════════════════════════════════════════════════════════════
print("\n── 3. Pulse duration error ──────────────────────────────────────────")
dur_eps = np.linspace(-0.15, 0.15, 21)
dur_sweep: dict[float, dict] = {}
for eps in dur_eps:
    seq = make_l1c(CHI_NOM, CHIP_NOM, KERR_NOM, dur_scale=1 + eps)
    dur_sweep[float(eps)] = eval_seq(seq)

results["dur_error_sweep"] = dur_sweep
F_dur = [dur_sweep[e]["F_proj"] for e in sorted(dur_sweep)]
print(f"  F at ε=-10%: {dur_sweep[min(dur_sweep, key=lambda e: abs(e+0.10))]['F_proj']:.5f}")
print(f"  F at ε=  0%: {dur_sweep[0.0]['F_proj']:.5f}")
print(f"  F at ε=+10%: {dur_sweep[min(dur_sweep, key=lambda e: abs(e-0.10))]['F_proj']:.5f}")


# ════════════════════════════════════════════════════════════════════════════
# 4. SQR phase offset (±0.2 rad)
# ════════════════════════════════════════════════════════════════════════════
print("\n── 4. SQR phase offset ──────────────────────────────────────────────")
phi_eps = np.linspace(-0.20, 0.20, 21)
phi_sweep: dict[float, float] = {}
for d_phi in phi_eps:
    seq = make_l1c(CHI_NOM, CHIP_NOM, KERR_NOM, phi_offset=d_phi)
    res = eval_seq(seq)
    phi_sweep[float(d_phi)] = res["F_proj"]

results["phi_offset_sweep"] = phi_sweep
print(f"  F at δφ=-0.1 rad: {phi_sweep[min(phi_sweep, key=lambda e: abs(e+0.10))]:.5f}")
print(f"  F at δφ=0:        {phi_sweep[min(phi_sweep, key=lambda e: abs(e))]:.5f}")
print(f"  F at δφ=+0.1 rad: {phi_sweep[min(phi_sweep, key=lambda e: abs(e-0.10))]:.5f}")


# ════════════════════════════════════════════════════════════════════════════
# 5. Sensitivity summary (RMS fidelity drop over ±5%)
# ════════════════════════════════════════════════════════════════════════════
def rms_drop(sweep: dict, band: float = 0.05) -> float:
    vals = [
        (v["F_proj"] if isinstance(v, dict) else v) - F_nom
        for e, v in sweep.items()
        if abs(e) <= band
    ]
    return float(np.sqrt(np.mean(np.array(vals) ** 2))) if vals else float("nan")

sens = {
    "chi_rms_drop":  rms_drop(chi_sweep),
    "amp_rms_drop":  rms_drop(amp_sweep),
    "dur_rms_drop":  rms_drop(dur_sweep),
    "phi_rms_drop":  rms_drop({e: {"F_proj": f} for e, f in phi_sweep.items()}, band=0.10),
}
print(f"\nRobustness (RMS fidelity drop over ±5% / ±0.1 rad for φ):")
for k, v in sens.items():
    print(f"  {k}: {v:.4f}")
results["sensitivity_summary"] = sens


# ════════════════════════════════════════════════════════════════════════════
# 6. Plots
# ════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(2, 2, figsize=(12, 9))

eps_vals = sorted(chi_sweep.keys())

# χ error
ax = axes[0, 0]
ax.plot([e * 100 for e in eps_vals], F_chi, "o-", color="steelblue", lw=2, markersize=5)
ax.axvline(0, ls="--", color="gray")
ax.axhline(F_nom, ls=":", color="gray", label=f"F₀={F_nom:.4f}")
ax.set_xlabel("χ error (%)", fontsize=11)
ax.set_ylabel("$F_{\\rm proj}$", fontsize=11)
ax.set_title("χ calibration error", fontsize=11)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_ylim(max(0, F_nom - 0.05), min(1.02, F_nom + 0.03))

# Amplitude error
ax = axes[0, 1]
ax.plot([e * 100 for e in sorted(amp_sweep)], F_amp, "s-", color="darkorange", lw=2, markersize=5)
ax.axvline(0, ls="--", color="gray")
ax.axhline(F_nom, ls=":", color="gray", label=f"F₀={F_nom:.4f}")
ax.set_xlabel("Amplitude error (%)", fontsize=11)
ax.set_ylabel("$F_{\\rm proj}$", fontsize=11)
ax.set_title("Amplitude calibration error", fontsize=11)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_ylim(max(0, F_nom - 0.05), min(1.02, F_nom + 0.03))

# Duration error
ax = axes[1, 0]
ax.plot([e * 100 for e in sorted(dur_sweep)], F_dur, "^-", color="forestgreen", lw=2, markersize=5)
ax.axvline(0, ls="--", color="gray")
ax.axhline(F_nom, ls=":", color="gray", label=f"F₀={F_nom:.4f}")
ax.set_xlabel("Duration error (%)", fontsize=11)
ax.set_ylabel("$F_{\\rm proj}$", fontsize=11)
ax.set_title("Duration calibration error", fontsize=11)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_ylim(max(0, F_nom - 0.05), min(1.02, F_nom + 0.03))

# Phase offset
ax = axes[1, 1]
ax.plot(sorted(phi_sweep.keys()),
        [phi_sweep[e] for e in sorted(phi_sweep.keys())],
        "D-", color="brown", lw=2, markersize=5)
ax.axvline(0, ls="--", color="gray")
ax.axhline(F_nom, ls=":", color="gray", label=f"F₀={F_nom:.4f}")
ax.set_xlabel("Phase offset δφ (rad)", fontsize=11)
ax.set_ylabel("$F_{\\rm proj}$", fontsize=11)
ax.set_title("SQR phase offset", fontsize=11)
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
ax.set_ylim(max(0, F_nom - 0.05), min(1.02, F_nom + 0.03))

fig.suptitle("Phase 6: Robustness of L1c (D+QubitRotation+SQR) optimised sequence", fontsize=13)
fig.tight_layout()
fig.savefig(FIG_DIR / "phase6_robustness.png", dpi=150)
fig.savefig(FIG_DIR / "phase6_robustness.pdf")
plt.close(fig)

with open(DATA_DIR / "phase6_results.json", "w") as f:
    json.dump(results, f, indent=2, default=str)

print(f"\nResults written to {DATA_DIR / 'phase6_results.json'}")
print("Phase 6 complete.")
