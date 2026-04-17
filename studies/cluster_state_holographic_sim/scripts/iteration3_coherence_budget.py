"""Iteration 3 — Experiment 5: Decoherence / coherence budget analysis.

Addresses:
  [P2 | MEDIUM] Add decoherence channels (analytical coherence budget)

Computes decoherence-limited fidelity for all winning strategies based on
their total gate times, using representative T1 and T2 values.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
STUDY_DIR = SCRIPT_DIR.parent
DATA_DIR = STUDY_DIR / "data"
FIG_DIR = STUDY_DIR / "figures"

DATA_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

# ── Physical decoherence parameters ──────────────────────────────────────────
T1_QUBIT = 30e-6      # qubit relaxation (s)
T2_QUBIT = 20e-6      # qubit dephasing (s) — T2 ≤ 2*T1
T1_CAVITY = 200e-6    # cavity relaxation (s)
KAPPA_STORAGE = 2 * np.pi * 10e3  # storage cavity decay rate (rad/s)

TWO_PI = 2 * np.pi
CHI = TWO_PI * (-2.84e6)
TAU_CZ = np.pi / abs(CHI)  # ≈ 176 ns


def decoherence_fidelity(t_total: float, t1: float = T1_QUBIT,
                          t2: float = T2_QUBIT) -> float:
    """Estimate fidelity loss due to relaxation and dephasing.

    Approximate formula for short gates (t_total << T1, T2):
    F_coh ≈ 1 - t_total/(2*T1) - t_total/T2_phi
    where 1/T2_phi = 1/T2 - 1/(2*T1)
    """
    gamma1 = 1 / t1
    gamma2 = 1 / t2
    gamma_phi = gamma2 - gamma1 / 2  # pure dephasing rate
    # Average fidelity loss for a 2-qubit-like gate
    F = np.exp(-(gamma1/2 + gamma_phi) * t_total)
    return float(F)


# ── Strategy gate times ──────────────────────────────────────────────────────
# From Iteration 2 decomposition_comparison.py:
# Strategy B: D+SQR+CP, 2 blocks = [D·SQR·CP·SQR]^2·D
#   5 Displacements (200ns ea) + 4 SQR (400ns ea) + 2 CP (200ns ea)
#   But optimizer can adjust durations — use nominal as upper bound
# Strategy D: D+R+FE, 2 blocks = [D·R·FE]^2·D·R
#   3 Displacements (200ns ea) + 3 Rotations (100ns ea) + 2 FE (~176ns ea)
# GRAPE: single pulse at specified duration

strategies = {
    "B_D+SQR+CP_2blk_ideal": {
        "label": "D+SQR+CP (2 blocks, ideal)",
        "gate_times_ns": {
            "D": [200, 200, 200, 200, 200],
            "SQR": [400, 400, 400, 400],
            "CP": [200, 200],
        },
    },
    "D_D+R+FE_2blk_ideal": {
        "label": "D+R+FE (2 blocks, ideal)",
        "gate_times_ns": {
            "D": [200, 200, 200],
            "R": [100, 100, 100],
            "FE": [176, 176],  # nominal χ-wait ≈ τ_CZ, will update from model-based
        },
    },
    "GRAPE_200ns": {"label": "GRAPE 200 ns", "total_ns": 200},
    "GRAPE_300ns": {"label": "GRAPE 300 ns", "total_ns": 300},
    "GRAPE_400ns": {"label": "GRAPE 400 ns", "total_ns": 400},
    "GRAPE_500ns": {"label": "GRAPE 500 ns", "total_ns": 500},
    "GRAPE_600ns": {"label": "GRAPE 600 ns", "total_ns": 600},
    "GRAPE_800ns": {"label": "GRAPE 800 ns", "total_ns": 800},
}

# Compute total gate times
for key, val in strategies.items():
    if "total_ns" not in val:
        total = 0
        for gate_type, times in val["gate_times_ns"].items():
            total += sum(times)
        val["total_ns"] = total

# ── Load model-based results for updated FE times if available ────────────────
try:
    mb_path = DATA_DIR / "iteration3_model_based.json"
    if mb_path.exists():
        mb_data = json.loads(mb_path.read_text(encoding="utf-8"))
        fe_analysis = mb_data.get("fe_wait_time_analysis", {})
        if "gates" in fe_analysis:
            fe_dur_ns = [g["duration_ns"] for g in fe_analysis["gates"]]
            strategies["D_D+R+FE_2blk_model"] = {
                "label": "D+R+FE (2 blocks, model)",
                "gate_times_ns": {
                    "D": [200, 200, 200],
                    "R": [100, 100, 100],
                    "FE": fe_dur_ns,
                },
            }
            total = sum([200]*3 + [100]*3 + fe_dur_ns)
            strategies["D_D+R+FE_2blk_model"]["total_ns"] = total
            print(f"Loaded model-based FE times: {fe_dur_ns}")
except Exception as e:
    print(f"Could not load model-based results: {e}")

# ── Prior fidelity results ────────────────────────────────────────────────────
ideal_fidelities = {
    "B_D+SQR+CP_2blk_ideal": 1.0000,
    "D_D+R+FE_2blk_ideal": 0.9999,
    "GRAPE_200ns": 0.9966,
    "GRAPE_300ns": 0.9957,
    "GRAPE_400ns": 0.9990,
}

# ═══════════════════════════════════════════════════════════════════════════════
# Compute coherence budget
# ═══════════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("  Iteration 3: Coherence Budget Analysis")
print("=" * 70)
print(f"  T1_qubit  = {T1_QUBIT*1e6:.0f} μs")
print(f"  T2_qubit  = {T2_QUBIT*1e6:.0f} μs")
print(f"  T1_cavity = {T1_CAVITY*1e6:.0f} μs")
print(f"  τ_CZ      = {TAU_CZ*1e9:.1f} ns")

budget = {}

for key, val in sorted(strategies.items()):
    t_total_ns = val["total_ns"]
    t_total_s = t_total_ns * 1e-9
    F_coh = decoherence_fidelity(t_total_s)
    F_ideal = ideal_fidelities.get(key)
    # Combined fidelity (multiplicative approximation)
    F_combined = F_coh * F_ideal if F_ideal is not None else None

    entry = {
        "label": val["label"],
        "total_time_ns": t_total_ns,
        "total_time_us": t_total_ns / 1000,
        "F_coherence": F_coh,
        "F_ideal": F_ideal,
        "F_combined": F_combined,
    }
    budget[key] = entry
    comb_str = f"  F_combined = {F_combined:.6f}" if F_combined else ""
    print(f"  {val['label']:35s}  t = {t_total_ns:6.0f} ns  "
          f"F_coh = {F_coh:.6f}{comb_str}")

# ═══════════════════════════════════════════════════════════════════════════════
# Figure: Coherence-limited fidelity vs gate time
# ═══════════════════════════════════════════════════════════════════════════════
print("\n── Generating coherence budget figure ──")

t_range_ns = np.linspace(50, 5000, 500)
F_coh_curve = np.array([decoherence_fidelity(t*1e-9) for t in t_range_ns])

fig, ax = plt.subplots(figsize=(7, 5))

# T1/T2 limit curve
ax.plot(t_range_ns, F_coh_curve, 'k-', lw=2, label=f'T1={T1_QUBIT*1e6:.0f} μs, T2={T2_QUBIT*1e6:.0f} μs')

# Plot strategies
strategy_colors = {
    "B_D+SQR+CP": '#EE6677',
    "D_D+R+FE": '#228833',
    "GRAPE": '#4477AA',
}
markers = {"B_D+SQR+CP": 's', "D_D+R+FE": 'D', "GRAPE": 'o'}

for key, val in sorted(budget.items()):
    color = '#999999'
    marker = 'o'
    for prefix, c in strategy_colors.items():
        if key.startswith(prefix):
            color = c
            marker = markers[prefix]
            break
    t_ns = val["total_time_ns"]
    # Plot coherence limit
    ax.plot(t_ns, val["F_coherence"], marker=marker, color=color, ms=10,
            markeredgecolor='black', markeredgewidth=0.5, zorder=5)
    # Plot combined fidelity if known
    if val["F_combined"] is not None:
        ax.plot(t_ns, val["F_combined"], marker=marker, color=color, ms=7,
                markeredgecolor='black', markeredgewidth=0.5, alpha=0.5, zorder=4)
        ax.annotate('', xy=(t_ns, val["F_combined"]),
                     xytext=(t_ns, val["F_coherence"]),
                     arrowprops=dict(arrowstyle='->', color=color, lw=1, alpha=0.5))

# Custom legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], color='k', lw=2, label='Decoherence limit'),
    Line2D([0], [0], marker='s', color='#EE6677', ms=10, lw=0,
           markeredgecolor='k', markeredgewidth=0.5, label='D+SQR+CP'),
    Line2D([0], [0], marker='D', color='#228833', ms=10, lw=0,
           markeredgecolor='k', markeredgewidth=0.5, label='D+R+FE'),
    Line2D([0], [0], marker='o', color='#4477AA', ms=10, lw=0,
           markeredgecolor='k', markeredgewidth=0.5, label='GRAPE'),
    Line2D([0], [0], marker='o', color='gray', ms=7, lw=0, alpha=0.5,
           markeredgecolor='k', markeredgewidth=0.5, label='Combined (ideal × coh)'),
]
ax.legend(handles=legend_elements, fontsize=9, loc='lower left')

ax.axhline(0.99, ls='--', color='red', alpha=0.4, label='99%')
ax.axhline(0.999, ls=':', color='darkred', alpha=0.4, label='99.9%')
ax.set_xlabel('Total gate time (ns)', fontsize=11)
ax.set_ylabel('Fidelity', fontsize=11)
ax.set_title('Coherence Budget: Gate Strategy Comparison', fontsize=12)
ax.set_ylim(0.96, 1.002)
ax.set_xlim(0, 3500)
ax.grid(True, alpha=0.3)
fig.tight_layout()

for fmt in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"coherence_budget.{fmt}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("  coherence_budget.{png,pdf} saved")

# ═══════════════════════════════════════════════════════════════════════════════
# Figure: T1/T2 sensitivity — how F_combined depends on T1, T2
# ═══════════════════════════════════════════════════════════════════════════════
print("\n── Generating T1/T2 sensitivity figure ──")

T1_range = np.linspace(10e-6, 100e-6, 50)
T2_range = np.linspace(5e-6, 100e-6, 50)

# For the three main strategies, plot F_combined vs T2 at T1=30μs
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Left: F_combined vs T2 for fixed T1
key_strategies = [
    ("B_D+SQR+CP_2blk_ideal", '#EE6677', 's'),
    ("D_D+R+FE_2blk_ideal", '#228833', 'D'),
    ("GRAPE_400ns", '#4477AA', 'o'),
    ("GRAPE_200ns", '#66CCEE', '^'),
]

ax = axes[0]
for key, color, marker in key_strategies:
    if key not in budget:
        continue
    t_s = budget[key]["total_time_ns"] * 1e-9
    F_ideal = budget[key]["F_ideal"]
    if F_ideal is None:
        continue
    F_vs_T2 = [F_ideal * decoherence_fidelity(t_s, T1_QUBIT, t2) for t2 in T2_range]
    ax.plot(T2_range * 1e6, F_vs_T2, color=color, lw=2, marker=marker,
            markevery=10, ms=6, label=budget[key]["label"])

ax.axhline(0.99, ls='--', color='red', alpha=0.4)
ax.axhline(0.999, ls=':', color='darkred', alpha=0.4)
ax.set_xlabel('T2 (μs)', fontsize=11)
ax.set_ylabel('Combined fidelity', fontsize=11)
ax.set_title(f'Fidelity vs T2 (T1 = {T1_QUBIT*1e6:.0f} μs)', fontsize=12)
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_ylim(0.95, 1.002)

# Right: F_combined vs T1 for fixed T2
ax = axes[1]
for key, color, marker in key_strategies:
    if key not in budget:
        continue
    t_s = budget[key]["total_time_ns"] * 1e-9
    F_ideal = budget[key]["F_ideal"]
    if F_ideal is None:
        continue
    F_vs_T1 = [F_ideal * decoherence_fidelity(t_s, t1, T2_QUBIT) for t1 in T1_range]
    ax.plot(T1_range * 1e6, F_vs_T1, color=color, lw=2, marker=marker,
            markevery=10, ms=6, label=budget[key]["label"])

ax.axhline(0.99, ls='--', color='red', alpha=0.4)
ax.axhline(0.999, ls=':', color='darkred', alpha=0.4)
ax.set_xlabel('T1 (μs)', fontsize=11)
ax.set_ylabel('Combined fidelity', fontsize=11)
ax.set_title(f'Fidelity vs T1 (T2 = {T2_QUBIT*1e6:.0f} μs)', fontsize=12)
ax.legend(fontsize=8)
ax.grid(True, alpha=0.3)
ax.set_ylim(0.95, 1.002)

fig.tight_layout()
for fmt in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"decoherence_sensitivity.{fmt}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("  decoherence_sensitivity.{png,pdf} saved")

# ═══════════════════════════════════════════════════════════════════════════════
# Save results
# ═══════════════════════════════════════════════════════════════════════════════
out = {
    "parameters": {
        "T1_qubit_us": T1_QUBIT * 1e6,
        "T2_qubit_us": T2_QUBIT * 1e6,
        "T1_cavity_us": T1_CAVITY * 1e6,
        "tau_cz_ns": TAU_CZ * 1e9,
    },
    "budget": budget,
}
out_path = DATA_DIR / "iteration3_coherence_budget.json"
out_path.write_text(json.dumps(out, indent=2, default=str), encoding="utf-8")
print(f"\nResults saved to {out_path}")
print("Done.")
