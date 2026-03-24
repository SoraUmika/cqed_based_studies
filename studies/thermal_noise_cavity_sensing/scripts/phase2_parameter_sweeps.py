"""
Phase 2: Parameter sweeps for cavity-only sensing feasibility.

Studies how the post-cavity steady-state photon number responds to the target bath
parameters and coupling fractions.  Identifies the sensing regime and quantifies
the sensitivity.

Sweeps:
  - n_target ∈ [0, 10] at fixed κ_target/κ_tot = 0.5
  - κ_target/κ_tot ∈ [0, 1] at fixed n_target = 2.0
  - 2D sensitivity heatmap in (n_target, κ_frac) space
  - Photon-number distributions P_n for representative conditions
  - Heating/cooling transients

Scientific questions answered:
  Q1. How strongly must the target device couple to the cavity for it to be
      distinguishable from the background floor?
  Q2. Is steady-state n̄ alone enough to infer n_target?
  Q3. What parameter combinations are degenerate at the cavity-only level?

Usage:
    python scripts/phase2_parameter_sweeps.py

Output:
    data/phase2_results.npz
    figures/phase2_nss_vs_ntarget.{png,pdf}
    figures/phase2_nss_vs_kappafrac.{png,pdf}
    figures/phase2_sensitivity_heatmap.{png,pdf}
    figures/phase2_transients.{png,pdf}
    figures/phase2_pn_conditions.{png,pdf}
"""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import runtime_compat  # noqa: F401 — must be before any qutip/cqed_sim import

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import qutip as qt

from common import (
    DATA_DIR, FIG_DIR, STYLE_PATH,
    KAPPA_TOT, N_CAV,
    ThermalBath, analytic_nbar_ss, analytic_kappa_tot,
    analytic_nbar_transient, thermal_pn,
    build_cavity_c_ops, check_truncation, default_baths,
    OMEGA_C, n_thermal, temperature_from_nbar,
)

if STYLE_PATH.exists():
    plt.style.use(STYLE_PATH)

COLORS = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377"]

print("=" * 60)
print("Phase 2: Parameter Sweeps — Sensing Feasibility")
print("=" * 60)

# ---------------------------------------------------------------------------
# Fixed background parameters
# ---------------------------------------------------------------------------

N_BG = 0.01                    # cold background bath occupation
N_INT = 0.00                   # internal loss (zero T)
KAPPA_FRAC_BG = 0.30
KAPPA_FRAC_INT = 0.20
KAPPA_FRAC_TARGET_DEFAULT = 0.50

# Background floor (signal when n_target = 0)
def n_floor(kappa_frac_target: float) -> float:
    """Steady-state n̄ when n_target = 0 (background only)."""
    baths_bg = [
        ThermalBath(kappa=(1.0 - kappa_frac_target) * KAPPA_FRAC_BG / (KAPPA_FRAC_BG + KAPPA_FRAC_INT) * KAPPA_TOT,
                    n_th=N_BG, label="bg"),
        ThermalBath(kappa=(1.0 - kappa_frac_target) * KAPPA_FRAC_INT / (KAPPA_FRAC_BG + KAPPA_FRAC_INT) * KAPPA_TOT,
                    n_th=N_INT, label="int"),
    ]
    return analytic_nbar_ss(baths_bg)


def baths_from_params(n_target: float, kappa_frac_target: float) -> list:
    """Build three-bath list from sensing study parameters."""
    kappa_residual = 1.0 - kappa_frac_target
    # Split residual between bg and int proportionally
    frac_bg = KAPPA_FRAC_BG / (KAPPA_FRAC_BG + KAPPA_FRAC_INT)
    frac_int = KAPPA_FRAC_INT / (KAPPA_FRAC_BG + KAPPA_FRAC_INT)
    return [
        ThermalBath(kappa=kappa_frac_target * KAPPA_TOT,          n_th=n_target, label="target"),
        ThermalBath(kappa=kappa_residual * frac_bg * KAPPA_TOT,   n_th=N_BG,     label="bg"),
        ThermalBath(kappa=kappa_residual * frac_int * KAPPA_TOT,  n_th=N_INT,    label="int"),
    ]


# ---------------------------------------------------------------------------
# Sweep 1: n̄_ss vs n_target
# ---------------------------------------------------------------------------
print("\n--- Sweep 1: n̄_ss vs n_target ---")

n_target_arr = np.linspace(0, 10, 60)
kf_vals = [0.1, 0.3, 0.5, 0.7, 0.9]

nss_vs_ntarget = np.zeros((len(kf_vals), len(n_target_arr)))

for i, kf in enumerate(kf_vals):
    for j, n_t in enumerate(n_target_arr):
        baths = baths_from_params(n_t, kf)
        nss_vs_ntarget[i, j] = analytic_nbar_ss(baths)

# Also compute analytic background floor for each kappa_frac_target
n_floor_arr = np.array([
    analytic_nbar_ss(baths_from_params(0.0, kf)) for kf in kf_vals
])

print(f"  Done. n_floor at kf=0.5: {n_floor_arr[2]:.4f}")

# ---------------------------------------------------------------------------
# Sweep 2: n̄_ss vs κ_target/κ_tot
# ---------------------------------------------------------------------------
print("--- Sweep 2: n̄_ss vs κ_target/κ_tot ---")

kappa_frac_arr = np.linspace(0, 1, 60)
n_target_fixed_vals = [0.0, 0.5, 1.0, 2.0, 5.0]

nss_vs_kappafrac = np.zeros((len(n_target_fixed_vals), len(kappa_frac_arr)))

for i, n_t in enumerate(n_target_fixed_vals):
    for j, kf in enumerate(kappa_frac_arr):
        baths = baths_from_params(n_t, kf)
        nss_vs_kappafrac[i, j] = analytic_nbar_ss(baths)

print(f"  Done.")

# ---------------------------------------------------------------------------
# Sweep 3: 2D sensitivity heatmap
# ---------------------------------------------------------------------------
print("--- Sweep 3: 2D sensitivity heatmap ---")

N_TARGET_GRID = np.linspace(0.01, 8, 40)
KAPPA_FRAC_GRID = np.linspace(0.02, 0.98, 40)
NN, KK = np.meshgrid(N_TARGET_GRID, KAPPA_FRAC_GRID)

# Signal = n̄_ss - n̄_floor (shift above background)
SIGNAL = np.zeros_like(NN)
for i in range(len(KAPPA_FRAC_GRID)):
    for j in range(len(N_TARGET_GRID)):
        kf = KAPPA_FRAC_GRID[i]
        n_t = N_TARGET_GRID[j]
        baths = baths_from_params(n_t, kf)
        n_bg_here = analytic_nbar_ss(baths_from_params(0.0, kf))
        SIGNAL[i, j] = analytic_nbar_ss(baths) - n_bg_here

print(f"  Done. Max signal: {SIGNAL.max():.3f}")

# Minimum detectable signal threshold (measurement resolution ~0.05 photons,
# corresponds to averaging ~400 shots with single-shot noise ~1 photon)
MIN_DETECTABLE = 0.05

# ---------------------------------------------------------------------------
# Sweep 4: Transient heating / cooling curves
# ---------------------------------------------------------------------------
print("--- Sweep 4: Transients ---")

# Heating: vacuum to thermal state
n_trans_vals = [0.5, 1.0, 2.0, 5.0]
kappa_tot = KAPPA_TOT
tau_max = 5.0 / kappa_tot
tlist = np.linspace(0, tau_max, 200)

transient_heat = {}
transient_cool = {}

for n_t in n_trans_vals:
    baths = baths_from_params(n_t, KAPPA_FRAC_TARGET_DEFAULT)
    n_ss = analytic_nbar_ss(baths)
    transient_heat[n_t] = analytic_nbar_transient(tlist, 0.0, baths)    # heat from 0
    transient_cool[n_t] = analytic_nbar_transient(tlist, n_ss * 2, baths)  # cool from hot

print(f"  Done.")

# ---------------------------------------------------------------------------
# Sweep 5: P_n distributions for representative conditions
# ---------------------------------------------------------------------------
print("--- Sweep 5: P_n distributions ---")

N = N_CAV
a = qt.destroy(N)
H = qt.qzero(N)
n_op = a.dag() * a

pn_conditions = [
    ("Cold cavity", 0.0, 0.5),   # n_target=0, kf=0.5 → near vacuum
    ("n_t=1, 50% coupling", 1.0, 0.5),
    ("n_t=2, 50% coupling", 2.0, 0.5),
    ("n_t=5, 50% coupling", 5.0, 0.5),
    ("n_t=2, 20% coupling", 2.0, 0.2),
]

pn_data = {}
for label, n_t, kf in pn_conditions:
    check_truncation(analytic_nbar_ss(baths_from_params(n_t, kf)), N)
    baths = baths_from_params(n_t, kf)
    c_ops = build_cavity_c_ops(a, baths)
    rho_ss = qt.steadystate(H, c_ops)
    pn_num = np.array([rho_ss[n, n].real for n in range(N)])
    n_ss = analytic_nbar_ss(baths)
    pn_data[label] = {"pn": pn_num, "n_ss": n_ss, "n_target": n_t, "kf": kf}
    print(f"  {label}: n̄_ss={n_ss:.4f}")

# Also compute equivalent temperatures
print("\n  Temperature interpretations (ω_c/2π = 5.241 GHz):")
for label, dat in pn_data.items():
    n_ss = dat["n_ss"]
    T_eff = temperature_from_nbar(OMEGA_C, max(n_ss, 1e-10))
    print(f"  {label}: n̄_ss={n_ss:.4f} → T_eff={T_eff*1e3:.1f} mK")

# ---------------------------------------------------------------------------
# Generate figures
# ---------------------------------------------------------------------------
print("\nGenerating figures...")

# Figure 1: n̄_ss vs n_target
fig, ax = plt.subplots(figsize=(3.375, 2.8))

for i, kf in enumerate(kf_vals):
    ax.plot(n_target_arr, nss_vs_ntarget[i], color=COLORS[i], lw=1.2,
            label=rf"$\kappa_t/\kappa_\mathrm{{tot}}={kf}$")

ax.axhline(n_floor_arr[2], color="gray", ls="--", lw=0.8,
           label=rf"Floor ($\kappa_t/\kappa_\mathrm{{tot}}=0.5$)")
ax.set_xlabel(r"Target bath occupation $n_\mathrm{target}$")
ax.set_ylabel(r"Cavity steady state $\bar{n}_\mathrm{ss}$")
ax.set_title(r"$\bar{n}_\mathrm{ss}$ vs target occupation")
ax.legend(fontsize=7, loc="upper left")
ax.set_xlim(0, 10)
ax.set_ylim(bottom=0)
plt.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"phase2_nss_vs_ntarget.{ext}")
plt.close(fig)
print("  Saved: phase2_nss_vs_ntarget")

# Figure 2: n̄_ss vs κ_target/κ_tot
fig, ax = plt.subplots(figsize=(3.375, 2.8))

for i, n_t in enumerate(n_target_fixed_vals):
    ax.plot(kappa_frac_arr, nss_vs_kappafrac[i], color=COLORS[i], lw=1.2,
            label=rf"$n_\mathrm{{target}}={n_t}$")

ax.set_xlabel(r"$\kappa_\mathrm{target}/\kappa_\mathrm{tot}$")
ax.set_ylabel(r"$\bar{n}_\mathrm{ss}$")
ax.set_title(r"$\bar{n}_\mathrm{ss}$ vs coupling fraction")
ax.legend(fontsize=7, loc="upper left")
ax.set_xlim(0, 1)
ax.set_ylim(bottom=0)
plt.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"phase2_nss_vs_kappafrac.{ext}")
plt.close(fig)
print("  Saved: phase2_nss_vs_kappafrac")

# Figure 3: 2D sensitivity heatmap
fig, ax = plt.subplots(figsize=(3.375, 2.8))

im = ax.pcolormesh(N_TARGET_GRID, KAPPA_FRAC_GRID, SIGNAL,
                   cmap="viridis", shading="auto")
cb = plt.colorbar(im, ax=ax)
cb.set_label(r"$\bar{n}_\mathrm{ss} - \bar{n}_\mathrm{floor}$", fontsize=8)
cb.ax.tick_params(labelsize=7)

# Mark minimum-detectable contour
cs = ax.contour(N_TARGET_GRID, KAPPA_FRAC_GRID, SIGNAL,
                levels=[MIN_DETECTABLE], colors=["white"], linewidths=1.0)
ax.clabel(cs, fmt=f"{MIN_DETECTABLE:.2f}", fontsize=7, colors="white")

ax.set_xlabel(r"$n_\mathrm{target}$")
ax.set_ylabel(r"$\kappa_\mathrm{target}/\kappa_\mathrm{tot}$")
ax.set_title("Sensing sensitivity map")
plt.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"phase2_sensitivity_heatmap.{ext}")
plt.close(fig)
print("  Saved: phase2_sensitivity_heatmap")

# Figure 4: Transients
fig, axes = plt.subplots(1, 2, figsize=(6.75, 2.5))
t_us = tlist * 1e6

ax = axes[0]
for i, n_t in enumerate(n_trans_vals):
    ax.plot(t_us, transient_heat[n_t], color=COLORS[i], lw=1.2,
            label=rf"$n_t={n_t}$")
ax.set_xlabel(r"$t$ ($\mu$s)")
ax.set_ylabel(r"$\bar{n}(t)$")
ax.set_title(r"Heating: $\bar{n}(0)=0$")
ax.legend(fontsize=7)
ax.set_ylim(bottom=0)

ax = axes[1]
for i, n_t in enumerate(n_trans_vals):
    n_ss_v = analytic_nbar_ss(baths_from_params(n_t, KAPPA_FRAC_TARGET_DEFAULT))
    ax.plot(t_us, transient_cool[n_t], color=COLORS[i], lw=1.2,
            label=rf"$n_t={n_t}$")
    ax.axhline(n_ss_v, color=COLORS[i], ls=":", lw=0.8, alpha=0.7)
ax.set_xlabel(r"$t$ ($\mu$s)")
ax.set_ylabel(r"$\bar{n}(t)$")
ax.set_title(r"Cooling: $\bar{n}(0)=2\bar{n}_\mathrm{ss}$")
ax.legend(fontsize=7)
ax.set_ylim(bottom=0)

plt.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"phase2_transients.{ext}")
plt.close(fig)
print("  Saved: phase2_transients")

# Figure 5: P_n distributions
fig, ax = plt.subplots(figsize=(3.375, 2.8))

n_show = 12
n_arr = np.arange(n_show)

for i, (label, dat) in enumerate(pn_data.items()):
    ax.semilogy(n_arr, np.maximum(dat["pn"][:n_show], 1e-12),
                "o-", ms=3, color=COLORS[i], lw=1.0, label=f"{label} (n̄={dat['n_ss']:.2f})")

ax.set_xlabel(r"Fock level $n$")
ax.set_ylabel(r"$P_n$")
ax.set_title("Photon-number distributions")
ax.legend(fontsize=6, loc="upper right")
ax.set_xlim(-0.5, n_show - 0.5)
ax.set_ylim(1e-7, 1.5)
plt.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"phase2_pn_conditions.{ext}")
plt.close(fig)
print("  Saved: phase2_pn_conditions")

# ---------------------------------------------------------------------------
# Answer scientific questions
# ---------------------------------------------------------------------------
print("\n--- Scientific questions ---")

# Q1: Minimum coupling fraction for measurable signal
print("\n  Q1: Min coupling fraction for Δn̄ > 0.05 at n_target = 2.0")
for kf in np.linspace(0.05, 0.5, 20):
    baths = baths_from_params(2.0, kf)
    n_ss = analytic_nbar_ss(baths)
    n_bg_v = analytic_nbar_ss(baths_from_params(0.0, kf))
    delta = n_ss - n_bg_v
    if delta > MIN_DETECTABLE:
        print(f"    κ_t/κ_tot = {kf:.3f}: Δn̄ = {delta:.4f} > {MIN_DETECTABLE}")
        break

# Q2: Degeneracy analysis (same n̄_ss from different (n_target, kappa_frac) pairs)
print("\n  Q2: Degeneracy — pairs with same n̄_ss ≈ 0.5")
TARGET_NSS = 0.5
deg_pairs = []
for n_t in np.linspace(0.1, 10, 30):
    for kf in np.linspace(0.01, 0.99, 30):
        baths = baths_from_params(n_t, kf)
        nss = analytic_nbar_ss(baths)
        if abs(nss - TARGET_NSS) < 0.01:
            deg_pairs.append((n_t, kf, nss))

if deg_pairs:
    print(f"    Found {len(deg_pairs)} parameter pairs giving n̄_ss ≈ {TARGET_NSS}:")
    for n_t, kf, nss in deg_pairs[:5]:
        print(f"    n_target={n_t:.2f}, κ_frac={kf:.2f} → n̄_ss={nss:.4f}")

print("\nQ3: κ_tot (from transient) partially breaks degeneracy — see Phase 4.")

# ---------------------------------------------------------------------------
# Save data
# ---------------------------------------------------------------------------
np.savez(
    DATA_DIR / "phase2_results.npz",
    n_target_arr=n_target_arr,
    kf_vals=np.array(kf_vals),
    nss_vs_ntarget=nss_vs_ntarget,
    n_floor_arr=n_floor_arr,
    kappa_frac_arr=kappa_frac_arr,
    n_target_fixed_vals=np.array(n_target_fixed_vals),
    nss_vs_kappafrac=nss_vs_kappafrac,
    N_TARGET_GRID=N_TARGET_GRID,
    KAPPA_FRAC_GRID=KAPPA_FRAC_GRID,
    SIGNAL=SIGNAL,
    tlist=tlist,
)
print("\nData saved: data/phase2_results.npz")

print("\n" + "=" * 60)
print("Phase 2 complete.")
print("=" * 60)
