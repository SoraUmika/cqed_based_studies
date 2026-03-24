# -*- coding: utf-8 -*-
"""
Phase 1: Minimal finite-temperature cavity model - implementation and validation.

Validates the multi-bath Lindblad master equation against three known analytic results:
  1. Steady-state mean photon number: n_ss = (sum kj nj) / (sum kj)
  2. Transient occupation: n(t) = n_ss + (n(0) - n_ss) exp(-kappa_tot t)
  3. Thermal photon-number distribution: P_n = nbar^n / (1+nbar)^{n+1}

Also includes a truncation-convergence study.

Usage:
    python scripts/phase1_cavity_model.py

Output:
    data/phase1_results.npz
    figures/phase1_steady_state_validation.{png,pdf}
    figures/phase1_transient_validation.{png,pdf}
    figures/phase1_pn_distribution.{png,pdf}
    figures/phase1_truncation_convergence.{png,pdf}
"""

import sys
import time
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
# Force UTF-8 output on Windows to avoid cp1252 encoding errors
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
)

# ---------------------------------------------------------------------------
# Load plot style
# ---------------------------------------------------------------------------
if STYLE_PATH.exists():
    plt.style.use(STYLE_PATH)

# Colour palette from style
COLORS = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377"]

# ---------------------------------------------------------------------------
# Parameters for this phase
# ---------------------------------------------------------------------------

N = 80             # Fock truncation for validation (larger than default for accuracy)
N_STEPS = 400      # time-evolution steps for transient study
TOL_SS = 1e-4      # tolerance for steady-state validation (some truncation error expected for n_th=5)
TOL_TRANSIENT = 5e-4
TOL_PN = 1e-4      # max absolute deviation in P_n

print("=" * 60)
print("Phase 1: Cavity Thermal Model — Validation")
print("=" * 60)

# ---------------------------------------------------------------------------
# Test 1: Single-bath steady state
# ---------------------------------------------------------------------------
print("\n--- Test 1: Single-bath steady state ---")

kappa_single = KAPPA_TOT
n_targets = [0.0, 0.5, 1.0, 2.0, 5.0]

a = qt.destroy(N)
n_op = a.dag() * a
H = qt.qzero(N)   # rotating frame, zero drive

ss_numeric = []
ss_analytic = []

for n_th in n_targets:
    bath = ThermalBath(kappa=kappa_single, n_th=n_th, label="single")
    c_ops = build_cavity_c_ops(a, [bath])
    rho_ss = qt.steadystate(H, c_ops)
    n_num = qt.expect(n_op, rho_ss)
    n_ana = n_th      # for single bath, n_ss = n_th
    ss_numeric.append(n_num)
    ss_analytic.append(n_ana)
    err = abs(n_num - n_ana)
    status = "PASS" if err < TOL_SS else "FAIL"
    print(f"  n_th={n_th:.1f}: numeric={n_num:.8f}, analytic={n_ana:.8f}, "
          f"err={err:.2e}  [{status}]")

ss_numeric = np.array(ss_numeric)
ss_analytic = np.array(ss_analytic)

# ---------------------------------------------------------------------------
# Test 2: Multi-bath steady state
# ---------------------------------------------------------------------------
print("\n--- Test 2: Multi-bath steady state ---")

# Three-bath scenario (target + background + internal)
multi_bath_cases = [
    dict(n_target=1.0, kf_t=0.5, n_bg=0.01, kf_bg=0.3, n_int=0.0, kf_int=0.2),
    dict(n_target=3.0, kf_t=0.4, n_bg=0.05, kf_bg=0.4, n_int=0.0, kf_int=0.2),
    dict(n_target=5.0, kf_t=0.6, n_bg=0.02, kf_bg=0.2, n_int=0.0, kf_int=0.2),
]

multi_ss_numeric = []
multi_ss_analytic = []

for case in multi_bath_cases:
    baths = [
        ThermalBath(kappa=case["kf_t"]   * KAPPA_TOT, n_th=case["n_target"], label="target"),
        ThermalBath(kappa=case["kf_bg"]  * KAPPA_TOT, n_th=case["n_bg"],     label="bg"),
        ThermalBath(kappa=case["kf_int"] * KAPPA_TOT, n_th=case["n_int"],    label="int"),
    ]
    c_ops = build_cavity_c_ops(a, baths)
    rho_ss = qt.steadystate(H, c_ops)
    n_num = qt.expect(n_op, rho_ss)
    n_ana = analytic_nbar_ss(baths)
    err = abs(n_num - n_ana)
    status = "PASS" if err < TOL_SS else "FAIL"
    multi_ss_numeric.append(n_num)
    multi_ss_analytic.append(n_ana)
    print(f"  n_t={case['n_target']}, kf_t={case['kf_t']}: "
          f"numeric={n_num:.8f}, analytic={n_ana:.8f}, err={err:.2e}  [{status}]")

# ---------------------------------------------------------------------------
# Test 3: Thermal photon-number distribution P_n
# ---------------------------------------------------------------------------
print("\n--- Test 3: Thermal P_n distribution ---")

pn_n_bars = [0.5, 1.0, 2.0, 5.0]
pn_results = {}

for n_bar in pn_n_bars:
    check_truncation(n_bar, N)
    bath = ThermalBath(kappa=KAPPA_TOT, n_th=n_bar, label="single")
    c_ops = build_cavity_c_ops(a, [bath])
    rho_ss = qt.steadystate(H, c_ops)

    # Extract numeric P_n from diagonal of density matrix
    pn_num = np.array([rho_ss[n, n].real for n in range(N)])
    pn_ana = thermal_pn(n_bar, N)

    max_dev = np.max(np.abs(pn_num - pn_ana))
    kl_div = np.sum(pn_ana[pn_ana > 1e-15] *
                    np.log(pn_ana[pn_ana > 1e-15] / np.maximum(pn_num[pn_ana > 1e-15], 1e-15)))
    status = "PASS" if max_dev < TOL_PN else "FAIL"
    pn_results[n_bar] = {"numeric": pn_num, "analytic": pn_ana}
    print(f"  n_bar={n_bar:.1f}: max|dP_n|={max_dev:.2e}, KL={kl_div:.2e}  [{status}]")

# ---------------------------------------------------------------------------
# Test 4: Transient occupation
# ---------------------------------------------------------------------------
print("\n--- Test 4: Transient occupation (single bath) ---")

n_th_trans = 2.0
kappa_trans = KAPPA_TOT
bath_trans = ThermalBath(kappa=kappa_trans, n_th=n_th_trans, label="single")
c_ops_trans = build_cavity_c_ops(a, [bath_trans])

# Start from vacuum |0>
psi0 = qt.fock(N, 0)
tau_max = 5.0 / kappa_trans
tlist = np.linspace(0, tau_max, N_STEPS)

t0 = time.time()
result = qt.mesolve(H, psi0, tlist, c_ops_trans, [n_op],
                   options={"nsteps": 20000, "atol": 1e-10, "rtol": 1e-8})
t_elapsed = time.time() - t0
n_bar_t = result.expect[0]
n_bar_ana_t = analytic_nbar_transient(tlist, 0.0, [bath_trans])

max_err_transient = np.max(np.abs(n_bar_t - n_bar_ana_t))
status = "PASS" if max_err_transient < TOL_TRANSIENT else "FAIL"
print(f"  n_th={n_th_trans}, start from vacuum: max err={max_err_transient:.2e}, "
      f"elapsed={t_elapsed:.2f}s  [{status}]")

# Also test cooling (start from hot state)
n_init_hot = 5.0
rho_hot = qt.thermal_dm(N, n_init_hot)  # thermal state with n_bar=5
bath_cold = ThermalBath(kappa=KAPPA_TOT, n_th=0.5, label="cold")
c_ops_cool = build_cavity_c_ops(a, [bath_cold])
result_cool = qt.mesolve(H, rho_hot, tlist, c_ops_cool, [n_op],
                        options={"nsteps": 20000, "atol": 1e-10, "rtol": 1e-8})
n_bar_cool = result_cool.expect[0]
n_bar_cool_ana = analytic_nbar_transient(tlist, n_init_hot, [bath_cold])
max_err_cool = np.max(np.abs(n_bar_cool - n_bar_cool_ana))
status_cool = "PASS" if max_err_cool < TOL_TRANSIENT else "FAIL"
print(f"  Cooling (hot->cold): max err={max_err_cool:.2e}  [{status_cool}]")

# ---------------------------------------------------------------------------
# Test 5: Truncation convergence
# ---------------------------------------------------------------------------
print("\n--- Test 5: Truncation convergence ---")

N_values = [5, 8, 10, 15, 20, 25, 30, 40, 50]
n_bar_trunc = 2.0   # challenging case for small N
bath_trunc = ThermalBath(kappa=KAPPA_TOT, n_th=n_bar_trunc, label="trunc")

trunc_nbar = []
for Nv in N_values:
    a_v = qt.destroy(Nv)
    n_v = a_v.dag() * a_v
    H_v = qt.qzero(Nv)
    c_v = build_cavity_c_ops(a_v, [bath_trunc])
    rho_v = qt.steadystate(H_v, c_v)
    trunc_nbar.append(qt.expect(n_v, rho_v))

trunc_nbar = np.array(trunc_nbar)
ref_nbar = trunc_nbar[-1]   # N=50 as reference
trunc_err = np.abs(trunc_nbar - ref_nbar)
print("  N      n_ss      err vs N=50")
for Nv, nb, err in zip(N_values, trunc_nbar, trunc_err):
    flag = " *" if err > 1e-3 else ""
    print(f"  {Nv:3d}    {nb:.6f}  {err:.2e}{flag}")

# ---------------------------------------------------------------------------
# Generate figures
# ---------------------------------------------------------------------------
print("\nGenerating figures...")

# Figure 1: Steady-state validation (single + multi bath)
fig, axes = plt.subplots(1, 2, figsize=(6.75, 2.5))

ax = axes[0]
ax.plot(ss_analytic, ss_analytic, "k--", lw=1.0, label="Analytic", zorder=0)
ax.plot(ss_analytic, ss_numeric, "o", ms=5, color=COLORS[0], label="Numeric (ss)")
ax.set_xlabel(r"$\bar{n}_\mathrm{analytic}$")
ax.set_ylabel(r"$\bar{n}_\mathrm{numeric}$")
ax.set_title("Single-bath steady state")
ax.legend()

ax = axes[1]
x = np.array(multi_ss_analytic)
ax.plot(x, x, "k--", lw=1.0, zorder=0)
ax.plot(x, np.array(multi_ss_numeric), "s", ms=5, color=COLORS[1], label="Numeric (ss)")
ax.set_xlabel(r"$\bar{n}_\mathrm{analytic}$")
ax.set_ylabel(r"$\bar{n}_\mathrm{numeric}$")
ax.set_title("Multi-bath steady state")
ax.legend()

plt.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"phase1_steady_state_validation.{ext}")
plt.close(fig)
print("  Saved: phase1_steady_state_validation")

# Figure 2: Transient validation
fig, axes = plt.subplots(1, 2, figsize=(6.75, 2.5))

kappa_khz = kappa_trans / (2 * np.pi * 1e3)
t_us = tlist * 1e6

ax = axes[0]
ax.plot(t_us, n_bar_ana_t, "--", color="k", lw=1.2, label="Analytic")
ax.plot(t_us, n_bar_t, color=COLORS[0], lw=1.2, label="Numeric")
ax.set_xlabel(r"$t$ ($\mu$s)")
ax.set_ylabel(r"$\bar{n}(t)$")
ax.set_title(rf"Heating: $n_0=0$, $n_{{th}}={n_th_trans}$")
ax.legend()

ax = axes[1]
ax.plot(t_us, n_bar_cool_ana, "--", color="k", lw=1.2, label="Analytic")
ax.plot(t_us, n_bar_cool, color=COLORS[1], lw=1.2, label="Numeric")
ax.set_xlabel(r"$t$ ($\mu$s)")
ax.set_ylabel(r"$\bar{n}(t)$")
ax.set_title(rf"Cooling: $n_0={n_init_hot}$, $n_{{th}}=0.5$")
ax.legend()

plt.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"phase1_transient_validation.{ext}")
plt.close(fig)
print("  Saved: phase1_transient_validation")

# Figure 3: Thermal P_n distribution
fig, axes = plt.subplots(1, 2, figsize=(6.75, 2.5))

n_show = 15
n_arr = np.arange(n_show)

ax = axes[0]
for i, (n_bar, dat) in enumerate(pn_results.items()):
    pn_a = dat["analytic"][:n_show]
    pn_n = dat["numeric"][:n_show]
    ax.semilogy(n_arr, pn_a, "--", color=COLORS[i], lw=1.0)
    ax.semilogy(n_arr, np.maximum(pn_n, 1e-15), "o", ms=3, color=COLORS[i],
                label=rf"$\bar{{n}}={n_bar}$")
ax.set_xlabel(r"Fock level $n$")
ax.set_ylabel(r"$P_n$")
ax.set_title("Thermal photon-number distributions")
ax.legend(fontsize=7)
ax.set_ylim(1e-8, 1.5)

ax = axes[1]
for i, (n_bar, dat) in enumerate(pn_results.items()):
    pn_a = dat["analytic"][:n_show]
    pn_n = dat["numeric"][:n_show]
    dev = np.abs(pn_n - pn_a)
    ax.semilogy(n_arr, np.maximum(dev, 1e-15), color=COLORS[i],
                lw=1.0, label=rf"$\bar{{n}}={n_bar}$")
ax.axhline(TOL_PN, color="gray", lw=0.8, ls="--", label=f"tol={TOL_PN}")
ax.set_xlabel(r"Fock level $n$")
ax.set_ylabel(r"$|P_n^\mathrm{num} - P_n^\mathrm{ana}|$")
ax.set_title("P_n deviation from analytic")
ax.legend(fontsize=7)

plt.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"phase1_pn_distribution.{ext}")
plt.close(fig)
print("  Saved: phase1_pn_distribution")

# Figure 4: Truncation convergence
fig, ax = plt.subplots(figsize=(3.375, 2.5))

ax.semilogy(N_values, np.maximum(trunc_err, 1e-12), "o-", color=COLORS[0], lw=1.2, ms=4)
ax.axhline(1e-3, color="gray", ls="--", lw=0.8, label="0.1% threshold")
ax.set_xlabel(r"Fock-space truncation $N$")
ax.set_ylabel(r"$|\bar{n}(N) - \bar{n}(N{=}40)|$")
ax.set_title(rf"Truncation convergence ($\bar{{n}}_{{th}}={n_bar_trunc}$)")
ax.legend()

plt.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"phase1_truncation_convergence.{ext}")
plt.close(fig)
print("  Saved: phase1_truncation_convergence")

# ---------------------------------------------------------------------------
# Save data
# ---------------------------------------------------------------------------
np.savez(
    DATA_DIR / "phase1_results.npz",
    # Test 1: single-bath ss
    ss_n_targets=np.array(n_targets),
    ss_numeric=ss_numeric,
    ss_analytic=ss_analytic,
    # Test 2: multi-bath ss
    multi_ss_numeric=np.array(multi_ss_numeric),
    multi_ss_analytic=np.array(multi_ss_analytic),
    # Test 3: P_n distributions
    pn_n_bars=np.array(pn_n_bars),
    # Test 4: transients
    tlist=tlist,
    n_bar_transient=n_bar_t,
    n_bar_transient_analytic=n_bar_ana_t,
    n_bar_cool=n_bar_cool,
    n_bar_cool_analytic=n_bar_cool_ana,
    # Test 5: truncation convergence
    trunc_N=np.array(N_values),
    trunc_nbar=trunc_nbar,
)
print("\nData saved: data/phase1_results.npz")

print("\n" + "=" * 60)
print("Phase 1 complete.")
print("=" * 60)
