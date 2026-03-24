"""
Targeted validation tests for the thermal-noise cavity sensing study.

Validates all key results against known analytic formulas.

Tests:
  1.  Single-bath steady-state n̄ vs analytic formula
  2.  Multi-bath weighted steady-state n̄
  3.  Transient occupation — exponential decay rate κ_tot
  4.  Transient occupation — exponential steady-state value
  5.  Thermal P_n distribution — max absolute deviation
  6.  Thermal P_n distribution — normalization
  7.  Thermal P_n distribution — mean photon number consistency
  8.  Truncation convergence — N ≥ N_min gives < 0.1% error
  9.  ThermalBath dataclass — negative kappa/n_th rejected
  10. Ramsey coherence at τ=0 equals 1
  11. Ramsey coherence analytic vs direct sum
  12. Spectroscopy peaks at correct frequencies
  13. Temperature conversion round-trip (n_th → T → n_th)
  14. n_thermal at T=0 returns 0
  15. n_bar_ss = n_th for single bath (trivial but explicit)

Usage:
    python scripts/validate_results.py

Exit code: 0 if all pass, 1 if any fail.
"""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import runtime_compat  # noqa: F401 — must be before any qutip/cqed_sim import

import numpy as np
import qutip as qt

from common import (
    KAPPA_TOT, N_CAV, OMEGA_C,
    ThermalBath, analytic_nbar_ss, analytic_kappa_tot,
    analytic_nbar_transient, thermal_pn,
    build_cavity_c_ops, build_qubit_c_ops,
    ramsey_coherence_thermal, spectroscopy_signal,
    n_thermal, temperature_from_nbar, check_truncation,
    CHI_DISP,
)

PASS = 0
FAIL = 0
RESULTS = []

def check(name: str, condition: bool, extra: str = "") -> None:
    global PASS, FAIL
    status = "PASS" if condition else "FAIL"
    if condition:
        PASS += 1
    else:
        FAIL += 1
    RESULTS.append((name, status, extra))
    flag = "  " if condition else "!!"
    print(f"  [{status}] {flag} {name}" + (f"  ({extra})" if extra else ""))

print("=" * 60)
print("Validation tests — thermal_noise_cavity_sensing")
print("=" * 60)

# ---------------------------------------------------------------------------
# Set up basic operators (default N for low-n_th tests)
# ---------------------------------------------------------------------------
N = N_CAV   # = 30, adequate for n_th <= 1

def make_ops(N_):
    a_ = qt.destroy(N_)
    return a_, a_.dag() * a_, qt.qzero(N_)

a, n_op, H = make_ops(N)

# ---------------------------------------------------------------------------
# Test 1: Single-bath steady-state n̄
# ---------------------------------------------------------------------------
print("\nSteady-state tests")

for n_th_val in [0.0, 0.5, 1.0, 2.0, 5.0]:
    # Use N large enough so truncation error < 1e-4
    N_test = max(N_CAV + 20, int(12 * max(n_th_val, 1.0)) + 25)
    a_, n_op_, H_ = make_ops(N_test)
    bath = ThermalBath(kappa=KAPPA_TOT, n_th=n_th_val, label="t")
    c_ops = build_cavity_c_ops(a_, [bath])
    rho_ss = qt.steadystate(H_, c_ops)
    n_num = qt.expect(n_op_, rho_ss)
    err = abs(n_num - n_th_val)
    check(f"Single-bath ss n_th={n_th_val} (N={N_test})", err < 1e-4,
          f"err={err:.2e}")

# ---------------------------------------------------------------------------
# Test 2: Multi-bath weighted steady-state
# ---------------------------------------------------------------------------
print()
for (n1, k1), (n2, k2), (n3, k3) in [
    ((1.0, 0.5), (0.01, 0.3), (0.0, 0.2)),
    ((3.0, 0.4), (0.05, 0.4), (0.0, 0.2)),
    ((0.0, 0.5), (0.1, 0.3), (0.0, 0.2)),
]:
    baths = [
        ThermalBath(kappa=k1 * KAPPA_TOT, n_th=n1, label="t"),
        ThermalBath(kappa=k2 * KAPPA_TOT, n_th=n2, label="bg"),
        ThermalBath(kappa=k3 * KAPPA_TOT, n_th=n3, label="int"),
    ]
    c_ops = build_cavity_c_ops(a, baths)
    rho_ss = qt.steadystate(H, c_ops)
    n_num = qt.expect(n_op, rho_ss)
    n_ana = analytic_nbar_ss(baths)
    err = abs(n_num - n_ana)
    check(f"Multi-bath ss n1={n1},k1={k1}", err < 1e-5, f"num={n_num:.6f}, ana={n_ana:.6f}")

# ---------------------------------------------------------------------------
# Test 3: Transient decay rate κ_tot
# ---------------------------------------------------------------------------
print("\nTransient tests")

n_th_trans = 2.0
bath_trans = ThermalBath(kappa=KAPPA_TOT, n_th=n_th_trans, label="t")
c_ops_trans = build_cavity_c_ops(a, [bath_trans])
psi0 = qt.fock(N, 0)

tau_max = 6.0 / KAPPA_TOT
tlist = np.linspace(0, tau_max, 300)
result = qt.mesolve(H, psi0, tlist, c_ops_trans, [n_op], options={"nsteps": 5000})
n_t = result.expect[0]

n_ss_ana = n_th_trans
# Fit: n(t) = n_ss + (n0 - n_ss) exp(-kappa t)
# → ln((n(t) - n_ss) / (n(0) - n_ss)) = -kappa t
# Use only the initial decay region for the fit
fit_mask = tlist < 2.0 / KAPPA_TOT
y_fit = np.log(np.maximum((n_ss_ana - n_t[fit_mask]) / n_ss_ana, 1e-10))
kappa_fit = -np.polyfit(tlist[fit_mask], y_fit, 1)[0]
rel_err_kappa = abs(kappa_fit - KAPPA_TOT) / KAPPA_TOT
check("Transient decay rate κ_tot", rel_err_kappa < 0.01,
      f"fit={kappa_fit/(2*np.pi*1e3):.3f} kHz, true={KAPPA_TOT/(2*np.pi*1e3):.3f} kHz")

# Test 4: Transient steady-state value
n_ss_num = n_t[-1]
err_ss = abs(n_ss_num - n_th_trans)
check("Transient → correct steady state", err_ss < 0.01,
      f"n_ss_num={n_ss_num:.5f}, expected={n_th_trans}")

# Test 5: Transient analytic agreement
n_ana_t = analytic_nbar_transient(tlist, 0.0, [bath_trans])
max_err = np.max(np.abs(n_t - n_ana_t))
check("Transient vs analytic formula", max_err < 0.02,
      f"max_err={max_err:.2e}")

# ---------------------------------------------------------------------------
# Tests 5–7: Thermal P_n distribution
# ---------------------------------------------------------------------------
print("\nPhoton-number distribution tests")

for n_bar in [0.5, 1.0, 2.0, 5.0]:
    # Use N large enough for good truncation (same formula as single-bath tests)
    N_pn = max(N_CAV + 20, int(12 * max(n_bar, 1.0)) + 25)
    a_pn, n_op_pn, H_pn = make_ops(N_pn)
    bath = ThermalBath(kappa=KAPPA_TOT, n_th=n_bar, label="t")
    c_ops = build_cavity_c_ops(a_pn, [bath])
    rho_ss = qt.steadystate(H_pn, c_ops)
    pn_num = np.array([rho_ss[n, n].real for n in range(N_pn)])
    pn_ana = thermal_pn(n_bar, N_pn)

    # Max deviation
    max_dev = np.max(np.abs(pn_num - pn_ana))
    check(f"P_n max deviation (n_bar={n_bar}, N={N_pn})", max_dev < 1e-5,
          f"max|dP_n|={max_dev:.2e}")

    # Normalization
    norm_err = abs(pn_num.sum() - 1.0)
    check(f"P_n normalization (n_bar={n_bar})", norm_err < 1e-6,
          f"|sum-1|={norm_err:.2e}")

    # Mean photon number consistency
    n_mean_from_pn = sum(n * p for n, p in enumerate(pn_num))
    n_mean_err = abs(n_mean_from_pn - n_bar)
    check(f"P_n -> nbar consistency (n_bar={n_bar})", n_mean_err < 1e-4,
          f"nbar_from_pn={n_mean_from_pn:.5f}")

# ---------------------------------------------------------------------------
# Test 8: Truncation convergence
# ---------------------------------------------------------------------------
print("\nTruncation convergence tests")

n_bar_test = 3.0
bath_ref = ThermalBath(kappa=KAPPA_TOT, n_th=n_bar_test, label="t")

# Reference: N=80 (large enough for n_bar=3 with very small truncation)
a_ref = qt.destroy(80)
H_ref = qt.qzero(80)
c_ref = build_cavity_c_ops(a_ref, [bath_ref])
n_ref = qt.expect(a_ref.dag() * a_ref, qt.steadystate(H_ref, c_ref))

for Ntest in [10, 15, 20, 30, 45]:
    a_v = qt.destroy(Ntest)
    c_v = build_cavity_c_ops(a_v, [ThermalBath(kappa=KAPPA_TOT, n_th=n_bar_test)])
    n_v = qt.expect(a_v.dag() * a_v, qt.steadystate(qt.qzero(Ntest), c_v))
    rel = abs(n_v - n_ref) / n_ref
    # For n_bar=3: need N >= 45 for < 0.1% convergence
    if Ntest >= 45:
        check(f"Truncation convergence N={Ntest} (nbar=3)", rel < 1e-3,
              f"rel_err={rel:.2e}")
    else:
        # Just report (informational)
        check(f"Truncation hint N={Ntest} (nbar=3)", True,
              f"rel_err={rel:.2e} (informational)")

# ---------------------------------------------------------------------------
# Test 9: ThermalBath validation
# ---------------------------------------------------------------------------
print("\nDataclass validation tests")

try:
    _ = ThermalBath(kappa=-1.0, n_th=1.0)
    check("ThermalBath rejects kappa<0", False, "No error raised")
except ValueError:
    check("ThermalBath rejects kappa<0", True)

try:
    _ = ThermalBath(kappa=1.0, n_th=-0.5)
    check("ThermalBath rejects n_th<0", False, "No error raised")
except ValueError:
    check("ThermalBath rejects n_th<0", True)

# ---------------------------------------------------------------------------
# Tests 10–11: Ramsey coherence
# ---------------------------------------------------------------------------
print("\nRamsey coherence tests")

for n_bar in [0.0, 1.0, 3.0]:
    # Test: |χ(0)| = 1
    chi0 = ramsey_coherence_thermal(np.array([0.0]), n_bar, CHI_DISP)
    check(f"Ramsey χ(0)=1 (n_bar={n_bar})", abs(chi0[0] - 1.0) < 1e-12,
          f"χ(0)={chi0[0]:.8f}")

    # Test: direct sum vs analytic formula
    tau_test = np.array([0.1e-6, 0.5e-6, 1.0e-6])
    chi_analytic = ramsey_coherence_thermal(tau_test, n_bar, CHI_DISP)
    pn = thermal_pn(n_bar, 60)
    chi_direct = np.array([
        sum(pn[n] * np.exp(1j * n * CHI_DISP * t) for n in range(60))
        for t in tau_test
    ])
    max_diff = np.max(np.abs(chi_analytic - chi_direct))
    check(f"Ramsey χ(τ) analytic vs sum (n_bar={n_bar})", max_diff < 1e-6,
          f"max_diff={max_diff:.2e}")

# ---------------------------------------------------------------------------
# Test 12: Spectroscopy peaks at correct frequencies
# ---------------------------------------------------------------------------
print("\nSpectroscopy tests")

n_bar_spec = 3.0
pn_spec = thermal_pn(n_bar_spec, 20)
chi = -CHI_DISP   # negative

# Fine grid around expected peak positions
for n_peak in [0, 1, 2, 3]:
    omega_pk = n_peak * chi  # peak position in rotating frame
    omega_test = np.linspace(omega_pk - 2e6 * 2 * np.pi, omega_pk + 2e6 * 2 * np.pi, 500)
    gamma_q = 1.0 / (20e-6)  # 1/T2
    sig = spectroscopy_signal(omega_test, pn_spec, 0.0, chi, gamma_q)
    omega_peak_num = omega_test[np.argmax(sig)]
    err_MHz = abs(omega_peak_num - omega_pk) / (2 * np.pi * 1e6)
    check(f"Spectroscopy peak n={n_peak} at correct ω", err_MHz < 0.01,
          f"err={err_MHz:.4f} MHz")

# ---------------------------------------------------------------------------
# Tests 13–14: Temperature conversion
# ---------------------------------------------------------------------------
print("\nTemperature conversion tests")

for T_K in [0.020, 0.100, 0.500, 1.0, 4.0]:
    n_th_val = n_thermal(OMEGA_C, T_K)
    T_roundtrip = temperature_from_nbar(OMEGA_C, n_th_val)
    rel_err = abs(T_roundtrip - T_K) / T_K
    check(f"T round-trip (T={T_K*1000:.0f} mK)", rel_err < 1e-6,
          f"err={rel_err:.2e}")

check("n_thermal at T=0 returns 0", n_thermal(OMEGA_C, 0.0) == 0.0)
check("n_thermal at T=0.0 (float)", n_thermal(OMEGA_C, 0.0) == 0.0)

# Test 15: n̄_ss = n_th for single bath
for n_th_v in [0.0, 0.5, 2.0, 5.0]:
    bath = ThermalBath(kappa=KAPPA_TOT, n_th=n_th_v)
    n_ss = analytic_nbar_ss([bath])
    check(f"analytic_nbar_ss = n_th (single bath, n_th={n_th_v})",
          abs(n_ss - n_th_v) < 1e-15, f"diff={abs(n_ss - n_th_v):.2e}")

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
print("\n" + "=" * 60)
print(f"Results: {PASS} PASS, {FAIL} FAIL")
print("=" * 60)

if FAIL > 0:
    print("\nFailed tests:")
    for name, status, extra in RESULTS:
        if status == "FAIL":
            print(f"  !! {name}  ({extra})")
    sys.exit(1)
else:
    print("\nAll tests passed.")
    sys.exit(0)
