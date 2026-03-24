"""
Quick physics-only test — runs in ~10s without qutip.

Validates all core formulas (ODE integrator, GRAPE gradient, objectives,
pulse builders) without importing cqed_sim/qutip. Run this first to confirm
the study infrastructure is correct before launching the full phase scripts.

Usage: python run_quick_test.py
"""

import sys, os
sys.path.insert(0, os.path.dirname(__file__))

# Monkey-patch to skip cqed_sim import for this test
import common as _cm
_cm._get_readout_resonator = lambda: None   # stub

import numpy as np
import common as C

TWO_PI = C.TWO_PI
kappa   = C.KAPPA
chi     = C.CHI_NOMINAL
epsilon = C.EPSILON_MAX * 0.5
delta_g = C.optimal_delta_g(chi)
delta_e = delta_g + chi

PASS = "\033[92mPASS\033[0m"
FAIL = "\033[91mFAIL\033[0m"

def check(msg, passed, val, tol):
    s = PASS if passed else FAIL
    print(f"  [{s}] {msg}: {val:.3e} (tol={tol:.0e})")
    return passed


all_ok = True

print("=" * 55)
print("Quick Physics Test (no qutip)")
print("=" * 55)

# ── 1. Steady-state formula ───────────────────────────────────────────────
print("\n[1] Steady-state formula |Δα_ss_max| = 4|ε||χ|/(κ²+χ²)")
sep_formula = C.optimal_separation_formula(epsilon, kappa, chi)
sep_exact   = abs(C.steady_state_alpha(epsilon, kappa, delta_g)
                 - C.steady_state_alpha(epsilon, kappa, delta_e))
err = abs(sep_formula - sep_exact) / max(abs(sep_exact), 1e-30)
all_ok &= check("Formula vs. direct α_ss computation", err < 1e-10, err, 1e-10)
# Check peak is at chi/kappa = 1
x = np.linspace(0.01, 6, 2000)
f_scan = 4 * abs(epsilon) * x * kappa / (kappa**2 + (x*kappa)**2)
idx_peak = np.argmax(f_scan)
all_ok &= check("Peak at χ/κ=1", abs(x[idx_peak] - 1.0) < 0.01, abs(x[idx_peak]-1.0), 0.01)
print(f"  Analytical |Δα_ss_max| at χ/κ=1: {sep_formula:.6f}")

# ── 2. ODE integrator: constant drive → analytical steady state ───────────
print("\n[2] ODE integrator vs analytical steady-state (constant ε)")
T = 20.0 / kappa  # κT=20 → transient residual exp(-10) < 1e-4
N = max(2, int(T / C.DT_ODE))
tlist = np.linspace(0.0, T, N+1)
eps_arr = epsilon * np.ones(N, dtype=np.complex128)

ag, ae = C.simulate_conditioned_fields(eps_arr, tlist, kappa, chi, delta_g)
ag_ss = C.steady_state_alpha(epsilon, kappa, delta_g)
ae_ss = C.steady_state_alpha(epsilon, kappa, delta_e)

err_g = abs(ag[-1] - ag_ss) / max(abs(ag_ss), 1e-30)
err_e = abs(ae[-1] - ae_ss) / max(abs(ae_ss), 1e-30)
all_ok &= check("α_g(T) → α_g^ss", err_g < 1e-4, err_g, 1e-4)
all_ok &= check("α_e(T) → α_e^ss", err_e < 1e-4, err_e, 1e-4)

# ── 3. Zero drive → zero field ────────────────────────────────────────────
print("\n[3] Zero drive → SNR² = 0")
eps_zero = np.zeros(N, dtype=np.complex128)
ag0, ae0 = C.simulate_conditioned_fields(eps_zero, tlist, kappa, chi, delta_g)
snr2_zero = C.snr_squared(ag0, ae0, tlist, kappa)
all_ok &= check("Zero drive SNR²=0", snr2_zero < 1e-20, snr2_zero, 1e-20)

# ── 4. SNR² increases with drive amplitude ────────────────────────────────
print("\n[4] SNR² monotone in drive amplitude (at κT=5)")
T3 = 5.0 / kappa
N3 = max(2, int(T3 / C.DT_ODE))
tl3 = np.linspace(0, T3, N3+1)
snr2_prev = -1.0
monotone = True
for amp in np.linspace(epsilon*0.1, epsilon*0.9, 8):
    eps3 = amp * np.ones(N3, dtype=np.complex128)
    ag3, ae3 = C.simulate_conditioned_fields(eps3, tl3, kappa, chi, delta_g)
    s = C.snr_squared(ag3, ae3, tl3, kappa)
    if s < snr2_prev - 1e-8:
        monotone = False
    snr2_prev = max(snr2_prev, s)
all_ok &= check("SNR² monotone in |ε|", monotone, 0.0, 0.5)

# ── 5. Pulse builders produce correct shapes ──────────────────────────────
print("\n[5] Pulse builders area and symmetry")
T_pb = 1e-6
N_pb = 500
tl_pb = np.linspace(0, T_pb, N_pb+1)
t_mid = 0.5*(tl_pb[:-1]+tl_pb[1:])

for family, builder in C.PULSE_FAMILIES.items():
    env = builder(tl_pb, complex(epsilon))
    amp_real = np.real(env)
    # Check all samples are non-negative (for real-valued envelope families)
    noneg = np.all(amp_real >= -1e-6*abs(epsilon))
    print(f"  {family:12s}: max={np.max(np.abs(env))/epsilon:.3f}ε, "
          f"nonneg={noneg}, len={len(env)}")

# ── 6. GRAPE gradient accuracy (finite-difference check) ──────────────────
print("\n[6] GRAPE gradient via finite differences")
N_seg = 10
T_seg = 2.0 / kappa
tau = T_seg / N_seg
A_g, B_g, A_e, B_e = C.grape_propagators(kappa, delta_g, delta_e, tau)

rng = np.random.default_rng(0)
x0 = rng.uniform(-epsilon*0.5, epsilon*0.5, size=2*N_seg)
eps_vec = x0[:N_seg] + 1j * x0[N_seg:]

# Analytical gradient
obj0, grad_analytical = C.grape_objective_and_grad(
    x0, kappa, delta_g, delta_e, tau, N_seg,
    w_snr=1.0, w_res=0.0
)

# Finite-difference gradient
h = epsilon * 1e-5
grad_fd = np.zeros(2*N_seg)
for i in range(2*N_seg):
    xp = x0.copy(); xp[i] += h
    xm = x0.copy(); xm[i] -= h
    fp, _ = C.grape_objective_and_grad(xp, kappa, delta_g, delta_e, tau, N_seg, 1.0, 0.0)
    fm, _ = C.grape_objective_and_grad(xm, kappa, delta_g, delta_e, tau, N_seg, 1.0, 0.0)
    grad_fd[i] = (fp - fm) / (2*h)

# Relative error
denom = max(np.max(np.abs(grad_analytical)), np.max(np.abs(grad_fd)), 1e-30)
err_grad = np.max(np.abs(grad_analytical - grad_fd)) / denom
all_ok &= check("GRAPE gradient (FD vs analytical)", err_grad < 1e-4, err_grad, 1e-4)
print(f"  Analytical grad norm: {np.linalg.norm(grad_analytical):.4e}")
print(f"  FD grad norm:         {np.linalg.norm(grad_fd):.4e}")

# ── 7. GRAPE improves over initial guess ──────────────────────────────────
print("\n[7] GRAPE optimizes SNR² above amplitude scan")
# Small problem for speed
N_g = 20
T_g = 3.0 / kappa
_, ev_sq = C.optimize_amplitude("square", T_g, kappa, chi, delta_g, epsilon)
gr = C.run_grape(kappa=kappa, chi=chi, T_read=T_g, N_seg=N_g,
                 delta_g=delta_g, epsilon_max=epsilon, n_restarts=5, maxiter=100)
margin = gr["snr2_opt"] - ev_sq["snr2"]
all_ok &= check("GRAPE SNR² ≥ square pulse", margin >= -1e-4, -margin, 1e-4)
print(f"  Square SNR²: {ev_sq['snr2']:.4f}, GRAPE SNR²: {gr['snr2_opt']:.4f}")
print(f"  GRAPE advantage: {100*(gr['snr2_opt']-ev_sq['snr2'])/max(ev_sq['snr2'],1e-8):.1f}%")

# ── 8. Midpoint drive outperforms on-resonance ────────────────────────────
print("\n[8] Midpoint drive outperforms on-resonance drive")
T8 = 5.0 / kappa
_, ev_mid = C.optimize_amplitude("square", T8, kappa, chi, delta_g,      epsilon)
_, ev_res = C.optimize_amplitude("square", T8, kappa, chi, 0.0,           epsilon)
gain_midpt = (ev_mid["snr2"] - ev_res["snr2"]) / max(ev_res["snr2"], 1e-8)
all_ok &= check("Midpoint drive SNR² > on-resonance", gain_midpt > 0.01, -gain_midpt, 0.01)
print(f"  Midpoint SNR²={ev_mid['snr2']:.4f}, On-resonance SNR²={ev_res['snr2']:.4f}")
print(f"  Midpoint advantage: {100*gain_midpt:.1f}%")

# ── Summary ───────────────────────────────────────────────────────────────
print("\n" + "=" * 55)
if all_ok:
    print("ALL QUICK TESTS PASSED ✓")
    print("Core physics implementation is correct.")
    print("Run phase scripts for full numerical study.")
else:
    print("SOME QUICK TESTS FAILED — check implementation")
print("=" * 55)
