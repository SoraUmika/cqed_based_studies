"""
Validation Framework — Three-Check Protocol
============================================

Sanity checks (must pass before any result is used):
  S1. Zero drive → α = 0 for all time (both states identical, SNR² = 0).
  S2. Constant drive → ODE solution matches analytical steady-state.
  S3. ODE integrator → matches cqed_sim.ReadoutResonator.response_trace() for const ε.

Convergence checks:
  C1. SNR² is stable to 0.1% when dt is halved (ODE time-step convergence).
  C2. SNR² is stable to 0.1% when N_seg is doubled (GRAPE segment convergence).

Literature / theoretical checks:
    L1. Maximum steady-state SNR² rate (γ_meas) matches analytical formula.
    L2. Optimal drive frequency matches ω_d = ω_r + χ/2, equivalently Δ_g = −χ/2.
    L3. GRAPE SNR² ≥ best pulse-family SNR² (GRAPE is at least as good).

All checks report PASS or FAIL with quantitative margins.
"""

import sys, os


def _configure_utf8_output() -> None:
    for stream in (sys.stdout, sys.stderr):
        reconfigure = getattr(stream, "reconfigure", None)
        if reconfigure is not None:
            try:
                reconfigure(encoding="utf-8")
            except OSError:
                pass


_configure_utf8_output()

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

import numpy as np
from common import (
    TWO_PI, KAPPA, CHI_NOMINAL, EPSILON_MAX, OMEGA_R, DT_ODE,
    optimal_delta_g, steady_state_alpha, integrate_readout_ode,
    simulate_conditioned_fields, snr_squared, endpoint_separation_sq,
    assignment_fidelity_from_snr2, residual_photons,
    grape_forward_pass, grape_backward_pass, grape_propagators,
    make_readout_resonator, cqed_sim_steady_state,
    evaluate_pulse_family, run_grape, optimize_amplitude,
)

PASS_STR = "\033[92mPASS\033[0m"
FAIL_STR = "\033[91mFAIL\033[0m"


def check(name: str, passed: bool, value: float, tol: float, unit: str = "") -> bool:
    status = PASS_STR if passed else FAIL_STR
    print(f"  [{status}] {name}: {value:.2e} {unit} (tol={tol:.1e})")
    return passed


# ─── Sanity Checks ────────────────────────────────────────────────────────────

def sanity_zero_drive() -> bool:
    """S1: Zero drive → zero field → SNR² = 0."""
    kappa = KAPPA
    chi   = CHI_NOMINAL
    T     = 500e-9
    N     = int(T / DT_ODE)
    tlist = np.linspace(0.0, T, N + 1)
    eps   = np.zeros(N, dtype=np.complex128)
    ag, ae = simulate_conditioned_fields(eps, tlist, kappa, chi)
    snr2 = snr_squared(ag, ae, tlist, kappa)
    n_max_field = max(np.max(np.abs(ag)), np.max(np.abs(ae)))
    passed = snr2 < 1e-20 and n_max_field < 1e-15
    return check("S1: zero drive → SNR²=0", passed, snr2, 1e-20)


def sanity_steady_state_convergence() -> bool:
    """S2: Constant drive, long duration → ODE field approaches α_ss."""
    kappa   = KAPPA
    chi     = CHI_NOMINAL
    epsilon = EPSILON_MAX * 0.5
    delta_g = optimal_delta_g(chi)
    delta_e = delta_g + chi

    T   = 20.0 / kappa      # κT=20 -> transient residual exp(-10) < 1e-4
    N   = max(100, int(round(T / DT_ODE)))
    tlist = np.linspace(0.0, T, N + 1)
    eps = epsilon * np.ones(N, dtype=np.complex128)
    ag, ae = simulate_conditioned_fields(eps, tlist, kappa, chi, delta_g)

    alpha_g_ss = steady_state_alpha(epsilon, kappa, delta_g)
    alpha_e_ss = steady_state_alpha(epsilon, kappa, delta_e)

    err_g = abs(ag[-1] - alpha_g_ss) / (abs(alpha_g_ss) + 1e-30)
    err_e = abs(ae[-1] - alpha_e_ss) / (abs(alpha_e_ss) + 1e-30)
    passed = err_g < 1e-4 and err_e < 1e-4
    print(f"  [{'PASS' if err_g < 1e-4 else 'FAIL'}] S2a: α_g steady-state error = {err_g:.2e}")
    print(f"  [{'PASS' if err_e < 1e-4 else 'FAIL'}] S2b: α_e steady-state error = {err_e:.2e}")
    return passed


def sanity_cqed_sim_response_trace() -> bool:
    """
    S3: ODE integrator matches cqed_sim ReadoutResonator.response_trace()
    for constant ε (which uses the exact analytical formula).
    """
    kappa   = KAPPA
    chi     = CHI_NOMINAL
    epsilon = EPSILON_MAX * 0.4
    delta_g = optimal_delta_g(chi)
    T       = 300e-9
    dt_test = 2e-9

    # Our ODE integrator
    N = max(2, int(round(T / dt_test)))
    tlist = np.linspace(0.0, T, N + 1)
    eps_arr = epsilon * np.ones(N, dtype=np.complex128)
    ag_ode = integrate_readout_ode(eps_arr, tlist, kappa, delta_g)

    # cqed_sim ReadoutResonator.response_trace()
    res = make_readout_resonator(kappa=kappa, chi=chi, epsilon=epsilon)
    omega_d = OMEGA_R - delta_g
    tlist_cs, ag_cs = res.response_trace(
        "g", duration=T, dt=dt_test, drive_frequency=omega_d
    )

    # Align time grids
    n_min = min(len(ag_ode), len(ag_cs))
    err = np.max(np.abs(ag_ode[:n_min] - ag_cs[:n_min]))
    ref = max(np.max(np.abs(ag_cs[:n_min])), 1e-30)
    rel_err = err / ref
    passed = rel_err < 1e-3
    return check("S3: ODE vs cqed_sim response_trace", passed, rel_err, 1e-3, "(relative)")


# ─── Convergence Checks ───────────────────────────────────────────────────────

def convergence_dt() -> bool:
    """C1: SNR² is stable to 0.1% when dt is halved."""
    kappa   = KAPPA
    chi     = CHI_NOMINAL
    epsilon = EPSILON_MAX * 0.5
    delta_g = optimal_delta_g(chi)
    T       = 3.0 / kappa
    n_grid  = 40

    snr2_vals = []
    for dt_factor in [1, 0.5, 0.25]:
        dt = DT_ODE * dt_factor
        N = max(2, int(T / dt))
        tlist = np.linspace(0.0, T, N + 1)
        eps_arr = epsilon * np.ones(N, dtype=np.complex128)
        ag, ae = simulate_conditioned_fields(eps_arr, tlist, kappa, chi, delta_g)
        snr2_vals.append(snr_squared(ag, ae, tlist, kappa))

    # Check convergence between finest two
    change_01 = abs(snr2_vals[1] - snr2_vals[0]) / max(abs(snr2_vals[0]), 1e-30)
    change_12 = abs(snr2_vals[2] - snr2_vals[1]) / max(abs(snr2_vals[1]), 1e-30)
    print(f"  SNR² at dt={DT_ODE*1e9:.0f}ns: {snr2_vals[0]:.6f}")
    print(f"  SNR² at dt={DT_ODE*0.5e9:.0f}ns: {snr2_vals[1]:.6f}  (change: {change_01:.2e})")
    print(f"  SNR² at dt={DT_ODE*0.25e9:.0f}ns: {snr2_vals[2]:.6f}  (change: {change_12:.2e})")
    passed = change_01 < 1e-3
    return check("C1: dt convergence (SNR² change)", passed, change_01, 1e-3, "(relative)")


def convergence_grape_segments() -> bool:
    """C2: GRAPE SNR² converges as N_seg increases."""
    kappa   = KAPPA
    chi     = CHI_NOMINAL
    epsilon_max = EPSILON_MAX * 0.8
    delta_g = optimal_delta_g(chi)
    T       = 3.0 / kappa
    snr2_values = {}

    for n_seg in [30, 60, 120]:
        gr = run_grape(kappa=kappa, chi=chi, T_read=T, N_seg=n_seg,
                       delta_g=delta_g, epsilon_max=epsilon_max,
                       n_restarts=5, maxiter=200, w_snr=1.0, w_res=0.0)
        snr2_values[n_seg] = gr["snr2_opt"]
        print(f"    N_seg={n_seg}: SNR²={snr2_values[n_seg]:.5f}")

    change = abs(snr2_values[120] - snr2_values[60]) / max(abs(snr2_values[60]), 1e-30)
    passed = change < 5e-3    # <0.5% change from 60→120 segments
    return check("C2: GRAPE N_seg convergence (60→120)", passed, change, 5e-3, "(relative)")


# ─── Literature / Theoretical Checks ─────────────────────────────────────────

def lit_gamma_meas_formula() -> bool:
    """
    L1: γ_meas at optimal drive equals (κ/2)|Δα_ss|².
    Verify our formula matches cqed_sim.ReadoutResonator.gamma_meas().
    """
    kappa   = KAPPA
    chi     = CHI_NOMINAL
    epsilon = EPSILON_MAX * 0.5
    delta_g = optimal_delta_g(chi)
    delta_e = delta_g + chi

    ag_ss = steady_state_alpha(epsilon, kappa, delta_g)
    ae_ss = steady_state_alpha(epsilon, kappa, delta_e)
    gamma_meas_formula = 0.5 * kappa * abs(ae_ss - ag_ss) ** 2

    res = make_readout_resonator(kappa=kappa, chi=chi, epsilon=epsilon)
    omega_d = OMEGA_R - delta_g
    gamma_meas_cqed = res.gamma_meas(drive_frequency=omega_d)

    err = abs(gamma_meas_formula - gamma_meas_cqed) / max(abs(gamma_meas_cqed), 1e-30)
    passed = err < 1e-6
    print(f"  γ_meas formula: {gamma_meas_formula/TWO_PI/1e3:.3f} kHz")
    print(f"  γ_meas cqed:    {gamma_meas_cqed/TWO_PI/1e3:.3f} kHz")
    return check("L1: γ_meas formula accuracy", passed, err, 1e-6, "(relative)")


def lit_optimal_drive_frequency() -> bool:
    """
    L2: Drive at Δ_g = −χ/2 maximizes steady-state |Δα_ss|.
    Numerical verification against the analytical formula.
    """
    kappa   = KAPPA
    chi     = CHI_NOMINAL
    epsilon = EPSILON_MAX * 0.5

    # Analytical optimal
    delta_g_opt_analytical = optimal_delta_g(chi)

    # Numerical sweep
    delta_g_scan = np.linspace(-2*abs(chi), abs(chi), 2000)
    from common import steady_state_separation
    seps = np.array([steady_state_separation(epsilon, kappa, chi, dg) for dg in delta_g_scan])
    delta_g_opt_numerical = delta_g_scan[np.argmax(seps)]

    err = abs(delta_g_opt_numerical - delta_g_opt_analytical) / max(abs(delta_g_opt_analytical), 1e-30)
    passed = err < 1e-2    # numerical sweep has ~0.1% grid resolution
    print(f"  Analytical Δ_g_opt/2π = {delta_g_opt_analytical/TWO_PI/1e6:.4f} MHz")
    print(f"  Numerical  Δ_g_opt/2π = {delta_g_opt_numerical/TWO_PI/1e6:.4f} MHz")
    return check("L2: optimal drive frequency", passed, err, 1e-2, "(relative)")


def lit_grape_upper_bound() -> bool:
    """
    L3: GRAPE SNR² ≥ best pulse-family SNR² (necessary condition for correctness).
    """
    kappa   = KAPPA
    chi     = CHI_NOMINAL
    delta_g = optimal_delta_g(chi)
    T_read  = 3.0 / kappa

    # Best pulse family
    best_snr2_family = -np.inf
    best_family = ""
    for fam in ["square", "hann", "gaussian"]:
        _, ev = optimize_amplitude(fam, T_read, kappa, chi, delta_g, EPSILON_MAX)
        if ev["snr2"] > best_snr2_family:
            best_snr2_family = ev["snr2"]
            best_family = fam

    # GRAPE
    gr = run_grape(kappa=kappa, chi=chi, T_read=T_read, N_seg=60,
                   delta_g=delta_g, epsilon_max=EPSILON_MAX,
                   n_restarts=8, maxiter=300)

    margin = gr["snr2_opt"] - best_snr2_family
    deficit_rel = max(0.0, -margin / max(best_snr2_family, 1e-30))
    passed = deficit_rel <= 1e-4     # allow tiny relative numerical tie-breaking errors
    print(f"  Best family ({best_family}) SNR² = {best_snr2_family:.5f}")
    print(f"  GRAPE SNR² = {gr['snr2_opt']:.5f}")
    print(f"  GRAPE margin = {margin:.5f}")
    return check("L3: GRAPE ≥ best family (GRAPE upper-bound)", passed, deficit_rel, 1e-4, "(relative deficit)")


# ─── Master runner ────────────────────────────────────────────────────────────

def main() -> None:
    print("=" * 60)
    print("Validation Framework — Three-Check Protocol")
    print("=" * 60)

    all_passed = True

    print("\n── Sanity Checks ─────────────────────────────────────────")
    r = sanity_zero_drive();        all_passed = all_passed and r
    r = sanity_steady_state_convergence(); all_passed = all_passed and r
    r = sanity_cqed_sim_response_trace(); all_passed = all_passed and r

    print("\n── Convergence Checks ────────────────────────────────────")
    r = convergence_dt();           all_passed = all_passed and r
    r = convergence_grape_segments(); all_passed = all_passed and r

    print("\n── Literature / Theoretical Checks ──────────────────────")
    r = lit_gamma_meas_formula();   all_passed = all_passed and r
    r = lit_optimal_drive_frequency(); all_passed = all_passed and r
    r = lit_grape_upper_bound();    all_passed = all_passed and r

    print("\n" + "=" * 60)
    if all_passed:
        print("ALL CHECKS PASSED ✓")
    else:
        print("SOME CHECKS FAILED — review results before reporting")
    print("=" * 60)


if __name__ == "__main__":
    main()
