"""
Phase 5 — Robust Optimization Study
=====================================

Goals:
  1. Evaluate sensitivity of the nominal-optimal pulse (square, Hann, GRAPE)
     to ±20% uncertainty in χ and κ.
  2. Implement robust GRAPE: optimize over an ensemble of (χ, κ) parameter samples.
  3. Compare nominal-optimal vs. robust-optimal on the full uncertainty distribution.
  4. Determine whether robust optimization is worth the additional cost.
  5. Study which pulse family is most naturally robust.

Robustness metric:
  For each pulse design, evaluate SNR² over a grid of (χ_err, κ_err) ∈ [−20%, +20%]²
  and report:
    - mean SNR² over the uncertainty distribution
    - minimum SNR² (worst-case)
    - 10th percentile SNR²
    - coefficient of variation (std/mean)

Robust GRAPE objective:
  L_robust = −(1/M) Σ_{m=1}^{M} SNR²(ε; λ_m)
  where λ_m = (χ_m, κ_m) are ensemble samples.

Outputs:
  data/phase5_results.npz
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

import numpy as np
from scipy.optimize import minimize
from common import (
    TWO_PI, KAPPA, CHI_NOMINAL, EPSILON_MAX,
    N_GRAPE_SEGMENTS, N_GRAPE_RESTARTS, GRAPE_MAXITER, DT_ODE,
    optimal_delta_g, run_grape, optimize_amplitude,
    snr_squared, assignment_fidelity_from_snr2, residual_photons,
    simulate_conditioned_fields, evaluate_pulse_family,
    grape_propagators, grape_forward_pass, grape_backward_pass,
    PULSE_FAMILIES,
)

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
os.makedirs(DATA_DIR, exist_ok=True)

KAPPA_T_MAIN = 5.0
N_SEG = N_GRAPE_SEGMENTS
UNCERTAINTY_FRAC = 0.20    # ±20% uncertainty in χ and κ


# ─── Sensitivity analysis ─────────────────────────────────────────────────────

def evaluate_pulse_under_errors(
    epsilon_t: np.ndarray,
    T_read: float,
    kappa_nom: float,
    chi_nom: float,
    delta_g_nom: float,
    chi_err_frac: float,
    kappa_err_frac: float,
    dt: float = DT_ODE,
) -> float:
    """
    Evaluate SNR² for a given pulse under perturbed parameters.

    The drive frequency (and hence delta_g) was designed for the nominal parameters.
    Under errors, the actual detuning shifts, reducing performance.

    Args:
        epsilon_t     : complex pulse envelope, length (N_pts,) or (N_seg,)
        T_read        : readout duration (s)
        kappa_nom     : nominal κ (rad/s)
        chi_nom       : nominal χ (rad/s)
        delta_g_nom   : Δ_g used when designing the pulse
        chi_err_frac  : fractional error in χ, i.e., χ_actual = χ_nom · (1 + χ_err_frac)
        kappa_err_frac: fractional error in κ

    Returns:
        SNR² under perturbed parameters
    """
    kappa_actual = kappa_nom * (1.0 + kappa_err_frac)
    chi_actual   = chi_nom   * (1.0 + chi_err_frac)
    # Drive frequency was optimized for nominal δg; actual δg shifts
    # Δ_g = ω_r − ω_d; since ω_d = ω_r − δg_nom, actual Δ_g = ω_r − ω_d = δg_nom
    # but ω_r(e) = ω_r + χ_actual so Δ_e = Δ_g_nom + χ_actual
    delta_g_actual = delta_g_nom    # drive frequency unchanged (hardware calibration)

    N = len(epsilon_t)
    if N < 2:
        N = max(2, int(round(T_read / dt)))
    tlist = np.linspace(0.0, T_read, N + 1)
    ag, ae = simulate_conditioned_fields(
        epsilon_t, tlist, kappa_actual, chi_actual, delta_g_actual
    )
    return snr_squared(ag, ae, tlist, kappa_actual)


def sensitivity_map_2d(
    epsilon_t: np.ndarray,
    T_read: float,
    kappa_nom: float,
    chi_nom: float,
    delta_g_nom: float,
    n_grid: int = 15,
    uncertainty: float = UNCERTAINTY_FRAC,
    dt: float = DT_ODE,
) -> dict:
    """
    2-D sensitivity map over (χ_err, κ_err) ∈ [−uncertainty, +uncertainty]².

    Returns:
        chi_err_grid, kappa_err_grid, snr2_map (2D array)
    """
    chi_errs   = np.linspace(-uncertainty, uncertainty, n_grid)
    kappa_errs = np.linspace(-uncertainty, uncertainty, n_grid)
    snr2_map = np.zeros((n_grid, n_grid))

    for i, ce in enumerate(chi_errs):
        for j, ke in enumerate(kappa_errs):
            snr2_map[i, j] = evaluate_pulse_under_errors(
                epsilon_t, T_read, kappa_nom, chi_nom, delta_g_nom, ce, ke, dt
            )
    return {
        "chi_err_grid":   chi_errs,
        "kappa_err_grid": kappa_errs,
        "snr2_map": snr2_map,
    }


def sensitivity_stats(snr2_map: np.ndarray, snr2_nominal: float) -> dict:
    """
    Compute robustness statistics from a 2D sensitivity map.
    """
    flat = snr2_map.ravel()
    return {
        "mean":    float(np.mean(flat)),
        "std":     float(np.std(flat)),
        "min":     float(np.min(flat)),
        "p10":     float(np.percentile(flat, 10)),
        "cv":      float(np.std(flat) / max(np.mean(flat), 1e-30)),
        "nominal": float(snr2_nominal),
        "mean_loss_frac": float((snr2_nominal - np.mean(flat)) / max(snr2_nominal, 1e-30)),
    }


# ─── Robust GRAPE ─────────────────────────────────────────────────────────────

def robust_grape_objective_and_grad(
    x: np.ndarray,
    kappa_nom: float,
    chi_nom: float,
    delta_g_nom: float,
    tau: float,
    N: int,
    ensemble: list[tuple[float, float]],
) -> tuple[float, np.ndarray]:
    """
    Robust GRAPE: objective = −mean(SNR²) over (χ, κ) ensemble.

    ensemble : list of (chi_actual, kappa_actual) tuples
    """
    epsilon_vec = x[:N] + 1j * x[N:]
    total_obj  = 0.0
    total_grad = np.zeros(2 * N)
    M = len(ensemble)

    for (chi_m, kappa_m) in ensemble:
        delta_e_m = delta_g_nom + chi_m
        A_g, B_g, A_e, B_e = grape_propagators(kappa_m, delta_g_nom, delta_e_m, tau)
        ag_m, ae_m = grape_forward_pass(epsilon_vec, kappa_m, delta_g_nom, delta_e_m, tau)
        tlist_m = np.arange(N + 1) * tau
        snr2_m = snr_squared(ag_m, ae_m, tlist_m, kappa_m)

        # Gradient for this ensemble member
        grad_w = grape_backward_pass(ag_m, ae_m, kappa_m, tau, A_g, A_e, B_g, B_e)
        grad_wirtinger_m = -1.0 * kappa_m * tau * (B_e * np.array([
            np.conj(ae_m[n+1] - ag_m[n+1]) if n < N else 0.0 for n in range(N)
        ]))
        # Use the proper backward pass result
        grad_wirtinger_m = grape_backward_pass(ag_m, ae_m, kappa_m, tau, A_g, A_e, B_g, B_e)
        grad_wirtinger_m = -grad_wirtinger_m   # negative for maximization

        grad_real_m = np.concatenate([
            2.0 * np.real(grad_wirtinger_m),
            2.0 * np.imag(grad_wirtinger_m),
        ])
        total_obj  += -snr2_m      # minimize −SNR²
        total_grad += grad_real_m

    return float(total_obj / M), total_grad / M


def run_robust_grape(
    kappa_nom: float = KAPPA,
    chi_nom: float = None,
    kappa_t: float = KAPPA_T_MAIN,
    epsilon_max: float = EPSILON_MAX,
    n_seg: int = N_SEG,
    n_restarts: int = N_GRAPE_RESTARTS,
    maxiter: int = GRAPE_MAXITER,
    uncertainty: float = UNCERTAINTY_FRAC,
    n_ensemble: int = 5,
    seed: int = 42,
) -> dict:
    """
    Run robust GRAPE using a Sobol-like grid ensemble of (χ, κ) perturbations.

    Returns the robust-optimal pulse and its nominal and mean SNR².
    """
    if chi_nom is None:
        chi_nom = CHI_NOMINAL
    delta_g_nom = optimal_delta_g(chi_nom)
    T_read = kappa_t / kappa_nom
    tau = T_read / n_seg
    bounds = [(-epsilon_max, epsilon_max)] * (2 * n_seg)

    # Build ensemble: uniform grid in (χ_err, κ_err) ∈ [−unc, +unc]
    errs = np.linspace(-uncertainty, uncertainty, n_ensemble)
    ensemble = []
    for ce in errs:
        for ke in errs:
            ensemble.append((chi_nom * (1.0 + ce), kappa_nom * (1.0 + ke)))

    rng = np.random.default_rng(seed)
    best_obj = np.inf
    best_x = None

    for _ in range(n_restarts):
        x0 = rng.uniform(-epsilon_max, epsilon_max, size=2 * n_seg)
        result = minimize(
            robust_grape_objective_and_grad,
            x0,
            jac=True,
            method="L-BFGS-B",
            bounds=bounds,
            args=(kappa_nom, chi_nom, delta_g_nom, tau, n_seg, ensemble),
            options={"maxiter": maxiter, "ftol": 1e-12, "gtol": 1e-8},
        )
        if result.fun < best_obj:
            best_obj = result.fun
            best_x = result.x

    eps_robust = best_x[:n_seg] + 1j * best_x[n_seg:]

    # Evaluate on nominal parameters
    tlist = np.linspace(0.0, T_read, n_seg + 1)
    delta_e_nom = delta_g_nom + chi_nom
    ag_nom, ae_nom = grape_forward_pass(eps_robust, kappa_nom, delta_g_nom, delta_e_nom, tau)
    snr2_nom = snr_squared(ag_nom, ae_nom, tlist, kappa_nom)

    return {
        "epsilon_robust": eps_robust,
        "snr2_nominal": snr2_nom,
        "F_assign_nominal": assignment_fidelity_from_snr2(snr2_nom),
        "n_res_nominal": residual_photons(ag_nom, ae_nom),
        "alpha_g": ag_nom,
        "alpha_e": ae_nom,
        "tlist": tlist,
        "ensemble_size": len(ensemble),
        "kappa_t": kappa_t,
        "n_seg": n_seg,
    }


def compare_nominal_vs_robust(
    kappa: float = KAPPA,
    chi: float = None,
    kappa_t: float = KAPPA_T_MAIN,
    epsilon_max: float = EPSILON_MAX,
    n_restarts: int = N_GRAPE_RESTARTS,
    n_grid_sens: int = 15,
    uncertainty: float = UNCERTAINTY_FRAC,
) -> dict:
    """
    Full comparison:
      - Nominal-optimal GRAPE (square, Hann, GRAPE)
      - Robust GRAPE
      Each evaluated on 2D sensitivity map.
    """
    if chi is None:
        chi = CHI_NOMINAL
    delta_g = optimal_delta_g(chi)
    T_read = kappa_t / kappa

    # Nominal GRAPE
    print("    Nominal GRAPE...")
    gr_nom = run_grape(kappa=kappa, chi=chi, T_read=T_read, N_seg=N_SEG,
                       delta_g=delta_g, epsilon_max=epsilon_max,
                       n_restarts=n_restarts, w_snr=1.0, w_res=0.0)

    # Robust GRAPE
    print("    Robust GRAPE (5×5 ensemble)...")
    gr_rob = run_robust_grape(kappa_nom=kappa, chi_nom=chi, kappa_t=kappa_t,
                              epsilon_max=epsilon_max, n_restarts=n_restarts,
                              n_ensemble=5)

    # Best pulse-family (square, optimized amplitude)
    print("    Best square pulse...")
    _, ev_sq = optimize_amplitude("square", T_read, kappa, chi, delta_g, epsilon_max)

    # Best Hann pulse
    _, ev_hann = optimize_amplitude("hann", T_read, kappa, chi, delta_g, epsilon_max)

    # Sensitivity maps for each
    N = N_SEG
    dt_grape = T_read / N
    tlist_grape = np.linspace(0.0, T_read, N + 1)

    # For pulse families, resample to dense grid
    dt_dense = DT_ODE
    N_dense = max(2, int(round(T_read / dt_dense)))

    print("    Computing sensitivity maps...")
    results = {}

    for label, epsilon_t, snr2_nom in [
        ("GRAPE_nominal", gr_nom["epsilon_opt"], gr_nom["snr2_opt"]),
        ("GRAPE_robust",  gr_rob["epsilon_robust"], gr_rob["snr2_nominal"]),
        ("square",        ev_sq["epsilon_t"],    ev_sq["snr2"]),
        ("hann",          ev_hann["epsilon_t"],  ev_hann["snr2"]),
    ]:
        print(f"      {label}...")
        smap = sensitivity_map_2d(
            epsilon_t, T_read, kappa, chi, delta_g,
            n_grid=n_grid_sens, uncertainty=uncertainty
        )
        stats = sensitivity_stats(smap["snr2_map"], snr2_nom)
        results[label] = {
            "snr2_map": smap["snr2_map"],
            "stats": stats,
            "snr2_nominal": snr2_nom,
        }
    results["chi_err_grid"]   = smap["chi_err_grid"]
    results["kappa_err_grid"] = smap["kappa_err_grid"]

    return {
        "comparison": results,
        "eps_grape_nominal": gr_nom["epsilon_opt"],
        "eps_grape_robust": gr_rob["epsilon_robust"],
        "eps_square": ev_sq["epsilon_t"],
        "eps_hann": ev_hann["epsilon_t"],
        "tlist": gr_nom["tlist"],
        "chi_err_grid": smap["chi_err_grid"],
        "kappa_err_grid": smap["kappa_err_grid"],
        "kappa_t": kappa_t,
    }


def main() -> None:
    print("=" * 60)
    print("Phase 5: Robust Optimization Study")
    print("=" * 60)

    kappa = KAPPA
    chi = CHI_NOMINAL

    print("\n[1/1] Nominal vs. robust comparison (all pulse types)...")
    comp = compare_nominal_vs_robust(kappa=kappa, chi=chi, n_restarts=N_GRAPE_RESTARTS)

    print("\n  Results summary:")
    print(f"  {'Design':20s} {'Nominal SNR²':12s} {'Mean SNR²':10s} {'Min SNR²':10s} {'CV':8s} {'Loss%':8s}")
    print(f"  {'-'*72}")
    for label in ["square", "hann", "GRAPE_nominal", "GRAPE_robust"]:
        if label in comp["comparison"]:
            st = comp["comparison"][label]["stats"]
            print(f"  {label:20s} {st['nominal']:12.3f} {st['mean']:10.3f} "
                  f"{st['min']:10.3f} {st['cv']:8.4f} {100*st['mean_loss_frac']:8.2f}%")

    # ── Save ─────────────────────────────────────────────────────────────────
    outpath = os.path.join(DATA_DIR, "phase5_results.npz")
    save_dict = {
        "chi_err_grid":   comp["chi_err_grid"],
        "kappa_err_grid": comp["kappa_err_grid"],
        "kappa_t":        comp["kappa_t"],
    }
    for label in ["GRAPE_nominal", "GRAPE_robust", "square", "hann"]:
        if label in comp["comparison"]:
            d = comp["comparison"][label]
            prefix = f"{label}_"
            save_dict[prefix + "snr2_map"]  = d["snr2_map"]
            save_dict[prefix + "snr2_nom"]  = d["stats"]["nominal"]
            save_dict[prefix + "mean_snr2"] = d["stats"]["mean"]
            save_dict[prefix + "min_snr2"]  = d["stats"]["min"]
            save_dict[prefix + "p10_snr2"]  = d["stats"]["p10"]
            save_dict[prefix + "cv"]        = d["stats"]["cv"]
            save_dict[prefix + "mean_loss_frac"] = d["stats"]["mean_loss_frac"]
    save_dict["eps_grape_nominal"] = comp["eps_grape_nominal"]
    save_dict["eps_grape_robust"]  = comp["eps_grape_robust"]
    save_dict["tlist"]             = comp["tlist"]
    save_dict["kappa_nom"]         = kappa
    save_dict["chi_nom"]           = chi

    np.savez(outpath, **save_dict)
    print(f"\nResults saved to {outpath}")

    print("\n" + "=" * 60)
    print("Phase 5 Summary:")
    for label in ["square", "hann", "GRAPE_nominal", "GRAPE_robust"]:
        if label in comp["comparison"]:
            st = comp["comparison"][label]["stats"]
            print(f"  {label}: CV={st['cv']:.4f}, mean loss={100*st['mean_loss_frac']:.1f}%")
    print("=" * 60)


if __name__ == "__main__":
    main()
