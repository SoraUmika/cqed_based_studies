"""
probe_library.py — Chi Ramsey probe generation and chi inference.

This module implements the multi-Fock dispersive Ramsey probe and the chi
inference procedure from noisy probe data.

API Gap Note
-----------
The multi-Fock dispersive Ramsey probe is NOT in cqed_sim.calibration_targets.
run_ramsey in cqed_sim provides single-frequency surrogate probing but does not
support the Fock-number-resolved version needed for chi_higher inference.
We implement the exact analytical forward model directly (valid in the dispersive
RWA limit, after ideal pi/2 pulses):

    P_e(t, n) = 0.5 * (1 + cos(chi_eff(n) * t) * exp(-t / T2*))
    chi_eff(n) = chi * n + chi_higher * n * (n - 1)

where n is the Fock number. SPAM is applied via the confusion matrix M:

    P_obs = M[1,0] + (M[1,1] - M[1,0]) * P_e

Shot noise is sampled binomially: counts_e / n_shots.

Physical sign convention: chi < 0 for typical cQED (qubit frequency redshifted
by photons). chi_eff(n) < 0 for n > 0, so cos(chi_eff(n) * t) oscillates at
|chi_eff(n)| rad/s. The inference is sign-sensitive.

Chi inference uses scipy.optimize.minimize (L-BFGS-B) to fit [chi, chi_higher,
log_T2, amplitude corrections A_n] from multi-Fock data simultaneously.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

from cqed_sim.sim.noise import NoiseSpec


# ---------------------------------------------------------------------------
# T2* from NoiseSpec
# ---------------------------------------------------------------------------


def t2_star_from_noise(noise_spec: NoiseSpec) -> float:
    """
    Compute the effective T2* dephasing time from a NoiseSpec.

    Uses the standard formula:
        1/T2* = 1/(2*T1) + 1/Tphi

    where T1 is the energy relaxation time and Tphi is the pure dephasing time.

    Parameters
    ----------
    noise_spec : NoiseSpec
        Noise specification with t1 and tphi fields.

    Returns
    -------
    float
        T2* in seconds. Returns a large value (1 s) if no dephasing is set.
    """
    rate = 0.0
    if noise_spec.t1 is not None and noise_spec.t1 > 0.0:
        rate += 1.0 / (2.0 * float(noise_spec.t1))
    if noise_spec.tphi is not None and noise_spec.tphi > 0.0:
        rate += 1.0 / float(noise_spec.tphi)
    if rate <= 0.0:
        return 1.0  # effectively infinite coherence
    return 1.0 / rate


# ---------------------------------------------------------------------------
# Analytical forward model
# ---------------------------------------------------------------------------


def _chi_eff(n: int, chi: float, chi_higher: float) -> float:
    """
    Compute the effective Ramsey oscillation frequency for Fock state n.

    chi_eff(n) = chi * n + chi_higher * n * (n - 1)

    For the standard dispersive Hamiltonian H_disp = chi * n_cav * sigma_z / 2,
    the qubit frequency shift at n photons is chi_eff(n). chi < 0 is typical.

    Parameters
    ----------
    n : int
        Fock number (photon count in cavity).
    chi : float
        Linear dispersive coupling (rad/s). Negative is typical.
    chi_higher : float
        Second-order coefficient (rad/s). chi_higher * n*(n-1) correction.

    Returns
    -------
    float
        Effective oscillation frequency (rad/s) at Fock state n.
    """
    return float(chi) * int(n) + float(chi_higher) * int(n) * (int(n) - 1)


def _p_e_ideal(delays: np.ndarray, n: int, chi: float, chi_higher: float, t2_star: float) -> np.ndarray:
    """
    Compute the ideal (noiseless) excited-state probability for Fock state n.

    P_e(t, n) = 0.5 * (1 + cos(chi_eff(n) * t) * exp(-t / T2*))

    For n=0: chi_eff(0) = 0, so P_e(t, 0) = 0.5 * (1 + exp(-t/T2*)).
    This is the pure dephasing decay without oscillation.

    Parameters
    ----------
    delays : np.ndarray
        Array of delay times in seconds.
    n : int
        Fock number.
    chi : float
        Linear dispersive coupling (rad/s).
    chi_higher : float
        Second-order chi coefficient (rad/s).
    t2_star : float
        T2* coherence time in seconds.

    Returns
    -------
    np.ndarray
        P_e values at each delay time.
    """
    t = np.asarray(delays, dtype=float)
    freq = _chi_eff(n, chi, chi_higher)
    return 0.5 * (1.0 + np.cos(freq * t) * np.exp(-t / float(t2_star)))


def _apply_spam(p_e: np.ndarray, confusion_matrix: np.ndarray) -> np.ndarray:
    """
    Apply SPAM (State Preparation And Measurement) via the confusion matrix.

    M[i,j] = P(observe state i | true state j)
    P_obs_e = M[1,0] * P_g + M[1,1] * P_e
            = M[1,0] * (1 - P_e) + M[1,1] * P_e
            = M[1,0] + (M[1,1] - M[1,0]) * P_e

    Parameters
    ----------
    p_e : np.ndarray
        True excited-state probability (before SPAM).
    confusion_matrix : np.ndarray, shape (2, 2)
        Readout confusion matrix. M[observed, true].

    Returns
    -------
    np.ndarray
        Observed P_e after SPAM corruption.
    """
    M = np.asarray(confusion_matrix, dtype=float)
    return M[1, 0] + (M[1, 1] - M[1, 0]) * np.asarray(p_e, dtype=float)


# ---------------------------------------------------------------------------
# Probe generation
# ---------------------------------------------------------------------------


def run_chi_ramsey_probe(
    chi_true: float,
    chi_higher_true: float,
    t2_star: float,
    confusion_matrix: np.ndarray,
    n_shots: int,
    delays_s: np.ndarray,
    fock_levels: list[int] | None = None,
    rng: np.random.Generator | None = None,
) -> dict:
    """
    Simulate a multi-Fock dispersive Ramsey experiment.

    For each Fock level n in fock_levels:
    1. Compute ideal P_e(t, n) using the true chi, chi_higher, and T2*.
    2. Apply SPAM via confusion_matrix.
    3. Sample binomially to get noisy shot counts.

    Note: This uses the analytical forward model, not cqed_sim.simulate_sequence.
    This is justified because:
    - The probe uses ideal pi/2 pulses (fast compared to chi evolution).
    - The dispersive RWA is exact for the qubit oscillation at the level of the
      chi measurement.
    - Adding a full Hamiltonian simulation would require specifying pi/2 pulse
      details that are not part of the chi inference setup.

    Parameters
    ----------
    chi_true : float
        True dispersive coupling (rad/s).
    chi_higher_true : float
        True second-order chi (rad/s).
    t2_star : float
        T2* coherence time (seconds).
    confusion_matrix : np.ndarray, shape (2, 2)
        SPAM confusion matrix. M[observed, true].
    n_shots : int
        Number of measurement shots per time point per Fock level.
    delays_s : np.ndarray
        Array of delay times in seconds.
    fock_levels : list[int], optional
        Fock levels to probe. Default: [1, 2, 3]. (n=0 has no chi oscillation.)
    rng : np.random.Generator, optional
        Random number generator. Created fresh if None.

    Returns
    -------
    dict with keys:
        'delays'       : np.ndarray of delay times (s)
        'n_list'       : list of Fock levels probed
        'p_e_observed' : np.ndarray of shape (n_fock, n_delays)
                         Noisy observed excited-state probability
        'p_e_ideal'    : np.ndarray of shape (n_fock, n_delays)
                         Noiseless ideal P_e (for reference / debug)
        'p_e_spam'     : np.ndarray of shape (n_fock, n_delays)
                         P_e after SPAM but before shot noise
        't2_star'      : float, T2* used
        'chi_true'     : float, chi_true used
        'chi_higher_true' : float, chi_higher_true used
    """
    if fock_levels is None:
        fock_levels = [1, 2, 3]
    if rng is None:
        rng = np.random.default_rng()

    delays = np.asarray(delays_s, dtype=float)
    n_delays = len(delays)
    n_fock = len(fock_levels)

    p_e_ideal = np.zeros((n_fock, n_delays), dtype=float)
    p_e_spam = np.zeros((n_fock, n_delays), dtype=float)
    p_e_observed = np.zeros((n_fock, n_delays), dtype=float)

    for i, n in enumerate(fock_levels):
        p_ideal = _p_e_ideal(delays, n, chi_true, chi_higher_true, t2_star)
        p_spam = _apply_spam(p_ideal, confusion_matrix)
        # Clip to valid probability range before binomial sampling
        p_spam_clipped = np.clip(p_spam, 0.0, 1.0)
        counts = rng.binomial(int(n_shots), p_spam_clipped)
        p_obs = counts.astype(float) / float(n_shots)

        p_e_ideal[i] = p_ideal
        p_e_spam[i] = p_spam
        p_e_observed[i] = p_obs

    return {
        "delays": delays,
        "n_list": list(fock_levels),
        "p_e_observed": p_e_observed,
        "p_e_ideal": p_e_ideal,
        "p_e_spam": p_e_spam,
        "t2_star": float(t2_star),
        "chi_true": float(chi_true),
        "chi_higher_true": float(chi_higher_true),
    }


# ---------------------------------------------------------------------------
# Chi inference
# ---------------------------------------------------------------------------


def _forward_model_spam(
    delays: np.ndarray,
    n: int,
    chi: float,
    chi_higher: float,
    t2: float,
    amplitude: float,
    confusion_matrix: np.ndarray,
) -> np.ndarray:
    """
    Forward model for P_e_obs with amplitude correction and SPAM.

    P_e_obs(t, n) = M[1,0] + (M[1,1] - M[1,0]) * 0.5 * (1 + A * cos(chi_eff(n)*t) * exp(-t/T2))

    The amplitude parameter A accounts for state preparation imperfections
    (e.g., imperfect pi/2 pulse), separately from the confusion matrix. Typically A ~ 1.

    Parameters
    ----------
    delays : np.ndarray
        Delay times in seconds.
    n : int
        Fock number.
    chi : float
        Dispersive coupling (rad/s).
    chi_higher : float
        Second-order chi (rad/s).
    t2 : float
        T2* (seconds). Must be positive.
    amplitude : float
        Contrast amplitude (0 to 1). Accounts for imperfect state prep.
    confusion_matrix : np.ndarray
        SPAM confusion matrix.

    Returns
    -------
    np.ndarray
        Predicted observed P_e at each delay.
    """
    t = np.asarray(delays, dtype=float)
    t2 = max(float(t2), 1e-12)
    freq = _chi_eff(n, chi, chi_higher)
    p_e = 0.5 * (1.0 + float(amplitude) * np.cos(freq * t) * np.exp(-t / t2))
    p_e = np.clip(p_e, 0.0, 1.0)
    return _apply_spam(p_e, confusion_matrix)


def _build_residuals(
    params: np.ndarray,
    delays: np.ndarray,
    fock_levels: list[int],
    p_e_observed: np.ndarray,
    confusion_matrix: np.ndarray,
    n_shots: int,
) -> float:
    """
    Chi-squared cost function for the chi inference fit.

    Parameters
    ----------
    params : np.ndarray
        Fit parameters: [chi, chi_higher, log_T2, A_1, A_2, ..., A_n_fock]
        where A_i is the amplitude correction for fock_levels[i].
    delays : np.ndarray
        Delay times (s).
    fock_levels : list[int]
        Fock levels probed.
    p_e_observed : np.ndarray, shape (n_fock, n_delays)
        Observed P_e data.
    confusion_matrix : np.ndarray
        SPAM confusion matrix.
    n_shots : int
        Number of shots per data point (for sigma^2 estimation).

    Returns
    -------
    float
        Sum of squared residuals (weighted by shot noise variance).
    """
    chi = params[0]
    chi_higher = params[1]
    log_t2 = params[2]
    t2 = np.exp(log_t2)
    amplitudes = params[3:]

    total = 0.0
    for i, n in enumerate(fock_levels):
        a = float(amplitudes[i]) if i < len(amplitudes) else 1.0
        p_model = _forward_model_spam(delays, n, chi, chi_higher, t2, a, confusion_matrix)
        p_obs = p_e_observed[i]
        # Binomial variance: sigma^2 = p*(1-p)/n_shots (avoid division by zero)
        sigma2 = np.clip(p_model * (1.0 - p_model), 1e-6, None) / float(n_shots)
        residuals = (p_obs - p_model) ** 2 / sigma2
        total += float(np.mean(residuals))

    return total


def infer_chi_from_probe(
    probe_data: dict,
    confusion_matrix: np.ndarray,
    n_shots: int,
    chi_initial: float,
    chi_higher_initial: float = 0.0,
) -> dict:
    """
    Infer chi and chi_higher from multi-Fock Ramsey probe data.

    Uses scipy.optimize.minimize (L-BFGS-B) to fit the forward model to noisy
    observed P_e data from multiple Fock levels simultaneously.

    Fit parameters (scaled for numerical stability):
        [chi_MHz, chi_higher_kHz, log_T2, A_n_1, A_n_2, ..., A_n_k]

    where chi_MHz = chi / (2*pi*1e6) and chi_higher_kHz = chi_higher / (2*pi*1e3),
    and A_n_i is a per-Fock amplitude correction.

    Scaling is essential for L-BFGS-B convergence: chi in rad/s is ~1e7 which
    creates ill-conditioned gradient steps. Scaling to MHz units ensures the
    parameter magnitudes are O(1-10) throughout optimization.

    Uncertainty estimation uses the diagonal Hessian approximation from the
    L-BFGS-B fit: sigma_i ~ 1/sqrt(H_ii) where H_ii is estimated numerically.

    Parameters
    ----------
    probe_data : dict
        Output from run_chi_ramsey_probe. Must contain keys:
        'delays', 'n_list', 'p_e_observed', 't2_star'.
    confusion_matrix : np.ndarray
        SPAM confusion matrix used in the probe.
    n_shots : int
        Number of shots per data point.
    chi_initial : float
        Initial guess for chi (rad/s). Should be the learner's prior.
    chi_higher_initial : float
        Initial guess for chi_higher (rad/s). Default 0.

    Returns
    -------
    dict with keys:
        'chi'                : float, inferred chi (rad/s)
        'chi_higher'         : float, inferred chi_higher (rad/s)
        'T2_star'            : float, inferred T2* (s)
        'amplitudes'         : np.ndarray, per-Fock amplitude corrections
        'chi_uncertainty'    : float, 1-sigma uncertainty on chi (rad/s)
        'chi_higher_uncertainty' : float, 1-sigma uncertainty on chi_higher (rad/s)
        'fit_cost'           : float, final cost function value
        'fit_success'        : bool, whether optimization converged
        'fit_message'        : str, optimizer message
        'n_fev'              : int, number of function evaluations
    """
    delays = np.asarray(probe_data["delays"], dtype=float)
    fock_levels = list(probe_data["n_list"])
    p_e_observed = np.asarray(probe_data["p_e_observed"], dtype=float)
    t2_initial = float(probe_data.get("t2_star", 50e-6))
    n_fock = len(fock_levels)

    # Scaling constants: work in (MHz, kHz) internally so parameters are O(1-10)
    # This is critical for L-BFGS-B convergence since chi in rad/s is ~1e7.
    CHI_SCALE = 2.0 * np.pi * 1e6     # 1 unit = 1 MHz (in rad/s)
    CHI_H_SCALE = 2.0 * np.pi * 1e3   # 1 unit = 1 kHz (in rad/s)

    # Cost function with scaled parameters [chi_MHz, chi_h_kHz, log_T2, A1,..,An]
    def _residuals_scaled(params_scaled):
        chi = params_scaled[0] * CHI_SCALE
        chi_h = params_scaled[1] * CHI_H_SCALE
        # params_scaled[2] is log_T2, passed directly to _build_residuals which expects log_T2
        amplitudes = params_scaled[3:]
        return _build_residuals(
            np.concatenate([[chi, chi_h, params_scaled[2]], amplitudes]),
            delays, fock_levels, p_e_observed, confusion_matrix, n_shots
        )

    # Bounds in scaled (MHz / kHz) units — covers ±15 MHz in chi
    chi_lo_scaled = -15.0    # -15 MHz lower bound
    chi_hi_scaled = 1.0      # +1 MHz upper bound (small positive chi allowed)
    chi_higher_range_scaled = 500.0  # ±500 kHz for chi_higher
    bounds_scaled = (
        [(chi_lo_scaled, chi_hi_scaled),
         (-chi_higher_range_scaled, chi_higher_range_scaled),
         (np.log(1e-9), np.log(1e-3))]  # T2: 1 ns to 1 ms
        + [(0.05, 1.1)] * n_fock
    )

    # Two-stage strategy: coarse 1D chi scan + L-BFGS-B refinement.
    # Rationale: the chi Ramsey cost landscape has many local minima (one per oscillation
    # period). L-BFGS-B from chi_initial can get stuck in a nearby local minimum that is
    # not the global one. A coarse scan identifies the global minimum basin, which is then
    # refined by L-BFGS-B. The scan costs ~120 evaluations with T2 fixed.
    n_scan = 120
    chi_scan_scaled = np.linspace(chi_lo_scaled, chi_hi_scaled, n_scan)
    t2_log_init = np.log(max(t2_initial, 1e-9))

    best_cost_scan = float("inf")
    best_chi_scan_scaled = float(chi_initial) / CHI_SCALE

    for chi_s in chi_scan_scaled:
        params_scan = np.concatenate([
            [chi_s, 0.0, t2_log_init],
            np.ones(n_fock) * 0.95,
        ])
        c = _residuals_scaled(params_scan)
        if c < best_cost_scan:
            best_cost_scan = c
            best_chi_scan_scaled = float(chi_s)

    # Refinement starting from best scan point
    x0 = np.concatenate([
        [best_chi_scan_scaled,
         float(chi_higher_initial) / CHI_H_SCALE,
         t2_log_init],
        np.ones(n_fock) * 0.95,
    ])

    result = minimize(
        _residuals_scaled,
        x0,
        method="L-BFGS-B",
        bounds=bounds_scaled,
        options={"maxiter": 10000, "ftol": 1e-14, "gtol": 1e-10},
    )

    chi_hat = float(result.x[0]) * CHI_SCALE
    chi_higher_hat = float(result.x[1]) * CHI_H_SCALE
    t2_hat = float(np.exp(result.x[2]))
    amplitudes_hat = result.x[3:]
    f0 = float(result.fun)
    total_nfev = n_scan + int(result.nfev)

    # Uncertainty estimation via finite-difference Hessian diagonal (in scaled units, then convert)
    eps_chi_s = 0.01   # 0.01 MHz
    eps_chi_h_s = 0.1  # 0.1 kHz

    def sigma_from_hessian_diag_scaled(x_opt, index, eps):
        x_plus = x_opt.copy()
        x_minus = x_opt.copy()
        x_plus[index] += eps
        x_minus[index] -= eps
        f_plus = _residuals_scaled(x_plus)
        f_minus = _residuals_scaled(x_minus)
        d2f = (f_plus - 2.0 * f0 + f_minus) / (eps ** 2)
        if d2f <= 1e-30:
            return float("inf")
        return 1.0 / np.sqrt(d2f)

    sigma_chi_scaled = sigma_from_hessian_diag_scaled(result.x, 0, eps_chi_s)
    sigma_chi_h_scaled = sigma_from_hessian_diag_scaled(result.x, 1, eps_chi_h_s)

    return {
        "chi": chi_hat,
        "chi_higher": chi_higher_hat,
        "T2_star": t2_hat,
        "amplitudes": amplitudes_hat,
        "chi_uncertainty": sigma_chi_scaled * CHI_SCALE,
        "chi_higher_uncertainty": sigma_chi_h_scaled * CHI_H_SCALE,
        "fit_cost": f0,
        "fit_success": bool(result.success),
        "fit_message": str(result.message),
        "n_fev": total_nfev,
        "best_chi_scan_MHz": best_chi_scan_scaled,
    }


# ---------------------------------------------------------------------------
# Convenience: delay array builder
# ---------------------------------------------------------------------------


def make_ramsey_delays(
    chi_estimate: float,
    n_periods: float = 3.0,
    n_points: int = 60,
    t_min_s: float = 0.0,
) -> np.ndarray:
    """
    Build a delay array that covers n_periods of chi Ramsey oscillations.

    The characteristic period is T = 2*pi / |chi_estimate|.

    Parameters
    ----------
    chi_estimate : float
        Estimated chi (rad/s). Used to set the oscillation period.
    n_periods : float
        Number of periods to cover in the delay array.
    n_points : int
        Number of delay points.
    t_min_s : float
        Minimum delay time (s). Default 0.

    Returns
    -------
    np.ndarray
        Delay times in seconds.
    """
    period = 2.0 * np.pi / max(abs(float(chi_estimate)), 1.0)
    t_max = float(n_periods) * period
    return np.linspace(float(t_min_s), t_max, int(n_points))


if __name__ == "__main__":
    import sys
    sys.path.insert(0, ".")
    from models import CHI_TRUE, CHI_HIGHER_TRUE, CHI_PRIOR, CONFUSION_MATRIX, N_SHOTS_PROBE, TRUTH_NOISE

    rng = np.random.default_rng(42)
    t2_star = t2_star_from_noise(TRUTH_NOISE)
    print(f"T2* = {t2_star * 1e6:.2f} us")

    delays = make_ramsey_delays(CHI_TRUE, n_periods=3.0, n_points=80)
    print(f"Delay range: 0 to {delays[-1]*1e6:.2f} us, {len(delays)} points")

    probe_data = run_chi_ramsey_probe(
        chi_true=CHI_TRUE,
        chi_higher_true=CHI_HIGHER_TRUE,
        t2_star=t2_star,
        confusion_matrix=CONFUSION_MATRIX,
        n_shots=N_SHOTS_PROBE,
        delays_s=delays,
        fock_levels=[1, 2, 3],
        rng=rng,
    )
    print(f"Probe data shape: {probe_data['p_e_observed'].shape}")

    infer_result = infer_chi_from_probe(
        probe_data=probe_data,
        confusion_matrix=CONFUSION_MATRIX,
        n_shots=N_SHOTS_PROBE,
        chi_initial=CHI_PRIOR,
        chi_higher_initial=0.0,
    )
    print(f"Inferred chi     = {infer_result['chi'] / (2*np.pi) / 1e6:.4f} MHz  (true: {CHI_TRUE/(2*np.pi)/1e6:.4f} MHz)")
    print(f"Inferred chi_h   = {infer_result['chi_higher'] / (2*np.pi) / 1e3:.2f} kHz  (true: {CHI_HIGHER_TRUE/(2*np.pi)/1e3:.2f} kHz)")
    print(f"Inferred T2*     = {infer_result['T2_star'] * 1e6:.2f} us  (true: {t2_star*1e6:.2f} us)")
    print(f"Chi error        = {abs(infer_result['chi'] - CHI_TRUE) / abs(CHI_TRUE) * 100:.2f}%")
    print(f"Fit converged    = {infer_result['fit_success']}")
    print(f"sigma_chi        = {infer_result['chi_uncertainty'] / (2*np.pi) / 1e3:.1f} kHz")
