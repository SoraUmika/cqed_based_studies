"""
Shared physical parameters, ODE integrator, objective functions, pulse envelope
builders, and GRAPE adjoint-gradient machinery for the readout pulse optimization study.

Physical model (linear dispersive, rotating frame at drive frequency ω_d):

    dα_j/dt = −(κ/2 + iΔ_j) α_j − i ε(t)      j ∈ {g, e}

where:
    α_j  : resonator coherent amplitude (qubit-state-conditioned)
    κ    : resonator decay rate (rad/s)
    Δ_g  = ω_r − ω_d        (ground state detuning)
    Δ_e  = ω_r + χ − ω_d    (excited state detuning)
    ε(t) : complex drive envelope (rad/s)

Output field (input-output theory):  s_j(t) = √κ · α_j(t)

Sign convention follows cqed_sim.measurement.ReadoutResonator:
    α_ss = −iε / (κ/2 + iΔ)    (steady-state solution)

cqed_sim gap documented in README: ReadoutResonator.response_trace() only supports
constant ε; the ODE integrators and GRAPE gradient below extend to arbitrary ε(t).
"""

from __future__ import annotations

import sys
import os
import math
import numpy as np

from runtime_compat import patch_windows_qutip_import

patch_windows_qutip_import()

# ── cqed_sim path (imports deferred; cqed_sim/qutip take minutes on first run) ─
_CQED_SIM_PATH = (
    r"C:\Users\dazzl\Box\Shyam Shankar Quantum Circuits Group"
    r"\Users\Users_JianJun\cQED_simulation"
)
if _CQED_SIM_PATH not in sys.path:
    sys.path.insert(0, _CQED_SIM_PATH)


def _get_readout_resonator():
    """Lazy import of ReadoutResonator to avoid slow qutip startup at module level."""
    from cqed_sim.measurement import ReadoutResonator
    return ReadoutResonator

TWO_PI = 2.0 * np.pi

# ── Nominal physical parameters ───────────────────────────────────────────────
#   Readout resonator
OMEGA_R     = TWO_PI * 7.5e9       # rad/s   resonator frequency (7.5 GHz)
KAPPA       = TWO_PI * 5.0e6       # rad/s   resonator linewidth (κ/2π = 5 MHz)
CHI_NOMINAL = TWO_PI * (-5.0e6)    # rad/s   dispersive shift (χ/2π = −5 MHz → χ/κ = 1)
G_COUP      = TWO_PI * 100.0e6     # rad/s   qubit–resonator coupling
OMEGA_Q     = TWO_PI * 5.5e9       # rad/s   qubit frequency

#   Drive amplitude
EPSILON_MAX = TWO_PI * 5.0e6       # rad/s   max drive amplitude (≤5 ss photons on res.)

#   Qubit decoherence
T1_QUBIT    = 30.0e-6              # s       qubit T1
T2_QUBIT    = 20.0e-6              # s       qubit T2
T_PHI_QUBIT = 1.0 / max(1.0/T2_QUBIT - 0.5/T1_QUBIT, 1e-30)  # s  pure dephasing time

#   Scan parameters
CHI_KAPPA_VALUES = np.array([0.25, 0.5, 1.0, 2.0, 4.0])     # χ/κ scan
CHI_VALUES = CHI_KAPPA_VALUES * KAPPA                         # |χ| values (rad/s)

#   Readout duration: dimensionless κT values to scan
KAPPA_T_VALUES = np.array([0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0])
T_READ_VALUES  = KAPPA_T_VALUES / KAPPA                        # seconds

#   GRAPE settings
N_GRAPE_SEGMENTS = 60      # piecewise-constant time segments
N_GRAPE_RESTARTS = 10      # random restarts for global coverage
GRAPE_MAXITER    = 300     # L-BFGS-B max iterations per restart

#   ODE integration
DT_ODE = 2.0e-9            # s   nominal time step (2 ns)


# ═══════════════════════════════════════════════════════════════════════════════
# 1. ODE INTEGRATORS (gap in cqed_sim: ReadoutResonator only supports const ε)
# ═══════════════════════════════════════════════════════════════════════════════

def integrate_readout_ode(
    epsilon_t: np.ndarray,
    tlist: np.ndarray,
    kappa: float,
    delta: float,
    alpha0: complex = 0.0,
) -> np.ndarray:
    """
    Integrate the linear dispersive resonator ODE with piecewise-constant drive ε(t):

        dα/dt = −(κ/2 + iΔ) α − i ε(t)

    Uses the exact piecewise-constant propagator (matrix exponential for scalar):
        α(t_{k+1}) = A · α(t_k) + B · ε_k
        A = exp(−λ τ_k),   B = −i(A − 1)/λ,   λ = κ/2 + iΔ

    This is exact (not an approximation) for each constant-ε segment, and
    is orders of magnitude faster than a Python-loop RK4 for long trajectories.

    Args:
        epsilon_t : drive envelope, one complex value per interval, length N.
        tlist     : time grid, length N+1 (may be non-uniform).
        kappa     : resonator decay rate (rad/s).
        delta     : resonator-drive detuning Δ = ω_r(state) − ω_d (rad/s).
        alpha0    : initial cavity field amplitude (default 0).

    Returns:
        alpha : complex field trajectory at tlist points, shape (N+1,).
    """
    lambda_eff = complex(0.5 * kappa + 1j * delta)
    N = len(tlist) - 1
    eps = np.asarray(epsilon_t, dtype=np.complex128)
    dt_arr = np.diff(tlist)           # shape (N,)

    # Precompute propagator coefficients for each segment
    A_arr = np.exp(-lambda_eff * dt_arr)         # shape (N,)
    if abs(lambda_eff) > 1e-30:
        B_arr = -1j * (A_arr - 1.0) / lambda_eff
    else:
        B_arr = -1j * dt_arr                     # limit λ→0

    # Fast path: constant-amplitude drive → use analytical solution
    # α(t) = α_ss + (α₀ − α_ss) exp(−λ t)   where α_ss = −iε/λ
    if eps.size > 0 and np.allclose(eps, eps[0], rtol=1e-10, atol=0):
        eps_c = complex(eps[0])
        if abs(lambda_eff) > 1e-30:
            alpha_ss = -1j * eps_c / lambda_eff
        else:
            alpha_ss = 0.0 + 0j
        t0 = float(tlist[0])
        return alpha_ss + (complex(alpha0) - alpha_ss) * np.exp(-lambda_eff * (tlist - t0))

    # General piecewise-constant path: step through segments
    alpha = np.empty(N + 1, dtype=np.complex128)
    alpha[0] = complex(alpha0)
    for k in range(N):
        alpha[k + 1] = A_arr[k] * alpha[k] + B_arr[k] * eps[k]

    return alpha


def simulate_conditioned_fields(
    epsilon_t: np.ndarray,
    tlist: np.ndarray,
    kappa: float,
    chi: float,
    delta_g: float = 0.0,
    alpha_g0: complex = 0.0,
    alpha_e0: complex = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Simulate qubit-state-conditioned resonator fields for the same drive ε(t).

    Args:
        epsilon_t : drive envelope, length len(tlist)−1.
        tlist     : time grid, length N+1.
        kappa     : resonator decay rate (rad/s).
        chi       : dispersive shift (rad/s).  Δ_e = Δ_g + χ.
        delta_g   : ground-state resonator-drive detuning (default 0 → drive on |g⟩).
        alpha_g0, alpha_e0 : initial fields.

    Returns:
        (alpha_g, alpha_e) : conditioned field trajectories at tlist points.
    """
    delta_e = delta_g + chi
    alpha_g = integrate_readout_ode(epsilon_t, tlist, kappa, delta_g, alpha_g0)
    alpha_e = integrate_readout_ode(epsilon_t, tlist, kappa, delta_e, alpha_e0)
    return alpha_g, alpha_e


# ═══════════════════════════════════════════════════════════════════════════════
# 2. OBJECTIVE FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def snr_squared(
    alpha_g: np.ndarray,
    alpha_e: np.ndarray,
    tlist: np.ndarray,
    kappa: float,
) -> float:
    """
    Heterodyne-detected SNR² for quantum-limited readout:

        SNR² = κ ∫₀ᵀ |α_e(t) − α_g(t)|² dt

    This equals |s_e − s_g|² integrated in time, where s_j = √κ α_j is the
    output field (input-output theory). Normalized to vacuum noise per unit time.

    Returns a non-negative float.
    """
    diff_sq = np.abs(alpha_e - alpha_g) ** 2
    return float(kappa * np.trapezoid(diff_sq, tlist))


def endpoint_separation_sq(alpha_g: np.ndarray, alpha_e: np.ndarray) -> float:
    """
    |α_e(T) − α_g(T)|²  — endpoint pointer-state separation.

    Captures only the final-time phase-space separation, ignoring the transient.
    """
    return float(abs(alpha_e[-1] - alpha_g[-1]) ** 2)


def assignment_fidelity_from_snr2(snr2: float) -> float:
    """
    Assignment fidelity for symmetric Gaussian readout noise:

        F_A = 1 − (1/2) erfc(√(SNR²/2)) = 1 − (1/2) erfc(SNR/√2)

    Valid for heterodyne detection in the quantum limit.
    """
    snr = math.sqrt(max(float(snr2), 0.0))
    return float(1.0 - 0.5 * math.erfc(snr / math.sqrt(2.0)))


def residual_photons(alpha_g: np.ndarray, alpha_e: np.ndarray) -> float:
    """
    Mean residual resonator photon number at end of pulse:
        n_res = (|α_e(T)|² + |α_g(T)|²) / 2
    """
    return float(0.5 * (abs(alpha_e[-1]) ** 2 + abs(alpha_g[-1]) ** 2))


def measurement_dephasing_factor(
    alpha_g: np.ndarray,
    alpha_e: np.ndarray,
    tlist: np.ndarray,
    kappa: float,
) -> float:
    """
    Measurement-induced dephasing factor for the qubit during readout.

    The qubit off-diagonal density matrix element decays as:
        ρ_{ge}(T) = ρ_{ge}(0) · exp(−∫₀ᵀ (κ/2)|α_e − α_g|² dt)

    Returns the exponent factor η_φ = exp(−Γ_φ T) ∈ [0, 1].

    Note: Γ_φ = (κ/2)|α_e − α_g|² is the instantaneous dephasing rate.
    This is distinct from the measurement rate γ_meas = (κ/2)|Δα_ss|² (steady-state).
    """
    integrand = (kappa / 2.0) * np.abs(alpha_e - alpha_g) ** 2
    total_dephasing = np.trapezoid(integrand, tlist)
    return float(np.exp(-total_dephasing))


def multi_objective(
    alpha_g: np.ndarray,
    alpha_e: np.ndarray,
    tlist: np.ndarray,
    kappa: float,
    w_snr: float = 1.0,
    w_res: float = 0.1,
    n_res_scale: float = 1.0,
) -> float:
    """
    Multi-objective readout cost (to be minimized):

        L = −w_snr · SNR² + w_res · n_res(T) / n_res_scale

    n_res_scale normalizes the residual-photon penalty relative to SNR².
    """
    s2 = snr_squared(alpha_g, alpha_e, tlist, kappa)
    n_res = residual_photons(alpha_g, alpha_e)
    return -w_snr * s2 + w_res * n_res / max(n_res_scale, 1e-30)


# ═══════════════════════════════════════════════════════════════════════════════
# 3. OPTIMAL DRIVE FREQUENCY AND STEADY-STATE ANALYSIS
# ═══════════════════════════════════════════════════════════════════════════════

def optimal_delta_g(chi: float) -> float:
    """
    Ground-state detuning Δ_g that maximizes steady-state |α_e^ss − α_g^ss|.

    At Δ_g = −χ/2 the drive sits at the midpoint between the two resonances,
    maximizing the pointer separation per unit drive amplitude.
    """
    return -0.5 * chi


def steady_state_alpha(epsilon: complex, kappa: float, delta: float) -> complex:
    """Steady-state amplitude: α_ss = −iε / (κ/2 + iΔ)."""
    return complex(-1j * epsilon / (0.5 * kappa + 1j * delta))


def steady_state_separation(epsilon: complex, kappa: float, chi: float, delta_g: float) -> float:
    """
    |α_e^ss − α_g^ss| for given drive detuning Δ_g.
    Δ_e = Δ_g + χ.
    """
    delta_e = delta_g + chi
    ag = steady_state_alpha(epsilon, kappa, delta_g)
    ae = steady_state_alpha(epsilon, kappa, delta_e)
    return float(abs(ae - ag))


def max_steady_state_separation(epsilon: complex, kappa: float, chi: float) -> float:
    """
    Maximum steady-state separation at optimal drive frequency (Δ_g = −χ/2).

    Analytical result:
        |Δα^ss_max| = 4|ε||χ| / (κ² + χ²) = 4|ε|/κ · |χ/κ| / (1 + (χ/κ)²)

    Maximized at |χ/κ| = 1 → |Δα^ss_max|_{χ=κ} = 2|ε|/κ.
    """
    delta_g = optimal_delta_g(chi)
    return steady_state_separation(epsilon, kappa, chi, delta_g)


def optimal_separation_formula(epsilon: float, kappa: float, chi: float) -> float:
    """
    Analytical formula for max steady-state separation (real form):
        |Δα^ss_max| = 4|ε||χ| / (κ² + χ²)
    """
    return 4.0 * abs(epsilon) * abs(chi) / (kappa**2 + chi**2)


# ═══════════════════════════════════════════════════════════════════════════════
# 4. PULSE ENVELOPE BUILDERS (midpoint samples for ODE integration)
# ═══════════════════════════════════════════════════════════════════════════════

def _midpoints(tlist: np.ndarray) -> np.ndarray:
    """Return midpoints of tlist intervals."""
    return 0.5 * (tlist[:-1] + tlist[1:])


def build_square(tlist: np.ndarray, amplitude: complex) -> np.ndarray:
    """Rectangular (constant) drive envelope."""
    return amplitude * np.ones(len(tlist) - 1, dtype=np.complex128)


def build_gaussian(
    tlist: np.ndarray, amplitude: complex, sigma_frac: float = 0.25
) -> np.ndarray:
    """
    Gaussian envelope centered at t = T/2 with σ = sigma_frac · T.

    Peak-normalized so the instantaneous drive never exceeds |amplitude|.
    """
    T = float(tlist[-1] - tlist[0])
    t0 = float(tlist[0])
    t_mid = _midpoints(tlist)
    sigma = sigma_frac * T
    env = np.exp(-0.5 * ((t_mid - (t0 + 0.5 * T)) / sigma) ** 2)
    peak = np.max(np.abs(env))
    if peak > 1e-15:
        env = env / peak
    phase_factor = amplitude / (abs(amplitude) + 1e-60)
    return (phase_factor * abs(amplitude) * env).astype(np.complex128)


def build_hann(tlist: np.ndarray, amplitude: complex) -> np.ndarray:
    """
    Hann (cosine-squared / von Hann) window: env(t) = sin²(π t/T).

    Peak-normalized so the instantaneous drive never exceeds |amplitude|.
    """
    T = float(tlist[-1] - tlist[0])
    t0 = float(tlist[0])
    t_mid = _midpoints(tlist)
    t_rel = (t_mid - t0) / T
    env = np.sin(np.pi * t_rel) ** 2   # area = T/2
    phase_factor = amplitude / (abs(amplitude) + 1e-60)
    return (phase_factor * abs(amplitude) * env).astype(np.complex128)


def build_cosine_rise(
    tlist: np.ndarray, amplitude: complex, rise_frac: float = 0.15
) -> np.ndarray:
    """
    Flat-top with cosine rise/fall edges (trapezoidal cosine).

        env(t) = 0.5 (1 − cos(πt/τ_rise))   for t < τ_rise
                 1                             for τ_rise ≤ t ≤ T − τ_rise
                 0.5 (1 − cos(π(T−t)/τ_rise)) for t > T − τ_rise

    Peak-normalized so the instantaneous drive never exceeds |amplitude|.
    """
    T = float(tlist[-1] - tlist[0])
    t0 = float(tlist[0])
    t_rise = rise_frac * T
    t_mid = _midpoints(tlist)
    t_rel = t_mid - t0
    env = np.ones(len(t_rel))
    left = t_rel < t_rise
    right = t_rel > (T - t_rise)
    if t_rise > 0:
        env[left] = 0.5 * (1.0 - np.cos(np.pi * t_rel[left] / t_rise))
        env[right] = 0.5 * (1.0 - np.cos(np.pi * (T - t_rel[right]) / t_rise))
    peak = np.max(np.abs(env))
    if peak > 1e-15:
        env = env / peak
    phase_factor = amplitude / (abs(amplitude) + 1e-60)
    return (phase_factor * abs(amplitude) * env).astype(np.complex128)


def build_spline(
    tlist: np.ndarray, amplitude: complex, n_knots: int = 4
) -> np.ndarray:
    """
    Cubic spline envelope with n_knots interior knots, clamped to zero at edges.

    Peak-normalized so the instantaneous drive never exceeds |amplitude|.
    """
    from scipy.interpolate import CubicSpline  # lazy import — scipy.special hangs if top-level
    T = float(tlist[-1] - tlist[0])
    t0 = float(tlist[0])
    # Knot positions in [0, 1] including boundary (clamped to 0 at both ends)
    x_knots = np.linspace(0.0, 1.0, n_knots + 2)
    # Default bell-shaped values at knots
    y_knots = np.sin(np.pi * x_knots) ** 2
    cs = CubicSpline(x_knots, y_knots, bc_type="clamped")
    t_mid = _midpoints(tlist)
    t_rel = (t_mid - t0) / T
    env = np.clip(np.real(cs(t_rel)), 0.0, None)
    peak = np.max(np.abs(env))
    if peak > 1e-15:
        env = env / peak
    phase_factor = amplitude / (abs(amplitude) + 1e-60)
    return (phase_factor * abs(amplitude) * env).astype(np.complex128)


def build_depletion_pulse(
    tlist: np.ndarray,
    alpha_at_t0: complex,
    kappa: float,
    delta_g: float,
    delta_e: float,
    n_steps: int = 1,
) -> np.ndarray:
    """
    Compute the optimal piecewise-constant depletion pulse that drives α → 0.

    For the linear model with constant depletion drive ε_dep, the resonator field
    evolves as:
        α(t) = α₀ exp(−λ t) − iε_dep (exp(−λ t) − 1) / λ

    Setting α(T_dep) = 0 gives ε_dep = iα₀ λ / (1 − exp(−λ T_dep)).

    Here we return a simple single-segment depletion drive for α_g(t0),
    which is a reasonable approximation for the ground state. The returned
    envelope drives α_g → 0; α_e will not be perfectly depleted (residual
    depends on the χ/κ ratio).

    Args:
        tlist     : time grid for the depletion segment, length N_dep+1.
        alpha_at_t0 : resonator field amplitude at start of depletion segment.
        kappa, delta_g, delta_e : resonator parameters.
        n_steps   : number of piecewise segments (1 = optimal single depletion).

    Returns:
        epsilon_dep : complex drive envelope, length len(tlist)−1.
    """
    T_dep = float(tlist[-1] - tlist[0])
    # Use mean of g and e detunings for depletion (approximation)
    delta_mean = 0.5 * (delta_g + delta_e)
    lambda_dep = 0.5 * kappa + 1j * delta_mean
    denom = 1.0 - np.exp(-lambda_dep * T_dep)
    if abs(denom) < 1e-14:
        epsilon_dep = 0.0 + 0j
    else:
        epsilon_dep = 1j * complex(alpha_at_t0) * lambda_dep / denom
    return epsilon_dep * np.ones(len(tlist) - 1, dtype=np.complex128)


PULSE_FAMILIES = {
    "square":      build_square,
    "gaussian":    build_gaussian,
    "hann":        build_hann,
    "cosine_rise": build_cosine_rise,
    "spline":      build_spline,
}


# ═══════════════════════════════════════════════════════════════════════════════
# 5. GRAPE ADJOINT GRADIENT (linear dispersive model)
# ═══════════════════════════════════════════════════════════════════════════════
#
# Piecewise-constant parameterization: ε(t) = ε_k for t ∈ [t_k, t_{k+1}),
# time step τ = T/N (uniform grid assumed).
#
# Propagator per segment:
#   α_j(t_{k+1}) = A_j · α_j(t_k) + B_j · ε_k
#   A_j = exp(−λ_j τ),   λ_j = κ/2 + iΔ_j
#   B_j = −i (A_j − 1) / λ_j
#
# SNR² objective (to maximize):
#   J = κτ Σ_{n=1}^{N} |Δα_n|²    (Δα_n = α_e(t_n) − α_g(t_n))
#
# Wirtinger gradient ∂J/∂ε_k*:
#   ∂J/∂ε_k* = κτ (B_e · Λ_e^k − B_g · Λ_g^k)
#
# where Λ_j^k = Σ_{n=k}^{N-1} A_j^{n-k} (Δα_{n+1})*   [backward recursion]
#
# Real-variable gradient (ε_k = x_k + i y_k):
#   ∂J/∂x_k = 2 Re[∂J/∂ε_k*]
#   ∂J/∂y_k = 2 Im[∂J/∂ε_k*]
#
# This is equivalent to the adjoint / co-state method.

def grape_propagators(
    kappa: float, delta_g: float, delta_e: float, tau: float
) -> tuple[complex, complex, complex, complex]:
    """
    Compute per-step propagator coefficients A_j, B_j for GRAPE.

    Returns:
        A_g, B_g, A_e, B_e
    """
    lambda_g = 0.5 * kappa + 1j * delta_g
    lambda_e = 0.5 * kappa + 1j * delta_e
    A_g = np.exp(-lambda_g * tau)
    A_e = np.exp(-lambda_e * tau)
    B_g = -1j * (A_g - 1.0) / lambda_g if abs(lambda_g) > 1e-30 else -1j * tau
    B_e = -1j * (A_e - 1.0) / lambda_e if abs(lambda_e) > 1e-30 else -1j * tau
    return complex(A_g), complex(B_g), complex(A_e), complex(B_e)


def grape_forward_pass(
    epsilon_vec: np.ndarray,
    kappa: float,
    delta_g: float,
    delta_e: float,
    tau: float,
    alpha_g0: complex = 0.0,
    alpha_e0: complex = 0.0,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Forward propagation for GRAPE (piecewise-constant ε).

    Args:
        epsilon_vec : complex drive array, shape (N,).
        kappa, delta_g, delta_e, tau : resonator parameters and step size.
        alpha_g0, alpha_e0 : initial fields (default 0).

    Returns:
        (alpha_g, alpha_e) : field trajectories at time points 0..N, shape (N+1,).
    """
    N = len(epsilon_vec)
    A_g, B_g, A_e, B_e = grape_propagators(kappa, delta_g, delta_e, tau)
    alpha_g = np.zeros(N + 1, dtype=np.complex128)
    alpha_e = np.zeros(N + 1, dtype=np.complex128)
    alpha_g[0] = complex(alpha_g0)
    alpha_e[0] = complex(alpha_e0)
    for k in range(N):
        eps = complex(epsilon_vec[k])
        alpha_g[k + 1] = A_g * alpha_g[k] + B_g * eps
        alpha_e[k + 1] = A_e * alpha_e[k] + B_e * eps
    return alpha_g, alpha_e


def grape_backward_pass(
    alpha_g: np.ndarray,
    alpha_e: np.ndarray,
    kappa: float,
    tau: float,
    A_g: complex,
    A_e: complex,
    B_g: complex,
    B_e: complex,
) -> np.ndarray:
    """
    Compute Wirtinger gradient ∂J/∂ε_k* using the adjoint (backward) pass.

    Returns:
        grad_wirtinger : complex array of shape (N,).
    """
    N = len(alpha_g) - 1
    delta_alpha = alpha_e - alpha_g    # shape (N+1,)

    # Backward recursion: Λ_j^k = Σ_{n=k}^{N-1} A_j^{n-k} (Δα_{n+1})*
    Lambda_g = np.zeros(N, dtype=np.complex128)
    Lambda_e = np.zeros(N, dtype=np.complex128)

    # Seed at k = N−1
    Lambda_g[N - 1] = np.conj(delta_alpha[N])
    Lambda_e[N - 1] = np.conj(delta_alpha[N])

    for k in range(N - 2, -1, -1):
        Lambda_g[k] = np.conj(delta_alpha[k + 1]) + A_g * Lambda_g[k + 1]
        Lambda_e[k] = np.conj(delta_alpha[k + 1]) + A_e * Lambda_e[k + 1]

    # ∂J/∂ε_k* = κτ (B_e Λ_e^k − B_g Λ_g^k)
    grad_wirtinger = kappa * tau * (B_e * Lambda_e - B_g * Lambda_g)
    return grad_wirtinger


def grape_objective_and_grad(
    x: np.ndarray,
    kappa: float,
    delta_g: float,
    delta_e: float,
    tau: float,
    N: int,
    w_snr: float = 1.0,
    w_res: float = 0.0,
    n_res_scale: float = 1.0,
) -> tuple[float, np.ndarray]:
    """
    Compute the GRAPE objective (negative SNR²) and its real gradient w.r.t. x.

    Args:
        x  : real parameter vector [Re(ε_0),...,Re(ε_{N-1}), Im(ε_0),...,Im(ε_{N-1})].
        w_snr : weight for SNR² term (positive → maximize SNR²).
        w_res : weight for residual photon penalty (positive → penalize n_res).

    Returns:
        (objective, gradient) where objective is to be minimized.
    """
    epsilon_vec = x[:N] + 1j * x[N:]
    A_g, B_g, A_e, B_e = grape_propagators(kappa, delta_g, delta_e, tau)

    # Forward pass
    alpha_g, alpha_e = grape_forward_pass(epsilon_vec, kappa, delta_g, delta_e, tau)
    delta_alpha = alpha_e - alpha_g
    tlist = np.arange(N + 1) * tau

    # SNR² = κτ Σ |Δα_n|²
    snr2 = float(kappa * tau * np.sum(np.abs(delta_alpha) ** 2))

    # Residual photons
    n_res = float(0.5 * (abs(alpha_e[-1]) ** 2 + abs(alpha_g[-1]) ** 2))

    # Objective: minimize L = −w_snr·SNR² + w_res·n_res/n_res_scale
    objective = -w_snr * snr2 + w_res * n_res / max(n_res_scale, 1e-30)

    # ── Gradient w.r.t. ε_k (Wirtinger, for SNR² part) ──────────────────────
    # ∂(−SNR²)/∂ε_k* = −κτ (B_e Λ_e^k − B_g Λ_g^k)
    grad_w = grape_backward_pass(alpha_g, alpha_e, kappa, tau, A_g, A_e, B_g, B_e)
    grad_w_snr = -w_snr * grad_w   # negate for minimization

    # ── Gradient for residual photon penalty ─────────────────────────────────
    # n_res = (|α_e(T)|² + |α_g(T)|²)/2
    # ∂n_res/∂ε_k* = Σ_j (A_j^{N-1-k} B_j α_j(T)*) / 2
    if w_res != 0.0:
        grad_w_res = np.zeros(N, dtype=np.complex128)
        for k in range(N):
            n_steps_remaining = N - 1 - k
            grad_w_res[k] = 0.5 * (
                A_e ** n_steps_remaining * B_e * np.conj(alpha_e[-1])
                + A_g ** n_steps_remaining * B_g * np.conj(alpha_g[-1])
            )
        grad_w_res = w_res / max(n_res_scale, 1e-30) * grad_w_res
    else:
        grad_w_res = np.zeros(N, dtype=np.complex128)

    grad_wirtinger = grad_w_snr + grad_w_res

    # Convert Wirtinger → real gradient:  ∂L/∂x_k = 2 Re[∂L/∂ε_k*]
    #                                      ∂L/∂y_k = -2 Im[∂L/∂ε_k*]
    grad_real = np.concatenate([
        2.0 * np.real(grad_wirtinger),
        -2.0 * np.imag(grad_wirtinger),
    ])
    return float(objective), grad_real


def run_grape(
    kappa: float,
    chi: float,
    T_read: float,
    N_seg: int,
    delta_g: float,
    epsilon_max: float,
    n_restarts: int = N_GRAPE_RESTARTS,
    maxiter: int = GRAPE_MAXITER,
    w_snr: float = 1.0,
    w_res: float = 0.0,
    n_res_scale: float = 1.0,
    seed: int = 42,
) -> dict:
    """
    Run GRAPE optimization for the linear dispersive readout model.

    Args:
        kappa       : resonator decay rate (rad/s).
        chi         : dispersive shift (rad/s).
        T_read      : readout duration (s).
        N_seg       : number of piecewise-constant segments.
        delta_g     : ground-state detuning (rad/s).
        epsilon_max : amplitude bound (rad/s).
        n_restarts  : number of random restarts.
        maxiter     : max L-BFGS-B iterations per restart.
        w_snr, w_res, n_res_scale : objective weights.
        seed        : NumPy RNG seed.

    Returns:
        dict with keys:
            'epsilon_opt'  : complex array of optimized pulse, shape (N_seg,)
            'snr2_opt'     : optimized SNR²
            'n_res_opt'    : residual photons
            'F_assign_opt' : assignment fidelity
            'alpha_g', 'alpha_e' : field trajectories (N_seg+1,)
            'objective_history': list of best objectives across restarts
    """
    tau = T_read / N_seg
    delta_e = delta_g + chi

    # Try scipy L-BFGS-B first; fall back to Adam (numpy-only) if scipy hangs/unavailable.
    try:
        import importlib
        _scipy_opt = importlib.import_module("scipy.optimize")
        _minimize = _scipy_opt.minimize
        _use_scipy = True
    except Exception:
        _use_scipy = False

    rng = np.random.default_rng(seed)
    best_obj = np.inf
    best_x = None
    obj_history = []

    initial_guesses = [
        np.concatenate([
            np.full(N_seg, float(epsilon_max), dtype=np.float64),
            np.zeros(N_seg, dtype=np.float64),
        ])
    ]
    initial_guesses.extend(
        rng.uniform(-epsilon_max, epsilon_max, size=2 * N_seg)
        for _ in range(n_restarts)
    )

    for x0 in initial_guesses:
        if _use_scipy:
            bounds = [(-epsilon_max, epsilon_max)] * (2 * N_seg)
            result = _minimize(
                grape_objective_and_grad,
                x0,
                jac=True,
                method="L-BFGS-B",
                bounds=bounds,
                args=(kappa, delta_g, delta_e, tau, N_seg, w_snr, w_res, n_res_scale),
                options={"maxiter": maxiter, "ftol": 1e-12, "gtol": 1e-8},
            )
            x_final = result.x
            obj_final = result.fun
        else:
            x_final, obj_final = _adam_minimize(
                grape_objective_and_grad,
                x0,
                lo=-epsilon_max,
                hi=epsilon_max,
                maxiter=max(maxiter * 3, 500),
                args=(kappa, delta_g, delta_e, tau, N_seg, w_snr, w_res, n_res_scale),
            )
        obj_history.append(obj_final)
        if obj_final < best_obj:
            best_obj = obj_final
            best_x = x_final

    eps_opt = best_x[:N_seg] + 1j * best_x[N_seg:]
    tlist = np.linspace(0.0, T_read, N_seg + 1)
    ag, ae = grape_forward_pass(eps_opt, kappa, delta_g, delta_e, tau)
    snr2_val = snr_squared(ag, ae, tlist, kappa)
    n_res_val = residual_photons(ag, ae)
    F_val = assignment_fidelity_from_snr2(snr2_val)

    return {
        "epsilon_opt": eps_opt,
        "snr2_opt": snr2_val,
        "n_res_opt": n_res_val,
        "F_assign_opt": F_val,
        "alpha_g": ag,
        "alpha_e": ae,
        "objective_history": obj_history,
        "tlist": tlist,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# 6. PULSE-FAMILY EVALUATION (for Phases 2 and 5)
# ═══════════════════════════════════════════════════════════════════════════════

def evaluate_pulse_family(
    family: str,
    amplitude: complex,
    T_read: float,
    kappa: float,
    chi: float,
    delta_g: float,
    dt: float = DT_ODE,
) -> dict:
    """
    Evaluate a named pulse family with given amplitude and duration.

    Returns a dict with SNR², F_assign, endpoint separation, residual photons,
    and measurement dephasing factor.
    """
    N = max(2, int(round(T_read / dt)))
    tlist = np.linspace(0.0, T_read, N + 1)
    builder = PULSE_FAMILIES[family]
    epsilon_t = builder(tlist, amplitude)
    ag, ae = simulate_conditioned_fields(epsilon_t, tlist, kappa, chi, delta_g)
    snr2_val = snr_squared(ag, ae, tlist, kappa)
    return {
        "snr2": snr2_val,
        "F_assign": assignment_fidelity_from_snr2(snr2_val),
        "endpoint_sep2": endpoint_separation_sq(ag, ae),
        "n_res": residual_photons(ag, ae),
        "dephasing_factor": measurement_dephasing_factor(ag, ae, tlist, kappa),
        "alpha_g": ag,
        "alpha_e": ae,
        "tlist": tlist,
        "epsilon_t": epsilon_t,
    }


def optimize_amplitude(
    family: str,
    T_read: float,
    kappa: float,
    chi: float,
    delta_g: float,
    epsilon_max: float = EPSILON_MAX,
    dt: float = DT_ODE,
    n_grid: int = 80,
) -> tuple[float, dict]:
    """
    1-D golden-section search over |amplitude| ∈ [0, epsilon_max] to maximize SNR².

    Returns:
        (epsilon_opt_abs, eval_dict)
    """
    eps_grid = np.linspace(epsilon_max * 0.01, epsilon_max, n_grid)
    best_snr2 = -np.inf
    best_eps = eps_grid[0]
    for eps in eps_grid:
        ev = evaluate_pulse_family(family, float(eps), T_read, kappa, chi, delta_g, dt)
        if ev["snr2"] > best_snr2:
            best_snr2 = ev["snr2"]
            best_eps = float(eps)

    # Fine refinement around best
    lo = max(epsilon_max * 0.01, best_eps * 0.8)
    hi = min(epsilon_max, best_eps * 1.2)
    eps_fine = np.linspace(lo, hi, 40)
    for eps in eps_fine:
        ev = evaluate_pulse_family(family, float(eps), T_read, kappa, chi, delta_g, dt)
        if ev["snr2"] > best_snr2:
            best_snr2 = ev["snr2"]
            best_eps = float(eps)

    best_eval = evaluate_pulse_family(family, best_eps, T_read, kappa, chi, delta_g, dt)
    return best_eps, best_eval


# ═══════════════════════════════════════════════════════════════════════════════
# 7. CQED_SIM BRIDGE (steady-state validation)
# ═══════════════════════════════════════════════════════════════════════════════

def make_readout_resonator(
    omega_r: float = OMEGA_R,
    kappa: float = KAPPA,
    chi: float = CHI_NOMINAL,
    epsilon: complex = EPSILON_MAX,
    drive_frequency: float | None = None,
):
    """
    Construct a cqed_sim ReadoutResonator for steady-state validation.
    Drive frequency defaults to the midpoint ω_d = ω_r + χ/2.
    (cqed_sim imported lazily to avoid slow qutip startup.)
    """
    ReadoutResonator = _get_readout_resonator()
    if drive_frequency is None:
        drive_frequency = omega_r + 0.5 * chi
    return ReadoutResonator(
        omega_r=omega_r,
        kappa=kappa,
        g=G_COUP,
        epsilon=epsilon,
        chi=chi,
        drive_frequency=drive_frequency,
    )


def cqed_sim_steady_state(
    kappa: float, chi: float, epsilon: complex, delta_g: float
) -> tuple[complex, complex]:
    """
    Compute steady-state amplitudes via cqed_sim.ReadoutResonator for cross-check.

    Uses the analytical formula α_ss = −iε/(κ/2 + iΔ) internally.
    (cqed_sim imported lazily to avoid slow qutip startup.)
    """
    ReadoutResonator = _get_readout_resonator()
    omega_r = OMEGA_R
    omega_d = omega_r - delta_g      # delta_g = omega_r - omega_d
    res = ReadoutResonator(
        omega_r=omega_r, kappa=kappa, g=G_COUP, epsilon=epsilon, chi=chi
    )
    amps = res.steady_state_amplitudes(drive_frequency=omega_d)
    return amps["g"], amps["e"]
