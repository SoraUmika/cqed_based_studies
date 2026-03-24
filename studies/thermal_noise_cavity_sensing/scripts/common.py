"""
Common physics definitions for the thermal-noise cavity sensing study.

Physical setup:
  - A post cavity (bosonic mode a, frequency omega_c) is coupled to one or more
    thermal baths, including a target device bath, a background cold bath, and
    an internal loss bath.
  - Optionally, a dispersively coupled ancilla qubit (transmon) allows ancilla-
    based readout of cavity photon number / distribution.

Documented cqed_sim gaps (see README):
  - Multi-bath Lindblad operators are built here manually (cqed_sim NoiseSpec
    supports only a single kappa/nth per cavity mode).
  - Steady-state solving uses qt.steadystate() directly.
  - Cavity-only model does not use DispersiveTransmonCavityModel.

Reference parameters match the sqr_pulse_waveform_design study for consistency.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

STUDY_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = STUDY_DIR / "data"
FIG_DIR = STUDY_DIR / "figures"
STYLE_PATH = (
    Path(__file__).resolve().parents[2]
    / ".github"
    / "skills"
    / "publication-figures"
    / "assets"
    / "cqed_style.mplstyle"
)

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------

HBAR = 1.0545718e-34  # J·s
KB = 1.380649e-23     # J/K

# ---------------------------------------------------------------------------
# Reference cqed_sim experimental parameters
# (matching sqr_pulse_waveform_design/scripts/common.py)
# ---------------------------------------------------------------------------

OMEGA_C = 2 * np.pi * 5.241e9   # rad/s  cavity frequency
OMEGA_Q = 2 * np.pi * 6.197e9   # rad/s  transmon frequency

# Dispersive coupling |chi| (magnitude; actual chi < 0 in cqed_sim convention).
# Qubit transition frequency in Fock state n: omega_q^(n) = omega_q - n*CHI_DISP
CHI_DISP = 2 * np.pi * 2.84e6   # rad/s

T1_Q = 20e-6                     # s   qubit T1
T2Q = 20e-6                      # s   qubit T2* (total, includes pure dephasing)

# Derived qubit rates
GAMMA_DOWN = 1.0 / T1_Q          # rad/s   relaxation rate
GAMMA_PHI = 1.0 / T2Q - 1.0 / (2 * T1_Q)   # rad/s   pure dephasing rate

# ---------------------------------------------------------------------------
# Default sensing study parameters
# ---------------------------------------------------------------------------

KAPPA_TOT = 2 * np.pi * 100e3   # rad/s  total cavity decay (100 kHz)
KAPPA_FRAC_TARGET = 0.50         # kappa_target / kappa_tot
KAPPA_FRAC_BG = 0.30             # kappa_background / kappa_tot
KAPPA_FRAC_INT = 0.20            # kappa_internal / kappa_tot

N_TH_BG = 0.01                   # background cold-bath occupation
N_TH_INT = 0.00                  # internal loss bath (zero temperature)

N_CAV = 30                       # default Fock-space truncation


# ---------------------------------------------------------------------------
# ThermalBath dataclass
# ---------------------------------------------------------------------------

@dataclass
class ThermalBath:
    """
    Specification for a single thermal bath coupled to a bosonic cavity mode.

    Attributes
    ----------
    kappa : float
        Coupling rate to this bath (rad/s). Must be >= 0.
    n_th : float
        Thermal occupation number of this bath (dimensionless). Must be >= 0.
    label : str
        Descriptive label for this bath (e.g. 'target', 'background').
    """

    kappa: float
    n_th: float
    label: str = ""

    def __post_init__(self) -> None:
        if self.kappa < 0:
            raise ValueError(f"kappa must be >= 0, got {self.kappa}")
        if self.n_th < 0:
            raise ValueError(f"n_th must be >= 0, got {self.n_th}")


# ---------------------------------------------------------------------------
# Thermal occupation conversions
# ---------------------------------------------------------------------------

def n_thermal(omega: float, T: float) -> float:
    """
    Bose-Einstein thermal occupation number at angular frequency omega
    and temperature T.

        n_th(omega, T) = 1 / (exp(hbar*omega / kB*T) - 1)

    Parameters
    ----------
    omega : float
        Angular frequency (rad/s).
    T : float
        Temperature (K). Returns 0 for T <= 0.
    """
    if T <= 0.0:
        return 0.0
    x = HBAR * omega / (KB * T)
    if x > 700.0:          # avoid overflow; n_th ~ exp(-x) -> 0
        return np.exp(-x)
    return 1.0 / np.expm1(x)


def temperature_from_nbar(omega: float, n_bar: float) -> float:
    """
    Effective temperature [K] corresponding to mean occupation n_bar
    at angular frequency omega.

        T = hbar*omega / (kB * ln(1 + 1/n_bar))

    Returns 0 for n_bar <= 0.
    """
    if n_bar <= 0.0:
        return 0.0
    return HBAR * omega / (KB * np.log1p(1.0 / n_bar))


# ---------------------------------------------------------------------------
# Analytic formulas for cavity thermal dynamics
# ---------------------------------------------------------------------------

def analytic_nbar_ss(baths: Sequence[ThermalBath]) -> float:
    """
    Analytic steady-state mean cavity photon number for multiple thermal baths.

        n̄_ss = (Σ_j κ_j n_j) / (Σ_j κ_j)

    This result follows from detailed balance applied to the multi-bath Lindblad
    equation.  The steady state is a thermal state with effective occupation n̄_ss,
    independent of the initial state.
    """
    kappa_tot = sum(b.kappa for b in baths)
    if kappa_tot <= 0.0:
        return 0.0
    return sum(b.kappa * b.n_th for b in baths) / kappa_tot


def analytic_kappa_tot(baths: Sequence[ThermalBath]) -> float:
    """Total cavity decay rate κ_tot = Σ_j κ_j (rad/s)."""
    return sum(b.kappa for b in baths)


def analytic_nbar_transient(
    t: np.ndarray, n_init: float, baths: Sequence[ThermalBath]
) -> np.ndarray:
    """
    Analytic transient mean cavity photon number starting from n_init.

        n̄(t) = n̄_ss + (n̄(0) − n̄_ss) * exp(−κ_tot * t)

    Parameters
    ----------
    t : np.ndarray
        Time array (s).
    n_init : float
        Initial mean photon number.
    baths : list of ThermalBath
        List of thermal baths.
    """
    n_ss = analytic_nbar_ss(baths)
    kappa_tot = analytic_kappa_tot(baths)
    return n_ss + (n_init - n_ss) * np.exp(-kappa_tot * np.asarray(t))


def thermal_pn(n_bar: float, n_max: int) -> np.ndarray:
    """
    Thermal photon-number distribution.

        P_n = n̄^n / (1 + n̄)^{n+1}

    Uses log-space arithmetic for numerical stability.

    Parameters
    ----------
    n_bar : float
        Mean photon number.
    n_max : int
        Number of Fock levels to compute (returns P_0 ... P_{n_max-1}).
    """
    n = np.arange(n_max, dtype=float)
    if n_bar <= 0.0:
        pn = np.zeros(n_max)
        pn[0] = 1.0
        return pn
    log_pn = n * np.log(n_bar) - (n + 1.0) * np.log1p(n_bar)
    return np.exp(log_pn)


def ramsey_coherence_thermal(
    tau: np.ndarray, n_bar: float, chi: float
) -> np.ndarray:
    """
    Qubit coherence function during Ramsey under a thermal cavity state
    (no qubit T2 envelope — apply separately).

    Exact analytic result for a thermal photon-number distribution:

        χ(τ) = Σ_n P_n e^{i n χ τ} = 1 / (1 + n̄(1 − e^{i χ τ}))

    Parameters
    ----------
    tau : np.ndarray
        Free-evolution time array (s).
    n_bar : float
        Mean cavity photon number.
    chi : float
        Dispersive coupling (rad/s).  Sign determines phase direction.
    """
    phase = chi * np.asarray(tau, dtype=float)
    return 1.0 / (1.0 + n_bar * (1.0 - np.exp(1j * phase)))


# ---------------------------------------------------------------------------
# Lindblad collapse-operator builders
# ---------------------------------------------------------------------------

def build_cavity_c_ops(a, baths: Sequence[ThermalBath]) -> list:
    """
    Build Lindblad collapse operators for a cavity coupled to multiple
    thermal baths.

    For each bath j:
        c_j^- = sqrt(κ_j * (n_j + 1)) * a      [emission into bath]
        c_j^+ = sqrt(κ_j * n_j)       * a†     [absorption from bath]

    The multi-bath Lindblad equation with these operators has the same
    steady state as a single effective bath with:
        κ_eff = Σ_j κ_j
        n_eff = Σ_j κ_j n_j / Σ_j κ_j
    but keeping them separate correctly handles the full dynamics.

    Parameters
    ----------
    a : qt.Qobj
        Cavity lowering operator.
    baths : list of ThermalBath
        List of thermal baths.
    """
    c_ops = []
    for bath in baths:
        if bath.kappa <= 0.0:
            continue
        # Emission (always present)
        c_ops.append(np.sqrt(bath.kappa * (bath.n_th + 1.0)) * a)
        # Absorption (only when bath is thermal)
        if bath.n_th > 0.0:
            c_ops.append(np.sqrt(bath.kappa * bath.n_th) * a.dag())
    return c_ops


def build_qubit_c_ops(sm, sz, gamma_down: float, gamma_phi: float,
                       n_th_q: float = 0.0) -> list:
    """
    Build Lindblad collapse operators for a two-level ancilla qubit.

    Operators:
        sqrt(Γ↓ (1 + n_q))  σ_-          [relaxation]
        sqrt(Γ↓ n_q)         σ_+          [thermal excitation]
        sqrt(Γ_φ / 2)        σ_z          [pure dephasing]

    The factor 1/2 in the dephasing operator follows the cqed_sim convention
    where σ_z has eigenvalues ±1 (factor from the Lindblad superoperator).

    Parameters
    ----------
    sm : qt.Qobj  lowering operator σ_-
    sz : qt.Qobj  Pauli-Z operator σ_z
    gamma_down : float  relaxation rate 1/T1 (rad/s)
    gamma_phi : float   pure dephasing rate (rad/s)
    n_th_q : float      qubit thermal occupation (negligible at 6 GHz / 20 mK)
    """
    c_ops = []
    if gamma_down > 0.0:
        c_ops.append(np.sqrt(gamma_down * (1.0 + n_th_q)) * sm)
        if n_th_q > 0.0:
            c_ops.append(np.sqrt(gamma_down * n_th_q) * sm.dag())
    if gamma_phi > 0.0:
        c_ops.append(np.sqrt(gamma_phi / 2.0) * sz)
    return c_ops


# ---------------------------------------------------------------------------
# Spectroscopy helpers
# ---------------------------------------------------------------------------

def spectroscopy_signal(
    omega_probe: np.ndarray,
    pn: np.ndarray,
    omega_q: float,
    chi: float,
    gamma_q: float,
) -> np.ndarray:
    """
    Synthetic number-selective spectroscopy signal.

    Models a weak continuous drive on the ancilla qubit.  In steady state,
    the qubit excitation probability is a sum of Lorentzians, one per Fock
    level n, weighted by P_n:

        S(ω) = Σ_n P_n * (γ_q/2)² / ((ω - ω_q^(n))² + (γ_q/2)²)

    where ω_q^(n) = ω_q + n * χ.

    This approximation holds in the weak-drive, number-resolved limit
    (χ >> γ_q), which is well satisfied here (χ T2 >> 1).

    Parameters
    ----------
    omega_probe : np.ndarray
        Probe frequencies (rad/s).
    pn : np.ndarray
        Photon-number distribution P_n for n = 0, 1, ..., len(pn)-1.
    omega_q : float
        Qubit transition frequency (rad/s).  Use 0 for rotating frame.
    chi : float
        Dispersive coupling (rad/s); sign sets direction of peak shift.
    gamma_q : float
        Qubit linewidth (rad/s) = 1/T2_q.
    """
    omega_probe = np.asarray(omega_probe)
    signal = np.zeros_like(omega_probe, dtype=float)
    half_gamma = gamma_q / 2.0
    for n, pn_val in enumerate(pn):
        if pn_val < 1e-12:
            continue
        omega_n = omega_q + n * chi
        signal += pn_val * half_gamma**2 / ((omega_probe - omega_n)**2 + half_gamma**2)
    return signal


# ---------------------------------------------------------------------------
# Truncation check
# ---------------------------------------------------------------------------

def check_truncation(n_bar: float, N: int, tol: float = 1e-4) -> bool:
    """
    Check whether Fock-space truncation N is adequate for mean occupation n_bar.

    Computes how much probability mass lies beyond n = N-1 and warns if it
    exceeds tol.

    Returns True if truncation is too tight, False if it is adequate.
    """
    pn = thermal_pn(n_bar, N)
    tail = max(0.0, 1.0 - pn.sum())
    if tail > tol:
        print(
            f"  [TRUNCATION WARNING] N={N} clips {100*tail:.4g}% of P_n "
            f"for n_bar={n_bar:.4f}.  Consider N >= {int(np.ceil(5 * n_bar)) + 10}."
        )
        return True
    return False


# ---------------------------------------------------------------------------
# Default three-bath configuration
# ---------------------------------------------------------------------------

def default_baths(
    n_target: float = 1.0,
    kappa_frac_target: float = KAPPA_FRAC_TARGET,
    n_bg: float = N_TH_BG,
    kappa_frac_bg: float = KAPPA_FRAC_BG,
    n_int: float = N_TH_INT,
    kappa_frac_int: float = KAPPA_FRAC_INT,
    kappa_tot: float = KAPPA_TOT,
) -> list:
    """
    Return the default three-bath configuration for the sensing study.

    Baths:
      - target:     target device (e.g. attenuator), variable occupation
      - background: cold background bath (dilution fridge environment)
      - internal:   internal cavity loss (zero temperature)
    """
    return [
        ThermalBath(kappa=kappa_frac_target * kappa_tot, n_th=n_target,
                    label="target"),
        ThermalBath(kappa=kappa_frac_bg * kappa_tot,     n_th=n_bg,
                    label="background"),
        ThermalBath(kappa=kappa_frac_int * kappa_tot,    n_th=n_int,
                    label="internal"),
    ]
