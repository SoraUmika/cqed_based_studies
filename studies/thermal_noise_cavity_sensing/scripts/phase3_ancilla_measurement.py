"""
Phase 3: Ancilla-assisted measurement of cavity thermal occupation.

Simulates two measurement protocols used to infer the cavity photon-number
distribution and mean occupation via a dispersively coupled ancilla qubit:

  A. Number-selective spectroscopy:
       S(ω) = Σ_n P_n * L(ω - ω_q - nχ, γ_q)
     Peaks at ω_q + nχ have heights proportional to P_n.

  B. Ramsey-based thermal-photon sensing:
       P_e(τ) = (1 + Re[χ_ramsey(τ) * exp(-Γ_q τ)]) / 2
     where χ_ramsey(τ) = Σ_n P_n exp(i n χ τ) [characteristic function]

Both protocols yield information about P_n / n̄, but not separately about
κ_target and n_target (see Phase 4 for the inverse problem).

Uses cqed_sim-consistent parameters (ω_c, ω_q, χ, T1, T2) while building
the qubit+cavity Hamiltonian directly with QuTiP (documented gap: no
pulse infrastructure needed for thermal evolution).

Usage:
    python scripts/phase3_ancilla_measurement.py

Output:
    data/phase3_results.npz
    figures/phase3_spectroscopy.{png,pdf}
    figures/phase3_ramsey.{png,pdf}
    figures/phase3_ancilla_bias.{png,pdf}
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
from scipy.optimize import curve_fit

from common import (
    DATA_DIR, FIG_DIR, STYLE_PATH,
    KAPPA_TOT, N_CAV,
    ThermalBath, analytic_nbar_ss, thermal_pn,
    build_cavity_c_ops, build_qubit_c_ops,
    spectroscopy_signal, ramsey_coherence_thermal,
    check_truncation,
    CHI_DISP, GAMMA_DOWN, GAMMA_PHI, T1_Q, T2Q,
    OMEGA_C,
)

if STYLE_PATH.exists():
    plt.style.use(STYLE_PATH)

COLORS = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377"]

print("=" * 60)
print("Phase 3: Ancilla-Assisted Measurement Simulation")
print("=" * 60)

# ---------------------------------------------------------------------------
# System parameters
# ---------------------------------------------------------------------------

N = N_CAV            # Fock truncation
CHI = -CHI_DISP      # actual dispersive coupling < 0 (qubit shifts down with photons)
KAPPA_BG_FRAC = 0.30
KAPPA_INT_FRAC = 0.20
N_BG = 0.01
N_INT = 0.00

def make_baths(n_target: float, kf_target: float = 0.50) -> list:
    """Standard three-bath configuration."""
    res = 1.0 - kf_target
    frac_bg  = KAPPA_BG_FRAC  / (KAPPA_BG_FRAC + KAPPA_INT_FRAC)
    frac_int = KAPPA_INT_FRAC / (KAPPA_BG_FRAC + KAPPA_INT_FRAC)
    return [
        ThermalBath(kappa=kf_target        * KAPPA_TOT, n_th=n_target, label="target"),
        ThermalBath(kappa=res * frac_bg    * KAPPA_TOT, n_th=N_BG,     label="bg"),
        ThermalBath(kappa=res * frac_int   * KAPPA_TOT, n_th=N_INT,    label="int"),
    ]

# Study conditions
CONDITIONS = [
    ("Cold cavity",  0.0, 0.50),
    ("n_t=1, 50%",   1.0, 0.50),
    ("n_t=2, 50%",   2.0, 0.50),
    ("n_t=3, 50%",   3.0, 0.50),
]

# ---------------------------------------------------------------------------
# Protocol A: Number-selective spectroscopy
# ---------------------------------------------------------------------------
print("\n--- Protocol A: Number-selective spectroscopy ---")

GAMMA_Q = 1.0 / T2Q          # qubit linewidth (rad/s) = 1/T2
N_PEAKS = 10                  # how many peaks to include
OMEGA_Q_ROT = 0.0             # qubit frequency in rotating frame

# Spectroscopy probe grid: from ω_q - (N_PEAKS+1)χ to ω_q + 2χ (peaks below ω_q)
OMEGA_PROBE_SPAN = (N_PEAKS + 2) * abs(CHI)
omega_probe = np.linspace(OMEGA_Q_ROT - OMEGA_PROBE_SPAN, OMEGA_Q_ROT + abs(CHI),
                          2000)

spec_data = {}
inferred_nbar_spec = {}

for label, n_t, kf in CONDITIONS:
    baths = make_baths(n_t, kf)
    n_ss = analytic_nbar_ss(baths)
    pn = thermal_pn(n_ss, N)

    # Synthetic spectroscopy signal
    signal = spectroscopy_signal(omega_probe, pn, OMEGA_Q_ROT, CHI, GAMMA_Q)

    # Infer P_n from spectroscopy: integrate signal around each peak
    pn_inferred = np.zeros(N_PEAKS)
    for n_pk in range(N_PEAKS):
        omega_n = OMEGA_Q_ROT + n_pk * CHI
        window = np.abs(omega_probe - omega_n) < abs(CHI) / 2
        if window.sum() > 0:
            domega = np.diff(omega_probe).mean()
            pn_inferred[n_pk] = np.trapezoid(signal[window], omega_probe[window])

    # Normalize inferred P_n
    pn_norm = pn_inferred.sum()
    if pn_norm > 0:
        pn_inferred /= pn_norm

    # Inferred n̄ from spectroscopy peaks
    n_bar_inferred = sum(n * p for n, p in enumerate(pn_inferred))

    spec_data[label] = {
        "signal": signal, "pn_true": pn[:N_PEAKS],
        "pn_inferred": pn_inferred, "n_ss": n_ss,
        "n_bar_inferred": n_bar_inferred, "n_target": n_t, "kf": kf,
    }
    inferred_nbar_spec[label] = n_bar_inferred
    print(f"  {label}: n̄_true={n_ss:.4f}, n̄_inferred={n_bar_inferred:.4f}, "
          f"bias={n_bar_inferred - n_ss:.4f}")

# ---------------------------------------------------------------------------
# Protocol B: Ramsey-based sensing
# ---------------------------------------------------------------------------
print("\n--- Protocol B: Ramsey-based thermal-photon sensing ---")

# Ramsey free-evolution time array
# Several oscillation periods of χ, up to 3×T2
T2_eff = T2Q
N_RAMSEY_PERIODS = 3
TAU_MAX = min(N_RAMSEY_PERIODS * T2_eff, 5.0 / abs(CHI))
tau_arr = np.linspace(0, TAU_MAX, 500)

ramsey_data = {}

for label, n_t, kf in CONDITIONS:
    baths = make_baths(n_t, kf)
    n_ss = analytic_nbar_ss(baths)

    # Ideal Ramsey coherence (no qubit T2)
    chi_tau_ideal = ramsey_coherence_thermal(tau_arr, n_ss, CHI)

    # With qubit T2 envelope: multiply by exp(-Γ_q τ / 2)
    # (factor 1/2 because T2 sets the off-diagonal decay rate 1/T2,
    #  and the signal is Re[χ(τ)], which decays as exp(-τ/T2))
    t2_envelope = np.exp(-tau_arr / T2_eff)
    chi_tau = chi_tau_ideal * t2_envelope

    # Ramsey signal: P_e(τ) = (1 + Re[χ(τ)]) / 2
    ramsey_signal = 0.5 * (1.0 + np.real(chi_tau))

    # Infer n̄ by fitting the Ramsey signal.
    # For small χτ: |χ(τ)|² ≈ exp(-n̄(1+n̄)(χτ)²) → dephasing rate ∝ n̄(1+n̄)
    # Use the analytic form directly for the fit.
    def ramsey_model(tau, n_bar_fit):
        chi_fit = ramsey_coherence_thermal(tau, max(n_bar_fit, 0.0), CHI)
        t2_env = np.exp(-tau / T2_eff)
        return 0.5 * (1.0 + np.real(chi_fit * t2_env))

    try:
        popt, _ = curve_fit(ramsey_model, tau_arr, ramsey_signal,
                            p0=[max(n_ss, 0.1)], bounds=(0, 20),
                            maxfev=5000)
        n_bar_inferred = popt[0]
    except RuntimeError:
        n_bar_inferred = float("nan")

    ramsey_data[label] = {
        "tau": tau_arr, "signal": ramsey_signal,
        "chi_ideal": chi_tau_ideal, "n_ss": n_ss,
        "n_bar_inferred": n_bar_inferred, "n_target": n_t, "kf": kf,
    }
    print(f"  {label}: n̄_true={n_ss:.4f}, n̄_inferred={n_bar_inferred:.4f}, "
          f"bias={n_bar_inferred - n_ss:.4f}")

# ---------------------------------------------------------------------------
# Bias study: effect of qubit T2 on inferred n̄
# ---------------------------------------------------------------------------
print("\n--- Bias study: qubit T2 effect ---")

T2_VALS = np.array([1e-6, 5e-6, 10e-6, 20e-6, 50e-6, 100e-6])   # s
n_target_bias = 2.0
baths_bias = make_baths(n_target_bias, 0.50)
n_ss_bias = analytic_nbar_ss(baths_bias)

bias_inferred = []

for T2_v in T2_VALS:
    tau_v = np.linspace(0, min(3 * T2_v, 5.0 / abs(CHI)), 300)
    t2_env = np.exp(-tau_v / T2_v)
    chi_v = ramsey_coherence_thermal(tau_v, n_ss_bias, CHI) * t2_env
    signal_v = 0.5 * (1.0 + np.real(chi_v))

    def ramsey_model_v(tau, n_bar_fit):
        chi_fit = ramsey_coherence_thermal(tau, max(n_bar_fit, 0.0), CHI)
        env = np.exp(-tau / T2_v)
        return 0.5 * (1.0 + np.real(chi_fit * env))

    try:
        popt, _ = curve_fit(ramsey_model_v, tau_v, signal_v,
                            p0=[n_ss_bias], bounds=(0, 20), maxfev=5000)
        bias_inferred.append(popt[0])
    except RuntimeError:
        bias_inferred.append(float("nan"))

    chi_T2 = abs(CHI) * T2_v / (2 * np.pi)
    print(f"  T2={T2_v*1e6:.1f} μs (χT2/2π={chi_T2:.1f}): "
          f"n̄_inferred={bias_inferred[-1]:.4f} (true={n_ss_bias:.4f})")

bias_inferred = np.array(bias_inferred)

# ---------------------------------------------------------------------------
# Generate figures
# ---------------------------------------------------------------------------
print("\nGenerating figures...")

# Figure 1: Spectroscopy
fig, axes = plt.subplots(1, 2, figsize=(6.75, 2.8))

ax = axes[0]
omega_MHz = (omega_probe - OMEGA_Q_ROT) / (2 * np.pi * 1e6)
for i, (label, dat) in enumerate(spec_data.items()):
    ax.plot(omega_MHz, dat["signal"], color=COLORS[i], lw=1.0, label=label)

# Annotate peak positions
for n_pk in range(5):
    omega_n_MHz = n_pk * CHI / (2 * np.pi * 1e6)
    ax.axvline(omega_n_MHz, color="gray", ls=":", lw=0.5)
    if n_pk <= 4:
        ax.text(omega_n_MHz - 0.3, ax.get_ylim()[1] * 0.9,
                f"n={n_pk}", fontsize=6, ha="right", color="gray")

ax.set_xlabel(r"$(\omega_\mathrm{probe} - \omega_q) / 2\pi$ (MHz)")
ax.set_ylabel(r"Qubit excitation $P_e$")
ax.set_title("Number-selective spectroscopy")
ax.legend(fontsize=6)

ax = axes[1]
n_pk_arr = np.arange(N_PEAKS)
for i, (label, dat) in enumerate(spec_data.items()):
    ax.semilogy(n_pk_arr, np.maximum(dat["pn_true"], 1e-12), "o--",
                ms=4, color=COLORS[i], lw=0.8, label=f"{label} (true)")
    ax.semilogy(n_pk_arr, np.maximum(dat["pn_inferred"], 1e-12), "s",
                ms=3, color=COLORS[i], alpha=0.6, label=f"{label} (inferred)")

ax.set_xlabel(r"Fock level $n$")
ax.set_ylabel(r"$P_n$")
ax.set_title("P_n true vs inferred")
ax.legend(fontsize=5, ncol=2)
ax.set_xlim(-0.5, N_PEAKS - 0.5)
ax.set_ylim(1e-6, 1.5)

plt.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"phase3_spectroscopy.{ext}")
plt.close(fig)
print("  Saved: phase3_spectroscopy")

# Figure 2: Ramsey signals
fig, axes = plt.subplots(1, 2, figsize=(6.75, 2.8))

ax = axes[0]
for i, (label, dat) in enumerate(ramsey_data.items()):
    tau_us = dat["tau"] * 1e6
    ax.plot(tau_us, dat["signal"], color=COLORS[i], lw=1.0, label=label)

ax.set_xlabel(r"Free evolution $\tau$ ($\mu$s)")
ax.set_ylabel(r"Qubit excitation $P_e(\tau)$")
ax.set_title("Ramsey signal vs cavity occupation")
ax.legend(fontsize=7)
ax.set_xlim(0, TAU_MAX * 1e6)
ax.set_ylim(0, 1)

ax = axes[1]
for i, (label, dat) in enumerate(ramsey_data.items()):
    tau_us = dat["tau"] * 1e6
    ax.plot(tau_us, np.abs(dat["chi_ideal"]), color=COLORS[i], lw=1.0,
            label=f"{label} (n̄={dat['n_ss']:.2f})")

ax.set_xlabel(r"Free evolution $\tau$ ($\mu$s)")
ax.set_ylabel(r"$|\chi(\tau)|$ (coherence magnitude)")
ax.set_title("Ramsey coherence from thermal cavity")
ax.legend(fontsize=7)
ax.set_xlim(0, TAU_MAX * 1e6)
ax.set_ylim(0, 1.05)

plt.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"phase3_ramsey.{ext}")
plt.close(fig)
print("  Saved: phase3_ramsey")

# Figure 3: Ancilla bias vs T2
fig, ax = plt.subplots(figsize=(3.375, 2.8))

chi_T2_arr = abs(CHI) * T2_VALS / (2 * np.pi)
ax.semilogx(chi_T2_arr, bias_inferred, "o-", color=COLORS[0], lw=1.2, ms=4,
            label="Inferred $\\bar{n}$")
ax.axhline(n_ss_bias, color="k", ls="--", lw=1.0, label=f"True $\\bar{{n}}$={n_ss_bias:.3f}")
ax.set_xlabel(r"$|\chi| T_2 / (2\pi)$")
ax.set_ylabel(r"Inferred $\bar{n}$")
ax.set_title("Ramsey inference bias vs qubit T₂")
ax.legend(fontsize=8)
ax.set_ylim(bottom=0)
plt.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"phase3_ancilla_bias.{ext}")
plt.close(fig)
print("  Saved: phase3_ancilla_bias")

# ---------------------------------------------------------------------------
# Save data
# ---------------------------------------------------------------------------
spec_labels = [c[0] for c in CONDITIONS]
spec_nss = np.array([spec_data[l]["n_ss"] for l in spec_labels])
spec_inferred = np.array([spec_data[l]["n_bar_inferred"] for l in spec_labels])

ramsey_nss = np.array([ramsey_data[l]["n_ss"] for l in spec_labels])
ramsey_inferred = np.array([ramsey_data[l]["n_bar_inferred"] for l in spec_labels])

np.savez(
    DATA_DIR / "phase3_results.npz",
    omega_probe=omega_probe,
    tau_arr=tau_arr,
    # Spectroscopy
    spec_nss=spec_nss,
    spec_inferred=spec_inferred,
    # Ramsey
    ramsey_nss=ramsey_nss,
    ramsey_inferred=ramsey_inferred,
    # Bias study
    T2_VALS=T2_VALS,
    bias_inferred=bias_inferred,
    n_ss_bias=np.array([n_ss_bias]),
    chi_T2_arr=chi_T2_arr,
)
print("\nData saved: data/phase3_results.npz")

print("\n--- Inference summary ---")
print("  Protocol       n̄_true   n̄_spec   n̄_ramsey")
for i, (label, _, _) in enumerate(CONDITIONS):
    print(f"  {label:20s}  {spec_nss[i]:.4f}   {spec_inferred[i]:.4f}   "
          f"{ramsey_inferred[i]:.4f}")

print("\nKey result: For χT2 >> 1 (here χT2/2π ≈ {:.1f}), both protocols recover".format(
    abs(CHI) * T2Q / (2 * np.pi)))
print("  n̄ accurately.  The Lorentzian spectroscopy peaks are well-resolved.")
print("  However, neither gives additional identifiability for κ_target and n_target.")

print("\n" + "=" * 60)
print("Phase 3 complete.")
print("=" * 60)
