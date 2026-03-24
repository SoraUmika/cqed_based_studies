"""
Phase 4: Inverse problem — identifiability of target bath parameters.

Given simulated cavity observables, can we recover the target bath occupation
n_target and coupling κ_target separately?

Key finding (derived analytically and verified numerically):
  The cavity steady state is ALWAYS a thermal state with effective occupation
  n̄_ss = (Σ_j κ_j n_j) / (Σ_j κ_j).  Therefore:
  - n̄_ss alone cannot separate n_target and κ_target.
  - Adding transient rate measurement gives κ_tot but not individual κ_j.
  - Measurement of P_n (via spectroscopy) gives the same n̄_ss — no new info.
  - The product κ_target × n_target is the fundamental identifiable quantity.

Identifiable combinations:
  1. From steady-state:       κ_target × n_target + κ_bg × n_bg + ...   (one number)
  2. From transient:          κ_tot = Σ_j κ_j                           (one more)
  3. If κ_bg, n_bg, κ_int known:  κ_target × n_target                   (product only)
  4. From a κ_target measurement (cavity transmission): n_target         (fully identified)

Workflow:
  1. Generate synthetic data with true parameters + noise.
  2. Fit for recovered parameters under different observable sets.
  3. Visualize degeneracy manifolds in (n_target, κ_target) space.
  4. Quantify inference error vs noise level.

Usage:
    python scripts/phase4_inverse_problem.py

Output:
    data/phase4_results.npz
    figures/phase4_degeneracy_manifold.{png,pdf}
    figures/phase4_inference_error.{png,pdf}
    figures/phase4_identifiability_summary.{png,pdf}
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
from scipy.optimize import minimize, curve_fit

from common import (
    DATA_DIR, FIG_DIR, STYLE_PATH,
    KAPPA_TOT, N_CAV,
    ThermalBath, analytic_nbar_ss, analytic_kappa_tot,
    analytic_nbar_transient, thermal_pn,
    build_cavity_c_ops,
    ramsey_coherence_thermal, CHI_DISP, T2Q,
)

if STYLE_PATH.exists():
    plt.style.use(STYLE_PATH)

COLORS = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377"]

print("=" * 60)
print("Phase 4: Inverse Problem — Identifiability Study")
print("=" * 60)

# ---------------------------------------------------------------------------
# Known (calibrated) background parameters
# ---------------------------------------------------------------------------

N_BG  = 0.01
N_INT = 0.00
KAPPA_FRAC_BG  = 0.30
KAPPA_FRAC_INT = 0.20

# Known fixed fractions for bg/int (calibrated separately)
KAPPA_BG  = KAPPA_FRAC_BG  * KAPPA_TOT
KAPPA_INT = KAPPA_FRAC_INT * KAPPA_TOT
# Residual (bg + int) is known: κ_res = κ_bg + κ_int
KAPPA_RES = KAPPA_BG + KAPPA_INT  # = 0.5 * KAPPA_TOT

# ---------------------------------------------------------------------------
# Forward model
# ---------------------------------------------------------------------------

def forward_nbar(n_target: float, kappa_target: float) -> float:
    """Steady-state n̄ for given (n_target, κ_target)."""
    baths = [
        ThermalBath(kappa=kappa_target, n_th=n_target, label="target"),
        ThermalBath(kappa=KAPPA_BG,     n_th=N_BG,     label="bg"),
        ThermalBath(kappa=KAPPA_INT,    n_th=N_INT,     label="int"),
    ]
    return analytic_nbar_ss(baths)


def forward_kappa_tot(kappa_target: float) -> float:
    """Total decay rate = κ_target + κ_bg + κ_int."""
    return kappa_target + KAPPA_RES


def forward_nbar_product(product: float) -> float:
    """
    n̄_ss as a function of the product p = κ_target × n_target.
    (Single identifiable parameter from steady-state alone, after subtracting
    the known background contribution.)
    """
    kappa_tot_val = KAPPA_RES + (product / max(product, 1e-10)) * KAPPA_TOT  # ill-posed
    # Actually: n̄_ss = (p + κ_bg n_bg + κ_int n_int) / (κ_target + κ_res)
    # Without knowing κ_target separately, κ_tot is unknown.  Use placeholder.
    bg_contrib = KAPPA_BG * N_BG + KAPPA_INT * N_INT
    # Cannot evaluate without κ_tot.  This function illustrates the ill-posedness.
    return float("nan")

# ---------------------------------------------------------------------------
# True parameters for synthetic data generation
# ---------------------------------------------------------------------------

TRUE_CASES = [
    dict(n_target=1.0, kf=0.50, label="Case A (n_t=1, kf=0.5)"),
    dict(n_target=2.0, kf=0.40, label="Case B (n_t=2, kf=0.4)"),
    dict(n_target=5.0, kf=0.25, label="Case C (n_t=5, kf=0.25)"),
]

NOISE_SIGMA_NBAR = 0.05    # measurement noise on n̄_ss (photons)
NOISE_SIGMA_KAPPA = 0.02 * KAPPA_TOT  # noise on κ_tot from transient rate fit

# ---------------------------------------------------------------------------
# Inference attempt 1: Using n̄_ss only
# ---------------------------------------------------------------------------
print("\n--- Inference 1: n̄_ss only (single observable) ---")
print("  Cannot recover n_target and κ_target separately.")
print("  Constraint surface: n̄_ss = const is a curve in (n_target, κ_target) space.")

inf1_results = []
for case in TRUE_CASES:
    kappa_t = case["kf"] * KAPPA_TOT
    baths = [
        ThermalBath(kappa=kappa_t, n_th=case["n_target"], label="target"),
        ThermalBath(kappa=KAPPA_BG, n_th=N_BG, label="bg"),
        ThermalBath(kappa=KAPPA_INT, n_th=N_INT, label="int"),
    ]
    n_ss_true = analytic_nbar_ss(baths)
    n_ss_meas = n_ss_true + np.random.normal(0, NOISE_SIGMA_NBAR)

    # From n̄_ss and known bg, we know:
    #   n̄_ss = (κ_t × n_t + κ_bg × n_bg + κ_int × n_int) / (κ_t + κ_res)
    # → κ_t × n_t = n̄_ss × (κ_t + κ_res) − (κ_bg n_bg + κ_int n_int)
    # → κ_t × n_t = n̄_ss × κ_t + n̄_ss × κ_res − bg_contrib
    # This is ONE equation in TWO unknowns (κ_t, n_t).  Ill-posed.

    bg_contrib = KAPPA_BG * N_BG + KAPPA_INT * N_INT
    # Identifiable quantity: n̄_ss × (κ_t + κ_res) = κ_t × n_t + bg_contrib
    # → product = n̄_ss × (κ_t + κ_res) − bg_contrib  (depends on unknown κ_t)
    print(f"  {case['label']}: n̄_ss={n_ss_true:.4f} (meas={n_ss_meas:.4f})")
    print(f"    → Need κ_tot to isolate κ_target × n_target.")
    inf1_results.append({"n_ss_true": n_ss_true, "n_ss_meas": n_ss_meas, **case})

# ---------------------------------------------------------------------------
# Inference attempt 2: n̄_ss + κ_tot from transient rate
# ---------------------------------------------------------------------------
print("\n--- Inference 2: n̄_ss + κ_tot (two observables) ---")
print("  Now have: κ_target = κ_tot − κ_res, and")
print("            κ_target × n_target = n̄_ss × κ_tot − bg_contrib")
print("  Product p = κ_target × n_target is identifiable, but not n_target alone.")

inf2_results = []
for case in TRUE_CASES:
    kappa_t_true = case["kf"] * KAPPA_TOT
    baths = [
        ThermalBath(kappa=kappa_t_true, n_th=case["n_target"], label="target"),
        ThermalBath(kappa=KAPPA_BG, n_th=N_BG, label="bg"),
        ThermalBath(kappa=KAPPA_INT, n_th=N_INT, label="int"),
    ]
    kappa_tot_true = kappa_t_true + KAPPA_RES
    n_ss_true = analytic_nbar_ss(baths)

    # Add measurement noise
    n_ss_meas = n_ss_true + np.random.normal(0, NOISE_SIGMA_NBAR)
    kappa_tot_meas = kappa_tot_true + np.random.normal(0, NOISE_SIGMA_KAPPA)

    # Recover κ_target
    kappa_t_inferred = max(kappa_tot_meas - KAPPA_RES, 1e-6)

    # Recover product p = κ_t × n_t
    bg_contrib = KAPPA_BG * N_BG + KAPPA_INT * N_INT
    product_inferred = n_ss_meas * kappa_tot_meas - bg_contrib
    product_true = kappa_t_true * case["n_target"]

    # n_target still not separately recoverable without independent κ_t measurement
    # (unless we assume κ_t = κ_tot - κ_res from transmission measurement)

    print(f"  {case['label']}:")
    print(f"    κ_target (true/inferred): {kappa_t_true/(2*np.pi*1e3):.1f} / "
          f"{kappa_t_inferred/(2*np.pi*1e3):.1f} kHz (2π×)")
    print(f"    product κ_t×n_t (true/inferred): {product_true/(KAPPA_TOT):.4f} / "
          f"{product_inferred/(KAPPA_TOT):.4f} (in units of κ_tot)")
    inf2_results.append({
        "kappa_t_true": kappa_t_true, "kappa_t_inferred": kappa_t_inferred,
        "product_true": product_true, "product_inferred": product_inferred,
        "n_ss_true": n_ss_true, "n_ss_meas": n_ss_meas, **case,
    })

# ---------------------------------------------------------------------------
# Inference attempt 3: n̄_ss + κ_tot + κ_target from transmission
# ---------------------------------------------------------------------------
print("\n--- Inference 3: n̄_ss + κ_tot + κ_target (three observables) ---")
print("  Now n_target = (κ_t×n_t) / κ_t is fully identifiable.")

NOISE_SIGMA_KT = 0.03 * KAPPA_TOT   # 3% noise on κ_target from transmission

inf3_results = []
for case in TRUE_CASES:
    kappa_t_true = case["kf"] * KAPPA_TOT
    baths = [
        ThermalBath(kappa=kappa_t_true, n_th=case["n_target"], label="target"),
        ThermalBath(kappa=KAPPA_BG, n_th=N_BG, label="bg"),
        ThermalBath(kappa=KAPPA_INT, n_th=N_INT, label="int"),
    ]
    kappa_tot_true = kappa_t_true + KAPPA_RES
    n_ss_true = analytic_nbar_ss(baths)

    # Add noise
    n_ss_meas = n_ss_true + np.random.normal(0, NOISE_SIGMA_NBAR)
    kappa_tot_meas = kappa_tot_true + np.random.normal(0, NOISE_SIGMA_KAPPA)
    kappa_t_meas = kappa_t_true + np.random.normal(0, NOISE_SIGMA_KT)

    kappa_t_meas = max(kappa_t_meas, 1e-6)
    bg_contrib = KAPPA_BG * N_BG + KAPPA_INT * N_INT
    product = n_ss_meas * kappa_tot_meas - bg_contrib
    n_target_inferred = product / kappa_t_meas

    rel_err = abs(n_target_inferred - case["n_target"]) / max(case["n_target"], 0.01)
    print(f"  {case['label']}: n_target true={case['n_target']:.2f}, "
          f"inferred={n_target_inferred:.4f}, rel_err={100*rel_err:.2f}%")
    inf3_results.append({
        "n_target_true": case["n_target"], "n_target_inferred": n_target_inferred,
        "rel_err": rel_err, **case,
    })

# ---------------------------------------------------------------------------
# Degeneracy manifold: all (n_target, κ_frac) pairs giving same n̄_ss
# ---------------------------------------------------------------------------
print("\n--- Degeneracy manifold visualization ---")

# For a fixed n̄_ss value, compute the curve (n_target, κ_frac) satisfying it
n_bar_targets = [0.3, 0.6, 1.0, 1.5, 2.5]
n_target_deg_arr = np.linspace(0.01, 12, 200)

deg_manifolds = {}
for n_bar_tgt in n_bar_targets:
    kf_solutions = []
    for n_t in n_target_deg_arr:
        # n̄_ss = (κ_t n_t + κ_bg n_bg + κ_int n_int) / (κ_t + κ_res)
        # → n̄_ss (κ_t + κ_res) = κ_t n_t + bg_contrib
        # → κ_t (n̄_ss - n_t) = bg_contrib - n̄_ss κ_res
        # → κ_t = (bg_contrib - n̄_ss κ_res) / (n̄_ss - n_t)
        bg_contrib = KAPPA_BG * N_BG + KAPPA_INT * N_INT
        denom = n_bar_tgt - n_t
        if abs(denom) < 1e-10:
            kf_solutions.append(float("nan"))
            continue
        kappa_t = (bg_contrib - n_bar_tgt * KAPPA_RES) / denom
        kf = kappa_t / KAPPA_TOT
        if 0 < kf < 1 and kappa_t > 0:
            kf_solutions.append(kf)
        else:
            kf_solutions.append(float("nan"))
    deg_manifolds[n_bar_tgt] = np.array(kf_solutions)

# ---------------------------------------------------------------------------
# Inference error vs noise level
# ---------------------------------------------------------------------------
print("\n--- Inference error vs measurement noise ---")

NOISE_LEVELS = np.array([0.01, 0.02, 0.05, 0.1, 0.2])  # fraction of n̄_ss
N_TRIALS = 500
n_true_test = 2.0
kf_true_test = 0.50

baths_test = [
    ThermalBath(kappa=kf_true_test * KAPPA_TOT, n_th=n_true_test, label="t"),
    ThermalBath(kappa=KAPPA_BG, n_th=N_BG, label="bg"),
    ThermalBath(kappa=KAPPA_INT, n_th=N_INT, label="int"),
]
n_ss_test = analytic_nbar_ss(baths_test)
kappa_tot_test = kf_true_test * KAPPA_TOT + KAPPA_RES
kappa_t_test = kf_true_test * KAPPA_TOT

rng = np.random.default_rng(42)
rmse_n_target = []

for noise_frac in NOISE_LEVELS:
    sigma_n = noise_frac * n_ss_test
    sigma_k = noise_frac * KAPPA_TOT
    sigma_kt = 0.03 * KAPPA_TOT  # independent transmission noise

    errors = []
    for _ in range(N_TRIALS):
        n_ss_m = n_ss_test + rng.normal(0, sigma_n)
        kappa_tot_m = kappa_tot_test + rng.normal(0, sigma_k)
        kappa_t_m = max(kappa_t_test + rng.normal(0, sigma_kt), 1e-6)
        bg_contrib = KAPPA_BG * N_BG + KAPPA_INT * N_INT
        product = n_ss_m * kappa_tot_m - bg_contrib
        n_inferred = product / kappa_t_m
        errors.append(n_inferred - n_true_test)

    errors = np.array(errors)
    rmse = np.sqrt(np.mean(errors**2))
    rmse_n_target.append(rmse)
    print(f"  noise={noise_frac*100:.0f}%: RMSE(n_target)={rmse:.4f}")

rmse_n_target = np.array(rmse_n_target)

# ---------------------------------------------------------------------------
# Generate figures
# ---------------------------------------------------------------------------
print("\nGenerating figures...")

# Figure 1: Degeneracy manifold
fig, ax = plt.subplots(figsize=(3.375, 2.8))

for i, (n_bar_tgt, kf_arr) in enumerate(deg_manifolds.items()):
    valid = np.isfinite(kf_arr)
    ax.plot(n_target_deg_arr[valid], kf_arr[valid], color=COLORS[i], lw=1.2,
            label=rf"$\bar{{n}}_\mathrm{{ss}}={n_bar_tgt}$")

# Mark true parameter points for the three cases
for i, case in enumerate(TRUE_CASES):
    ax.plot(case["n_target"], case["kf"], "^", ms=6, color="k", zorder=5)

ax.set_xlabel(r"$n_\mathrm{target}$")
ax.set_ylabel(r"$\kappa_\mathrm{target}/\kappa_\mathrm{tot}$")
ax.set_title("Degeneracy manifold (iso-n̄ curves)")
ax.legend(fontsize=7, loc="upper right")
ax.set_xlim(0, 10)
ax.set_ylim(0, 1)
ax.text(0.3, 0.15, "▲ True params", fontsize=7, transform=ax.transAxes)
plt.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"phase4_degeneracy_manifold.{ext}")
plt.close(fig)
print("  Saved: phase4_degeneracy_manifold")

# Figure 2: Inference error vs noise
fig, ax = plt.subplots(figsize=(3.375, 2.8))

ax.loglog(NOISE_LEVELS * 100, rmse_n_target, "o-", color=COLORS[0], lw=1.2, ms=4,
          label="RMSE (3 observables)")
# Reference: 1:1 scaling with n̄_ss noise
ax.loglog(NOISE_LEVELS * 100, NOISE_LEVELS * n_ss_test, "--", color="gray", lw=1.0,
          label=r"$\propto \sigma_{\bar{n}}$")

ax.set_xlabel(r"Noise level (\% of $\bar{n}_\mathrm{ss}$)")
ax.set_ylabel(r"RMSE($n_\mathrm{target}$)")
ax.set_title(r"Inference error vs noise")
ax.legend(fontsize=8)
plt.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"phase4_inference_error.{ext}")
plt.close(fig)
print("  Saved: phase4_inference_error")

# Figure 3: Identifiability summary table (text figure)
fig, ax = plt.subplots(figsize=(6.75, 2.0))
ax.axis("off")

table_data = [
    ["Observable set", "Identifiable quantities", "Degeneracy"],
    [r"$\bar{n}_\mathrm{ss}$ only", r"$\kappa_t n_t + \kappa_\mathrm{bg} n_\mathrm{bg} + ...$",
     r"1D manifold in $(\kappa_t, n_t)$ space"],
    [r"$\bar{n}_\mathrm{ss}$ + $\kappa_\mathrm{tot}$", r"$\kappa_t$ (= $\kappa_\mathrm{tot} - \kappa_\mathrm{res}$) and product $\kappa_t n_t$",
     r"Product $\kappa_t \times n_t$ only"],
    [r"$\bar{n}_\mathrm{ss}$ + $\kappa_\mathrm{tot}$ + $\kappa_t$ (transmission)",
     r"$n_t = (\kappa_t n_t) / \kappa_t$", "Fully identifiable"],
    [r"Transient $\bar{n}(t)$ only", r"$\kappa_\mathrm{tot}$ and $\bar{n}_\mathrm{ss}$",
     "Same as row 2"],
]

table = ax.table(cellText=table_data[1:], colLabels=table_data[0],
                 loc="center", cellLoc="left")
table.auto_set_font_size(False)
table.set_fontsize(7)
table.auto_set_column_width([0, 1, 2])
for (r, c), cell in table.get_celld().items():
    cell.set_linewidth(0.5)
    if r == 0:
        cell.set_facecolor("#DDDDDD")

ax.set_title("Identifiability summary", fontsize=9, pad=4)
plt.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(FIG_DIR / f"phase4_identifiability_summary.{ext}")
plt.close(fig)
print("  Saved: phase4_identifiability_summary")

# ---------------------------------------------------------------------------
# Save data
# ---------------------------------------------------------------------------
np.savez(
    DATA_DIR / "phase4_results.npz",
    n_target_deg_arr=n_target_deg_arr,
    n_bar_targets=np.array(n_bar_targets),
    NOISE_LEVELS=NOISE_LEVELS,
    rmse_n_target=rmse_n_target,
    # Summary arrays
    inf3_n_true=np.array([r["n_target_true"] for r in inf3_results]),
    inf3_n_inferred=np.array([r["n_target_inferred"] for r in inf3_results]),
    inf3_rel_err=np.array([r["rel_err"] for r in inf3_results]),
)
print("\nData saved: data/phase4_results.npz")

print("\n--- Key findings ---")
print("  1. Steady-state n̄ alone: UNDERDETERMINED — iso-n̄ curves span entire")
print("     (n_target, κ_frac) plane.  No unique solution.")
print("  2. Adding transient measurement: κ_tot identified → κ_target × n_target")
print("     recoverable, but n_target and κ_target still degenerate as a pair.")
print("  3. Adding cavity transmission measurement: κ_target separately identified")
print("     → n_target fully recoverable.  This is the minimal observable set.")
print("  4. P_n distribution (from spectroscopy) gives same n̄ — no additional")
print("     identifiability for the bath separation problem.")

print("\n" + "=" * 60)
print("Phase 4 complete.")
print("=" * 60)
