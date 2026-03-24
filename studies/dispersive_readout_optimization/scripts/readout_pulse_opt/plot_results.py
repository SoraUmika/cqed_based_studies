"""
Plot Results — Publication-Quality Figures for All Phases
==========================================================

Generates 8 main figures:
  Fig 1: Steady-state separation |Δα_ss| vs χ/κ and drive detuning (Phase 1)
  Fig 2: SNR² vs κT — ODE vs steady-state approx, midpoint vs resonance drive (Phase 1)
  Fig 3: Objective comparison — endpoint sep vs integrated SNR² (Phase 1)
  Fig 4: Pulse-family comparison — SNR² vs χ/κ at κT=5 (Phase 2)
  Fig 5: Readout+depletion — residual photons vs depletion time (Phase 2)
  Fig 6: GRAPE vs baselines — SNR² across χ/κ and κT (Phase 3)
  Fig 7: Net fidelity vs readout duration with T1 effects (Phase 4)
  Fig 8: Robustness — sensitivity maps and stats comparison (Phase 5)

Uses colorblind-friendly Tol Bright palette consistent with existing studies.
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

from common import TWO_PI, assignment_fidelity_from_snr2

DATA_DIR    = os.path.join(os.path.dirname(__file__), "..", "data")
FIGURES_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)

# ── Tol Bright colorblind-safe palette ───────────────────────────────────────
TOL = {
    "blue":   "#4477AA",
    "cyan":   "#66CCEE",
    "green":  "#228833",
    "yellow": "#CCBB44",
    "red":    "#EE6677",
    "purple": "#AA3377",
    "grey":   "#BBBBBB",
    "black":  "#000000",
}
FAMILY_COLORS = {
    "square":      TOL["blue"],
    "gaussian":    TOL["green"],
    "hann":        TOL["red"],
    "cosine_rise": TOL["yellow"],
    "spline":      TOL["purple"],
}
FAMILY_LABELS = {
    "square":      "Square",
    "gaussian":    "Gaussian",
    "hann":        "Hann",
    "cosine_rise": "Cos-rise",
    "spline":      "Spline",
}
FAMILY_MARKERS = {
    "square": "o", "gaussian": "s", "hann": "^", "cosine_rise": "D", "spline": "v"
}

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "legend.fontsize": 10,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "axes.linewidth": 1.2,
    "figure.dpi": 150,
})


def savefig(name: str) -> None:
    for ext in ["png", "pdf"]:
        fpath = os.path.join(FIGURES_DIR, f"{name}.{ext}")
        plt.savefig(fpath, bbox_inches="tight")
    print(f"  Saved: {name}.{{png,pdf}}")
    plt.close()


# ─── Figure 1: Steady-state separation ───────────────────────────────────────
def fig_steady_state(d1: dict) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    # Panel (a): |Δα_ss_max| vs χ/κ
    ax = axes[0]
    ck = d1["chi_kappa_chi_kappa_scan"]
    sep = d1["chi_kappa_sep_max_analytical"]
    ax.plot(ck, sep, color=TOL["blue"], lw=2)
    ax.axvline(1.0, ls="--", color=TOL["grey"], lw=1.2, label=r"$\chi/\kappa = 1$")
    ax.set_xlabel(r"$\chi/\kappa$")
    ax.set_ylabel(r"$|\Delta\alpha^\mathrm{ss}_\mathrm{max}|$ (normalized to $|\varepsilon|/\kappa$)")
    ax.set_title(r"(a) Steady-state pointer separation")
    ax.legend()
    # Normalize y axis by ε/κ (the peak at χ/κ=1 is 2ε/κ)
    ax.set_xlim(0, 6)
    ax.grid(True, alpha=0.3)
    ax.annotate("peak at $\\chi/\\kappa=1$", xy=(1.0, sep[np.argmin(np.abs(ck - 1.0))]),
                xytext=(2.0, sep[np.argmin(np.abs(ck - 1.0))]*0.9),
                arrowprops=dict(arrowstyle="->", color="black"), fontsize=10)

    # Panel (b): |Δα_ss| vs drive detuning Δ_g/χ
    ax = axes[1]
    dg = d1["drive_freq_delta_g_values"] / abs(d1["kappa"] * 1.0)  # just plot in units of Δ_g/2π
    # Actually plot in MHz
    dg_MHz = d1["drive_freq_delta_g_values"] / TWO_PI / 1e6
    sep_dg = d1["drive_freq_sep"]
    ax.plot(dg_MHz, sep_dg, color=TOL["red"], lw=2, label=r"$|\Delta\alpha^{\rm ss}|$")
    dg_opt = d1["delta_g_opt"] / TWO_PI / 1e6
    ax.axvline(dg_opt, ls="--", color=TOL["green"], lw=1.5, label=f"Opt $\\Delta_g/2\\pi = {dg_opt:.1f}$ MHz")
    ax.set_xlabel(r"$\Delta_g / 2\pi$ (MHz)")
    ax.set_ylabel(r"$|\Delta\alpha^{\rm ss}(\Delta_g)|$")
    ax.set_title(r"(b) Drive frequency optimization ($\chi/\kappa = 1$)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    savefig("fig1_steady_state")


# ─── Figure 2: SNR vs κT ─────────────────────────────────────────────────────
def fig_snr_vs_kappa_t(d1: dict) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    kT = d1["snr_kT_kappa_t_scan"]
    ax = axes[0]
    ax.plot(kT, d1["snr_kT_snr2_ode_opt"],  color=TOL["blue"],  lw=2,  label="Midpoint drive (ODE)")
    ax.plot(kT, d1["snr_kT_snr2_ode_res"],  color=TOL["red"],   lw=2,  ls="--", label="On-resonance (ODE)")
    ax.plot(kT, d1["snr_kT_snr2_ss_approx"], color=TOL["grey"], lw=1.5, ls=":", label=r"SS approx: $\kappa T|\Delta\alpha^{\rm ss}|^2$")
    ax.set_xlabel(r"$\kappa T$")
    ax.set_ylabel(r"$\mathrm{SNR}^2$")
    ax.set_title(r"(a) Integrated SNR$^2$ vs readout duration")
    ax.legend(fontsize=9)
    ax.set_xlim(0, 20)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.plot(kT, d1["snr_kT_F_assign_ode_opt"], color=TOL["blue"], lw=2,    label="Midpoint drive")
    ax.plot(kT, d1["snr_kT_F_assign_ode_res"], color=TOL["red"],  lw=2, ls="--", label="On-resonance")
    ax.axhline(0.99, ls=":", color=TOL["grey"], lw=1.2)
    ax.axhline(0.999, ls=":", color=TOL["grey"], lw=1.2)
    ax.set_xlabel(r"$\kappa T$")
    ax.set_ylabel(r"Assignment fidelity $F_A$")
    ax.set_title(r"(b) Assignment fidelity vs readout duration")
    ax.set_ylim(0.5, 1.0)
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    savefig("fig2_snr_vs_kappa_t")


# ─── Figure 3: Objective comparison ──────────────────────────────────────────
def fig_objective_comparison(d1: dict) -> None:
    fig, ax = plt.subplots(figsize=(6, 4))
    kT = d1["obj_cmp_kappa_t_scan"]
    snr2 = d1["obj_cmp_snr2_arr"]
    endsep = d1["obj_cmp_endsep_arr"]
    snr2_ss = d1["obj_cmp_snr2_ss_approx"]
    eff = d1["obj_cmp_efficiency"]
    endsep_ss = float(d1["obj_cmp_endsep_ss"])

    ax2 = ax.twinx()
    l1, = ax.plot(kT, snr2,   color=TOL["blue"],  lw=2,  label=r"Integrated SNR$^2$ (ODE)")
    l2, = ax.plot(kT, snr2_ss, color=TOL["blue"], lw=1.5, ls=":", label=r"SS approx SNR$^2$")
    l3, = ax2.plot(kT, endsep, color=TOL["red"],   lw=2,  label=r"Endpoint $|\Delta\alpha(T)|^2$")
    ax2.axhline(endsep_ss, color=TOL["red"], lw=1, ls="--")

    ax.set_xlabel(r"$\kappa T$")
    ax.set_ylabel(r"$\mathrm{SNR}^2$", color=TOL["blue"])
    ax2.set_ylabel(r"Endpoint sep $|\Delta\alpha(T)|^2$", color=TOL["red"])
    ax.tick_params(axis="y", labelcolor=TOL["blue"])
    ax2.tick_params(axis="y", labelcolor=TOL["red"])
    ax.set_title("Endpoint separation vs integrated SNR² as objectives")
    ax.legend(handles=[l1, l2, l3], loc="upper left", fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 20)
    plt.tight_layout()
    savefig("fig3_objective_comparison")


# ─── Figure 4: Pulse family comparison ───────────────────────────────────────
def fig_pulse_families(d2: dict) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    chi_kappa = d2["fam_chi_kappa_vals"]
    families = list(d2.get("fam_families", ["square", "gaussian", "hann", "cosine_rise", "spline"]))
    snr2 = d2["fam_snr2"]
    F    = d2["fam_F_assign"]

    ax = axes[0]
    for i, fam in enumerate(families):
        ax.plot(chi_kappa, snr2[i], marker=FAMILY_MARKERS.get(fam, "o"),
                color=FAMILY_COLORS.get(fam, "black"),
                label=FAMILY_LABELS.get(fam, fam), lw=1.8, ms=6)
    ax.set_xlabel(r"$|\chi|/\kappa$")
    ax.set_ylabel(r"$\mathrm{SNR}^2$")
    ax.set_title(r"(a) SNR² vs $\chi/\kappa$  ($\kappa T = 5$)")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    dur_kT = d2["dur_kappa_t_vals"]
    dur_snr2 = d2["dur_snr2"]
    for i, fam in enumerate(families):
        ax.plot(dur_kT, dur_snr2[i], color=FAMILY_COLORS.get(fam, "black"),
                label=FAMILY_LABELS.get(fam, fam), lw=1.8)
    ax.set_xlabel(r"$\kappa T$")
    ax.set_ylabel(r"$\mathrm{SNR}^2$")
    ax.set_title(r"(b) SNR² vs $\kappa T$ at $\chi/\kappa = 1$")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    savefig("fig4_pulse_families")


# ─── Figure 5: Depletion study ────────────────────────────────────────────────
def fig_depletion(d2: dict) -> None:
    fig, ax = plt.subplots(figsize=(6, 4))
    kT_dep = d2["dep_kT_dep"]
    n_nat  = d2["dep_n_res_natural"]
    n_act  = d2["dep_n_res_active"]
    n_read = float(d2["dep_n_res_read_only"])

    ax.axhline(n_read, color=TOL["grey"], lw=2, ls=":", label="Read-only (no depletion)")
    ax.plot(kT_dep, n_nat, color=TOL["red"],  lw=2, label="Natural ring-down")
    ax.plot(kT_dep, n_act, color=TOL["blue"], lw=2, label="Active depletion")
    ax.set_xlabel(r"Additional depletion time $\kappa T_\mathrm{dep}$")
    ax.set_ylabel(r"Residual photons $n_\mathrm{res}$")
    ax.set_title("Readout depletion: natural vs. active")
    ax.legend()
    ax.set_yscale("log")
    ax.grid(True, alpha=0.3, which="both")
    plt.tight_layout()
    savefig("fig5_depletion")


# ─── Figure 6: GRAPE vs families ─────────────────────────────────────────────
def fig_grape_comparison(d3: dict) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    families = list(d3.get("families", ["square", "gaussian", "hann", "cosine_rise", "spline"]))

    # Panel (a): SNR² vs χ/κ
    ax = axes[0]
    chi_kappa = d3["grape_chi_kappa_vals"]
    snr2_grape   = d3["grape_snr2_grape"]
    snr2_families = d3["grape_snr2_family"]

    for i, fam in enumerate(families):
        ax.plot(chi_kappa, snr2_families[i], marker=FAMILY_MARKERS.get(fam, "o"),
                color=FAMILY_COLORS.get(fam, "grey"), lw=1.5, ms=5, alpha=0.7,
                label=FAMILY_LABELS.get(fam, fam))
    ax.plot(chi_kappa, snr2_grape, color=TOL["black"], marker="*", ms=10, lw=2.5,
            label="GRAPE", zorder=5)
    ax.set_xlabel(r"$|\chi|/\kappa$")
    ax.set_ylabel(r"$\mathrm{SNR}^2$")
    ax.set_title(r"(a) GRAPE vs. pulse families ($\kappa T = 5$)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Panel (b): GRAPE vs κT
    ax = axes[1]
    kT_vals = d3["kT_vals"]
    snr2_g = d3["kT_snr2_grape"]
    snr2_sq = d3["kT_snr2_square"]
    snr2_h  = d3["kT_snr2_hann"]
    ax.plot(kT_vals, snr2_g,  color=TOL["black"], lw=2.5, marker="*", ms=9, label="GRAPE")
    ax.plot(kT_vals, snr2_sq, color=TOL["blue"],  lw=2,   marker="o", ms=7, label="Square")
    ax.plot(kT_vals, snr2_h,  color=TOL["red"],   lw=2,   marker="^", ms=7, label="Hann")
    ax.set_xlabel(r"$\kappa T$")
    ax.set_ylabel(r"$\mathrm{SNR}^2$")
    ax.set_title(r"(b) GRAPE vs. baselines across $\kappa T$ ($\chi/\kappa = 1$)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    savefig("fig6_grape_comparison")


# ─── Figure 7: Net fidelity with open-system effects ─────────────────────────
def fig_open_system(d4: dict) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    # Panel (a): F_net vs κT
    ax = axes[0]
    kT = d4["opt_kT_vals"]
    F_net = d4["opt_F_net"]
    snr2  = d4["opt_snr2"]
    eta   = d4["opt_eta_phi"]

    ax.plot(kT, F_net, color=TOL["blue"], lw=2.5, label=r"$F_\mathrm{net}$ (with $T_1$)")
    ax2 = ax.twinx()
    ax2.plot(kT, eta, color=TOL["red"], lw=2, ls="--", label=r"$\eta_\varphi$ (dephasing)")
    ax.set_xlabel(r"$\kappa T$")
    ax.set_ylabel(r"Net assignment fidelity $F_\mathrm{net}$", color=TOL["blue"])
    ax2.set_ylabel(r"Dephasing factor $\eta_\varphi$", color=TOL["red"])
    ax.tick_params(axis="y", labelcolor=TOL["blue"])
    ax2.tick_params(axis="y", labelcolor=TOL["red"])
    ax.set_title(r"(a) Net readout fidelity vs $\kappa T$ with $T_1$ and dephasing")
    kT_opt = float(d4["opt_kT_optimal"])
    ax.axvline(kT_opt, ls=":", color=TOL["green"], lw=1.5,
               label=f"Optimal $\\kappa T = {kT_opt:.1f}$")
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(True, alpha=0.3)

    # Panel (b): SNR² and dephasing vs drive amplitude
    ax = axes[1]
    eps_kappa = d4["damp_eps_kappa"]
    snr2_amp  = d4["damp_snr2"]
    eta_amp   = d4["damp_eta_phi"]
    n_ss      = d4["damp_n_ss"]
    ax.plot(eps_kappa, snr2_amp, color=TOL["blue"], lw=2, label=r"$\mathrm{SNR}^2$")
    ax2 = ax.twinx()
    ax2.plot(eps_kappa, eta_amp, color=TOL["red"], lw=2, ls="--", label=r"$\eta_\varphi$")
    ax.set_xlabel(r"$|\varepsilon|/\kappa$ (drive amplitude)")
    ax.set_ylabel(r"$\mathrm{SNR}^2$", color=TOL["blue"])
    ax2.set_ylabel(r"Dephasing factor $\eta_\varphi$", color=TOL["red"])
    ax.tick_params(axis="y", labelcolor=TOL["blue"])
    ax2.tick_params(axis="y", labelcolor=TOL["red"])
    ax.set_title(r"(b) SNR$^2$ and dephasing vs drive amplitude")
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    savefig("fig7_open_system")


# ─── Figure 8: Robustness ─────────────────────────────────────────────────────
def fig_robustness(d5: dict) -> None:
    labels_ordered = ["square", "hann", "GRAPE_nominal", "GRAPE_robust"]
    nice_labels = {"square": "Square", "hann": "Hann",
                   "GRAPE_nominal": "Nom. GRAPE", "GRAPE_robust": "Rob. GRAPE"}
    colors = [TOL["blue"], TOL["red"], TOL["green"], TOL["purple"]]

    chi_err   = d5["chi_err_grid"]
    kappa_err = d5["kappa_err_grid"]

    fig = plt.figure(figsize=(14, 5))
    gs = gridspec.GridSpec(1, len(labels_ordered) + 1, width_ratios=[1]*len(labels_ordered) + [0.08])

    vmin = np.inf; vmax = -np.inf
    for lbl in labels_ordered:
        key = lbl + "_snr2_map"
        if key in d5:
            vmin = min(vmin, np.min(d5[key]))
            vmax = max(vmax, np.max(d5[key]))

    for k, lbl in enumerate(labels_ordered):
        ax = fig.add_subplot(gs[k])
        key = lbl + "_snr2_map"
        if key in d5:
            im = ax.pcolormesh(
                kappa_err * 100, chi_err * 100, d5[key],
                cmap="RdYlGn", vmin=vmin, vmax=vmax
            )
        nom_snr2 = float(d5.get(lbl + "_snr2_nom", 0))
        mean_snr2 = float(d5.get(lbl + "_mean_snr2", 0))
        cv = float(d5.get(lbl + "_cv", 0))
        ax.set_title(f"{nice_labels[lbl]}\nnominal={nom_snr2:.2f}, CV={cv:.3f}", fontsize=10)
        ax.set_xlabel(r"$\kappa$ error (%)")
        if k == 0:
            ax.set_ylabel(r"$\chi$ error (%)")
        else:
            ax.set_yticklabels([])
        ax.set_xticks([-20, -10, 0, 10, 20])
        ax.set_yticks([-20, -10, 0, 10, 20])

    cax = fig.add_subplot(gs[-1])
    sm = ScalarMappable(cmap="RdYlGn", norm=Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    plt.colorbar(sm, cax=cax, label=r"$\mathrm{SNR}^2$")

    plt.suptitle(r"Robustness: SNR$^2$ under $\pm 20\%$ parameter uncertainty", fontsize=12)
    plt.tight_layout()
    savefig("fig8_robustness")


# ─── Figure 9: GRAPE waveform shape ──────────────────────────────────────────
def fig_waveform_shapes(d3: dict) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    tlist = d3["wf_tlist"]
    t_ns = tlist * 1e9
    t_sq_ns = d3.get("wf_tlist_sq", tlist) * 1e9
    t_hann_ns = d3.get("wf_tlist_hann", tlist) * 1e9

    # Drive amplitude (GRAPE)
    ax = axes[0, 0]
    eps_grape = d3["wf_eps_grape"]
    t_seg = np.linspace(t_ns[0], t_ns[-1], len(eps_grape))
    ax.step(t_seg, np.real(eps_grape) / TWO_PI / 1e6, where="mid",
            color=TOL["blue"], lw=2, label="I (Re ε)")
    ax.step(t_seg, np.imag(eps_grape) / TWO_PI / 1e6, where="mid",
            color=TOL["red"], lw=2, ls="--", label="Q (Im ε)")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel(r"Drive $\varepsilon(t)/2\pi$ (MHz)")
    ax.set_title("(a) GRAPE waveform")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # Phase-space trajectories
    ax = axes[0, 1]
    ag_grape = d3["wf_ag_grape"]
    ae_grape = d3["wf_ae_grape"]
    ag_sq = d3["wf_ag_sq"]
    ae_sq = d3["wf_ae_sq"]
    ax.plot(np.real(ag_grape), np.imag(ag_grape), color=TOL["blue"],  lw=2, label="GRAPE |g⟩")
    ax.plot(np.real(ae_grape), np.imag(ae_grape), color=TOL["red"],   lw=2, label="GRAPE |e⟩")
    ax.plot(np.real(ag_sq),   np.imag(ag_sq),   color=TOL["blue"],  lw=1.5, ls="--", label="Square |g⟩")
    ax.plot(np.real(ae_sq),   np.imag(ae_sq),   color=TOL["red"],   lw=1.5, ls="--", label="Square |e⟩")
    ax.plot(0, 0, "ko", ms=8, zorder=5)
    ax.set_xlabel(r"Re $\alpha$")
    ax.set_ylabel(r"Im $\alpha$")
    ax.set_title("(b) Phase-space trajectories")
    ax.legend(fontsize=8)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)

    # Separation |Δα(t)|
    ax = axes[1, 0]
    delta_grape = np.abs(ae_grape - ag_grape)
    delta_sq    = np.abs(ae_sq - ag_sq)
    ag_hann = d3["wf_ag_hann"]
    ae_hann = d3["wf_ae_hann"]
    delta_hann = np.abs(ae_hann - ag_hann)
    ax.plot(t_ns, delta_grape, color=TOL["black"], lw=2, label=f"GRAPE (SNR²={d3['wf_snr2_grape']:.2f})")
    ax.plot(t_sq_ns, delta_sq,   color=TOL["blue"],  lw=2, ls="--", label=f"Square (SNR²={d3['wf_snr2_sq']:.2f})")
    ax.plot(t_hann_ns, delta_hann, color=TOL["red"],   lw=2, ls="-.", label=f"Hann (SNR²={d3['wf_snr2_hann']:.2f})")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel(r"$|\alpha_e(t) - \alpha_g(t)|$")
    ax.set_title("(c) Pointer separation vs time")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # Multi-objective Pareto
    ax = axes[1, 1]
    mo_snr2 = d3.get("mo_snr2", None)
    mo_nres = d3.get("mo_n_res", None)
    if mo_snr2 is not None and mo_nres is not None:
        ax.scatter(mo_nres, mo_snr2, c=range(len(mo_snr2)), cmap="viridis", s=80, zorder=5)
        ax.plot(mo_nres, mo_snr2, color=TOL["grey"], lw=1)
        ax.set_xlabel(r"Residual photons $n_\mathrm{res}$")
        ax.set_ylabel(r"$\mathrm{SNR}^2$")
        ax.set_title("(d) Pareto frontier: SNR² vs residual photons")
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    savefig("fig9_waveforms")


def load_results() -> dict:
    data = {}
    for phase in range(1, 6):
        fname = os.path.join(DATA_DIR, f"phase{phase}_results.npz")
        if os.path.exists(fname):
            arr = np.load(fname, allow_pickle=True)
            data[phase] = dict(arr)
            print(f"  Loaded phase{phase}_results.npz ({len(arr.files)} arrays)")
        else:
            print(f"  WARNING: {fname} not found — skipping phase {phase}")
            data[phase] = {}
    return data


def main() -> None:
    print("=" * 60)
    print("Generating publication-quality figures")
    print("=" * 60)

    data = load_results()
    d1 = data.get(1, {})
    d2 = data.get(2, {})
    d3 = data.get(3, {})
    d4 = data.get(4, {})
    d5 = data.get(5, {})

    if d1:
        print("\n[Fig 1] Steady-state separation...")
        try:
            fig_steady_state(d1)
        except Exception as e:
            print(f"  Error: {e}")

        print("\n[Fig 2] SNR² vs κT...")
        try:
            fig_snr_vs_kappa_t(d1)
        except Exception as e:
            print(f"  Error: {e}")

        print("\n[Fig 3] Objective comparison...")
        try:
            fig_objective_comparison(d1)
        except Exception as e:
            print(f"  Error: {e}")

    if d2:
        print("\n[Fig 4] Pulse family comparison...")
        try:
            fig_pulse_families(d2)
        except Exception as e:
            print(f"  Error: {e}")

        print("\n[Fig 5] Depletion study...")
        try:
            fig_depletion(d2)
        except Exception as e:
            print(f"  Error: {e}")

    if d3:
        print("\n[Fig 6] GRAPE vs families...")
        try:
            fig_grape_comparison(d3)
        except Exception as e:
            print(f"  Error: {e}")

        print("\n[Fig 9] Waveform shapes...")
        try:
            # Merge multi-objective data into d3 for this figure
            d3_ext = dict(d3)
            d3_ext["mo_snr2"] = d3.get("mo_snr2", None)
            d3_ext["mo_n_res"] = d3.get("mo_n_res", None)
            fig_waveform_shapes(d3_ext)
        except Exception as e:
            print(f"  Error: {e}")

    if d4:
        print("\n[Fig 7] Open-system effects...")
        try:
            fig_open_system(d4)
        except Exception as e:
            print(f"  Error: {e}")

    if d5:
        print("\n[Fig 8] Robustness sensitivity maps...")
        try:
            fig_robustness(d5)
        except Exception as e:
            print(f"  Error: {e}")

    print("\nAll figures saved to:", FIGURES_DIR)


if __name__ == "__main__":
    main()
