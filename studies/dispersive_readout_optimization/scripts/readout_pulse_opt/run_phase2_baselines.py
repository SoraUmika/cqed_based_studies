"""
Phase 2 — Baseline Pulse-Family Comparison
==========================================

Goals:
  1. Compare five pulse families (square, Gaussian, Hann, cosine-rise, spline)
     across χ/κ ∈ {0.25, 0.5, 1.0, 2.0, 4.0} at fixed κT = 5.
  2. Optimize drive amplitude independently for each family and each χ/κ.
  3. Compare against readout-only at midpoint vs on-resonance drive.
  4. Study readout-only vs readout+natural-ringdown vs readout+active-depletion,
     measuring residual photons and total cycle time penalty.
  5. Sweep κT at fixed χ/κ = 1 to find the optimal readout duration per family.

Outputs:
  data/phase2_results.npz
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

import numpy as np
from common import (
    TWO_PI, KAPPA, CHI_NOMINAL, EPSILON_MAX,
    CHI_KAPPA_VALUES, KAPPA_T_VALUES,
    DT_ODE, N_GRAPE_SEGMENTS,
    optimal_delta_g, simulate_conditioned_fields,
    snr_squared, endpoint_separation_sq, assignment_fidelity_from_snr2,
    residual_photons, measurement_dephasing_factor,
    evaluate_pulse_family, optimize_amplitude,
    PULSE_FAMILIES, build_square, build_depletion_pulse,
)

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
os.makedirs(DATA_DIR, exist_ok=True)

# Fixed κT for family comparison sweep
KAPPA_T_FIXED = 5.0
FAMILIES = list(PULSE_FAMILIES.keys())


def compare_pulse_families_vs_chi_kappa(
    kappa: float = KAPPA,
    kappa_t: float = KAPPA_T_FIXED,
    epsilon_max: float = EPSILON_MAX,
    dt: float = DT_ODE,
) -> dict:
    """
    For each χ/κ value and each pulse family, find the best amplitude and
    evaluate SNR², F_assign, endpoint separation, residual photons.
    """
    T_read = kappa_t / kappa
    chi_kappa_vals = CHI_KAPPA_VALUES
    chi_vals = chi_kappa_vals * kappa

    snr2       = np.zeros((len(FAMILIES), len(chi_kappa_vals)))
    F_assign   = np.zeros_like(snr2)
    endsep     = np.zeros_like(snr2)
    n_res      = np.zeros_like(snr2)
    dephase    = np.zeros_like(snr2)
    eps_opt    = np.zeros_like(snr2)

    for j, chi in enumerate(chi_vals):
        delta_g = optimal_delta_g(chi)
        print(f"  χ/κ = {chi_kappa_vals[j]:.2f} | Δ_g/2π = {delta_g/TWO_PI*1e-6:.2f} MHz")
        for i, family in enumerate(FAMILIES):
            best_eps, ev = optimize_amplitude(
                family, T_read, kappa, chi, delta_g, epsilon_max, dt
            )
            snr2[i, j]     = ev["snr2"]
            F_assign[i, j] = ev["F_assign"]
            endsep[i, j]   = ev["endpoint_sep2"]
            n_res[i, j]    = ev["n_res"]
            dephase[i, j]  = ev["dephasing_factor"]
            eps_opt[i, j]  = best_eps

    return {
        "families": FAMILIES,
        "chi_kappa_vals": chi_kappa_vals,
        "snr2": snr2,
        "F_assign": F_assign,
        "endsep": endsep,
        "n_res": n_res,
        "dephase": dephase,
        "eps_opt": eps_opt,
        "kappa_t": kappa_t,
    }


def compare_drive_frequencies(
    kappa: float = KAPPA,
    chi: float = None,
    kappa_t: float = KAPPA_T_FIXED,
    epsilon_max: float = EPSILON_MAX,
    dt: float = DT_ODE,
) -> dict:
    """
    Compare three drive strategies for the square pulse:
    (a) On resonance with |g⟩: Δ_g = 0
    (b) Midpoint (optimal): Δ_g = −χ/2
    (c) On resonance with |e⟩: Δ_g = −χ

    Shows that midpoint drive dominates for all κT.
    """
    if chi is None:
        chi = CHI_NOMINAL   # χ/κ = 1

    T_read = kappa_t / kappa
    strategies = {
        "on_g":     0.0,
        "midpoint": optimal_delta_g(chi),
        "on_e":     -chi,
    }
    kappa_t_scan = np.linspace(0.3, 20.0, 100)
    T_scan = kappa_t_scan / kappa

    snr2_results = {name: np.zeros(len(T_scan)) for name in strategies}
    F_assign_results = {name: np.zeros(len(T_scan)) for name in strategies}

    for k, (name, delta_g) in enumerate(strategies.items()):
        for i, T in enumerate(T_scan):
            N = max(2, int(round(T / dt)))
            tlist = np.linspace(0.0, T, N + 1)
            best_eps, ev = optimize_amplitude(
                "square", T, kappa, chi, delta_g, epsilon_max, dt, n_grid=40
            )
            snr2_results[name][i]    = ev["snr2"]
            F_assign_results[name][i] = ev["F_assign"]

    return {
        "kappa_t_scan": kappa_t_scan,
        "strategies": list(strategies.keys()),
        "delta_g_values": list(strategies.values()),
        "snr2": {k: v.tolist() for k, v in snr2_results.items()},
        "F_assign": {k: v.tolist() for k, v in F_assign_results.items()},
    }


def optimal_duration_per_family(
    kappa: float = KAPPA,
    chi: float = None,
    epsilon_max: float = EPSILON_MAX,
    dt: float = DT_ODE,
) -> dict:
    """
    For each pulse family, sweep κT ∈ KAPPA_T_VALUES and record SNR² and F_assign
    at χ/κ = 1 (nominal).  This identifies the saturation behavior per family.
    """
    if chi is None:
        chi = CHI_NOMINAL
    delta_g = optimal_delta_g(chi)
    kT_vals = np.linspace(0.5, 20.0, 60)
    T_vals = kT_vals / kappa

    snr2_mat   = np.zeros((len(FAMILIES), len(kT_vals)))
    F_assign_mat = np.zeros_like(snr2_mat)
    n_res_mat  = np.zeros_like(snr2_mat)

    for i, family in enumerate(FAMILIES):
        for j, T in enumerate(T_vals):
            best_eps, ev = optimize_amplitude(
                family, T, kappa, chi, delta_g, epsilon_max, dt, n_grid=40
            )
            snr2_mat[i, j]    = ev["snr2"]
            F_assign_mat[i, j] = ev["F_assign"]
            n_res_mat[i, j]   = ev["n_res"]

    return {
        "families": FAMILIES,
        "kappa_t_vals": kT_vals,
        "snr2": snr2_mat,
        "F_assign": F_assign_mat,
        "n_res": n_res_mat,
        "chi_kappa": abs(chi) / kappa,
    }


def study_depletion(
    kappa: float = KAPPA,
    chi: float = None,
    kappa_t_read: float = 3.0,
    epsilon_max: float = EPSILON_MAX,
    dt: float = DT_ODE,
) -> dict:
    """
    Compare three depletion strategies after a square readout pulse:
      (a) readout only — no depletion, photons remain.
      (b) natural ring-down — wait additional κT_dep = 3 with no drive.
      (c) active depletion — apply optimal single-segment depletion pulse.

    Metrics: residual photons n_res, total time κ(T_read+T_dep), SNR² from readout.
    """
    if chi is None:
        chi = CHI_NOMINAL
    delta_g = optimal_delta_g(chi)
    T_read = kappa_t_read / kappa

    # Optimal amplitude for square pulse
    best_eps, ev_read = optimize_amplitude(
        "square", T_read, kappa, chi, delta_g, epsilon_max, dt
    )

    # State of resonator at end of readout
    ag = ev_read["alpha_g"]
    ae = ev_read["alpha_e"]
    alpha_g_end = ag[-1]
    alpha_e_end = ae[-1]
    snr2_read = ev_read["snr2"]

    # ── (a) readout only ─────────────────────────────────────────────────────
    res_read_only = {
        "n_res": ev_read["n_res"],
        "kappa_T_total": kappa_t_read,
        "snr2": snr2_read,
    }

    # ── (b) natural ring-down ────────────────────────────────────────────────
    kT_dep_values = np.linspace(0.5, 10.0, 60)
    T_dep_values = kT_dep_values / kappa
    n_res_nat = np.zeros(len(T_dep_values))
    for k, T_dep in enumerate(T_dep_values):
        N_dep = max(2, int(round(T_dep / dt)))
        tlist_dep = np.linspace(0.0, T_dep, N_dep + 1)
        eps_zero = np.zeros(N_dep, dtype=np.complex128)
        from common import integrate_readout_ode
        ag_dep = integrate_readout_ode(eps_zero, tlist_dep, kappa, delta_g,      alpha_g_end)
        ae_dep = integrate_readout_ode(eps_zero, tlist_dep, kappa, delta_g + chi, alpha_e_end)
        n_res_nat[k] = residual_photons(ag_dep, ae_dep)

    # ── (c) active depletion ─────────────────────────────────────────────────
    n_res_active = np.zeros(len(T_dep_values))
    for k, T_dep in enumerate(T_dep_values):
        N_dep = max(2, int(round(T_dep / dt)))
        tlist_dep = np.linspace(0.0, T_dep, N_dep + 1)
        # Build depletion drive targeting α_g → 0
        eps_dep = build_depletion_pulse(
            tlist_dep, alpha_g_end, kappa, delta_g, delta_g + chi
        )
        from common import integrate_readout_ode
        ag_dep = integrate_readout_ode(eps_dep, tlist_dep, kappa, delta_g,       alpha_g_end)
        ae_dep = integrate_readout_ode(eps_dep, tlist_dep, kappa, delta_g + chi, alpha_e_end)
        n_res_active[k] = residual_photons(ag_dep, ae_dep)

    return {
        "kappa_t_read": kappa_t_read,
        "snr2_read": snr2_read,
        "alpha_g_end": alpha_g_end,
        "alpha_e_end": alpha_e_end,
        "n_res_read_only": ev_read["n_res"],
        "kT_dep_values": kT_dep_values,
        "n_res_natural": n_res_nat,
        "n_res_active": n_res_active,
        "F_assign": ev_read["F_assign"],
        "best_eps": best_eps,
    }


def waveform_examples(
    kappa: float = KAPPA,
    chi: float = None,
    kappa_t: float = 5.0,
    epsilon: float = None,
    dt: float = DT_ODE,
) -> dict:
    """
    Evaluate all pulse families at the same amplitude and store field trajectories
    for waveform shape comparison figure.
    """
    if chi is None:
        chi = CHI_NOMINAL
    if epsilon is None:
        epsilon = EPSILON_MAX * 0.6
    T_read = kappa_t / kappa
    delta_g = optimal_delta_g(chi)
    N = max(2, int(round(T_read / dt)))
    tlist = np.linspace(0.0, T_read, N + 1)
    out = {"tlist": tlist, "families": FAMILIES}
    for family in FAMILIES:
        ev = evaluate_pulse_family(family, epsilon, T_read, kappa, chi, delta_g, dt)
        out[f"alpha_g_{family}"] = ev["alpha_g"]
        out[f"alpha_e_{family}"] = ev["alpha_e"]
        out[f"epsilon_{family}"] = ev["epsilon_t"]
        out[f"snr2_{family}"] = ev["snr2"]
        out[f"n_res_{family}"] = ev["n_res"]
    return out


def main() -> None:
    print("=" * 60)
    print("Phase 2: Baseline Pulse-Family Comparison")
    print("=" * 60)

    kappa = KAPPA

    # ── 1. Pulse family vs χ/κ comparison ───────────────────────────────────
    print("\n[1/4] Comparing pulse families vs χ/κ at κT = 5...")
    family_data = compare_pulse_families_vs_chi_kappa(kappa)
    for i, fam in enumerate(FAMILIES):
        snr2_max = np.max(family_data["snr2"][i])
        best_chi_kappa = family_data["chi_kappa_vals"][np.argmax(family_data["snr2"][i])]
        print(f"  {fam:12s}:  max SNR²={snr2_max:.2f}  at χ/κ={best_chi_kappa:.2f}")

    # ── 2. Drive frequency comparison (square pulse) ─────────────────────────
    print("\n[2/4] Comparing drive strategies (square pulse, χ/κ=1)...")
    drive_data = compare_drive_frequencies(kappa, CHI_NOMINAL)
    idx5 = np.argmin(np.abs(np.array(drive_data["kappa_t_scan"]) - 5.0))
    for name in drive_data["strategies"]:
        print(f"  {name:10s}: SNR²(κT=5)={drive_data['snr2'][name][idx5]:.3f}, "
              f"F_assign={drive_data['F_assign'][name][idx5]:.4f}")

    # ── 3. Optimal duration per family ──────────────────────────────────────
    print("\n[3/4] Finding optimal κT per family at χ/κ=1...")
    duration_data = optimal_duration_per_family(kappa)
    for i, fam in enumerate(FAMILIES):
        kT_best = duration_data["kappa_t_vals"][np.argmax(duration_data["snr2"][i])]
        snr2_best = np.max(duration_data["snr2"][i])
        print(f"  {fam:12s}: best κT = {kT_best:.1f}, SNR² = {snr2_best:.2f}")

    # ── 4. Depletion study ───────────────────────────────────────────────────
    print("\n[4/4] Studying readout+depletion strategies...")
    dep_data = study_depletion(kappa)
    idx_1 = np.argmin(np.abs(dep_data["kT_dep_values"] - 1.0))
    idx_3 = np.argmin(np.abs(dep_data["kT_dep_values"] - 3.0))
    print(f"  SNR² (readout only): {dep_data['snr2_read']:.3f}")
    print(f"  n_res (readout only): {dep_data['n_res_read_only']:.4f}")
    print(f"  n_res after 1/κ natural decay:    {dep_data['n_res_natural'][idx_1]:.4f}")
    print(f"  n_res after 1/κ active depletion: {dep_data['n_res_active'][idx_1]:.4f}")
    print(f"  n_res after 3/κ natural decay:    {dep_data['n_res_natural'][idx_3]:.4f}")
    print(f"  n_res after 3/κ active depletion: {dep_data['n_res_active'][idx_3]:.4f}")

    # ── Waveform examples ────────────────────────────────────────────────────
    print("\n[extra] Generating waveform trajectory examples...")
    wf_data = waveform_examples(kappa)

    # ── Save ─────────────────────────────────────────────────────────────────
    outpath = os.path.join(DATA_DIR, "phase2_results.npz")

    # Flatten dicts for npz storage
    save_dict = {}
    for k, v in family_data.items():
        if isinstance(v, (list, np.ndarray)):
            save_dict[f"fam_{k}"] = np.asarray(v) if not isinstance(v[0], str) else np.array(v, dtype=object)
        else:
            save_dict[f"fam_{k}"] = v

    save_dict.update({
        "drive_kappa_t_scan": drive_data["kappa_t_scan"],
        "drive_snr2_on_g":    drive_data["snr2"]["on_g"],
        "drive_snr2_midpt":   drive_data["snr2"]["midpoint"],
        "drive_snr2_on_e":    drive_data["snr2"]["on_e"],
        "drive_F_on_g":       drive_data["F_assign"]["on_g"],
        "drive_F_midpt":      drive_data["F_assign"]["midpoint"],
        "drive_F_on_e":       drive_data["F_assign"]["on_e"],
    })
    save_dict.update({
        "dur_kappa_t_vals": duration_data["kappa_t_vals"],
        "dur_snr2":         duration_data["snr2"],
        "dur_F_assign":     duration_data["F_assign"],
        "dur_n_res":        duration_data["n_res"],
    })
    save_dict.update({
        "dep_kT_dep": dep_data["kT_dep_values"],
        "dep_n_res_read_only": dep_data["n_res_read_only"],
        "dep_n_res_natural":   dep_data["n_res_natural"],
        "dep_n_res_active":    dep_data["n_res_active"],
        "dep_snr2_read":       dep_data["snr2_read"],
        "dep_kappa_t_read":    dep_data["kappa_t_read"],
    })
    save_dict.update({
        "wf_tlist":  wf_data["tlist"],
        **{f"wf_ag_{fam}":  wf_data[f"alpha_g_{fam}"] for fam in FAMILIES},
        **{f"wf_ae_{fam}":  wf_data[f"alpha_e_{fam}"] for fam in FAMILIES},
        **{f"wf_eps_{fam}": wf_data[f"epsilon_{fam}"]  for fam in FAMILIES},
    })
    save_dict["families"] = np.array(FAMILIES, dtype=object)

    np.savez(outpath, **save_dict)
    print(f"\nResults saved to {outpath}")

    print("\n" + "=" * 60)
    print("Phase 2 Summary:")
    # Best family at χ/κ = 1
    idx_ck1 = np.argmin(np.abs(family_data["chi_kappa_vals"] - 1.0))
    snr2_at_ck1 = family_data["snr2"][:, idx_ck1]
    best_idx = np.argmax(snr2_at_ck1)
    print(f"  Best family at χ/κ=1: {FAMILIES[best_idx]} (SNR²={snr2_at_ck1[best_idx]:.3f})")
    print(f"  Square pulse SNR²: {snr2_at_ck1[0]:.3f}")
    print(f"  Improvement over square: {snr2_at_ck1[best_idx]/snr2_at_ck1[0]:.3f}x")
    print("=" * 60)


if __name__ == "__main__":
    main()
