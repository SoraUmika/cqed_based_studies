"""
Phase 3 — Piecewise-Constant and GRAPE-Like Control
====================================================

Goals:
  1. Run GRAPE optimization for N_seg ∈ {20, 40, 60, 100} segments and assess
     convergence of SNR² vs. number of segments.
  2. Compare GRAPE vs. best pulse-family baseline at χ/κ ∈ {0.25, 0.5, 1.0, 2.0, 4.0}.
  3. Study GRAPE with multi-objective (SNR² + residual photon penalty): joint readout
     + depletion in a single N_seg-segment pulse.
  4. Show GRAPE waveform shapes at χ/κ = 1 for single-objective and multi-objective.
  5. Determine at which κT and χ/κ regime GRAPE clearly outperforms pulse families.

Outputs:
  data/phase3_results.npz
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

import numpy as np
from common import (
    TWO_PI, KAPPA, CHI_NOMINAL, EPSILON_MAX,
    CHI_KAPPA_VALUES, KAPPA_T_VALUES,
    N_GRAPE_SEGMENTS, N_GRAPE_RESTARTS, GRAPE_MAXITER, DT_ODE,
    optimal_delta_g, run_grape, optimize_amplitude,
    assignment_fidelity_from_snr2, residual_photons, snr_squared,
    simulate_conditioned_fields, evaluate_pulse_family,
    PULSE_FAMILIES,
)

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
os.makedirs(DATA_DIR, exist_ok=True)

KAPPA_T_MAIN = 5.0     # main κT for GRAPE comparisons
FAMILIES = list(PULSE_FAMILIES.keys())


def grape_vs_families_chi_kappa(
    kappa: float = KAPPA,
    kappa_t: float = KAPPA_T_MAIN,
    epsilon_max: float = EPSILON_MAX,
    n_seg: int = N_GRAPE_SEGMENTS,
    n_restarts: int = N_GRAPE_RESTARTS,
) -> dict:
    """
    For each χ/κ, run GRAPE and the best pulse-family baseline.
    Returns SNR², F_assign, n_res for GRAPE and each family.
    """
    chi_kappa_vals = CHI_KAPPA_VALUES
    T_read = kappa_t / kappa

    # GRAPE results
    snr2_grape  = np.zeros(len(chi_kappa_vals))
    F_grape     = np.zeros(len(chi_kappa_vals))
    n_res_grape = np.zeros(len(chi_kappa_vals))

    # Family results (best family at each χ/κ)
    snr2_family  = np.zeros((len(FAMILIES), len(chi_kappa_vals)))
    F_family     = np.zeros_like(snr2_family)
    n_res_family = np.zeros_like(snr2_family)

    for j, chi_kappa in enumerate(chi_kappa_vals):
        chi = chi_kappa * kappa
        delta_g = optimal_delta_g(chi)
        print(f"\n  χ/κ = {chi_kappa:.2f}")

        # GRAPE
        print(f"    Running GRAPE (N_seg={n_seg}, {n_restarts} restarts)...")
        gr = run_grape(
            kappa=kappa, chi=chi, T_read=T_read, N_seg=n_seg,
            delta_g=delta_g, epsilon_max=epsilon_max,
            n_restarts=n_restarts, maxiter=GRAPE_MAXITER,
            w_snr=1.0, w_res=0.0,
        )
        snr2_grape[j]  = gr["snr2_opt"]
        F_grape[j]     = gr["F_assign_opt"]
        n_res_grape[j] = gr["n_res_opt"]
        print(f"    GRAPE: SNR²={gr['snr2_opt']:.3f}, F={gr['F_assign_opt']:.4f}")

        # Pulse families
        for i, fam in enumerate(FAMILIES):
            best_eps, ev = optimize_amplitude(fam, T_read, kappa, chi, delta_g, epsilon_max)
            snr2_family[i, j]  = ev["snr2"]
            F_family[i, j]     = ev["F_assign"]
            n_res_family[i, j] = ev["n_res"]

    return {
        "chi_kappa_vals": chi_kappa_vals,
        "snr2_grape": snr2_grape,
        "F_grape": F_grape,
        "n_res_grape": n_res_grape,
        "snr2_family": snr2_family,
        "F_family": F_family,
        "n_res_family": n_res_family,
        "families": FAMILIES,
        "kappa_t": kappa_t,
        "n_seg": n_seg,
    }


def grape_segment_convergence(
    kappa: float = KAPPA,
    chi: float = None,
    kappa_t: float = KAPPA_T_MAIN,
    epsilon_max: float = EPSILON_MAX,
    n_restarts: int = N_GRAPE_RESTARTS,
) -> dict:
    """
    Study convergence of GRAPE SNR² as a function of the number of segments N_seg.
    """
    if chi is None:
        chi = CHI_NOMINAL
    delta_g = optimal_delta_g(chi)
    T_read = kappa_t / kappa
    seg_values = [10, 20, 30, 40, 60, 80, 100, 150]

    snr2_vals = []
    F_vals = []
    n_res_vals = []

    for n_seg in seg_values:
        print(f"  N_seg = {n_seg}...")
        gr = run_grape(
            kappa=kappa, chi=chi, T_read=T_read, N_seg=n_seg,
            delta_g=delta_g, epsilon_max=epsilon_max,
            n_restarts=n_restarts, maxiter=GRAPE_MAXITER,
            w_snr=1.0, w_res=0.0,
        )
        snr2_vals.append(gr["snr2_opt"])
        F_vals.append(gr["F_assign_opt"])
        n_res_vals.append(gr["n_res_opt"])
        print(f"    SNR²={gr['snr2_opt']:.4f}, F={gr['F_assign_opt']:.5f}")

    return {
        "seg_values": np.array(seg_values),
        "snr2": np.array(snr2_vals),
        "F_assign": np.array(F_vals),
        "n_res": np.array(n_res_vals),
        "chi_kappa": abs(chi) / kappa,
        "kappa_t": kappa_t,
    }


def grape_kappa_t_sweep(
    kappa: float = KAPPA,
    chi: float = None,
    epsilon_max: float = EPSILON_MAX,
    n_seg: int = N_GRAPE_SEGMENTS,
    n_restarts: int = N_GRAPE_RESTARTS,
) -> dict:
    """
    Sweep κT at fixed χ/κ = 1, compare GRAPE vs. square and Hann baselines.
    """
    if chi is None:
        chi = CHI_NOMINAL
    delta_g = optimal_delta_g(chi)
    kT_vals = np.array([1.0, 2.0, 3.0, 5.0, 8.0, 12.0])
    T_vals = kT_vals / kappa

    snr2_grape = np.zeros(len(kT_vals))
    F_grape    = np.zeros(len(kT_vals))
    snr2_sq    = np.zeros(len(kT_vals))
    snr2_hann  = np.zeros(len(kT_vals))

    for i, (kT, T) in enumerate(zip(kT_vals, T_vals)):
        print(f"  κT = {kT:.1f}...")
        gr = run_grape(
            kappa=kappa, chi=chi, T_read=T, N_seg=n_seg,
            delta_g=delta_g, epsilon_max=epsilon_max,
            n_restarts=n_restarts, maxiter=GRAPE_MAXITER,
        )
        snr2_grape[i] = gr["snr2_opt"]
        F_grape[i]    = gr["F_assign_opt"]

        _, ev_sq   = optimize_amplitude("square", T, kappa, chi, delta_g, epsilon_max)
        _, ev_hann = optimize_amplitude("hann",   T, kappa, chi, delta_g, epsilon_max)
        snr2_sq[i]   = ev_sq["snr2"]
        snr2_hann[i] = ev_hann["snr2"]
        print(f"    GRAPE={snr2_grape[i]:.3f}, Square={snr2_sq[i]:.3f}, Hann={snr2_hann[i]:.3f}")

    return {
        "kappa_t_vals": kT_vals,
        "snr2_grape": snr2_grape,
        "F_grape": F_grape,
        "snr2_square": snr2_sq,
        "snr2_hann": snr2_hann,
        "chi_kappa": abs(chi) / kappa,
    }


def grape_multi_objective(
    kappa: float = KAPPA,
    chi: float = None,
    kappa_t: float = KAPPA_T_MAIN,
    epsilon_max: float = EPSILON_MAX,
    n_seg: int = N_GRAPE_SEGMENTS,
    n_restarts: int = N_GRAPE_RESTARTS,
) -> dict:
    """
    Multi-objective GRAPE: optimize SNR² and minimize residual photons jointly.

    Study the Pareto tradeoff between SNR² and n_res by sweeping w_res.
    """
    if chi is None:
        chi = CHI_NOMINAL
    delta_g = optimal_delta_g(chi)
    T_read = kappa_t / kappa

    # Single-objective GRAPE (baseline for normalization)
    gr_snr = run_grape(kappa=kappa, chi=chi, T_read=T_read, N_seg=n_seg,
                       delta_g=delta_g, epsilon_max=epsilon_max,
                       n_restarts=n_restarts, w_snr=1.0, w_res=0.0)
    snr2_scale = gr_snr["snr2_opt"]    # use for normalization in multi-obj

    # Sweep weight w_res (penalty on residual photons)
    w_res_values = np.array([0.0, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0])
    snr2_pareto  = np.zeros(len(w_res_values))
    n_res_pareto = np.zeros(len(w_res_values))
    F_pareto     = np.zeros(len(w_res_values))

    for k, w_res in enumerate(w_res_values):
        print(f"  w_res = {w_res:.2f}...")
        gr = run_grape(kappa=kappa, chi=chi, T_read=T_read, N_seg=n_seg,
                       delta_g=delta_g, epsilon_max=epsilon_max,
                       n_restarts=n_restarts, w_snr=1.0, w_res=w_res,
                       n_res_scale=max(snr2_scale, 1.0))
        snr2_pareto[k]  = gr["snr2_opt"]
        n_res_pareto[k] = gr["n_res_opt"]
        F_pareto[k]     = gr["F_assign_opt"]
        print(f"    SNR²={snr2_pareto[k]:.3f}, n_res={n_res_pareto[k]:.4f}")

    # Store the full optimized waveform for w_res=0 and w_res=0.5
    gr_snr_only = gr_snr
    gr_balanced = run_grape(kappa=kappa, chi=chi, T_read=T_read, N_seg=n_seg,
                            delta_g=delta_g, epsilon_max=epsilon_max,
                            n_restarts=n_restarts, w_snr=1.0, w_res=0.5,
                            n_res_scale=max(snr2_scale, 1.0))

    return {
        "w_res_values": w_res_values,
        "snr2_pareto": snr2_pareto,
        "n_res_pareto": n_res_pareto,
        "F_pareto": F_pareto,
        "snr2_nominal": snr2_scale,
        # Waveforms
        "eps_snr_only": gr_snr_only["epsilon_opt"],
        "alpha_g_snr_only": gr_snr_only["alpha_g"],
        "alpha_e_snr_only": gr_snr_only["alpha_e"],
        "tlist_snr_only": gr_snr_only["tlist"],
        "eps_balanced": gr_balanced["epsilon_opt"],
        "alpha_g_balanced": gr_balanced["alpha_g"],
        "alpha_e_balanced": gr_balanced["alpha_e"],
        "snr2_balanced": gr_balanced["snr2_opt"],
        "n_res_balanced": gr_balanced["n_res_opt"],
        "n_seg": n_seg,
        "kappa_t": kappa_t,
    }


def grape_waveform_detail(
    kappa: float = KAPPA,
    chi: float = None,
    kappa_t: float = KAPPA_T_MAIN,
    epsilon_max: float = EPSILON_MAX,
    n_seg: int = N_GRAPE_SEGMENTS,
    n_restarts: int = N_GRAPE_RESTARTS,
) -> dict:
    """
    Return detailed waveform data for χ/κ=1, single-objective GRAPE.
    Used for the waveform-shape figure.
    """
    if chi is None:
        chi = CHI_NOMINAL
    delta_g = optimal_delta_g(chi)
    T_read = kappa_t / kappa
    gr = run_grape(kappa=kappa, chi=chi, T_read=T_read, N_seg=n_seg,
                   delta_g=delta_g, epsilon_max=epsilon_max,
                   n_restarts=n_restarts, w_snr=1.0, w_res=0.0)

    # Compare against square and Hann at the same T
    _, ev_sq   = optimize_amplitude("square", T_read, kappa, chi, delta_g, epsilon_max)
    _, ev_hann = optimize_amplitude("hann",   T_read, kappa, chi, delta_g, epsilon_max)

    return {
        "tlist": gr["tlist"],
        "eps_grape": gr["epsilon_opt"],
        "alpha_g_grape": gr["alpha_g"],
        "alpha_e_grape": gr["alpha_e"],
        "snr2_grape": gr["snr2_opt"],
        "n_res_grape": gr["n_res_opt"],
        # Square baseline
        "tlist_sq": ev_sq["tlist"],
        "alpha_g_sq": ev_sq["alpha_g"],
        "alpha_e_sq": ev_sq["alpha_e"],
        "snr2_sq": ev_sq["snr2"],
        # Hann baseline
        "tlist_hann": ev_hann["tlist"],
        "alpha_g_hann": ev_hann["alpha_g"],
        "alpha_e_hann": ev_hann["alpha_e"],
        "snr2_hann": ev_hann["snr2"],
        "n_seg": n_seg,
        "kappa_t": kappa_t,
    }


def main() -> None:
    print("=" * 60)
    print("Phase 3: GRAPE-Like Optimization")
    print("=" * 60)

    kappa = KAPPA

    # ── 1. GRAPE vs families across χ/κ ─────────────────────────────────────
    print("\n[1/5] GRAPE vs pulse families across χ/κ...")
    comp_data = grape_vs_families_chi_kappa(kappa)

    for j, ck in enumerate(comp_data["chi_kappa_vals"]):
        best_fam_idx = np.argmax(comp_data["snr2_family"][:, j])
        best_fam_snr2 = comp_data["snr2_family"][best_fam_idx, j]
        gain = (comp_data["snr2_grape"][j] - best_fam_snr2) / max(best_fam_snr2, 1e-30)
        print(f"  χ/κ={ck:.2f}: GRAPE SNR²={comp_data['snr2_grape'][j]:.3f}, "
              f"best-fam ({FAMILIES[best_fam_idx]})={best_fam_snr2:.3f}, "
              f"GRAPE gain={100*gain:.1f}%")

    # ── 2. Segment convergence ───────────────────────────────────────────────
    print("\n[2/5] GRAPE segment convergence at χ/κ=1...")
    seg_data = grape_segment_convergence(kappa)

    # ── 3. GRAPE κT sweep ────────────────────────────────────────────────────
    print("\n[3/5] GRAPE vs square/Hann across κT at χ/κ=1...")
    kT_data = grape_kappa_t_sweep(kappa)
    for k, kT in enumerate(kT_data["kappa_t_vals"]):
        gain_sq = (kT_data["snr2_grape"][k] - kT_data["snr2_square"][k]) / max(kT_data["snr2_square"][k], 1e-30)
        print(f"  κT={kT:.1f}: GRAPE gain over square = {100*gain_sq:.1f}%")

    # ── 4. Multi-objective GRAPE ─────────────────────────────────────────────
    print("\n[4/5] Multi-objective GRAPE: SNR² vs n_res tradeoff...")
    multi_data = grape_multi_objective(kappa)
    print(f"  SNR² (w_res=0):    {multi_data['snr2_pareto'][0]:.3f}")
    print(f"  n_res (w_res=0):   {multi_data['n_res_pareto'][0]:.4f}")
    print(f"  SNR² (w_res=0.5):  {multi_data['snr2_balanced']:.3f}")
    print(f"  n_res (w_res=0.5): {multi_data['n_res_balanced']:.4f}")

    # ── 5. Waveform detail ───────────────────────────────────────────────────
    print("\n[5/5] Computing waveform detail for figure...")
    wf_data = grape_waveform_detail(kappa)
    print(f"  GRAPE SNR²={wf_data['snr2_grape']:.3f}, Square={wf_data['snr2_sq']:.3f}, "
          f"Hann={wf_data['snr2_hann']:.3f}")

    # ── Save ─────────────────────────────────────────────────────────────────
    outpath = os.path.join(DATA_DIR, "phase3_results.npz")
    np.savez(
        outpath,
        # Family comparison
        grape_chi_kappa_vals=comp_data["chi_kappa_vals"],
        grape_snr2_grape=comp_data["snr2_grape"],
        grape_F_grape=comp_data["F_grape"],
        grape_n_res_grape=comp_data["n_res_grape"],
        grape_snr2_family=comp_data["snr2_family"],
        grape_F_family=comp_data["F_family"],
        grape_n_res_family=comp_data["n_res_family"],
        # Segment convergence
        seg_values=seg_data["seg_values"],
        seg_snr2=seg_data["snr2"],
        seg_F_assign=seg_data["F_assign"],
        # κT sweep
        kT_vals=kT_data["kappa_t_vals"],
        kT_snr2_grape=kT_data["snr2_grape"],
        kT_snr2_square=kT_data["snr2_square"],
        kT_snr2_hann=kT_data["snr2_hann"],
        # Multi-objective
        mo_w_res=multi_data["w_res_values"],
        mo_snr2=multi_data["snr2_pareto"],
        mo_n_res=multi_data["n_res_pareto"],
        mo_F=multi_data["F_pareto"],
        mo_eps_snr_only=multi_data["eps_snr_only"],
        mo_eps_balanced=multi_data["eps_balanced"],
        mo_alpha_g_snr_only=multi_data["alpha_g_snr_only"],
        mo_alpha_e_snr_only=multi_data["alpha_e_snr_only"],
        mo_alpha_g_balanced=multi_data["alpha_g_balanced"],
        mo_alpha_e_balanced=multi_data["alpha_e_balanced"],
        mo_tlist=multi_data["tlist_snr_only"],
        # Waveform detail
        wf_tlist=wf_data["tlist"],
        wf_eps_grape=wf_data["eps_grape"],
        wf_ag_grape=wf_data["alpha_g_grape"],
        wf_ae_grape=wf_data["alpha_e_grape"],
        wf_snr2_grape=wf_data["snr2_grape"],
        wf_tlist_sq=wf_data["tlist_sq"],
        wf_ag_sq=wf_data["alpha_g_sq"],
        wf_ae_sq=wf_data["alpha_e_sq"],
        wf_snr2_sq=wf_data["snr2_sq"],
        wf_tlist_hann=wf_data["tlist_hann"],
        wf_ag_hann=wf_data["alpha_g_hann"],
        wf_ae_hann=wf_data["alpha_e_hann"],
        wf_snr2_hann=wf_data["snr2_hann"],
        families=np.array(FAMILIES, dtype=object),
    )
    print(f"\nResults saved to {outpath}")

    print("\n" + "=" * 60)
    print("Phase 3 Summary:")
    best_j = np.argmax(comp_data["snr2_grape"])
    print(f"  GRAPE best at χ/κ={comp_data['chi_kappa_vals'][best_j]:.2f}: "
          f"SNR²={comp_data['snr2_grape'][best_j]:.3f}")
    print(f"  Segment convergence: {seg_data['snr2'][-1]:.4f} vs {seg_data['snr2'][-2]:.4f} "
          f"(last two N_seg: {seg_data['seg_values'][-1]}, {seg_data['seg_values'][-2]})")
    print("=" * 60)


if __name__ == "__main__":
    main()
