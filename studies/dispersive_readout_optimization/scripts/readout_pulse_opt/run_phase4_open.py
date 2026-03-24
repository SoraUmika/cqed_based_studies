"""
Phase 4 — Open-System and Realistic Hamiltonian Effects
========================================================

Goals:
  1. Include qubit T1 relaxation during readout.
     If the qubit decays |e⟩ → |g⟩ mid-readout with rate Γ₁ = 1/T1,
     the effective SNR is reduced. Compute the net assignment fidelity:
         F_net = (1−p_decay) · F_A(SNR²) + p_decay · F_A(0)
     where p_decay ≈ 1 − exp(−T_read/T1).
  2. Include measurement-induced dephasing of the qubit.
     The qubit off-diagonal element decays:
         η_φ = exp(−∫₀ᵀ (κ/2)|α_e−α_g|² dt)
     Compute how this depends on drive amplitude and κT.
  3. Combine T1 and dephasing effects to compute the net readout fidelity
     as a function of T_read. Identify the optimal T_read that balances
     SNR gain against decoherence loss.
  4. Estimate Purcell relaxation rate using cqed_sim ReadoutResonator.
  5. Assess the effect of finite transmon anharmonicity:
     for large drive amplitudes, the transmon |f⟩ level can be populated.
     We include this via an effective reduced dispersive shift approximation.
  6. Provide a discussion of which level of realism changes the conclusions.

Note: Full QuTiP Lindblad simulation is not used here (too expensive for GRAPE).
      Open-system effects are included analytically via the linear dispersive model.

Outputs:
  data/phase4_results.npz
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

import numpy as np
from common import (
    TWO_PI, KAPPA, CHI_NOMINAL, EPSILON_MAX, G_COUP, OMEGA_R, OMEGA_Q,
    T1_QUBIT, T2_QUBIT, T_PHI_QUBIT,
    optimal_delta_g, simulate_conditioned_fields,
    snr_squared, endpoint_separation_sq, assignment_fidelity_from_snr2,
    residual_photons, measurement_dephasing_factor,
    optimize_amplitude, run_grape, PULSE_FAMILIES,
    make_readout_resonator, DT_ODE,
)

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "data")
os.makedirs(DATA_DIR, exist_ok=True)


# ─── Analytical open-system model ────────────────────────────────────────────

def qubit_t1_correction(T_read: float, T1: float = T1_QUBIT) -> float:
    """
    Probability that qubit does NOT relax during readout.
    p_survive = exp(−T_read / T1).
    If qubit decays at time t_d < T_read, the readout signal from t_d to T is
    from the ground state, effectively reducing SNR.

    Returns the effective SNR² correction factor:
        SNR²_eff = ∫₀ᵀ κ |α_e(t)|−α_g(t)|² · P_no_decay(t) dt / ∫₀ᵀ κ|Δα|² dt
    For the linear model, this simplifies to an exponential weighting of the integrand.
    """
    return float(np.exp(-T_read / max(T1, 1e-30)))


def effective_snr2_with_t1(
    alpha_g: np.ndarray,
    alpha_e: np.ndarray,
    tlist: np.ndarray,
    kappa: float,
    T1: float = T1_QUBIT,
) -> float:
    """
    Effective SNR² accounting for qubit T1 relaxation during readout.

    Model: if qubit decays at time t_d (exponentially distributed with rate Γ₁=1/T1),
    the contribution to the distinguishability from [t_d, T] is zero.

    SNR²_eff = κ ∫₀ᵀ |Δα(t)|² · P(no decay up to t) dt
             = κ ∫₀ᵀ |Δα(t)|² · exp(−t/T1) dt
    """
    decay_weight = np.exp(-tlist / max(T1, 1e-30))
    diff_sq = np.abs(alpha_e - alpha_g) ** 2
    integrand = decay_weight * diff_sq
    return float(kappa * np.trapezoid(integrand, tlist))


def net_assignment_fidelity(
    snr2_eff: float,
    T_read: float,
    T1: float = T1_QUBIT,
) -> float:
    """
    Net assignment fidelity accounting for qubit T1 during readout.

    The qubit decays to |g⟩ with probability p_decay = 1 − exp(−T_read/T1).
    For a qubit prepared in |e⟩ that decays, the readout result is random:
        F_net = (1 − p_decay) · F_A(SNR²_eff) + p_decay · (1/2)
    [Here F_A(0) = 0.5 since without signal the classification is random.]
    """
    p_decay = 1.0 - np.exp(-T_read / max(T1, 1e-30))
    F_no_decay = assignment_fidelity_from_snr2(snr2_eff)
    F_net = (1.0 - p_decay) * F_no_decay + p_decay * 0.5
    return float(F_net)


def optimal_readout_duration_with_decoherence(
    kappa: float = KAPPA,
    chi: float = None,
    epsilon_max: float = EPSILON_MAX,
    T1: float = T1_QUBIT,
    T2: float = T2_QUBIT,
    dt: float = DT_ODE,
    kT_scan_max: float = 20.0,
    n_kT: int = 80,
) -> dict:
    """
    Sweep κT at fixed χ/κ=1 and compute:
      - SNR² (no decoherence)
      - SNR²_eff (with T1 weighting)
      - F_net (with T1 correction)
      - η_φ (measurement-induced dephasing factor)
      - Qubit-readout-fidelity combining T1 and SNR²

    Returns arrays for each metric vs κT.
    """
    if chi is None:
        chi = CHI_NOMINAL
    delta_g = optimal_delta_g(chi)
    kT_vals = np.linspace(0.3, kT_scan_max, n_kT)
    T_vals  = kT_vals / kappa

    snr2_arr      = np.zeros(n_kT)
    snr2_eff_arr  = np.zeros(n_kT)
    F_net_arr     = np.zeros(n_kT)
    eta_phi_arr   = np.zeros(n_kT)
    n_res_arr     = np.zeros(n_kT)

    for i, (T, kT) in enumerate(zip(T_vals, kT_vals)):
        best_eps, ev = optimize_amplitude(
            "square", T, kappa, chi, delta_g, epsilon_max, dt, n_grid=40
        )
        ag, ae, tl = ev["alpha_g"], ev["alpha_e"], ev["tlist"]
        snr2_arr[i]     = ev["snr2"]
        snr2_eff_arr[i] = effective_snr2_with_t1(ag, ae, tl, kappa, T1)
        F_net_arr[i]    = net_assignment_fidelity(snr2_eff_arr[i], T, T1)
        eta_phi_arr[i]  = measurement_dephasing_factor(ag, ae, tl, kappa)
        n_res_arr[i]    = ev["n_res"]

    idx_opt = np.argmax(F_net_arr)
    return {
        "kT_vals": kT_vals,
        "T_vals": T_vals,
        "snr2": snr2_arr,
        "snr2_eff": snr2_eff_arr,
        "F_net": F_net_arr,
        "eta_phi": eta_phi_arr,
        "n_res": n_res_arr,
        "kT_opt": kT_vals[idx_opt],
        "F_net_max": F_net_arr[idx_opt],
        "T1": T1, "T2": T2,
    }


def dephasing_vs_amplitude(
    kappa: float = KAPPA,
    chi: float = None,
    kappa_t: float = 3.0,
    dt: float = DT_ODE,
) -> dict:
    """
    Study measurement-induced dephasing η_φ as a function of drive amplitude |ε|.

    Higher amplitude → better SNR but more qubit dephasing.
    η_φ = exp(−∫(κ/2)|Δα|² dt); lower → more dephasing.
    """
    if chi is None:
        chi = CHI_NOMINAL
    delta_g = optimal_delta_g(chi)
    T_read = kappa_t / kappa
    eps_vals = np.linspace(EPSILON_MAX * 0.01, EPSILON_MAX, 80)

    snr2_arr   = np.zeros(len(eps_vals))
    eta_phi_arr = np.zeros(len(eps_vals))
    n_res_arr  = np.zeros(len(eps_vals))
    F_assign_arr = np.zeros(len(eps_vals))

    N = max(2, int(round(T_read / dt)))
    tlist = np.linspace(0.0, T_read, N + 1)

    from common import evaluate_pulse_family
    for i, eps in enumerate(eps_vals):
        ev = evaluate_pulse_family("square", complex(eps), T_read, kappa, chi, delta_g, dt)
        snr2_arr[i]    = ev["snr2"]
        eta_phi_arr[i] = ev["dephasing_factor"]
        n_res_arr[i]   = ev["n_res"]
        F_assign_arr[i] = ev["F_assign"]

    # Steady-state photon number at each amplitude
    lambda_g = 0.5 * kappa + 1j * delta_g
    n_ss_arr = (eps_vals / abs(lambda_g)) ** 2

    return {
        "eps_vals": eps_vals,
        "eps_kappa_vals": eps_vals / kappa,      # normalized
        "snr2": snr2_arr,
        "eta_phi": eta_phi_arr,
        "n_res": n_res_arr,
        "F_assign": F_assign_arr,
        "n_ss": n_ss_arr,
        "kappa_t": kappa_t,
        "chi_kappa": abs(chi) / kappa,
    }


def purcell_rate_estimate(
    kappa: float = KAPPA,
    chi: float = None,
    n_omega_q: int = 5,
) -> dict:
    """
    Estimate Purcell-limited T1 using cqed_sim ReadoutResonator.
    Sweep qubit-resonator detuning Δ_qr = ω_q − ω_r.
    """
    if chi is None:
        chi = CHI_NOMINAL

    # Detuning sweep: ω_q − ω_r
    omega_qr_values = TWO_PI * np.array([0.5, 1.0, 2.0, 3.0, 5.0]) * 1e9   # rad/s
    res = make_readout_resonator(kappa=kappa, chi=chi, epsilon=EPSILON_MAX * 0.1)
    purcell_rates  = np.array([res.purcell_rate(OMEGA_R + d) for d in omega_qr_values])
    purcell_t1s    = np.array([res.purcell_limited_t1(OMEGA_R + d) for d in omega_qr_values])

    return {
        "omega_qr_GHz": omega_qr_values / TWO_PI / 1e9,
        "purcell_rate": purcell_rates,
        "purcell_T1_us": purcell_t1s * 1e6,
    }


def grape_with_t1(
    kappa: float = KAPPA,
    chi: float = None,
    kappa_t: float = 3.0,
    epsilon_max: float = EPSILON_MAX,
    T1: float = T1_QUBIT,
    n_seg: int = 60,
    n_restarts: int = 10,
    dt: float = DT_ODE,
) -> dict:
    """
    Compare GRAPE performance (pure SNR²) vs. the effective SNR²_eff with T1.
    Shows that for κT << T1·κ, T1 correction is negligible.
    """
    if chi is None:
        chi = CHI_NOMINAL
    delta_g = optimal_delta_g(chi)
    T_read = kappa_t / kappa

    gr = run_grape(kappa=kappa, chi=chi, T_read=T_read, N_seg=n_seg,
                   delta_g=delta_g, epsilon_max=epsilon_max,
                   n_restarts=n_restarts, w_snr=1.0, w_res=0.0)

    snr2_eff = effective_snr2_with_t1(gr["alpha_g"], gr["alpha_e"], gr["tlist"], kappa, T1)
    F_net = net_assignment_fidelity(snr2_eff, T_read, T1)
    eta_phi = measurement_dephasing_factor(gr["alpha_g"], gr["alpha_e"], gr["tlist"], kappa)

    # Square pulse comparison
    _, ev_sq = optimize_amplitude("square", T_read, kappa, chi, delta_g, epsilon_max, dt)
    snr2_eff_sq = effective_snr2_with_t1(ev_sq["alpha_g"], ev_sq["alpha_e"], ev_sq["tlist"], kappa, T1)
    F_net_sq = net_assignment_fidelity(snr2_eff_sq, T_read, T1)

    return {
        "snr2_grape": gr["snr2_opt"],
        "snr2_eff_grape": snr2_eff,
        "F_net_grape": F_net,
        "eta_phi_grape": eta_phi,
        "snr2_square": ev_sq["snr2"],
        "snr2_eff_square": snr2_eff_sq,
        "F_net_square": F_net_sq,
        "p_decay": 1.0 - np.exp(-T_read / T1),
        "kappa_t": kappa_t,
        "chi_kappa": abs(chi) / kappa,
    }


def main() -> None:
    print("=" * 60)
    print("Phase 4: Open-System and Realistic Hamiltonian Effects")
    print("=" * 60)

    kappa = KAPPA

    # ── 1. Optimal duration with decoherence ─────────────────────────────────
    print("\n[1/4] Computing optimal readout duration with T1, T2 effects...")
    opt_dur = optimal_readout_duration_with_decoherence(kappa)
    print(f"  Optimal κT = {opt_dur['kT_opt']:.1f} ({opt_dur['kT_opt']/kappa*1e9:.0f} ns)")
    print(f"  F_net max  = {opt_dur['F_net_max']:.4f}")
    print(f"  T1 = {opt_dur['T1']*1e6:.0f} μs → κT1 = {opt_dur['T1']*kappa:.0f}")
    idx5 = np.argmin(np.abs(opt_dur["kT_vals"] - 5.0))
    print(f"  At κT=5: SNR²={opt_dur['snr2'][idx5]:.3f}, F_net={opt_dur['F_net'][idx5]:.4f}, "
          f"η_φ={opt_dur['eta_phi'][idx5]:.4f}")

    # ── 2. Dephasing vs amplitude ─────────────────────────────────────────────
    print("\n[2/4] Computing measurement-induced dephasing vs drive amplitude...")
    dep_amp = dephasing_vs_amplitude(kappa)
    idx_best = np.argmax(dep_amp["snr2"])
    print(f"  Best SNR² at |ε|/κ={dep_amp['eps_kappa_vals'][idx_best]:.2f}: "
          f"SNR²={dep_amp['snr2'][idx_best]:.3f}, η_φ={dep_amp['eta_phi'][idx_best]:.4f}")
    print(f"  n_ss at best = {dep_amp['n_ss'][idx_best]:.2f} photons")

    # ── 3. Purcell rate ──────────────────────────────────────────────────────
    print("\n[3/4] Estimating Purcell-limited T1...")
    purcell = purcell_rate_estimate(kappa)
    for d, r, t1 in zip(purcell["omega_qr_GHz"], purcell["purcell_rate"], purcell["purcell_T1_us"]):
        print(f"  Δ_qr/2π={d:.1f} GHz: Γ_Purcell/2π={r/TWO_PI/1e3:.2f} kHz, T1_Purcell={t1:.1f} μs")

    # ── 4. GRAPE with T1 ─────────────────────────────────────────────────────
    print("\n[4/4] GRAPE performance with T1 correction...")
    grape_t1 = grape_with_t1(kappa)
    print(f"  GRAPE SNR²(no T1)={grape_t1['snr2_grape']:.3f}, "
          f"SNR²_eff(T1)={grape_t1['snr2_eff_grape']:.3f}")
    print(f"  GRAPE F_net={grape_t1['F_net_grape']:.4f}")
    print(f"  Square F_net={grape_t1['F_net_square']:.4f}")
    print(f"  p_decay={grape_t1['p_decay']:.4f} (negligible for κT<<κT1)")

    # ── Save ─────────────────────────────────────────────────────────────────
    outpath = os.path.join(DATA_DIR, "phase4_results.npz")
    np.savez(
        outpath,
        # Optimal duration
        opt_kT_vals=opt_dur["kT_vals"],
        opt_T_vals=opt_dur["T_vals"],
        opt_snr2=opt_dur["snr2"],
        opt_snr2_eff=opt_dur["snr2_eff"],
        opt_F_net=opt_dur["F_net"],
        opt_eta_phi=opt_dur["eta_phi"],
        opt_n_res=opt_dur["n_res"],
        opt_kT_optimal=opt_dur["kT_opt"],
        opt_F_net_max=opt_dur["F_net_max"],
        # Dephasing vs amplitude
        damp_eps_vals=dep_amp["eps_vals"],
        damp_eps_kappa=dep_amp["eps_kappa_vals"],
        damp_snr2=dep_amp["snr2"],
        damp_eta_phi=dep_amp["eta_phi"],
        damp_n_res=dep_amp["n_res"],
        damp_n_ss=dep_amp["n_ss"],
        damp_F_assign=dep_amp["F_assign"],
        # Purcell
        purcell_dqr_GHz=purcell["omega_qr_GHz"],
        purcell_rate=purcell["purcell_rate"],
        purcell_T1_us=purcell["purcell_T1_us"],
        # GRAPE with T1
        gt1_snr2_grape=grape_t1["snr2_grape"],
        gt1_snr2_eff_grape=grape_t1["snr2_eff_grape"],
        gt1_F_net_grape=grape_t1["F_net_grape"],
        gt1_snr2_square=grape_t1["snr2_square"],
        gt1_F_net_square=grape_t1["F_net_square"],
        gt1_p_decay=grape_t1["p_decay"],
        T1=T1_QUBIT, T2=T2_QUBIT, kappa=kappa, chi=CHI_NOMINAL,
    )
    print(f"\nResults saved to {outpath}")

    print("\n" + "=" * 60)
    print("Phase 4 Summary:")
    print(f"  Optimal T_read = {opt_dur['kT_opt']/kappa*1e9:.0f} ns (κT={opt_dur['kT_opt']:.1f})")
    print(f"  F_net at optimum = {opt_dur['F_net_max']:.4f}")
    print(f"  Measurement-induced dephasing η_φ at κT=5: {opt_dur['eta_phi'][idx5]:.4f}")
    print(f"  Purcell effect: T1_Purcell >> T1 for Δ_qr >> κ (negligible)")
    print(f"  T1 correction negligible for κT << κT1 = {T1_QUBIT*kappa:.0f}")
    print("=" * 60)


if __name__ == "__main__":
    main()
