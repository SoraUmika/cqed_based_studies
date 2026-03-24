"""
Phase 5: Open-system dynamics with T1, T2 decoherence.

Model includes higher-order corrections (χ', K) from Phase 4 plus
Lindblad decoherence with T1 = T2 = 20 μs.

For density-matrix outputs, fidelity is computed as average state fidelity
against closed-system ideal final states (from Phase-4 model).

Usage:
    python scripts/run_phase5_open_system.py

Output:
    data/phase5_results.npz
"""

import sys
import time
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import runtime_compat  # noqa: F401

from common import (
    CHI, CHI_PRIME, KERR, CHI_T_VALUES, DT, N_CAV, N_FOCK, N_TR,
    THETA_TARGET, PHI_TARGET, TARGET_N0, T1, T2,
    build_frame, build_model, build_noise_spec, duration_from_chi_t,
    state_fidelity_dm, z_corrected_target_fidelity, spectator_z_fidelity,
    target_qubit_unitary, extract_branch_unitaries,
)
from run_phase1_phase2 import (
    build_single_tone_gaussian, build_multitone_pulse,
    build_square_pulse, build_cosine_squared_pulse,
)


def simulate_open_system(model, frame, pulses, drive_ops, duration, noise_spec):
    """Simulate with Lindblad noise, returning density matrices.

    Returns (ideal_final_kets, noisy_final_dms, ideal_unitaries).
    """
    from cqed_sim.sequence import SequenceCompiler
    from cqed_sim.sim import SimulationConfig, prepare_simulation

    compiler = SequenceCompiler(dt=DT)
    compiled = compiler.compile(pulses, t_end=duration + 4 * DT)
    config = SimulationConfig(frame=frame, store_states=False)

    # Closed-system session (ideal reference)
    session_ideal = prepare_simulation(
        model, compiled, drive_ops, config=config, e_ops={})
    # Open-system session (with noise)
    session_noisy = prepare_simulation(
        model, compiled, drive_ops, config=config, noise=noise_spec, e_ops={})

    ideal_kets = []
    noisy_dms = []
    for n in range(N_FOCK):
        for q in [0, 1]:
            psi0 = model.basis_state(q, n)
            # Ideal (closed-system)
            res_ideal = session_ideal.run(psi0)
            ideal_kets.append(res_ideal.final_state)
            # Noisy (open-system → density matrix)
            res_noisy = session_noisy.run(psi0)
            noisy_dms.append(res_noisy.final_state)

    ideal_unitaries = extract_branch_unitaries(ideal_kets, model, N_FOCK)
    return ideal_kets, noisy_dms, ideal_unitaries


def compute_open_system_metrics(model, ideal_kets, noisy_dms, ideal_unitaries, n0, theta, phi):
    """Compute state fidelities for open-system outputs.

    Computes three metrics per branch:
      1. F_deco: state fidelity(noisy, closed-system) — decoherence contribution
      2. F_net:  state fidelity(noisy, ideal_cphase_target) — total performance
      3. F_ideal: closed-system cphase SQR fidelity — reference

    For computational-basis inputs, F_net reduces to population fidelity:
      - Target branch: checks population transfer |g⟩→|e⟩
      - Spectator branches: checks population preservation

    Returns dict with branch-resolved and subspace-averaged values.
    """
    U_target = target_qubit_unitary(theta, phi)

    # Get Z-corrections from ideal unitaries (closed-system)
    _, alpha_opt = z_corrected_target_fidelity(ideal_unitaries[n0], U_target)

    branch_fid_deco = np.zeros(N_FOCK)
    branch_fid_net = np.zeros(N_FOCK)
    branch_fid_ideal = np.zeros(N_FOCK)

    for n in range(N_FOCK):
        idx_g = 2 * n
        idx_e = 2 * n + 1
        ideal_g = ideal_kets[idx_g]
        ideal_e = ideal_kets[idx_e]
        rho_g = noisy_dms[idx_g]
        rho_e = noisy_dms[idx_e]

        # --- Decoherence fidelity: noisy vs closed-system ---
        F_deco_g = state_fidelity_dm(ideal_g, rho_g)
        F_deco_e = state_fidelity_dm(ideal_e, rho_e)
        branch_fid_deco[n] = (F_deco_g + F_deco_e) / 2.0

        # --- Net fidelity: noisy vs ideal cphase target ---
        # For target branch (R_X(π)): |g,n0⟩→|e,n0⟩, |e,n0⟩→|g,n0⟩
        # For spectator (identity): |g,n⟩→|g,n⟩, |e,n⟩→|e,n⟩
        if n == n0:
            target_g = model.basis_state(1, n)  # R_X(π)|g⟩ → |e⟩
            target_e = model.basis_state(0, n)  # R_X(π)|e⟩ → |g⟩
        else:
            target_g = model.basis_state(0, n)  # identity
            target_e = model.basis_state(1, n)
        F_net_g = state_fidelity_dm(target_g, rho_g)
        F_net_e = state_fidelity_dm(target_e, rho_e)
        branch_fid_net[n] = (F_net_g + F_net_e) / 2.0

        # --- Ideal cphase fidelity (closed-system reference) ---
        if n == n0:
            fid_ideal, _ = z_corrected_target_fidelity(ideal_unitaries[n], U_target)
        else:
            fid_ideal, _ = spectator_z_fidelity(ideal_unitaries[n])
        branch_fid_ideal[n] = fid_ideal

    return {
        "cphase_fid_ideal": np.mean(branch_fid_ideal),
        "cphase_fid_net": np.mean(branch_fid_net),
        "deco_fid": np.mean(branch_fid_deco),
        "branch_fid_deco": branch_fid_deco,
        "branch_fid_net": branch_fid_net,
        "branch_fid_ideal": branch_fid_ideal,
    }


def run_scan():
    """Run the full χT scan with higher-order corrections and decoherence."""
    model_p5 = build_model(chi_prime=CHI_PRIME, kerr=KERR)
    frame_p5 = build_frame(model_p5)
    noise_spec = build_noise_spec()

    n_chi_t = len(CHI_T_VALUES)
    n_families = 4
    family_names = [
        "single_tone_gaussian", "square",
        "cosine_squared", "multitone_one_segment",
    ]

    cphase_fid_net = np.zeros((n_families, n_chi_t))
    cphase_fid_ideal = np.zeros((n_families, n_chi_t))
    deco_fid = np.zeros((n_families, n_chi_t))
    branch_fid_net_all = np.zeros((n_families, n_chi_t, N_FOCK))
    branch_fid_deco_all = np.zeros((n_families, n_chi_t, N_FOCK))
    branch_fid_ideal_all = np.zeros((n_families, n_chi_t, N_FOCK))

    print(f"Phase 5: {n_families} families × {n_chi_t} χT/(2π) values (open-system)")
    print(f"Model: χ'=2π×{CHI_PRIME/(2*np.pi)/1e3:.0f} kHz, K=2π×{KERR/(2*np.pi)/1e3:.0f} kHz")
    print(f"Noise: T1={T1*1e6:.0f} μs, T2={T2*1e6:.0f} μs")
    print(f"Target: n0={TARGET_N0}, θ=π, φ=0")
    print()

    t_start = time.time()

    for i_fam, fam_name in enumerate(family_names):
        print(f"  Family: {fam_name}")

        for i_ct, chi_t in enumerate(CHI_T_VALUES):
            T = duration_from_chi_t(chi_t)

            if fam_name == "single_tone_gaussian":
                pulses, drive_ops = build_single_tone_gaussian(
                    model_p5, frame_p5, TARGET_N0, THETA_TARGET, PHI_TARGET, T)
            elif fam_name == "square":
                pulses, drive_ops = build_square_pulse(
                    model_p5, frame_p5, TARGET_N0, THETA_TARGET, PHI_TARGET, T)
            elif fam_name == "cosine_squared":
                pulses, drive_ops = build_cosine_squared_pulse(
                    model_p5, frame_p5, TARGET_N0, THETA_TARGET, PHI_TARGET, T)
            elif fam_name == "multitone_one_segment":
                pulses, drive_ops = build_multitone_pulse(
                    model_p5, frame_p5, N_FOCK, TARGET_N0, THETA_TARGET, PHI_TARGET, T)

            ideal_kets, noisy_dms, ideal_unitaries = simulate_open_system(
                model_p5, frame_p5, pulses, drive_ops, T, noise_spec)

            metrics = compute_open_system_metrics(
                model_p5, ideal_kets, noisy_dms, ideal_unitaries, TARGET_N0,
                THETA_TARGET, PHI_TARGET)

            cphase_fid_net[i_fam, i_ct] = metrics["cphase_fid_net"]
            cphase_fid_ideal[i_fam, i_ct] = metrics["cphase_fid_ideal"]
            deco_fid[i_fam, i_ct] = metrics["deco_fid"]
            branch_fid_net_all[i_fam, i_ct] = metrics["branch_fid_net"]
            branch_fid_deco_all[i_fam, i_ct] = metrics["branch_fid_deco"]
            branch_fid_ideal_all[i_fam, i_ct] = metrics["branch_fid_ideal"]

            print(f"    χT/(2π)={chi_t:5.2f}  T={T*1e6:6.2f} μs  "
                  f"F_ideal={metrics['cphase_fid_ideal']:.6f}  "
                  f"F_net={metrics['cphase_fid_net']:.6f}  "
                  f"F_deco={metrics['deco_fid']:.6f}")

    elapsed = time.time() - t_start
    print(f"\nPhase 5 scan complete in {elapsed:.1f} s")

    data_dir = SCRIPT_DIR.parent / "data"
    np.savez(
        data_dir / "phase5_results.npz",
        chi_t_values=CHI_T_VALUES,
        family_names=family_names,
        cphase_fid_net=cphase_fid_net,
        cphase_fid_ideal=cphase_fid_ideal,
        deco_fid=deco_fid,
        branch_fid_net=branch_fid_net_all,
        branch_fid_deco=branch_fid_deco_all,
        branch_fid_ideal=branch_fid_ideal_all,
        n_fock=N_FOCK,
        target_n0=TARGET_N0,
        theta_target=THETA_TARGET,
        phi_target=PHI_TARGET,
        chi_rad_s=CHI,
        chi_prime_rad_s=CHI_PRIME,
        kerr_rad_s=KERR,
        t1_s=T1,
        t2_s=T2,
    )
    print(f"Results saved to {data_dir / 'phase5_results.npz'}")


if __name__ == "__main__":
    run_scan()
