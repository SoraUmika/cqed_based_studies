"""
Phase 4: Higher-order corrections (χ', cavity self-Kerr K).

Same scan as Phase 1+2 but with the full dispersive model including:
  χ' = 2π × (-21 kHz)   — second-order dispersive shift
  K  = 2π × (-28 kHz)   — cavity self-Kerr

Pulse frequencies are recalibrated to the corrected branch frequencies.

Usage:
    python scripts/run_phase4_higher_order.py

Output:
    data/phase4_results.npz
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
    THETA_TARGET, PHI_TARGET, TARGET_N0, build_frame, build_model,
    branch_frequencies, duration_from_chi_t, extract_leakage,
)
from run_phase1_phase2 import (
    build_single_tone_gaussian, build_multitone_pulse,
    build_square_pulse, build_cosine_squared_pulse,
    simulate_waveform, compute_metrics,
)

# ---------------------------------------------------------------------------
# Modified simulation using Phase-4 model
# ---------------------------------------------------------------------------

def simulate_waveform_p4(model, frame, pulses, drive_ops, duration):
    """Same as simulate_waveform but accepts pre-built model/frame."""
    from cqed_sim.sequence import SequenceCompiler
    from cqed_sim.sim import SimulationConfig, prepare_simulation

    compiler = SequenceCompiler(dt=DT)
    compiled = compiler.compile(pulses, t_end=duration + 4 * DT)
    config = SimulationConfig(frame=frame, store_states=False)
    session = prepare_simulation(model, compiled, drive_ops, config=config, e_ops={})

    initial_states = []
    for n in range(N_FOCK):
        initial_states.append(model.basis_state(0, n))
        initial_states.append(model.basis_state(1, n))

    final_states = []
    for psi0 in initial_states:
        result = session.run(psi0)
        final_states.append(result.final_state)

    from common import extract_branch_unitaries, extract_leakage
    unitaries = extract_branch_unitaries(final_states, model, N_FOCK)
    return final_states, unitaries


def run_scan():
    """Run the full χT scan with higher-order corrections."""
    model_p4 = build_model(chi_prime=CHI_PRIME, kerr=KERR)
    frame_p4 = build_frame(model_p4)

    # Show branch frequency shifts
    model_ideal = build_model()
    frame_ideal = build_frame(model_ideal)
    freqs_ideal = branch_frequencies(model_ideal, frame_ideal, N_FOCK)
    freqs_p4 = branch_frequencies(model_p4, frame_p4, N_FOCK)
    print("Branch frequency shifts from χ', K:")
    for n in range(N_FOCK):
        delta = (freqs_p4[n] - freqs_ideal[n]) / (2 * np.pi * 1e3)
        print(f"  n={n}: Δf = {delta:+.1f} kHz")
    print()

    n_chi_t = len(CHI_T_VALUES)
    n_families = 4
    family_names = [
        "single_tone_gaussian", "square",
        "cosine_squared", "multitone_one_segment",
    ]

    true_sqr_fid = np.zeros((n_families, n_chi_t))
    cphase_sqr_fid = np.zeros((n_families, n_chi_t))
    target_fid = np.zeros((n_families, n_chi_t))
    target_z_phase = np.zeros((n_families, n_chi_t))
    spec_phase_spread = np.zeros((n_families, n_chi_t))
    spec_max_transverse = np.zeros((n_families, n_chi_t))
    branch_fid_true_all = np.zeros((n_families, n_chi_t, N_FOCK))
    branch_fid_cphase_all = np.zeros((n_families, n_chi_t, N_FOCK))
    spec_phases_all = np.zeros((n_families, n_chi_t, N_FOCK))
    leakage_all = np.zeros((n_families, n_chi_t, N_FOCK))

    print(f"Phase 4: {n_families} families × {n_chi_t} χT/(2π) values")
    print(f"Model: χ'=2π×{CHI_PRIME/(2*np.pi)/1e3:.0f} kHz, K=2π×{KERR/(2*np.pi)/1e3:.0f} kHz")
    print(f"Target: n0={TARGET_N0}, θ=π, φ=0")
    print()

    t_start = time.time()

    for i_fam, fam_name in enumerate(family_names):
        print(f"  Family: {fam_name}")

        for i_ct, chi_t in enumerate(CHI_T_VALUES):
            T = duration_from_chi_t(chi_t)

            if fam_name == "single_tone_gaussian":
                pulses, drive_ops = build_single_tone_gaussian(
                    model_p4, frame_p4, TARGET_N0, THETA_TARGET, PHI_TARGET, T)
            elif fam_name == "square":
                pulses, drive_ops = build_square_pulse(
                    model_p4, frame_p4, TARGET_N0, THETA_TARGET, PHI_TARGET, T)
            elif fam_name == "cosine_squared":
                pulses, drive_ops = build_cosine_squared_pulse(
                    model_p4, frame_p4, TARGET_N0, THETA_TARGET, PHI_TARGET, T)
            elif fam_name == "multitone_one_segment":
                pulses, drive_ops = build_multitone_pulse(
                    model_p4, frame_p4, N_FOCK, TARGET_N0, THETA_TARGET, PHI_TARGET, T)

            final_states, unitaries = simulate_waveform_p4(model_p4, frame_p4, pulses, drive_ops, T)
            leakage_all[i_fam, i_ct] = extract_leakage(final_states, model_p4, N_FOCK)
            metrics = compute_metrics(unitaries, TARGET_N0, THETA_TARGET, PHI_TARGET)

            true_sqr_fid[i_fam, i_ct] = metrics["true_sqr_fidelity"]
            cphase_sqr_fid[i_fam, i_ct] = metrics["cphase_sqr_fidelity"]
            target_fid[i_fam, i_ct] = metrics["target_fidelity"]
            target_z_phase[i_fam, i_ct] = metrics["target_z_phase"]
            branch_fid_true_all[i_fam, i_ct] = metrics["branch_fidelities_true"]
            branch_fid_cphase_all[i_fam, i_ct] = metrics["branch_fidelities_cphase"]
            spec_phases_all[i_fam, i_ct] = metrics["spectator_phases"]

            sp = metrics["spectator_phases"]
            spec_mask = np.arange(N_FOCK) != TARGET_N0
            if np.sum(spec_mask) > 1:
                spec_phase_spread[i_fam, i_ct] = np.ptp(sp[spec_mask])
            spec_max_transverse[i_fam, i_ct] = np.max(
                metrics["spectator_transverse_errors"][spec_mask])

            print(f"    χT/(2π)={chi_t:5.2f}  T={T*1e9:8.1f} ns  "
                  f"F_true={metrics['true_sqr_fidelity']:.6f}  "
                  f"F_cphase={metrics['cphase_sqr_fidelity']:.6f}  "
                  f"F_target={metrics['target_fidelity']:.6f}")

    elapsed = time.time() - t_start
    print(f"\nPhase 4 scan complete in {elapsed:.1f} s")

    data_dir = SCRIPT_DIR.parent / "data"
    np.savez(
        data_dir / "phase4_results.npz",
        chi_t_values=CHI_T_VALUES,
        family_names=family_names,
        true_sqr_fidelity=true_sqr_fid,
        cphase_sqr_fidelity=cphase_sqr_fid,
        target_fidelity=target_fid,
        target_z_phase=target_z_phase,
        spectator_phase_spread=spec_phase_spread,
        spectator_max_transverse=spec_max_transverse,
        branch_fidelities_true=branch_fid_true_all,
        branch_fidelities_cphase=branch_fid_cphase_all,
        spectator_phases=spec_phases_all,
        leakage=leakage_all,
        n_fock=N_FOCK,
        target_n0=TARGET_N0,
        theta_target=THETA_TARGET,
        phi_target=PHI_TARGET,
        chi_rad_s=CHI,
        chi_prime_rad_s=CHI_PRIME,
        kerr_rad_s=KERR,
    )
    print(f"Results saved to {data_dir / 'phase4_results.npz'}")


if __name__ == "__main__":
    run_scan()
