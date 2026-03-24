"""
Phase 1+2: Baseline single-tone Gaussian SQR and multitone comparison.

Scans χT for the four waveform families and two target types (true SQR,
conditional-phase SQR), computing branch-resolved process fidelity, spectator
phase, and overall subspace fidelity.

Usage:
    python scripts/run_phase1_phase2.py

Output:
    data/phase1_phase2_results.npz
    figures/fidelity_vs_chi_t.{png,pdf}
    figures/spectator_phase_vs_chi_t.{png,pdf}
    figures/branch_fidelity_heatmap.{png,pdf}
"""

import sys
import time
from pathlib import Path

import numpy as np

# Ensure the scripts directory is importable
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import runtime_compat  # noqa: F401

from common import (
    CHI, CHI_T_VALUES, DT, N_CAV, N_FOCK, N_TR, THETA_TARGET, PHI_TARGET,
    TARGET_N0, TOL_BRIGHT, build_frame, build_model, branch_frequencies,
    conditional_process_fidelity, duration_from_chi_t,
    extract_branch_unitaries, extract_leakage,
    identity_fidelity_with_z, spectator_transverse_error,
    spectator_z_fidelity, target_qubit_unitary, z_corrected_target_fidelity,
)

from cqed_sim.core import FrameSpec
from cqed_sim.core.frequencies import (
    carrier_for_transition_frequency,
    manifold_transition_frequency,
)
from cqed_sim.pulses import Pulse
from cqed_sim.pulses.envelopes import (
    gaussian_envelope,
    normalized_gaussian,
    gaussian_area_fraction,
    multitone_gaussian_envelope,
    MultitoneTone,
)
from cqed_sim.sequence import SequenceCompiler
from cqed_sim.sim import SimulationConfig, simulate_sequence, prepare_simulation

# ---------------------------------------------------------------------------
# Helpers: Pulse builders for each waveform family
# ---------------------------------------------------------------------------

SIGMA_FRAC = 1.0 / 6.0
GAUSS_AREA = gaussian_area_fraction(SIGMA_FRAC)


def _make_gaussian_envelope(sigma_frac=SIGMA_FRAC):
    """Return a normalized Gaussian envelope callable with given sigma fraction."""
    def env(t_rel):
        return normalized_gaussian(t_rel, sigma_fraction=sigma_frac)
    return env


def build_single_tone_gaussian(model, frame, n0, theta, phi, duration):
    """Build a single-tone Gaussian selective pulse targeting branch n0.

    Returns (pulse_list, drive_ops).
    """
    omega_n0 = manifold_transition_frequency(model, n0, frame)
    carrier = carrier_for_transition_frequency(omega_n0)
    # RWA: rotation angle θ = 2 * amp * T * area_fraction
    # With normalized_gaussian (unit area), θ = 2 * amp * T
    amp = theta / (2 * duration)

    pulse = Pulse(
        channel="q",
        t0=0.0,
        duration=duration,
        envelope=_make_gaussian_envelope(),
        carrier=carrier,
        phase=phi,
        amp=amp,
    )
    return [pulse], {"q": "qubit"}


def build_multitone_pulse(model, frame, n_fock, n0, theta, phi, duration):
    """Build a one-segment common-envelope multitone pulse.

    Active tone on branch n0 with amplitude θ/(2T·area); spectator tones with
    zero amplitude (baseline). Uses the library's multitone Gaussian envelope.

    The tone carrier frequencies are set to carrier_for_transition_frequency(ω_n)
    to match the exp(+iωt) sign convention used by the Pulse class.

    Returns (pulse_list, drive_ops).
    """
    sigma_frac = 1.0 / 6.0
    tone_specs = []

    for n in range(n_fock):
        omega_n = manifold_transition_frequency(model, n, frame)
        # Use carrier frequency (−ω_n) to match Pulse sign convention
        tone_carrier = carrier_for_transition_frequency(omega_n)
        if n == n0:
            # multitone_gaussian_envelope uses normalized_gaussian internally
            # (unit area), so no GAUSS_AREA correction needed
            amp_n = theta / (2 * duration)
        else:
            amp_n = 0.0
        tone_specs.append(MultitoneTone(
            manifold=n,
            omega_rad_s=tone_carrier,
            amp_rad_s=amp_n,
            phase_rad=phi if n == n0 else 0.0,
        ))

    def env(t_rel):
        return multitone_gaussian_envelope(
            t_rel, duration_s=duration, sigma_fraction=sigma_frac,
            tone_specs=tone_specs,
        )

    pulse = Pulse(
        channel="q",
        t0=0.0,
        duration=duration,
        envelope=env,
        carrier=0.0,
        phase=0.0,
        amp=1.0,
    )
    return [pulse], {"q": "qubit"}


def build_square_pulse(model, frame, n0, theta, phi, duration):
    """Build a square (rectangular) selective pulse targeting branch n0.

    The constant envelope has sinc-like spectral response, with stronger
    sidelobes than a Gaussian. This probes the effect of spectral leakage on
    spectator branches.

    Returns (pulse_list, drive_ops).
    """
    omega_n0 = manifold_transition_frequency(model, n0, frame)
    carrier = carrier_for_transition_frequency(omega_n0)
    # Constant envelope has unit area: ∫₀¹ 1 dt_rel = 1, so θ = 2·amp·T
    amp = theta / (2 * duration)

    def env_square(t_rel):
        return np.ones_like(np.asarray(t_rel, dtype=float))

    pulse = Pulse(
        channel="q",
        t0=0.0,
        duration=duration,
        envelope=env_square,
        carrier=carrier,
        phase=phi,
        amp=amp,
    )
    return [pulse], {"q": "qubit"}


def build_cosine_squared_pulse(model, frame, n0, theta, phi, duration):
    """Build a cosine-squared (Hann) selective pulse targeting branch n0.

    env(t_rel) = 2·cos²(π(t_rel − 0.5)) for t_rel ∈ [0, 1].
    This smooth envelope has very low spectral sidelobes, offering better
    spectator isolation than a Gaussian at the cost of slightly broader mainlobe.
    Area ∫₀¹ env dt_rel = 1, so θ = 2·amp·T.

    Returns (pulse_list, drive_ops).
    """
    omega_n0 = manifold_transition_frequency(model, n0, frame)
    carrier = carrier_for_transition_frequency(omega_n0)
    amp = theta / (2 * duration)

    def env_cos2(t_rel):
        t = np.asarray(t_rel, dtype=float)
        return 2.0 * np.cos(np.pi * (t - 0.5))**2

    pulse = Pulse(
        channel="q",
        t0=0.0,
        duration=duration,
        envelope=env_cos2,
        carrier=carrier,
        phase=phi,
        amp=amp,
    )
    return [pulse], {"q": "qubit"}


# ---------------------------------------------------------------------------
# Core simulation: run one waveform at one χT value
# ---------------------------------------------------------------------------

def simulate_waveform(model, frame, pulses, drive_ops, duration):
    """Simulate the pulse with all 2*N_FOCK initial basis states.

    Returns list of final states and the branch-resolved 2×2 unitaries.
    """
    compiler = SequenceCompiler(dt=DT)
    compiled = compiler.compile(pulses, t_end=duration + 4 * DT)
    config = SimulationConfig(frame=frame, store_states=False)

    session = prepare_simulation(
        model, compiled, drive_ops,
        config=config, e_ops={},
    )

    initial_states = []
    for n in range(N_FOCK):
        initial_states.append(model.basis_state(0, n))  # |g, n⟩
        initial_states.append(model.basis_state(1, n))  # |e, n⟩

    final_states = []
    for psi0 in initial_states:
        result = session.run(psi0)
        final_states.append(result.final_state)

    unitaries = extract_branch_unitaries(final_states, model, N_FOCK)
    return final_states, unitaries


# ---------------------------------------------------------------------------
# Metrics for one waveform at one χT
# ---------------------------------------------------------------------------

def compute_metrics(unitaries, n0, theta, phi):
    """Compute all fidelity metrics from branch-resolved unitaries.

    The target branch fidelity is Z-corrected: we find the optimal virtual-Z
    frame update α that maximizes F_proc(Z(α) U_n0, R(θ,φ)).

    True SQR fidelity applies the *same* global Z(α) to all branches, then
    measures spectator identity fidelity F_proc(Z(α) U_n, I).

    Conditional-phase SQR allows each spectator branch its own Z-correction,
    so spectator fidelity is max_φ F_proc(Z(φ) U_n, I) per branch.

    Returns dict with:
      - target_fidelity: Z-corrected process fidelity on the target branch
      - target_z_phase: optimal virtual-Z phase α
      - true_sqr_fidelity: average branch fidelity for true SQR target
      - cphase_sqr_fidelity: average branch fidelity for conditional-phase SQR
      - spectator_phases: array of best-fit Z-phases for spectator branches
      - spectator_transverse_errors: array of transverse errors for spectators
      - branch_fidelities_true: per-branch fidelity for true SQR
      - branch_fidelities_cphase: per-branch fidelity for conditional-phase SQR
    """
    U_target = target_qubit_unitary(theta, phi)

    # Z-corrected target branch fidelity & optimal frame update
    target_fid, alpha_opt = z_corrected_target_fidelity(unitaries[n0], U_target)

    # Per-branch metrics
    branch_fid_true = np.zeros(N_FOCK)
    branch_fid_cphase = np.zeros(N_FOCK)
    spec_phases = np.zeros(N_FOCK)
    spec_transverse = np.zeros(N_FOCK)

    for n in range(N_FOCK):
        if n == n0:
            branch_fid_true[n] = target_fid
            branch_fid_cphase[n] = target_fid
            spec_phases[n] = alpha_opt
            spec_transverse[n] = 0.0
        else:
            # True SQR: same global Z(alpha_opt) applied to spectator, vs identity
            branch_fid_true[n] = identity_fidelity_with_z(unitaries[n], alpha_opt)
            # Conditional-phase SQR: each spectator gets its own optimal Z
            z_fid, z_phase = spectator_z_fidelity(unitaries[n])
            branch_fid_cphase[n] = z_fid
            spec_phases[n] = z_phase
            spec_transverse[n] = spectator_transverse_error(unitaries[n])

    # Subspace fidelities (average over branches, equally weighted)
    true_sqr_fid = np.mean(branch_fid_true)
    cphase_sqr_fid = np.mean(branch_fid_cphase)

    return {
        "target_fidelity": target_fid,
        "target_z_phase": alpha_opt,
        "true_sqr_fidelity": true_sqr_fid,
        "cphase_sqr_fidelity": cphase_sqr_fid,
        "spectator_phases": spec_phases,
        "spectator_transverse_errors": spec_transverse,
        "branch_fidelities_true": branch_fid_true,
        "branch_fidelities_cphase": branch_fid_cphase,
    }


# ---------------------------------------------------------------------------
# Main scan
# ---------------------------------------------------------------------------

WAVEFORM_FAMILIES = {
    "single_tone_gaussian": build_single_tone_gaussian,
    "square": None,              # special handling
    "cosine_squared": None,      # special handling
    "multitone_one_segment": None,  # special handling
}


def run_scan():
    """Run the full χT scan for all waveform families."""
    model = build_model()
    frame = build_frame(model)

    n_chi_t = len(CHI_T_VALUES)
    n_families = 4
    family_names = [
        "single_tone_gaussian",
        "square",
        "cosine_squared",
        "multitone_one_segment",
    ]

    # Result arrays
    true_sqr_fid = np.zeros((n_families, n_chi_t))
    cphase_sqr_fid = np.zeros((n_families, n_chi_t))
    target_fid = np.zeros((n_families, n_chi_t))
    target_z_phase = np.zeros((n_families, n_chi_t))
    spec_phase_spread = np.zeros((n_families, n_chi_t))
    spec_max_transverse = np.zeros((n_families, n_chi_t))
    branch_fid_true_all = np.zeros((n_families, n_chi_t, N_FOCK))
    branch_fid_cphase_all = np.zeros((n_families, n_chi_t, N_FOCK))
    spec_phases_all = np.zeros((n_families, n_chi_t, N_FOCK))

    # Leakage tracking
    leakage_all = np.zeros((n_families, n_chi_t, N_FOCK))

    print(f"Running scan: {n_families} families × {n_chi_t} χT/(2π) values")
    print(f"Target: n0={TARGET_N0}, θ=π, φ=0 (X_π)")
    print(f"Truncated Fock space: n=0..{N_FOCK-1}, n_tr={N_TR}")
    print()

    t_start = time.time()

    for i_fam, fam_name in enumerate(family_names):
        print(f"  Family: {fam_name}")

        for i_ct, chi_t in enumerate(CHI_T_VALUES):
            T = duration_from_chi_t(chi_t)

            # Build pulse for this family
            if fam_name == "single_tone_gaussian":
                pulses, drive_ops = build_single_tone_gaussian(
                    model, frame, TARGET_N0, THETA_TARGET, PHI_TARGET, T,
                )
            elif fam_name == "square":
                pulses, drive_ops = build_square_pulse(
                    model, frame, TARGET_N0, THETA_TARGET, PHI_TARGET, T,
                )
            elif fam_name == "cosine_squared":
                pulses, drive_ops = build_cosine_squared_pulse(
                    model, frame, TARGET_N0, THETA_TARGET, PHI_TARGET, T,
                )
            elif fam_name == "multitone_one_segment":
                pulses, drive_ops = build_multitone_pulse(
                    model, frame, N_FOCK, TARGET_N0, THETA_TARGET, PHI_TARGET, T,
                )
            else:
                raise ValueError(f"Unknown family: {fam_name}")

            # Simulate
            final_states, unitaries = simulate_waveform(model, frame, pulses, drive_ops, T)

            # Leakage
            leakage_all[i_fam, i_ct] = extract_leakage(final_states, model, N_FOCK)

            # Metrics
            metrics = compute_metrics(unitaries, TARGET_N0, THETA_TARGET, PHI_TARGET)

            true_sqr_fid[i_fam, i_ct] = metrics["true_sqr_fidelity"]
            cphase_sqr_fid[i_fam, i_ct] = metrics["cphase_sqr_fidelity"]
            target_fid[i_fam, i_ct] = metrics["target_fidelity"]
            target_z_phase[i_fam, i_ct] = metrics["target_z_phase"]
            branch_fid_true_all[i_fam, i_ct] = metrics["branch_fidelities_true"]
            branch_fid_cphase_all[i_fam, i_ct] = metrics["branch_fidelities_cphase"]
            spec_phases_all[i_fam, i_ct] = metrics["spectator_phases"]

            # Spectator phase spread (max - min among non-target branches)
            sp = metrics["spectator_phases"]
            spec_mask = np.arange(N_FOCK) != TARGET_N0
            if np.sum(spec_mask) > 1:
                spec_phase_spread[i_fam, i_ct] = np.ptp(sp[spec_mask])
            else:
                spec_phase_spread[i_fam, i_ct] = 0.0
            spec_max_transverse[i_fam, i_ct] = np.max(
                metrics["spectator_transverse_errors"][spec_mask]
            )

            leak_max = np.max(leakage_all[i_fam, i_ct])
            print(f"    χT/(2π)={chi_t:5.2f}  T={T*1e9:8.1f} ns  "
                  f"F_true={metrics['true_sqr_fidelity']:.6f}  "
                  f"F_cphase={metrics['cphase_sqr_fidelity']:.6f}  "
                  f"F_target={metrics['target_fidelity']:.6f}  "
                  f"leak={leak_max:.2e}")

    elapsed = time.time() - t_start
    print(f"\nScan complete in {elapsed:.1f} s")

    # Save results
    data_dir = SCRIPT_DIR.parent / "data"
    np.savez(
        data_dir / "phase1_phase2_results.npz",
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
    )
    print(f"Results saved to {data_dir / 'phase1_phase2_results.npz'}")

    return {
        "chi_t_values": CHI_T_VALUES,
        "family_names": family_names,
        "true_sqr_fidelity": true_sqr_fid,
        "cphase_sqr_fidelity": cphase_sqr_fid,
        "target_fidelity": target_fid,
        "target_z_phase": target_z_phase,
        "spectator_phase_spread": spec_phase_spread,
        "spectator_max_transverse": spec_max_transverse,
        "branch_fidelities_true": branch_fid_true_all,
        "branch_fidelities_cphase": branch_fid_cphase_all,
        "spectator_phases": spec_phases_all,
    }


if __name__ == "__main__":
    run_scan()
