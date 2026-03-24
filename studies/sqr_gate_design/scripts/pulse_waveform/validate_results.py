"""
Validation script for the SQR pulse-waveform design study.

Runs the AGENTS.md three-check validation gate:
1. Sanity checks (zero drive, unitarity, symmetry)
2. Convergence tests (Hilbert space dimension, time step)
3. Literature comparison (N/A for this study — no REP class)

Usage:
    python scripts/validate_results.py
"""

import sys
from pathlib import Path
import numpy as np


def _configure_utf8_output() -> None:
    for stream in (sys.stdout, sys.stderr):
        reconfigure = getattr(stream, "reconfigure", None)
        if reconfigure is not None:
            try:
                reconfigure(encoding="utf-8")
            except OSError:
                pass


_configure_utf8_output()

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR.parents[1] / "data" / SCRIPT_DIR.name
sys.path.insert(0, str(SCRIPT_DIR))

import runtime_compat  # noqa: F401

from common import (
    build_model, build_frame, duration_from_chi_t, extract_branch_unitaries,
    z_corrected_target_fidelity, spectator_z_fidelity, target_qubit_unitary,
    N_FOCK, TARGET_N0, THETA_TARGET, PHI_TARGET, DT, N_CAV, N_TR,
)
from run_phase1_phase2 import (
    build_single_tone_gaussian, build_cosine_squared_pulse, simulate_waveform,
    compute_metrics,
)

PASS = "\033[92mPASS\033[0m"
FAIL = "\033[91mFAIL\033[0m"


def check_sanity():
    """Check 1: Sanity checks."""
    print("=" * 60)
    print("CHECK 1: SANITY CHECKS")
    print("=" * 60)
    all_pass = True
    model = build_model()
    frame = build_frame(model)
    T = duration_from_chi_t(10)

    # (a) Zero drive → no population transfer (diagonal propagator)
    # Note: In the rotating frame, the qubit accumulates Z-phase from the
    # dispersive interaction even without drive. We check for zero POPULATION
    # TRANSFER (off-diagonal elements), not identity.
    print("\n(a) Zero drive → no population transfer:")
    from cqed_sim.pulses import Pulse
    from cqed_sim.pulses.envelopes import normalized_gaussian
    pulse = Pulse(
        channel="q", t0=0.0, duration=T,
        envelope=lambda t_rel: normalized_gaussian(t_rel, sigma_fraction=1.0/6.0),
        carrier=0.0, phase=0.0, amp=0.0,  # zero amplitude
    )
    _, unitaries = simulate_waveform(model, frame, [pulse], {"q": "qubit"}, T)
    for n in range(N_FOCK):
        off_diag = abs(unitaries[n][0, 1])**2 + abs(unitaries[n][1, 0])**2
        status = PASS if off_diag < 1e-6 else FAIL
        if off_diag >= 1e-6:
            all_pass = False
        print(f"  Branch n={n}: off-diag pop = {off_diag:.2e}  [{status}]")

    # (b) Unitarity of branch operators
    print("\n(b) Unitarity of branch operators (|det| = 1):")
    pulses, drive_ops = build_single_tone_gaussian(
        model, frame, TARGET_N0, THETA_TARGET, PHI_TARGET, T,
    )
    _, unitaries = simulate_waveform(model, frame, pulses, drive_ops, T)
    for n in range(N_FOCK):
        det = np.linalg.det(unitaries[n])
        unitarity_err = abs(abs(det) - 1.0)
        status = PASS if unitarity_err < 1e-6 else FAIL
        if unitarity_err >= 1e-6:
            all_pass = False
        print(f"  Branch n={n}: |det(U_n)| = {abs(det):.10f}, err = {unitarity_err:.2e}  [{status}]")

    # (c) Target branch rotation angle is π
    print("\n(c) Target branch rotation angle = π (full population inversion):")
    U_target_branch = unitaries[TARGET_N0]
    pop_inversion = abs(U_target_branch[0, 0])**2
    status = PASS if pop_inversion < 1e-6 else FAIL
    if pop_inversion >= 1e-6:
        all_pass = False
    print(f"  |U[0,0]|² = {pop_inversion:.10f} (should be ~0)  [{status}]")
    off_diag_mag = abs(U_target_branch[0, 1])**2
    status2 = PASS if abs(off_diag_mag - 1.0) < 1e-6 else FAIL
    if abs(off_diag_mag - 1.0) >= 1e-6:
        all_pass = False
    print(f"  |U[0,1]|² = {off_diag_mag:.10f} (should be ~1)  [{status2}]")

    # (d) Gaussian and multitone give identical results
    print("\n(d) Gaussian ≡ multitone (single active tone):")
    from run_phase1_phase2 import build_multitone_pulse
    pulses_mt, drive_ops_mt = build_multitone_pulse(
        model, frame, N_FOCK, TARGET_N0, THETA_TARGET, PHI_TARGET, T,
    )
    _, unitaries_mt = simulate_waveform(model, frame, pulses_mt, drive_ops_mt, T)
    max_diff = 0
    for n in range(N_FOCK):
        diff = np.max(np.abs(unitaries[n] - unitaries_mt[n]))
        max_diff = max(max_diff, diff)
    status = PASS if max_diff < 1e-4 else FAIL
    if max_diff >= 1e-4:
        all_pass = False
    print(f"  max|U_gauss - U_multi| = {max_diff:.2e}  [{status}]")

    return all_pass


def check_convergence():
    """Check 2: Convergence with respect to numerical parameters."""
    print("\n" + "=" * 60)
    print("CHECK 2: CONVERGENCE")
    print("=" * 60)
    all_pass = True

    chi_t = 10
    T = duration_from_chi_t(chi_t)

    # (a) Hilbert space dimension: N_CAV = 6 vs 8
    print("\n(a) Cavity dimension convergence (N_cav = 6 vs 8):")
    model_6 = build_model(n_cav=6)
    model_8 = build_model(n_cav=8)
    frame_6 = build_frame(model_6)
    frame_8 = build_frame(model_8)

    pulses_6, ops_6 = build_single_tone_gaussian(
        model_6, frame_6, TARGET_N0, THETA_TARGET, PHI_TARGET, T,
    )
    pulses_8, ops_8 = build_single_tone_gaussian(
        model_8, frame_8, TARGET_N0, THETA_TARGET, PHI_TARGET, T,
    )
    _, U6 = simulate_waveform(model_6, frame_6, pulses_6, ops_6, T)
    metrics_6 = compute_metrics(U6, TARGET_N0, THETA_TARGET, PHI_TARGET)

    # Need to extract using n_fock=4 from the 8-dim model
    from cqed_sim.sequence import SequenceCompiler
    from cqed_sim.sim import SimulationConfig, simulate_sequence, prepare_simulation
    compiler = SequenceCompiler(dt=DT)
    compiled = compiler.compile(pulses_8, t_end=T + 4 * DT)
    config = SimulationConfig(frame=frame_8, store_states=False)
    session = prepare_simulation(model_8, compiled, ops_8, config=config, e_ops={})
    final_states_8 = []
    for n in range(N_FOCK):
        for q in [0, 1]:
            psi0 = model_8.basis_state(q, n)
            result = session.run(psi0)
            final_states_8.append(result.final_state)
    U8 = extract_branch_unitaries(final_states_8, model_8, N_FOCK)
    metrics_8 = compute_metrics(U8, TARGET_N0, THETA_TARGET, PHI_TARGET)

    delta_cphase = abs(metrics_8["cphase_sqr_fidelity"] - metrics_6["cphase_sqr_fidelity"])
    status = PASS if delta_cphase < 1e-4 else FAIL
    if delta_cphase >= 1e-4:
        all_pass = False
    print(f"  F_cphase(N=6) = {metrics_6['cphase_sqr_fidelity']:.8f}")
    print(f"  F_cphase(N=8) = {metrics_8['cphase_sqr_fidelity']:.8f}")
    print(f"  ΔF = {delta_cphase:.2e}  [{status}]")

    # (b) Time step convergence: dt = 2ns vs 1ns
    print("\n(b) Time step convergence (dt = 2 ns vs 1 ns):")
    model = build_model()
    frame = build_frame(model)
    pulses, ops = build_cosine_squared_pulse(
        model, frame, TARGET_N0, THETA_TARGET, PHI_TARGET, T,
    )

    # dt = 2ns (baseline)
    compiler_2 = SequenceCompiler(dt=2e-9)
    compiled_2 = compiler_2.compile(pulses, t_end=T + 4 * 2e-9)
    config_2 = SimulationConfig(frame=frame, store_states=False)
    session_2 = prepare_simulation(model, compiled_2, ops, config=config_2, e_ops={})
    fs_2 = []
    for n in range(N_FOCK):
        for q in [0, 1]:
            fs_2.append(session_2.run(model.basis_state(q, n)).final_state)
    U_dt2 = extract_branch_unitaries(fs_2, model, N_FOCK)
    m_dt2 = compute_metrics(U_dt2, TARGET_N0, THETA_TARGET, PHI_TARGET)

    # dt = 1ns
    compiler_1 = SequenceCompiler(dt=1e-9)
    compiled_1 = compiler_1.compile(pulses, t_end=T + 4 * 1e-9)
    config_1 = SimulationConfig(frame=frame, store_states=False)
    session_1 = prepare_simulation(model, compiled_1, ops, config=config_1, e_ops={})
    fs_1 = []
    for n in range(N_FOCK):
        for q in [0, 1]:
            fs_1.append(session_1.run(model.basis_state(q, n)).final_state)
    U_dt1 = extract_branch_unitaries(fs_1, model, N_FOCK)
    m_dt1 = compute_metrics(U_dt1, TARGET_N0, THETA_TARGET, PHI_TARGET)

    delta_dt = abs(m_dt1["cphase_sqr_fidelity"] - m_dt2["cphase_sqr_fidelity"])
    status = PASS if delta_dt < 1e-4 else FAIL
    if delta_dt >= 1e-4:
        all_pass = False
    print(f"  F_cphase(dt=2ns) = {m_dt2['cphase_sqr_fidelity']:.8f}")
    print(f"  F_cphase(dt=1ns) = {m_dt1['cphase_sqr_fidelity']:.8f}")
    print(f"  ΔF = {delta_dt:.2e}  [{status}]")

    return all_pass


def check_literature():
    """Check 3: Literature comparison."""
    print("\n" + "=" * 60)
    print("CHECK 3: LITERATURE COMPARISON")
    print("=" * 60)
    print("\n  This is an OPT/ANA/DES class study — no published benchmarks")
    print("  to reproduce.  Key internal consistency checks:")
    print("  - Gaussian and multitone produce identical results (checked above)")
    print("  - F_cphase > F_true at all χT (confirmed by scan)")
    print("  - Spectator transverse error decreases monotonically and is strongly suppressed at large χT")

    # Verify scaling of spectator transverse error
    data = np.load(DATA_DIR / "phase1_phase2_results.npz")
    chi_t = data["chi_t_values"]
    # Use single-tone Gaussian (family 0)
    spec_max_trans = data["spectator_max_transverse"][0]

    diffs = np.diff(spec_max_trans)
    monotone = bool(np.all(diffs <= 1.0e-9))
    tail_mask = chi_t >= 5
    tail_max = float(np.max(spec_max_trans[tail_mask]))

    monotone_status = PASS if monotone else FAIL
    tail_status = PASS if tail_max < 1.0e-6 else FAIL

    print(f"\n  Max spectator transverse error for χT/(2π) ≥ 5: {tail_max:.2e}  [{tail_status}]")
    print(f"  Spectator transverse error monotone with χT: {monotone_status}")

    return monotone and (tail_max < 1.0e-6)


def check_phase4():
    """Check Phase 4: higher-order corrections produce small fidelity shift."""
    print("\n" + "=" * 60)
    print("CHECK 4: PHASE 4 — HIGHER-ORDER CORRECTIONS")
    print("=" * 60)
    all_pass = True

    p12 = np.load(DATA_DIR / "phase1_phase2_results.npz")
    p4 = np.load(DATA_DIR / "phase4_results.npz")

    # At χT=10, cphase fidelity should change by <1% due to χ', K
    chi_t = p12["chi_t_values"]
    idx_10 = np.argmin(np.abs(chi_t - 10))

    print("\n  cphase SQR fidelity at χT=10 (Gaussian):")
    f_p12 = p12["cphase_sqr_fidelity"][0, idx_10]
    f_p4 = p4["cphase_sqr_fidelity"][0, idx_10]
    delta = abs(f_p12 - f_p4)
    status = PASS if delta < 0.01 else FAIL
    if delta >= 0.01:
        all_pass = False
    print(f"  Phase 1-2: {f_p12:.6f}")
    print(f"  Phase 4:   {f_p4:.6f}")
    print(f"  ΔF = {delta:.2e}  [{status}]")

    # Branch frequencies shift: n=2 should shift by ~42 kHz (χ' = -21 kHz)
    from common import build_model, build_frame, branch_frequencies, CHI_PRIME, KERR
    m_ideal = build_model()
    m_p4 = build_model(chi_prime=CHI_PRIME, kerr=KERR)
    f_ideal = build_frame(m_ideal)
    f_p4f = build_frame(m_p4)
    freqs_ideal = branch_frequencies(m_ideal, f_ideal, N_FOCK)
    freqs_p4 = branch_frequencies(m_p4, f_p4f, N_FOCK)
    delta_n2 = (freqs_p4[2] - freqs_ideal[2]) / (2 * np.pi)
    print(f"\n  Branch n=2 frequency shift: {delta_n2/1e3:.1f} kHz")
    # Expected: ~42 kHz = 2 × chi_prime_hz
    expected_shift = 2 * CHI_PRIME / (2 * np.pi)
    status = PASS if abs(delta_n2 - expected_shift) < 1e3 else FAIL
    if abs(delta_n2 - expected_shift) >= 1e3:
        all_pass = False
    print(f"  Expected ~ {expected_shift/1e3:.1f} kHz  [{status}]")

    return all_pass


def check_phase5():
    """Check Phase 5: decoherence fidelity scales as ~1 - T/(2T1)."""
    print("\n" + "=" * 60)
    print("CHECK 5: PHASE 5 — OPEN-SYSTEM DECOHERENCE")
    print("=" * 60)
    all_pass = True

    p5 = np.load(DATA_DIR / "phase5_results.npz")
    chi_t = p5["chi_t_values"]
    chi_rad = float(p5["chi_rad_s"])
    t1_s = float(p5["t1_s"])

    # Decoherence fidelity should be envelope-independent (common T/T1 limit)
    idx_max = np.argmax(chi_t)  # use largest available χT/(2π)
    chi_t_ref = chi_t[idx_max]
    print(f"\n  Decoherence fidelity envelope-independence (chiT/(2pi)={chi_t_ref:.0f}):")
    f_deco_all = p5["deco_fid"][:, idx_max]
    spread = np.ptp(f_deco_all)
    status = PASS if spread < 5e-3 else FAIL
    if spread >= 5e-3:
        all_pass = False
    print(f"  F_deco range: [{f_deco_all.min():.6f}, {f_deco_all.max():.6f}]")
    print(f"  Spread = {spread:.2e}  [{status}]")

    # Scaling: F_deco ≈ 1 - T/(2T1) at large χT
    # chi_t is χT/(2π), so T = chi_t / f_chi where f_chi = |χ|/(2π)
    f_chi = abs(chi_rad) / (2 * np.pi)
    print(f"\n  Decoherence scaling (Gaussian, chiT/(2pi)={chi_t_ref:.0f}):")
    T_s = chi_t[idx_max] / f_chi
    F_measured = float(p5["deco_fid"][0, idx_max])
    F_analytic = 1 - T_s / (2 * t1_s)
    delta = abs(F_measured - F_analytic)
    status = PASS if delta < 0.02 else FAIL
    if delta >= 0.02:
        all_pass = False
    print(f"  F_deco (measured)  = {F_measured:.6f}")
    print(f"  F_deco (1-T/2T1) = {F_analytic:.6f}")
    print(f"  T = {T_s*1e6:.2f} us, T/(2T1) = {T_s/(2*t1_s):.4f}")
    print(f"  dF = {delta:.4f}  [{status}]")

    # Net fidelity should have a peak (optimal χT)
    print("\n  Net fidelity peak detection (Gaussian):")
    f_net = p5["cphase_fid_net"][0]
    idx_peak = np.argmax(f_net)
    chi_t_peak = chi_t[idx_peak]
    f_peak = f_net[idx_peak]
    # Peak should be in range χT/(2π) = 1–10
    status = PASS if 1 <= chi_t_peak <= 10 else FAIL
    if not (1 <= chi_t_peak <= 10):
        all_pass = False
    print(f"  Peak F_net = {f_peak:.6f} at chiT/(2pi) = {chi_t_peak}")
    print(f"  Expected range: chiT/(2pi) = 1-10  [{status}]")

    return all_pass


def main():
    print("SQR Pulse-Waveform Design Study — Validation (Extended)")
    print("=" * 60)

    s1 = check_sanity()
    s2 = check_convergence()
    s3 = check_literature()
    s4 = check_phase4()
    s5 = check_phase5()

    print("\n" + "=" * 60)
    print("VALIDATION SUMMARY")
    print("=" * 60)
    print(f"  Sanity checks:     {PASS if s1 else FAIL}")
    print(f"  Convergence:       {PASS if s2 else FAIL}")
    print(f"  Literature/scaling: {PASS if s3 else FAIL}")
    print(f"  Phase 4 (χ', K):   {PASS if s4 else FAIL}")
    print(f"  Phase 5 (T1, T2):  {PASS if s5 else FAIL}")

    if s1 and s2 and s3 and s4 and s5:
        print(f"\n  ALL CHECKS PASSED — ready for report")
    else:
        print(f"\n  SOME CHECKS FAILED — review before reporting")


if __name__ == "__main__":
    main()
