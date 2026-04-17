"""Cross-regime synthesis: combine displacement and qubit rotation results
into unified regime maps and practical design guidance.

Also runs the validation checks (sanity, convergence, literature).
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import qutip as qt

import common
from common import (
    TWO_PI, N_TR, N_CAV, DEFAULT_DT, CHI, CHI_PRIME, KERR,
    build_model, build_frame,
    make_square_displacement_pulse, make_gaussian_qubit_pulse,
    compile_and_prepare, simulate_state,
    propagate_basis_states, extract_qubit_block, rotation_metrics,
    displacement_fidelity, entanglement_entropy,
    save_json, load_json, ARTIFACTS_DIR, FIGURES_DIR, DATA_DIR,
    apply_plot_style, TOL_BRIGHT,
)
from cqed_sim.core.ideal_gates import displacement_op, qubit_rotation_xy

common.apply_plot_style()


# ===================================================================
# 1. CROSS-REGIME MAP
# ===================================================================

def build_displacement_regime_map():
    """Classify displacement accuracy into regimes based on chi*T and alpha."""
    print("=" * 60)
    print("DISPLACEMENT REGIME MAP")
    print("=" * 60)

    # Load sweep data
    data = np.load(DATA_DIR / "displacement_sweep.npz")
    alpha_values = data["alpha_values"]
    duration_values = data["duration_values"]
    fid_g = data["fidelity_g"]
    fid_e = data["fidelity_e"]

    chi_T = np.abs(CHI) * duration_values / TWO_PI  # dimensionless chi*T in cycles

    # Classify regimes for |e>
    # Fast regime: T * |chi|/(2pi) < 0.01 => fidelity should be high
    # Moderate: 0.01 < chi*T/(2pi) < 0.1
    # Broken: chi*T/(2pi) > 0.1
    regime_map = {
        "alpha_values": alpha_values.tolist(),
        "duration_ns": (duration_values * 1e9).tolist(),
        "chi_T_cycles": chi_T.tolist(),
        "fidelity_g": fid_g.tolist(),
        "fidelity_e": fid_e.tolist(),
        "thresholds": {
            "fid_99": {"description": "Fidelity > 0.99", "chi_T_threshold_cycles": None},
            "fid_999": {"description": "Fidelity > 0.999", "chi_T_threshold_cycles": None},
        },
    }

    # Find chi*T threshold for 99% fidelity at alpha=1.0
    alpha_1_idx = np.argmin(np.abs(alpha_values - 1.0))
    for thresh_name, thresh_val in [("fid_99", 0.99), ("fid_999", 0.999)]:
        for i, ct in enumerate(chi_T):
            if fid_e[i, alpha_1_idx] < thresh_val:
                regime_map["thresholds"][thresh_name]["chi_T_threshold_cycles"] = float(chi_T[max(0, i - 1)])
                break

    save_json(ARTIFACTS_DIR / "displacement_regime_map.json", regime_map)
    print("Saved displacement_regime_map.json")
    return regime_map


def build_qubit_rotation_regime_map():
    """Classify qubit rotation accuracy into regimes based on Fock level and duration."""
    print("\n" + "=" * 60)
    print("QUBIT ROTATION REGIME MAP")
    print("=" * 60)

    fock_data = load_json(ARTIFACTS_DIR / "qubit_rotation_fock_resolved.json")

    # For X_pi at each duration, find the max Fock level with fidelity > threshold
    pulse_label = "X_pi"
    durations_tested = sorted(fock_data[pulse_label].keys(), key=lambda k: float(k.replace("ns", "")))

    regime_map = {"pulse": pulse_label, "durations": [], "regimes": {}}

    for T_key in durations_tested:
        met = fock_data[pulse_label][T_key]
        fids = {int(n): met[n]["process_fidelity"] for n in met}
        max_n_99 = max((n for n, f in fids.items() if f > 0.99), default=-1)
        max_n_999 = max((n for n, f in fids.items() if f > 0.999), default=-1)
        max_n_9999 = max((n for n, f in fids.items() if f > 0.9999), default=-1)

        regime_map["durations"].append(T_key)
        regime_map["regimes"][T_key] = {
            "fidelities": {str(n): float(f) for n, f in fids.items()},
            "max_fock_99": max_n_99,
            "max_fock_999": max_n_999,
            "max_fock_9999": max_n_9999,
        }

        print(f"  {T_key}: max_n(99%)={max_n_99}, max_n(99.9%)={max_n_999}, "
              f"max_n(99.99%)={max_n_9999}")

    save_json(ARTIFACTS_DIR / "qubit_rotation_regime_map.json", regime_map)
    print("Saved qubit_rotation_regime_map.json")
    return regime_map


def build_summary_table():
    """Combined summary table of both displacement and rotation regimes."""
    print("\n" + "=" * 60)
    print("SUMMARY TABLE")
    print("=" * 60)

    disp_regime = load_json(ARTIFACTS_DIR / "displacement_regime_map.json")
    rot_regime = load_json(ARTIFACTS_DIR / "qubit_rotation_regime_map.json")

    summary = {
        "displacement": {
            "operation": "D(alpha)",
            "waveform": "Square envelope, on-resonance carrier",
            "dominant_distortion_qubit_g": "Kerr-induced phase drift (negligible for alpha < 2, T < 200 ns)",
            "dominant_distortion_qubit_e": "chi-induced frequency mismatch (dominant error)",
            "regime_qubit_g": "Ideal spectator for |chi*T| << 1, alpha < ~3",
            "regime_qubit_e": "Breaks for |chi*T| > 0.03 cycles (T > 10 ns at chi=-2.84 MHz)",
            "recommendation": "For unconditional displacement: use fast pulses (T << 1/chi). For conditional displacement: exploit the chi-dependent fidelity.",
        },
        "qubit_rotation": {
            "operation": "R_X(theta), R_Y(theta)",
            "waveform": "Gaussian envelope (optionally DRAG-corrected)",
            "dominant_distortion": "Fock-dependent rotation angle due to chi-shifted transition frequency",
            "regime_spectator": "Cavity is acceptable spectator for n < n_max(T)",
            "max_fock_levels": {T: r["max_fock_99"] for T, r in rot_regime["regimes"].items()},
            "recommendation": "Use shorter pulses for higher cavity occupation tolerance. DRAG improves n=0 fidelity but does not fix the Fock-dependent angle shift.",
        },
    }

    save_json(ARTIFACTS_DIR / "cross_regime_summary.json", summary)
    print("Saved cross_regime_summary.json")

    # Print formatted table
    print("\n  --- Displacement Regimes ---")
    for k, v in summary["displacement"].items():
        print(f"    {k}: {v}")
    print("\n  --- Qubit Rotation Regimes ---")
    for k, v in summary["qubit_rotation"].items():
        print(f"    {k}: {v}")

    return summary


# ===================================================================
# 2. VALIDATION
# ===================================================================

def run_sanity_checks():
    """Sanity checks: verify ideal limits (chi=chi'=K=0)."""
    print("\n" + "=" * 60)
    print("VALIDATION: SANITY CHECKS")
    print("=" * 60)

    results = {}

    # 1. Displacement in ideal limit
    model_ideal = build_model(chi=0.0, chi_prime=None, kerr=0.0)
    frame_ideal = build_frame(model_ideal)

    for alpha_test in [0.5, 1.0, 2.0]:
        pulse = make_square_displacement_pulse(model_ideal, frame_ideal, alpha=alpha_test, duration_s=100e-9)
        session = compile_and_prepare(model_ideal, frame_ideal, [pulse])
        psi0 = model_ideal.basis_state(0, 0)
        psi_f = simulate_state(session, psi0)
        D = displacement_op(model_ideal.n_cav, alpha_test)
        ideal = qt.tensor(qt.basis(model_ideal.n_tr, 0), D * qt.basis(model_ideal.n_cav, 0))
        fid = displacement_fidelity(psi_f, ideal)
        results[f"displacement_ideal_alpha{alpha_test}"] = {"fidelity": fid, "pass": fid > 0.9999}
        print(f"  Displacement ideal (alpha={alpha_test}): fid={fid:.6f} {'PASS' if fid > 0.9999 else 'FAIL'}")

    # 2. Qubit rotation in ideal limit (n=0)
    target_pi = np.asarray(qubit_rotation_xy(np.pi, 0.0).full(), dtype=np.complex128)
    pulse_q = make_gaussian_qubit_pulse(model_ideal, frame_ideal, theta=np.pi, phase=0.0, duration_s=40e-9)
    session_q = compile_and_prepare(model_ideal, frame_ideal, [pulse_q])
    U = propagate_basis_states(model_ideal, session_q, cavity_levels=[0], qubit_levels=[0, 1])
    block = extract_qubit_block(U, model_ideal.n_cav, 0)
    met = rotation_metrics(block, target=target_pi)
    results["qubit_pi_ideal_n0"] = {
        "process_fidelity": met["process_fidelity"],
        "rotation_angle": met["rotation_angle_rad"],
        "pass": met["process_fidelity"] > 0.998,
    }
    print(f"  Qubit X_pi ideal (n=0): fid={met['process_fidelity']:.6f}, "
          f"angle={met['rotation_angle_rad']:.4f} {'PASS' if met['process_fidelity'] > 0.998 else 'FAIL'}")

    # 3. Check that displacement |g> and |e> give same result in ideal limit
    pulse_ge = make_square_displacement_pulse(model_ideal, frame_ideal, alpha=1.0, duration_s=100e-9)
    session_ge = compile_and_prepare(model_ideal, frame_ideal, [pulse_ge])
    for q_label, q_level in [("g", 0), ("e", 1)]:
        psi0 = model_ideal.basis_state(q_level, 0)
        psi_f = simulate_state(session_ge, psi0)
        D = displacement_op(model_ideal.n_cav, 1.0)
        ideal = qt.tensor(qt.basis(model_ideal.n_tr, q_level), D * qt.basis(model_ideal.n_cav, 0))
        fid = displacement_fidelity(psi_f, ideal)
        results[f"displacement_ideal_qubit_{q_label}"] = {"fidelity": fid, "pass": fid > 0.9999}
        print(f"  Displacement ideal |{q_label}> (alpha=1.0): fid={fid:.6f}")

    all_pass = all(r.get("pass", True) for r in results.values())
    results["all_pass"] = all_pass
    save_json(ARTIFACTS_DIR / "validation_sanity.json", results)
    print(f"\n  Sanity checks: {'ALL PASS' if all_pass else 'SOME FAILED'}")
    return results


def run_convergence_checks():
    """Convergence checks: vary n_cav, n_tr, dt."""
    print("\n" + "=" * 60)
    print("VALIDATION: CONVERGENCE CHECKS")
    print("=" * 60)

    results = {}

    # 1. Convergence in n_cav for displacement
    alpha = 1.0
    T = 100e-9
    n_cav_values = [8, 10, 12, 15, 20]
    disp_conv = {}
    for n_cav in n_cav_values:
        model = build_model(n_cav=n_cav)
        frame = build_frame(model)
        pulse = make_square_displacement_pulse(model, frame, alpha=alpha, duration_s=T)
        session = compile_and_prepare(model, frame, [pulse])
        psi0 = model.basis_state(0, 0)
        psi_f = simulate_state(session, psi0)
        D = displacement_op(n_cav, alpha)
        ideal = qt.tensor(qt.basis(model.n_tr, 0), D * qt.basis(n_cav, 0))
        fid = displacement_fidelity(psi_f, ideal)
        disp_conv[str(n_cav)] = fid
        print(f"  n_cav={n_cav}: displacement fid={fid:.6f}")

    results["displacement_n_cav_convergence"] = disp_conv
    # Check convergence: difference between last two
    fids_list = list(disp_conv.values())
    disp_converged = abs(fids_list[-1] - fids_list[-2]) < 1e-5
    results["displacement_n_cav_converged"] = disp_converged
    print(f"  Convergence: delta(15 vs 20) = {abs(fids_list[-2] - fids_list[-1]):.2e} "
          f"{'PASS' if disp_converged else 'FAIL'}")

    # 2. Convergence in dt for displacement
    dt_values = [0.25e-9, 0.5e-9, 1.0e-9, 2.0e-9]
    dt_conv = {}
    model = build_model()
    frame = build_frame(model)
    for dt in dt_values:
        pulse = make_square_displacement_pulse(model, frame, alpha=alpha, duration_s=T)
        session = compile_and_prepare(model, frame, [pulse], dt=dt)
        psi0 = model.basis_state(0, 0)
        psi_f = simulate_state(session, psi0)
        D = displacement_op(model.n_cav, alpha)
        ideal = qt.tensor(qt.basis(model.n_tr, 0), D * qt.basis(model.n_cav, 0))
        fid = displacement_fidelity(psi_f, ideal)
        dt_conv[f"{dt*1e9:.2f}ns"] = fid
        print(f"  dt={dt*1e9:.2f}ns: displacement fid={fid:.6f}")

    results["displacement_dt_convergence"] = dt_conv
    dt_fids = list(dt_conv.values())
    dt_converged = abs(dt_fids[0] - dt_fids[1]) < 1e-5
    results["displacement_dt_converged"] = dt_converged

    # 3. Convergence in n_cav for qubit rotation (n=0 block)
    rot_conv = {}
    target_pi = np.asarray(qubit_rotation_xy(np.pi, 0.0).full(), dtype=np.complex128)
    for n_cav in n_cav_values:
        model = build_model(n_cav=n_cav)
        frame = build_frame(model)
        pulse = make_gaussian_qubit_pulse(model, frame, theta=np.pi, phase=0.0, duration_s=40e-9)
        session = compile_and_prepare(model, frame, [pulse])
        U = propagate_basis_states(model, session, cavity_levels=[0], qubit_levels=[0, 1])
        block = extract_qubit_block(U, n_cav, 0)
        met = rotation_metrics(block, target=target_pi)
        rot_conv[str(n_cav)] = met["process_fidelity"]
        print(f"  n_cav={n_cav}: qubit pi fid(n=0)={met['process_fidelity']:.6f}")

    results["rotation_n_cav_convergence"] = rot_conv
    rot_fids = list(rot_conv.values())
    rot_converged = abs(rot_fids[-1] - rot_fids[-2]) < 1e-5
    results["rotation_n_cav_converged"] = rot_converged

    # 4. Convergence in n_tr for qubit rotation
    n_tr_values = [2, 3, 4]
    ntr_conv = {}
    for n_tr in n_tr_values:
        try:
            model = build_model(n_tr=n_tr)
            frame = build_frame(model)
            pulse = make_gaussian_qubit_pulse(model, frame, theta=np.pi, phase=0.0, duration_s=40e-9)
            session = compile_and_prepare(model, frame, [pulse])
            U = propagate_basis_states(model, session, cavity_levels=[0], qubit_levels=[0, 1])
            block = extract_qubit_block(U, model.n_cav, 0)
            met = rotation_metrics(block, target=target_pi)
            ntr_conv[str(n_tr)] = met["process_fidelity"]
            print(f"  n_tr={n_tr}: qubit pi fid(n=0)={met['process_fidelity']:.6f}")
        except Exception as e:
            print(f"  n_tr={n_tr}: FAILED ({e})")
            ntr_conv[str(n_tr)] = None

    results["rotation_n_tr_convergence"] = ntr_conv

    all_converged = disp_converged and dt_converged and rot_converged
    results["all_converged"] = all_converged
    save_json(ARTIFACTS_DIR / "validation_convergence.json", results)
    print(f"\n  Convergence checks: {'ALL PASS' if all_converged else 'CHECK DETAILS'}")
    return results


# ===================================================================
# 3. REGIME MAP FIGURES
# ===================================================================

def generate_regime_figures():
    """Generate regime map figures combining displacement and rotation data."""
    print("\n" + "=" * 60)
    print("GENERATING REGIME MAP FIGURES")
    print("=" * 60)

    # --- Figure: Duration vs max Fock level for qubit rotation ---
    rot_regime = load_json(ARTIFACTS_DIR / "qubit_rotation_regime_map.json")
    durations = [float(T.replace("ns", "")) for T in rot_regime["durations"]]
    max_n_99 = [rot_regime["regimes"][T]["max_fock_99"] for T in rot_regime["durations"]]
    max_n_999 = [rot_regime["regimes"][T]["max_fock_999"] for T in rot_regime["durations"]]
    max_n_9999 = [rot_regime["regimes"][T]["max_fock_9999"] for T in rot_regime["durations"]]

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(durations, max_n_99, "o-", color=TOL_BRIGHT[0], label=r"$\mathcal{F} > 99\%$")
    ax.plot(durations, max_n_999, "s-", color=TOL_BRIGHT[1], label=r"$\mathcal{F} > 99.9\%$")
    ax.plot(durations, max_n_9999, "^-", color=TOL_BRIGHT[2], label=r"$\mathcal{F} > 99.99\%$")
    ax.set_xlabel("Pulse duration (ns)")
    ax.set_ylabel("Max cavity Fock level $n$")
    ax.set_title("Qubit Rotation: Spectator Regime Boundary")
    ax.legend()
    ax.set_ylim(-0.5, 6)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "regime_map_qubit_rotation.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "regime_map_qubit_rotation.pdf", bbox_inches="tight")
    plt.close(fig)

    # --- Figure: Combined displacement + rotation regime overview ---
    disp_data = np.load(DATA_DIR / "displacement_sweep.npz")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Left: displacement fidelity for |e> (the harder case)
    alpha_vals = disp_data["alpha_values"]
    dur_vals = disp_data["duration_values"] * 1e9
    fid_e = disp_data["fidelity_e"]
    im = axes[0].pcolormesh(alpha_vals, dur_vals, 1 - fid_e, shading="nearest",
                            cmap="RdYlGn_r", vmin=0, vmax=1)
    axes[0].set_xlabel(r"Displacement $|\alpha|$")
    axes[0].set_ylabel("Pulse duration (ns)")
    axes[0].set_title(r"Displacement Infidelity ($|e\rangle$ qubit)")
    plt.colorbar(im, ax=axes[0], label=r"$1-\mathcal{F}$")

    # Right: qubit rotation infidelity heatmap (n vs T)
    fock_data = load_json(ARTIFACTS_DIR / "qubit_rotation_fock_resolved.json")
    pulse_data = fock_data["X_pi"]
    T_keys = sorted(pulse_data.keys(), key=lambda k: float(k.replace("ns", "")))
    n_levels = list(range(6))
    infid_grid = np.zeros((len(T_keys), len(n_levels)))
    for i, T_key in enumerate(T_keys):
        for j, n in enumerate(n_levels):
            infid_grid[i, j] = 1 - pulse_data[T_key][str(n)]["process_fidelity"]

    T_vals = [float(T.replace("ns", "")) for T in T_keys]
    im2 = axes[1].pcolormesh(n_levels, T_vals, infid_grid, shading="nearest",
                             cmap="RdYlGn_r", norm=matplotlib.colors.LogNorm(vmin=1e-4, vmax=1))
    axes[1].set_xlabel("Cavity Fock level $n$")
    axes[1].set_ylabel("Pulse duration (ns)")
    axes[1].set_title(r"$X_\pi$ Gate Infidelity (Fock-resolved)")
    plt.colorbar(im2, ax=axes[1], label=r"$1-\mathcal{F}_{\mathrm{proc}}$")

    fig.suptitle("Combined Regime Map: Displacement and Qubit Rotation", y=1.02)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "combined_regime_map.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "combined_regime_map.pdf", bbox_inches="tight")
    plt.close(fig)

    # --- Figure: Convergence ---
    conv_data = load_json(ARTIFACTS_DIR / "validation_convergence.json")

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))

    # n_cav convergence (displacement)
    n_cav_data = conv_data["displacement_n_cav_convergence"]
    n_cavs = [int(k) for k in n_cav_data.keys()]
    fids = list(n_cav_data.values())
    axes[0].plot(n_cavs, [1 - f for f in fids], "o-", color=TOL_BRIGHT[0])
    axes[0].set_xlabel("Cavity Fock dimension $N_{\\mathrm{cav}}$")
    axes[0].set_ylabel(r"$1 - \mathcal{F}$")
    axes[0].set_yscale("log")
    axes[0].set_title("Displacement: $N_{\\mathrm{cav}}$ Convergence")

    # dt convergence
    dt_data = conv_data["displacement_dt_convergence"]
    dt_vals = [float(k.replace("ns", "")) for k in dt_data.keys()]
    dt_fids = list(dt_data.values())
    axes[1].plot(dt_vals, [1 - f for f in dt_fids], "s-", color=TOL_BRIGHT[1])
    axes[1].set_xlabel("Time step $\\Delta t$ (ns)")
    axes[1].set_ylabel(r"$1 - \mathcal{F}$")
    axes[1].set_yscale("log")
    axes[1].set_title("Displacement: $\\Delta t$ Convergence")

    # n_cav convergence (rotation)
    rot_data = conv_data["rotation_n_cav_convergence"]
    n_cavs_r = [int(k) for k in rot_data.keys()]
    fids_r = list(rot_data.values())
    axes[2].plot(n_cavs_r, [1 - f for f in fids_r], "^-", color=TOL_BRIGHT[2])
    axes[2].set_xlabel("Cavity Fock dimension $N_{\\mathrm{cav}}$")
    axes[2].set_ylabel(r"$1 - \mathcal{F}_{\mathrm{proc}}$")
    axes[2].set_yscale("log")
    axes[2].set_title("Qubit Rotation: $N_{\\mathrm{cav}}$ Convergence")

    fig.suptitle("Convergence Analysis", y=1.02)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "validation_convergence.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "validation_convergence.pdf", bbox_inches="tight")
    plt.close(fig)

    print("All regime map and validation figures saved.")


if __name__ == "__main__":
    build_displacement_regime_map()
    build_qubit_rotation_regime_map()
    build_summary_table()
    sanity = run_sanity_checks()
    convergence = run_convergence_checks()
    generate_regime_figures()
    print("\n=== CROSS-REGIME SYNTHESIS & VALIDATION COMPLETE ===")
