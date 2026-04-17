"""Qubit rotation pulse study: characterize how Gaussian and DRAG-Gaussian
pulses behave when the cavity is populated in a dispersive qubit+cavity system.

Extracts Fock-resolved effective SU(2) blocks, rotation angles, axes,
conditional phases, leakage, and entanglement metrics.
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
    TWO_PI, N_TR, N_CAV, DEFAULT_DT, CHI, CHI_PRIME, KERR, SIGMA_FRACTION,
    build_model, build_frame,
    make_gaussian_qubit_pulse,
    compile_and_prepare, simulate_state,
    propagate_basis_states, extract_qubit_block,
    nearest_su2, su2_to_rotation, rotation_metrics,
    entanglement_entropy, qubit_purity,
    qubit_state_from_joint, bloch_vector,
    save_json, ARTIFACTS_DIR, FIGURES_DIR, DATA_DIR,
    apply_plot_style, TOL_BRIGHT, PAULI_X, PAULI_Y, PAULI_Z,
)
from cqed_sim.core.ideal_gates import qubit_rotation_xy

common.apply_plot_style()

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
# Pulse types: (theta, phase, label)
PULSE_TYPES = [
    (np.pi, 0.0, "X_pi"),
    (np.pi / 2, 0.0, "X_pi2"),
    (np.pi, np.pi / 2, "Y_pi"),
    (np.pi / 2, np.pi / 2, "Y_pi2"),
]

DURATION_VALUES = np.array([20, 30, 40, 60, 80, 100]) * 1e-9
DRAG_VALUES = np.array([0.0, 0.5e-9, 1.0e-9, 2.0e-9, 5.0e-9])
CAVITY_LEVELS_PROBE = list(range(6))  # Fock 0 through 5
N_CAV_STUDY = 10  # Smaller Hilbert space since we only probe low Fock states


def run_fock_resolved_rotation_study():
    """Extract Fock-resolved qubit rotation parameters for each pulse type."""
    print("=" * 60)
    print("FOCK-RESOLVED QUBIT ROTATION STUDY")
    print("=" * 60)
    t_start = time.time()

    model = build_model(n_cav=N_CAV_STUDY)
    frame = build_frame(model)

    all_results = {}

    for theta, phase, pulse_label in PULSE_TYPES:
        print(f"\n--- Pulse: {pulse_label} (theta={theta:.3f}, phi={phase:.3f}) ---")

        target_2x2 = np.asarray(qubit_rotation_xy(theta, phase).full(), dtype=np.complex128)

        pulse_results = {}
        for T in DURATION_VALUES:
            pulse = make_gaussian_qubit_pulse(
                model, frame, theta=theta, phase=phase, duration_s=T,
                manifold_level=0, drag=0.0,
            )
            session = compile_and_prepare(model, frame, [pulse])
            U = propagate_basis_states(
                model, session,
                cavity_levels=CAVITY_LEVELS_PROBE,
                qubit_levels=[0, 1],
            )

            fock_metrics = {}
            for n in CAVITY_LEVELS_PROBE:
                block = extract_qubit_block(U, model.n_cav, n)
                met = rotation_metrics(block, target=target_2x2)
                fock_metrics[str(n)] = met

            pulse_results[f"{T*1e9:.0f}ns"] = fock_metrics
            fids = [fock_metrics[str(n)]["process_fidelity"] for n in CAVITY_LEVELS_PROBE]
            print(f"  T={T*1e9:.0f}ns: fid(n=0)={fids[0]:.6f}, fid(n=5)={fids[5]:.6f}")

        all_results[pulse_label] = pulse_results

    elapsed = time.time() - t_start
    print(f"\nFock-resolved study completed in {elapsed:.1f} s")

    save_json(ARTIFACTS_DIR / "qubit_rotation_fock_resolved.json", all_results)
    print("Saved qubit_rotation_fock_resolved.json")
    return all_results


def run_drag_sweep():
    """Sweep DRAG coefficient for the X_pi pulse at fixed duration."""
    print("\n" + "=" * 60)
    print("DRAG COEFFICIENT SWEEP")
    print("=" * 60)
    t_start = time.time()

    model = build_model(n_cav=N_CAV_STUDY)
    frame = build_frame(model)

    theta = np.pi
    phase = 0.0
    T = 40e-9
    target_2x2 = np.asarray(qubit_rotation_xy(theta, phase).full(), dtype=np.complex128)

    drag_results = {}
    for drag_val in DRAG_VALUES:
        pulse = make_gaussian_qubit_pulse(
            model, frame, theta=theta, phase=phase, duration_s=T,
            drag=drag_val,
        )
        session = compile_and_prepare(model, frame, [pulse])
        U = propagate_basis_states(
            model, session,
            cavity_levels=CAVITY_LEVELS_PROBE,
            qubit_levels=[0, 1],
        )

        fock_met = {}
        for n in CAVITY_LEVELS_PROBE:
            block = extract_qubit_block(U, model.n_cav, n)
            met = rotation_metrics(block, target=target_2x2)
            fock_met[str(n)] = met

        # Also check leakage to |f>  
        psi0_g0 = model.basis_state(0, 0)
        psi_f = simulate_state(session, psi0_g0)
        full = np.asarray(psi_f.full()).ravel()
        leakage = sum(abs(full[2 * model.n_cav + n]) ** 2 for n in range(model.n_cav))

        drag_results[f"drag_{drag_val*1e9:.1f}ns"] = {
            "fock_metrics": fock_met,
            "leakage_to_f": float(leakage),
        }
        print(f"  DRAG={drag_val*1e9:.1f}ns: fid(n=0)={fock_met['0']['process_fidelity']:.6f}, "
              f"leakage={leakage:.6f}")

    save_json(ARTIFACTS_DIR / "qubit_rotation_drag_sweep.json", drag_results)
    elapsed = time.time() - t_start
    print(f"DRAG sweep completed in {elapsed:.1f} s")
    return drag_results


def run_entanglement_with_cavity():
    """Measure qubit-cavity entanglement after qubit pulse for different cavity states."""
    print("\n" + "=" * 60)
    print("ENTANGLEMENT WITH CAVITY OCCUPATION")
    print("=" * 60)

    model = build_model(n_cav=N_CAV_STUDY)
    frame = build_frame(model)

    theta = np.pi
    phase = 0.0
    T = 40e-9

    pulse = make_gaussian_qubit_pulse(model, frame, theta=theta, phase=phase, duration_s=T)
    session = compile_and_prepare(model, frame, [pulse])

    # Test coherent state initial conditions
    alpha_values = [0.0, 0.5, 1.0, 1.5, 2.0]
    ent_results = {}
    for alpha in alpha_values:
        if alpha == 0.0:
            cav_state = qt.basis(model.n_cav, 0)
        else:
            cav_state = qt.coherent(model.n_cav, alpha)
        psi0 = qt.tensor(qt.basis(model.n_tr, 0), cav_state)
        psi_f = simulate_state(session, psi0)

        ent = entanglement_entropy(psi_f, model.n_tr, model.n_cav)
        pur = qubit_purity(psi_f)

        # Ideal: R tensor I on cavity
        R_ideal = qubit_rotation_xy(theta, phase)
        R_full = qt.tensor(R_ideal, qt.qeye(model.n_cav))
        # Extend R to 3 levels (embed in n_tr)
        R_full_3 = qt.Qobj(
            np.eye(model.n_tr * model.n_cav, dtype=complex),
            dims=[[model.n_tr, model.n_cav], [model.n_tr, model.n_cav]]
        )
        full_mat = np.asarray(R_full_3.full())
        r2 = np.asarray(R_ideal.full())
        for n in range(model.n_cav):
            idx_g = n
            idx_e = model.n_cav + n
            full_mat[np.ix_([idx_g, idx_e], [idx_g, idx_e])] = r2
        ideal_out = qt.Qobj(full_mat, dims=R_full_3.dims) * psi0

        state_fid = float(abs(ideal_out.dag() * psi_f) ** 2)

        ent_results[f"alpha_{alpha:.1f}"] = {
            "entanglement_entropy": ent,
            "qubit_purity": pur,
            "state_fidelity": state_fid,
            "mean_photon": float(alpha ** 2),
        }
        print(f"  alpha={alpha:.1f}: ent={ent:.4f}, purity={pur:.4f}, fid={state_fid:.6f}")

    save_json(ARTIFACTS_DIR / "qubit_rotation_entanglement_vs_alpha.json", ent_results)
    return ent_results


def generate_figures(fock_results, drag_results, ent_results):
    """Generate all qubit rotation study figures."""
    print("\n" + "=" * 60)
    print("GENERATING FIGURES")
    print("=" * 60)

    # --- Figure 1: Fock-resolved rotation angle vs cavity level ---
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    for ax, (pulse_label, pulse_data) in zip(axes.ravel(), fock_results.items()):
        for i, T_key in enumerate(sorted(pulse_data.keys(), key=lambda k: float(k.replace("ns", "")))):
            met = pulse_data[T_key]
            angles = [met[str(n)]["rotation_angle_rad"] for n in CAVITY_LEVELS_PROBE]
            ax.plot(CAVITY_LEVELS_PROBE, angles, "o-",
                    color=TOL_BRIGHT[i % len(TOL_BRIGHT)], label=T_key, markersize=5)

        if "pi2" in pulse_label:
            ax.axhline(np.pi / 2, color="gray", linestyle="--", alpha=0.5)
        else:
            ax.axhline(np.pi, color="gray", linestyle="--", alpha=0.5)

        ax.set_xlabel("Cavity Fock level $n$")
        ax.set_ylabel("Rotation angle (rad)")
        ax.set_title(pulse_label.replace("_", " "))
        ax.legend(fontsize=7)

    fig.suptitle("Fock-Resolved Qubit Rotation Angle", y=1.02)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "qubit_rotation_fock_resolved_angle.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "qubit_rotation_fock_resolved_angle.pdf", bbox_inches="tight")
    plt.close(fig)

    # --- Figure 2: Fock-resolved process fidelity ---
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    for ax, (pulse_label, pulse_data) in zip(axes.ravel(), fock_results.items()):
        for i, T_key in enumerate(sorted(pulse_data.keys(), key=lambda k: float(k.replace("ns", "")))):
            met = pulse_data[T_key]
            fids = [met[str(n)]["process_fidelity"] for n in CAVITY_LEVELS_PROBE]
            ax.semilogy(CAVITY_LEVELS_PROBE, [1 - f for f in fids], "o-",
                        color=TOL_BRIGHT[i % len(TOL_BRIGHT)], label=T_key, markersize=5)

        ax.set_xlabel("Cavity Fock level $n$")
        ax.set_ylabel(r"$1 - \mathcal{F}_{\mathrm{proc}}$")
        ax.set_title(pulse_label.replace("_", " "))
        ax.legend(fontsize=7)

    fig.suptitle("Fock-Resolved Qubit Rotation Infidelity", y=1.02)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "qubit_rotation_fock_resolved_infidelity.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "qubit_rotation_fock_resolved_infidelity.pdf", bbox_inches="tight")
    plt.close(fig)

    # --- Figure 3: DRAG sweep ---
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    drag_ns = [float(k.split("_")[1].replace("ns", "")) for k in sorted(drag_results.keys())]
    fids_n0 = [drag_results[k]["fock_metrics"]["0"]["process_fidelity"] for k in sorted(drag_results.keys())]
    leakages = [drag_results[k]["leakage_to_f"] for k in sorted(drag_results.keys())]

    axes[0].plot(drag_ns, [1 - f for f in fids_n0], "o-", color=TOL_BRIGHT[0])
    axes[0].set_xlabel("DRAG coefficient (ns)")
    axes[0].set_ylabel(r"$1 - \mathcal{F}_{\mathrm{proc}}$ (n=0)")
    axes[0].set_yscale("log")
    axes[0].set_title("Rotation Infidelity vs. DRAG")

    axes[1].plot(drag_ns, leakages, "s-", color=TOL_BRIGHT[1])
    axes[1].set_xlabel("DRAG coefficient (ns)")
    axes[1].set_ylabel("Leakage to $|f\\rangle$")
    axes[1].set_yscale("log")
    axes[1].set_title("Leakage vs. DRAG")

    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "qubit_rotation_drag_sweep.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "qubit_rotation_drag_sweep.pdf", bbox_inches="tight")
    plt.close(fig)

    # --- Figure 4: Entanglement vs cavity occupation ---
    fig, ax = plt.subplots(figsize=(6, 4))
    alphas = sorted(ent_results.keys(), key=lambda k: float(k.split("_")[1]))
    alpha_vals = [float(k.split("_")[1]) for k in alphas]
    ents = [ent_results[k]["entanglement_entropy"] for k in alphas]
    fids = [ent_results[k]["state_fidelity"] for k in alphas]

    ax2 = ax.twinx()
    l1, = ax.plot(alpha_vals, ents, "o-", color=TOL_BRIGHT[0], label="Entanglement")
    l2, = ax2.plot(alpha_vals, [1 - f for f in fids], "s--", color=TOL_BRIGHT[1], label="Infidelity")

    ax.set_xlabel(r"Initial coherent state $|\alpha|$")
    ax.set_ylabel("Entanglement entropy (bits)", color=TOL_BRIGHT[0])
    ax2.set_ylabel(r"$1 - \mathcal{F}$", color=TOL_BRIGHT[1])
    ax2.set_yscale("log")
    ax.legend(handles=[l1, l2], loc="upper left")
    ax.set_title(r"Qubit $X_\pi$ Pulse: Entanglement with Cavity")

    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "qubit_rotation_entanglement_vs_cavity.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "qubit_rotation_entanglement_vs_cavity.pdf", bbox_inches="tight")
    plt.close(fig)

    # --- Figure 5: Error budget ---
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    T_key = "40ns"  # Representative duration
    for ax, (pulse_label, pulse_data) in zip(axes.ravel(), fock_results.items()):
        if T_key not in pulse_data:
            continue
        met = pulse_data[T_key]
        residual_z = [met[str(n)].get("residual_z_rad", 0.0) for n in CAVITY_LEVELS_PROBE]
        transverse = [met[str(n)].get("transverse_error_rad", 0.0) for n in CAVITY_LEVELS_PROBE]
        angle_err = [met[str(n)].get("angle_error_rad", 0.0) for n in CAVITY_LEVELS_PROBE]

        x = np.arange(len(CAVITY_LEVELS_PROBE))
        w = 0.25
        ax.bar(x - w, residual_z, w, label="Residual Z", color=TOL_BRIGHT[0])
        ax.bar(x, transverse, w, label="Transverse error", color=TOL_BRIGHT[1])
        ax.bar(x + w, angle_err, w, label="Angle error", color=TOL_BRIGHT[2])
        ax.set_xlabel("Cavity Fock level $n$")
        ax.set_ylabel("Error (rad)")
        ax.set_title(f"{pulse_label.replace('_', ' ')} (T={T_key})")
        ax.set_xticks(x)
        ax.set_xticklabels(CAVITY_LEVELS_PROBE)
        ax.legend(fontsize=7)

    fig.suptitle("Error Budget: Fock-Resolved Rotation Errors", y=1.02)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "qubit_rotation_error_budget.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "qubit_rotation_error_budget.pdf", bbox_inches="tight")
    plt.close(fig)

    print("All qubit rotation figures saved.")


if __name__ == "__main__":
    fock_results = run_fock_resolved_rotation_study()
    drag_results = run_drag_sweep()
    ent_results = run_entanglement_with_cavity()
    generate_figures(fock_results, drag_results, ent_results)
    print("\n=== QUBIT ROTATION STUDY COMPLETE ===")
