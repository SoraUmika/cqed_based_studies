"""Displacement pulse study: characterize how square cavity drive pulses depart
from ideal D(alpha) in the dispersive qubit+cavity system.

Sweeps over amplitude, duration, initial qubit state, and coupling ablations.
Saves all results and generates figures.
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
    make_square_displacement_pulse,
    compile_and_prepare, simulate_state,
    propagate_basis_states, extract_qubit_block, extract_cavity_block,
    displacement_fidelity, cavity_state_from_joint, entanglement_entropy,
    qubit_purity, fock_populations,
    save_json, ARTIFACTS_DIR, FIGURES_DIR, DATA_DIR,
    apply_plot_style, TOL_BRIGHT,
)
from cqed_sim.core.ideal_gates import displacement_op

common.apply_plot_style()

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
ALPHA_VALUES = np.linspace(0.1, 2.0, 20)
DURATION_VALUES = np.array([20, 50, 100, 200, 500]) * 1e-9  # ns -> s
QUBIT_STATES = {"g": 0, "e": 1}

# Ablation configs: (chi, chi_prime, kerr, label)
ABLATION_CONFIGS = [
    (0.0, None, 0.0, "ideal"),
    (CHI, None, 0.0, "chi_only"),
    (CHI, CHI_PRIME, 0.0, "chi_chiprime"),
    (CHI, None, KERR, "chi_kerr"),
    (CHI, CHI_PRIME, KERR, "full"),
]


def run_displacement_sweep():
    """Main displacement amplitude x duration sweep with full model."""
    print("=" * 60)
    print("DISPLACEMENT PULSE STUDY")
    print("=" * 60)
    t_start = time.time()

    model = build_model()
    frame = build_frame(model)

    results = {}
    for q_label, q_level in QUBIT_STATES.items():
        print(f"\n--- Qubit initial state: |{q_label}> ---")
        fid_grid = np.zeros((len(DURATION_VALUES), len(ALPHA_VALUES)))
        ent_grid = np.zeros_like(fid_grid)
        purity_grid = np.zeros_like(fid_grid)

        for i, T in enumerate(DURATION_VALUES):
            for j, alpha_mag in enumerate(ALPHA_VALUES):
                alpha = complex(alpha_mag)
                pulse = make_square_displacement_pulse(model, frame, alpha=alpha, duration_s=T)
                session = compile_and_prepare(model, frame, [pulse])

                psi0 = model.basis_state(q_level, 0)
                psi_f = simulate_state(session, psi0)

                # Ideal: qubit state unchanged, cavity displaced
                D = displacement_op(model.n_cav, alpha)
                cav_ideal = D * qt.basis(model.n_cav, 0)
                ideal = qt.tensor(qt.basis(model.n_tr, q_level), cav_ideal)

                fid = displacement_fidelity(psi_f, ideal)
                ent = entanglement_entropy(psi_f, model.n_tr, model.n_cav)
                pur = qubit_purity(psi_f)

                fid_grid[i, j] = fid
                ent_grid[i, j] = ent
                purity_grid[i, j] = pur

            print(f"  T={T*1e9:.0f}ns: min fid={fid_grid[i].min():.4f}, max ent={ent_grid[i].max():.4f}")

        results[q_label] = {
            "fidelity": fid_grid,
            "entanglement_entropy": ent_grid,
            "qubit_purity": purity_grid,
        }

    elapsed = time.time() - t_start
    print(f"\nSweep completed in {elapsed:.1f} s")

    # Save data
    np.savez(
        DATA_DIR / "displacement_sweep.npz",
        alpha_values=ALPHA_VALUES,
        duration_values=DURATION_VALUES,
        fidelity_g=results["g"]["fidelity"],
        fidelity_e=results["e"]["fidelity"],
        entanglement_g=results["g"]["entanglement_entropy"],
        entanglement_e=results["e"]["entanglement_entropy"],
        purity_g=results["g"]["qubit_purity"],
        purity_e=results["e"]["qubit_purity"],
    )
    print("Saved displacement_sweep.npz")
    return results


def run_ablation_study():
    """Compare displacement fidelity across coupling configurations."""
    print("\n" + "=" * 60)
    print("ABLATION STUDY (chi, chi', K)")
    print("=" * 60)
    t_start = time.time()

    alpha = 1.0
    T = 100e-9
    ablation_results = {}

    for chi_val, chip_val, kerr_val, label in ABLATION_CONFIGS:
        model = build_model(chi=chi_val, chi_prime=chip_val, kerr=kerr_val)
        frame = build_frame(model)
        pulse = make_square_displacement_pulse(model, frame, alpha=alpha, duration_s=T)
        session = compile_and_prepare(model, frame, [pulse])

        fids = {}
        ents = {}
        for q_label, q_level in QUBIT_STATES.items():
            psi0 = model.basis_state(q_level, 0)
            psi_f = simulate_state(session, psi0)

            D = displacement_op(model.n_cav, alpha)
            cav_ideal = D * qt.basis(model.n_cav, 0)
            ideal = qt.tensor(qt.basis(model.n_tr, q_level), cav_ideal)

            fids[q_label] = displacement_fidelity(psi_f, ideal)
            ents[q_label] = entanglement_entropy(psi_f, model.n_tr, model.n_cav)

        ablation_results[label] = {"fidelity": fids, "entanglement": ents}
        print(f"  {label:15s}: fid_g={fids['g']:.6f}, fid_e={fids['e']:.6f}, "
              f"ent_g={ents['g']:.4f}, ent_e={ents['e']:.4f}")

    save_json(ARTIFACTS_DIR / "displacement_ablation.json", ablation_results)
    elapsed = time.time() - t_start
    print(f"Ablation completed in {elapsed:.1f} s")
    return ablation_results


def run_time_resolved_trajectory():
    """Time-resolved phase-space trajectory for representative cases."""
    print("\n" + "=" * 60)
    print("TIME-RESOLVED TRAJECTORY")
    print("=" * 60)

    model = build_model()
    frame = build_frame(model)
    alpha = 1.0
    T = 200e-9

    pulse = make_square_displacement_pulse(model, frame, alpha=alpha, duration_s=T)
    session = compile_and_prepare(model, frame, [pulse], store_states=True)

    trajectories = {}
    for q_label, q_level in QUBIT_STATES.items():
        psi0 = model.basis_state(q_level, 0)
        result = session.run(psi0)

        times = np.asarray(result.solver_result.times)
        states = result.states

        # Sample at ~20 time points
        n_samples = min(20, len(states))
        sample_idx = np.linspace(0, len(states) - 1, n_samples, dtype=int)

        x_mean = []
        p_mean = []
        ent_traj = []
        fid_traj = []

        a = model.operators()["a"]
        x_op = (a + a.dag()) / np.sqrt(2)
        p_op = -1j * (a - a.dag()) / np.sqrt(2)

        for idx in sample_idx:
            st = states[idx]
            x_mean.append(float(np.real(qt.expect(x_op, st))))
            p_mean.append(float(np.real(qt.expect(p_op, st))))
            ent_traj.append(entanglement_entropy(st, model.n_tr, model.n_cav))

            # Instantaneous ideal target
            t_frac = times[idx] / T
            alpha_t = alpha * t_frac if times[idx] < T else alpha
            D_t = displacement_op(model.n_cav, alpha_t)
            cav_t = D_t * qt.basis(model.n_cav, 0)
            ideal_t = qt.tensor(qt.basis(model.n_tr, q_level), cav_t)
            fid_traj.append(displacement_fidelity(st, ideal_t))

        trajectories[q_label] = {
            "times_ns": (times[sample_idx] * 1e9).tolist(),
            "x_mean": x_mean,
            "p_mean": p_mean,
            "entanglement": ent_traj,
            "fidelity": fid_traj,
        }
        print(f"  |{q_label}>: final ent={ent_traj[-1]:.4f}, final fid={fid_traj[-1]:.4f}")

    save_json(ARTIFACTS_DIR / "displacement_trajectory.json", trajectories)
    print("Saved displacement_trajectory.json")
    return trajectories


def generate_figures(sweep_results, ablation_results, trajectory_results):
    """Generate all displacement study figures."""
    print("\n" + "=" * 60)
    print("GENERATING FIGURES")
    print("=" * 60)

    # --- Figure 1: Fidelity heatmaps ---
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, (q_label, data) in zip(axes, sweep_results.items()):
        fid = data["fidelity"]
        im = ax.pcolormesh(
            ALPHA_VALUES, DURATION_VALUES * 1e9, 1 - fid,
            shading="nearest", cmap="viridis",
        )
        ax.set_xlabel(r"Displacement $|\alpha|$")
        ax.set_ylabel("Pulse duration (ns)")
        ax.set_title(f"Qubit in $|{q_label}\\rangle$")
        plt.colorbar(im, ax=ax, label=r"$1 - \mathcal{F}$")
    fig.suptitle("Displacement Infidelity vs. Amplitude and Duration", y=1.02)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "displacement_infidelity_heatmap.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "displacement_infidelity_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)

    # --- Figure 2: Entanglement heatmaps ---
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, (q_label, data) in zip(axes, sweep_results.items()):
        ent = data["entanglement_entropy"]
        im = ax.pcolormesh(
            ALPHA_VALUES, DURATION_VALUES * 1e9, ent,
            shading="nearest", cmap="inferno",
        )
        ax.set_xlabel(r"Displacement $|\alpha|$")
        ax.set_ylabel("Pulse duration (ns)")
        ax.set_title(f"Qubit in $|{q_label}\\rangle$")
        plt.colorbar(im, ax=ax, label="Entanglement entropy (bits)")
    fig.suptitle("Qubit-Cavity Entanglement During Displacement", y=1.02)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "displacement_entanglement_heatmap.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "displacement_entanglement_heatmap.pdf", bbox_inches="tight")
    plt.close(fig)

    # --- Figure 3: Ablation bar chart ---
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    labels = [cfg[3] for cfg in ABLATION_CONFIGS]
    x = np.arange(len(labels))
    width = 0.35

    fid_g = [ablation_results[l]["fidelity"]["g"] for l in labels]
    fid_e = [ablation_results[l]["fidelity"]["e"] for l in labels]
    axes[0].bar(x - width / 2, [1 - f for f in fid_g], width, label=r"$|g\rangle$", color=TOL_BRIGHT[0])
    axes[0].bar(x + width / 2, [1 - f for f in fid_e], width, label=r"$|e\rangle$", color=TOL_BRIGHT[1])
    axes[0].set_ylabel(r"$1 - \mathcal{F}$")
    axes[0].set_xticks(x)
    axes[0].set_xticklabels(labels, rotation=30, ha="right")
    axes[0].set_yscale("log")
    axes[0].legend()
    axes[0].set_title("Displacement Infidelity by Coupling")

    ent_g = [ablation_results[l]["entanglement"]["g"] for l in labels]
    ent_e = [ablation_results[l]["entanglement"]["e"] for l in labels]
    axes[1].bar(x - width / 2, ent_g, width, label=r"$|g\rangle$", color=TOL_BRIGHT[0])
    axes[1].bar(x + width / 2, ent_e, width, label=r"$|e\rangle$", color=TOL_BRIGHT[1])
    axes[1].set_ylabel("Entanglement entropy (bits)")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=30, ha="right")
    axes[1].legend()
    axes[1].set_title("Entanglement by Coupling")

    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "displacement_ablation.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "displacement_ablation.pdf", bbox_inches="tight")
    plt.close(fig)

    # --- Figure 4: Phase-space trajectory ---
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, (q_label, traj) in zip(axes, trajectory_results.items()):
        ax.plot(traj["x_mean"], traj["p_mean"], "o-", color=TOL_BRIGHT[0], markersize=4)
        ax.plot(traj["x_mean"][0], traj["p_mean"][0], "s", color=TOL_BRIGHT[2], markersize=8, label="Start")
        ax.plot(traj["x_mean"][-1], traj["p_mean"][-1], "*", color=TOL_BRIGHT[1], markersize=12, label="End")
        ax.set_xlabel(r"$\langle X \rangle / \sqrt{2}$")
        ax.set_ylabel(r"$\langle P \rangle / \sqrt{2}$")
        ax.set_title(f"Qubit in $|{q_label}\\rangle$")
        ax.legend()
        ax.set_aspect("equal")
    fig.suptitle(r"Phase-Space Trajectory During Displacement ($\alpha=1$, $T=200$ ns)", y=1.02)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "displacement_phase_space.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "displacement_phase_space.pdf", bbox_inches="tight")
    plt.close(fig)

    # --- Figure 5: Fidelity vs alpha line cuts ---
    fig, ax = plt.subplots(figsize=(6, 4))
    for i, T in enumerate(DURATION_VALUES):
        ax.semilogy(ALPHA_VALUES, 1 - sweep_results["g"]["fidelity"][i],
                     label=f"$|g\\rangle$, T={T*1e9:.0f} ns", linestyle="-", color=TOL_BRIGHT[i % len(TOL_BRIGHT)])
        ax.semilogy(ALPHA_VALUES, 1 - sweep_results["e"]["fidelity"][i],
                     label=f"$|e\\rangle$, T={T*1e9:.0f} ns", linestyle="--", color=TOL_BRIGHT[i % len(TOL_BRIGHT)])
    ax.set_xlabel(r"Displacement $|\alpha|$")
    ax.set_ylabel(r"$1 - \mathcal{F}$")
    ax.set_title("Displacement Infidelity vs. Amplitude")
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "displacement_fidelity_vs_alpha.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "displacement_fidelity_vs_alpha.pdf", bbox_inches="tight")
    plt.close(fig)

    print("All displacement figures saved.")


if __name__ == "__main__":
    sweep_results = run_displacement_sweep()
    ablation_results = run_ablation_study()
    trajectory_results = run_time_resolved_trajectory()
    generate_figures(sweep_results, ablation_results, trajectory_results)
    print("\n=== DISPLACEMENT STUDY COMPLETE ===")
