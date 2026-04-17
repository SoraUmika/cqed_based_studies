"""Echoed displacement study: D(alpha/2) -> X_pi -> D(alpha/2) -> X_pi.

The echo scheme aims to refocus the chi-induced frequency mismatch between
qubit branches during cavity displacement. When the qubit starts in a
superposition |+x>, the first half-displacement accrues a chi-dependent
phase error; the X_pi flip swaps the qubit state so the second half accrues
the opposite error, canceling to leading order O(chi*T).

This script compares four protocols:
  1. Bare displacement D(alpha) with a single square pulse
  2. Echoed displacement D(alpha/2)->X_pi->D(alpha/2)->X_pi
  3. Reference |g>-only displacement (best-case baseline)
  4. Reference |e>-only displacement (worst-case baseline)

For the echo, two pi-pulse variants are tested: Gaussian and DRAG.
Two fairness conventions are used:
  (a) Fixed displacement duration: T_disp is the same for bare and echo;
      echo total time = T_disp + 2*T_pi.
  (b) Fixed total duration: echo uses T_disp - 2*T_pi per displacement half;
      bare uses T_disp. Only valid when T_disp > 2*T_pi.

Units: rad/s and seconds throughout.
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
    displacement_fidelity, entanglement_entropy, qubit_purity,
    save_json, ARTIFACTS_DIR, FIGURES_DIR, DATA_DIR,
    apply_plot_style, TOL_BRIGHT,
)
from cqed_sim.core.ideal_gates import displacement_op
from cqed_sim.pulses.pulse import Pulse

common.apply_plot_style()

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
ALPHA_VALUES = np.array([0.5, 1.0, 2.0])
DISP_DURATIONS_NS = np.array([10, 20, 50, 100, 200])
PI_PULSE_DURATIONS_NS = np.array([20, 40])

# Initial qubit states to test
QUBIT_LABELS = {
    "g": 0,
    "e": 1,
    "plus_x": "superposition",
}

DRAG_COEFF_OPTIMAL = 0.5e-9  # 0.5 ns, from DRAG sweep result


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def make_pi_pulse(model, frame, *, t0, duration_s, drag=0.0, label="Xpi"):
    """Construct an X_pi qubit pulse at the given start time."""
    return make_gaussian_qubit_pulse(
        model, frame,
        theta=np.pi,
        phase=0.0,
        duration_s=duration_s,
        manifold_level=0,
        drag=drag,
        t0=t0,
        label=label,
    )


def make_initial_state(model, qubit_label):
    """Prepare the initial joint state |qubit> x |0_cav>."""
    if qubit_label == "g":
        return model.basis_state(0, 0)
    elif qubit_label == "e":
        return model.basis_state(1, 0)
    elif qubit_label == "plus_x":
        psi_g = model.basis_state(0, 0)
        psi_e = model.basis_state(1, 0)
        psi = (psi_g + psi_e).unit()
        return psi
    else:
        raise ValueError(f"Unknown qubit label: {qubit_label}")


def ideal_displaced_state(model, alpha, qubit_label):
    """Ideal target state: |qubit> x D(alpha)|0>.

    For the echo sequence D(a/2)->Xpi->D(a/2)->Xpi, if the pi pulses are
    ideal, the net qubit operation is identity (two pi rotations = 2pi = I),
    while each half-displacement should combine to D(alpha) on the cavity.
    So the target is the same as for a bare displacement: |qubit> x D(alpha)|0>.
    """
    d_op = displacement_op(model.n_cav, complex(alpha))
    cav_displaced = d_op * qt.basis(model.n_cav, 0)
    if qubit_label == "g":
        q_state = qt.basis(model.n_tr, 0)
    elif qubit_label == "e":
        q_state = qt.basis(model.n_tr, 1)
    elif qubit_label == "plus_x":
        q_state = (qt.basis(model.n_tr, 0) + qt.basis(model.n_tr, 1)).unit()
    else:
        raise ValueError(f"Unknown qubit label: {qubit_label}")
    return qt.tensor(q_state, cav_displaced)


def run_bare_displacement(model, frame, alpha, duration_s, initial_state):
    """Single square displacement pulse D(alpha)."""
    pulse = make_square_displacement_pulse(
        model, frame, alpha=complex(alpha), duration_s=duration_s
    )
    session = compile_and_prepare(model, frame, [pulse])
    return simulate_state(session, initial_state)


def run_echoed_displacement(
    model, frame, alpha, disp_duration_s, pi_duration_s, initial_state,
    *, drag=0.0,
):
    """Echoed sequence: D(alpha/2) -> X_pi -> D(alpha/2) -> X_pi.

    Each displacement half has duration disp_duration_s / 2.
    Total sequence time = disp_duration_s + 2 * pi_duration_s.
    """
    half_dur = disp_duration_s / 2.0
    half_alpha = complex(alpha) / 2.0

    t0_d1 = 0.0
    t0_pi1 = t0_d1 + half_dur
    t0_d2 = t0_pi1 + pi_duration_s
    t0_pi2 = t0_d2 + half_dur

    d1 = make_square_displacement_pulse(
        model, frame, alpha=half_alpha, duration_s=half_dur, t0=t0_d1, label="D1"
    )
    pi1 = make_pi_pulse(
        model, frame, t0=t0_pi1, duration_s=pi_duration_s, drag=drag, label="Xpi_1"
    )
    d2 = make_square_displacement_pulse(
        model, frame, alpha=half_alpha, duration_s=half_dur, t0=t0_d2, label="D2"
    )
    pi2 = make_pi_pulse(
        model, frame, t0=t0_pi2, duration_s=pi_duration_s, drag=drag, label="Xpi_2"
    )

    pulses = [d1, pi1, d2, pi2]
    session = compile_and_prepare(model, frame, pulses)
    return simulate_state(session, initial_state)


def evaluate_protocol(final_state, target_state, model):
    """Compute fidelity, entanglement, and purity metrics."""
    fid = displacement_fidelity(final_state, target_state)
    ent = entanglement_entropy(final_state, model.n_tr, model.n_cav)
    pur = qubit_purity(final_state)
    return {
        "fidelity": fid,
        "infidelity": 1.0 - fid,
        "entanglement_entropy_bits": ent,
        "qubit_purity": pur,
    }


# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------
def run_echo_study():
    print("=" * 60)
    print("ECHOED DISPLACEMENT STUDY")
    print("=" * 60)
    t_start = time.time()

    model = build_model()
    frame = build_frame(model)

    results = {}

    for alpha_val in ALPHA_VALUES:
        alpha_key = f"alpha_{alpha_val:.1f}"
        results[alpha_key] = {}

        for T_disp_ns in DISP_DURATIONS_NS:
            T_disp_s = T_disp_ns * 1e-9
            dur_key = f"{int(T_disp_ns)}ns"
            results[alpha_key][dur_key] = {}

            for q_label in QUBIT_LABELS:
                psi0 = make_initial_state(model, q_label)
                target = ideal_displaced_state(model, alpha_val, q_label)

                # --- Bare displacement ---
                psi_bare = run_bare_displacement(model, frame, alpha_val, T_disp_s, psi0)
                bare_metrics = evaluate_protocol(psi_bare, target, model)
                bare_metrics["protocol"] = "bare"
                bare_metrics["total_time_ns"] = T_disp_ns

                entry = {"bare": bare_metrics}

                # --- Echoed displacement with various pi-pulse settings ---
                for T_pi_ns in PI_PULSE_DURATIONS_NS:
                    T_pi_s = T_pi_ns * 1e-9

                    for drag_label, drag_val in [("gaussian", 0.0), ("drag", DRAG_COEFF_OPTIMAL)]:
                        echo_key = f"echo_{drag_label}_{int(T_pi_ns)}ns"

                        # Convention A: fixed displacement duration
                        psi_echo = run_echoed_displacement(
                            model, frame, alpha_val, T_disp_s, T_pi_s, psi0,
                            drag=drag_val,
                        )
                        echo_metrics = evaluate_protocol(psi_echo, target, model)
                        echo_metrics["protocol"] = f"echo_{drag_label}"
                        echo_metrics["pi_pulse_duration_ns"] = T_pi_ns
                        echo_metrics["total_time_ns"] = T_disp_ns + 2 * T_pi_ns
                        echo_metrics["convention"] = "fixed_disp_duration"
                        entry[echo_key + "_fixeddisp"] = echo_metrics

                        # Convention B: fixed total duration
                        # The echo must fit inside T_disp total, so each
                        # displacement half = (T_disp - 2*T_pi) / 2
                        if T_disp_ns > 2 * T_pi_ns:
                            effective_disp_s = (T_disp_s - 2 * T_pi_s)
                            # effective alpha per half to achieve same total displacement
                            # with reduced time: amplitude scales inversely
                            psi_echo_ft = run_echoed_displacement(
                                model, frame, alpha_val, effective_disp_s, T_pi_s, psi0,
                                drag=drag_val,
                            )
                            echo_ft_metrics = evaluate_protocol(psi_echo_ft, target, model)
                            echo_ft_metrics["protocol"] = f"echo_{drag_label}_fixedtotal"
                            echo_ft_metrics["pi_pulse_duration_ns"] = T_pi_ns
                            echo_ft_metrics["total_time_ns"] = T_disp_ns
                            echo_ft_metrics["effective_disp_duration_ns"] = (T_disp_ns - 2 * T_pi_ns)
                            echo_ft_metrics["convention"] = "fixed_total_duration"
                            entry[echo_key + "_fixedtotal"] = echo_ft_metrics

                results[alpha_key][dur_key][q_label] = entry
                print(
                    f"  alpha={alpha_val:.1f}, T={T_disp_ns}ns, qubit={q_label}: "
                    f"bare_fid={bare_metrics['fidelity']:.6f}, "
                    f"echo_gauss_20ns={entry.get('echo_gaussian_20ns_fixeddisp', {}).get('fidelity', 'N/A')}"
                )

    elapsed = time.time() - t_start
    print(f"\nSweep completed in {elapsed:.1f}s")

    # Save results
    save_json(
        ARTIFACTS_DIR / "echoed_displacement_results.json",
        results,
        description="Echoed displacement study: bare vs echo D(a/2)->Xpi->D(a/2)->Xpi comparison",
        load_instructions="Use common.load_json() to load. Keys: alpha_{val}/duration/qubit_label/protocol.",
    )

    return results


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def generate_figures(results):
    """Generate comparison figures for the echoed displacement study."""
    print("\nGenerating figures...")

    # --- Figure 1: Fidelity comparison bar chart for key cases ---
    _plot_fidelity_comparison(results)

    # --- Figure 2: Fidelity improvement ratio (echo / bare) heatmap ---
    _plot_improvement_heatmap(results)

    # --- Figure 3: Infidelity vs displacement duration for fixed alpha ---
    _plot_infidelity_vs_duration(results)

    # --- Figure 4: Effect of pi-pulse variant on echo fidelity ---
    _plot_pi_pulse_comparison(results)

    print("All figures saved.")


def _plot_fidelity_comparison(results):
    """Bar chart: bare vs echo fidelity for select (alpha, T) pairs, qubit in |+x>."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharey=True)

    for ax_idx, alpha_val in enumerate(ALPHA_VALUES):
        ax = axes[ax_idx]
        alpha_key = f"alpha_{alpha_val:.1f}"

        durations = []
        fid_bare = []
        fid_echo_gauss = []
        fid_echo_drag = []

        for T_ns in DISP_DURATIONS_NS:
            dur_key = f"{int(T_ns)}ns"
            entry = results[alpha_key][dur_key].get("plus_x", {})
            if not entry:
                continue
            durations.append(T_ns)
            fid_bare.append(entry.get("bare", {}).get("fidelity", 0))

            # Pick echo with 20ns pi-pulse, fixed-disp convention
            eg = entry.get("echo_gaussian_20ns_fixeddisp", {})
            ed = entry.get("echo_drag_20ns_fixeddisp", {})
            fid_echo_gauss.append(eg.get("fidelity", 0))
            fid_echo_drag.append(ed.get("fidelity", 0))

        x = np.arange(len(durations))
        w = 0.25
        ax.bar(x - w, fid_bare, w, label="Bare", color=TOL_BRIGHT[0])
        ax.bar(x, fid_echo_gauss, w, label=r"Echo (Gauss $\pi$)", color=TOL_BRIGHT[1])
        ax.bar(x + w, fid_echo_drag, w, label=r"Echo (DRAG $\pi$)", color=TOL_BRIGHT[2])
        ax.set_xticks(x)
        ax.set_xticklabels([f"{int(d)}" for d in durations])
        ax.set_xlabel("Displacement duration (ns)")
        ax.set_title(fr"$|\alpha| = {alpha_val:.1f}$")
        if ax_idx == 0:
            ax.set_ylabel("Fidelity")
        ax.set_ylim(0, 1.05)
        ax.axhline(0.99, color="gray", ls="--", lw=0.8, alpha=0.5)
        ax.legend(fontsize=7, loc="lower left")

    fig.suptitle(r"Echoed vs Bare Displacement — Qubit in $|{+}x\rangle$ (20 ns $\pi$-pulse, fixed disp. dur.)", fontsize=12)
    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(str(FIGURES_DIR / f"echo_fidelity_comparison.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Saved echo_fidelity_comparison")


def _plot_improvement_heatmap(results):
    """Heatmap of fidelity improvement ratio for |+x> with echo_gaussian_20ns."""
    fid_mat_bare = np.zeros((len(DISP_DURATIONS_NS), len(ALPHA_VALUES)))
    fid_mat_echo = np.zeros_like(fid_mat_bare)

    for j, alpha_val in enumerate(ALPHA_VALUES):
        alpha_key = f"alpha_{alpha_val:.1f}"
        for i, T_ns in enumerate(DISP_DURATIONS_NS):
            dur_key = f"{int(T_ns)}ns"
            entry = results[alpha_key][dur_key].get("plus_x", {})
            fid_mat_bare[i, j] = entry.get("bare", {}).get("fidelity", 0)
            fid_mat_echo[i, j] = entry.get("echo_gaussian_20ns_fixeddisp", {}).get("fidelity", 0)

    # Improvement = echo infidelity / bare infidelity (< 1 means echo is better)
    infid_bare = 1.0 - fid_mat_bare
    infid_echo = 1.0 - fid_mat_echo
    ratio = np.where(infid_bare > 1e-10, infid_echo / infid_bare, 1.0)

    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(
        ratio,
        aspect="auto",
        origin="lower",
        extent=[ALPHA_VALUES[0], ALPHA_VALUES[-1], DISP_DURATIONS_NS[0], DISP_DURATIONS_NS[-1]],
        cmap="RdYlGn_r",
        vmin=0, vmax=2,
    )
    ax.set_xlabel(r"$|\alpha|$")
    ax.set_ylabel("Displacement duration (ns)")
    ax.set_title(r"Echo Infidelity / Bare Infidelity ($|{+}x\rangle$, Gauss 20 ns $\pi$)")
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Infidelity ratio (< 1 = echo better)")
    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(str(FIGURES_DIR / f"echo_improvement_heatmap.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Saved echo_improvement_heatmap")


def _plot_infidelity_vs_duration(results):
    """Infidelity vs duration for each alpha, comparing bare and echo on |+x> and |e>."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 5), sharey=True)

    for ax_idx, alpha_val in enumerate(ALPHA_VALUES):
        ax = axes[ax_idx]
        alpha_key = f"alpha_{alpha_val:.1f}"

        for q_label, ls in [("plus_x", "-"), ("e", "--")]:
            infid_bare = []
            infid_echo = []
            durations = []
            for T_ns in DISP_DURATIONS_NS:
                dur_key = f"{int(T_ns)}ns"
                entry = results[alpha_key][dur_key].get(q_label, {})
                if not entry:
                    continue
                durations.append(T_ns)
                infid_bare.append(1.0 - entry.get("bare", {}).get("fidelity", 0))
                infid_echo.append(1.0 - entry.get("echo_drag_20ns_fixeddisp", {}).get("fidelity", 0))

            label_suffix = r"$|{+}x\rangle$" if q_label == "plus_x" else r"$|e\rangle$"
            ax.semilogy(durations, infid_bare, f"o{ls}", color=TOL_BRIGHT[0],
                        label=f"Bare {label_suffix}", markersize=5)
            ax.semilogy(durations, infid_echo, f"s{ls}", color=TOL_BRIGHT[2],
                        label=f"Echo {label_suffix}", markersize=5)

        ax.set_xlabel("Displacement duration (ns)")
        ax.set_title(fr"$|\alpha| = {alpha_val:.1f}$")
        if ax_idx == 0:
            ax.set_ylabel(r"Infidelity $1 - \mathcal{F}$")
        ax.axhline(0.01, color="gray", ls=":", lw=0.8)
        ax.legend(fontsize=7)

    fig.suptitle(r"Infidelity: Bare vs Echo (DRAG 20 ns $\pi$, fixed disp. dur.)", fontsize=12)
    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(str(FIGURES_DIR / f"echo_infidelity_vs_duration.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Saved echo_infidelity_vs_duration")


def _plot_pi_pulse_comparison(results):
    """Compare echo fidelity across pi-pulse variants for |+x> at alpha=1.0."""
    alpha_key = "alpha_1.0"
    q_label = "plus_x"

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    for ax_idx, T_pi_ns in enumerate(PI_PULSE_DURATIONS_NS):
        ax = axes[ax_idx]

        for drag_label, color, marker in [
            ("gaussian", TOL_BRIGHT[1], "o"),
            ("drag", TOL_BRIGHT[2], "s"),
        ]:
            echo_key_fd = f"echo_{drag_label}_{int(T_pi_ns)}ns_fixeddisp"
            echo_key_ft = f"echo_{drag_label}_{int(T_pi_ns)}ns_fixedtotal"

            durations = []
            infid_fd = []
            infid_ft = []
            for T_ns in DISP_DURATIONS_NS:
                dur_key = f"{int(T_ns)}ns"
                entry = results[alpha_key][dur_key].get(q_label, {})
                if not entry:
                    continue
                durations.append(T_ns)
                fd = entry.get(echo_key_fd, {})
                ft = entry.get(echo_key_ft, {})
                infid_fd.append(1.0 - fd.get("fidelity", 0) if fd else np.nan)
                infid_ft.append(1.0 - ft.get("fidelity", 0) if ft else np.nan)

            label_env = "Gauss" if drag_label == "gaussian" else "DRAG"
            ax.semilogy(durations, infid_fd, marker + "-", color=color,
                        label=f"{label_env} (fixed disp)", markersize=5)
            ax.semilogy(durations, infid_ft, marker + "--", color=color, alpha=0.6,
                        label=f"{label_env} (fixed total)", markersize=5)

        # Also plot bare for reference
        bare_infid = []
        dur_list = []
        for T_ns in DISP_DURATIONS_NS:
            dur_key = f"{int(T_ns)}ns"
            entry = results[alpha_key][dur_key].get(q_label, {})
            if entry:
                dur_list.append(T_ns)
                bare_infid.append(1.0 - entry.get("bare", {}).get("fidelity", 0))
        ax.semilogy(dur_list, bare_infid, "d-", color=TOL_BRIGHT[0], label="Bare", markersize=5)

        ax.set_xlabel("Displacement duration (ns)")
        ax.set_ylabel(r"Infidelity $1 - \mathcal{F}$")
        ax.set_title(fr"$\pi$-pulse = {int(T_pi_ns)} ns, $|\alpha|=1.0$, $|{{+}}x\rangle$")
        ax.axhline(0.01, color="gray", ls=":", lw=0.8)
        ax.legend(fontsize=7)

    fig.suptitle(r"Echo $\pi$-Pulse Variant Comparison", fontsize=12)
    plt.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(str(FIGURES_DIR / f"echo_pi_pulse_comparison.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("  Saved echo_pi_pulse_comparison")


# ---------------------------------------------------------------------------
# Analytic estimate
# ---------------------------------------------------------------------------
def analytic_estimate():
    """Quick toggling-frame analysis of the echo scheme.

    In the dispersive frame, the cavity picks up a qubit-state-dependent
    frequency shift chi. For a bare displacement of duration T:
      |g>: cavity driven on resonance -> D(alpha)
      |e>: cavity detuned by chi -> D(alpha * sinc(chi*T/2) * exp(-i*chi*T/2))

    For the echo D(a/2)->Xpi->D(a/2)->Xpi with ideal pi-pulses:
      First half (T/2): qubit in original state
      After Xpi: qubit flipped
      Second half (T/2): qubit in flipped state
      After Xpi: qubit restored

    The cavity phase from chi averages out to leading order,
    leaving a residual ~ O((chi*T/2)^2) error.
    """
    chi_val = CHI / TWO_PI  # Hz
    print("\n" + "=" * 60)
    print("ANALYTIC PILOT: TOGGLING-FRAME ECHO ESTIMATE")
    print("=" * 60)

    analytic_results = {}
    for T_ns in [10, 20, 50, 100, 200]:
        T = T_ns * 1e-9
        chi_T = abs(CHI) * T
        # Bare |e> infidelity scales as ~ (chi*T)^2 for small chi*T
        bare_error = chi_T ** 2
        # Echo residual scales as ~ (chi*T/2)^4 / 4 for symmetric echo
        echo_error = (chi_T / 2) ** 4 / 4.0
        improvement = bare_error / echo_error if echo_error > 0 else float("inf")
        print(f"  T={T_ns}ns: chi*T={chi_T:.3f} rad, "
              f"bare~{bare_error:.2e}, echo~{echo_error:.2e}, "
              f"improvement~{improvement:.1f}x")
        analytic_results[f"{T_ns}ns"] = {
            "chi_T_rad": float(chi_T),
            "bare_error_estimate": float(bare_error),
            "echo_error_estimate": float(echo_error),
            "improvement_factor": float(improvement),
        }

    save_json(
        ARTIFACTS_DIR / "echo_analytic_estimate.json",
        analytic_results,
        description="Analytic toggling-frame estimate for echo displacement error scaling",
    )

    return analytic_results


# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------
def print_summary(results):
    """Print a summary table of key results."""
    print("\n" + "=" * 60)
    print("SUMMARY: ECHOED DISPLACEMENT RESULTS")
    print("=" * 60)
    print(f"{'alpha':>6} {'T(ns)':>6} {'qubit':>7} {'bare':>10} {'echo_g20':>10} {'echo_d20':>10} {'improvement':>12}")
    print("-" * 70)

    for alpha_val in ALPHA_VALUES:
        alpha_key = f"alpha_{alpha_val:.1f}"
        for T_ns in DISP_DURATIONS_NS:
            dur_key = f"{int(T_ns)}ns"
            for q_label in ["plus_x", "e"]:
                entry = results[alpha_key][dur_key].get(q_label, {})
                if not entry:
                    continue
                fid_bare = entry.get("bare", {}).get("fidelity", 0)
                fid_eg = entry.get("echo_gaussian_20ns_fixeddisp", {}).get("fidelity", 0)
                fid_ed = entry.get("echo_drag_20ns_fixeddisp", {}).get("fidelity", 0)
                infid_bare = 1.0 - fid_bare
                infid_echo = 1.0 - max(fid_eg, fid_ed)
                ratio = infid_bare / infid_echo if infid_echo > 1e-10 else float("inf")
                print(f"{alpha_val:>6.1f} {T_ns:>6.0f} {q_label:>7} "
                      f"{fid_bare:>10.6f} {fid_eg:>10.6f} {fid_ed:>10.6f} {ratio:>12.1f}x")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    analytic_estimate()
    results = run_echo_study()
    generate_figures(results)
    print_summary(results)
    print("\nDone.")
