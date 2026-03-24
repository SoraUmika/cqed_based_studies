"""
study_phase5.py — Systematic parameter sweeps for the gray-box adaptive control study.

This script runs six sweep studies:

  Phase 5.1: Chi mismatch sweep (extended range)
  Phase 5.2: Noise strength sweep (eval with realistic decoherence)
  Phase 5.3: Readout imperfection sweep
  Phase 5.4: Probe budget (n_shots) sweep
  Phase 5.5: Chi drift study
  Phase 5.6: Hamiltonian omission study (chi_higher and Kerr effects)

All results saved to data/ as .npz files.

Usage (from study directory):
    python scripts/study_phase5.py --phase 5.1   # run only one phase
    python scripts/study_phase5.py                # run all phases

Runtime: 1-4 hours for full sweep. Each sub-phase can be run independently.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

STUDY_DIR = Path(__file__).parent.parent
sys.path.insert(0, str(STUDY_DIR / "scripts"))

from models import (
    make_truth_model,
    make_learner_model,
    make_frame,
    make_grape_subspace,
    make_target_matrix,
    CHI_TRUE,
    CHI_PRIOR,
    CHI_HIGHER_TRUE,
    KERR_TRUE,
    N_SHOTS_PROBE,
    CONFUSION_MATRIX,
    TRUTH_NOISE,
)
from cqed_sim import NoiseSpec
from probe_library import (
    run_chi_ramsey_probe,
    infer_chi_from_probe,
    make_ramsey_delays,
    t2_star_from_noise,
)
from control import run_grape, run_grape_multistart, eval_on_model

DATA_DIR = STUDY_DIR / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)

# Shared GRAPE settings
GRAPE_MAXITER = 200
GRAPE_N_STEPS = 16
GRAPE_DT_S = 10e-9
GRAPE_SEED = 42  # fallback single seed
# Multi-start seeds for phase5 helpers (2 seeds for tractability)
GRAPE_SEEDS = [2, 9]

# Fixed mismatch for phases 5.2-5.6
FIXED_MISMATCH = 0.30
CHI_PRIOR_FIXED = CHI_TRUE * (1.0 - FIXED_MISMATCH)

# Probe settings
PROBE_N_DELAYS = 80
PROBE_N_PERIODS = 3.0


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def nominal_fidelity_for_chi_prior(
    chi_prior, truth_model, truth_frame, subspace, target, verbose=False
):
    """Run multi-start GRAPE on learner(chi_prior) and eval on truth (noiseless)."""
    learner = make_learner_model(chi=chi_prior)
    frame = make_frame(learner)
    gr, prob, _, _ = run_grape_multistart(learner, frame, subspace, target,
                                          seeds=GRAPE_SEEDS,
                                          n_steps=GRAPE_N_STEPS, dt_s=GRAPE_DT_S,
                                          maxiter=GRAPE_MAXITER)
    avg_f, _ = eval_on_model(gr, prob, truth_model, truth_frame, eval_noise=None)
    return avg_f


def perfect_fidelity(truth_model, truth_frame, subspace, target):
    """Run multi-start GRAPE on perfect model (chi=chi_true) and eval on truth."""
    learner = make_learner_model(chi=CHI_TRUE)
    frame = make_frame(learner)
    gr, prob, _, _ = run_grape_multistart(learner, frame, subspace, target,
                                          seeds=GRAPE_SEEDS,
                                          n_steps=GRAPE_N_STEPS, dt_s=GRAPE_DT_S,
                                          maxiter=GRAPE_MAXITER)
    avg_f, _ = eval_on_model(gr, prob, truth_model, truth_frame, eval_noise=None)
    return avg_f


def gray_box_fidelity(
    chi_prior, truth_model, truth_frame, subspace, target,
    confusion_mat, n_shots, rng, t2_star, noise_for_eval=None,
):
    """
    Run full gray-box pipeline and return fidelity on truth model.

    Parameters
    ----------
    chi_prior : float
        Learner's initial chi estimate.
    truth_model : DispersiveTransmonCavityModel
        Truth model (probed and evaluated).
    truth_frame : FrameSpec
    subspace : Subspace
    target : np.ndarray
    confusion_mat : np.ndarray
        SPAM confusion matrix for probe.
    n_shots : int
        Shots per probe data point.
    rng : np.random.Generator
    t2_star : float
        T2* used in probe forward model.
    noise_for_eval : NoiseSpec or None
        Noise for GRAPE evaluation.

    Returns
    -------
    tuple: (fidelity, chi_hat, chi_error_frac)
    """
    delays = make_ramsey_delays(chi_prior, n_periods=PROBE_N_PERIODS, n_points=PROBE_N_DELAYS)
    probe_data = run_chi_ramsey_probe(
        chi_true=truth_model.chi,
        chi_higher_true=float(truth_model.chi_higher[0]) if len(truth_model.chi_higher) > 0 else 0.0,
        t2_star=t2_star,
        confusion_matrix=confusion_mat,
        n_shots=int(n_shots),
        delays_s=delays,
        fock_levels=[1, 2, 3],
        rng=rng,
    )
    infer_result = infer_chi_from_probe(
        probe_data=probe_data,
        confusion_matrix=confusion_mat,
        n_shots=int(n_shots),
        chi_initial=chi_prior,
    )
    chi_hat = infer_result["chi"]
    chi_error = abs(chi_hat - truth_model.chi) / abs(truth_model.chi)

    learner_corr = make_learner_model(chi=chi_hat)
    frame_corr = make_frame(learner_corr)
    gr, prob, _, _ = run_grape_multistart(learner_corr, frame_corr, subspace, target,
                                          seeds=GRAPE_SEEDS,
                                          n_steps=GRAPE_N_STEPS, dt_s=GRAPE_DT_S,
                                          maxiter=GRAPE_MAXITER)
    avg_f, _ = eval_on_model(gr, prob, truth_model, truth_frame, eval_noise=noise_for_eval)
    return avg_f, chi_hat, chi_error


# ---------------------------------------------------------------------------
# Phase 5.1 — Chi mismatch sweep (extended range)
# ---------------------------------------------------------------------------


def run_phase_5_1(verbose=True):
    """
    Sweep chi_prior / chi_true ratios over extended range.

    chi_prior = chi_true * ratio  for ratio in [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.3]
    """
    print("\n=== Phase 5.1: Chi mismatch sweep ===")
    ratios = np.array([1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.3])
    chi_priors = CHI_TRUE * ratios   # both negative, ratio reduces magnitude
    truth_model = make_truth_model()
    truth_frame = make_frame(truth_model)
    subspace = make_grape_subspace()
    target = make_target_matrix()
    rng = np.random.default_rng(111)
    t2_star = t2_star_from_noise(TRUTH_NOISE)

    nominal_fs = np.zeros(len(ratios))
    gray_box_fs = np.zeros(len(ratios))
    perfect_fs = np.zeros(len(ratios))
    chi_hat_vals = np.zeros(len(ratios))
    chi_errors = np.zeros(len(ratios))

    for i, (ratio, chi_p) in enumerate(zip(ratios, chi_priors)):
        if verbose:
            print(f"  ratio={ratio:.1f}  chi_prior={chi_p/(2*np.pi)/1e6:.3f} MHz")

        nominal_fs[i] = nominal_fidelity_for_chi_prior(
            chi_p, truth_model, truth_frame, subspace, target
        )
        gb_f, chi_hat, chi_err = gray_box_fidelity(
            chi_p, truth_model, truth_frame, subspace, target,
            CONFUSION_MATRIX, N_SHOTS_PROBE, rng, t2_star
        )
        gray_box_fs[i] = gb_f
        chi_hat_vals[i] = chi_hat
        chi_errors[i] = chi_err
        perfect_fs[i] = perfect_fidelity(truth_model, truth_frame, subspace, target)

        if verbose:
            print(f"    nominal={nominal_fs[i]:.5f}  gray_box={gb_f:.5f}  perfect={perfect_fs[i]:.5f}")

    np.savez(DATA_DIR / "phase5_1_chi_mismatch.npz",
             ratios=ratios, chi_priors=chi_priors,
             nominal_fidelities=nominal_fs, gray_box_fidelities=gray_box_fs,
             perfect_fidelities=perfect_fs,
             chi_hat_values=chi_hat_vals, chi_error_fracs=chi_errors,
             chi_true=CHI_TRUE)
    print(f"  Saved: {DATA_DIR / 'phase5_1_chi_mismatch.npz'}")


# ---------------------------------------------------------------------------
# Phase 5.2 — Noise strength sweep
# ---------------------------------------------------------------------------


def run_phase_5_2(verbose=True):
    """
    Sweep noise strength (T1 values) for nominal, gray-box, and perfect methods.

    chi mismatch fixed at 30%.
    T1 in {inf, 100us, 50us, 20us}.
    """
    print("\n=== Phase 5.2: Noise strength sweep ===")
    t1_values = [None, 100e-6, 50e-6, 20e-6]
    t1_labels = ["inf", "100us", "50us", "20us"]

    truth_model = make_truth_model()
    truth_frame = make_frame(truth_model)
    subspace = make_grape_subspace()
    target = make_target_matrix()
    rng = np.random.default_rng(222)
    t2_star = t2_star_from_noise(TRUTH_NOISE)

    nominal_fs_noisy = np.zeros(len(t1_values))
    gray_box_fs_noisy = np.zeros(len(t1_values))
    perfect_fs_noisy = np.zeros(len(t1_values))

    chi_prior = CHI_PRIOR_FIXED

    # First run GRAPE once for each method (noiseless optimization)
    # Then evaluate with noise varying
    if verbose:
        print(f"  chi_prior = {chi_prior/(2*np.pi)/1e6:.3f} MHz ({FIXED_MISMATCH*100:.0f}% mismatch)")

    # Build GRAPE results once (noiseless) using multi-start
    learner_nom = make_learner_model(chi=chi_prior)
    frame_nom = make_frame(learner_nom)
    gr_nom, prob_nom, _, _ = run_grape_multistart(learner_nom, frame_nom, subspace, target,
                                                   seeds=GRAPE_SEEDS, maxiter=GRAPE_MAXITER)

    # Gray-box corrected
    delays = make_ramsey_delays(chi_prior, n_periods=PROBE_N_PERIODS, n_points=PROBE_N_DELAYS)
    probe_data = run_chi_ramsey_probe(
        chi_true=truth_model.chi,
        chi_higher_true=float(truth_model.chi_higher[0]) if len(truth_model.chi_higher) > 0 else 0.0,
        t2_star=t2_star,
        confusion_matrix=CONFUSION_MATRIX,
        n_shots=N_SHOTS_PROBE,
        delays_s=delays,
        fock_levels=[1, 2, 3],
        rng=rng,
    )
    infer_result = infer_chi_from_probe(probe_data, CONFUSION_MATRIX, N_SHOTS_PROBE, chi_prior)
    chi_hat = infer_result["chi"]
    learner_gb = make_learner_model(chi=chi_hat)
    frame_gb = make_frame(learner_gb)
    gr_gb, prob_gb, _, _ = run_grape_multistart(learner_gb, frame_gb, subspace, target,
                                                 seeds=GRAPE_SEEDS, maxiter=GRAPE_MAXITER)

    learner_perf = make_learner_model(chi=CHI_TRUE)
    frame_perf = make_frame(learner_perf)
    gr_perf, prob_perf, _, _ = run_grape_multistart(learner_perf, frame_perf, subspace, target,
                                                     seeds=GRAPE_SEEDS, maxiter=GRAPE_MAXITER)

    for i, (t1, label) in enumerate(zip(t1_values, t1_labels)):
        if t1 is None:
            noise = None
        else:
            tphi = TRUTH_NOISE.tphi
            noise = NoiseSpec(t1=float(t1), tphi=float(tphi) if tphi else None,
                              kappa=TRUTH_NOISE.kappa)

        nom_f, _ = eval_on_model(gr_nom, prob_nom, truth_model, truth_frame, eval_noise=noise)
        gb_f, _ = eval_on_model(gr_gb, prob_gb, truth_model, truth_frame, eval_noise=noise)
        perf_f, _ = eval_on_model(gr_perf, prob_perf, truth_model, truth_frame, eval_noise=noise)

        nominal_fs_noisy[i] = nom_f
        gray_box_fs_noisy[i] = gb_f
        perfect_fs_noisy[i] = perf_f

        if verbose:
            print(f"  T1={label}: nominal={nom_f:.5f}  gray_box={gb_f:.5f}  perfect={perf_f:.5f}")

    t1_array = np.array([float("inf") if t is None else t for t in t1_values])
    np.savez(DATA_DIR / "phase5_2_noise_sweep.npz",
             t1_values=t1_array, t1_labels=np.array(t1_labels),
             nominal_fidelities=nominal_fs_noisy,
             gray_box_fidelities=gray_box_fs_noisy,
             perfect_fidelities=perfect_fs_noisy,
             chi_hat=np.array(chi_hat), chi_prior=np.array(chi_prior))
    print(f"  Saved: {DATA_DIR / 'phase5_2_noise_sweep.npz'}")


# ---------------------------------------------------------------------------
# Phase 5.3 — Readout imperfection sweep
# ---------------------------------------------------------------------------


def run_phase_5_3(verbose=True):
    """
    Sweep readout error rate (confusion_01) for chi inference quality and gate fidelity.

    confusion_01: probability of measuring |e> when state is |g> (off-diagonal M[1,0]).
    confusion_01 in {0.0, 0.02, 0.05, 0.10, 0.15}.
    """
    print("\n=== Phase 5.3: Readout imperfection sweep ===")
    confusion_01_values = np.array([0.0, 0.02, 0.05, 0.10, 0.15])
    # Perfect complement: M[0,1] = 2 * confusion_01 (assume symmetric-ish)
    # Actually, fix M[0,0]=1-c01, M[1,0]=c01, M[0,1]=2*c01, M[1,1]=1-2*c01
    # (consistent with increasing readout degradation)

    truth_model = make_truth_model()
    truth_frame = make_frame(truth_model)
    subspace = make_grape_subspace()
    target = make_target_matrix()
    t2_star = t2_star_from_noise(TRUTH_NOISE)
    rng = np.random.default_rng(333)

    chi_prior = CHI_PRIOR_FIXED
    chi_errors = np.zeros(len(confusion_01_values))
    chi_hat_values = np.zeros(len(confusion_01_values))
    nominal_fs = np.zeros(len(confusion_01_values))
    gray_box_fs = np.zeros(len(confusion_01_values))
    perfect_fs = np.zeros(len(confusion_01_values))

    # Nominal and perfect GRAPE are confusion-independent
    nom_f = nominal_fidelity_for_chi_prior(chi_prior, truth_model, truth_frame, subspace, target)
    perf_f = perfect_fidelity(truth_model, truth_frame, subspace, target)
    nominal_fs[:] = nom_f
    perfect_fs[:] = perf_f

    for i, c01 in enumerate(confusion_01_values):
        # Build confusion matrix with this error rate
        c10 = 2.0 * float(c01)  # assume complementary readout errors
        conf_mat = np.array([
            [1.0 - float(c01), float(c10)],
            [float(c01), 1.0 - float(c10)],
        ])
        conf_mat = np.clip(conf_mat, 0.0, 1.0)

        gb_f, chi_hat, chi_err = gray_box_fidelity(
            chi_prior, truth_model, truth_frame, subspace, target,
            conf_mat, N_SHOTS_PROBE, rng, t2_star
        )
        gray_box_fs[i] = gb_f
        chi_hat_values[i] = chi_hat
        chi_errors[i] = chi_err

        if verbose:
            print(f"  c01={c01:.3f}: chi_err={chi_err*100:.2f}%  "
                  f"gray_box={gb_f:.5f}  nominal={nom_f:.5f}  perfect={perf_f:.5f}")

    np.savez(DATA_DIR / "phase5_3_readout_sweep.npz",
             confusion_01_values=confusion_01_values,
             chi_errors=chi_errors, chi_hat_values=chi_hat_values,
             nominal_fidelities=nominal_fs, gray_box_fidelities=gray_box_fs,
             perfect_fidelities=perfect_fs, chi_prior=np.array(chi_prior))
    print(f"  Saved: {DATA_DIR / 'phase5_3_readout_sweep.npz'}")


# ---------------------------------------------------------------------------
# Phase 5.4 — Probe budget sweep
# ---------------------------------------------------------------------------


def run_phase_5_4(verbose=True):
    """
    Sweep number of measurement shots per probe point.

    n_shots in {50, 100, 300, 500, 1000, 3000}.
    """
    print("\n=== Phase 5.4: Probe budget sweep ===")
    n_shots_values = np.array([50, 100, 300, 500, 1000, 3000])

    truth_model = make_truth_model()
    truth_frame = make_frame(truth_model)
    subspace = make_grape_subspace()
    target = make_target_matrix()
    t2_star = t2_star_from_noise(TRUTH_NOISE)
    rng = np.random.default_rng(444)

    chi_prior = CHI_PRIOR_FIXED
    chi_errors = np.zeros(len(n_shots_values))
    chi_hat_values = np.zeros(len(n_shots_values))
    chi_uncertainties = np.zeros(len(n_shots_values))
    gray_box_fs = np.zeros(len(n_shots_values))

    # Nominal and perfect only once
    nom_f = nominal_fidelity_for_chi_prior(chi_prior, truth_model, truth_frame, subspace, target)
    perf_f = perfect_fidelity(truth_model, truth_frame, subspace, target)

    for i, n_shots in enumerate(n_shots_values):
        delays = make_ramsey_delays(chi_prior, n_periods=PROBE_N_PERIODS, n_points=PROBE_N_DELAYS)
        probe_data = run_chi_ramsey_probe(
            chi_true=truth_model.chi,
            chi_higher_true=float(truth_model.chi_higher[0]) if len(truth_model.chi_higher) > 0 else 0.0,
            t2_star=t2_star,
            confusion_matrix=CONFUSION_MATRIX,
            n_shots=int(n_shots),
            delays_s=delays,
            fock_levels=[1, 2, 3],
            rng=rng,
        )
        infer_result = infer_chi_from_probe(probe_data, CONFUSION_MATRIX, int(n_shots), chi_prior)
        chi_hat = infer_result["chi"]
        chi_err = abs(chi_hat - truth_model.chi) / abs(truth_model.chi)
        chi_errors[i] = chi_err
        chi_hat_values[i] = chi_hat
        chi_uncertainties[i] = infer_result["chi_uncertainty"]

        learner = make_learner_model(chi=chi_hat)
        frame = make_frame(learner)
        gr, prob, _, _ = run_grape_multistart(learner, frame, subspace, target,
                                              seeds=GRAPE_SEEDS, maxiter=GRAPE_MAXITER)
        gb_f, _ = eval_on_model(gr, prob, truth_model, truth_frame, eval_noise=None)
        gray_box_fs[i] = gb_f

        if verbose:
            print(f"  n_shots={n_shots}: chi_err={chi_err*100:.2f}%  gray_box={gb_f:.5f}")

    np.savez(DATA_DIR / "phase5_4_probe_budget.npz",
             n_shots_values=n_shots_values,
             chi_errors=chi_errors, chi_hat_values=chi_hat_values,
             chi_uncertainties=chi_uncertainties,
             gray_box_fidelities=gray_box_fs,
             nominal_fidelity=np.array(nom_f),
             perfect_fidelity=np.array(perf_f),
             chi_prior=np.array(chi_prior))
    print(f"  Saved: {DATA_DIR / 'phase5_4_probe_budget.npz'}")


# ---------------------------------------------------------------------------
# Phase 5.5 — Chi drift study
# ---------------------------------------------------------------------------


def run_phase_5_5(verbose=True):
    """
    Simulate slow chi drift and compare recalibration strategies.

    chi_true_t = chi_true * (1 + drift_rate * t) where t is the experiment cycle index.
    drift_rates in {0, 0.1%, 0.5%, 1%} of chi_true per cycle.
    Recalibration strategies: never, every 5 cycles, every 20 cycles.
    Total cycles: 100.
    """
    print("\n=== Phase 5.5: Chi drift study ===")
    N_CYCLES = 100
    drift_rates = [0.0, 0.001, 0.005, 0.01]
    drift_labels = ["0%", "0.1%", "0.5%", "1%"]
    recal_intervals = [None, 5, 20]   # None = never
    recal_labels = ["never", "every_5", "every_20"]

    subspace = make_grape_subspace()
    target = make_target_matrix()
    rng = np.random.default_rng(555)
    t2_star = t2_star_from_noise(TRUTH_NOISE)

    # Shape: (n_drift_rates, n_recal_strategies, n_cycles)
    fidelity_traces = np.zeros((len(drift_rates), len(recal_intervals), N_CYCLES))

    for di, (drift_rate, d_label) in enumerate(zip(drift_rates, drift_labels)):
        if verbose:
            print(f"  drift_rate = {d_label}")

        for ri, (recal_interval, r_label) in enumerate(zip(recal_intervals, recal_labels)):
            # Current effective chi at start (prior)
            chi_current = CHI_PRIOR_FIXED  # start from learner's initial prior

            # Build initial GRAPE pulse (prior model) — single seed for drift study tractability
            learner_init = make_learner_model(chi=chi_current)
            frame_init = make_frame(learner_init)
            gr_current, prob_current, _ = run_grape(
                learner_init, frame_init, subspace, target,
                maxiter=GRAPE_MAXITER, seed=2
            )

            for cycle in range(N_CYCLES):
                # Current true chi at this cycle
                chi_true_t = CHI_TRUE * (1.0 + drift_rate * cycle)

                # Build truth model at current chi
                truth_model_t = make_learner_model(
                    chi=chi_true_t,
                    chi_higher_val=float(CHI_HIGHER_TRUE),
                    kerr_val=float(KERR_TRUE),
                )
                truth_frame_t = make_frame(truth_model_t)

                # Evaluate current pulse on current truth model
                avg_f, _ = eval_on_model(gr_current, prob_current, truth_model_t, truth_frame_t)
                fidelity_traces[di, ri, cycle] = avg_f

                # Recalibration trigger
                should_recal = (
                    recal_interval is not None
                    and (cycle + 1) % recal_interval == 0
                )
                if should_recal:
                    # Run gray-box probe and update
                    chi_prior_t = chi_current  # use current estimate as prior
                    delays = make_ramsey_delays(chi_prior_t, n_periods=PROBE_N_PERIODS, n_points=PROBE_N_DELAYS)
                    probe_data = run_chi_ramsey_probe(
                        chi_true=chi_true_t,
                        chi_higher_true=float(CHI_HIGHER_TRUE),
                        t2_star=t2_star,
                        confusion_matrix=CONFUSION_MATRIX,
                        n_shots=N_SHOTS_PROBE,
                        delays_s=delays,
                        rng=rng,
                    )
                    infer_result = infer_chi_from_probe(
                        probe_data, CONFUSION_MATRIX, N_SHOTS_PROBE, chi_prior_t
                    )
                    chi_current = infer_result["chi"]

                    # Re-run GRAPE with updated chi — single seed for tractability
                    learner_new = make_learner_model(chi=chi_current)
                    frame_new = make_frame(learner_new)
                    gr_current, prob_current, _ = run_grape(
                        learner_new, frame_new, subspace, target,
                        maxiter=GRAPE_MAXITER, seed=2
                    )

            if verbose:
                final_f = fidelity_traces[di, ri, -1]
                mean_f = np.mean(fidelity_traces[di, ri, :])
                print(f"    {r_label}: mean={mean_f:.5f}  final={final_f:.5f}")

    np.savez(DATA_DIR / "phase5_5_drift.npz",
             fidelity_traces=fidelity_traces,
             drift_rates=np.array(drift_rates),
             drift_labels=np.array(drift_labels),
             recal_intervals=np.array([r if r is not None else -1 for r in recal_intervals]),
             recal_labels=np.array(recal_labels),
             n_cycles=np.array(N_CYCLES),
             chi_true=np.array(CHI_TRUE))
    print(f"  Saved: {DATA_DIR / 'phase5_5_drift.npz'}")


# ---------------------------------------------------------------------------
# Phase 5.6 — Hamiltonian omission study
# ---------------------------------------------------------------------------


def run_phase_5_6(verbose=True):
    """
    Study the effect of omitting chi_higher from the learner's model.

    The truth model always has the full chi_higher_true and kerr_true.
    The learner's GRAPE uses chi=chi_true but varying chi_higher knowledge.

    chi_higher_multiplier in {0, 0.5, 1, 2} of chi_higher_true
    (0 = omitted, 1 = perfectly known, 2 = overestimated).
    """
    print("\n=== Phase 5.6: Hamiltonian omission study ===")
    chi_higher_multipliers = np.array([0.0, 0.5, 1.0, 2.0])

    truth_model = make_truth_model()
    truth_frame = make_frame(truth_model)
    subspace = make_grape_subspace()
    target = make_target_matrix()

    fidelities_no_corr = np.zeros(len(chi_higher_multipliers))   # chi correct, chi_h varied
    fidelities_chi_also_wrong = np.zeros(len(chi_higher_multipliers))  # chi_prior + chi_h varied

    for i, mult in enumerate(chi_higher_multipliers):
        chi_higher_used = CHI_HIGHER_TRUE * mult

        # Scenario A: learner knows chi_true but not chi_higher
        learner_a = make_learner_model(chi=CHI_TRUE, chi_higher_val=chi_higher_used)
        frame_a = make_frame(learner_a)
        gr_a, prob_a, _, _ = run_grape_multistart(learner_a, frame_a, subspace, target,
                                                   seeds=GRAPE_SEEDS, maxiter=GRAPE_MAXITER)
        f_a, _ = eval_on_model(gr_a, prob_a, truth_model, truth_frame)
        fidelities_no_corr[i] = f_a

        # Scenario B: learner has chi_prior (30% off) AND wrong chi_higher
        learner_b = make_learner_model(chi=CHI_PRIOR_FIXED, chi_higher_val=chi_higher_used)
        frame_b = make_frame(learner_b)
        gr_b, prob_b, _, _ = run_grape_multistart(learner_b, frame_b, subspace, target,
                                                   seeds=GRAPE_SEEDS, maxiter=GRAPE_MAXITER)
        f_b, _ = eval_on_model(gr_b, prob_b, truth_model, truth_frame)
        fidelities_chi_also_wrong[i] = f_b

        if verbose:
            print(f"  mult={mult:.1f}: chi_h={chi_higher_used/(2*np.pi)/1e3:.1f} kHz  "
                  f"f_chi_correct={f_a:.5f}  f_chi_wrong={f_b:.5f}")

    np.savez(DATA_DIR / "phase5_6_omission.npz",
             chi_higher_multipliers=chi_higher_multipliers,
             chi_higher_values=CHI_HIGHER_TRUE * chi_higher_multipliers,
             fidelities_chi_correct=fidelities_no_corr,
             fidelities_chi_wrong=fidelities_chi_also_wrong,
             chi_true=np.array(CHI_TRUE),
             chi_higher_true=np.array(CHI_HIGHER_TRUE),
             kerr_true=np.array(KERR_TRUE))
    print(f"  Saved: {DATA_DIR / 'phase5_6_omission.npz'}")


# ---------------------------------------------------------------------------
# CLI entrypoint
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Study Phase 5: systematic parameter sweeps for gray-box adaptive control"
    )
    parser.add_argument(
        "--phase", type=str, default="all",
        help="Which phase to run: 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, or 'all'",
    )
    parser.add_argument("--quiet", action="store_true", help="Suppress progress output")
    args = parser.parse_args()

    verbose = not args.quiet
    t0 = time.perf_counter()

    runners = {
        "5.1": run_phase_5_1,
        "5.2": run_phase_5_2,
        "5.3": run_phase_5_3,
        "5.4": run_phase_5_4,
        "5.5": run_phase_5_5,
        "5.6": run_phase_5_6,
    }

    if args.phase == "all":
        for key, runner in runners.items():
            runner(verbose=verbose)
    elif args.phase in runners:
        runners[args.phase](verbose=verbose)
    else:
        print(f"Unknown phase '{args.phase}'. Choose from: {', '.join(runners.keys())} or 'all'.")
        sys.exit(1)

    print(f"\nTotal wall time: {time.perf_counter() - t0:.1f} s")


if __name__ == "__main__":
    main()
