"""Extensive simultaneous multitone SQR study based on cqed_sim.

This script answers three separate questions:

1. When does one common multitone waveform already implement multiple SQR-like
   rotations with acceptable branch-local fidelity?
2. When does that same waveform survive a strict logical targeted-subspace
   benchmark?
3. Can simple corrections rescue the hard large-angle regime?

Outputs
-------
data/simultaneous_multitone_sqr_results.npz
data/simultaneous_multitone_sqr_summary.json
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from common import (
    CASE_LABELS,
    CASE_SPECS,
    CHI_T_VALUES,
    DATA_DIR,
    DT_S,
    LOGICAL_LEVELS,
    N_TR_QUTRIT,
    OBJECTIVE_WEIGHTS,
    THETA_VALUES,
    build_frame,
    build_model,
    build_targets,
    duration_from_chi_t,
    logical_indices,
    process_fidelity,
    run_config_for_chi_t,
    target_and_spectator_means,
)

from cqed_sim.calibration import (
    ConditionedMultitoneCorrections,
    ConditionedOptimizationConfig,
    build_block_rotation_target_operator,
    build_spanning_state_transfer_set,
    optimize_conditioned_multitone,
    run_conditioned_multitone_validation,
    run_targeted_subspace_multitone_validation,
)
from cqed_sim.core.ideal_gates import logical_block_phase_op
from cqed_sim.sim import SimulationConfig, prepare_simulation


OUT_NPZ = DATA_DIR / "simultaneous_multitone_sqr_results.npz"
OUT_JSON = DATA_DIR / "simultaneous_multitone_sqr_summary.json"


def _target_operator(target_levels: tuple[int, ...], theta: float) -> tuple:
    targets = build_targets(target_levels, theta)
    operator = build_block_rotation_target_operator(targets, logical_levels=LOGICAL_LEVELS)
    transfer = build_spanning_state_transfer_set(operator)
    return targets, operator, transfer


def _baseline_bundle(
    model,
    run_config,
    target_levels: tuple[int, ...],
    theta: float,
):
    targets, operator, transfer = _target_operator(target_levels, theta)
    reduced = run_conditioned_multitone_validation(
        model,
        targets,
        run_config,
        simulation_mode="reduced",
    )
    strict = run_targeted_subspace_multitone_validation(
        model,
        targets,
        run_config,
        logical_levels=LOGICAL_LEVELS,
        target_operator=operator,
        transfer_set=transfer,
        objective_weights=OBJECTIVE_WEIGHTS,
    )
    return targets, operator, reduced, strict


def run_baseline_grid() -> dict[str, object]:
    """Run the baseline simultaneous multitone scan."""
    case_names = tuple(CASE_SPECS.keys())
    n_cases = len(case_names)
    n_chi = len(CHI_T_VALUES)
    n_theta = len(THETA_VALUES)
    n_levels = len(LOGICAL_LEVELS)

    reduced_fidelity = np.zeros((n_cases, n_chi, n_theta), dtype=float)
    reduced_target_mean = np.zeros_like(reduced_fidelity)
    reduced_spectator_mean = np.zeros_like(reduced_fidelity)
    reduced_target_theta_mean = np.zeros_like(reduced_fidelity)
    reduced_target_theta_ratio = np.zeros_like(reduced_fidelity)
    reduced_sector_fidelity = np.zeros((n_cases, n_chi, n_theta, n_levels), dtype=float)

    strict_fidelity = np.zeros_like(reduced_fidelity)
    compiled_fidelity = np.zeros_like(reduced_fidelity)
    state_transfer_mean = np.zeros_like(reduced_fidelity)
    state_transfer_min = np.zeros_like(reduced_fidelity)
    same_block_mean = np.zeros_like(reduced_fidelity)
    leakage_mean = np.zeros_like(reduced_fidelity)
    block_phase_rms = np.zeros_like(reduced_fidelity)
    strict_sector_fidelity = np.zeros((n_cases, n_chi, n_theta, n_levels), dtype=float)
    best_fit_phase = np.zeros((n_cases, n_chi, n_theta, n_levels), dtype=float)

    model = build_model()
    for case_index, case_name in enumerate(case_names):
        target_levels = CASE_SPECS[case_name]
        print(f"[baseline] {case_name}")
        for chi_index, chi_t in enumerate(CHI_T_VALUES):
            run_config = run_config_for_chi_t(model, float(chi_t))
            for theta_index, theta in enumerate(THETA_VALUES):
                _, _, reduced, strict = _baseline_bundle(model, run_config, target_levels, float(theta))
                reduced_fidelity[case_index, chi_index, theta_index] = float(reduced.weighted_mean_fidelity)
                target_mean, spectator_mean = target_and_spectator_means(reduced.sector_metrics, target_levels)
                reduced_target_mean[case_index, chi_index, theta_index] = target_mean
                reduced_spectator_mean[case_index, chi_index, theta_index] = spectator_mean
                target_theta_values = [
                    float(metric.theta_simulated_rad)
                    for metric in reduced.sector_metrics
                    if int(metric.n) in set(target_levels)
                ]
                theta_target = float(theta)
                theta_mean = float(np.mean(target_theta_values)) if target_theta_values else 0.0
                reduced_target_theta_mean[case_index, chi_index, theta_index] = theta_mean
                reduced_target_theta_ratio[case_index, chi_index, theta_index] = (
                    theta_mean / theta_target if abs(theta_target) > 1.0e-15 else 0.0
                )
                reduced_sector_fidelity[case_index, chi_index, theta_index, :] = [
                    float(metric.fidelity) for metric in reduced.sector_metrics
                ]

                strict_fidelity[case_index, chi_index, theta_index] = float(strict.restricted_process_fidelity)
                compiled_fidelity[case_index, chi_index, theta_index] = float(strict.best_fit_restricted_process_fidelity)
                state_transfer_mean[case_index, chi_index, theta_index] = float(strict.state_transfer_fidelity_mean)
                state_transfer_min[case_index, chi_index, theta_index] = float(strict.state_transfer_fidelity_min)
                same_block_mean[case_index, chi_index, theta_index] = float(strict.same_block_population_mean)
                leakage_mean[case_index, chi_index, theta_index] = float(strict.leakage_outside_target_mean)
                block_phase_rms[case_index, chi_index, theta_index] = float(
                    strict.block_phase_diagnostics.rms_block_phase_error_rad
                    if strict.block_phase_diagnostics is not None
                    else np.nan
                )
                strict_sector_fidelity[case_index, chi_index, theta_index, :] = [
                    float(metric.fidelity) for metric in strict.conditioned_sector_metrics
                ]
                best_fit_phase[case_index, chi_index, theta_index, :] = np.asarray(
                    strict.best_fit_logical_block_phase.phases_for_levels(LOGICAL_LEVELS),
                    dtype=float,
                )

    return {
        "case_names": case_names,
        "reduced_fidelity": reduced_fidelity,
        "reduced_target_mean": reduced_target_mean,
        "reduced_spectator_mean": reduced_spectator_mean,
        "reduced_target_theta_mean": reduced_target_theta_mean,
        "reduced_target_theta_ratio": reduced_target_theta_ratio,
        "reduced_sector_fidelity": reduced_sector_fidelity,
        "strict_fidelity": strict_fidelity,
        "compiled_fidelity": compiled_fidelity,
        "state_transfer_mean": state_transfer_mean,
        "state_transfer_min": state_transfer_min,
        "same_block_mean": same_block_mean,
        "leakage_mean": leakage_mean,
        "block_phase_rms": block_phase_rms,
        "strict_sector_fidelity": strict_sector_fidelity,
        "best_fit_phase": best_fit_phase,
    }


def run_target_amplitude_scan() -> dict[str, object]:
    """Scan a simple target-tone amplitude correction family for a hard case."""
    case_name = "pair_adjacent"
    target_levels = CASE_SPECS[case_name]
    theta = float(np.pi)
    chi_t = 3.0
    scales = np.linspace(-1.5, 3.0, 19, dtype=float)

    model = build_model()
    run_config = run_config_for_chi_t(model, chi_t)
    targets = build_targets(target_levels, theta)

    reduced_fidelity = np.zeros_like(scales)
    target_mean = np.zeros_like(scales)
    spectator_mean = np.zeros_like(scales)
    per_sector = np.zeros((scales.size, len(LOGICAL_LEVELS)), dtype=float)

    for index, scale in enumerate(scales):
        corrections = ConditionedMultitoneCorrections(
            d_lambda=(float(scale), float(scale), 0.0, 0.0),
            d_alpha=(0.0, 0.0, 0.0, 0.0),
            d_omega_rad_s=(0.0, 0.0, 0.0, 0.0),
        )
        reduced = run_conditioned_multitone_validation(
            model,
            targets,
            run_config,
            corrections=corrections,
            simulation_mode="reduced",
        )
        reduced_fidelity[index] = float(reduced.weighted_mean_fidelity)
        t_mean, s_mean = target_and_spectator_means(reduced.sector_metrics, target_levels)
        target_mean[index] = t_mean
        spectator_mean[index] = s_mean
        per_sector[index, :] = [float(metric.fidelity) for metric in reduced.sector_metrics]

    return {
        "case_name": case_name,
        "chi_t": chi_t,
        "theta_over_pi": theta / np.pi,
        "scales": scales,
        "reduced_fidelity": reduced_fidelity,
        "target_mean": target_mean,
        "spectator_mean": spectator_mean,
        "per_sector": per_sector,
    }


def _random_initial_corrections(seed: int) -> ConditionedMultitoneCorrections:
    rng = np.random.default_rng(seed)
    return ConditionedMultitoneCorrections(
        d_lambda=tuple(float(x) for x in rng.uniform(-3.0, 3.0, size=len(LOGICAL_LEVELS))),
        d_alpha=tuple(float(x) for x in rng.uniform(-np.pi, np.pi, size=len(LOGICAL_LEVELS))),
        d_omega_rad_s=tuple(float(2.0 * np.pi * x) for x in rng.uniform(-1.5e6, 1.5e6, size=len(LOGICAL_LEVELS))),
    )


def run_multistart_correction_case() -> dict[str, object]:
    """Check whether broader waveform-side correction rescues a hard case."""
    case_name = "pair_adjacent"
    target_levels = CASE_SPECS[case_name]
    theta = float(np.pi)
    chi_t = 3.0
    seeds = (101, 202, 303)

    model = build_model()
    run_config = run_config_for_chi_t(model, chi_t)
    targets, operator, transfer = _target_operator(target_levels, theta)
    baseline_reduced = run_conditioned_multitone_validation(
        model,
        targets,
        run_config,
        simulation_mode="reduced",
    )
    baseline_strict = run_targeted_subspace_multitone_validation(
        model,
        targets,
        run_config,
        logical_levels=LOGICAL_LEVELS,
        target_operator=operator,
        transfer_set=transfer,
        objective_weights=OBJECTIVE_WEIGHTS,
    )

    optimization_config = ConditionedOptimizationConfig(
        parameters=("d_lambda", "d_alpha", "d_omega"),
        maxiter_stage1=4,
        maxiter_stage2=6,
        d_lambda_bounds=(-4.0, 4.0),
        d_alpha_bounds=(-np.pi, np.pi),
        d_omega_hz_bounds=(-2.0e6, 2.0e6),
    )

    seed_rows: list[dict[str, object]] = []
    best_result = None
    best_opt = None
    for seed in seeds:
        print(f"[multistart] pair_adjacent chiT=3 theta=pi seed={seed}")
        initial_corrections = _random_initial_corrections(seed)
        optimization = optimize_conditioned_multitone(
            model,
            targets,
            run_config,
            initial_corrections=initial_corrections,
            optimization_config=optimization_config,
            simulation_mode="reduced",
        )
        optimized = optimization.optimized_result
        row = {
            "seed": int(seed),
            "initial_weighted_mean_fidelity": float(optimization.initial_result.weighted_mean_fidelity),
            "optimized_weighted_mean_fidelity": float(optimized.weighted_mean_fidelity),
            "optimized_sector_fidelity": [float(metric.fidelity) for metric in optimized.sector_metrics],
            "history_length": int(len(optimization.history)),
            "optimized_d_lambda": [float(x) for x in optimization.optimized_corrections.d_lambda],
            "optimized_d_alpha": [float(x) for x in optimization.optimized_corrections.d_alpha],
            "optimized_d_omega_hz": [
                float(x / (2.0 * np.pi)) for x in optimization.optimized_corrections.d_omega_rad_s
            ],
        }
        seed_rows.append(row)
        score = float(optimized.weighted_mean_fidelity)
        if best_result is None or score > best_result:
            best_result = score
            best_opt = optimization

    assert best_opt is not None
    best_strict = run_targeted_subspace_multitone_validation(
        model,
        targets,
        run_config,
        corrections=best_opt.optimized_corrections,
        logical_levels=LOGICAL_LEVELS,
        target_operator=operator,
        transfer_set=transfer,
        objective_weights=OBJECTIVE_WEIGHTS,
    )

    return {
        "case_name": case_name,
        "chi_t": chi_t,
        "theta_over_pi": theta / np.pi,
        "baseline_reduced_fidelity": float(baseline_reduced.weighted_mean_fidelity),
        "baseline_strict_fidelity": float(baseline_strict.restricted_process_fidelity),
        "baseline_compiled_fidelity": float(baseline_strict.best_fit_restricted_process_fidelity),
        "seed_rows": seed_rows,
        "best_reduced_fidelity": float(best_opt.optimized_result.weighted_mean_fidelity),
        "best_strict_fidelity": float(best_strict.restricted_process_fidelity),
        "best_compiled_fidelity": float(best_strict.best_fit_restricted_process_fidelity),
        "best_state_transfer_mean": float(best_strict.state_transfer_fidelity_mean),
        "best_state_transfer_min": float(best_strict.state_transfer_fidelity_min),
    }


def _mean_f_leakage(compiled, drive_ops, chi_t: float) -> tuple[float, float]:
    model = build_model(n_tr=N_TR_QUTRIT)
    frame = build_frame(model)
    session = prepare_simulation(
        model,
        compiled,
        drive_ops,
        config=SimulationConfig(frame=frame, max_step=DT_S, store_states=False),
        e_ops={},
    )
    n_cav = int(model.n_cav)
    leakages: list[float] = []
    for qubit_level in (0, 1):
        for cavity_level in LOGICAL_LEVELS:
            result = session.run(model.basis_state(qubit_level, int(cavity_level)))
            vec = np.asarray(result.final_state.full(), dtype=np.complex128).reshape(-1)
            f_population = float(np.sum(np.abs(vec[2 * n_cav : 3 * n_cav]) ** 2))
            leakages.append(f_population)
    return float(np.mean(leakages)), float(np.max(leakages))


def run_qutrit_leakage_spotcheck() -> dict[str, object]:
    """Replay representative pulses on an n_tr=3 model."""
    specs = (
        ("pair_adjacent", 3.0, np.pi / 8.0),
        ("pair_adjacent", 3.0, np.pi / 2.0),
        ("pair_adjacent", 3.0, np.pi),
        ("triple_low", 3.0, np.pi),
    )
    model = build_model()
    rows: list[dict[str, object]] = []
    for case_name, chi_t, theta in specs:
        target_levels = CASE_SPECS[case_name]
        run_config = run_config_for_chi_t(model, chi_t)
        _, _, _, strict = _baseline_bundle(model, run_config, target_levels, float(theta))
        leak_mean, leak_max = _mean_f_leakage(strict.compiled, strict.waveform.drive_ops, chi_t)
        rows.append(
            {
                "case_name": case_name,
                "chi_t": float(chi_t),
                "theta_over_pi": float(theta / np.pi),
                "strict_fidelity": float(strict.restricted_process_fidelity),
                "compiled_fidelity": float(strict.best_fit_restricted_process_fidelity),
                "mean_f_leakage": leak_mean,
                "max_f_leakage": leak_max,
            }
        )
    return {"rows": rows}


def run_segmented_small_angle_check() -> dict[str, object]:
    """Check whether repeated small-angle steps can replace a bad direct pi gate."""
    chi_t = 3.0
    step_theta = float(np.pi / 8.0)
    target_theta = float(np.pi)
    n_steps = 8
    model = build_model()
    run_config = run_config_for_chi_t(model, chi_t)
    restricted_indices = np.asarray(logical_indices(model), dtype=int)

    rows: list[dict[str, object]] = []
    for case_name in ("pair_adjacent", "triple_low"):
        target_levels = CASE_SPECS[case_name]
        _, _, _, step_result = _baseline_bundle(model, run_config, target_levels, step_theta)
        phases = step_result.best_fit_logical_block_phase.phases_for_levels(LOGICAL_LEVELS)
        phase_matrix = np.asarray(
            logical_block_phase_op(
                phases,
                fock_levels=LOGICAL_LEVELS,
                cavity_dim=int(model.n_cav),
                qubit_dim=int(model.n_tr),
            ).full(),
            dtype=np.complex128,
        )
        corrected_full = phase_matrix @ np.asarray(step_result.full_operator, dtype=np.complex128)
        corrected_step = corrected_full[np.ix_(restricted_indices, restricted_indices)]
        segmented_operator = np.linalg.matrix_power(corrected_step, n_steps)

        _, target_operator, _, direct_result = _baseline_bundle(model, run_config, target_levels, target_theta)
        segmented_fidelity = process_fidelity(target_operator, segmented_operator)
        rows.append(
            {
                "case_name": case_name,
                "chi_t": chi_t,
                "step_theta_over_pi": step_theta / np.pi,
                "n_steps": int(n_steps),
                "small_step_strict_fidelity": float(step_result.restricted_process_fidelity),
                "small_step_compiled_fidelity": float(step_result.best_fit_restricted_process_fidelity),
                "direct_pi_strict_fidelity": float(direct_result.restricted_process_fidelity),
                "direct_pi_compiled_fidelity": float(direct_result.best_fit_restricted_process_fidelity),
                "segmented_operator_fidelity": float(segmented_fidelity),
                "note": "Segmented fidelity is computed by composing the corrected restricted logical operator and does not include additional qutrit leakage accumulation.",
            }
        )
    return {"rows": rows}


def build_summary(
    baseline: dict[str, object],
    amplitude_scan: dict[str, object],
    multistart: dict[str, object],
    qutrit: dict[str, object],
    segmented: dict[str, object],
) -> dict[str, object]:
    """Produce a compact human-readable summary."""
    case_names = list(baseline["case_names"])
    strict = np.asarray(baseline["strict_fidelity"], dtype=float)
    compiled = np.asarray(baseline["compiled_fidelity"], dtype=float)
    reduced = np.asarray(baseline["reduced_fidelity"], dtype=float)

    representative_rows = []
    chi_index = int(np.where(np.isclose(CHI_T_VALUES, 3.0))[0][0])
    for case_index, case_name in enumerate(case_names):
        representative_rows.append(
            {
                "case_name": case_name,
                "label": CASE_LABELS[case_name],
                "chi_t": 3.0,
                "theta_over_pi": [float(x) for x in THETA_VALUES / np.pi],
                "reduced_fidelity": [float(x) for x in reduced[case_index, chi_index, :]],
                "strict_fidelity": [float(x) for x in strict[case_index, chi_index, :]],
                "compiled_fidelity": [float(x) for x in compiled[case_index, chi_index, :]],
                "target_theta_sim_over_pi": [
                    float(x / np.pi) for x in baseline["reduced_target_theta_mean"][case_index, chi_index, :]
                ],
                "target_theta_ratio": [
                    float(x) for x in baseline["reduced_target_theta_ratio"][case_index, chi_index, :]
                ],
            }
        )

    max_gain_index = np.unravel_index(np.argmax(compiled - strict), strict.shape)
    max_gain = float((compiled - strict)[max_gain_index])
    best_case = case_names[int(max_gain_index[0])]

    return {
        "study": "simultaneous_multitone_sqr_design",
        "angles_over_pi": [float(x) for x in THETA_VALUES / np.pi],
        "chi_t_values": [float(x) for x in CHI_T_VALUES],
        "representative_rows": representative_rows,
        "max_compiled_gain": {
            "gain": max_gain,
            "case_name": best_case,
            "chi_t": float(CHI_T_VALUES[int(max_gain_index[1])]),
            "theta_over_pi": float(THETA_VALUES[int(max_gain_index[2])] / np.pi),
        },
        "amplitude_scan_peak_reduced_fidelity": float(np.max(np.asarray(amplitude_scan["reduced_fidelity"], dtype=float))),
        "multistart_best_reduced_fidelity": float(multistart["best_reduced_fidelity"]),
        "multistart_best_strict_fidelity": float(multistart["best_strict_fidelity"]),
        "qutrit_spotcheck_max_leakage": float(
            np.max([float(row["max_f_leakage"]) for row in qutrit["rows"]], initial=0.0)
        ),
        "segmented_rows": segmented["rows"],
        "main_takeaway": (
            "For the common Gaussian multitone family studied here, the simultaneous waveform acts almost like identity on each branch: "
            "target-branch simulated angles stay at O(10^-4) rad even when the requested rotation is O(1) rad. "
            "Small-angle fidelity is therefore mostly trivial identity overlap, not successful simultaneous control. "
            "Cavity-only compilation gives only modest gains, and simple waveform corrections do not rescue the hard regime."
        ),
    }


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    baseline = run_baseline_grid()
    amplitude_scan = run_target_amplitude_scan()
    multistart = run_multistart_correction_case()
    qutrit = run_qutrit_leakage_spotcheck()
    segmented = run_segmented_small_angle_check()
    summary = build_summary(baseline, amplitude_scan, multistart, qutrit, segmented)

    np.savez_compressed(
        OUT_NPZ,
        case_names=np.asarray(baseline["case_names"], dtype=object),
        case_labels=np.asarray([CASE_LABELS[name] for name in baseline["case_names"]], dtype=object),
        chi_t_values=np.asarray(CHI_T_VALUES, dtype=float),
        theta_values=np.asarray(THETA_VALUES, dtype=float),
        theta_over_pi=np.asarray(THETA_VALUES / np.pi, dtype=float),
        reduced_fidelity=np.asarray(baseline["reduced_fidelity"], dtype=float),
        reduced_target_mean=np.asarray(baseline["reduced_target_mean"], dtype=float),
        reduced_spectator_mean=np.asarray(baseline["reduced_spectator_mean"], dtype=float),
        reduced_target_theta_mean=np.asarray(baseline["reduced_target_theta_mean"], dtype=float),
        reduced_target_theta_ratio=np.asarray(baseline["reduced_target_theta_ratio"], dtype=float),
        reduced_sector_fidelity=np.asarray(baseline["reduced_sector_fidelity"], dtype=float),
        strict_fidelity=np.asarray(baseline["strict_fidelity"], dtype=float),
        compiled_fidelity=np.asarray(baseline["compiled_fidelity"], dtype=float),
        compiled_gain=np.asarray(baseline["compiled_fidelity"], dtype=float)
        - np.asarray(baseline["strict_fidelity"], dtype=float),
        state_transfer_mean=np.asarray(baseline["state_transfer_mean"], dtype=float),
        state_transfer_min=np.asarray(baseline["state_transfer_min"], dtype=float),
        same_block_mean=np.asarray(baseline["same_block_mean"], dtype=float),
        leakage_mean=np.asarray(baseline["leakage_mean"], dtype=float),
        block_phase_rms=np.asarray(baseline["block_phase_rms"], dtype=float),
        strict_sector_fidelity=np.asarray(baseline["strict_sector_fidelity"], dtype=float),
        best_fit_phase=np.asarray(baseline["best_fit_phase"], dtype=float),
        amplitude_scan_scales=np.asarray(amplitude_scan["scales"], dtype=float),
        amplitude_scan_reduced_fidelity=np.asarray(amplitude_scan["reduced_fidelity"], dtype=float),
        amplitude_scan_target_mean=np.asarray(amplitude_scan["target_mean"], dtype=float),
        amplitude_scan_spectator_mean=np.asarray(amplitude_scan["spectator_mean"], dtype=float),
        amplitude_scan_per_sector=np.asarray(amplitude_scan["per_sector"], dtype=float),
    )

    json_payload = {
        "summary": summary,
        "amplitude_scan": {
            "case_name": amplitude_scan["case_name"],
            "chi_t": float(amplitude_scan["chi_t"]),
            "theta_over_pi": float(amplitude_scan["theta_over_pi"]),
            "scales": [float(x) for x in amplitude_scan["scales"]],
            "reduced_fidelity": [float(x) for x in amplitude_scan["reduced_fidelity"]],
            "target_mean": [float(x) for x in amplitude_scan["target_mean"]],
            "spectator_mean": [float(x) for x in amplitude_scan["spectator_mean"]],
        },
        "multistart": multistart,
        "qutrit_spotcheck": qutrit,
        "segmented_check": segmented,
    }
    OUT_JSON.write_text(json.dumps(json_payload, indent=2), encoding="utf-8")
    print(f"Saved {OUT_NPZ}")
    print(f"Saved {OUT_JSON}")


if __name__ == "__main__":
    main()
