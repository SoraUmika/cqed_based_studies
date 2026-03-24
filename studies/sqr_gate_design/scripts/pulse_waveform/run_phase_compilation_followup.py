"""
Phase-compilation follow-up for the SQR waveform-design study.

This extension targets the central post-report question:
how much of the strict-logical gap is recoverable by a cavity-only
Fock-dependent phase layer, as opposed to branch-local qubit-Z structure
or genuinely non-diagonal control error?

The script runs two complementary analyses:

1. Single-target SQR on an enlarged logical window (N=8) for Gaussian and
   cosine-squared baselines. This exposes the Fock-phase profile over many
   spectator branches and tests exact / polynomial cavity-phase correction.
2. A representative structured all-branch multitone case (Phase-B NM
   optimizer at short gate times), where branch-local Z freedom is strong
   but cavity-only phase compilation may fail.

Outputs:
    data/phase_compilation_results.npz
    data/phase_compilation_summary.json
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
STUDY_DIR = SCRIPT_DIR.parent
sys.path.insert(0, str(SCRIPT_DIR))

import runtime_compat  # noqa: F401

from common import N_FOCK, N_TR, build_frame, build_model, duration_from_chi_t
from phase_compilation_common import (
    allbranch_blocks,
    compiled_phase_metrics,
    full_operator_from_basis_outputs,
    single_target_blocks,
)
from run_allbranch_multitone import (
    build_allbranch_multitone_pulse,
    optimize_allbranch_multitone,
    simulate_and_extract as simulate_allbranch_basis,
)
from run_followup_multitone import (
    build_cosine_squared_pulse,
    build_single_tone_gaussian,
    simulate_and_extract as simulate_single_target,
)

DATA_DIR = STUDY_DIR / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)
SAVE_PATH = DATA_DIR / "phase_compilation_results.npz"
SUMMARY_PATH = DATA_DIR / "phase_compilation_summary.json"

SINGLE_TARGET_LOGICAL_N = 8
SINGLE_TARGET_BRANCH = 1
SINGLE_TARGET_THETA = np.pi
SINGLE_TARGET_PHI = 0.0
SINGLE_TARGET_CHI_T = np.array([1.0, 2.0, 3.0, 5.0], dtype=float)
SINGLE_TARGET_FAMILIES = (
    ("single_tone_gaussian", build_single_tone_gaussian),
    ("cosine_squared", build_cosine_squared_pulse),
)

ALLBRANCH_CHI_T = np.array([0.5, 1.0], dtype=float)
ALLBRANCH_THETA = np.pi
ALLBRANCH_PHI = 0.0


def round_float(value: float, digits: int = 6) -> float:
    return round(float(value), digits)


def save_results(payload: dict[str, np.ndarray | float | int | object]) -> None:
    np.savez_compressed(str(SAVE_PATH), **payload)
    print(f"  [saved {SAVE_PATH.name}]", flush=True)


def load_allbranch_grape_reference() -> dict[float, float]:
    """Load previously computed all-branch GRAPE fidelities for matching chiT."""
    path = DATA_DIR / "allbranch_multitone_results.npz"
    if not path.exists():
        return {}
    data = np.load(path, allow_pickle=True)
    return {
        float(chi_t): float(fid)
        for chi_t, fid in zip(data["chi_t_values"], data["grape_fid"], strict=True)
    }


def init_payload() -> dict[str, np.ndarray | float | int | object]:
    single_shape = (len(SINGLE_TARGET_FAMILIES), len(SINGLE_TARGET_CHI_T))
    single_phase_shape = single_shape + (SINGLE_TARGET_LOGICAL_N,)
    allbranch_shape = (len(ALLBRANCH_CHI_T),)
    allbranch_block_shape = (len(ALLBRANCH_CHI_T), N_FOCK)

    payload: dict[str, np.ndarray | float | int | object] = {
        "single_family_names": np.array([name for name, _ in SINGLE_TARGET_FAMILIES], dtype=object),
        "single_chi_t_values": SINGLE_TARGET_CHI_T.copy(),
        "single_logical_n": np.array(SINGLE_TARGET_LOGICAL_N, dtype=int),
        "single_target_branch": np.array(SINGLE_TARGET_BRANCH, dtype=int),
        "single_target_theta": np.array(SINGLE_TARGET_THETA, dtype=float),
        "single_target_phi": np.array(SINGLE_TARGET_PHI, dtype=float),
        "allbranch_chi_t_values": ALLBRANCH_CHI_T.copy(),
        "allbranch_logical_n": np.array(N_FOCK, dtype=int),
    }

    for key in (
        "single_raw_strict_fid",
        "single_raw_global_z_fid",
        "single_branch_local_z_relaxed_fid",
        "single_exact_cavity_compiled_fid",
        "single_linear_cavity_compiled_fid",
        "single_quadratic_cavity_compiled_fid",
        "single_cubic_cavity_compiled_fid",
        "single_linear_phase_fit_rms",
        "single_quadratic_phase_fit_rms",
        "single_cubic_phase_fit_rms",
        "single_same_block_population_mean",
        "single_same_block_population_min",
        "single_leakage_mean",
        "single_leakage_max",
        "single_pair_superposition_raw_mean",
        "single_pair_superposition_raw_min",
        "single_pair_superposition_compiled_mean",
        "single_pair_superposition_compiled_min",
    ):
        payload[key] = np.full(single_shape, np.nan, dtype=float)

    for key in (
        "single_exact_cavity_phases",
        "single_linear_phase_profile",
        "single_quadratic_phase_profile",
        "single_cubic_phase_profile",
        "single_branch_local_z_phases",
        "single_branch_local_z_fids",
    ):
        payload[key] = np.full(single_phase_shape, np.nan, dtype=float)

    for key in (
        "allbranch_raw_strict_fid",
        "allbranch_raw_global_z_fid",
        "allbranch_branch_local_z_relaxed_fid",
        "allbranch_exact_cavity_compiled_fid",
        "allbranch_linear_cavity_compiled_fid",
        "allbranch_quadratic_cavity_compiled_fid",
        "allbranch_linear_phase_fit_rms",
        "allbranch_same_block_population_mean",
        "allbranch_same_block_population_min",
        "allbranch_leakage_mean",
        "allbranch_leakage_max",
        "allbranch_grape_reference_fid",
    ):
        payload[key] = np.full(allbranch_shape, np.nan, dtype=float)

    for key in (
        "allbranch_exact_cavity_phases",
        "allbranch_linear_phase_profile",
        "allbranch_branch_local_z_phases",
        "allbranch_branch_local_z_fids",
        "allbranch_opt_amp_ratios",
        "allbranch_opt_drive_phases",
    ):
        payload[key] = np.full(allbranch_block_shape, np.nan, dtype=float)

    return payload


def run_single_target_scan(payload: dict[str, np.ndarray | float | int | object]) -> None:
    print("=" * 72)
    print("Single-target enlarged-window cavity-phase scan (logical N = 8)")
    print("=" * 72)

    model = build_model(n_cav=SINGLE_TARGET_LOGICAL_N + 2, n_tr=N_TR)
    frame = build_frame(model)
    target_blocks = single_target_blocks(
        SINGLE_TARGET_LOGICAL_N,
        SINGLE_TARGET_BRANCH,
        SINGLE_TARGET_THETA,
        SINGLE_TARGET_PHI,
    )

    for family_index, (family_name, builder) in enumerate(SINGLE_TARGET_FAMILIES):
        print(f"\nFamily: {family_name}", flush=True)
        for chi_index, chi_t_value in enumerate(SINGLE_TARGET_CHI_T):
            duration_s = duration_from_chi_t(float(chi_t_value))
            pulses, drive_ops, total_duration_s = builder(
                model,
                frame,
                SINGLE_TARGET_BRANCH,
                SINGLE_TARGET_THETA,
                SINGLE_TARGET_PHI,
                duration_s,
            )
            full_operator, final_states = simulate_single_target(
                model,
                frame,
                pulses,
                drive_ops,
                SINGLE_TARGET_LOGICAL_N,
                total_duration_s,
            )
            metrics = compiled_phase_metrics(
                full_operator,
                final_states,
                model,
                SINGLE_TARGET_LOGICAL_N,
                target_blocks,
                coherence_stats=True,
            )

            for key, source in (
                ("single_raw_strict_fid", "raw_strict_fid"),
                ("single_raw_global_z_fid", "raw_global_z_fid"),
                ("single_branch_local_z_relaxed_fid", "branch_local_z_relaxed_fid"),
                ("single_exact_cavity_compiled_fid", "exact_cavity_compiled_fid"),
                ("single_linear_cavity_compiled_fid", "linear_cavity_compiled_fid"),
                ("single_quadratic_cavity_compiled_fid", "quadratic_cavity_compiled_fid"),
                ("single_cubic_cavity_compiled_fid", "cubic_cavity_compiled_fid"),
                ("single_linear_phase_fit_rms", "linear_phase_fit_rms"),
                ("single_quadratic_phase_fit_rms", "quadratic_phase_fit_rms"),
                ("single_cubic_phase_fit_rms", "cubic_phase_fit_rms"),
                ("single_same_block_population_mean", "same_block_population_mean"),
                ("single_same_block_population_min", "same_block_population_min"),
                ("single_leakage_mean", "leakage_mean"),
                ("single_leakage_max", "leakage_max"),
                ("single_pair_superposition_raw_mean", "pair_superposition_raw_mean"),
                ("single_pair_superposition_raw_min", "pair_superposition_raw_min"),
                ("single_pair_superposition_compiled_mean", "pair_superposition_compiled_mean"),
                ("single_pair_superposition_compiled_min", "pair_superposition_compiled_min"),
            ):
                payload[key][family_index, chi_index] = float(metrics[source])

            for key, source in (
                ("single_exact_cavity_phases", "exact_cavity_phases"),
                ("single_linear_phase_profile", "linear_phase_profile"),
                ("single_quadratic_phase_profile", "quadratic_phase_profile"),
                ("single_cubic_phase_profile", "cubic_phase_profile"),
                ("single_branch_local_z_phases", "branch_local_z_phases"),
                ("single_branch_local_z_fids", "branch_local_z_fids"),
            ):
                payload[key][family_index, chi_index] = np.asarray(metrics[source], dtype=float)

            print(
                f"  chiT/2pi={chi_t_value:>3.1f}: "
                f"F_raw={metrics['raw_global_z_fid']:.4f}  "
                f"F_diag={metrics['linear_cavity_compiled_fid']:.4f}  "
                f"F_branchZ={metrics['branch_local_z_relaxed_fid']:.4f}  "
                f"phase_rms={metrics['linear_phase_fit_rms']:.2e}",
                flush=True,
            )

        save_results(payload)


def run_allbranch_short_gate_cases(payload: dict[str, np.ndarray | float | int | object]) -> None:
    print("\n" + "=" * 72)
    print("All-branch structured multitone representative cases (Phase-B NM)")
    print("=" * 72)

    model = build_model()
    frame = build_frame(model)
    target_blocks = allbranch_blocks(N_FOCK, ALLBRANCH_THETA, ALLBRANCH_PHI)
    grape_reference = load_allbranch_grape_reference()

    for chi_index, chi_t_value in enumerate(ALLBRANCH_CHI_T):
        print(f"\nchiT/2pi = {chi_t_value:.1f}", flush=True)
        start = time.time()
        duration_s = duration_from_chi_t(float(chi_t_value))
        opt = optimize_allbranch_multitone(
            model,
            frame,
            duration_s,
            use_detuning=False,
            use_de=False,
            verbose=True,
        )
        pulses, drive_ops, total_duration_s = build_allbranch_multitone_pulse(
            model,
            frame,
            duration_s,
            opt["amp_ratios"],
            opt["phases"],
        )
        final_states = simulate_allbranch_basis(model, frame, pulses, drive_ops, total_duration_s)
        full_operator = full_operator_from_basis_outputs(final_states, model, N_FOCK)
        metrics = compiled_phase_metrics(
            full_operator,
            final_states,
            model,
            N_FOCK,
            target_blocks,
            coherence_stats=False,
        )

        for key, source in (
            ("allbranch_raw_strict_fid", "raw_strict_fid"),
            ("allbranch_raw_global_z_fid", "raw_global_z_fid"),
            ("allbranch_branch_local_z_relaxed_fid", "branch_local_z_relaxed_fid"),
            ("allbranch_exact_cavity_compiled_fid", "exact_cavity_compiled_fid"),
            ("allbranch_linear_cavity_compiled_fid", "linear_cavity_compiled_fid"),
            ("allbranch_quadratic_cavity_compiled_fid", "quadratic_cavity_compiled_fid"),
            ("allbranch_linear_phase_fit_rms", "linear_phase_fit_rms"),
            ("allbranch_same_block_population_mean", "same_block_population_mean"),
            ("allbranch_same_block_population_min", "same_block_population_min"),
            ("allbranch_leakage_mean", "leakage_mean"),
            ("allbranch_leakage_max", "leakage_max"),
        ):
            payload[key][chi_index] = float(metrics[source])

        for key, source in (
            ("allbranch_exact_cavity_phases", "exact_cavity_phases"),
            ("allbranch_linear_phase_profile", "linear_phase_profile"),
            ("allbranch_branch_local_z_phases", "branch_local_z_phases"),
            ("allbranch_branch_local_z_fids", "branch_local_z_fids"),
        ):
            payload[key][chi_index] = np.asarray(metrics[source], dtype=float)

        payload["allbranch_opt_amp_ratios"][chi_index] = np.asarray(opt["amp_ratios"], dtype=float)
        payload["allbranch_opt_drive_phases"][chi_index] = np.asarray(opt["phases"], dtype=float)
        payload["allbranch_grape_reference_fid"][chi_index] = grape_reference.get(float(chi_t_value), np.nan)

        print(
            f"  F_raw={metrics['raw_global_z_fid']:.4f}  "
            f"F_diag={metrics['exact_cavity_compiled_fid']:.4f}  "
            f"F_branchZ={metrics['branch_local_z_relaxed_fid']:.4f}  "
            f"F_GRAPE={payload['allbranch_grape_reference_fid'][chi_index]:.4f}  "
            f"elapsed={time.time() - start:.1f}s",
            flush=True,
        )
        save_results(payload)


def write_summary(payload: dict[str, np.ndarray | float | int | object]) -> None:
    """Write a compact JSON summary for README/report reuse."""
    family_names = list(payload["single_family_names"])
    chi_values = payload["single_chi_t_values"]
    rep_index = int(np.where(np.isclose(chi_values, 3.0))[0][0])
    summary = {
        "single_target_window": {
            "logical_n": int(payload["single_logical_n"]),
            "chi_t_values": [round_float(value, 3) for value in chi_values],
            "families": {},
        },
        "allbranch_structured_short_gate": {},
    }

    for family_index, family_name in enumerate(family_names):
        summary["single_target_window"]["families"][str(family_name)] = {
            "raw_global_z_fidelity": [
                round_float(value) for value in payload["single_raw_global_z_fid"][family_index]
            ],
            "linear_cavity_compiled_fidelity": [
                round_float(value) for value in payload["single_linear_cavity_compiled_fid"][family_index]
            ],
            "branch_local_z_relaxed_fidelity": [
                round_float(value) for value in payload["single_branch_local_z_relaxed_fid"][family_index]
            ],
            "linear_phase_rms": [
                round_float(value, 10) for value in payload["single_linear_phase_fit_rms"][family_index]
            ],
            "representative_chi_t_3": {
                "exact_cavity_phases": [
                    round_float(value) for value in payload["single_exact_cavity_phases"][family_index, rep_index]
                ],
                "pair_superposition_raw_mean": round_float(
                    payload["single_pair_superposition_raw_mean"][family_index, rep_index]
                ),
                "pair_superposition_compiled_mean": round_float(
                    payload["single_pair_superposition_compiled_mean"][family_index, rep_index]
                ),
                "pair_superposition_raw_min": round_float(
                    payload["single_pair_superposition_raw_min"][family_index, rep_index]
                ),
                "pair_superposition_compiled_min": round_float(
                    payload["single_pair_superposition_compiled_min"][family_index, rep_index]
                ),
            },
        }

    for chi_index, chi_t_value in enumerate(payload["allbranch_chi_t_values"]):
        summary["allbranch_structured_short_gate"][str(round_float(chi_t_value, 3))] = {
            "raw_global_z_fidelity": round_float(payload["allbranch_raw_global_z_fid"][chi_index]),
            "exact_cavity_compiled_fidelity": round_float(
                payload["allbranch_exact_cavity_compiled_fid"][chi_index]
            ),
            "branch_local_z_relaxed_fidelity": round_float(
                payload["allbranch_branch_local_z_relaxed_fid"][chi_index]
            ),
            "grape_reference_fidelity": round_float(payload["allbranch_grape_reference_fid"][chi_index]),
            "exact_cavity_phases": [
                round_float(value) for value in payload["allbranch_exact_cavity_phases"][chi_index]
            ],
            "branch_local_z_phases": [
                round_float(value) for value in payload["allbranch_branch_local_z_phases"][chi_index]
            ],
        }

    with SUMMARY_PATH.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
    print(f"\nSummary written to {SUMMARY_PATH.name}", flush=True)


def main() -> None:
    payload = init_payload()
    start = time.time()

    run_single_target_scan(payload)
    run_allbranch_short_gate_cases(payload)
    write_summary(payload)

    total = time.time() - start
    print(f"\nPhase-compilation follow-up complete in {total:.1f}s")
    print(f"Results saved to {SAVE_PATH}")


if __name__ == "__main__":
    main()
