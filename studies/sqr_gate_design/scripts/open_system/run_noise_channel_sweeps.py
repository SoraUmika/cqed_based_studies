"""Run the A2 realistic-noise sweeps on the reusable two-mode SQR baselines.

This script intentionally separates two effects that were collapsed together in
the completed SQR study:

1. explicit multilevel transmon relaxation, and
2. storage thermal occupation.

The outputs are intended to feed A2.3 execution and A2.4 reporting directly.
"""

from __future__ import annotations

from dataclasses import dataclass
import sys
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import runtime_compat  # noqa: F401

from common import (
    CHI_T_VALUES,
    DATA_DIR,
    FAMILY_NAMES,
    KAPPA_STORAGE_DEFAULT,
    N_STORAGE_LOGICAL,
    TARGET_STORAGE_LEVEL,
    branch_average,
    build_basis_initial_states,
    build_family_pulse,
    build_frame,
    build_legacy_noise_spec,
    build_multilevel_noise_spec,
    build_session,
    build_target_state,
    build_two_mode_model,
    duration_from_chi_t,
    run_session_over_states,
    state_fidelity_pure_target,
)

MULTILEVEL_T1_CASES_US = np.array(
    [
        [20.0, 5.0],
        [30.0, 10.0],
        [50.0, 15.0],
    ],
    dtype=float,
)
THERMAL_OCCUPATIONS = np.array([0.00, 0.01, 0.02, 0.05, 0.10], dtype=float)
OUTPUT_PATH = DATA_DIR / "noise_channel_sweeps.npz"
CHECKPOINT_PATH = DATA_DIR / "noise_channel_sweeps.partial.npz"


@dataclass
class FamilyCase:
    model: object
    frame: object
    pulses: list
    drive_ops: dict[str, str]
    labels: list[tuple[int, int]]
    initial_states: list
    ideal_final: list
    target_states: list
    duration: float


def prepare_family_case(family_name: str, chi_t_2pi: float) -> FamilyCase:
    """Build and cache the ideal replay objects for one family and chiT point."""
    model = build_two_mode_model()
    frame = build_frame(model)
    duration = duration_from_chi_t(float(chi_t_2pi))
    pulses, drive_ops = build_family_pulse(
        model,
        frame,
        family_name,
        duration=duration,
        target_storage_level=TARGET_STORAGE_LEVEL,
        n_storage_levels=N_STORAGE_LOGICAL,
    )
    labels, initial_states = build_basis_initial_states(model, n_storage_levels=N_STORAGE_LOGICAL)
    ideal_session = build_session(model, frame, pulses, drive_ops, duration=duration, noise=None)
    ideal_final = run_session_over_states(ideal_session, initial_states)
    target_states = [
        build_target_state(
            model,
            storage_level=storage_level,
            qubit_level=qubit_level,
            target_storage_level=TARGET_STORAGE_LEVEL,
        )
        for qubit_level, storage_level in labels
    ]
    return FamilyCase(
        model=model,
        frame=frame,
        pulses=list(pulses),
        drive_ops=dict(drive_ops),
        labels=labels,
        initial_states=initial_states,
        ideal_final=ideal_final,
        target_states=target_states,
        duration=duration,
    )


def save_checkpoint(
    family_names: list[str],
    chi_t_values: np.ndarray,
    legacy_to_ideal: np.ndarray,
    legacy_to_target: np.ndarray,
    legacy_branch_to_ideal: np.ndarray,
    legacy_branch_to_target: np.ndarray,
    multilevel_to_ideal: np.ndarray,
    multilevel_to_target: np.ndarray,
    thermal_to_ideal: np.ndarray,
    thermal_to_target: np.ndarray,
    legacy_done: np.ndarray,
    multilevel_done: np.ndarray,
    thermal_done: np.ndarray,
) -> None:
    np.savez(
        CHECKPOINT_PATH,
        family_names=np.array(family_names, dtype=object),
        chi_t_values=chi_t_values,
        multilevel_t1_cases_us=MULTILEVEL_T1_CASES_US,
        thermal_occupations=THERMAL_OCCUPATIONS,
        kappa_storage_rad_s=KAPPA_STORAGE_DEFAULT,
        legacy_fidelity_to_ideal=legacy_to_ideal,
        legacy_fidelity_to_target=legacy_to_target,
        legacy_branch_to_ideal=legacy_branch_to_ideal,
        legacy_branch_to_target=legacy_branch_to_target,
        multilevel_fidelity_to_ideal=multilevel_to_ideal,
        multilevel_fidelity_to_target=multilevel_to_target,
        thermal_fidelity_to_ideal=thermal_to_ideal,
        thermal_fidelity_to_target=thermal_to_target,
        legacy_done=legacy_done,
        multilevel_done=multilevel_done,
        thermal_done=thermal_done,
    )


def load_checkpoint(
    family_names: list[str],
    chi_t_values: np.ndarray,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    legacy_shape = (len(family_names), len(chi_t_values))
    legacy_branch_shape = (len(family_names), len(chi_t_values), N_STORAGE_LOGICAL)
    multilevel_shape = (len(MULTILEVEL_T1_CASES_US), len(family_names), len(chi_t_values))
    thermal_shape = (len(THERMAL_OCCUPATIONS), len(family_names), len(chi_t_values))

    legacy_to_ideal = np.zeros(legacy_shape, dtype=float)
    legacy_to_target = np.zeros(legacy_shape, dtype=float)
    legacy_branch_to_ideal = np.zeros(legacy_branch_shape, dtype=float)
    legacy_branch_to_target = np.zeros(legacy_branch_shape, dtype=float)
    multilevel_to_ideal = np.zeros(multilevel_shape, dtype=float)
    multilevel_to_target = np.zeros(multilevel_shape, dtype=float)
    thermal_to_ideal = np.zeros(thermal_shape, dtype=float)
    thermal_to_target = np.zeros(thermal_shape, dtype=float)
    legacy_done = np.zeros(legacy_shape, dtype=bool)
    multilevel_done = np.zeros(multilevel_shape, dtype=bool)
    thermal_done = np.zeros(thermal_shape, dtype=bool)

    if not CHECKPOINT_PATH.exists():
        return (
            legacy_to_ideal,
            legacy_to_target,
            legacy_branch_to_ideal,
            legacy_branch_to_target,
            multilevel_to_ideal,
            multilevel_to_target,
            thermal_to_ideal,
            thermal_to_target,
            legacy_done,
            multilevel_done,
            thermal_done,
            np.zeros(0, dtype=float),
        )

    payload = np.load(CHECKPOINT_PATH, allow_pickle=True)
    if tuple(payload["legacy_fidelity_to_ideal"].shape) != legacy_shape:
        return (
            legacy_to_ideal,
            legacy_to_target,
            legacy_branch_to_ideal,
            legacy_branch_to_target,
            multilevel_to_ideal,
            multilevel_to_target,
            thermal_to_ideal,
            thermal_to_target,
            legacy_done,
            multilevel_done,
            thermal_done,
            np.zeros(0, dtype=float),
        )

    legacy_to_ideal[...] = payload["legacy_fidelity_to_ideal"]
    legacy_to_target[...] = payload["legacy_fidelity_to_target"]
    legacy_branch_to_ideal[...] = payload["legacy_branch_to_ideal"]
    legacy_branch_to_target[...] = payload["legacy_branch_to_target"]
    multilevel_to_ideal[...] = payload["multilevel_fidelity_to_ideal"]
    multilevel_to_target[...] = payload["multilevel_fidelity_to_target"]
    thermal_to_ideal[...] = payload["thermal_fidelity_to_ideal"]
    thermal_to_target[...] = payload["thermal_fidelity_to_target"]
    legacy_done[...] = payload["legacy_done"]
    multilevel_done[...] = payload["multilevel_done"]
    thermal_done[...] = payload["thermal_done"]
    return (
        legacy_to_ideal,
        legacy_to_target,
        legacy_branch_to_ideal,
        legacy_branch_to_target,
        multilevel_to_ideal,
        multilevel_to_target,
        thermal_to_ideal,
        thermal_to_target,
        legacy_done,
        multilevel_done,
        thermal_done,
        np.array(payload.files, dtype=object),
    )


def evaluate_prepared_case(case: FamilyCase, noise_spec) -> dict[str, np.ndarray | float]:
    """Replay one prepared pulse family at one chiT point with a selected noise model."""
    noisy_session = build_session(
        case.model,
        case.frame,
        case.pulses,
        case.drive_ops,
        duration=case.duration,
        noise=noise_spec,
    )
    noisy_final = run_session_over_states(noisy_session, case.initial_states)

    ideal_ref_fid = np.zeros(len(case.labels), dtype=float)
    target_fid = np.zeros(len(case.labels), dtype=float)
    for idx, (ideal_state, target_state, noisy_state) in enumerate(
        zip(case.ideal_final, case.target_states, noisy_final)
    ):
        ideal_ref_fid[idx] = state_fidelity_pure_target(ideal_state, noisy_state)
        target_fid[idx] = state_fidelity_pure_target(target_state, noisy_state)

    return {
        "mean_fidelity_to_ideal": float(np.mean(ideal_ref_fid)),
        "mean_fidelity_to_target": float(np.mean(target_fid)),
        "branch_fidelity_to_ideal": branch_average(ideal_ref_fid, N_STORAGE_LOGICAL),
        "branch_fidelity_to_target": branch_average(target_fid, N_STORAGE_LOGICAL),
    }


def main() -> None:
    family_names = list(FAMILY_NAMES)
    chi_t_values = np.asarray(CHI_T_VALUES, dtype=float)

    (
        legacy_to_ideal,
        legacy_to_target,
        legacy_branch_to_ideal,
        legacy_branch_to_target,
        multilevel_to_ideal,
        multilevel_to_target,
        thermal_to_ideal,
        thermal_to_target,
        legacy_done,
        multilevel_done,
        thermal_done,
        resumed_files,
    ) = load_checkpoint(family_names, chi_t_values)
    if resumed_files.size > 0:
        print(f"Resuming from {CHECKPOINT_PATH}")

    case_cache: dict[tuple[str, float], FamilyCase] = {}

    def get_case(family_name: str, chi_t_value: float) -> FamilyCase:
        key = (family_name, float(chi_t_value))
        case = case_cache.get(key)
        if case is None:
            case = prepare_family_case(family_name, float(chi_t_value))
            case_cache[key] = case
        return case

    print("=" * 68)
    print("A2 realistic-noise sweeps")
    print("=" * 68)

    legacy_noise = build_legacy_noise_spec(
        kappa_storage=KAPPA_STORAGE_DEFAULT,
        nth_storage=0.0,
    )

    for family_index, family_name in enumerate(family_names):
        print(f"\n[legacy] family={family_name}")
        for chi_index, chi_t_value in enumerate(chi_t_values):
            if legacy_done[family_index, chi_index]:
                print(f"  chiT/2pi={chi_t_value:4.1f}  resumed")
                continue
            metrics = evaluate_prepared_case(get_case(family_name, float(chi_t_value)), legacy_noise)
            legacy_to_ideal[family_index, chi_index] = metrics["mean_fidelity_to_ideal"]
            legacy_to_target[family_index, chi_index] = metrics["mean_fidelity_to_target"]
            legacy_branch_to_ideal[family_index, chi_index] = metrics["branch_fidelity_to_ideal"]
            legacy_branch_to_target[family_index, chi_index] = metrics["branch_fidelity_to_target"]
            legacy_done[family_index, chi_index] = True
            save_checkpoint(
                family_names,
                chi_t_values,
                legacy_to_ideal,
                legacy_to_target,
                legacy_branch_to_ideal,
                legacy_branch_to_target,
                multilevel_to_ideal,
                multilevel_to_target,
                thermal_to_ideal,
                thermal_to_target,
                legacy_done,
                multilevel_done,
                thermal_done,
            )
            print(
                f"  chiT/2pi={chi_t_value:4.1f}  "
                f"F_to_ideal={metrics['mean_fidelity_to_ideal']:.6f}  "
                f"F_to_target={metrics['mean_fidelity_to_target']:.6f}"
            )

    for case_index, t1_case_us in enumerate(MULTILEVEL_T1_CASES_US):
        transmon_t1 = tuple(float(value) * 1.0e-6 for value in t1_case_us)
        noise_spec = build_multilevel_noise_spec(
            transmon_t1=transmon_t1,
            kappa_storage=KAPPA_STORAGE_DEFAULT,
            nth_storage=0.0,
            for_three_mode=False,
        )
        print(f"\n[multilevel] T1_ge={t1_case_us[0]:.1f} us, T1_fe={t1_case_us[1]:.1f} us")
        for family_index, family_name in enumerate(family_names):
            for chi_index, chi_t_value in enumerate(chi_t_values):
                if multilevel_done[case_index, family_index, chi_index]:
                    continue
                metrics = evaluate_prepared_case(get_case(family_name, float(chi_t_value)), noise_spec)
                multilevel_to_ideal[case_index, family_index, chi_index] = metrics["mean_fidelity_to_ideal"]
                multilevel_to_target[case_index, family_index, chi_index] = metrics["mean_fidelity_to_target"]
                multilevel_done[case_index, family_index, chi_index] = True
                save_checkpoint(
                    family_names,
                    chi_t_values,
                    legacy_to_ideal,
                    legacy_to_target,
                    legacy_branch_to_ideal,
                    legacy_branch_to_target,
                    multilevel_to_ideal,
                    multilevel_to_target,
                    thermal_to_ideal,
                    thermal_to_target,
                    legacy_done,
                    multilevel_done,
                    thermal_done,
                )

    for thermal_index, nth_storage in enumerate(THERMAL_OCCUPATIONS):
        noise_spec = build_multilevel_noise_spec(
            transmon_t1=(30.0e-6, 10.0e-6),
            kappa_storage=KAPPA_STORAGE_DEFAULT,
            nth_storage=float(nth_storage),
            for_three_mode=False,
        )
        print(f"\n[thermal] nth_storage={nth_storage:.2f}")
        for family_index, family_name in enumerate(family_names):
            for chi_index, chi_t_value in enumerate(chi_t_values):
                if thermal_done[thermal_index, family_index, chi_index]:
                    continue
                metrics = evaluate_prepared_case(get_case(family_name, float(chi_t_value)), noise_spec)
                thermal_to_ideal[thermal_index, family_index, chi_index] = metrics["mean_fidelity_to_ideal"]
                thermal_to_target[thermal_index, family_index, chi_index] = metrics["mean_fidelity_to_target"]
                thermal_done[thermal_index, family_index, chi_index] = True
                save_checkpoint(
                    family_names,
                    chi_t_values,
                    legacy_to_ideal,
                    legacy_to_target,
                    legacy_branch_to_ideal,
                    legacy_branch_to_target,
                    multilevel_to_ideal,
                    multilevel_to_target,
                    thermal_to_ideal,
                    thermal_to_target,
                    legacy_done,
                    multilevel_done,
                    thermal_done,
                )

    np.savez(
        OUTPUT_PATH,
        family_names=np.array(family_names, dtype=object),
        chi_t_values=chi_t_values,
        multilevel_t1_cases_us=MULTILEVEL_T1_CASES_US,
        thermal_occupations=THERMAL_OCCUPATIONS,
        kappa_storage_rad_s=KAPPA_STORAGE_DEFAULT,
        legacy_fidelity_to_ideal=legacy_to_ideal,
        legacy_fidelity_to_target=legacy_to_target,
        legacy_branch_to_ideal=legacy_branch_to_ideal,
        legacy_branch_to_target=legacy_branch_to_target,
        multilevel_fidelity_to_ideal=multilevel_to_ideal,
        multilevel_fidelity_to_target=multilevel_to_target,
        thermal_fidelity_to_ideal=thermal_to_ideal,
        thermal_fidelity_to_target=thermal_to_target,
    )
    CHECKPOINT_PATH.unlink(missing_ok=True)
    print(f"\nSaved {OUTPUT_PATH}")


if __name__ == "__main__":
    main()