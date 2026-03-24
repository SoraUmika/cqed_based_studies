"""Validation checks for the open-system SQR deep-dive study."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import qutip as qt

from common import (
    DT,
    KAPPA_R,
    KAPPA_STORAGE_DEFAULT,
    N_CAV_TWO_MODE,
    N_READOUT_SIM,
    N_STORAGE_LOGICAL,
    N_STORAGE_SIM,
    N_TR,
    TARGET_STORAGE_LEVEL,
    TPHI_READOUT_DEFAULT,
    TPHI_STORAGE_DEFAULT,
    build_basis_initial_states,
    build_family_pulse,
    build_frame,
    build_multilevel_noise_spec,
    build_readout_square_pulse,
    build_session,
    build_storage_superposition_state,
    build_target_state,
    build_three_mode_model,
    build_two_mode_model,
    duration_from_chi_t,
    reduce_qubit_storage,
    run_session_over_states,
    state_fidelity_pure_target,
    storage_coherence,
)
from run_grape_noisy_replay import (
    MAXITER,
    N_SLICES,
    REPLAY_SUBSTEPS_PER_SLICE,
    extract_reported_fidelity,
    replay_archived_command_values,
    solve_grape_for_duration,
)


STUDY_DIR = Path(__file__).resolve().parents[2]
DATA_DIR = STUDY_DIR / "data" / Path(__file__).resolve().parent.name
CONVERGENCE_DIR = DATA_DIR / "convergence"
NOISE_PATH = DATA_DIR / "noise_channel_sweeps.npz"
PURCELL_PATH = DATA_DIR / "purcell_and_backaction.npz"
THREE_MODE_PATH = DATA_DIR / "three_mode_readout_effects.npz"
GRAPE_PATH = DATA_DIR / "grape_noisy_replay.npz"
CONVERGENCE_PATH = CONVERGENCE_DIR / "a2_convergence.npz"
CONVERGENCE_PARTIAL_PATH = CONVERGENCE_DIR / "a2_convergence.partial.npz"
CONVERGENCE_TARGET = 5.0e-4

THREE_MODE_PARTIAL_KEYS = (
    "three_mode_base_reduced",
    "three_mode_base_coherence",
    "three_mode_dt_reduced",
    "three_mode_dt_coherence",
    "three_mode_dims_reduced",
    "three_mode_dims_coherence",
)
GRAPE_COMPUTE_PARTIAL_KEYS = (
    "grape_base_noisy",
    "grape_dt_noisy",
    "grape_iter_objective",
)
GRAPE_METADATA_PARTIAL_KEYS = (
    "grape_archive_substeps",
    "grape_archive_noisy_ref",
    "grape_command_l2",
)
GRAPE_PARTIAL_KEYS = GRAPE_COMPUTE_PARTIAL_KEYS + GRAPE_METADATA_PARTIAL_KEYS

REPRESENTATIVE_T1_US = np.array([30.0, 10.0], dtype=float)
REPRESENTATIVE_TWO_MODE_FAMILY = "square"
REPRESENTATIVE_TWO_MODE_CHI_T = 2.0
REPRESENTATIVE_THREE_MODE_FAMILY = "square"
REPRESENTATIVE_THREE_MODE_CHI_T = 3.0
REPRESENTATIVE_READOUT_AMP_MHZ = 1.0
REPRESENTATIVE_GRAPE_CHI_T = 2.0


def load_payload(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"{path} not found.")
    return np.load(path, allow_pickle=True)


def load_partial_results(path: Path) -> dict[str, float]:
    if not path.exists():
        return {}
    payload = np.load(path, allow_pickle=True)
    return {key: float(np.asarray(payload[key], dtype=float)) for key in payload.files}


def save_partial_results(path: Path, partial_results: dict[str, float]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path, **partial_results)


def persist_partial_results(partial_results: dict[str, float]) -> None:
    if partial_results:
        save_partial_results(CONVERGENCE_PARTIAL_PATH, partial_results)
    elif CONVERGENCE_PARTIAL_PATH.exists():
        CONVERGENCE_PARTIAL_PATH.unlink()


def drop_partial_keys(partial_results: dict[str, float], keys: tuple[str, ...]) -> bool:
    changed = False
    for key in keys:
        if key in partial_results:
            partial_results.pop(key)
            changed = True
    return changed


def sanitize_partial_results(
    partial_results: dict[str, float],
    *,
    grape_archive_noisy: float,
    grape_archive_substeps: int,
    grape_command_values: np.ndarray,
) -> dict[str, float]:
    changed = False

    if any(partial_results.get(key, 0.0) > 1.05 for key in ("three_mode_base_coherence", "three_mode_dt_coherence", "three_mode_dims_coherence")):
        changed |= drop_partial_keys(partial_results, THREE_MODE_PARTIAL_KEYS)

    expected_command_l2 = float(np.linalg.norm(np.asarray(grape_command_values, dtype=float)))
    expected_grape_metadata = {
        "grape_archive_substeps": float(grape_archive_substeps),
        "grape_archive_noisy_ref": float(grape_archive_noisy),
        "grape_command_l2": expected_command_l2,
    }
    grape_metadata_stale = any(
        key in partial_results and abs(partial_results[key] - value) > 1.0e-12
        for key, value in expected_grape_metadata.items()
    )
    grape_metadata_missing = any(key not in partial_results for key in expected_grape_metadata)
    if grape_metadata_stale or grape_metadata_missing:
        changed |= drop_partial_keys(partial_results, GRAPE_PARTIAL_KEYS)
        partial_results.update(expected_grape_metadata)
        changed = True

    if changed:
        persist_partial_results(partial_results)
    return partial_results


def get_or_compute_scalar(
    partial_results: dict[str, float],
    key: str,
    description: str,
    compute,
) -> float:
    if key in partial_results:
        print(f"Loaded partial {description}: {partial_results[key]:.6e}")
        return float(partial_results[key])
    print(f"Running {description}...")
    value = float(compute())
    partial_results[key] = value
    save_partial_results(CONVERGENCE_PARTIAL_PATH, partial_results)
    print(f"Saved partial convergence state to {CONVERGENCE_PARTIAL_PATH}")
    return value


def assert_monotone_nonincreasing(name: str, values: np.ndarray, *, tol: float = 1.0e-6) -> None:
    diffs = np.diff(np.asarray(values, dtype=float))
    if np.any(diffs > tol):
        raise AssertionError(f"{name} is not monotone non-increasing within tolerance {tol}.")


def qs_target_state(model, qubit_level: int, storage_level: int) -> qt.Qobj:
    target_qubit = 1 - int(qubit_level) if int(storage_level) == TARGET_STORAGE_LEVEL else int(qubit_level)
    target_vec = np.kron(
        np.eye(model.n_tr, dtype=np.complex128)[:, target_qubit],
        np.eye(model.n_storage, dtype=np.complex128)[:, storage_level],
    )
    return qt.Qobj(target_vec.reshape(-1, 1), dims=[[model.n_tr, model.n_storage], [1, 1]])


def mean_two_mode_target_fidelity(
    family_name: str,
    chi_t_value: float,
    *,
    noise_spec,
    n_cav: int = N_CAV_TWO_MODE,
    n_tr: int = N_TR,
    dt: float = DT,
) -> float:
    model = build_two_mode_model(n_cav=n_cav, n_tr=n_tr)
    frame = build_frame(model)
    duration = duration_from_chi_t(float(chi_t_value))
    pulses, drive_ops = build_family_pulse(
        model,
        frame,
        family_name,
        duration=duration,
        target_storage_level=TARGET_STORAGE_LEVEL,
        n_storage_levels=N_STORAGE_LOGICAL,
    )
    labels, initial_states = build_basis_initial_states(model, n_storage_levels=N_STORAGE_LOGICAL)
    session = build_session(
        model,
        frame,
        pulses,
        drive_ops,
        duration=duration,
        noise=noise_spec,
        dt=dt,
    )
    final_states = run_session_over_states(session, initial_states)

    fidelity = np.zeros(len(labels), dtype=float)
    for idx, ((qubit_level, storage_level), final_state) in enumerate(zip(labels, final_states)):
        target_state = build_target_state(
            model,
            storage_level=storage_level,
            qubit_level=qubit_level,
            target_storage_level=TARGET_STORAGE_LEVEL,
        )
        fidelity[idx] = state_fidelity_pure_target(target_state, final_state)
    return float(np.mean(fidelity))


def three_mode_metrics(
    family_name: str,
    chi_t_value: float,
    readout_amp_mhz: float,
    *,
    n_storage: int = N_STORAGE_SIM,
    n_readout: int = N_READOUT_SIM,
    n_tr: int = N_TR,
    dt: float = DT,
) -> dict[str, float]:
    model = build_three_mode_model(n_storage=n_storage, n_readout=n_readout, n_tr=n_tr)
    frame = build_frame(model)
    duration = duration_from_chi_t(float(chi_t_value))
    labels, initial_states = build_basis_initial_states(model, n_storage_levels=N_STORAGE_LOGICAL, readout_level=0)
    sqr_pulses, sqr_drive_ops = build_family_pulse(
        model,
        frame,
        family_name,
        duration=duration,
        target_storage_level=TARGET_STORAGE_LEVEL,
        n_storage_levels=N_STORAGE_LOGICAL,
        readout_level=0,
    )
    readout_pulses, readout_drive_ops = build_readout_square_pulse(
        model,
        frame,
        duration=duration,
        amplitude=2.0 * np.pi * float(readout_amp_mhz) * 1.0e6,
        storage_level=0,
        readout_level=0,
    )
    noise_spec = build_multilevel_noise_spec(
        transmon_t1=(30.0e-6, 10.0e-6),
        kappa_storage=KAPPA_STORAGE_DEFAULT,
        nth_storage=0.0,
        tphi_storage=TPHI_STORAGE_DEFAULT,
        kappa_readout=KAPPA_R,
        nth_readout=0.0,
        tphi_readout=TPHI_READOUT_DEFAULT,
        for_three_mode=True,
    )

    reference_superposition = build_storage_superposition_state(model, readout_level=0)
    reference_coherence = abs(storage_coherence(reference_superposition))

    pulses = list(sqr_pulses) + list(readout_pulses)
    drive_ops = dict(sqr_drive_ops)
    drive_ops.update(readout_drive_ops)
    session = build_session(
        model,
        frame,
        pulses,
        drive_ops,
        duration=duration,
        noise=noise_spec,
        dt=dt,
    )
    coherence_session = build_session(
        model,
        frame,
        readout_pulses,
        readout_drive_ops,
        duration=duration,
        noise=noise_spec,
        dt=dt,
    )
    final_states = run_session_over_states(session, initial_states)
    reduced_target_fidelity = np.zeros(len(labels), dtype=float)
    for idx, ((qubit_level, storage_level), final_state) in enumerate(zip(labels, final_states)):
        reduced_state = reduce_qubit_storage(final_state)
        reduced_target_fidelity[idx] = state_fidelity_pure_target(
            qs_target_state(model, qubit_level, storage_level),
            reduced_state,
        )

    noisy_superposition = coherence_session.run(reference_superposition).final_state
    coherence_ratio = 0.0
    if reference_coherence > 0.0:
        coherence_ratio = abs(storage_coherence(noisy_superposition)) / reference_coherence
    return {
        "reduced_target_fidelity": float(np.mean(reduced_target_fidelity)),
        "coherence_ratio": float(coherence_ratio),
    }


def grape_objective_metric(
    chi_t_value: float,
    *,
    n_slices: int = N_SLICES,
    maxiter: int = MAXITER,
 ) -> float:
    duration = duration_from_chi_t(float(chi_t_value))
    _, _, _, result = solve_grape_for_duration(
        duration,
        n_slices=n_slices,
        maxiter=maxiter,
    )
    return float(extract_reported_fidelity(result))


def main() -> None:
    noise = load_payload(NOISE_PATH)
    purcell = load_payload(PURCELL_PATH)
    three_mode = load_payload(THREE_MODE_PATH)
    grape = load_payload(GRAPE_PATH)
    partial_results = load_partial_results(CONVERGENCE_PARTIAL_PATH)

    family_names = [str(name) for name in noise["family_names"]]
    chi_t_values = np.asarray(noise["chi_t_values"], dtype=float)
    legacy_to_ideal = np.asarray(noise["legacy_fidelity_to_ideal"], dtype=float)
    thermal_to_target = np.asarray(noise["thermal_fidelity_to_target"], dtype=float)
    thermal_occupations = np.asarray(noise["thermal_occupations"], dtype=float)

    print("A2 validation")
    print("=" * 40)

    for family_index, family_name in enumerate(family_names):
        assert_monotone_nonincreasing(
            f"legacy fidelity to ideal for {family_name}",
            legacy_to_ideal[family_index],
            tol=5.0e-4,
        )
    print("Legacy fidelity-to-ideal traces are monotone non-increasing in chiT.")

    if not np.all(purcell["purcell_t1_with_filter_s"] > purcell["purcell_t1_no_filter_s"]):
        raise AssertionError("Purcell filter did not improve the Purcell-limited T1 at every detuning.")
    print("Purcell filter improves the Purcell-limited T1 across the detuning sweep.")

    for family_index, family_name in enumerate(family_names):
        for chi_index, chi_t_value in enumerate(chi_t_values):
            assert_monotone_nonincreasing(
                f"thermal target fidelity for {family_name}, chiT={chi_t_value}",
                thermal_to_target[:, family_index, chi_index],
                tol=3.0e-3,
            )
    print("Thermal-occupation sweeps are monotone non-increasing within tolerance.")

    mean_three_mode_fidelity = np.mean(np.asarray(three_mode["reduced_target_fidelity"], dtype=float), axis=(0, 1))
    mean_three_mode_coherence = np.mean(np.asarray(three_mode["storage_coherence_ratio"], dtype=float), axis=(0, 1))
    assert_monotone_nonincreasing("mean three-mode reduced fidelity vs readout amplitude", mean_three_mode_fidelity, tol=5.0e-3)
    assert_monotone_nonincreasing("mean three-mode coherence ratio vs readout amplitude", mean_three_mode_coherence, tol=5.0e-3)
    print("Average three-mode readout backaction worsens monotonically with readout amplitude.")

    objective_fidelity = np.asarray(grape["objective_fidelity"], dtype=float)
    noisy_target_fidelity = np.asarray(grape["noisy_fidelity_to_target"], dtype=float)
    if np.any(objective_fidelity + 1.0e-6 < noisy_target_fidelity):
        raise AssertionError("Noisy GRAPE replay exceeded the reported closed-system objective fidelity.")
    print("GRAPE closed-system objective fidelities upper-bound the noisy replay fidelities.")

    grape_chi_t_values = np.asarray(grape["chi_t_values"], dtype=float)
    grape_index = int(np.argmin(np.abs(grape_chi_t_values - REPRESENTATIVE_GRAPE_CHI_T)))
    grape_archive_objective = float(objective_fidelity[grape_index])
    grape_archive_noisy = float(noisy_target_fidelity[grape_index])
    if "command_values" not in grape.files:
        raise AssertionError("GRAPE archive is missing command_values. Rerun run_grape_noisy_replay.py with the current script.")
    archived_command_values = np.asarray(grape["command_values"], dtype=float)
    grape_command_values = np.asarray(archived_command_values[grape_index], dtype=float)
    grape_archive_substeps = int(np.asarray(grape["replay_substeps_per_slice"], dtype=int))
    partial_results = sanitize_partial_results(
        partial_results,
        grape_archive_noisy=grape_archive_noisy,
        grape_archive_substeps=grape_archive_substeps,
        grape_command_values=grape_command_values,
    )

    t1_case_index = int(np.argmin(np.sum(np.abs(np.asarray(noise["multilevel_t1_cases_us"], dtype=float) - REPRESENTATIVE_T1_US), axis=1)))
    rep_family_index = family_names.index(REPRESENTATIVE_TWO_MODE_FAMILY)
    rep_chi_index = int(np.argmin(np.abs(chi_t_values - REPRESENTATIVE_TWO_MODE_CHI_T)))
    rep_noise = build_multilevel_noise_spec(
        transmon_t1=(30.0e-6, 10.0e-6),
        kappa_storage=KAPPA_STORAGE_DEFAULT,
        nth_storage=0.0,
        for_three_mode=False,
    )

    two_mode_archive = float(np.asarray(noise["multilevel_fidelity_to_target"], dtype=float)[t1_case_index, rep_family_index, rep_chi_index])
    two_mode_base = get_or_compute_scalar(
        partial_results,
        "two_mode_multilevel_base",
        "two-mode representative baseline",
        lambda: mean_two_mode_target_fidelity(
            REPRESENTATIVE_TWO_MODE_FAMILY,
            REPRESENTATIVE_TWO_MODE_CHI_T,
            noise_spec=rep_noise,
        ),
    )
    two_mode_dt = get_or_compute_scalar(
        partial_results,
        "two_mode_multilevel_dt",
        "two-mode dt refinement",
        lambda: mean_two_mode_target_fidelity(
            REPRESENTATIVE_TWO_MODE_FAMILY,
            REPRESENTATIVE_TWO_MODE_CHI_T,
            noise_spec=rep_noise,
            dt=DT / 2.0,
        ),
    )
    two_mode_dims = get_or_compute_scalar(
        partial_results,
        "two_mode_multilevel_dims",
        "two-mode truncation refinement",
        lambda: mean_two_mode_target_fidelity(
            REPRESENTATIVE_TWO_MODE_FAMILY,
            REPRESENTATIVE_TWO_MODE_CHI_T,
            noise_spec=rep_noise,
            n_cav=N_CAV_TWO_MODE + 2,
            n_tr=N_TR + 1,
        ),
    )

    three_mode_family_names = [str(name) for name in three_mode["family_names"]]
    three_mode_chi_t_values = np.asarray(three_mode["chi_t_values"], dtype=float)
    three_mode_amp_values = np.asarray(three_mode["readout_amplitudes_mhz"], dtype=float)
    three_mode_family_index = three_mode_family_names.index(REPRESENTATIVE_THREE_MODE_FAMILY)
    three_mode_chi_index = int(np.argmin(np.abs(three_mode_chi_t_values - REPRESENTATIVE_THREE_MODE_CHI_T)))
    three_mode_amp_index = int(np.argmin(np.abs(three_mode_amp_values - REPRESENTATIVE_READOUT_AMP_MHZ)))
    three_mode_archive = float(np.asarray(three_mode["reduced_target_fidelity"], dtype=float)[three_mode_family_index, three_mode_chi_index, three_mode_amp_index])
    if "three_mode_base_reduced" in partial_results and "three_mode_base_coherence" in partial_results:
        print(
            "Loaded partial three-mode representative baseline: "
            f"reduced={partial_results['three_mode_base_reduced']:.6e}, "
            f"coherence={partial_results['three_mode_base_coherence']:.6e}"
        )
        three_mode_base = {
            "reduced_target_fidelity": float(partial_results["three_mode_base_reduced"]),
            "coherence_ratio": float(partial_results["three_mode_base_coherence"]),
        }
    else:
        print("Running three-mode representative baseline...")
        three_mode_base = three_mode_metrics(
            REPRESENTATIVE_THREE_MODE_FAMILY,
            REPRESENTATIVE_THREE_MODE_CHI_T,
            REPRESENTATIVE_READOUT_AMP_MHZ,
        )
        partial_results["three_mode_base_reduced"] = float(three_mode_base["reduced_target_fidelity"])
        partial_results["three_mode_base_coherence"] = float(three_mode_base["coherence_ratio"])
        save_partial_results(CONVERGENCE_PARTIAL_PATH, partial_results)
        print(f"Saved partial convergence state to {CONVERGENCE_PARTIAL_PATH}")

    if "three_mode_dt_reduced" in partial_results and "three_mode_dt_coherence" in partial_results:
        print(
            "Loaded partial three-mode dt refinement: "
            f"reduced={partial_results['three_mode_dt_reduced']:.6e}, "
            f"coherence={partial_results['three_mode_dt_coherence']:.6e}"
        )
        three_mode_dt = {
            "reduced_target_fidelity": float(partial_results["three_mode_dt_reduced"]),
            "coherence_ratio": float(partial_results["three_mode_dt_coherence"]),
        }
    else:
        print("Running three-mode dt refinement...")
        three_mode_dt = three_mode_metrics(
            REPRESENTATIVE_THREE_MODE_FAMILY,
            REPRESENTATIVE_THREE_MODE_CHI_T,
            REPRESENTATIVE_READOUT_AMP_MHZ,
            dt=DT / 2.0,
        )
        partial_results["three_mode_dt_reduced"] = float(three_mode_dt["reduced_target_fidelity"])
        partial_results["three_mode_dt_coherence"] = float(three_mode_dt["coherence_ratio"])
        save_partial_results(CONVERGENCE_PARTIAL_PATH, partial_results)
        print(f"Saved partial convergence state to {CONVERGENCE_PARTIAL_PATH}")

    if "three_mode_dims_reduced" in partial_results and "three_mode_dims_coherence" in partial_results:
        print(
            "Loaded partial three-mode truncation refinement: "
            f"reduced={partial_results['three_mode_dims_reduced']:.6e}, "
            f"coherence={partial_results['three_mode_dims_coherence']:.6e}"
        )
        three_mode_dims = {
            "reduced_target_fidelity": float(partial_results["three_mode_dims_reduced"]),
            "coherence_ratio": float(partial_results["three_mode_dims_coherence"]),
        }
    else:
        print("Running three-mode truncation refinement...")
        three_mode_dims = three_mode_metrics(
            REPRESENTATIVE_THREE_MODE_FAMILY,
            REPRESENTATIVE_THREE_MODE_CHI_T,
            REPRESENTATIVE_READOUT_AMP_MHZ,
            n_storage=N_STORAGE_SIM + 2,
            n_readout=N_READOUT_SIM + 2,
        )
        partial_results["three_mode_dims_reduced"] = float(three_mode_dims["reduced_target_fidelity"])
        partial_results["three_mode_dims_coherence"] = float(three_mode_dims["coherence_ratio"])
        save_partial_results(CONVERGENCE_PARTIAL_PATH, partial_results)
        print(f"Saved partial convergence state to {CONVERGENCE_PARTIAL_PATH}")

    if "grape_base_noisy" in partial_results:
        print(
            "Loaded partial GRAPE baseline replay: "
            f"noisy={partial_results['grape_base_noisy']:.6e}"
        )
        grape_base = {
            "objective_fidelity": grape_archive_objective,
            "noisy_target_fidelity": float(partial_results["grape_base_noisy"]),
        }
    else:
        print("Running GRAPE baseline replay from archived control...")
        grape_base = {
            "objective_fidelity": grape_archive_objective,
            "noisy_target_fidelity": replay_archived_command_values(
                REPRESENTATIVE_GRAPE_CHI_T,
                grape_command_values,
                n_slices=N_SLICES,
                replay_substeps_per_slice=REPLAY_SUBSTEPS_PER_SLICE,
            ),
        }
        partial_results["grape_base_noisy"] = float(grape_base["noisy_target_fidelity"])
        persist_partial_results(partial_results)
        print(f"Saved partial convergence state to {CONVERGENCE_PARTIAL_PATH}")

    if "grape_dt_noisy" in partial_results:
        print(f"Loaded partial GRAPE dt refinement: noisy={partial_results['grape_dt_noisy']:.6e}")
        grape_dt = {
            "noisy_target_fidelity": float(partial_results["grape_dt_noisy"]),
        }
    else:
        print("Running GRAPE dt refinement replay from archived control...")
        grape_dt = {
            "noisy_target_fidelity": replay_archived_command_values(
                REPRESENTATIVE_GRAPE_CHI_T,
                grape_command_values,
                n_slices=N_SLICES,
                replay_substeps_per_slice=REPLAY_SUBSTEPS_PER_SLICE * 2,
            ),
        }
        partial_results["grape_dt_noisy"] = float(grape_dt["noisy_target_fidelity"])
        persist_partial_results(partial_results)
        print(f"Saved partial convergence state to {CONVERGENCE_PARTIAL_PATH}")

    if "grape_iter_objective" in partial_results:
        print(f"Loaded partial GRAPE iteration refinement: objective={partial_results['grape_iter_objective']:.6e}")
        grape_iter = {
            "objective_fidelity": float(partial_results["grape_iter_objective"]),
        }
    else:
        print("Running GRAPE iteration refinement...")
        grape_iter = {
            "objective_fidelity": grape_objective_metric(
                REPRESENTATIVE_GRAPE_CHI_T,
                n_slices=N_SLICES,
                maxiter=MAXITER + 100,
            ),
        }
        partial_results["grape_iter_objective"] = float(grape_iter["objective_fidelity"])
        persist_partial_results(partial_results)
        print(f"Saved partial convergence state to {CONVERGENCE_PARTIAL_PATH}")

    two_mode_archive_delta = abs(two_mode_base - two_mode_archive)
    two_mode_dt_delta = abs(two_mode_dt - two_mode_base)
    two_mode_dims_delta = abs(two_mode_dims - two_mode_base)
    three_mode_archive_delta = abs(three_mode_base["reduced_target_fidelity"] - three_mode_archive)
    three_mode_dt_delta = abs(three_mode_dt["reduced_target_fidelity"] - three_mode_base["reduced_target_fidelity"])
    three_mode_dims_delta = abs(three_mode_dims["reduced_target_fidelity"] - three_mode_base["reduced_target_fidelity"])
    grape_archive_objective_delta = 0.0
    grape_archive_noisy_delta = abs(grape_base["noisy_target_fidelity"] - grape_archive_noisy)
    grape_objective_delta = abs(grape_iter["objective_fidelity"] - grape_archive_objective)
    grape_noisy_delta = abs(grape_dt["noisy_target_fidelity"] - grape_base["noisy_target_fidelity"])

    print(f"Two-mode archive delta: {two_mode_archive_delta:.6e}")
    print(f"Two-mode dt delta: {two_mode_dt_delta:.6e}")
    print(f"Two-mode truncation delta: {two_mode_dims_delta:.6e}")
    print(f"Three-mode archive delta: {three_mode_archive_delta:.6e}")
    print(f"Three-mode dt delta: {three_mode_dt_delta:.6e}")
    print(f"Three-mode truncation delta: {three_mode_dims_delta:.6e}")
    print(f"GRAPE archive objective delta: {grape_archive_objective_delta:.6e}")
    print(f"GRAPE archive noisy replay delta: {grape_archive_noisy_delta:.6e}")
    print(f"GRAPE objective delta: {grape_objective_delta:.6e}")
    print(f"GRAPE noisy replay delta: {grape_noisy_delta:.6e}")

    if two_mode_archive_delta > CONVERGENCE_TARGET or two_mode_dt_delta > CONVERGENCE_TARGET or two_mode_dims_delta > CONVERGENCE_TARGET:
        raise AssertionError(f"Two-mode convergence check exceeded {CONVERGENCE_TARGET:.0e}.")
    if three_mode_archive_delta > CONVERGENCE_TARGET or three_mode_dt_delta > CONVERGENCE_TARGET or three_mode_dims_delta > CONVERGENCE_TARGET:
        raise AssertionError(f"Three-mode convergence check exceeded {CONVERGENCE_TARGET:.0e}.")
    if grape_archive_objective_delta > CONVERGENCE_TARGET or grape_archive_noisy_delta > CONVERGENCE_TARGET:
        raise AssertionError(f"GRAPE archive-consistency check exceeded {CONVERGENCE_TARGET:.0e}.")
    if grape_objective_delta > CONVERGENCE_TARGET or grape_noisy_delta > CONVERGENCE_TARGET:
        raise AssertionError(f"GRAPE convergence check exceeded {CONVERGENCE_TARGET:.0e}.")

    CONVERGENCE_DIR.mkdir(parents=True, exist_ok=True)
    np.savez(
        CONVERGENCE_PATH,
        two_mode_archive=two_mode_archive,
        two_mode_base=two_mode_base,
        two_mode_dt=two_mode_dt,
        two_mode_dims=two_mode_dims,
        two_mode_archive_delta=two_mode_archive_delta,
        two_mode_dt_delta=two_mode_dt_delta,
        two_mode_dims_delta=two_mode_dims_delta,
        three_mode_archive=three_mode_archive,
        three_mode_base_reduced=three_mode_base["reduced_target_fidelity"],
        three_mode_dt_reduced=three_mode_dt["reduced_target_fidelity"],
        three_mode_dims_reduced=three_mode_dims["reduced_target_fidelity"],
        three_mode_base_coherence=three_mode_base["coherence_ratio"],
        three_mode_dt_coherence=three_mode_dt["coherence_ratio"],
        three_mode_dims_coherence=three_mode_dims["coherence_ratio"],
        three_mode_archive_delta=three_mode_archive_delta,
        three_mode_dt_delta=three_mode_dt_delta,
        three_mode_dims_delta=three_mode_dims_delta,
        grape_archive_objective=grape_archive_objective,
        grape_archive_noisy=grape_archive_noisy,
        grape_base_objective=grape_base["objective_fidelity"],
        grape_base_noisy=grape_base["noisy_target_fidelity"],
        grape_dt_noisy=grape_dt["noisy_target_fidelity"],
        grape_iter_objective=grape_iter["objective_fidelity"],
        grape_archive_objective_delta=grape_archive_objective_delta,
        grape_archive_noisy_delta=grape_archive_noisy_delta,
        grape_objective_delta=grape_objective_delta,
        grape_noisy_delta=grape_noisy_delta,
    )
    if CONVERGENCE_PARTIAL_PATH.exists():
        CONVERGENCE_PARTIAL_PATH.unlink()
    print(f"Saved {CONVERGENCE_PATH}")
    print("Validation passed.")


if __name__ == "__main__":
    main()
