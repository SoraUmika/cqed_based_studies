"""Simulate concurrent SQR and readout in the three-mode cqed_sim model."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import qutip as qt

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import runtime_compat  # noqa: F401

from common import (
    CHI_T_VALUES,
    DATA_DIR,
    FAMILY_NAMES,
    KAPPA_R,
    KAPPA_STORAGE_DEFAULT,
    N_STORAGE_LOGICAL,
    TARGET_STORAGE_LEVEL,
    TPHI_READOUT_DEFAULT,
    TPHI_STORAGE_DEFAULT,
    branch_average,
    build_basis_initial_states,
    build_family_pulse,
    build_frame,
    build_multilevel_noise_spec,
    build_readout_square_pulse,
    build_session,
    build_storage_superposition_state,
    build_target_state,
    build_three_mode_model,
    duration_from_chi_t,
    reduce_qubit_storage,
    run_session_over_states,
    state_fidelity_pure_target,
    storage_coherence,
)

REPRESENTATIVE_FAMILIES = ("square", "cosine_squared", "multitone_one_segment")
REPRESENTATIVE_CHI_T = np.array([1.5, 3.0, 5.0], dtype=float)
READOUT_AMPLITUDES_MHZ = np.array([0.5, 1.0, 2.5], dtype=float)
OUTPUT_PATH = DATA_DIR / "three_mode_readout_effects.npz"
CHECKPOINT_PATH = DATA_DIR / "three_mode_readout_effects.partial.npz"


def save_checkpoint(
    reduced_target_fidelity: np.ndarray,
    strict_full_fidelity: np.ndarray,
    storage_coherence_ratio: np.ndarray,
    branch_reduced_fidelity: np.ndarray,
    done_mask: np.ndarray,
) -> None:
    np.savez(
        CHECKPOINT_PATH,
        family_names=np.array(REPRESENTATIVE_FAMILIES, dtype=object),
        chi_t_values=REPRESENTATIVE_CHI_T,
        readout_amplitudes_mhz=READOUT_AMPLITUDES_MHZ,
        reduced_target_fidelity=reduced_target_fidelity,
        strict_full_fidelity=strict_full_fidelity,
        branch_reduced_fidelity=branch_reduced_fidelity,
        storage_coherence_ratio=storage_coherence_ratio,
        done_mask=done_mask,
        kappa_storage_rad_s=KAPPA_STORAGE_DEFAULT,
        kappa_readout_rad_s=KAPPA_R,
    )


def load_checkpoint() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    shape = (len(REPRESENTATIVE_FAMILIES), len(REPRESENTATIVE_CHI_T), len(READOUT_AMPLITUDES_MHZ))
    reduced_target_fidelity = np.zeros(shape, dtype=float)
    strict_full_fidelity = np.zeros(shape, dtype=float)
    storage_coherence_ratio = np.zeros(shape, dtype=float)
    branch_reduced_fidelity = np.zeros(shape + (N_STORAGE_LOGICAL,), dtype=float)
    done_mask = np.zeros(shape, dtype=bool)
    if not CHECKPOINT_PATH.exists():
        return reduced_target_fidelity, strict_full_fidelity, storage_coherence_ratio, branch_reduced_fidelity, done_mask

    payload = np.load(CHECKPOINT_PATH, allow_pickle=True)
    if tuple(payload["reduced_target_fidelity"].shape) != shape:
        return reduced_target_fidelity, strict_full_fidelity, storage_coherence_ratio, branch_reduced_fidelity, done_mask

    reduced_target_fidelity[...] = payload["reduced_target_fidelity"]
    strict_full_fidelity[...] = payload["strict_full_fidelity"]
    storage_coherence_ratio[...] = payload["storage_coherence_ratio"]
    branch_reduced_fidelity[...] = payload["branch_reduced_fidelity"]
    done_mask[...] = payload["done_mask"]
    return reduced_target_fidelity, strict_full_fidelity, storage_coherence_ratio, branch_reduced_fidelity, done_mask


def qs_target_state_ket(model, qubit_level: int, storage_level: int):
    target_qubit = 1 - int(qubit_level) if int(storage_level) == TARGET_STORAGE_LEVEL else int(qubit_level)
    return np.kron(
        np.eye(model.n_tr, dtype=np.complex128)[:, target_qubit],
        np.eye(model.n_storage, dtype=np.complex128)[:, storage_level],
    )


def main() -> None:
    model = build_three_mode_model()
    frame = build_frame(model)
    labels, initial_states = build_basis_initial_states(model, n_storage_levels=N_STORAGE_LOGICAL, readout_level=0)

    (
        reduced_target_fidelity,
        strict_full_fidelity,
        storage_coherence_ratio,
        branch_reduced_fidelity,
        done_mask,
    ) = load_checkpoint()
    if CHECKPOINT_PATH.exists():
        print(f"Resuming from {CHECKPOINT_PATH}")

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

    print("=" * 68)
    print("A2 three-mode concurrent-readout study")
    print("=" * 68)

    for family_index, family_name in enumerate(REPRESENTATIVE_FAMILIES):
        print(f"\nfamily={family_name}")
        for chi_index, chi_t_value in enumerate(REPRESENTATIVE_CHI_T):
            duration = duration_from_chi_t(float(chi_t_value))
            sqr_pulses, sqr_drive_ops = build_family_pulse(
                model,
                frame,
                family_name,
                duration=duration,
                target_storage_level=TARGET_STORAGE_LEVEL,
                n_storage_levels=N_STORAGE_LOGICAL,
                readout_level=0,
            )
            reference_superposition = build_storage_superposition_state(model, readout_level=0)
            reference_coherence = abs(storage_coherence(reference_superposition))

            for amp_index, readout_amp_mhz in enumerate(READOUT_AMPLITUDES_MHZ):
                if done_mask[family_index, chi_index, amp_index]:
                    print(f"  chiT/2pi={chi_t_value:3.1f}  amp={readout_amp_mhz:3.1f} MHz  resumed")
                    continue
                readout_pulses, readout_drive_ops = build_readout_square_pulse(
                    model,
                    frame,
                    duration=duration,
                    amplitude=2.0 * np.pi * float(readout_amp_mhz) * 1.0e6,
                    storage_level=0,
                    readout_level=0,
                )
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
                )
                coherence_session = build_session(
                    model,
                    frame,
                    readout_pulses,
                    readout_drive_ops,
                    duration=duration,
                    noise=noise_spec,
                )
                final_states = run_session_over_states(session, initial_states)

                fidelity_to_qs_target = np.zeros(len(labels), dtype=float)
                fidelity_to_full_target = np.zeros(len(labels), dtype=float)
                for state_index, ((qubit_level, storage_level), final_state) in enumerate(zip(labels, final_states)):
                    full_target = build_target_state(
                        model,
                        storage_level=storage_level,
                        qubit_level=qubit_level,
                        target_storage_level=TARGET_STORAGE_LEVEL,
                        readout_level=0,
                    )
                    fidelity_to_full_target[state_index] = state_fidelity_pure_target(full_target, final_state)

                    reduced_state = reduce_qubit_storage(final_state)
                    qs_target_vec = qs_target_state_ket(model, qubit_level, storage_level)
                    qs_target = qt.Qobj(
                        qs_target_vec.reshape(-1, 1),
                        dims=[[model.n_tr, model.n_storage], [1, 1]],
                    )
                    fidelity_to_qs_target[state_index] = state_fidelity_pure_target(qs_target, reduced_state)

                reduced_target_fidelity[family_index, chi_index, amp_index] = float(np.mean(fidelity_to_qs_target))
                strict_full_fidelity[family_index, chi_index, amp_index] = float(np.mean(fidelity_to_full_target))
                branch_reduced_fidelity[family_index, chi_index, amp_index] = branch_average(
                    fidelity_to_qs_target,
                    N_STORAGE_LOGICAL,
                )

                noisy_superposition = coherence_session.run(reference_superposition).final_state
                coherence_ratio = 0.0
                if reference_coherence > 0.0:
                    coherence_ratio = abs(storage_coherence(noisy_superposition)) / reference_coherence
                storage_coherence_ratio[family_index, chi_index, amp_index] = float(coherence_ratio)
                done_mask[family_index, chi_index, amp_index] = True
                save_checkpoint(
                    reduced_target_fidelity,
                    strict_full_fidelity,
                    storage_coherence_ratio,
                    branch_reduced_fidelity,
                    done_mask,
                )

                print(
                    f"  chiT/2pi={chi_t_value:3.1f}  amp={readout_amp_mhz:3.1f} MHz  "
                    f"F_qs={reduced_target_fidelity[family_index, chi_index, amp_index]:.6f}  "
                    f"F_full={strict_full_fidelity[family_index, chi_index, amp_index]:.6f}  "
                    f"coh_ratio={storage_coherence_ratio[family_index, chi_index, amp_index]:.6f}"
                )

    np.savez(
        OUTPUT_PATH,
        family_names=np.array(REPRESENTATIVE_FAMILIES, dtype=object),
        chi_t_values=REPRESENTATIVE_CHI_T,
        readout_amplitudes_mhz=READOUT_AMPLITUDES_MHZ,
        reduced_target_fidelity=reduced_target_fidelity,
        strict_full_fidelity=strict_full_fidelity,
        branch_reduced_fidelity=branch_reduced_fidelity,
        storage_coherence_ratio=storage_coherence_ratio,
        kappa_storage_rad_s=KAPPA_STORAGE_DEFAULT,
        kappa_readout_rad_s=KAPPA_R,
    )
    CHECKPOINT_PATH.unlink(missing_ok=True)
    print(f"\nSaved {OUTPUT_PATH}")


if __name__ == "__main__":
    main()