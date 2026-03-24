"""Validation checks for the simultaneous multitone SQR study."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from common import CASE_SPECS, CHI_T_VALUES, LOGICAL_LEVELS, THETA_VALUES, OBJECTIVE_WEIGHTS, build_model, run_config_for_chi_t
from cqed_sim.calibration import (
    build_block_rotation_target_operator,
    build_spanning_state_transfer_set,
    run_targeted_subspace_multitone_validation,
)


STUDY_DIR = Path(__file__).resolve().parents[2]
DATA_DIR = STUDY_DIR / "data" / Path(__file__).resolve().parent.name
DATA_PATH = DATA_DIR / "simultaneous_multitone_sqr_results.npz"
SUMMARY_PATH = DATA_DIR / "simultaneous_multitone_sqr_summary.json"


def load_payload():
    if not DATA_PATH.exists():
        raise FileNotFoundError(f"{DATA_PATH} not found. Run run_multitone_simultaneous_sqr_study.py first.")
    if not SUMMARY_PATH.exists():
        raise FileNotFoundError(f"{SUMMARY_PATH} not found. Run run_multitone_simultaneous_sqr_study.py first.")
    data = np.load(DATA_PATH, allow_pickle=True)
    summary = json.loads(SUMMARY_PATH.read_text(encoding="utf-8"))
    return data, summary


def assert_monotone_nonincreasing(name: str, values: np.ndarray, *, tol: float = 5.0e-4) -> None:
    diffs = np.diff(np.asarray(values, dtype=float))
    if np.any(diffs > tol):
        raise AssertionError(f"{name} is not monotone non-increasing within tolerance {tol}.")


def representative_validation(case_name: str, theta: float, *, dt_s: float, n_cav: int) -> dict[str, float]:
    target_levels = CASE_SPECS[case_name]
    model = build_model(n_cav=n_cav)
    run_config = run_config_for_chi_t(model, 3.0, dt_s=dt_s)
    from common import build_targets  # local import keeps this script lightweight

    targets = build_targets(target_levels, theta)
    operator = build_block_rotation_target_operator(targets, logical_levels=LOGICAL_LEVELS)
    transfer = build_spanning_state_transfer_set(operator)
    result = run_targeted_subspace_multitone_validation(
        model,
        targets,
        run_config,
        logical_levels=LOGICAL_LEVELS,
        target_operator=operator,
        transfer_set=transfer,
        objective_weights=OBJECTIVE_WEIGHTS,
    )
    return {
        "strict": float(result.restricted_process_fidelity),
        "compiled": float(result.best_fit_restricted_process_fidelity),
        "state_mean": float(result.state_transfer_fidelity_mean),
    }


def main() -> None:
    data, summary = load_payload()

    case_names = [str(x) for x in data["case_names"]]
    strict = np.asarray(data["strict_fidelity"], dtype=float)
    compiled = np.asarray(data["compiled_fidelity"], dtype=float)
    reduced = np.asarray(data["reduced_fidelity"], dtype=float)
    same_block = np.asarray(data["same_block_mean"], dtype=float)
    compiled_gain = np.asarray(data["compiled_gain"], dtype=float)

    print("Simultaneous multitone SQR validation")
    print("=" * 42)

    min_same_block = float(np.min(same_block))
    print(f"Minimum same-block population mean: {min_same_block:.9f}")
    if min_same_block < 0.99999:
        raise AssertionError("Same-block population dropped below 0.99999.")

    min_compiled_gain = float(np.min(compiled_gain))
    print(f"Minimum compiled gain: {min_compiled_gain:.6e}")
    if min_compiled_gain < -1.0e-5:
        raise AssertionError("Best-fit compiled fidelity fell below raw strict fidelity by more than numerical tolerance.")

    for case_index, case_name in enumerate(case_names):
        for chi_index, chi_t in enumerate(CHI_T_VALUES):
            assert_monotone_nonincreasing(
                f"reduced[{case_name}, chiT={chi_t}]",
                reduced[case_index, chi_index, :],
            )
            assert_monotone_nonincreasing(
                f"strict[{case_name}, chiT={chi_t}]",
                strict[case_index, chi_index, :],
                tol=1.0e-3,
            )
    print("Angle sweeps are monotone non-increasing within tolerance.")

    pair_base = representative_validation("pair_adjacent", float(np.pi / 4.0), dt_s=4.0e-9, n_cav=6)
    pair_dt = representative_validation("pair_adjacent", float(np.pi / 4.0), dt_s=2.0e-9, n_cav=6)
    pair_ncav = representative_validation("pair_adjacent", float(np.pi / 4.0), dt_s=4.0e-9, n_cav=7)

    triple_base = representative_validation("triple_low", float(np.pi), dt_s=4.0e-9, n_cav=6)
    triple_dt = representative_validation("triple_low", float(np.pi), dt_s=2.0e-9, n_cav=6)
    triple_ncav = representative_validation("triple_low", float(np.pi), dt_s=4.0e-9, n_cav=7)

    def max_delta(ref: dict[str, float], test: dict[str, float]) -> float:
        return max(abs(float(test[key]) - float(ref[key])) for key in ref)

    pair_dt_delta = max_delta(pair_base, pair_dt)
    pair_ncav_delta = max_delta(pair_base, pair_ncav)
    triple_dt_delta = max_delta(triple_base, triple_dt)
    triple_ncav_delta = max_delta(triple_base, triple_ncav)

    print(f"Pair chiT=3, theta=pi/4 dt refinement delta: {pair_dt_delta:.6f}")
    print(f"Pair chiT=3, theta=pi/4 n_cav refinement delta: {pair_ncav_delta:.6f}")
    print(f"Triple chiT=3, theta=pi dt refinement delta: {triple_dt_delta:.6f}")
    print(f"Triple chiT=3, theta=pi n_cav refinement delta: {triple_ncav_delta:.6f}")

    if pair_dt_delta > 5.0e-3 or triple_dt_delta > 5.0e-3:
        raise AssertionError("Time-step convergence check exceeded 5e-3.")
    if pair_ncav_delta > 5.0e-4 or triple_ncav_delta > 5.0e-4:
        raise AssertionError("Cavity truncation convergence check exceeded 5e-4.")

    qutrit_rows = summary["qutrit_spotcheck"]["rows"]
    max_qutrit_leakage = max(float(row["max_f_leakage"]) for row in qutrit_rows)
    print(f"Maximum qutrit replay leakage: {max_qutrit_leakage:.3e}")
    if max_qutrit_leakage > 1.0e-6:
        raise AssertionError("Representative qutrit replay leakage exceeded 1e-6.")

    print("Validation passed.")


if __name__ == "__main__":
    main()
