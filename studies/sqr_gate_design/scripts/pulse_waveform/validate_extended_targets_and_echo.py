"""Data validation for the extended generalized-target and echoed-SQR study."""

from __future__ import annotations

from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_PATH = SCRIPT_DIR.parent / "data" / "extended_targets_results.npz"

PASS = "PASS"
FAIL = "FAIL"


def print_check(name: str, passed: bool, detail: str):
    status = PASS if passed else FAIL
    print(f"{name}: {status}")
    print(f"  {detail}")


def main():
    data = np.load(DATA_PATH, allow_pickle=True)

    prefixes = ["axis_scan", "branch_scan", "trunc_scan", "representative_scan"]
    all_ok = True

    print("Extended SQR study validation")
    print("=" * 60)

    # Check 1: no NaN/inf in all key metrics.
    finite_ok = True
    finite_count = 0
    for prefix in prefixes:
        for metric_name in (
            "branch_true_mean",
            "branch_cphase_mean",
            "joint_strict_fidelity",
            "joint_best_block_phase_fidelity",
            "state_transfer_mean",
            "state_transfer_min",
            "same_block_population_mean",
            "same_block_population_min",
            "leakage_mean",
            "leakage_max",
            "spectator_phase_spread",
            "spectator_max_transverse",
            "block_global_phase_spread",
        ):
            values = np.asarray(data[f"{prefix}_{metric_name}"])
            finite_ok &= np.all(np.isfinite(values))
            finite_count += values.size
    print_check("Check 1", finite_ok, f"Verified finiteness across {finite_count} stored metric values.")
    all_ok &= finite_ok

    # Check 2: cavity-block-phase-relaxed logical fidelity must exceed the strict logical fidelity.
    relaxed_ok = True
    worst_gap = 0.0
    for prefix in prefixes:
        strict = np.asarray(data[f"{prefix}_joint_strict_fidelity"], dtype=float)
        relaxed = np.asarray(data[f"{prefix}_joint_best_block_phase_fidelity"], dtype=float)
        gap = relaxed - strict
        worst_gap = min(worst_gap, float(np.min(gap)))
        relaxed_ok &= np.all(gap >= -1.0e-10)
    print_check("Check 2", relaxed_ok, f"Worst relaxed-minus-strict gap = {worst_gap:.2e}.")
    all_ok &= relaxed_ok

    # Check 3: conditional-phase branch metric must not be worse than the true-SQR branch metric.
    cphase_ok = True
    worst_cphase_gap = 0.0
    for prefix in prefixes:
        true_values = np.asarray(data[f"{prefix}_branch_true_mean"], dtype=float)
        cphase_values = np.asarray(data[f"{prefix}_branch_cphase_mean"], dtype=float)
        gap = cphase_values - true_values
        worst_cphase_gap = min(worst_cphase_gap, float(np.min(gap)))
        cphase_ok &= np.all(gap >= -1.0e-10)
    print_check("Check 3", cphase_ok, f"Worst cphase-minus-true gap = {worst_cphase_gap:.2e}.")
    all_ok &= cphase_ok

    # Check 4: ordinary single-tone should preserve x/y symmetry in the non-echo generalized scan.
    families = list(data["family_names"])
    family_index = families.index("single_tone_gaussian")
    axis_labels = list(data["axis_labels"])
    x_index = axis_labels.index("X")
    y_index = axis_labels.index("Y")
    x_values = np.asarray(data["axis_scan_joint_strict_fidelity"][family_index, :, x_index, :], dtype=float)
    y_values = np.asarray(data["axis_scan_joint_strict_fidelity"][family_index, :, y_index, :], dtype=float)
    symmetry_delta = np.max(np.abs(x_values - y_values))
    symmetry_ok = symmetry_delta < 5.0e-8
    print_check("Check 4", symmetry_ok, f"Max X-vs-Y logical-fidelity delta = {symmetry_delta:.2e}.")
    all_ok &= symmetry_ok

    # Check 5: same-block population must remain high in the representative four-family comparison.
    same_block = np.asarray(data["representative_scan_same_block_population_mean"], dtype=float)
    same_block_min = float(np.min(same_block))
    same_block_ok = same_block_min > 0.95
    print_check("Check 5", same_block_ok, f"Minimum representative same-block population mean = {same_block_min:.6f}.")
    all_ok &= same_block_ok

    # Check 6: representative arrays must be populated for each family and chiT.
    representative_phase_shape = tuple(np.asarray(data["representative_best_fit_block_phases"]).shape)
    expected_shape = (
        len(families),
        len(np.asarray(data["chi_t_values"])),
        int(np.asarray(data["representative_logical_n"]).item()),
    )
    representative_ok = representative_phase_shape == expected_shape
    print_check("Check 6", representative_ok, f"Representative block-phase array shape = {representative_phase_shape}; expected {expected_shape}.")
    all_ok &= representative_ok

    print("=" * 60)
    if all_ok:
        print("ALL EXTENDED VALIDATION CHECKS PASSED")
    else:
        print("SOME EXTENDED VALIDATION CHECKS FAILED")


if __name__ == "__main__":
    main()
