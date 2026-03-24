"""
Validation for the SQR phase-compilation follow-up.

This is the Step-4 validation gate for the new extension. It focuses on the
new claims introduced by the follow-up rather than re-running the original
baseline validation suite.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
STUDY_DIR = SCRIPT_DIR.parent
DATA_PATH = STUDY_DIR / "data" / "phase_compilation_results.npz"

PASS = "\033[92mPASS\033[0m"
FAIL = "\033[91mFAIL\033[0m"


def status(ok: bool) -> str:
    return PASS if ok else FAIL


def check_single_target_phase_compression(data) -> bool:
    print("=" * 64)
    print("CHECK 1: SINGLE-TARGET PHASE PROFILE COMPRESSION")
    print("=" * 64)
    rms = np.asarray(data["single_linear_phase_fit_rms"], dtype=float)
    worst = float(np.max(rms))
    ok = worst < 1e-3
    print(f"  Worst linear-fit RMS phase error = {worst:.3e} rad  [{status(ok)}]")
    return ok


def check_single_target_compilation_gain(data) -> bool:
    print("\n" + "=" * 64)
    print("CHECK 2: SINGLE-TARGET CAVITY COMPILATION IMPROVES STRICT LOGICAL FIDELITY")
    print("=" * 64)
    chi_t = np.asarray(data["single_chi_t_values"], dtype=float)
    raw = np.asarray(data["single_raw_global_z_fid"], dtype=float)
    compiled = np.asarray(data["single_linear_cavity_compiled_fid"], dtype=float)
    gain = compiled - raw
    min_gain = float(np.min(gain))
    practical_min_gain = float(np.min(gain[:, chi_t >= 2.0]))
    ok = min_gain > 0.015 and practical_min_gain > 0.02
    print(f"  Minimum gain over all families/chiT      = {min_gain:.4f}")
    print(f"  Minimum gain for practical chiT >= 2.0 = {practical_min_gain:.4f}  [{status(ok)}]")
    return ok


def check_single_target_coherence_gain(data) -> bool:
    print("\n" + "=" * 64)
    print("CHECK 3: SUPERPOSITION BENCHMARK IMPROVES AFTER COMPILATION")
    print("=" * 64)
    chi_t = np.asarray(data["single_chi_t_values"], dtype=float)
    rep_index = int(np.where(np.isclose(chi_t, 3.0))[0][0])
    raw_mean = np.asarray(data["single_pair_superposition_raw_mean"], dtype=float)[:, rep_index]
    compiled_mean = np.asarray(data["single_pair_superposition_compiled_mean"], dtype=float)[:, rep_index]
    raw_min = np.asarray(data["single_pair_superposition_raw_min"], dtype=float)[:, rep_index]
    compiled_min = np.asarray(data["single_pair_superposition_compiled_min"], dtype=float)[:, rep_index]
    mean_gain = float(np.min(compiled_mean - raw_mean))
    min_gain = float(np.min(compiled_min - raw_min))
    ok = mean_gain > 0.01 and min_gain > 0.15
    print(f"  Worst-case mean-fidelity gain at chiT=3 = {mean_gain:.4f}")
    print(f"  Worst-case min-fidelity gain at chiT=3  = {min_gain:.4f}  [{status(ok)}]")
    return ok


def check_allbranch_failure_mode(data) -> bool:
    print("\n" + "=" * 64)
    print("CHECK 4: ALL-BRANCH SHORT-GATE CASE IS NOT FIXED BY CAVITY-ONLY COMPILATION")
    print("=" * 64)
    cavity = np.asarray(data["allbranch_exact_cavity_compiled_fid"], dtype=float)
    branchz = np.asarray(data["allbranch_branch_local_z_relaxed_fid"], dtype=float)
    gap = branchz - cavity
    min_gap = float(np.min(gap))
    ok = min_gap > 0.4
    print(f"  Minimum gap F_branchZ - F_cavity = {min_gap:.4f}  [{status(ok)}]")
    return ok


def check_leakage(data) -> bool:
    print("\n" + "=" * 64)
    print("CHECK 5: SINGLE-TARGET ENLARGED-WINDOW LEAKAGE REMAINS NEGLIGIBLE")
    print("=" * 64)
    leakage = np.asarray(data["single_leakage_max"], dtype=float)
    worst = float(np.max(leakage))
    ok = worst < 1e-6
    print(f"  Worst leakage = {worst:.3e}  [{status(ok)}]")
    return ok


def main() -> None:
    if not DATA_PATH.exists():
        raise FileNotFoundError(f"Missing data file: {DATA_PATH}")

    data = np.load(DATA_PATH, allow_pickle=True)

    checks = [
        check_single_target_phase_compression(data),
        check_single_target_compilation_gain(data),
        check_single_target_coherence_gain(data),
        check_allbranch_failure_mode(data),
        check_leakage(data),
    ]

    print("\n" + "=" * 64)
    print("FOLLOW-UP VALIDATION SUMMARY")
    print("=" * 64)
    for idx, ok in enumerate(checks, start=1):
        print(f"  Check {idx}: {status(ok)}")
    if all(checks):
        print("\n  ALL FOLLOW-UP CHECKS PASSED")
    else:
        print("\n  FOLLOW-UP VALIDATION FAILED")


if __name__ == "__main__":
    main()
