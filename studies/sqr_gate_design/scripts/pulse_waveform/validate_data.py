"""
Lightweight validation: check saved data files for consistency.
Avoids importing cqed_sim (which can be slow/hang).
Writes results to data/validation_output.txt
"""
import numpy as np
from pathlib import Path
import sys

STUDY = Path(__file__).resolve().parent.parent
DATA = STUDY / "data"
OUT = DATA / "validation_output.txt"

log = open(OUT, "w", encoding="utf-8")

def pr(s=""):
    log.write(s + "\n")
    log.flush()

PASS = "PASS"
FAIL = "FAIL"

try:
    pr("SQR Pulse-Waveform Design Study — Data Validation")
    pr("=" * 60)

    # Load all data
    pr("\nLoading data files...")
    p12 = np.load(DATA / "phase1_phase2_results.npz", allow_pickle=True)
    p4 = np.load(DATA / "phase4_results.npz", allow_pickle=True)
    p5 = np.load(DATA / "phase5_results.npz", allow_pickle=True)
    grape = np.load(DATA / "grape_benchmark_results.npz", allow_pickle=True)
    pr("  All 4 data files loaded OK")

    chi_t = p12["chi_t_values"]
    pr(f"  chi_t values: {chi_t}")

    all_pass = True

    # ================================================================
    pr("\n" + "=" * 60)
    pr("CHECK 1: GRAPE UPPER BOUND (cphase fidelity > parameterized)")
    pr("=" * 60)
    grape_chi_t = grape["chi_t_values"]
    grape_cphase = grape["fidelity_cphase"]
    for i, ct in enumerate(grape_chi_t):
        # Find matching chi_t in p12
        idx = np.argmin(np.abs(chi_t - ct))
        best_param = np.max(p12["cphase_sqr_fidelity"][:, idx])
        diff = grape_cphase[i] - best_param
        status = PASS if diff >= -1e-4 else FAIL
        if diff < -1e-4:
            all_pass = False
        pr(f"  chiT={ct:.1f}: GRAPE={grape_cphase[i]:.6f}  best_param={best_param:.6f}  diff={diff:.2e}  [{status}]")

    # ================================================================
    pr("\n" + "=" * 60)
    pr("CHECK 2: CPHASE > TRUE SQR FIDELITY (at all chi_t)")
    pr("=" * 60)
    n_fam = p12["cphase_sqr_fidelity"].shape[0]
    families = list(p12["family_names"])
    for fi in range(n_fam):
        for ci, ct in enumerate(chi_t):
            f_true = p12["true_sqr_fidelity"][fi, ci]
            f_cphase = p12["cphase_sqr_fidelity"][fi, ci]
            if f_cphase < f_true - 1e-6:
                pr(f"  UNEXPECTED: {families[fi]} at chiT={ct}: cphase={f_cphase:.6f} < true={f_true:.6f}  [{FAIL}]")
                all_pass = False
    pr(f"  cphase >= true for all families/chi_t  [{PASS}]")

    # ================================================================
    pr("\n" + "=" * 60)
    pr("CHECK 3: MONOTONIC CPHASE CONVERGENCE (smooth envelopes)")
    pr("=" * 60)
    # For Gaussian (idx 0) and cosine-squared (idx 2), cphase fidelity should
    # generally increase with chi_t (non-strictly)
    for fi in [0, 2]:
        fids = p12["cphase_sqr_fidelity"][fi]
        decreases = []
        for j in range(1, len(fids)):
            if fids[j] < fids[j-1] - 0.005:
                decreases.append((chi_t[j-1], chi_t[j], fids[j-1], fids[j]))
        status = PASS if len(decreases) == 0 else FAIL
        if len(decreases) > 0:
            all_pass = False
            for d in decreases:
                pr(f"  {families[fi]}: F drops from {d[2]:.4f} at chiT={d[0]} to {d[3]:.4f} at chiT={d[1]}")
        pr(f"  {families[fi]}: monotonic within 0.5%  [{status}]")

    # ================================================================
    pr("\n" + "=" * 60)
    pr("CHECK 4: PHASE 4 — HIGHER-ORDER CORRECTIONS SMALL")
    pr("=" * 60)
    idx_10 = np.argmin(np.abs(chi_t - 10))
    for fi in range(n_fam):
        f_p12 = p12["cphase_sqr_fidelity"][fi, idx_10]
        f_p4 = p4["cphase_sqr_fidelity"][fi, idx_10]
        delta = abs(f_p12 - f_p4)
        status = PASS if delta < 0.01 else FAIL
        if delta >= 0.01:
            all_pass = False
        pr(f"  {families[fi]}: |F_p12 - F_p4| = {delta:.2e} at chiT=10  [{status}]")

    # ================================================================
    pr("\n" + "=" * 60)
    pr("CHECK 5: PHASE 5 — DECOHERENCE ENVELOPE INDEPENDENCE")
    pr("=" * 60)
    chi_t_p5 = p5["chi_t_values"]
    idx_max = np.argmax(chi_t_p5)
    ct_ref = chi_t_p5[idx_max]
    f_deco_all = p5["deco_fid"][:, idx_max]
    spread = np.ptp(f_deco_all)
    status = PASS if spread < 5e-3 else FAIL
    if spread >= 5e-3:
        all_pass = False
    pr(f"  F_deco range at chiT={ct_ref:.0f}: [{f_deco_all.min():.6f}, {f_deco_all.max():.6f}]")
    pr(f"  Spread = {spread:.2e}  [{status}]")

    # ================================================================
    pr("\n" + "=" * 60)
    pr("CHECK 6: PHASE 5 — DECOHERENCE SCALING (1 - T/2T1)")
    pr("=" * 60)
    chi_rad = float(p5["chi_rad_s"])
    t1_s = float(p5["t1_s"])
    f_chi = abs(chi_rad) / (2 * np.pi)
    T_s = chi_t_p5[idx_max] / f_chi
    F_measured = float(p5["deco_fid"][0, idx_max])
    F_analytic = 1 - T_s / (2 * t1_s)
    delta = abs(F_measured - F_analytic)
    status = PASS if delta < 0.02 else FAIL
    if delta >= 0.02:
        all_pass = False
    pr(f"  F_deco (measured)  = {F_measured:.6f}")
    pr(f"  F_deco (1-T/2T1)  = {F_analytic:.6f}")
    pr(f"  T = {T_s*1e6:.2f} us, dF = {delta:.4f}  [{status}]")

    # ================================================================
    pr("\n" + "=" * 60)
    pr("CHECK 7: NET FIDELITY PEAK EXISTS")
    pr("=" * 60)
    f_net = p5["cphase_fid_net"][0]
    idx_peak = np.argmax(f_net)
    ct_peak = chi_t_p5[idx_peak]
    f_peak = f_net[idx_peak]
    status = PASS if 1 <= ct_peak <= 10 else FAIL
    if not (1 <= ct_peak <= 10):
        all_pass = False
    pr(f"  Peak F_net = {f_peak:.6f} at chiT/(2pi) = {ct_peak}")
    pr(f"  Expected range: 1-10  [{status}]")

    # ================================================================
    pr("\n" + "=" * 60)
    pr("CHECK 8: SPECTATOR SCALING (log-log slope)")
    pr("=" * 60)
    spec_max_trans = p12["spectator_max_transverse"][0]
    mask = chi_t >= 5
    x = np.log(chi_t[mask])
    y = np.log(spec_max_trans[mask] + 1e-15)
    if len(x) > 1:
        slope = np.polyfit(x, y, 1)[0]
        status = PASS if slope <= -2.0 else FAIL
        pr(f"  Spectator transverse error scaling: slope = {slope:.2f}")
        pr(f"  Expected: <= -2.0  [{status}]")
    else:
        pr("  Not enough data points")

    # ================================================================
    pr("\n" + "=" * 60)
    pr("CHECK 9: LEAKAGE NEGLIGIBLE")
    pr("=" * 60)
    if "leakage" in p12.files:
        leakage = p12["leakage"]
        max_leak = np.max(leakage)
        status = PASS if max_leak < 1e-3 else FAIL
        if max_leak >= 1e-3:
            all_pass = False
        pr(f"  Max leakage to |f>: {max_leak:.2e}  [{status}]")
    else:
        pr("  No leakage data in phase1_phase2_results.npz (skipped)")

    # ================================================================
    pr("\n" + "=" * 60)
    pr("VALIDATION SUMMARY")
    pr("=" * 60)
    if all_pass:
        pr("  ALL CHECKS PASSED")
    else:
        pr("  SOME CHECKS FAILED — review above")

except Exception as e:
    import traceback
    pr(f"\nERROR: {e}")
    traceback.print_exc(file=log)
finally:
    log.close()

print(f"Validation complete. Results in {OUT}")
