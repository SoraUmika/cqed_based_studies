"""Iteration 3 — Experiment 2: Extended GRAPE sweep (500, 600, 800 ns).

Addresses:
  [P2 | LOW] Extend GRAPE sweep to longer durations
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
STUDY_DIR = SCRIPT_DIR.parent
DATA_DIR = STUDY_DIR / "data"

CQED_SIM = Path(
    r"C:\Users\dazzl\Box\Shyam Shankar Quantum Circuits Group"
    r"\Users\Users_JianJun\cQED_simulation"
)
if str(CQED_SIM) not in sys.path:
    sys.path.insert(0, str(CQED_SIM))

from cqed_sim import (
    DispersiveTransmonCavityModel, FrameSpec,
    ModelControlChannelSpec, PiecewiseConstantTimeGrid,
)
from cqed_sim.optimal_control import (
    GrapeConfig,
    GrapeSolver,
    LeakagePenalty as OCLeakagePenalty,
    UnitaryObjective as OCUnitaryObjective,
    build_control_problem_from_model,
)
from cqed_sim.unitary_synthesis import Subspace, subspace_unitary_fidelity
from cqed_sim.unitary_synthesis.targets import make_target

# ── Device parameters ─────────────────────────────────────────────────────────
TWO_PI = 2 * np.pi
OMEGA_Q = TWO_PI * 6.150e9
OMEGA_C = TWO_PI * 5.241e9
ALPHA = TWO_PI * (-255e6)
CHI = TWO_PI * (-2.84e6)
CHI_P = TWO_PI * (-21e3)
KERR = TWO_PI * (-28e3)

N_CAV = 8
N_TR = 2
FULL_DIM = N_TR * N_CAV

GRAPE_AMP_BOUND = TWO_PI * 50e6
GRAPE_MAXITER = 300
GRAPE_SEEDS = [17, 42, 73]

# ── Target and subspace ──────────────────────────────────────────────────────
U_target = make_target("cluster", n_match=1)
subspace = Subspace.custom(FULL_DIM, [0, 1, N_CAV, N_CAV + 1],
                           ["|g,0>", "|g,1>", "|e,0>", "|e,1>"])

# ── Build physical model ──────────────────────────────────────────────────────
model = DispersiveTransmonCavityModel(
    omega_q=OMEGA_Q, alpha=ALPHA,
    omega_c=OMEGA_C, chi=CHI, chi_higher=(CHI_P,), kerr=KERR,
    n_cav=N_CAV, n_tr=N_TR,
)
frame = FrameSpec(omega_q_frame=OMEGA_Q, omega_c_frame=OMEGA_C)

DATA_DIR.mkdir(parents=True, exist_ok=True)

# Prior results for context
prior_grape = {50: 0.6337, 100: 0.9494, 150: 0.9561, 200: 0.9966, 300: 0.9957, 400: 0.9990}

# ═══════════════════════════════════════════════════════════════════════════════
# Extended GRAPE sweep: 500, 600, 800 ns
# ═══════════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("  Iteration 3: Extended GRAPE Sweep (500, 600, 800 ns)")
print("=" * 70)

new_durations_ns = [500, 600, 800]
new_results = {}

for dur_ns in new_durations_ns:
    dur_s = dur_ns * 1e-9
    n_slices = max(10, int(dur_ns / 4))  # 4 ns per slice
    dt_s = dur_s / n_slices

    print(f"\n── GRAPE {dur_ns} ns ({n_slices} slices, dt={dt_s*1e9:.1f} ns) ──")

    best_fid = 0.0
    best_info = None
    t0 = time.perf_counter()

    for seed in GRAPE_SEEDS:
        try:
            problem = build_control_problem_from_model(
                model, frame=frame,
                time_grid=PiecewiseConstantTimeGrid.uniform(steps=n_slices, dt_s=dt_s),
                channel_specs=(
                    ModelControlChannelSpec(
                        name="storage", target="storage",
                        quadratures=("I", "Q"),
                        amplitude_bounds=(-GRAPE_AMP_BOUND, GRAPE_AMP_BOUND),
                        export_channel="storage",
                    ),
                    ModelControlChannelSpec(
                        name="qubit", target="qubit",
                        quadratures=("I", "Q"),
                        amplitude_bounds=(-GRAPE_AMP_BOUND, GRAPE_AMP_BOUND),
                        export_channel="qubit",
                    ),
                ),
                objectives=(
                    OCUnitaryObjective(
                        target_operator=U_target,
                        subspace=subspace,
                        ignore_global_phase=True,
                        name=f"cluster_{dur_ns}ns",
                    ),
                ),
                penalties=(
                    OCLeakagePenalty(weight=0.02, subspace=subspace),
                ),
            )
            result = GrapeSolver(GrapeConfig(
                maxiter=GRAPE_MAXITER, seed=seed, random_scale=0.3,
            )).solve(problem)

            fid = float(result.metrics.get("nominal_fidelity",
                        result.metrics.get("fidelity", 0.0)))
            print(f"    Seed {seed}: F={fid:.6f}")

            if fid > best_fid:
                best_fid = fid
                best_info = {
                    "seed": seed,
                    "converged": bool(result.metrics.get("converged", False)),
                    "iterations": int(result.metrics.get("iterations",
                                     result.metrics.get("n_iter", 0))),
                }
        except Exception as exc:
            print(f"    Seed {seed} failed: {exc}")

    elapsed = time.perf_counter() - t0
    entry = {
        "duration_ns": dur_ns,
        "n_slices": n_slices,
        "fidelity": best_fid,
        "elapsed_s": elapsed,
        "best_info": best_info,
    }
    new_results[f"GRAPE_{dur_ns}ns"] = entry
    print(f"  Best: F={best_fid:.6f} ({elapsed:.1f}s)")

# Combine with prior results
all_grape = dict(prior_grape)
for key, val in new_results.items():
    all_grape[val["duration_ns"]] = val["fidelity"]

# ═══════════════════════════════════════════════════════════════════════════════
# Save results
# ═══════════════════════════════════════════════════════════════════════════════
out = {
    "new_results": new_results,
    "all_grape_fidelities": {str(k): v for k, v in sorted(all_grape.items())},
    "prior_results": {str(k): v for k, v in prior_grape.items()},
}
out_path = DATA_DIR / "iteration3_grape_extension.json"
out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
print(f"\nResults saved to {out_path}")

print("\n" + "=" * 70)
print("  FULL GRAPE SWEEP SUMMARY")
print("=" * 70)
for dur_ns, fid in sorted(all_grape.items()):
    marker = " (NEW)" if dur_ns in [500, 600, 800] else ""
    print(f"  {dur_ns:4d} ns: F = {fid:.6f}{marker}")
print("=" * 70)
print("Done.")
