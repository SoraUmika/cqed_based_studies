"""Validation checks for the procedural readout study."""

from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

STUDY_DIR = Path(__file__).resolve().parents[2]
DATA_DIR = STUDY_DIR / "data" / Path(__file__).resolve().parent.name
if str(Path(__file__).resolve().parent) not in sys.path:
    sys.path.insert(0, str(Path(__file__).resolve().parent))

from readout_opt import DEFAULT_CONFIG, set_nulling_tail_kappa, solve_linear_response
from readout_opt.pulse_families import get_family
from readout_opt.simulate import evaluate_full_design


def assert_small(name: str, value: float, tol: float) -> None:
    status = "PASS" if value <= tol else "FAIL"
    print(f"[{status}] {name}: {value:.4e} (tol={tol:.4e})")
    if value > tol:
        raise SystemExit(1)


def main() -> None:
    summary_path = DATA_DIR / "study_summary.json"
    if not summary_path.exists():
        raise SystemExit("Run scripts/run_study.py first.")
    summary = json.loads(summary_path.read_text())

    set_nulling_tail_kappa(DEFAULT_CONFIG.kappa)
    params = np.array([0.20, 0.0, 0.15, 0.0, 1.0, 2.0, 1.0, 1.0, 0.0], dtype=float)
    design = get_family("nulling_tail").builder(
        params,
        DEFAULT_CONFIG.representative_duration,
        DEFAULT_CONFIG.dt,
        DEFAULT_CONFIG.amp_max,
        DEFAULT_CONFIG.chi,
    )
    tlist = np.arange(len(design.waveform) + 1, dtype=float) * DEFAULT_CONFIG.dt
    linear = solve_linear_response(design.waveform, tlist, kappa=DEFAULT_CONFIG.kappa, chi=DEFAULT_CONFIG.chi, delta_g=design.delta_g)
    assert_small("Nulling-tail final alpha_g", abs(linear.alpha_g[-1]), 1.0e-9)
    assert_small("Nulling-tail final alpha_e", abs(linear.alpha_e[-1]), 1.0e-9)

    convergence = summary["convergence"]
    assert_small("Fine-dt fidelity delta", float(convergence["fidelity_dt_delta"]), 3.0e-2)
    assert_small("Truncation fidelity delta", float(convergence["fidelity_trunc_delta"]), 3.0e-2)

    square_params = np.array(summary["representative"]["square"]["design"]["params"], dtype=float)
    square_design = get_family("square").builder(square_params, DEFAULT_CONFIG.representative_duration, DEFAULT_CONFIG.dt, DEFAULT_CONFIG.amp_max, DEFAULT_CONFIG.chi)
    zero_design = type(square_design)(
        family=square_design.family,
        params=square_design.params.copy(),
        waveform=np.zeros_like(square_design.waveform),
        dt=square_design.dt,
        duration=square_design.duration,
        delta_g=square_design.delta_g,
        metadata=dict(square_design.metadata),
    )
    zero_eval = evaluate_full_design(zero_design, DEFAULT_CONFIG)
    assert_small("Zero-drive residual photons", float(zero_eval.metrics.residual_photons), 1.0e-6)
    total_to_second = zero_design.duration + float(zero_eval.metadata["wait_time"])
    expected_defect = 0.5 * (1.0 - np.exp(-total_to_second / DEFAULT_CONFIG.t1))
    assert_small(
        "Zero-drive QND defect vs T1 floor",
        float(abs((1.0 - zero_eval.metrics.qnd_preservation) - expected_defect)),
        1.0e-4,
    )

    print("All validation checks passed.")


if __name__ == "__main__":
    main()
