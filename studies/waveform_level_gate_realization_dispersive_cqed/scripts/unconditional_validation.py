"""Spot-check validation for the unconditional displacement study."""

from __future__ import annotations

import time

import common
from common import ARTIFACTS_DIR, build_frame, build_model, compile_and_prepare, save_json
from unconditional_displacement_study import (
    branch_vacuum_metrics,
    single_tone_pulse,
    variant_model,
)


def evaluate_square_case(*, model, frame, alpha: complex, duration_s: float, dt: float) -> dict:
    pulses, _meta = single_tone_pulse(
        model=model,
        frame=frame,
        alpha=alpha,
        duration_s=duration_s,
        family="square",
        filter_bw_mhz=None,
    )
    session = compile_and_prepare(model, frame, pulses, dt=dt)
    return branch_vacuum_metrics(model, session, alpha_target=alpha)


def main() -> None:
    started = time.time()
    alpha = 1.0
    duration_s = 80.0e-9

    ideal_model = build_model(chi=0.0, chi_prime=None, kerr=0.0, n_cav=15, n_tr=3)
    ideal_frame = build_frame(ideal_model)
    ideal_metrics = evaluate_square_case(
        model=ideal_model,
        frame=ideal_frame,
        alpha=alpha,
        duration_s=duration_s,
        dt=common.DEFAULT_DT,
    )

    ncav_rows = []
    for n_cav in (12, 15, 18):
        model = variant_model("full", n_cav=n_cav, n_tr=3)
        frame = build_frame(model)
        metrics = evaluate_square_case(
            model=model,
            frame=frame,
            alpha=alpha,
            duration_s=duration_s,
            dt=common.DEFAULT_DT,
        )
        ncav_rows.append({"n_cav": n_cav, **metrics})

    dt_rows = []
    for dt_ns in (0.25, 0.5, 1.0):
        model = variant_model("full", n_cav=15, n_tr=3)
        frame = build_frame(model)
        metrics = evaluate_square_case(
            model=model,
            frame=frame,
            alpha=alpha,
            duration_s=duration_s,
            dt=dt_ns * 1.0e-9,
        )
        dt_rows.append({"dt_ns": dt_ns, **metrics})

    ntr_rows = []
    for n_tr in (3, 4):
        model = variant_model("full", n_cav=15, n_tr=n_tr)
        frame = build_frame(model)
        metrics = evaluate_square_case(
            model=model,
            frame=frame,
            alpha=alpha,
            duration_s=duration_s,
            dt=common.DEFAULT_DT,
        )
        ntr_rows.append({"n_tr": n_tr, **metrics})

    payload = {
        "representative_case": {"protocol": "naive square", "alpha_target": alpha, "duration_ns": 80.0},
        "ideal_limit": ideal_metrics,
        "n_cav_sweep": ncav_rows,
        "dt_sweep": dt_rows,
        "n_tr_sweep": ntr_rows,
        "wall_time_s": time.time() - started,
    }
    save_json(ARTIFACTS_DIR / "unconditional_validation_spotcheck.json", payload)


if __name__ == "__main__":
    main()
