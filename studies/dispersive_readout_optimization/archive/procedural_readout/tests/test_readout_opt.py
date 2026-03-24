from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from readout_opt import DEFAULT_CONFIG, assignment_fidelity_from_snr2, classification_probability, set_nulling_tail_kappa, solve_linear_response
from readout_opt.pulse_families import get_family


def test_classification_probability_matches_template_limit() -> None:
    tlist = np.linspace(0.0, 1.0, 5)
    signal_g = np.zeros_like(tlist, dtype=np.complex128)
    signal_e = np.ones_like(tlist, dtype=np.complex128)
    prob_e = classification_probability(signal_e, signal_g, signal_e, tlist, eta=1.0, target="e")
    expected = assignment_fidelity_from_snr2(float(np.trapezoid(np.ones_like(tlist), tlist)))
    assert abs(prob_e - expected) < 1.0e-12


def test_nulling_tail_exactly_zeros_linear_response() -> None:
    cfg = DEFAULT_CONFIG
    set_nulling_tail_kappa(cfg.kappa)
    params = np.array([0.20, 0.0, 0.15, 0.0, 1.0, 2.0, 1.0, 1.0, 0.0], dtype=float)
    design = get_family("nulling_tail").builder(params, cfg.representative_duration, cfg.dt, cfg.amp_max, cfg.chi)
    tlist = np.arange(len(design.waveform) + 1, dtype=float) * cfg.dt
    response = solve_linear_response(design.waveform, tlist, kappa=cfg.kappa, chi=cfg.chi, delta_g=design.delta_g)
    assert abs(response.alpha_g[-1]) < 1.0e-9
    assert abs(response.alpha_e[-1]) < 1.0e-9


def test_square_family_respects_amplitude_cap() -> None:
    cfg = DEFAULT_CONFIG
    design = get_family("square").builder(np.array([1.0, 0.0, 0.0]), cfg.representative_duration, cfg.dt, cfg.amp_max, cfg.chi)
    assert np.max(np.abs(design.waveform)) <= cfg.amp_max + 1.0e-12
