"""Simulation helpers for linear and full-model readout replay."""

from __future__ import annotations

from dataclasses import dataclass, replace
import sys
from typing import Any

import numpy as np

from runtime_compat import patch_windows_qutip_import

from .bounds import solve_linear_response
from .config import CQED_SIM_PATH, ReadoutStudyConfig
from .metrics import ReadoutMetrics, build_metrics, cavity_expectation_from_quadratures, repeated_consistency_from_signals
from .pulse_families import PulseDesign

patch_windows_qutip_import()
if str(CQED_SIM_PATH) not in sys.path:
    sys.path.insert(0, str(CQED_SIM_PATH))

from cqed_sim.core import DispersiveTransmonCavityModel, FrameSpec  # noqa: E402
from cqed_sim.pulses import Pulse  # noqa: E402
from cqed_sim.sequence import SequenceCompiler  # noqa: E402
from cqed_sim.sim import SimulationConfig, prepare_simulation  # noqa: E402
from cqed_sim.sim.noise import NoiseSpec  # noqa: E402


@dataclass(frozen=True)
class ReplayEvaluation:
    """Evaluation bundle returned by either the linear or full-model replay path."""

    regime: str
    design: PulseDesign
    t_pulse: np.ndarray
    first_signal_g: np.ndarray
    first_signal_e: np.ndarray
    second_signal_g: np.ndarray
    second_signal_e: np.ndarray
    alpha_g: np.ndarray
    alpha_e: np.ndarray
    metrics: ReadoutMetrics
    metadata: dict[str, Any]


def build_model(cfg: ReadoutStudyConfig) -> DispersiveTransmonCavityModel:
    return DispersiveTransmonCavityModel(
        omega_c=cfg.omega_c,
        omega_q=cfg.omega_q,
        alpha=cfg.alpha,
        chi=cfg.chi,
        kerr=cfg.kerr,
        n_cav=cfg.n_cav,
        n_tr=cfg.n_tr,
    )


def make_noise(cfg: ReadoutStudyConfig) -> NoiseSpec:
    return NoiseSpec(t1=cfg.t1, tphi=cfg.tphi, kappa=cfg.kappa)


def pulse_tlist(design: PulseDesign) -> np.ndarray:
    return np.arange(len(design.waveform) + 1, dtype=float) * design.dt


def two_pulse_waveform(design: PulseDesign, wait_time: float) -> tuple[list[Pulse], float]:
    """Build the two-pulse protocol used for repeated-readout diagnostics."""
    sample_rate = 1.0 / design.dt
    pulse_1 = Pulse("cavity", 0.0, design.duration, design.waveform, amp=1.0, sample_rate=sample_rate, label="readout_1")
    pulse_2 = Pulse(
        "cavity",
        design.duration + wait_time,
        design.duration,
        design.waveform,
        amp=1.0,
        sample_rate=sample_rate,
        label="readout_2",
    )
    return [pulse_1, pulse_2], float(2.0 * design.duration + wait_time)


def build_frame(cfg: ReadoutStudyConfig, delta_g: float) -> FrameSpec:
    return FrameSpec(
        omega_q_frame=cfg.omega_q,
        omega_c_frame=cfg.drive_frequency(delta_g),
    )


def apply_overrides(cfg: ReadoutStudyConfig, **overrides: float | int) -> ReadoutStudyConfig:
    return replace(cfg, **overrides)


def evaluate_linear_design(
    design: PulseDesign,
    cfg: ReadoutStudyConfig,
    *,
    eta: float | None = None,
    wait_time: float | None = None,
) -> ReplayEvaluation:
    """Evaluate a design in the ideal linear-dispersive model."""
    eta = cfg.eta_nominal if eta is None else float(eta)
    wait_time = cfg.wait_time(design.duration) if wait_time is None else float(wait_time)
    n_wait = max(1, int(round(wait_time / design.dt)))
    wait_time = n_wait * design.dt
    t_first = pulse_tlist(design)
    response = solve_linear_response(design.waveform, t_first, kappa=cfg.kappa, chi=cfg.chi, delta_g=design.delta_g)

    eps_protocol = np.concatenate(
        [
            design.waveform,
            np.zeros(n_wait, dtype=np.complex128),
            design.waveform,
        ]
    )
    t_protocol = np.arange(len(eps_protocol) + 1, dtype=float) * design.dt
    repeated = solve_linear_response(eps_protocol, t_protocol, kappa=cfg.kappa, chi=cfg.chi, delta_g=design.delta_g)
    idx_start_2 = len(design.waveform) + n_wait
    idx_end_1 = len(design.waveform)
    first_signal_g = repeated.output_g[: idx_end_1 + 1]
    first_signal_e = repeated.output_e[: idx_end_1 + 1]
    second_signal_g = repeated.output_g[idx_start_2 : idx_start_2 + len(design.waveform) + 1]
    second_signal_e = repeated.output_e[idx_start_2 : idx_start_2 + len(design.waveform) + 1]
    t_axis = np.arange(len(first_signal_g), dtype=float) * design.dt

    metrics = build_metrics(
        output_g=first_signal_g,
        output_e=first_signal_e,
        tlist=t_axis,
        eta=eta,
        t1=cfg.t1,
        alpha_g=response.alpha_g,
        alpha_e=response.alpha_e,
        residual_after_wait_value=0.5 * (abs(repeated.alpha_g[idx_start_2]) ** 2 + abs(repeated.alpha_e[idx_start_2]) ** 2),
        qnd_preservation=1.0,
        leakage=0.0,
        repeat_consistency=repeated_consistency_from_signals(
            first_g=first_signal_g,
            first_e=first_signal_e,
            second_g=second_signal_g,
            second_e=second_signal_e,
            tlist=t_axis,
            eta=eta,
        ),
        n_crit=cfg.n_crit,
    )
    return ReplayEvaluation(
        regime="linear",
        design=design,
        t_pulse=t_axis,
        first_signal_g=first_signal_g,
        first_signal_e=first_signal_e,
        second_signal_g=second_signal_g,
        second_signal_e=second_signal_e,
        alpha_g=response.alpha_g,
        alpha_e=response.alpha_e,
        metrics=metrics,
        metadata={"wait_time": wait_time, "n_wait_steps": n_wait},
    )


def evaluate_full_design(
    design: PulseDesign,
    cfg: ReadoutStudyConfig,
    *,
    eta: float | None = None,
    wait_time: float | None = None,
) -> ReplayEvaluation:
    """Replay a design through `cqed_sim` with cavity damping and qubit decoherence."""
    eta = cfg.eta_nominal if eta is None else float(eta)
    wait_time = cfg.wait_time(design.duration) if wait_time is None else float(wait_time)
    n_wait = max(1, int(round(wait_time / design.dt)))
    wait_time = n_wait * design.dt

    model = build_model(cfg)
    pulses, t_end = two_pulse_waveform(design, wait_time)
    compiler = SequenceCompiler(dt=design.dt)
    compiled = compiler.compile(pulses, t_end=t_end)
    session = prepare_simulation(
        model,
        compiled,
        {"cavity": "cavity"},
        config=SimulationConfig(frame=build_frame(cfg, design.delta_g)),
        noise=make_noise(cfg),
    )
    state_g = model.basis_state(0, 0)
    state_e = model.basis_state(1, 0)
    result_g, result_e = session.run_many([state_g, state_e], max_workers=1)

    n_steps = len(design.waveform)
    idx_end_1 = n_steps
    idx_start_2 = n_steps + n_wait

    alpha_trace_g = cavity_expectation_from_quadratures(result_g.expectations["x_c"], result_g.expectations["p_c"])
    alpha_trace_e = cavity_expectation_from_quadratures(result_e.expectations["x_c"], result_e.expectations["p_c"])
    root_kappa = float(np.sqrt(cfg.kappa))
    output_trace_g = root_kappa * alpha_trace_g
    output_trace_e = root_kappa * alpha_trace_e

    first_signal_g = output_trace_g[: idx_end_1 + 1]
    first_signal_e = output_trace_e[: idx_end_1 + 1]
    second_signal_g = output_trace_g[idx_start_2 : idx_start_2 + n_steps + 1]
    second_signal_e = output_trace_e[idx_start_2 : idx_start_2 + n_steps + 1]
    alpha_first_g = alpha_trace_g[: idx_end_1 + 1]
    alpha_first_e = alpha_trace_e[: idx_end_1 + 1]
    t_axis = np.arange(len(first_signal_g), dtype=float) * design.dt

    p_g_preserve = float(result_g.expectations["P_g"][idx_start_2])
    p_e_preserve = float(result_e.expectations["P_e"][idx_start_2])
    leakage_g = max(0.0, 1.0 - float(result_g.expectations["P_g"][idx_start_2]) - float(result_g.expectations["P_e"][idx_start_2]))
    leakage_e = max(0.0, 1.0 - float(result_e.expectations["P_g"][idx_start_2]) - float(result_e.expectations["P_e"][idx_start_2]))
    qnd_preservation = 0.5 * (p_g_preserve + p_e_preserve)
    leakage = 0.5 * (leakage_g + leakage_e)

    metrics = build_metrics(
        output_g=first_signal_g,
        output_e=first_signal_e,
        tlist=t_axis,
        eta=eta,
        t1=cfg.t1,
        alpha_g=alpha_first_g,
        alpha_e=alpha_first_e,
        residual_after_wait_value=0.5 * (
            abs(alpha_trace_g[idx_start_2]) ** 2 + abs(alpha_trace_e[idx_start_2]) ** 2
        ),
        qnd_preservation=qnd_preservation,
        leakage=leakage,
        repeat_consistency=repeated_consistency_from_signals(
            first_g=first_signal_g,
            first_e=first_signal_e,
            second_g=second_signal_g,
            second_e=second_signal_e,
            tlist=t_axis,
            eta=eta,
        ),
        n_crit=cfg.n_crit,
    )
    return ReplayEvaluation(
        regime="full",
        design=design,
        t_pulse=t_axis,
        first_signal_g=first_signal_g,
        first_signal_e=first_signal_e,
        second_signal_g=second_signal_g,
        second_signal_e=second_signal_e,
        alpha_g=alpha_first_g,
        alpha_e=alpha_first_e,
        metrics=metrics,
        metadata={
            "wait_time": wait_time,
            "n_wait_steps": n_wait,
            "qnd_g": p_g_preserve,
            "qnd_e": p_e_preserve,
            "leakage_g": leakage_g,
            "leakage_e": leakage_e,
        },
    )
