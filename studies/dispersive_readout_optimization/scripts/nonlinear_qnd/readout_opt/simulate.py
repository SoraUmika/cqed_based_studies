"""Simulation helpers for linear, multilevel, hardware-aware, and nonlinear readout replay."""

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
from cqed_sim.core.drive_targets import TransmonTransitionDriveSpec  # noqa: E402
from cqed_sim.core.frequencies import carrier_for_transition_frequency, transmon_transition_frequency  # noqa: E402
from cqed_sim.pulses import Pulse  # noqa: E402
from cqed_sim.pulses.hardware import HardwareConfig  # noqa: E402
from cqed_sim.sequence import SequenceCompiler  # noqa: E402
from cqed_sim.sim import SimulationConfig, prepare_simulation  # noqa: E402
from cqed_sim.sim.noise import NoiseSpec  # noqa: E402


@dataclass(frozen=True)
class ReplayEvaluation:
    """Evaluation bundle returned by either the linear or physical replay path."""

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


def make_noise(cfg: ReadoutStudyConfig, *, multilevel: bool) -> NoiseSpec:
    if multilevel:
        return NoiseSpec(
            transmon_t1=cfg.transmon_t1,
            tphi=cfg.tphi,
            kappa=cfg.kappa,
            nth=cfg.nth,
        )
    return NoiseSpec(t1=cfg.t1, tphi=cfg.tphi, kappa=cfg.kappa)


def pulse_tlist(design: PulseDesign) -> np.ndarray:
    return np.arange(len(design.waveform) + 1, dtype=float) * design.dt


def build_frame(cfg: ReadoutStudyConfig, delta_g: float) -> FrameSpec:
    return FrameSpec(
        omega_q_frame=cfg.omega_q,
        omega_c_frame=cfg.drive_frequency(delta_g),
    )


def apply_overrides(cfg: ReadoutStudyConfig, **overrides: float | int | object) -> ReadoutStudyConfig:
    return replace(cfg, **overrides)


def _two_pulse_wait_steps(design: PulseDesign, wait_time: float | None, cfg: ReadoutStudyConfig) -> tuple[int, float]:
    wait = cfg.wait_time(design.duration) if wait_time is None else float(wait_time)
    n_wait = max(1, int(round(wait / design.dt)))
    return n_wait, float(n_wait * design.dt)


def _phase_slew_limit(waveform: np.ndarray, dt: float, max_slew: float | None) -> np.ndarray:
    if max_slew is None or max_slew <= 0.0:
        return np.asarray(waveform, dtype=np.complex128).copy()
    wf = np.asarray(waveform, dtype=np.complex128)
    amp = np.abs(wf)
    raw_phase = np.angle(wf)
    phase = np.empty_like(raw_phase)
    phase[0] = raw_phase[0]
    for idx in range(1, len(raw_phase)):
        candidate = raw_phase[idx]
        if amp[idx] <= 1.0e-15:
            candidate = phase[idx - 1]
        delta = np.angle(np.exp(1j * (candidate - phase[idx - 1])))
        max_step = float(max_slew) * float(dt)
        phase[idx] = phase[idx - 1] + np.clip(delta, -max_step, max_step)
    return amp * np.exp(1j * phase)


def _transport_waveform(design: PulseDesign, cfg: ReadoutStudyConfig, *, use_hardware: bool) -> np.ndarray:
    waveform = np.asarray(design.waveform, dtype=np.complex128)
    if not use_hardware:
        return waveform.copy()
    return _phase_slew_limit(waveform, design.dt, cfg.hardware.phase_slew_rad_per_s)


def _mean_linear_occupancy(waveform: np.ndarray, design: PulseDesign, cfg: ReadoutStudyConfig) -> tuple[np.ndarray, float]:
    t_axis = np.arange(len(waveform) + 1, dtype=float) * design.dt
    response = solve_linear_response(waveform, t_axis, kappa=cfg.kappa, chi=cfg.chi, delta_g=design.delta_g)
    n_g = 0.5 * (np.abs(response.alpha_g[:-1]) ** 2 + np.abs(response.alpha_g[1:]) ** 2)
    n_e = 0.5 * (np.abs(response.alpha_e[:-1]) ** 2 + np.abs(response.alpha_e[1:]) ** 2)
    n_mean = 0.5 * (n_g + n_e)
    return n_mean, float(np.max(n_mean) if n_mean.size else 0.0)


def _mixing_waveforms(waveform: np.ndarray, design: PulseDesign, cfg: ReadoutStudyConfig) -> tuple[np.ndarray, np.ndarray, float]:
    n_mean, n_peak = _mean_linear_occupancy(waveform, design, cfg)
    if not np.isfinite(cfg.n_crit) or cfg.n_crit <= 0.0:
        activation = np.zeros_like(n_mean, dtype=float)
    else:
        ratio = n_mean / cfg.n_crit
        onset = cfg.mixing.onset_ratio
        activation = np.clip((ratio - onset) / max(onset, 1.0e-12), 0.0, 3.0)
        activation = activation ** cfg.mixing.occupancy_exponent
    slew = np.concatenate([np.zeros(1, dtype=np.complex128), np.diff(waveform)])
    phase_term = np.exp(1j * cfg.mixing.phase_lag_rad)
    slew_phase_term = np.exp(1j * cfg.mixing.slew_phase_rad)
    ge = (
        cfg.mixing.ge_scale * activation * waveform * phase_term
        + cfg.mixing.slew_ge_scale * activation * slew * slew_phase_term
    ).astype(np.complex128)
    ef = (
        cfg.mixing.ef_scale * activation * waveform * phase_term
        + cfg.mixing.slew_ef_scale * activation * slew * slew_phase_term
    ).astype(np.complex128)
    return ge, ef, n_peak


def _hardware_dict(cfg: ReadoutStudyConfig) -> dict[str, HardwareConfig]:
    hw = cfg.hardware
    return {
        "cavity": HardwareConfig(
            gain_i=hw.gain_i,
            gain_q=hw.gain_q,
            quadrature_skew=hw.quadrature_skew_rad,
            dc_i=hw.dc_i,
            dc_q=hw.dc_q,
            image_leakage=hw.image_leakage,
            channel_gain=hw.channel_gain,
            zoh_samples=hw.zoh_samples,
            lowpass_bw=hw.lowpass_bw_hz,
            detuning=hw.detuning,
            timing_quantum=hw.timing_quantum,
            amplitude_bits=hw.amplitude_bits,
        )
    }


def _build_protocol(
    design: PulseDesign,
    cfg: ReadoutStudyConfig,
    *,
    use_hardware: bool,
    use_effective_mixing: bool,
    wait_time: float | None,
) -> dict[str, Any]:
    n_wait, wait = _two_pulse_wait_steps(design, wait_time, cfg)
    sample_rate = 1.0 / design.dt
    transport_waveform = _transport_waveform(design, cfg, use_hardware=use_hardware)
    model = build_model(cfg)
    frame = build_frame(cfg, design.delta_g)

    pulses: list[Pulse] = [
        Pulse("cavity", 0.0, design.duration, transport_waveform, amp=1.0, sample_rate=sample_rate, label="readout_1"),
        Pulse(
            "cavity",
            design.duration + wait,
            design.duration,
            transport_waveform,
            amp=1.0,
            sample_rate=sample_rate,
            label="readout_2",
        ),
    ]
    drive_ops: dict[str, str | TransmonTransitionDriveSpec] = {"cavity": "cavity"}
    mixing_peak_est = 0.0
    mix_ge_peak = 0.0
    mix_ef_peak = 0.0

    if use_effective_mixing:
        mix_ge, mix_ef, mixing_peak_est = _mixing_waveforms(transport_waveform, design, cfg)
        ge_carrier = carrier_for_transition_frequency(
            transmon_transition_frequency(model, cavity_level=0, lower_level=0, upper_level=1, frame=frame)
        )
        pulses.extend(
            [
                Pulse(
                    "mix_ge",
                    0.0,
                    design.duration,
                    mix_ge,
                    carrier=ge_carrier,
                    amp=1.0,
                    sample_rate=sample_rate,
                    label="mix_ge_1",
                ),
                Pulse(
                    "mix_ge",
                    design.duration + wait,
                    design.duration,
                    mix_ge,
                    carrier=ge_carrier,
                    amp=1.0,
                    sample_rate=sample_rate,
                    label="mix_ge_2",
                ),
            ]
        )
        drive_ops["mix_ge"] = TransmonTransitionDriveSpec(lower_level=0, upper_level=1)
        mix_ge_peak = float(np.max(np.abs(mix_ge)) if mix_ge.size else 0.0)
        if cfg.n_tr >= 3:
            ef_carrier = carrier_for_transition_frequency(
                transmon_transition_frequency(model, cavity_level=0, lower_level=1, upper_level=2, frame=frame)
            )
            pulses.extend(
                [
                    Pulse(
                        "mix_ef",
                        0.0,
                        design.duration,
                        mix_ef,
                        carrier=ef_carrier,
                        amp=1.0,
                        sample_rate=sample_rate,
                        label="mix_ef_1",
                    ),
                    Pulse(
                        "mix_ef",
                        design.duration + wait,
                        design.duration,
                        mix_ef,
                        carrier=ef_carrier,
                        amp=1.0,
                        sample_rate=sample_rate,
                        label="mix_ef_2",
                    ),
                ]
            )
            drive_ops["mix_ef"] = TransmonTransitionDriveSpec(lower_level=1, upper_level=2)
            mix_ef_peak = float(np.max(np.abs(mix_ef)) if mix_ef.size else 0.0)

    compiler = SequenceCompiler(dt=design.dt, hardware=_hardware_dict(cfg) if use_hardware else None)
    t_end = float(2.0 * design.duration + wait)
    compiled = compiler.compile(pulses, t_end=t_end)
    n_steps = len(transport_waveform)
    distorted_first = np.asarray(compiled.channels["cavity"].distorted[:n_steps], dtype=np.complex128)
    baseband_first = np.asarray(compiled.channels["cavity"].baseband[:n_steps], dtype=np.complex128)
    diff_norm = np.linalg.norm(distorted_first - np.asarray(design.waveform, dtype=np.complex128))
    ref_norm = max(np.linalg.norm(np.asarray(design.waveform, dtype=np.complex128)), 1.0e-15)

    return {
        "model": model,
        "frame": frame,
        "compiled": compiled,
        "drive_ops": drive_ops,
        "n_wait": n_wait,
        "wait_time": wait,
        "transport_waveform": transport_waveform,
        "baseband_first": baseband_first,
        "distorted_first": distorted_first,
        "hardware_rms_error": float(diff_norm / ref_norm),
        "mixing_peak_est": float(mixing_peak_est),
        "mix_ge_peak": float(mix_ge_peak),
        "mix_ef_peak": float(mix_ef_peak),
    }


def transport_analysis(
    design: PulseDesign,
    cfg: ReadoutStudyConfig,
    *,
    regime: str,
) -> dict[str, np.ndarray | float]:
    use_hardware = regime in {"hardware", "rich"}
    use_effective_mixing = regime in {"nonlinear", "rich"}
    protocol = _build_protocol(
        design,
        cfg,
        use_hardware=use_hardware,
        use_effective_mixing=use_effective_mixing,
        wait_time=None,
    )
    return {
        "program_waveform": np.asarray(design.waveform, dtype=np.complex128),
        "transport_waveform": np.asarray(protocol["transport_waveform"], dtype=np.complex128),
        "distorted_waveform": np.asarray(protocol["distorted_first"], dtype=np.complex128),
        "hardware_rms_error": float(protocol["hardware_rms_error"]),
        "mixing_peak_est": float(protocol["mixing_peak_est"]),
        "mix_ge_peak": float(protocol["mix_ge_peak"]),
        "mix_ef_peak": float(protocol["mix_ef_peak"]),
    }


def evaluate_linear_design(
    design: PulseDesign,
    cfg: ReadoutStudyConfig,
    *,
    eta: float | None = None,
    wait_time: float | None = None,
) -> ReplayEvaluation:
    """Evaluate a design in the ideal linear-dispersive model."""
    eta = cfg.eta_nominal if eta is None else float(eta)
    n_wait, wait = _two_pulse_wait_steps(design, wait_time, cfg)
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
        baseline_transition=0.0,
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
        metadata={"wait_time": wait, "n_wait_steps": n_wait, "hardware_rms_error": 0.0},
    )


def _evaluate_physical_design(
    design: PulseDesign,
    cfg: ReadoutStudyConfig,
    *,
    regime: str,
    eta: float | None,
    wait_time: float | None,
    use_multilevel: bool,
    use_hardware: bool,
    use_effective_mixing: bool,
) -> ReplayEvaluation:
    eta = cfg.eta_nominal if eta is None else float(eta)
    protocol = _build_protocol(
        design,
        cfg,
        use_hardware=use_hardware,
        use_effective_mixing=use_effective_mixing,
        wait_time=wait_time,
    )
    model = protocol["model"]
    wait = float(protocol["wait_time"])
    n_wait = int(protocol["n_wait"])
    session = prepare_simulation(
        model,
        protocol["compiled"],
        protocol["drive_ops"],
        config=SimulationConfig(frame=protocol["frame"]),
        noise=make_noise(cfg, multilevel=use_multilevel),
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
    baseline_transition = 0.5 * (1.0 - np.exp(-(design.duration + wait) / cfg.t1))

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
        baseline_transition=baseline_transition,
    )
    return ReplayEvaluation(
        regime=regime,
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
            "wait_time": wait,
            "n_wait_steps": n_wait,
            "qnd_g": p_g_preserve,
            "qnd_e": p_e_preserve,
            "leakage_g": leakage_g,
            "leakage_e": leakage_e,
            "hardware_rms_error": float(protocol["hardware_rms_error"]),
            "mixing_peak_est": float(protocol["mixing_peak_est"]),
            "mix_ge_peak": float(protocol["mix_ge_peak"]),
            "mix_ef_peak": float(protocol["mix_ef_peak"]),
            "baseline_transition": float(baseline_transition),
        },
    )


def evaluate_full_design(
    design: PulseDesign,
    cfg: ReadoutStudyConfig,
    *,
    eta: float | None = None,
    wait_time: float | None = None,
) -> ReplayEvaluation:
    """Legacy two-mode replay with aggregate T1/Tphi and no transport or extra mixing."""
    return _evaluate_physical_design(
        design,
        cfg,
        regime="full",
        eta=eta,
        wait_time=wait_time,
        use_multilevel=False,
        use_hardware=False,
        use_effective_mixing=False,
    )


def evaluate_multilevel_design(
    design: PulseDesign,
    cfg: ReadoutStudyConfig,
    *,
    eta: float | None = None,
    wait_time: float | None = None,
) -> ReplayEvaluation:
    """Replay with multilevel transmon relaxation and thermal occupation, but no hardware distortion."""
    return _evaluate_physical_design(
        design,
        cfg,
        regime="multilevel",
        eta=eta,
        wait_time=wait_time,
        use_multilevel=True,
        use_hardware=False,
        use_effective_mixing=False,
    )


def evaluate_hardware_design(
    design: PulseDesign,
    cfg: ReadoutStudyConfig,
    *,
    eta: float | None = None,
    wait_time: float | None = None,
) -> ReplayEvaluation:
    """Replay with multilevel noise plus native hardware transport distortion."""
    return _evaluate_physical_design(
        design,
        cfg,
        regime="hardware",
        eta=eta,
        wait_time=wait_time,
        use_multilevel=True,
        use_hardware=True,
        use_effective_mixing=False,
    )


def evaluate_nonlinear_design(
    design: PulseDesign,
    cfg: ReadoutStudyConfig,
    *,
    eta: float | None = None,
    wait_time: float | None = None,
) -> ReplayEvaluation:
    """Replay with multilevel noise plus occupancy- and slew-driven effective mixing."""
    return _evaluate_physical_design(
        design,
        cfg,
        regime="nonlinear",
        eta=eta,
        wait_time=wait_time,
        use_multilevel=True,
        use_hardware=False,
        use_effective_mixing=True,
    )


def evaluate_rich_design(
    design: PulseDesign,
    cfg: ReadoutStudyConfig,
    *,
    eta: float | None = None,
    wait_time: float | None = None,
) -> ReplayEvaluation:
    """Replay with multilevel noise, native hardware distortion, and effective strong-drive mixing."""
    return _evaluate_physical_design(
        design,
        cfg,
        regime="rich",
        eta=eta,
        wait_time=wait_time,
        use_multilevel=True,
        use_hardware=True,
        use_effective_mixing=True,
    )
