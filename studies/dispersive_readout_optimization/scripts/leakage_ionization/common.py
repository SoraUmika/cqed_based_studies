"""Shared helpers for the measurement-induced leakage and ionization study."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from math import erf, sqrt
from pathlib import Path
import json
import sys

import numpy as np

try:
    from .runtime_compat import patch_windows_qutip_import
except ImportError:  # pragma: no cover - direct script execution path
    from runtime_compat import patch_windows_qutip_import

patch_windows_qutip_import()

SCRIPT_DIR = Path(__file__).resolve().parent
STUDY_DIR = SCRIPT_DIR.parent
DATA_DIR = STUDY_DIR / "data"
FIGURES_DIR = STUDY_DIR / "figures"
REPORT_DIR = STUDY_DIR / "report"
CQED_SIM_PATH = Path(
    r"C:\Users\dazzl\Box\Shyam Shankar Quantum Circuits Group\Users\Users_JianJun\cQED_simulation"
)

if str(CQED_SIM_PATH) not in sys.path:
    sys.path.insert(0, str(CQED_SIM_PATH))

from cqed_sim.core import (  # noqa: E402
    DispersiveTransmonCavityModel,
    FrameSpec,
    carrier_for_transition_frequency,
    transmon_transition_frequency,
)
from cqed_sim.measurement import (  # noqa: E402
    AmplifierChain,
    ContinuousReadoutSpec,
    ReadoutChain,
    ReadoutResonator,
    StrongReadoutMixingSpec,
    build_strong_readout_disturbance,
    integrate_measurement_record,
    simulate_continuous_readout,
    strong_readout_drive_targets,
)
from cqed_sim.pulses import HardwareConfig, Pulse  # noqa: E402
from cqed_sim.sequence import SequenceCompiler  # noqa: E402
from cqed_sim.sim import NoiseSpec, SimulationConfig, prepare_simulation  # noqa: E402

TWO_PI = 2.0 * np.pi


@dataclass(frozen=True)
class StudyConfig:
    """Nominal physical, numerical, and sweep settings for the study."""

    omega_q: float = TWO_PI * 6.150e9
    omega_c: float = TWO_PI * 8.597e9
    alpha: float = TWO_PI * (-255.0e6)
    chi: float = TWO_PI * (-2.84e6)
    kappa: float = TWO_PI * 2.4e6
    kerr: float = TWO_PI * (-28.0e3)
    n_tr: int = 5
    n_cav: int = 18
    t1_ge: float = 30.0e-6
    t1_ef: float = 12.0e-6
    t1_high: float = 8.0e-6
    t2: float = 20.0e-6
    nth: float = 0.01
    dt: float = 4.0e-9
    post_buffer: float = 80.0e-9
    ion_level: int = 3
    amplitude_values_mhz: tuple[float, ...] = (1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5)
    duration_values_ns: tuple[float, ...] = (120.0, 240.0, 360.0, 480.0)
    detuning_values_mhz: tuple[float, ...] = (-1.0, -0.5, 0.0, 0.5, 1.0)
    bandwidth_values_mhz: tuple[float, ...] = (35.0, 75.0, 150.0)
    mitigation_shapes: tuple[str, ...] = ("square", "cosine", "gaussian")
    baseline_shape: str = "cosine"
    baseline_bandwidth_mhz: float = 75.0
    target_assignment: float = 0.95
    repeated_wait_factor: float = 2.5
    readout_noise_temperature: float = 4.0
    readout_gain: float = 12.0
    ntraj_representative: int = 48
    mixing_onset_ratio: float = 0.10
    mixing_occupancy_exponent: float = 1.20
    mixing_ge_scale: float = 6.0
    mixing_ef_scale: float = 3.6
    mixing_slew_ge_scale: float = 3.0
    mixing_slew_ef_scale: float = 1.5
    mixing_phase_lag_deg: float = -60.0
    mixing_slew_phase_deg: float = 35.0
    higher_ladder_scales: tuple[float, ...] = (0.45, 0.20, 0.10)
    seed: int = 7

    @property
    def tphi(self) -> float:
        inv_tphi = max(0.0, 1.0 / self.t2 - 1.0 / (2.0 * self.t1_ge))
        return float(np.inf if inv_tphi <= 0.0 else 1.0 / inv_tphi)

    @property
    def delta_qc(self) -> float:
        return float(self.omega_q - self.omega_c)

    @property
    def g_est(self) -> float:
        numerator = abs(self.chi) * abs(self.delta_qc) * abs(self.delta_qc + self.alpha)
        return float(np.sqrt(max(numerator / abs(self.alpha), 0.0)))

    @property
    def n_crit(self) -> float:
        return float((self.delta_qc / (2.0 * self.g_est)) ** 2)

    @property
    def amplitude_values(self) -> np.ndarray:
        return TWO_PI * 1.0e6 * np.asarray(self.amplitude_values_mhz, dtype=float)

    @property
    def duration_values(self) -> np.ndarray:
        return 1.0e-9 * np.asarray(self.duration_values_ns, dtype=float)

    @property
    def detuning_values(self) -> np.ndarray:
        return TWO_PI * 1.0e6 * np.asarray(self.detuning_values_mhz, dtype=float)

    @property
    def bandwidth_values(self) -> np.ndarray:
        return 1.0e6 * np.asarray(self.bandwidth_values_mhz, dtype=float)

    def wait_time(self, duration: float) -> float:
        return float(max(0.4 * duration, self.repeated_wait_factor / self.kappa))

    def mixing_spec(self) -> StrongReadoutMixingSpec:
        return StrongReadoutMixingSpec(
            n_crit=self.n_crit,
            onset_ratio=self.mixing_onset_ratio,
            occupancy_exponent=self.mixing_occupancy_exponent,
            ge_scale=self.mixing_ge_scale,
            ef_scale=self.mixing_ef_scale,
            slew_ge_scale=self.mixing_slew_ge_scale,
            slew_ef_scale=self.mixing_slew_ef_scale,
            phase_lag=np.deg2rad(self.mixing_phase_lag_deg),
            slew_phase=np.deg2rad(self.mixing_slew_phase_deg),
            higher_ladder_scales=self.higher_ladder_scales,
        )

    def as_dict(self) -> dict[str, object]:
        payload = asdict(self)
        payload["tphi"] = self.tphi
        payload["g_est"] = self.g_est
        payload["n_crit"] = self.n_crit
        return payload


@dataclass
class ProtocolBundle:
    model: DispersiveTransmonCavityModel
    frame: FrameSpec
    compiled: object
    drive_ops: dict[str, object]
    noise: NoiseSpec
    readout_chain: ReadoutChain
    distorted_waveform: np.ndarray
    disturbance: object
    duration: float
    detuning: float
    lowpass_bw_hz: float
    hardware_rms_error: float
    shape: str
    amplitude: float


def ensure_directories() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)


def save_json(path: Path, payload: object) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)


def _transmon_t1_ladder(cfg: StudyConfig, n_tr: int) -> tuple[float, ...]:
    ladder: list[float] = []
    for upper in range(1, n_tr):
        if upper == 1:
            ladder.append(cfg.t1_ge)
        elif upper == 2:
            ladder.append(cfg.t1_ef)
        else:
            ladder.append(max(cfg.t1_high / (1.0 + 0.2 * (upper - 3)), 2.0e-6))
    return tuple(float(value) for value in ladder)


def make_model(cfg: StudyConfig, *, n_tr: int | None = None, n_cav: int | None = None) -> DispersiveTransmonCavityModel:
    return DispersiveTransmonCavityModel(
        omega_c=cfg.omega_c,
        omega_q=cfg.omega_q,
        alpha=cfg.alpha,
        chi=cfg.chi,
        kerr=cfg.kerr,
        n_cav=cfg.n_cav if n_cav is None else int(n_cav),
        n_tr=cfg.n_tr if n_tr is None else int(n_tr),
    )


def make_noise(cfg: StudyConfig, *, n_tr: int | None = None) -> NoiseSpec:
    levels = cfg.n_tr if n_tr is None else int(n_tr)
    return NoiseSpec(
        transmon_t1=_transmon_t1_ladder(cfg, levels),
        tphi=cfg.tphi,
        kappa=cfg.kappa,
        nth=cfg.nth,
    )


def hardware_config(lowpass_bw_hz: float | None) -> dict[str, HardwareConfig] | None:
    if lowpass_bw_hz is None:
        return None
    return {
        "cavity": HardwareConfig(
            gain_i=1.0,
            gain_q=0.97,
            quadrature_skew=np.deg2rad(4.0),
            image_leakage=0.02,
            channel_gain=0.98,
            zoh_samples=2,
            lowpass_bw=float(lowpass_bw_hz),
            timing_quantum=8.0e-9,
            amplitude_bits=10,
        )
    }


def _raised_cosine_window(n_steps: int, rise_steps: int) -> np.ndarray:
    if rise_steps <= 0 or 2 * rise_steps >= n_steps:
        return np.ones(n_steps, dtype=float)
    ramp = np.sin(np.linspace(0.0, 0.5 * np.pi, rise_steps, endpoint=False)) ** 2
    flat = np.ones(n_steps - 2 * rise_steps, dtype=float)
    return np.concatenate([ramp, flat, ramp[::-1]])


def make_envelope(shape: str, amplitude: float, duration: float, dt: float) -> tuple[np.ndarray, float]:
    n_steps = max(1, int(np.round(float(duration) / float(dt))))
    actual_duration = float(n_steps * dt)
    rise_steps = max(2, int(round(0.18 * n_steps)))
    if shape == "square":
        window = np.ones(n_steps, dtype=float)
    elif shape == "cosine":
        window = _raised_cosine_window(n_steps, rise_steps)
    elif shape == "gaussian":
        axis = np.linspace(-1.0, 1.0, n_steps)
        window = np.exp(-0.5 * (axis / 0.35) ** 2)
        window /= np.max(window)
    else:
        raise ValueError(f"Unsupported shape '{shape}'.")
    return float(amplitude) * window.astype(np.complex128), actual_duration


def build_protocol(
    cfg: StudyConfig,
    *,
    amplitude: float,
    duration: float,
    detuning: float,
    shape: str,
    lowpass_bw_hz: float | None,
    n_tr: int | None = None,
    n_cav: int | None = None,
) -> ProtocolBundle:
    model = make_model(cfg, n_tr=n_tr, n_cav=n_cav)
    noise = make_noise(cfg, n_tr=model.n_tr)
    drive_frequency = float(cfg.omega_c + detuning)
    frame = FrameSpec(omega_q_frame=cfg.omega_q, omega_c_frame=drive_frequency)
    readout_waveform, actual_duration = make_envelope(shape, amplitude, duration, cfg.dt)
    compiler = SequenceCompiler(dt=cfg.dt, hardware=hardware_config(lowpass_bw_hz))
    base_pulse = Pulse(
        "cavity",
        0.0,
        actual_duration,
        readout_waveform,
        amp=1.0,
        sample_rate=1.0 / cfg.dt,
        label="readout",
    )
    cavity_only = compiler.compile([base_pulse], t_end=actual_duration + cfg.post_buffer)
    n_samples = readout_waveform.size
    distorted_waveform = np.asarray(cavity_only.channels["cavity"].distorted[:n_samples], dtype=np.complex128)
    reference_norm = max(np.linalg.norm(readout_waveform), 1.0e-15)
    hardware_rms_error = float(np.linalg.norm(distorted_waveform - readout_waveform) / reference_norm)

    resonator = ReadoutResonator(
        omega_r=cfg.omega_c,
        kappa=cfg.kappa,
        g=cfg.g_est,
        epsilon=0.0,
        chi=cfg.chi,
        drive_frequency=drive_frequency,
    )
    readout_chain = ReadoutChain(
        resonator=resonator,
        amplifier=AmplifierChain(
            noise_temperature=cfg.readout_noise_temperature,
            gain=cfg.readout_gain,
        ),
        integration_time=actual_duration,
        dt=cfg.dt,
    )
    disturbance = build_strong_readout_disturbance(
        resonator,
        distorted_waveform,
        dt=cfg.dt,
        spec=cfg.mixing_spec(),
        drive_frequency=drive_frequency,
    )

    pulses: list[Pulse] = [base_pulse]
    drive_ops: dict[str, object] = {"cavity": "cavity"}
    targets = strong_readout_drive_targets(cfg.mixing_spec(), max_transmon_level=model.n_tr)

    if "mix_ge" in targets and np.max(np.abs(disturbance.ge_envelope)) > 0.0:
        ge_carrier = carrier_for_transition_frequency(
            transmon_transition_frequency(model, cavity_level=0, lower_level=0, upper_level=1, frame=frame)
        )
        pulses.append(
            Pulse(
                "mix_ge",
                0.0,
                actual_duration,
                disturbance.ge_envelope,
                carrier=ge_carrier,
                amp=1.0,
                sample_rate=1.0 / cfg.dt,
                label="mix_ge",
            )
        )
        drive_ops["mix_ge"] = targets["mix_ge"]

    if model.n_tr >= 3 and "mix_ef" in targets and np.max(np.abs(disturbance.ef_envelope)) > 0.0:
        ef_carrier = carrier_for_transition_frequency(
            transmon_transition_frequency(model, cavity_level=0, lower_level=1, upper_level=2, frame=frame)
        )
        pulses.append(
            Pulse(
                "mix_ef",
                0.0,
                actual_duration,
                disturbance.ef_envelope,
                carrier=ef_carrier,
                amp=1.0,
                sample_rate=1.0 / cfg.dt,
                label="mix_ef",
            )
        )
        drive_ops["mix_ef"] = targets["mix_ef"]

        for channel, drive_spec in targets.items():
            if not channel.startswith(f"{cfg.mixing_spec().higher_channel_prefix}_"):
                continue
            envelope = disturbance.higher_envelopes.get(channel)
            if envelope is None or np.max(np.abs(envelope)) <= 0.0:
                continue
            carrier = carrier_for_transition_frequency(
                transmon_transition_frequency(
                    model,
                    cavity_level=0,
                    lower_level=drive_spec.lower_level,
                    upper_level=drive_spec.upper_level,
                    frame=frame,
                )
            )
            pulses.append(
                Pulse(
                    channel,
                    0.0,
                    actual_duration,
                    envelope,
                    carrier=carrier,
                    amp=1.0,
                    sample_rate=1.0 / cfg.dt,
                    label=channel,
                )
            )
            drive_ops[channel] = drive_spec

    compiled = compiler.compile(pulses, t_end=actual_duration + cfg.post_buffer)
    return ProtocolBundle(
        model=model,
        frame=frame,
        compiled=compiled,
        drive_ops=drive_ops,
        noise=noise,
        readout_chain=readout_chain,
        distorted_waveform=distorted_waveform,
        disturbance=disturbance,
        duration=actual_duration,
        detuning=float(detuning),
        lowpass_bw_hz=float("nan") if lowpass_bw_hz is None else float(lowpass_bw_hz),
        hardware_rms_error=hardware_rms_error,
        shape=shape,
        amplitude=float(amplitude),
    )


def assignment_fidelity_proxy(protocol: ProtocolBundle) -> float:
    trace_g = protocol.readout_chain.simulate_waveform(
        "g",
        protocol.distorted_waveform,
        dt=protocol.readout_chain.dt,
        include_noise=False,
    )
    trace_e = protocol.readout_chain.simulate_waveform(
        "e",
        protocol.distorted_waveform,
        dt=protocol.readout_chain.dt,
        include_noise=False,
    )
    sigma = protocol.readout_chain.integrated_noise_sigma(duration=protocol.duration, dt=protocol.readout_chain.dt)
    if sigma <= 0.0:
        return 1.0
    distance = float(np.linalg.norm(trace_e.iq_sample - trace_g.iq_sample))
    return float(0.5 * (1.0 + erf(distance / (2.0 * sqrt(2.0) * sigma))))


def _final_population(expectations: dict[str, np.ndarray], label: str) -> float:
    values = np.asarray(expectations.get(label, [0.0]), dtype=float)
    return float(values[-1])


def summarize_expectations(expectations: dict[str, np.ndarray], *, ion_level: int) -> dict[str, float]:
    populations = {key: _final_population(expectations, key) for key in expectations if key.startswith("P_q")}
    p_g = _final_population(expectations, "P_g")
    p_e = _final_population(expectations, "P_e")
    p_f = _final_population(expectations, "P_f")
    p_ion = sum(value for key, value in populations.items() if int(key[3:]) >= ion_level)
    return {
        "p_g": p_g,
        "p_e": p_e,
        "p_f": p_f,
        "p_leak": float(max(0.0, 1.0 - p_g - p_e)),
        "p_ion": float(max(0.0, p_ion)),
        "peak_n_c": float(np.max(np.asarray(expectations["n_c"], dtype=float))),
        "residual_n_c": float(np.asarray(expectations["n_c"], dtype=float)[-1]),
    }


def simulate_point(protocol: ProtocolBundle, cfg: StudyConfig, *, max_step: float | None = None) -> dict[str, object]:
    session = prepare_simulation(
        protocol.model,
        protocol.compiled,
        protocol.drive_ops,
        config=SimulationConfig(frame=protocol.frame, max_step=cfg.dt if max_step is None else float(max_step)),
        noise=protocol.noise,
    )
    g_result = session.run(protocol.model.basis_state(0, 0))
    e_result = session.run(protocol.model.basis_state(1, 0))
    g_repeat = session.run(g_result.final_state)
    e_repeat = session.run(e_result.final_state)

    g_summary = summarize_expectations(g_result.expectations, ion_level=cfg.ion_level)
    e_summary = summarize_expectations(e_result.expectations, ion_level=cfg.ion_level)
    qnd_consistency = 0.5 * (
        _final_population(g_repeat.expectations, "P_g") + _final_population(e_repeat.expectations, "P_e")
    )

    return {
        "amplitude": float(protocol.amplitude),
        "amplitude_mhz": float(protocol.amplitude / (TWO_PI * 1.0e6)),
        "duration": float(protocol.duration),
        "duration_ns": float(protocol.duration * 1.0e9),
        "detuning": float(protocol.detuning),
        "detuning_mhz": float(protocol.detuning / (TWO_PI * 1.0e6)),
        "lowpass_bw_hz": float(protocol.lowpass_bw_hz),
        "lowpass_bw_mhz": float(protocol.lowpass_bw_hz / 1.0e6),
        "shape": protocol.shape,
        "hardware_rms_error": float(protocol.hardware_rms_error),
        "peak_activation": float(protocol.disturbance.peak_activation),
        "peak_mean_occupancy": float(protocol.disturbance.peak_mean_occupancy),
        "assignment_proxy": assignment_fidelity_proxy(protocol),
        "qnd_consistency": float(qnd_consistency),
        "qnd_defect": float(max(0.0, 1.0 - qnd_consistency)),
        "g": g_summary,
        "e": e_summary,
        "mean_p_f": float(0.5 * (g_summary["p_f"] + e_summary["p_f"])),
        "mean_p_leak": float(0.5 * (g_summary["p_leak"] + e_summary["p_leak"])),
        "mean_p_ion": float(0.5 * (g_summary["p_ion"] + e_summary["p_ion"])),
        "mean_peak_n_c": float(0.5 * (g_summary["peak_n_c"] + e_summary["peak_n_c"])),
        "mean_residual_n_c": float(0.5 * (g_summary["residual_n_c"] + e_summary["residual_n_c"])),
    }


def run_continuous_replay(
    protocol: ProtocolBundle,
    cfg: StudyConfig,
    *,
    initial_level: int,
    ntraj: int | None = None,
) -> np.ndarray:
    replay = simulate_continuous_readout(
        protocol.model,
        protocol.compiled,
        protocol.model.basis_state(initial_level, 0),
        protocol.drive_ops,
        noise=protocol.noise,
        e_ops={},
        spec=ContinuousReadoutSpec(
            frame=protocol.frame,
            monitored_subsystem="cavity",
            ntraj=cfg.ntraj_representative if ntraj is None else int(ntraj),
            max_step=cfg.dt,
            keep_runs_results=True,
            store_measurement="end",
        ),
    )
    integrated = []
    for record in replay.measurement_records:
        value = integrate_measurement_record(record, dt=cfg.dt)
        integrated.append(float(np.ravel(np.asarray(value, dtype=float))[0]))
    return np.asarray(integrated, dtype=float)
