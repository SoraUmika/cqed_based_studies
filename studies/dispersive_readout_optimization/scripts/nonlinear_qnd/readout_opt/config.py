"""Shared configuration and physical constants for the nonlinear-QND readout study."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

TWO_PI = 2.0 * np.pi

STUDY_DIR = Path(__file__).resolve().parents[2]
DATA_DIR = STUDY_DIR / "data"
FIG_DIR = STUDY_DIR / "figures"
REPORT_DIR = STUDY_DIR / "report"
CQED_SIM_PATH = Path(
    r"C:\Users\dazzl\Box\Shyam Shankar Quantum Circuits Group\Users\Users_JianJun\cQED_simulation"
)


@dataclass(frozen=True)
class HardwareProfile:
    """Hardware transport limits applied through `SequenceCompiler` and study-local slew limiting."""

    gain_i: float = 1.0
    gain_q: float = 0.94
    quadrature_skew_deg: float = 5.0
    dc_i: float = 0.0
    dc_q: float = 0.0
    image_leakage: float = 0.03
    channel_gain: float = 0.97
    zoh_samples: int = 2
    lowpass_bw_hz: float | None = 75.0e6
    detuning_mhz: float = 0.0
    timing_quantum: float | None = 8.0e-9
    amplitude_bits: int | None = 9
    phase_slew_deg_per_ns: float | None = 6.0

    @property
    def quadrature_skew_rad(self) -> float:
        return float(np.deg2rad(self.quadrature_skew_deg))

    @property
    def detuning(self) -> float:
        return float(TWO_PI * self.detuning_mhz * 1.0e6)

    @property
    def phase_slew_rad_per_s(self) -> float | None:
        if self.phase_slew_deg_per_ns is None:
            return None
        return float(np.deg2rad(self.phase_slew_deg_per_ns) * 1.0e9)


@dataclass(frozen=True)
class EffectiveMixingConfig:
    """Phenomenological strong-drive mixing model tied to readout occupancy and waveform slew."""

    onset_ratio: float = 0.10
    occupancy_exponent: float = 1.20
    ge_scale: float = 0.090
    ef_scale: float = 0.050
    slew_ge_scale: float = 0.060
    slew_ef_scale: float = 0.030
    phase_lag_deg: float = -60.0
    slew_phase_deg: float = 35.0

    @property
    def phase_lag_rad(self) -> float:
        return float(np.deg2rad(self.phase_lag_deg))

    @property
    def slew_phase_rad(self) -> float:
        return float(np.deg2rad(self.slew_phase_deg))


@dataclass(frozen=True)
class ReadoutStudyConfig:
    """Nominal physical, transport, and optimization settings for the study."""

    omega_q: float = TWO_PI * 6.150e9
    omega_c: float = TWO_PI * 8.597e9
    alpha: float = TWO_PI * (-255.0e6)
    chi: float = TWO_PI * (-2.84e6)
    kappa: float = TWO_PI * 2.4e6
    kerr: float = TWO_PI * (-28.0e3)
    n_tr: int = 3
    n_cav: int = 14
    t1: float = 30.0e-6
    t1_fe: float = 12.0e-6
    t2: float = 20.0e-6
    nth: float = 0.015
    eta_nominal: float = 0.35
    eta_ideal: float = 1.0
    dt: float = 4.0e-9
    amp_max: float = TWO_PI * 4.0e6
    duration_grid_ns: tuple[float, ...] = (96.0, 240.0, 496.0)
    representative_duration_ns: float = 240.0
    repeated_wait_factor: float = 2.5
    reference_segments: int = 10
    benchmark_qnd_min: float = 0.985
    truncation_probe: tuple[int, ...] = (14, 18)
    chi_probe_scale: tuple[float, ...] = (0.9, 1.0, 1.1)
    kappa_probe_scale: tuple[float, ...] = (0.9, 1.0, 1.1)
    amp_probe_scale: tuple[float, ...] = (0.85, 1.0, 1.15, 1.30)
    detuning_probe_mhz: tuple[float, ...] = (-0.5, 0.0, 0.5)
    phase_probe_deg: tuple[float, ...] = (-12.0, 0.0, 12.0)
    hardware_bandwidth_probe_hz: tuple[float, ...] = (50.0e6, 75.0e6, 120.0e6)
    hardware_skew_probe_deg: tuple[float, ...] = (0.0, 5.0, 9.0)
    hardware_gain_q_probe: tuple[float, ...] = (0.90, 0.94, 0.98)
    hardware_bits_probe: tuple[int, ...] = (7, 9, 11)
    seed: int = 7
    hardware: HardwareProfile = field(default_factory=HardwareProfile)
    mixing: EffectiveMixingConfig = field(default_factory=EffectiveMixingConfig)

    @property
    def tphi(self) -> float:
        inv_tphi = max(0.0, 1.0 / self.t2 - 1.0 / (2.0 * self.t1))
        return float(np.inf if inv_tphi <= 0.0 else 1.0 / inv_tphi)

    @property
    def transmon_t1(self) -> tuple[float, float]:
        return (float(self.t1), float(self.t1_fe))

    @property
    def delta_qc(self) -> float:
        return float(self.omega_q - self.omega_c)

    @property
    def g_est(self) -> float:
        numerator = abs(self.chi) * abs(self.delta_qc) * abs(self.delta_qc + self.alpha)
        denom = abs(self.alpha) if abs(self.alpha) > 1.0e-30 else 1.0
        return float(np.sqrt(max(numerator / denom, 0.0)))

    @property
    def n_crit(self) -> float:
        g = self.g_est
        if g <= 1.0e-30:
            return float(np.inf)
        return float((self.delta_qc / (2.0 * g)) ** 2)

    @property
    def midpoint_delta_g(self) -> float:
        return float(-0.5 * self.chi)

    @property
    def duration_grid(self) -> np.ndarray:
        return 1.0e-9 * np.asarray(self.duration_grid_ns, dtype=float)

    @property
    def representative_duration(self) -> float:
        return float(self.representative_duration_ns * 1.0e-9)

    def drive_frequency(self, delta_g: float) -> float:
        return float(self.omega_c - delta_g)

    def wait_time(self, duration: float | None = None) -> float:
        duration = self.representative_duration if duration is None else float(duration)
        return float(max(duration * 0.4, self.repeated_wait_factor / self.kappa))

    def as_dict(self) -> dict[str, float | int | tuple[float, ...] | dict[str, float | int | None]]:
        payload = asdict(self)
        payload["tphi"] = self.tphi
        payload["transmon_t1"] = list(self.transmon_t1)
        payload["g_est"] = self.g_est
        payload["n_crit"] = self.n_crit
        payload["midpoint_delta_g"] = self.midpoint_delta_g
        payload["hardware"]["quadrature_skew_rad"] = self.hardware.quadrature_skew_rad
        payload["hardware"]["detuning"] = self.hardware.detuning
        payload["hardware"]["phase_slew_rad_per_s"] = self.hardware.phase_slew_rad_per_s
        payload["mixing"]["phase_lag_rad"] = self.mixing.phase_lag_rad
        payload["mixing"]["slew_phase_rad"] = self.mixing.slew_phase_rad
        return payload


DEFAULT_CONFIG = ReadoutStudyConfig()
