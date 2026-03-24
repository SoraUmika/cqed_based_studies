"""Shared configuration and physical constants for the readout study."""

from __future__ import annotations

from dataclasses import asdict, dataclass
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
class ReadoutStudyConfig:
    """Nominal physical and numerical settings for the study."""

    omega_q: float = TWO_PI * 6.150e9
    omega_c: float = TWO_PI * 8.597e9
    alpha: float = TWO_PI * (-255.0e6)
    chi: float = TWO_PI * (-2.84e6)
    kappa: float = TWO_PI * 2.4e6
    kerr: float = TWO_PI * (-28.0e3)
    n_tr: int = 3
    n_cav: int = 14
    t1: float = 30.0e-6
    t2: float = 20.0e-6
    eta_nominal: float = 0.35
    eta_ideal: float = 1.0
    dt: float = 4.0e-9
    amp_max: float = TWO_PI * 4.0e6
    duration_grid_ns: tuple[float, ...] = (96.0, 160.0, 240.0, 336.0, 496.0, 720.0)
    representative_duration_ns: float = 240.0
    repeated_wait_factor: float = 2.5
    truncation_probe: tuple[int, ...] = (14, 18)
    chi_probe_scale: tuple[float, ...] = (0.9, 1.0, 1.1)
    kappa_probe_scale: tuple[float, ...] = (0.9, 1.0, 1.1)
    amp_probe_scale: tuple[float, ...] = (0.95, 1.0, 1.05)
    detuning_probe_mhz: tuple[float, ...] = (-0.5, 0.0, 0.5)
    phase_probe_deg: tuple[float, ...] = (-12.0, 0.0, 12.0)
    seed: int = 7

    @property
    def tphi(self) -> float:
        inv_tphi = max(0.0, 1.0 / self.t2 - 1.0 / (2.0 * self.t1))
        return float(np.inf if inv_tphi <= 0.0 else 1.0 / inv_tphi)

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

    def as_dict(self) -> dict[str, float | int | tuple[float, ...]]:
        payload = asdict(self)
        payload["tphi"] = self.tphi
        payload["g_est"] = self.g_est
        payload["n_crit"] = self.n_crit
        payload["midpoint_delta_g"] = self.midpoint_delta_g
        return payload


DEFAULT_CONFIG = ReadoutStudyConfig()
