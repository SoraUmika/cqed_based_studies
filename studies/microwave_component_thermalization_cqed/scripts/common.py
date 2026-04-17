from __future__ import annotations

import json
import math
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import qutip as qt
from scipy.optimize import curve_fit

from cqed_sim.core import (
    BosonicModeSpec,
    DispersiveReadoutTransmonStorageModel,
    DispersiveTransmonCavityModel,
    ExchangeSpec,
    FrameSpec,
    TransmonModeSpec,
    UniversalCQEDModel,
)
from cqed_sim.observables import reduced_cavity_state, reduced_qubit_state
from cqed_sim.sequence import SequenceCompiler
from cqed_sim.sim import SimulationConfig, simulate_sequence
from cqed_sim.sim.noise import NoiseSpec, collapse_operators


STUDY_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = STUDY_ROOT / "data"
FIGURES_DIR = STUDY_ROOT / "figures"
ARTIFACTS_DIR = STUDY_ROOT / "artifacts"
REPORT_DIR = STUDY_ROOT / "report"

HBAR = 1.054_571_817e-34
KB = 1.380_649e-23
TWO_PI = 2.0 * math.pi


def ghz(value: float) -> float:
    return TWO_PI * float(value) * 1.0e9


def mhz(value: float) -> float:
    return TWO_PI * float(value) * 1.0e6


def bose_occupation(temperature_k: float | np.ndarray, omega_rad_s: float) -> np.ndarray:
    temperature = np.asarray(temperature_k, dtype=float)
    occupation = np.zeros_like(temperature, dtype=float)
    positive = temperature > 0.0
    if np.any(positive):
        exponent = HBAR * float(omega_rad_s) / (KB * temperature[positive])
        occupation[positive] = 1.0 / np.expm1(exponent)
    return occupation


def temperature_from_nth(nth: float | np.ndarray, omega_rad_s: float) -> np.ndarray:
    occupation = np.asarray(nth, dtype=float)
    out = np.zeros_like(occupation, dtype=float)
    positive = occupation > 0.0
    if np.any(positive):
        out[positive] = HBAR * float(omega_rad_s) / (KB * np.log1p(1.0 / occupation[positive]))
    return out


def thermal_variance(nbar: float | np.ndarray) -> np.ndarray:
    arr = np.asarray(nbar, dtype=float)
    return arr * (arr + 1.0)


def ensure_directories() -> None:
    for path in (DATA_DIR, FIGURES_DIR, ARTIFACTS_DIR, REPORT_DIR):
        path.mkdir(parents=True, exist_ok=True)


def save_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)


def save_figure(figure: plt.Figure, stem: str) -> None:
    figure.savefig(FIGURES_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    figure.savefig(FIGURES_DIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(figure)


def configure_matplotlib() -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "legend.fontsize": 9,
            "figure.titlesize": 12,
            "axes.prop_cycle": plt.cycler(color=["#1b5e20", "#0d47a1", "#e65100", "#6a1b9a", "#37474f"]),
        }
    )


@dataclass(frozen=True)
class SingleModeStudyConfig:
    omega_q: float = ghz(5.2)
    omega_r: float = ghz(7.0)
    alpha: float = ghz(-0.22)
    chi: float = mhz(-0.25)
    exchange_g: float = mhz(30.0)
    kappa_readout: float = 1.0 / 150.0e-9
    t1_qubit: float = 60.0e-6
    tphi_qubit: float = 120.0e-6
    n_cav_steady: int = 30
    n_cav_dynamic: int = 12
    n_tr: int = 3
    temperature_grid: tuple[float, ...] = (0.03, 0.05, 0.07, 0.10, 0.15, 0.20, 0.30, 0.50, 1.0, 2.0)
    dephasing_temperatures: tuple[float, ...] = (0.03, 0.05, 0.07, 0.10, 0.15, 0.20, 0.30, 0.50)
    cold_temperature: float = 0.05
    transient_temperatures: tuple[float, ...] = (0.20, 0.50)
    coherence_target: float = 20.0e-6


@dataclass(frozen=True)
class MultimodeStudyConfig:
    omega_q: float = ghz(5.2)
    omega_r: float = ghz(6.95)
    omega_storage_center: float = ghz(7.20)
    alpha: float = ghz(-0.22)
    chi_r: float = mhz(-0.25)
    chi_s: float = mhz(-0.08)
    qubit_readout_exchange: float = mhz(30.0)
    kappa_readout: float = 1.0 / 180.0e-9
    kappa_storage: float = 1.0 / 120.0e-9
    kappa_storage_sensitivity_ns: tuple[float, ...] = (300.0, 120.0, 60.0)
    n_tr: int = 3
    n_readout: int = 6
    n_storage: int = 6
    hot_storage_temperature: float = 0.35
    cold_readout_temperature: float = 0.05
    detuning_grid_mhz: tuple[float, ...] = (-240.0, -180.0, -120.0, -60.0, 0.0, 60.0, 120.0, 180.0, 240.0)
    coupling_grid_mhz: tuple[float, ...] = (0.0, 1.5, 3.0, 4.5, 6.0, 7.5, 9.0, 10.5, 12.0)
    safe_readout_threshold: float = 0.05
    safe_excitation_threshold: float = 0.01


def single_mode_frame(config: SingleModeStudyConfig) -> FrameSpec:
    return FrameSpec(omega_q_frame=config.omega_q, omega_r_frame=config.omega_r)


def multimode_frame(storage_omega: float, config: MultimodeStudyConfig) -> FrameSpec:
    return FrameSpec(omega_c_frame=storage_omega, omega_q_frame=config.omega_q, omega_r_frame=config.omega_r)


def make_dispersive_model(config: SingleModeStudyConfig, *, n_cav: int | None = None, n_tr: int | None = None):
    model = DispersiveTransmonCavityModel(
        omega_c=config.omega_r,
        omega_q=config.omega_q,
        alpha=config.alpha,
        chi=config.chi,
        n_cav=config.n_cav_steady if n_cav is None else int(n_cav),
        n_tr=config.n_tr if n_tr is None else int(n_tr),
    )
    return model, single_mode_frame(config)


def make_dressed_single_mode_model(
    config: SingleModeStudyConfig,
    *,
    n_cav: int | None = None,
    n_tr: int | None = None,
):
    model = UniversalCQEDModel(
        transmon=TransmonModeSpec(
            omega=config.omega_q,
            dim=config.n_tr if n_tr is None else int(n_tr),
            alpha=config.alpha,
        ),
        bosonic_modes=(
            BosonicModeSpec(
                label="readout",
                omega=config.omega_r,
                dim=config.n_cav_steady if n_cav is None else int(n_cav),
                aliases=("readout", "cavity"),
                frame_channel="r",
            ),
        ),
        exchange_terms=(ExchangeSpec("qubit", "readout", config.exchange_g),),
    )
    return model, single_mode_frame(config)


def make_multimode_dispersive_model(
    config: MultimodeStudyConfig,
    *,
    detuning_mhz: float,
    coupling_mhz: float,
    n_storage: int | None = None,
    n_readout: int | None = None,
    n_tr: int | None = None,
):
    storage_omega = config.omega_storage_center + mhz(detuning_mhz)
    model = DispersiveReadoutTransmonStorageModel(
        omega_s=storage_omega,
        omega_r=config.omega_r,
        omega_q=config.omega_q,
        alpha=config.alpha,
        chi_s=config.chi_s,
        chi_r=config.chi_r,
        n_storage=config.n_storage if n_storage is None else int(n_storage),
        n_readout=config.n_readout if n_readout is None else int(n_readout),
        n_tr=config.n_tr if n_tr is None else int(n_tr),
        exchange_terms=(ExchangeSpec("storage", "readout", mhz(coupling_mhz)),),
    )
    return model, multimode_frame(storage_omega, config), storage_omega


def make_multimode_dressed_model(
    config: MultimodeStudyConfig,
    *,
    detuning_mhz: float,
    coupling_mhz: float,
    n_storage: int | None = None,
    n_readout: int | None = None,
):
    storage_omega = config.omega_storage_center + mhz(detuning_mhz)
    model = DispersiveReadoutTransmonStorageModel(
        omega_s=storage_omega,
        omega_r=config.omega_r,
        omega_q=config.omega_q,
        alpha=0.0,
        chi_s=0.0,
        chi_r=0.0,
        n_storage=config.n_storage if n_storage is None else int(n_storage),
        n_readout=config.n_readout if n_readout is None else int(n_readout),
        n_tr=2,
        exchange_terms=(
            ExchangeSpec("storage", "readout", mhz(coupling_mhz)),
            ExchangeSpec("qubit", "readout", config.qubit_readout_exchange),
        ),
    )
    return model, multimode_frame(storage_omega, config), storage_omega


def make_single_mode_noise(config: SingleModeStudyConfig, nth: float, *, dressed: bool = False) -> NoiseSpec:
    if dressed:
        return NoiseSpec(
            kappa_readout=config.kappa_readout,
            nth_readout=float(nth),
            t1=config.t1_qubit,
            tphi=config.tphi_qubit,
        )
    return NoiseSpec(
        kappa=config.kappa_readout,
        nth=float(nth),
        t1=config.t1_qubit,
        tphi=config.tphi_qubit,
    )


def make_multimode_noise(
    config: MultimodeStudyConfig,
    *,
    nth_storage: float,
    kappa_storage: float | None = None,
    nth_readout: float | None = None,
) -> NoiseSpec:
    return NoiseSpec(
        kappa_storage=config.kappa_storage if kappa_storage is None else float(kappa_storage),
        nth_storage=float(nth_storage),
        kappa_readout=config.kappa_readout,
        nth_readout=0.0 if nth_readout is None else float(nth_readout),
        t1=60.0e-6,
        tphi=120.0e-6,
    )


def qubit_populations(state: qt.Qobj, *, n_tr: int) -> dict[str, float]:
    rho_q = reduced_qubit_state(state)
    diagonal = np.real(np.diag(rho_q.full()))
    populations = {"p_g": float(diagonal[0]), "p_e": float(diagonal[1] if len(diagonal) > 1 else 0.0)}
    populations["p_f"] = float(diagonal[2] if n_tr > 2 and len(diagonal) > 2 else 0.0)
    return populations


def qubit_coherence(state: qt.Qobj) -> float:
    rho_q = reduced_qubit_state(state).full()
    return float(abs(rho_q[0, 1]))


def pure_dephasing_ratio(state: qt.Qobj) -> float:
    rho_q = reduced_qubit_state(state).full()
    p_g = max(float(np.real(rho_q[0, 0])), 0.0)
    p_e = max(float(np.real(rho_q[1, 1])), 0.0)
    denom = math.sqrt(max(p_g * p_e, 1.0e-15))
    return float(abs(rho_q[0, 1]) / denom)


def mode_occupation_metrics(state: qt.Qobj, *, alias: str) -> dict[str, float]:
    if alias == "cavity":
        rho_mode = reduced_cavity_state(state)
    elif alias == "readout":
        rho_mode = qt.ptrace(state, 1 if len(state.dims[0]) == 2 else 2)
    elif alias == "storage":
        rho_mode = qt.ptrace(state, 1)
    else:
        raise ValueError(f"Unsupported mode alias '{alias}'.")
    diagonal = np.real(np.diag(rho_mode.full()))
    levels = np.arange(len(diagonal), dtype=float)
    mean_n = float(np.dot(levels, diagonal))
    variance_n = float(np.dot((levels - mean_n) ** 2, diagonal))
    return {
        "mean_n": mean_n,
        "variance_n": variance_n,
        "fock_probabilities": diagonal.tolist(),
    }


def spectroscopy_shift_mhz(chi: float, mean_n: float) -> float:
    return float(abs(chi) * mean_n / TWO_PI / 1.0e6)


def spectroscopy_broadening_mhz(chi: float, variance_n: float) -> float:
    return float(abs(chi) * math.sqrt(max(variance_n, 0.0)) / TWO_PI / 1.0e6)


def readout_noise_penalty(mean_n: float) -> float:
    return float(1.0 / math.sqrt(1.0 + 2.0 * max(mean_n, 0.0)))


def qubit_temperature_resolution_proxy(temperature: np.ndarray, excited_population: np.ndarray) -> np.ndarray:
    derivative = np.gradient(excited_population, temperature, edge_order=1)
    variance = np.clip(excited_population * (1.0 - excited_population), 1.0e-12, None)
    return np.sqrt(variance) / np.clip(np.abs(derivative), 1.0e-12, None)


def cavity_temperature_resolution_proxy(temperature: np.ndarray, mean_n: np.ndarray) -> np.ndarray:
    derivative = np.gradient(mean_n, temperature, edge_order=1)
    variance = np.clip(mean_n * (mean_n + 1.0), 1.0e-12, None)
    return np.sqrt(variance) / np.clip(np.abs(derivative), 1.0e-12, None)


def steady_state(model: Any, *, noise: NoiseSpec, frame: FrameSpec) -> qt.Qobj:
    hamiltonian = model.hamiltonian(frame=frame)
    c_ops = collapse_operators(model, noise)
    return qt.steadystate(hamiltonian, c_ops, method="direct")


def simulate_idle(
    model: Any,
    *,
    duration: float,
    sample_dt: float,
    frame: FrameSpec,
    initial_state: qt.Qobj,
    noise: NoiseSpec,
    max_step: float | None = None,
) -> tuple[np.ndarray, list[qt.Qobj]]:
    compiled = SequenceCompiler(dt=sample_dt).compile([], t_end=duration)
    result = simulate_sequence(
        model,
        compiled,
        initial_state,
        {},
        config=SimulationConfig(frame=frame, store_states=True, max_step=max_step),
        noise=noise,
    )
    states = [] if result.states is None else list(result.states)
    return np.asarray(compiled.tlist, dtype=float), states


def fit_decay_rate(times: np.ndarray, values: np.ndarray) -> dict[str, float]:
    time_axis = np.asarray(times, dtype=float)
    signal = np.asarray(values, dtype=float)
    mask = np.isfinite(signal) & (signal > 1.0e-6)
    if int(np.sum(mask)) < 5:
        return {"gamma": float("nan"), "tau": float("nan"), "amplitude": float("nan")}

    fit_times = time_axis[mask]
    fit_signal = signal[mask]

    def model_fn(t: np.ndarray, amplitude: float, gamma: float) -> np.ndarray:
        return amplitude * np.exp(-gamma * t)

    amplitude0 = float(np.clip(fit_signal[0], 1.0e-6, None))
    gamma0 = float(max(1.0 / max(fit_times[-1], 1.0e-9), 1.0e3))
    try:
        params, _ = curve_fit(
            model_fn,
            fit_times,
            fit_signal,
            p0=(amplitude0, gamma0),
            bounds=([0.0, 0.0], [2.0, np.inf]),
            maxfev=20_000,
        )
        gamma = float(params[1])
        amplitude = float(params[0])
    except RuntimeError:
        slope, intercept = np.polyfit(fit_times, np.log(fit_signal), deg=1)
        gamma = float(max(-slope, 0.0))
        amplitude = float(np.exp(intercept))
    tau = float(np.inf if gamma <= 0.0 else 1.0 / gamma)
    return {"gamma": gamma, "tau": tau, "amplitude": amplitude}


def fit_step_response(times: np.ndarray, values: np.ndarray) -> dict[str, float]:
    time_axis = np.asarray(times, dtype=float)
    signal = np.asarray(values, dtype=float)

    def model_fn(t: np.ndarray, y_inf: float, amplitude: float, tau: float) -> np.ndarray:
        return y_inf - amplitude * np.exp(-t / tau)

    params, _ = curve_fit(
        model_fn,
        time_axis,
        signal,
        p0=(float(signal[-1]), float(signal[-1] - signal[0]), float(max(time_axis[-1] / 5.0, 1.0e-9))),
        bounds=([-np.inf, -np.inf, 1.0e-9], [np.inf, np.inf, np.inf]),
        maxfev=20_000,
    )
    return {"y_inf": float(params[0]), "amplitude": float(params[1]), "tau": float(params[2])}


def elapsed_s(start_time: float) -> float:
    return float(time.perf_counter() - start_time)


def payload_with_runtime(payload: dict[str, Any], *, runtime_s: float) -> dict[str, Any]:
    enriched = dict(payload)
    enriched["runtime_s"] = float(runtime_s)
    return enriched


def format_seconds(value: float) -> str:
    if value >= 1.0e-3:
        return f"{1.0e3 * value:.2f} ms"
    return f"{1.0e6 * value:.2f} us"


def dataclass_payload(instance: Any) -> dict[str, Any]:
    return json.loads(json.dumps(asdict(instance)))
