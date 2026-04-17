"""Shared helpers for the storage gf-sideband active-cooling study."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import importlib.util
import inspect
import json
import math
from pathlib import Path
import sys
from typing import Callable, Iterable

import numpy as np

from runtime_compat import patch_windows_qutip_import

patch_windows_qutip_import()

import matplotlib.pyplot as plt
import qutip as qt

from physics_and_conventions.conventions import from_internal_units, to_internal_units

from cqed_sim.core.drive_targets import SidebandDriveSpec, TransmonTransitionDriveSpec
from cqed_sim.core.frame import FrameSpec
from cqed_sim.core.frequencies import carrier_for_transition_frequency
from cqed_sim.core.readout_model import DispersiveReadoutTransmonStorageModel
from cqed_sim.pulses.pulse import Pulse
from cqed_sim.sequence.scheduler import SequenceCompiler
from cqed_sim.sim.extractors import (
    cavity_wigner,
    readout_photon_number,
    reduced_storage_state,
    storage_photon_number,
    subsystem_level_population,
    transmon_level_populations,
)
from cqed_sim.sim.noise import NoiseSpec, pure_dephasing_time_from_t1_t2
from cqed_sim.sim.runner import SimulationConfig, simulate_sequence

SCRIPT_DIR = Path(__file__).resolve().parent
STUDY_DIR = SCRIPT_DIR.parent
DATA_DIR = STUDY_DIR / "data"
FIGURES_DIR = STUDY_DIR / "figures"
ARTIFACTS_DIR = STUDY_DIR / "artifacts"
REPORT_DIR = STUDY_DIR / "report"

for directory in (DATA_DIR, FIGURES_DIR, ARTIFACTS_DIR, REPORT_DIR):
    directory.mkdir(parents=True, exist_ok=True)

TWO_PI = 2.0 * np.pi
DEFAULT_DT_S = 0.25e-9
DEFAULT_RINGDOWN_MULTIPLE = 4.0


@dataclass(frozen=True)
class DeviceParameters:
    """Exact device values available in the local editable cqed_sim environment."""

    parameter_source: str
    readout_frequency_hz: float
    qubit_frequency_hz: float
    storage_frequency_hz: float
    readout_kappa_hz: float
    qubit_anharmonicity_hz: float
    chi_storage_hz: float
    chi_readout_hz: float
    storage_gf_sideband_nominal_hz: float
    storage_t1_s: float
    storage_t2_ramsey_s: float
    chi_storage_readout_hz: float = 0.0
    storage_kerr_hz: float = 0.0
    readout_kerr_hz: float = 0.0


def _fallback_device() -> DeviceParameters:
    return DeviceParameters(
        parameter_source="cQED_simulation/examples/sequential_sideband_reset.py",
        readout_frequency_hz=8_596_222_556.078796,
        qubit_frequency_hz=6_150_358_764.4830475,
        storage_frequency_hz=5_240_932_800.0,
        readout_kappa_hz=4.156e6,
        qubit_anharmonicity_hz=-255_669_694.5244608,
        chi_storage_hz=-2_840_421.0,
        chi_readout_hz=-3.0e6,
        storage_gf_sideband_nominal_hz=6_803_533_628.0,
        storage_t1_s=250.0e-6,
        storage_t2_ramsey_s=150.0e-6,
    )


def _load_device_from_local_example() -> DeviceParameters:
    """Load the exact sideband-reset device tuple from the editable cqed_sim tree."""

    try:
        import cqed_sim

        sim_root = Path(inspect.getfile(cqed_sim)).resolve().parents[1]
        example_path = sim_root / "examples" / "sequential_sideband_reset.py"
        if not example_path.exists():
            return _fallback_device()
        spec = importlib.util.spec_from_file_location("cqed_sim_sideband_reset_example", example_path)
        if spec is None or spec.loader is None:
            return _fallback_device()
        module = importlib.util.module_from_spec(spec)
        if str(sim_root) not in sys.path:
            sys.path.insert(0, str(sim_root))
        spec.loader.exec_module(module)
        if not hasattr(module, "_device"):
            return _fallback_device()
        raw_device = module._device()
        return DeviceParameters(
            parameter_source=str(example_path),
            readout_frequency_hz=float(raw_device.readout_frequency_hz),
            qubit_frequency_hz=float(raw_device.qubit_frequency_hz),
            storage_frequency_hz=float(raw_device.storage_frequency_hz),
            readout_kappa_hz=float(raw_device.readout_kappa_hz),
            qubit_anharmonicity_hz=float(raw_device.qubit_anharmonicity_hz),
            chi_storage_hz=float(raw_device.chi_storage_hz),
            chi_readout_hz=float(raw_device.chi_readout_hz),
            storage_gf_sideband_nominal_hz=float(raw_device.storage_gf_sideband_frequency_hz or 0.0),
            storage_t1_s=float(raw_device.storage_t1_s or 0.0),
            storage_t2_ramsey_s=float(raw_device.storage_t2_ramsey_s or 0.0),
            chi_storage_readout_hz=float(raw_device.chi_storage_readout_hz),
            storage_kerr_hz=float(raw_device.storage_kerr_hz),
            readout_kerr_hz=float(raw_device.readout_kerr_hz),
        )
    except Exception:
        return _fallback_device()


DEVICE = _load_device_from_local_example()

N_TR = 4
N_STORAGE = 7
N_READOUT = 3


EnvelopeFunc = Callable[[np.ndarray], np.ndarray]


def _normalized(values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    area = float(np.trapezoid(np.abs(values), grid))
    if abs(area) < 1.0e-15:
        return np.asarray(values, dtype=np.complex128)
    return np.asarray(values / area, dtype=np.complex128)


def square_envelope(t_rel: np.ndarray) -> np.ndarray:
    return np.ones_like(np.asarray(t_rel, dtype=float), dtype=np.complex128)


def gaussian_envelope(sigma_fraction: float = 0.18) -> EnvelopeFunc:
    sigma = float(sigma_fraction)
    grid = np.linspace(0.0, 1.0, 4001)
    base = np.exp(-0.5 * ((grid - 0.5) / sigma) ** 2)
    normalized = _normalized(base, grid)

    def envelope(t_rel: np.ndarray) -> np.ndarray:
        t = np.asarray(t_rel, dtype=float)
        return np.interp(t, grid, np.real(normalized)).astype(np.complex128)

    return envelope


def cosine_squared_envelope() -> EnvelopeFunc:
    grid = np.linspace(0.0, 1.0, 4001)
    base = np.cos(np.pi * (grid - 0.5)) ** 2
    normalized = _normalized(base, grid)

    def envelope(t_rel: np.ndarray) -> np.ndarray:
        t = np.asarray(t_rel, dtype=float)
        return np.interp(t, grid, np.real(normalized)).astype(np.complex128)

    return envelope


def flat_top_gaussian_envelope(ramp_fraction: float = 0.18) -> EnvelopeFunc:
    ramp = max(0.02, min(float(ramp_fraction), 0.45))
    sigma = ramp / 3.0
    grid = np.linspace(0.0, 1.0, 4001)
    base = np.ones_like(grid)
    left = grid < ramp
    right = grid > (1.0 - ramp)
    base[left] = np.exp(-0.5 * ((grid[left] - ramp) / sigma) ** 2)
    base[right] = np.exp(-0.5 * ((grid[right] - (1.0 - ramp)) / sigma) ** 2)
    normalized = _normalized(base, grid)

    def envelope(t_rel: np.ndarray) -> np.ndarray:
        t = np.asarray(t_rel, dtype=float)
        return np.interp(t, grid, np.real(normalized)).astype(np.complex128)

    return envelope


def bump_envelope() -> EnvelopeFunc:
    grid = np.linspace(0.0, 1.0, 4001)
    safe = np.clip(grid, 1.0e-4, 1.0 - 1.0e-4)
    base = np.exp(-1.0 / (safe * (1.0 - safe)))
    base[(grid <= 0.0) | (grid >= 1.0)] = 0.0
    normalized = _normalized(base, grid)

    def envelope(t_rel: np.ndarray) -> np.ndarray:
        t = np.asarray(t_rel, dtype=float)
        return np.interp(t, grid, np.real(normalized)).astype(np.complex128)

    return envelope


def phase_modulated_bump_envelope(modulation_index: float = 0.85) -> EnvelopeFunc:
    base_env = bump_envelope()
    grid = np.linspace(0.0, 1.0, 4001)
    base = base_env(grid) * np.exp(1j * float(modulation_index) * np.sin(2.0 * np.pi * (grid - 0.5)))
    normalized = _normalized(base, grid)

    def envelope(t_rel: np.ndarray) -> np.ndarray:
        t = np.asarray(t_rel, dtype=float)
        real_part = np.interp(t, grid, np.real(normalized))
        imag_part = np.interp(t, grid, np.imag(normalized))
        return (real_part + 1j * imag_part).astype(np.complex128)

    return envelope


def bb1_like_envelope() -> EnvelopeFunc:
    """A fixed-phase composite envelope with the standard BB1 phase pattern."""

    phi = float(np.arccos(-0.25))
    weights = np.asarray([1.0, 1.0, 2.0, 1.0], dtype=float)
    phases = np.asarray([0.0, phi, 3.0 * phi, phi], dtype=float)
    cuts = np.cumsum(weights) / np.sum(weights)
    grid = np.linspace(0.0, 1.0, 4001)
    base = np.zeros_like(grid, dtype=np.complex128)
    start = 0.0
    for stop, phase in zip(cuts, phases, strict=True):
        mask = (grid >= start) & (grid < stop)
        base[mask] = np.exp(1j * phase)
        start = float(stop)
    base[grid >= cuts[-1]] = np.exp(1j * phases[-1])
    normalized = _normalized(base, grid)

    def envelope(t_rel: np.ndarray) -> np.ndarray:
        t = np.asarray(t_rel, dtype=float)
        real_part = np.interp(t, grid, np.real(normalized))
        imag_part = np.interp(t, grid, np.imag(normalized))
        return (real_part + 1j * imag_part).astype(np.complex128)

    return envelope


ENVELOPE_BUILDERS: dict[str, Callable[[], EnvelopeFunc]] = {
    "square": lambda: square_envelope,
    "gaussian": lambda: gaussian_envelope(0.18),
    "cosine_squared": cosine_squared_envelope,
    "flat_top_gaussian": lambda: flat_top_gaussian_envelope(0.2),
    "bump": bump_envelope,
    "phase_modulated_bump": lambda: phase_modulated_bump_envelope(0.85),
    "bb1_like": bb1_like_envelope,
}


def build_model(
    *,
    device: DeviceParameters = DEVICE,
    n_storage: int = N_STORAGE,
    n_readout: int = N_READOUT,
    n_tr: int = N_TR,
) -> DispersiveReadoutTransmonStorageModel:
    return DispersiveReadoutTransmonStorageModel(
        omega_s=to_internal_units(device.storage_frequency_hz),
        omega_r=to_internal_units(device.readout_frequency_hz),
        omega_q=to_internal_units(device.qubit_frequency_hz),
        alpha=to_internal_units(device.qubit_anharmonicity_hz),
        chi_s=to_internal_units(device.chi_storage_hz),
        chi_r=to_internal_units(device.chi_readout_hz),
        chi_sr=to_internal_units(device.chi_storage_readout_hz),
        kerr_s=to_internal_units(device.storage_kerr_hz),
        kerr_r=to_internal_units(device.readout_kerr_hz),
        n_storage=int(n_storage),
        n_readout=int(n_readout),
        n_tr=int(n_tr),
    )


def build_frame(model: DispersiveReadoutTransmonStorageModel) -> FrameSpec:
    return FrameSpec(
        omega_c_frame=float(model.omega_s),
        omega_q_frame=float(model.omega_q),
        omega_r_frame=float(model.omega_r),
    )


def build_noise(device: DeviceParameters = DEVICE) -> NoiseSpec:
    tphi_storage = pure_dephasing_time_from_t1_t2(t1_s=device.storage_t1_s, t2_s=device.storage_t2_ramsey_s)
    return NoiseSpec(
        kappa_storage=None if device.storage_t1_s <= 0.0 else 1.0 / float(device.storage_t1_s),
        tphi_storage=tphi_storage,
        kappa_readout=float(device.readout_kappa_hz),
        nth_storage=0.0,
        nth_readout=0.0,
    )


def json_dump(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def csv_dump(path: Path, rows: Iterable[dict[str, object]]) -> None:
    rows = list(rows)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    header = list(rows[0].keys())
    lines = [",".join(header)]
    for row in rows:
        lines.append(",".join(str(row[key]) for key in header))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def hz(value_rad_s: float) -> float:
    return float(from_internal_units(float(value_rad_s)))


def mhz(value_rad_s: float) -> float:
    return hz(value_rad_s) / 1.0e6


def ghz(value_rad_s: float) -> float:
    return hz(value_rad_s) / 1.0e9


def basis_label(q_level: int, storage_level: int, readout_level: int = 0) -> str:
    labels = {0: "g", 1: "e", 2: "f", 3: "h"}
    q_text = labels.get(int(q_level), f"q{q_level}")
    # Internally the tensor order is (transmon, storage, readout), but the study
    # reports states in the user-requested |q, n_r, n_s> convention.
    return f"|{q_text},{int(readout_level)},{int(storage_level)}>"


def state_population(state: qt.Qobj, basis_state: qt.Qobj) -> float:
    projector = basis_state.proj()
    if state.isket:
        return float(abs(basis_state.overlap(state)) ** 2)
    return float(np.real((projector * state).tr()))


def transmon_excited_population(state: qt.Qobj) -> float:
    populations = transmon_level_populations(state)
    return float(sum(population for level, population in populations.items() if int(level) >= 1))


def sorted_basis_populations(
    state: qt.Qobj,
    model: DispersiveReadoutTransmonStorageModel,
    *,
    cutoff: int = 8,
) -> list[tuple[str, float]]:
    entries: list[tuple[str, float]] = []
    for q_level in range(model.n_tr):
        for storage_level in range(model.n_storage):
            for readout_level in range(model.n_readout):
                ket = model.basis_state(q_level, storage_level, readout_level)
                pop = state_population(state, ket)
                if pop > 1.0e-8:
                    entries.append((basis_label(q_level, storage_level, readout_level), pop))
    entries.sort(key=lambda item: item[1], reverse=True)
    return entries[:cutoff]


def storage_sideband_lab_frequency(model: DispersiveReadoutTransmonStorageModel, n: int) -> float:
    return float(model.basis_energy(2, int(n) - 1, 0) - model.basis_energy(0, int(n), 0))


def readout_dump_lab_frequency(model: DispersiveReadoutTransmonStorageModel, n: int) -> float:
    return float(model.basis_energy(2, int(n) - 1, 0) - model.basis_energy(0, int(n) - 1, 1))


def storage_sideband_rotating_frequency(
    model: DispersiveReadoutTransmonStorageModel,
    frame: FrameSpec,
    n: int,
) -> float:
    return float(
        model.sideband_transition_frequency(
            mode="storage",
            storage_level=int(n) - 1,
            readout_level=0,
            lower_level=0,
            upper_level=2,
            sideband="red",
            frame=frame,
        )
    )


def readout_dump_rotating_frequency(
    model: DispersiveReadoutTransmonStorageModel,
    frame: FrameSpec,
    n: int,
) -> float:
    return float(
        model.sideband_transition_frequency(
            mode="readout",
            storage_level=int(n) - 1,
            readout_level=0,
            lower_level=0,
            upper_level=2,
            sideband="red",
            frame=frame,
        )
    )


def sideband_matrix_element(model: DispersiveReadoutTransmonStorageModel, *, mode: str, n: int) -> float:
    _, lowering = model.sideband_drive_operators(mode=mode, lower_level=0, upper_level=2, sideband="red")
    if mode == "storage":
        source = model.basis_state(2, int(n) - 1, 0)
        target = model.basis_state(0, int(n), 0)
    elif mode == "readout":
        source = model.basis_state(2, int(n) - 1, 0)
        target = model.basis_state(0, int(n) - 1, 1)
    else:
        raise ValueError(f"Unsupported mode '{mode}'.")
    return float(abs(target.overlap(lowering * source)))


def expected_pi_duration_s(*, amplitude_hz: float, matrix_element: float) -> float:
    amplitude_rad_s = to_internal_units(float(amplitude_hz))
    return float(np.pi / (2.0 * amplitude_rad_s * float(matrix_element)))


def make_pulse(
    *,
    channel: str,
    carrier_rad_s: float,
    duration_s: float,
    amplitude_hz: float,
    envelope_name: str,
    phase_rad: float = 0.0,
    t0_s: float = 0.0,
    drag: float = 0.0,
    label: str | None = None,
) -> Pulse:
    if envelope_name not in ENVELOPE_BUILDERS:
        raise ValueError(f"Unsupported envelope family '{envelope_name}'.")
    envelope = ENVELOPE_BUILDERS[envelope_name]()
    return Pulse(
        channel=str(channel),
        t0=float(t0_s),
        duration=float(duration_s),
        envelope=envelope,
        amp=float(to_internal_units(amplitude_hz)),
        carrier=float(carrier_for_transition_frequency(carrier_rad_s)),
        phase=float(phase_rad),
        drag=float(drag),
        label=label,
    )


def compile_single_pulse(pulse: Pulse | None, *, duration_s: float, dt_s: float = DEFAULT_DT_S):
    pulses = [] if pulse is None else [pulse]
    t_end = max(float(duration_s), float(dt_s))
    return SequenceCompiler(dt=float(dt_s)).compile(pulses, t_end=t_end)


def simulate_single_stage(
    model: DispersiveReadoutTransmonStorageModel,
    initial_state: qt.Qobj,
    *,
    pulse: Pulse | None,
    duration_s: float,
    drive_ops: dict[str, object],
    frame: FrameSpec,
    noise: NoiseSpec | None = None,
    dt_s: float = DEFAULT_DT_S,
    store_states: bool = True,
):
    compiled = compile_single_pulse(pulse, duration_s=duration_s, dt_s=dt_s)
    result = simulate_sequence(
        model,
        compiled,
        initial_state,
        drive_ops,
        SimulationConfig(frame=frame, store_states=store_states, max_step=dt_s),
        noise=noise,
    )
    return compiled, result


def ladder_sequence_times(noise: NoiseSpec) -> float:
    if noise.kappa_readout is None or noise.kappa_readout <= 0.0:
        return 0.0
    return float(DEFAULT_RINGDOWN_MULTIPLE / float(noise.kappa_readout))


def plot_save(fig: plt.Figure, stem: str) -> None:
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def write_device_manifest(path: Path) -> None:
    json_dump(
        path,
        {
            "device_parameters": asdict(DEVICE),
            "internal_units": {
                "frequencies": "rad/s",
                "times": "s",
                "noise_rates": "1/s",
            },
            "tensor_ordering": ["transmon", "storage", "readout"],
            "report_state_ordering": ["transmon", "readout", "storage"],
        },
    )
