"""Shared helpers for the gf-sideband waveform optimization study."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import csv
import importlib.util
import inspect
import json
import math
from pathlib import Path
import platform
import sys
from typing import Callable

from runtime_compat import patch_windows_qutip_import

patch_windows_qutip_import()

import matplotlib.pyplot as plt
import numpy as np
import qutip as qt

from cqed_sim.core.drive_targets import SidebandDriveSpec
from cqed_sim.core.frame import FrameSpec
from cqed_sim.core.frequencies import carrier_for_transition_frequency
from cqed_sim.core.readout_model import DispersiveReadoutTransmonStorageModel
from cqed_sim.pulses.pulse import Pulse
from cqed_sim.sequence.scheduler import SequenceCompiler
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
DEFAULT_SWEEP_DT_S = 0.5e-9
DEFAULT_FINAL_DT_S = 0.25e-9
MODES = ("storage", "readout")
N_VALUES = (1, 2, 3)
Q_LABELS = {0: "g", 1: "e", 2: "f", 3: "h"}


def to_internal_units(value_hz: float) -> float:
    return TWO_PI * float(value_hz)


def from_internal_units(value_rad_s: float) -> float:
    return float(value_rad_s) / TWO_PI


def hz(value_rad_s: float) -> float:
    return from_internal_units(value_rad_s)


def mhz(value_rad_s: float) -> float:
    return hz(value_rad_s) / 1.0e6


def ghz(value_rad_s: float) -> float:
    return hz(value_rad_s) / 1.0e9


def trapezoid(values: np.ndarray, grid: np.ndarray) -> float:
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(values, grid))
    return float(np.trapz(values, grid))


@dataclass(frozen=True)
class DeviceParameters:
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


@dataclass(frozen=True)
class TransmonCoherenceReference:
    parameter_source: str
    qubit_t1_s: float | None
    qubit_t2_ramsey_s: float | None
    qubit_t2_echo_s: float | None = None
    qubit_tphi_ramsey_s: float | None = None
    notes: str = ""


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


def load_device_from_local_example() -> DeviceParameters:
    """Load the exact local device tuple from the editable cqed_sim tree."""
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


DEVICE = load_device_from_local_example()


def _fallback_transmon_reference() -> TransmonCoherenceReference:
    qubit_t1_s = 9812.873848245112e-9
    qubit_t2_ramsey_s = 6324.73112712837e-9
    qubit_t2_echo_s = 8381.0e-9
    return TransmonCoherenceReference(
        parameter_source="cQED_simulation/examples/workflows/simulate_fock_tomo_and_sqr_calibration.py",
        qubit_t1_s=qubit_t1_s,
        qubit_t2_ramsey_s=qubit_t2_ramsey_s,
        qubit_t2_echo_s=qubit_t2_echo_s,
        qubit_tphi_ramsey_s=pure_dephasing_time_from_t1_t2(t1_s=qubit_t1_s, t2_s=qubit_t2_ramsey_s),
        notes=(
            "Matched local tomography workflow sharing the same storage/readout frequencies and "
            "transmon anharmonicity scale as the sideband-reset example. The sideband-reset device "
            "itself does not include transmon coherence values, so this reference is used as an "
            "explicitly labeled sensitivity anchor rather than as a claimed ground-truth measurement."
        ),
    )


def load_transmon_reference_from_local_example() -> TransmonCoherenceReference:
    try:
        import cqed_sim

        sim_root = Path(inspect.getfile(cqed_sim)).resolve().parents[1]
        reference_path = sim_root / "examples" / "workflows" / "simulate_fock_tomo_and_sqr_calibration.py"
        if not reference_path.exists():
            return _fallback_transmon_reference()
        spec = importlib.util.spec_from_file_location("cqed_sim_transmon_reference_example", reference_path)
        if spec is None or spec.loader is None:
            return _fallback_transmon_reference()
        module = importlib.util.module_from_spec(spec)
        if str(sim_root) not in sys.path:
            sys.path.insert(0, str(sim_root))
        spec.loader.exec_module(module)
        if not hasattr(module, "DEVICE_CONFIG"):
            return _fallback_transmon_reference()
        config = module.DEVICE_CONFIG
        qubit_t1_s = float(config["qb_T1_relax"]) * 1.0e-9 if config.get("qb_T1_relax") is not None else None
        qubit_t2_ramsey_s = float(config["qb_T2_ramsey"]) * 1.0e-9 if config.get("qb_T2_ramsey") is not None else None
        qubit_t2_echo_s = float(config["qb_T2_echo"]) * 1.0e-9 if config.get("qb_T2_echo") is not None else None
        return TransmonCoherenceReference(
            parameter_source=str(reference_path),
            qubit_t1_s=qubit_t1_s,
            qubit_t2_ramsey_s=qubit_t2_ramsey_s,
            qubit_t2_echo_s=qubit_t2_echo_s,
            qubit_tphi_ramsey_s=pure_dephasing_time_from_t1_t2(t1_s=qubit_t1_s, t2_s=qubit_t2_ramsey_s),
            notes=(
                "Loaded from the local tomography workflow that reuses the same device-frequency scale. "
                "Used here as a transmon-decoherence reference because the sideband-reset example exports "
                "storage/readout decay but no transmon T1/T2 values."
            ),
        )
    except Exception:
        return _fallback_transmon_reference()


TRANSMON_REFERENCE = load_transmon_reference_from_local_example()


@dataclass(frozen=True)
class FamilyVariant:
    family: str
    variant: str
    envelope_name: str
    parameter: float | None = None
    description: str = ""


FAMILY_VARIANTS: tuple[FamilyVariant, ...] = (
    FamilyVariant("square", "square", "square", description="Constant envelope"),
    FamilyVariant("gaussian", "sigma_0p12", "gaussian", 0.12, "Narrow Gaussian"),
    FamilyVariant("gaussian", "sigma_0p18", "gaussian", 0.18, "Reference Gaussian"),
    FamilyVariant("gaussian", "sigma_0p24", "gaussian", 0.24, "Wide Gaussian"),
    FamilyVariant("cosine", "raised_cosine", "cosine", description="Single-lobe cosine"),
    FamilyVariant("flat_top_cosine", "ramp_0p15", "flat_top_cosine", 0.15, "Short cosine ramp"),
    FamilyVariant("flat_top_cosine", "ramp_0p25", "flat_top_cosine", 0.25, "Reference cosine ramp"),
    FamilyVariant("flat_top_cosine", "ramp_0p35", "flat_top_cosine", 0.35, "Long cosine ramp"),
    FamilyVariant("flat_top_gaussian", "ramp_0p15", "flat_top_gaussian", 0.15, "Short Gaussian ramp"),
    FamilyVariant("flat_top_gaussian", "ramp_0p25", "flat_top_gaussian", 0.25, "Reference Gaussian ramp"),
    FamilyVariant("flat_top_gaussian", "ramp_0p35", "flat_top_gaussian", 0.35, "Long Gaussian ramp"),
    FamilyVariant("smooth_bump", "compact_bump", "smooth_bump", description="Compact-support smooth ramp"),
    FamilyVariant("blackman", "blackman", "blackman", description="Blackman window"),
)

FAMILY_ORDER: tuple[str, ...] = (
    "square",
    "gaussian",
    "cosine",
    "flat_top_cosine",
    "flat_top_gaussian",
    "smooth_bump",
    "blackman",
)

FAMILY_DESCRIPTIONS: dict[str, str] = {
    "square": "Hard-edged baseline with maximal bandwidth.",
    "gaussian": "Smooth localized pulse with tunable sigma.",
    "cosine": "Raised-cosine pulse with smooth zero-valued edges and no flat top.",
    "flat_top_cosine": "Flat-top pulse with cosine ramps.",
    "flat_top_gaussian": "Flat-top pulse with Gaussian ramps.",
    "smooth_bump": "Compact-support smooth pulse with vanishing derivatives at both edges.",
    "blackman": "Windowed pulse that strongly suppresses sidelobes.",
}


def variant_key(mode: str, n: int, family: str, duration_ns: float) -> str:
    return f"{mode}_n{n}_{family}_{duration_ns:.1f}ns"


def basis_label(q_level: int, storage_level: int, readout_level: int) -> str:
    q_text = Q_LABELS.get(int(q_level), f"q{q_level}")
    return f"|{q_text},{int(readout_level)},{int(storage_level)}>"


def basis_state_for_mode(model: DispersiveReadoutTransmonStorageModel, mode: str, n: int) -> tuple[qt.Qobj, qt.Qobj]:
    if mode == "storage":
        return model.basis_state(2, int(n) - 1, 0), model.basis_state(0, int(n), 0)
    if mode == "readout":
        return model.basis_state(2, 0, int(n) - 1), model.basis_state(0, 0, int(n))
    raise ValueError(f"Unsupported mode '{mode}'.")


def sideband_lab_frequency(model: DispersiveReadoutTransmonStorageModel, mode: str, n: int) -> float:
    source, target = basis_state_for_mode(model, mode, n)
    source_indices = tuple(int(index) for index in source.dims[0])  # pragma: no cover
    target_indices = tuple(int(index) for index in target.dims[0])  # pragma: no cover
    del source_indices, target_indices
    if mode == "storage":
        return float(model.basis_energy(2, int(n) - 1, 0) - model.basis_energy(0, int(n), 0))
    return float(model.basis_energy(2, 0, int(n) - 1) - model.basis_energy(0, 0, int(n)))


def sideband_rotating_frequency(
    model: DispersiveReadoutTransmonStorageModel,
    frame: FrameSpec,
    mode: str,
    n: int,
) -> float:
    if mode == "storage":
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
    if mode == "readout":
        return float(
            model.sideband_transition_frequency(
                mode="readout",
                storage_level=0,
                readout_level=int(n) - 1,
                lower_level=0,
                upper_level=2,
                sideband="red",
                frame=frame,
            )
        )
    raise ValueError(f"Unsupported mode '{mode}'.")


def analytic_rotating_sideband_frequency(mode: str, n: int, device: DeviceParameters = DEVICE) -> float:
    if mode == "storage":
        return float(device.qubit_anharmonicity_hz + 2.0 * device.chi_storage_hz * (int(n) - 1) - device.storage_kerr_hz * (int(n) - 1))
    if mode == "readout":
        return float(device.qubit_anharmonicity_hz + 2.0 * device.chi_readout_hz * (int(n) - 1) - device.readout_kerr_hz * (int(n) - 1))
    raise ValueError(f"Unsupported mode '{mode}'.")


def sideband_matrix_element(model: DispersiveReadoutTransmonStorageModel, mode: str, n: int) -> float:
    _, lowering = model.sideband_drive_operators(mode=mode, lower_level=0, upper_level=2, sideband="red")
    source, target = basis_state_for_mode(model, mode, n)
    return float(abs(target.overlap(lowering * source)))


def expected_pi_amplitude_hz(duration_s: float, matrix_element: float) -> float:
    return float(1.0 / (4.0 * float(duration_s) * float(matrix_element)))


def expected_pi_duration_s(amplitude_hz: float, matrix_element: float) -> float:
    amplitude_rad_s = to_internal_units(float(amplitude_hz))
    return float(np.pi / (2.0 * amplitude_rad_s * float(matrix_element)))


def _normalize_base_envelope(values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    area = trapezoid(np.abs(values), grid)
    if abs(area) < 1.0e-15:
        return np.asarray(values, dtype=np.complex128)
    return np.asarray(values / area, dtype=np.complex128)


def _interp_real(grid: np.ndarray, values: np.ndarray, t_rel: np.ndarray) -> np.ndarray:
    return np.interp(np.asarray(t_rel, dtype=float), grid, np.asarray(values, dtype=float)).astype(np.complex128)


def _build_square() -> Callable[[np.ndarray], np.ndarray]:
    return lambda t_rel: np.ones_like(np.asarray(t_rel, dtype=float), dtype=np.complex128)


def _build_gaussian(sigma_fraction: float) -> Callable[[np.ndarray], np.ndarray]:
    grid = np.linspace(0.0, 1.0, 4001)
    base = np.exp(-0.5 * ((grid - 0.5) / float(sigma_fraction)) ** 2)
    normalized = _normalize_base_envelope(base, grid)
    return lambda t_rel: _interp_real(grid, np.real(normalized), t_rel)


def _build_cosine() -> Callable[[np.ndarray], np.ndarray]:
    grid = np.linspace(0.0, 1.0, 4001)
    base = np.sin(np.pi * grid)
    normalized = _normalize_base_envelope(base, grid)
    return lambda t_rel: _interp_real(grid, np.real(normalized), t_rel)


def _build_flat_top_cosine(ramp_fraction: float) -> Callable[[np.ndarray], np.ndarray]:
    ramp = max(0.02, min(float(ramp_fraction), 0.45))
    grid = np.linspace(0.0, 1.0, 4001)
    base = np.ones_like(grid)
    left = grid < ramp
    right = grid > (1.0 - ramp)
    base[left] = 0.5 * (1.0 - np.cos(np.pi * grid[left] / ramp))
    base[right] = 0.5 * (1.0 - np.cos(np.pi * (1.0 - grid[right]) / ramp))
    normalized = _normalize_base_envelope(base, grid)
    return lambda t_rel: _interp_real(grid, np.real(normalized), t_rel)


def _build_flat_top_gaussian(ramp_fraction: float) -> Callable[[np.ndarray], np.ndarray]:
    ramp = max(0.02, min(float(ramp_fraction), 0.45))
    sigma = ramp / 3.0
    grid = np.linspace(0.0, 1.0, 4001)
    base = np.ones_like(grid)
    left = grid < ramp
    right = grid > (1.0 - ramp)
    base[left] = np.exp(-0.5 * ((grid[left] - ramp) / sigma) ** 2)
    base[right] = np.exp(-0.5 * ((grid[right] - (1.0 - ramp)) / sigma) ** 2)
    normalized = _normalize_base_envelope(base, grid)
    return lambda t_rel: _interp_real(grid, np.real(normalized), t_rel)


def _build_smooth_bump() -> Callable[[np.ndarray], np.ndarray]:
    grid = np.linspace(0.0, 1.0, 4001)
    safe = np.clip(grid, 1.0e-4, 1.0 - 1.0e-4)
    base = np.exp(-1.0 / (safe * (1.0 - safe)))
    base[(grid <= 0.0) | (grid >= 1.0)] = 0.0
    normalized = _normalize_base_envelope(base, grid)
    return lambda t_rel: _interp_real(grid, np.real(normalized), t_rel)


def _build_blackman() -> Callable[[np.ndarray], np.ndarray]:
    grid = np.linspace(0.0, 1.0, 4001)
    base = np.blackman(grid.size)
    base = np.clip(base, 0.0, None)
    normalized = _normalize_base_envelope(base, grid)
    return lambda t_rel: _interp_real(grid, np.real(normalized), t_rel)


def build_envelope(variant: FamilyVariant) -> Callable[[np.ndarray], np.ndarray]:
    if variant.envelope_name == "square":
        return _build_square()
    if variant.envelope_name == "gaussian":
        assert variant.parameter is not None
        return _build_gaussian(float(variant.parameter))
    if variant.envelope_name == "cosine":
        return _build_cosine()
    if variant.envelope_name == "flat_top_cosine":
        assert variant.parameter is not None
        return _build_flat_top_cosine(float(variant.parameter))
    if variant.envelope_name == "flat_top_gaussian":
        assert variant.parameter is not None
        return _build_flat_top_gaussian(float(variant.parameter))
    if variant.envelope_name == "smooth_bump":
        return _build_smooth_bump()
    if variant.envelope_name == "blackman":
        return _build_blackman()
    raise ValueError(f"Unsupported envelope '{variant.envelope_name}'.")


def build_model(
    *,
    device: DeviceParameters = DEVICE,
    n_storage: int = 5,
    n_readout: int = 5,
    n_tr: int = 4,
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


def build_noise(
    device: DeviceParameters = DEVICE,
    *,
    transmon_t1_s: float | None = None,
    transmon_t2_ramsey_s: float | None = None,
    transmon_tphi_s: float | None = None,
) -> NoiseSpec:
    tphi_storage = pure_dephasing_time_from_t1_t2(t1_s=device.storage_t1_s, t2_s=device.storage_t2_ramsey_s)
    if transmon_tphi_s is None and transmon_t2_ramsey_s is not None:
        transmon_tphi_s = pure_dephasing_time_from_t1_t2(t1_s=transmon_t1_s, t2_s=transmon_t2_ramsey_s)
    return NoiseSpec(
        t1=transmon_t1_s,
        tphi=transmon_tphi_s,
        kappa_storage=None if device.storage_t1_s <= 0.0 else 1.0 / float(device.storage_t1_s),
        tphi_storage=tphi_storage,
        kappa_readout=float(device.readout_kappa_hz),
        nth_storage=0.0,
        nth_readout=0.0,
    )


def device_parameter_rows(
    device: DeviceParameters = DEVICE,
    transmon_reference: TransmonCoherenceReference = TRANSMON_REFERENCE,
) -> list[dict[str, object]]:
    storage_tphi_s = pure_dephasing_time_from_t1_t2(t1_s=device.storage_t1_s, t2_s=device.storage_t2_ramsey_s)
    return [
        {
            "category": "device",
            "parameter": "readout_frequency",
            "value": float(device.readout_frequency_hz),
            "unit": "Hz",
            "source": str(device.parameter_source),
            "note": "Readout-mode bare frequency.",
        },
        {
            "category": "device",
            "parameter": "qubit_frequency",
            "value": float(device.qubit_frequency_hz),
            "unit": "Hz",
            "source": str(device.parameter_source),
            "note": "Transmon bare frequency used in the dispersive model.",
        },
        {
            "category": "device",
            "parameter": "storage_frequency",
            "value": float(device.storage_frequency_hz),
            "unit": "Hz",
            "source": str(device.parameter_source),
            "note": "Storage-mode bare frequency.",
        },
        {
            "category": "device",
            "parameter": "readout_kappa",
            "value": float(device.readout_kappa_hz),
            "unit": "1/s",
            "source": str(device.parameter_source),
            "note": "Readout linewidth used directly in the open-system replay.",
        },
        {
            "category": "device",
            "parameter": "qubit_anharmonicity",
            "value": float(device.qubit_anharmonicity_hz),
            "unit": "Hz",
            "source": str(device.parameter_source),
            "note": "Transmon anharmonicity alpha.",
        },
        {
            "category": "device",
            "parameter": "chi_storage",
            "value": float(device.chi_storage_hz),
            "unit": "Hz",
            "source": str(device.parameter_source),
            "note": "Storage dispersive shift per transmon excitation.",
        },
        {
            "category": "device",
            "parameter": "chi_readout",
            "value": float(device.chi_readout_hz),
            "unit": "Hz",
            "source": str(device.parameter_source),
            "note": "Readout dispersive shift per transmon excitation.",
        },
        {
            "category": "device",
            "parameter": "storage_sideband_nominal_frequency",
            "value": float(device.storage_gf_sideband_nominal_hz),
            "unit": "Hz",
            "source": str(device.parameter_source),
            "note": "Nominal storage gf sideband frequency exported by the local example.",
        },
        {
            "category": "device",
            "parameter": "storage_t1",
            "value": float(device.storage_t1_s),
            "unit": "s",
            "source": str(device.parameter_source),
            "note": "Storage energy-relaxation time used in the available open-system model.",
        },
        {
            "category": "device",
            "parameter": "storage_t2_ramsey",
            "value": float(device.storage_t2_ramsey_s),
            "unit": "s",
            "source": str(device.parameter_source),
            "note": "Storage Ramsey coherence time from the local example.",
        },
        {
            "category": "derived",
            "parameter": "storage_tphi_ramsey",
            "value": None if storage_tphi_s is None else float(storage_tphi_s),
            "unit": "s",
            "source": "derived from storage_t1 and storage_t2_ramsey",
            "note": "Pure dephasing time inferred by cqed_sim convention.",
        },
        {
            "category": "reference",
            "parameter": "transmon_t1_reference",
            "value": None if transmon_reference.qubit_t1_s is None else float(transmon_reference.qubit_t1_s),
            "unit": "s",
            "source": str(transmon_reference.parameter_source),
            "note": "Matched local transmon coherence reference used for sensitivity analysis.",
        },
        {
            "category": "reference",
            "parameter": "transmon_t2_ramsey_reference",
            "value": None if transmon_reference.qubit_t2_ramsey_s is None else float(transmon_reference.qubit_t2_ramsey_s),
            "unit": "s",
            "source": str(transmon_reference.parameter_source),
            "note": "Matched local Ramsey coherence reference for the transmon.",
        },
        {
            "category": "derived",
            "parameter": "transmon_tphi_ramsey_reference",
            "value": None if transmon_reference.qubit_tphi_ramsey_s is None else float(transmon_reference.qubit_tphi_ramsey_s),
            "unit": "s",
            "source": "derived from transmon_t1_reference and transmon_t2_ramsey_reference",
            "note": "Pure dephasing time used when transmon dephasing is enabled.",
        },
    ]


def state_population(state: qt.Qobj, basis_state: qt.Qobj) -> float:
    projector = basis_state.proj()
    if state.isket:
        return float(abs(basis_state.overlap(state)) ** 2)
    return float(np.real((projector * state).tr()))


def sorted_basis_populations(
    state: qt.Qobj,
    model: DispersiveReadoutTransmonStorageModel,
    *,
    cutoff: int = 6,
) -> list[list[float | str]]:
    entries: list[list[float | str]] = []
    for q_level in range(model.n_tr):
        for storage_level in range(model.n_storage):
            for readout_level in range(model.n_readout):
                basis = model.basis_state(q_level, storage_level, readout_level)
                pop = state_population(state, basis)
                if pop > 1.0e-8:
                    entries.append([basis_label(q_level, storage_level, readout_level), float(pop)])
    entries.sort(key=lambda row: float(row[1]), reverse=True)
    return entries[:cutoff]


def pulse_timeseries(
    variant: FamilyVariant,
    duration_s: float,
    amplitude_hz: float,
    *,
    n_points: int = 4097,
) -> tuple[np.ndarray, np.ndarray]:
    grid = np.linspace(0.0, 1.0, int(n_points))
    envelope = build_envelope(variant)(grid)
    times_ns = grid * float(duration_s) * 1.0e9
    amplitude_mhz = envelope * float(amplitude_hz) / 1.0e6
    return times_ns, amplitude_mhz


def spectrum_magnitude(
    variant: FamilyVariant,
    duration_s: float,
    *,
    n_points: int = 4096,
) -> tuple[np.ndarray, np.ndarray]:
    times_ns, amplitude_mhz = pulse_timeseries(variant, duration_s, amplitude_hz=1.0e6, n_points=n_points + 1)
    dt_s = float(duration_s) / float(n_points)
    spectrum = np.fft.fftshift(np.fft.fft(np.asarray(amplitude_mhz, dtype=np.complex128)))
    freqs_hz = np.fft.fftshift(np.fft.fftfreq(amplitude_mhz.size, d=dt_s))
    magnitude = np.abs(spectrum)
    magnitude = magnitude / max(np.max(magnitude), 1.0e-15)
    return freqs_hz / 1.0e6, magnitude


def make_pulse(
    *,
    channel: str,
    carrier_rad_s: float,
    duration_s: float,
    amplitude_hz: float,
    variant: FamilyVariant,
    t0_s: float = 0.0,
    label: str | None = None,
) -> Pulse:
    return Pulse(
        channel=str(channel),
        t0=float(t0_s),
        duration=float(duration_s),
        envelope=build_envelope(variant),
        amp=float(to_internal_units(amplitude_hz)),
        carrier=float(carrier_for_transition_frequency(carrier_rad_s)),
        phase=0.0,
        drag=0.0,
        label=label,
    )


def compile_single_pulse(pulse: Pulse, *, duration_s: float, dt_s: float) -> object:
    return SequenceCompiler(dt=float(dt_s)).compile([pulse], t_end=float(max(duration_s, dt_s)))


def simulate_single_pulse(
    model: DispersiveReadoutTransmonStorageModel,
    initial_state: qt.Qobj,
    *,
    pulse: Pulse,
    duration_s: float,
    drive_target: SidebandDriveSpec,
    frame: FrameSpec,
    dt_s: float,
    noise: NoiseSpec | None = None,
    store_states: bool = False,
) -> tuple[object, object]:
    compiled = compile_single_pulse(pulse, duration_s=duration_s, dt_s=dt_s)
    result = simulate_sequence(
        model,
        compiled,
        initial_state,
        {str(pulse.channel): drive_target},
        SimulationConfig(frame=frame, store_states=store_states, max_step=float(dt_s)),
        noise=noise,
    )
    return compiled, result


def sideband_drive_target(mode: str) -> SidebandDriveSpec:
    return SidebandDriveSpec(mode=mode, lower_level=0, upper_level=2, sideband="red")


def wrap_phase(angle_rad: float) -> float:
    return float(np.angle(np.exp(1j * float(angle_rad))))


def projected_swap_metrics(
    model: DispersiveReadoutTransmonStorageModel,
    frame: FrameSpec,
    *,
    mode: str,
    n: int,
    pulse: Pulse,
    duration_s: float,
    dt_s: float,
) -> dict[str, float]:
    source, target = basis_state_for_mode(model, mode, n)
    drive_target = sideband_drive_target(mode)
    _, result_source = simulate_single_pulse(
        model,
        source,
        pulse=pulse,
        duration_s=duration_s,
        drive_target=drive_target,
        frame=frame,
        dt_s=dt_s,
        store_states=False,
    )
    _, result_target = simulate_single_pulse(
        model,
        target,
        pulse=pulse,
        duration_s=duration_s,
        drive_target=drive_target,
        frame=frame,
        dt_s=dt_s,
        store_states=False,
    )

    psi_source = result_source.final_state
    psi_target = result_target.final_state
    u_proj = np.array(
        [
            [source.overlap(psi_source), source.overlap(psi_target)],
            [target.overlap(psi_source), target.overlap(psi_target)],
        ],
        dtype=np.complex128,
    )
    target_swap = np.array([[0.0, -1.0j], [-1.0j, 0.0]], dtype=np.complex128)
    projected_fidelity = float(abs(np.trace(target_swap.conj().T @ u_proj)) ** 2 / 4.0)
    offdiag_a = u_proj[0, 1]
    offdiag_b = u_proj[1, 0]
    phase_asymmetry = wrap_phase(np.angle(offdiag_b) - np.angle(offdiag_a))
    return {
        "projected_swap_fidelity": projected_fidelity,
        "projected_diagonal_norm": float(np.linalg.norm(np.diag(u_proj))),
        "offdiag_phase_asymmetry_rad": phase_asymmetry,
        "source_column_subspace_norm": float(np.sum(np.abs(u_proj[:, 0]) ** 2)),
        "target_column_subspace_norm": float(np.sum(np.abs(u_proj[:, 1]) ** 2)),
    }


def candidate_target_transfer(
    final_state: qt.Qobj,
    source: qt.Qobj,
    target: qt.Qobj,
) -> dict[str, float]:
    p_target = state_population(final_state, target)
    p_source = state_population(final_state, source)
    leakage = max(0.0, 1.0 - p_target - p_source)
    return {
        "target_probability": float(p_target),
        "source_probability": float(p_source),
        "leakage_probability": float(leakage),
        "subspace_population": float(p_target + p_source),
    }


def rank_cases(cases: list[dict[str, object]]) -> list[dict[str, object]]:
    return sorted(
        cases,
        key=lambda row: (
            -float(row["target_probability"]),
            float(row["leakage_probability"]),
            float(row["duration_ns"]),
            abs(float(row["amplitude_scale"]) - 1.0),
        ),
    )


def sanitize_for_json(value: object) -> object:
    if isinstance(value, dict):
        return {str(key): sanitize_for_json(inner) for key, inner in value.items()}
    if isinstance(value, (list, tuple)):
        return [sanitize_for_json(inner) for inner in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, complex):
        return {"real": float(np.real(value)), "imag": float(np.imag(value))}
    return value


def json_dump(path: Path, payload: object) -> None:
    path.write_text(json.dumps(sanitize_for_json(payload), indent=2, sort_keys=True), encoding="utf-8")


def csv_dump(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    header = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=header)
        writer.writeheader()
        for row in rows:
            serializable = {}
            for key in header:
                value = row.get(key)
                if isinstance(value, (dict, list, tuple)):
                    serializable[key] = json.dumps(sanitize_for_json(value), sort_keys=True)
                else:
                    serializable[key] = sanitize_for_json(value)
            writer.writerow(serializable)


def plot_save(fig: plt.Figure, stem: str) -> None:
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def write_device_manifest(path: Path) -> None:
    import cqed_sim

    json_dump(
        path,
        {
            "device_parameters": asdict(DEVICE),
            "transmon_reference": asdict(TRANSMON_REFERENCE),
            "device_parameter_rows": device_parameter_rows(),
            "cqed_sim_file": str(inspect.getfile(cqed_sim)),
            "python_version": sys.version,
            "platform": platform.platform(),
            "internal_units": {
                "frequencies": "rad/s",
                "times": "s",
                "noise_rates": "1/s",
            },
            "tensor_ordering": ["transmon", "storage", "readout"],
            "reported_state_ordering": ["transmon", "readout", "storage"],
        },
    )
