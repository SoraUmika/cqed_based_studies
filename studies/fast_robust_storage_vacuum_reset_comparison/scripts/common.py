"""Shared helpers for the fast robust storage vacuum-reset comparison study."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import csv
import importlib.util
import inspect
import json
import math
from pathlib import Path
import sys
from typing import Callable, Iterable

try:
    from .runtime_compat import patch_windows_qutip_import
except ImportError:
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

REPO_ROOT = STUDY_DIR.parent.parent
TWO_PI = 2.0 * np.pi

DEFAULT_SWEEP_DT_S = 1.0e-9
DEFAULT_TRAJECTORY_DT_S = 0.5e-9
DEFAULT_CONTINUOUS_DURATION_S = 1.5e-6
DEFAULT_PULSED_RINGDOWN_MULTIPLES = (1.5, 2.0, 3.0, 4.0)

PULSED_RECOMMENDATION_PATH = REPO_ROOT / "studies" / "storage_active_cooling_gf_sideband" / "data" / "recommendation_table.csv"
WAVEFORM_TWO_TONE_CASES_PATH = REPO_ROOT / "studies" / "gf_sideband_waveform_optimization" / "data" / "two_tone_selected_cases.csv"

Q_LABELS = {0: "g", 1: "e", 2: "f", 3: "h", 4: "i", 5: "j"}


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


@dataclass(frozen=True)
class SchemeSummary:
    scheme_key: str
    scheme_label: str
    class_label: str
    auxiliary_only: bool
    mechanism: str
    uses_real_transmon_population: bool
    transmon_virtuality: str
    control_complexity: str
    experimental_realism: str
    recommended_duration_ns: float
    baseline_protocol_duration_ns: float
    baseline_storage_n_final: float
    baseline_transmon_excited_final: float
    baseline_ground_vacuum_final: float
    baseline_storage_n_steady: float
    baseline_transmon_excited_steady: float
    baseline_ground_vacuum_steady: float
    baseline_time_to_threshold_ns: float | None
    baseline_e_fold_time_ns: float | None
    max_hplus_population: float
    coherence_pass_fraction: float
    calibration_pass_fraction: float
    robustness_score: float
    notes: str


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


DEVICE = load_device_from_local_example()
TRANSMON_REFERENCE = load_transmon_reference_from_local_example()


EnvelopeFunc = Callable[[np.ndarray], np.ndarray]


def _normalize(values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    area = np.trapezoid(np.abs(values), grid) if hasattr(np, "trapezoid") else np.trapz(np.abs(values), grid)
    if abs(float(area)) < 1.0e-15:
        return np.asarray(values, dtype=np.complex128)
    return np.asarray(values / float(area), dtype=np.complex128)


def square_envelope() -> EnvelopeFunc:
    def envelope(t_rel: np.ndarray) -> np.ndarray:
        return np.ones_like(np.asarray(t_rel, dtype=float), dtype=np.complex128)

    return envelope


def gaussian_envelope(sigma_fraction: float = 0.18) -> EnvelopeFunc:
    sigma = max(0.02, float(sigma_fraction))
    grid = np.linspace(0.0, 1.0, 4001)
    base = np.exp(-0.5 * ((grid - 0.5) / sigma) ** 2)
    normalized = _normalize(base, grid)

    def envelope(t_rel: np.ndarray) -> np.ndarray:
        t = np.asarray(t_rel, dtype=float)
        return np.interp(t, grid, np.real(normalized)).astype(np.complex128)

    return envelope


def cosine_squared_envelope() -> EnvelopeFunc:
    grid = np.linspace(0.0, 1.0, 4001)
    base = np.cos(np.pi * (grid - 0.5)) ** 2
    normalized = _normalize(base, grid)

    def envelope(t_rel: np.ndarray) -> np.ndarray:
        t = np.asarray(t_rel, dtype=float)
        return np.interp(t, grid, np.real(normalized)).astype(np.complex128)

    return envelope


def bump_envelope() -> EnvelopeFunc:
    grid = np.linspace(0.0, 1.0, 4001)
    safe = np.clip(grid, 1.0e-4, 1.0 - 1.0e-4)
    base = np.exp(-1.0 / (safe * (1.0 - safe)))
    base[(grid <= 0.0) | (grid >= 1.0)] = 0.0
    normalized = _normalize(base, grid)

    def envelope(t_rel: np.ndarray) -> np.ndarray:
        t = np.asarray(t_rel, dtype=float)
        return np.interp(t, grid, np.real(normalized)).astype(np.complex128)

    return envelope


ENVELOPE_BUILDERS: dict[str, Callable[[], EnvelopeFunc]] = {
    "square": square_envelope,
    "gaussian": lambda: gaussian_envelope(0.18),
    "cosine_squared": cosine_squared_envelope,
    "bump": bump_envelope,
}


def build_model(
    *,
    device: DeviceParameters = DEVICE,
    n_storage: int = 6,
    n_readout: int = 4,
    n_tr: int = 5,
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
    *,
    device: DeviceParameters = DEVICE,
    transmon_t1_s: float | None = None,
    transmon_t2_ramsey_s: float | None = None,
    transmon_tphi_s: float | None = None,
    readout_kappa_scale: float = 1.0,
    nth_storage: float = 0.0,
    nth_readout: float = 0.0,
) -> NoiseSpec:
    storage_tphi_s = pure_dephasing_time_from_t1_t2(t1_s=device.storage_t1_s, t2_s=device.storage_t2_ramsey_s)
    if transmon_tphi_s is None and transmon_t2_ramsey_s is not None:
        transmon_tphi_s = pure_dephasing_time_from_t1_t2(t1_s=transmon_t1_s, t2_s=transmon_t2_ramsey_s)
    return NoiseSpec(
        t1=transmon_t1_s,
        tphi=transmon_tphi_s,
        kappa_storage=None if device.storage_t1_s <= 0.0 else 1.0 / float(device.storage_t1_s),
        tphi_storage=storage_tphi_s,
        kappa_readout=float(device.readout_kappa_hz) * float(readout_kappa_scale),
        nth_storage=float(nth_storage),
        nth_readout=float(nth_readout),
    )


def baseline_noise(
    *,
    t1_scale: float = 1.0,
    t2_scale: float = 1.0,
    readout_kappa_scale: float = 1.0,
    nth_storage: float = 0.0,
    nth_readout: float = 0.0,
) -> NoiseSpec:
    t1 = None if TRANSMON_REFERENCE.qubit_t1_s is None else float(TRANSMON_REFERENCE.qubit_t1_s) * float(t1_scale)
    t2 = None if TRANSMON_REFERENCE.qubit_t2_ramsey_s is None else float(TRANSMON_REFERENCE.qubit_t2_ramsey_s) * float(t2_scale)
    tphi = None if t1 is None or t2 is None else pure_dephasing_time_from_t1_t2(t1_s=t1, t2_s=t2)
    return build_noise(
        transmon_t1_s=t1,
        transmon_t2_ramsey_s=t2,
        transmon_tphi_s=tphi,
        readout_kappa_scale=readout_kappa_scale,
        nth_storage=nth_storage,
        nth_readout=nth_readout,
    )


def basis_label(q_level: int, storage_level: int, readout_level: int) -> str:
    q = Q_LABELS.get(int(q_level), f"q{int(q_level)}")
    return f"|{q},{int(readout_level)},{int(storage_level)}>"


def pulsed_recommendations() -> dict[int, dict[str, object]]:
    rows: dict[int, dict[str, object]] = {}
    with PULSED_RECOMMENDATION_PATH.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows[int(row["n"])] = {
                "n": int(row["n"]),
                "step_a_frequency_GHz": float(row["step_a_frequency_GHz"]),
                "step_b_frequency_GHz": float(row["step_b_frequency_GHz"]),
                "step_a_family": str(row["step_a_family"]),
                "step_b_family": str(row["step_b_family"]),
                "step_a_duration_ns": float(row["step_a_duration_ns"]),
                "step_b_duration_ns": float(row["step_b_duration_ns"]),
                "step_a_amplitude_MHz": float(row["step_a_amplitude_MHz"]),
                "step_b_amplitude_MHz": float(row["step_b_amplitude_MHz"]),
                "single_cycle_success": float(row["single_cycle_success"]),
            }
    return rows


def waveform_selected_two_tone_cases() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with WAVEFORM_TWO_TONE_CASES_PATH.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            rows.append(
                {
                    "n": int(row["n"]),
                    "mechanism": str(row["mechanism"]),
                    "target_coupling_MHz": float(row["target_coupling_MHz"]),
                    "common_detuning_MHz": float(row["common_detuning_MHz"]),
                    "storage_amplitude_MHz": float(row["storage_amplitude_MHz"]),
                    "readout_amplitude_MHz": float(row["readout_amplitude_MHz"]),
                    "case_role": str(row["case_role"]),
                    "selection_status": str(row["selection_status"]),
                }
            )
    return rows


def storage_sideband_rotating_frequency(model: DispersiveReadoutTransmonStorageModel, frame: FrameSpec, n: int) -> float:
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


def readout_sideband_rotating_frequency(model: DispersiveReadoutTransmonStorageModel, frame: FrameSpec, n: int) -> float:
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
    _, operator = model.sideband_drive_operators(mode=str(mode), lower_level=0, upper_level=2, sideband="red")
    if str(mode) == "storage":
        source = model.basis_state(0, int(n), 0)
        intermediate = model.basis_state(2, int(n) - 1, 0)
        return float(abs(source.overlap(operator * intermediate)))
    elif str(mode) == "readout":
        source = model.basis_state(2, int(n) - 1, 0)
        target = model.basis_state(0, int(n) - 1, 1)
        return float(abs(target.overlap(operator * source)))
    else:
        raise ValueError(f"Unsupported mode '{mode}'.")


def make_pulse(
    *,
    channel: str,
    carrier_rad_s: float,
    duration_s: float,
    amplitude_hz: float,
    envelope_name: str,
    t0_s: float = 0.0,
    label: str | None = None,
) -> Pulse:
    if envelope_name not in ENVELOPE_BUILDERS:
        raise ValueError(f"Unsupported envelope family '{envelope_name}'.")
    return Pulse(
        channel=str(channel),
        t0=float(t0_s),
        duration=float(duration_s),
        envelope=ENVELOPE_BUILDERS[envelope_name](),
        amp=float(to_internal_units(amplitude_hz)),
        carrier=float(carrier_for_transition_frequency(carrier_rad_s)),
        phase=0.0,
        drag=0.0,
        label=label,
    )


def run_sequence(
    model: DispersiveReadoutTransmonStorageModel,
    initial_state: qt.Qobj,
    *,
    pulses: list[Pulse],
    duration_s: float,
    drive_ops: dict[str, object],
    frame: FrameSpec,
    noise: NoiseSpec | None,
    dt_s: float,
    e_ops: dict[str, qt.Qobj] | None = None,
    store_states: bool = False,
):
    compiled = SequenceCompiler(dt=float(dt_s)).compile(pulses, t_end=float(max(duration_s, dt_s)))
    result = simulate_sequence(
        model,
        compiled,
        initial_state,
        drive_ops,
        SimulationConfig(frame=frame, store_states=store_states, max_step=float(dt_s)),
        noise=noise,
        e_ops=e_ops,
    )
    return compiled, result


def tensor_product(*factors: qt.Qobj) -> qt.Qobj:
    return qt.tensor(*factors)


def _as_dm(obj: qt.Qobj) -> qt.Qobj:
    return obj if obj.isoper else obj.proj()


def fock_state(model: DispersiveReadoutTransmonStorageModel, *, q_level: int = 0, storage_level: int = 0, readout_level: int = 0) -> qt.Qobj:
    dims = tuple(int(dim) for dim in model.subsystem_dims)
    return qt.tensor(
        qt.basis(dims[0], int(q_level)),
        qt.basis(dims[1], int(storage_level)),
        qt.basis(dims[2], int(readout_level)),
    )


def coherent_storage_state(
    model: DispersiveReadoutTransmonStorageModel,
    *,
    alpha: complex,
    q_level: int = 0,
    readout_level: int = 0,
) -> qt.Qobj:
    dims = tuple(int(dim) for dim in model.subsystem_dims)
    return qt.tensor(
        qt.basis(dims[0], int(q_level)),
        qt.coherent(dims[1], complex(alpha)),
        qt.basis(dims[2], int(readout_level)),
    )


def thermal_storage_state(
    model: DispersiveReadoutTransmonStorageModel,
    *,
    nbar: float,
    excited_reset_prob: float = 0.0,
) -> qt.Qobj:
    dims = tuple(int(dim) for dim in model.subsystem_dims)
    rho_q = (1.0 - float(excited_reset_prob)) * qt.basis(dims[0], 0).proj() + float(excited_reset_prob) * qt.basis(dims[0], 1).proj()
    rho_s = qt.thermal_dm(dims[1], float(nbar))
    rho_r = qt.basis(dims[2], 0).proj()
    return qt.tensor(rho_q, rho_s, rho_r)


def imperfect_reset_state(
    model: DispersiveReadoutTransmonStorageModel,
    *,
    storage_state: qt.Qobj,
    excited_reset_prob: float,
) -> qt.Qobj:
    dims = tuple(int(dim) for dim in model.subsystem_dims)
    rho_q = (1.0 - float(excited_reset_prob)) * qt.basis(dims[0], 0).proj() + float(excited_reset_prob) * qt.basis(dims[0], 1).proj()
    storage_dm = _as_dm(storage_state)
    if len(storage_dm.dims[0]) == 3:
        storage_factor = storage_dm.ptrace(1)
        readout_factor = storage_dm.ptrace(2)
    elif len(storage_dm.dims[0]) == 1:
        storage_factor = storage_dm
        readout_factor = qt.basis(dims[2], 0).proj()
    else:
        raise ValueError("Unsupported storage_state dimensions for imperfect reset construction.")
    return qt.tensor(rho_q, storage_factor, readout_factor)


def expectation_operators(model: DispersiveReadoutTransmonStorageModel) -> dict[str, qt.Qobj]:
    dims = tuple(int(dim) for dim in model.subsystem_dims)
    ops = model.operators()
    e_ops: dict[str, qt.Qobj] = {
        "n_storage": ops["n_s"],
        "n_readout": ops["n_r"],
        "P_g": qt.tensor(qt.basis(dims[0], 0).proj(), qt.qeye(dims[1]), qt.qeye(dims[2])),
        "P_e": qt.tensor(qt.basis(dims[0], 1).proj(), qt.qeye(dims[1]), qt.qeye(dims[2])),
        "P_f": qt.tensor(qt.basis(dims[0], 2).proj(), qt.qeye(dims[1]), qt.qeye(dims[2])),
        "P_g00": model.basis_state(0, 0, 0).proj(),
    }
    if dims[0] >= 4:
        e_ops["P_hplus"] = sum(
            qt.tensor(qt.basis(dims[0], level).proj(), qt.qeye(dims[1]), qt.qeye(dims[2]))
            for level in range(3, dims[0])
        )
    else:
        e_ops["P_hplus"] = 0.0 * e_ops["P_g"]
    return e_ops


def state_population(state: qt.Qobj, basis_state: qt.Qobj) -> float:
    projector = basis_state if basis_state.isoper else basis_state.proj()
    return float(np.real((projector * _as_dm(state)).tr()))


def storage_mean_photon_number(state: qt.Qobj, model: DispersiveReadoutTransmonStorageModel) -> float:
    return float(np.real((model.operators()["n_s"] * _as_dm(state)).tr()))


def readout_mean_photon_number(state: qt.Qobj, model: DispersiveReadoutTransmonStorageModel) -> float:
    return float(np.real((model.operators()["n_r"] * _as_dm(state)).tr()))


def transmon_excited_population(state: qt.Qobj, model: DispersiveReadoutTransmonStorageModel) -> float:
    dims = tuple(int(dim) for dim in model.subsystem_dims)
    rho = _as_dm(state).ptrace(0)
    return float(sum(np.real(rho[level, level]) for level in range(1, dims[0])))


def hplus_population(state: qt.Qobj, model: DispersiveReadoutTransmonStorageModel) -> float:
    dims = tuple(int(dim) for dim in model.subsystem_dims)
    if dims[0] <= 3:
        return 0.0
    rho = _as_dm(state).ptrace(0)
    return float(sum(np.real(rho[level, level]) for level in range(3, dims[0])))


def ground_vacuum_population(state: qt.Qobj, model: DispersiveReadoutTransmonStorageModel) -> float:
    return state_population(state, model.basis_state(0, 0, 0))


def linear_entropy(state: qt.Qobj) -> float:
    rho = _as_dm(state)
    return float(max(0.0, 1.0 - np.real((rho * rho).tr())))


def first_threshold_time_ns(times_s: np.ndarray, storage_curve: np.ndarray, transmon_curve: np.ndarray, *, storage_threshold: float = 0.01, transmon_threshold: float = 1.0e-3) -> float | None:
    mask = np.logical_and(np.asarray(storage_curve) <= float(storage_threshold), np.asarray(transmon_curve) <= float(transmon_threshold))
    if not np.any(mask):
        return None
    index = int(np.argmax(mask))
    return float(np.asarray(times_s)[index] * 1.0e9)


def e_fold_time_ns(times_s: np.ndarray, curve: np.ndarray) -> float | None:
    values = np.asarray(curve, dtype=float)
    if values.size == 0:
        return None
    initial = float(values[0])
    if initial <= 0.0:
        return None
    target = initial / math.e
    mask = values <= target
    if not np.any(mask):
        return None
    index = int(np.argmax(mask))
    return float(np.asarray(times_s, dtype=float)[index] * 1.0e9)


def estimate_decay_rate_hz(times_s: np.ndarray, storage_curve: np.ndarray) -> float | None:
    t = np.asarray(times_s, dtype=float)
    y = np.asarray(storage_curve, dtype=float)
    mask = np.logical_and(y > 0.05 * max(float(y[0]), 1.0e-12), y < 0.95 * max(float(y[0]), 1.0e-12))
    if np.count_nonzero(mask) < 4:
        return None
    coeffs = np.polyfit(t[mask], np.log(y[mask]), 1)
    rate = -float(coeffs[0])
    if not np.isfinite(rate) or rate <= 0.0:
        return None
    return rate


def sanitize_for_json(value: object) -> object:
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return [sanitize_for_json(item) for item in value.tolist()]
    if isinstance(value, dict):
        return {str(key): sanitize_for_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [sanitize_for_json(item) for item in value]
    if hasattr(value, "__dataclass_fields__"):
        return sanitize_for_json(asdict(value))
    return str(value)


def json_dump(path: Path, payload: object) -> None:
    path.write_text(json.dumps(sanitize_for_json(payload), indent=2), encoding="utf-8")


def csv_dump(path: Path, rows: Iterable[dict[str, object]]) -> None:
    rows = list(rows)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row.keys():
            if key not in seen:
                seen.add(key)
                fieldnames.append(str(key))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: sanitize_for_json(row.get(key)) for key in fieldnames})


def plot_save(fig: plt.Figure, stem: str) -> None:
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / f"{stem}.pdf", bbox_inches="tight")


def write_device_manifest(path: Path) -> None:
    payload = {
        "device": asdict(DEVICE),
        "transmon_reference": asdict(TRANSMON_REFERENCE),
        "pulsed_recommendation_source": str(PULSED_RECOMMENDATION_PATH),
        "two_tone_case_source": str(WAVEFORM_TWO_TONE_CASES_PATH),
    }
    json_dump(path, payload)
