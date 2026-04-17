"""Common helpers for the Fock-resolved black-box SQR inference study (v2)."""

from __future__ import annotations

import json
import os
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import matplotlib

matplotlib.use("Agg")

import numpy as np
import qutip as qt

import runtime_compat  # noqa: F401

from cqed_sim.calibration.conditioned_multitone import (
    ConditionedMultitoneCorrections,
    ConditionedMultitoneRunConfig,
    ConditionedQubitTargets,
    build_conditioned_multitone_tones,
    build_conditioned_multitone_waveform,
)
from cqed_sim.core import DispersiveTransmonCavityModel, FrameSpec, displacement_op
from cqed_sim.core.ideal_gates import qubit_rotation_axis, qubit_rotation_xy, sqr_op
from cqed_sim.sequence import SequenceCompiler
from cqed_sim.sim import SimulationConfig, prepare_simulation
from cqed_sim.sim.extractors import (
    conditioned_population,
    conditioned_qubit_state,
    truncate_to_qubit_subspace,
)


STUDY_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = STUDY_DIR / "data"
FIGURES_DIR = STUDY_DIR / "figures"
ARTIFACTS_DIR = STUDY_DIR / "artifacts"
REPORT_DIR = STUDY_DIR / "report"

for _path in (DATA_DIR, FIGURES_DIR, ARTIFACTS_DIR, REPORT_DIR):
    _path.mkdir(parents=True, exist_ok=True)


STYLE_PATH = (
    STUDY_DIR.parent.parent
    / ".github"
    / "skills"
    / "publication-figures"
    / "assets"
    / "cqed_style.mplstyle"
)

TWO_PI = 2.0 * np.pi

OMEGA_Q = TWO_PI * 6.150e9
OMEGA_C = TWO_PI * 5.241e9
ALPHA = TWO_PI * (-255.0e6)
CHI = TWO_PI * (-2.84e6)
CHI_PRIME = TWO_PI * (-21.0e3)
KERR = TWO_PI * (-28.0e3)

# Default cavity sector populations and qubit rotation angles for analytic cases.
DEFAULT_POPULATIONS = np.asarray([0.41, 0.27, 0.19, 0.13], dtype=float)
DEFAULT_THETA = np.asarray([np.pi / 3.0, np.pi / 2.0, 0.9 * np.pi, np.pi / 4.0], dtype=float)
DEFAULT_PHI = np.asarray([0.0, np.pi / 3.0, np.pi / 2.0, -np.pi / 4.0], dtype=float)
DEFAULT_CPSQR_PHASES = np.asarray([0.0, 0.35 * np.pi, -0.45 * np.pi, 0.7 * np.pi], dtype=float)

N_ACTIVE = 4
SHORT_DURATION_S = 0.35e-6   # intentionally short / under-optimized
LONG_DURATION_S = 1.0e-6     # near-ideal optimized duration
DEFAULT_DT_S = 4.0e-9
DEFAULT_SIGMA_FRACTION = 1.0 / 6.0

PAULI_X = qt.sigmax()
PAULI_Y = qt.sigmay()
PAULI_Z = qt.sigmaz()
PAULI_I = qt.qeye(2)


@dataclass(frozen=True)
class SectorSummary:
    level: int
    population: float
    leakage: float
    x: float
    y: float
    z: float

    @property
    def bloch_vector(self) -> np.ndarray:
        return np.asarray([self.x, self.y, self.z], dtype=float)

    def as_dict(self) -> dict[str, float | int]:
        return {
            "level": int(self.level),
            "population": float(self.population),
            "leakage": float(self.leakage),
            "x": float(self.x),
            "y": float(self.y),
            "z": float(self.z),
        }


def apply_plot_style() -> None:
    import matplotlib.pyplot as plt

    if STYLE_PATH.exists():
        plt.style.use(str(STYLE_PATH))


def json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        if np.iscomplexobj(value):
            return {
                "real": np.real(value).tolist(),
                "imag": np.imag(value).tolist(),
                "shape": list(value.shape),
            }
        return value.tolist()
    if isinstance(value, complex):
        return {"real": float(value.real), "imag": float(value.imag)}
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def save_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(json_ready(payload), indent=2, sort_keys=True)
    temp_path = path.with_suffix(path.suffix + ".tmp")
    last_error: OSError | None = None
    for _ in range(5):
        try:
            temp_path.write_text(text, encoding="utf-8")
            os.replace(str(temp_path), str(path))
            return
        except OSError as exc:
            last_error = exc
            time.sleep(0.25)
    if last_error is not None:
        raise last_error


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def save_figure(fig: Any, stem: str) -> dict[str, str]:
    png_path = FIGURES_DIR / f"{stem}.png"
    pdf_path = FIGURES_DIR / f"{stem}.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    return {"png": str(png_path), "pdf": str(pdf_path)}


def normalize_probabilities(values: Sequence[float]) -> np.ndarray:
    arr = np.asarray(values, dtype=float).reshape(-1)
    total = float(np.sum(arr))
    if total <= 0.0:
        raise ValueError("Probability weights must sum to a positive value.")
    return arr / total


def build_model(*, n_cav: int = N_ACTIVE, n_tr: int = 2) -> DispersiveTransmonCavityModel:
    return DispersiveTransmonCavityModel(
        omega_q=OMEGA_Q,
        omega_c=OMEGA_C,
        alpha=ALPHA,
        chi=CHI,
        chi_higher=(CHI_PRIME,),
        kerr=KERR,
        n_cav=int(n_cav),
        n_tr=int(n_tr),
    )


def build_frame(model: DispersiveTransmonCavityModel) -> FrameSpec:
    return FrameSpec(omega_q_frame=float(model.omega_q), omega_c_frame=float(model.omega_c))


def qubit_ground_dm() -> qt.Qobj:
    return qt.basis(2, 0).proj()


def qubit_excited_dm() -> qt.Qobj:
    return qt.basis(2, 1).proj()


def qubit_plus_dm() -> qt.Qobj:
    ket = (qt.basis(2, 0) + qt.basis(2, 1)).unit()
    return ket.proj()


def qubit_plus_y_dm() -> qt.Qobj:
    ket = (qt.basis(2, 0) + 1j * qt.basis(2, 1)).unit()
    return ket.proj()


def qubit_density_from_bloch(x: float, y: float, z: float) -> qt.Qobj:
    return 0.5 * (PAULI_I + float(x) * PAULI_X + float(y) * PAULI_Y + float(z) * PAULI_Z)


def bloch_from_density(rho_q: qt.Qobj) -> tuple[float, float, float]:
    return (
        float(np.real((rho_q * PAULI_X).tr())),
        float(np.real((rho_q * PAULI_Y).tr())),
        float(np.real((rho_q * PAULI_Z).tr())),
    )


def density_matrix_fidelity(rho_a: qt.Qobj, rho_b: qt.Qobj) -> float:
    """Bures fidelity F = (Tr sqrt(sqrt(rho_a) rho_b sqrt(rho_a)))^2."""
    return float(qt.fidelity(rho_a, rho_b) ** 2)


def trace_distance(rho_a: qt.Qobj, rho_b: qt.Qobj) -> float:
    return float(qt.tracedist(rho_a, rho_b))


def displacement_probability_matrix(n_cav: int, alpha: complex) -> np.ndarray:
    dmat = np.asarray(displacement_op(int(n_cav), complex(alpha)).full(), dtype=np.complex128)
    return np.abs(dmat) ** 2


def qubit_block_rotation(theta: float, phi: float) -> qt.Qobj:
    return qubit_rotation_xy(float(theta), float(phi))


def qubit_block_phase(phase: float) -> qt.Qobj:
    return qubit_rotation_axis(float(phase), "z")


def ideal_sqr_operator(theta_values: Sequence[float], phi_values: Sequence[float]) -> qt.Qobj:
    return sqr_op(np.asarray(theta_values, dtype=float), np.asarray(phi_values, dtype=float))


def cpsqr_like_operator(phases: Sequence[float]) -> qt.Qobj:
    """Fock-conditioned qubit phase gate: |n><n| ⊗ Rz(phi_n)."""
    phase_array = np.asarray(phases, dtype=float).reshape(-1)
    n_cav = int(phase_array.size)
    operator = 0 * qt.tensor(qt.qeye(2), qt.qeye(n_cav))
    for level, phase in enumerate(phase_array):
        projector = qt.basis(n_cav, level) * qt.basis(n_cav, level).dag()
        operator += qt.tensor(qubit_block_phase(float(phase)), projector)
    return operator


def cavity_mixture_state(probabilities: Sequence[float]) -> qt.Qobj:
    probs = normalize_probabilities(probabilities)
    dim = int(probs.size)
    rho = 0 * qt.qeye(dim)
    for level, prob in enumerate(probs):
        rho += float(prob) * qt.basis(dim, level).proj()
    return rho


def cavity_superposition_state(coefficients: Sequence[complex]) -> qt.Qobj:
    coeffs = np.asarray(coefficients, dtype=np.complex128).reshape(-1)
    norm = float(np.linalg.norm(coeffs))
    if norm <= 0.0:
        raise ValueError("Superposition coefficients must not all be zero.")
    coeffs = coeffs / norm
    ket = sum(coeffs[level] * qt.basis(coeffs.size, level) for level in range(coeffs.size))
    return ket


def coherent_state(alpha: complex, *, n_cav: int = N_ACTIVE) -> qt.Qobj:
    return qt.coherent(int(n_cav), complex(alpha))


def as_dm(state: qt.Qobj) -> qt.Qobj:
    return state if state.isoper else state.proj()


def sector_summary_from_state(state: qt.Qobj, level: int) -> SectorSummary:
    rho_q_full, population, valid = conditioned_qubit_state(state, n=int(level), fallback="zero")
    if not valid:
        rho_q = qt.Qobj(np.zeros((2, 2), dtype=np.complex128), dims=[[2], [2]])
        leakage = 0.0
        x = y = z = 0.0
    else:
        rho_q, leakage = truncate_to_qubit_subspace(rho_q_full)
        x, y, z = bloch_from_density(rho_q)
    return SectorSummary(
        level=int(level),
        population=float(population),
        leakage=float(leakage),
        x=float(x),
        y=float(y),
        z=float(z),
    )


def sector_summaries_from_state(state: qt.Qobj, levels: Iterable[int]) -> list[SectorSummary]:
    return [sector_summary_from_state(state, level) for level in levels]


def simulate_compiled_on_states(
    model: DispersiveTransmonCavityModel,
    compiled: Any,
    *,
    frame: FrameSpec,
    drive_ops: dict[str, str],
    initial_states: Sequence[qt.Qobj],
) -> list[qt.Qobj]:
    session = prepare_simulation(
        model,
        compiled,
        drive_ops,
        config=SimulationConfig(frame=frame, store_states=False),
    )
    results = session.run_many(initial_states)
    return [result.final_state for result in results]


def conditioned_multitone_targets(
    theta_values: Sequence[float] = DEFAULT_THETA,
    phi_values: Sequence[float] = DEFAULT_PHI,
) -> ConditionedQubitTargets:
    pairs = list(zip(theta_values, phi_values, strict=True))
    return ConditionedQubitTargets.from_spec(pairs, n_levels=len(pairs))


def conditioned_run_config(
    model: DispersiveTransmonCavityModel,
    *,
    duration_s: float,
    dt_s: float = DEFAULT_DT_S,
    sigma_fraction: float = DEFAULT_SIGMA_FRACTION,
) -> ConditionedMultitoneRunConfig:
    frame = build_frame(model)
    return ConditionedMultitoneRunConfig(
        frame=frame,
        duration_s=float(duration_s),
        dt_s=float(dt_s),
        sigma_fraction=float(sigma_fraction),
        tone_cutoff=1.0e-12,
        include_all_levels=False,
        max_step_s=float(dt_s),
    )


def build_multitone_waveform(
    model: DispersiveTransmonCavityModel,
    *,
    duration_s: float,
    theta_values: Sequence[float] = DEFAULT_THETA,
    phi_values: Sequence[float] = DEFAULT_PHI,
    corrections: ConditionedMultitoneCorrections | None = None,
) -> tuple[Any, ConditionedQubitTargets, ConditionedMultitoneRunConfig]:
    targets = conditioned_multitone_targets(theta_values=theta_values, phi_values=phi_values)
    run_config = conditioned_run_config(model, duration_s=float(duration_s))
    corr = ConditionedMultitoneCorrections.zeros(targets.n_levels) if corrections is None else corrections
    tones = build_conditioned_multitone_tones(model, targets, run_config, corrections=corr)
    waveform = build_conditioned_multitone_waveform(
        tones,
        run_config,
        channel="qubit",
        drive_target="qubit",
        label="study_multitone_waveform",
    )
    return waveform, targets, run_config


def compile_waveform(waveform: Any, run_config: ConditionedMultitoneRunConfig) -> Any:
    compiler = SequenceCompiler(dt=float(run_config.dt_s))
    return compiler.compile([waveform.pulse], t_end=float(run_config.duration_s + run_config.dt_s))
