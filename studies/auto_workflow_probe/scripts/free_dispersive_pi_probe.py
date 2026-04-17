"""Minimal model-level probe of free dispersive pi-phase accumulation.

The script keeps the original closed-form helper confirmation and adds an
independent static-Hamiltonian evolution cross-check within cqed_sim.
"""

from __future__ import annotations

import csv
import json
import math
import time
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from cqed_sim.core import FrameSpec
from cqed_sim.core.model import DispersiveTransmonCavityModel
from cqed_sim.gates.coupled import dispersive_phase


STUDY_NAME = "auto_workflow_probe"
OMEGA_Q_HZ = 6.150e9
OMEGA_C_HZ = 5.241e9
ALPHA_HZ = -255.0e6
CHI_HZ = -2.84e6
KERR_HZ = -28.0e3
QUBIT_DIM = 2
CAVITY_DIMS = (2, 3)
HAMILTONIAN_CROSS_CHECK_DIM = 3
NUM_SAMPLES = 51
TWO_PI = 2.0 * math.pi


def _timestamp_utc() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _find_repo_root() -> Path:
    for candidate in Path(__file__).resolve().parents:
        if (candidate / ".github").exists():
            return candidate
    msg = "Could not locate repository root from study script path."
    raise RuntimeError(msg)


REPO_ROOT = _find_repo_root()
STUDY_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = STUDY_ROOT / "data"
FIGURES_DIR = STUDY_ROOT / "figures"
ARTIFACTS_DIR = STUDY_ROOT / "artifacts"
STYLE_PATH = (
    REPO_ROOT
    / ".github"
    / "skills"
    / "publication-figures"
    / "assets"
    / "cqed_style.mplstyle"
)
DATA_PATH = DATA_DIR / "phase_difference_samples.csv"
ARTIFACT_PATH = ARTIFACTS_DIR / "free_dispersive_pi_probe_summary.json"
FIGURE_STEM = FIGURES_DIR / "phase_difference_vs_idle_time"
CROSS_CHECK_DATA_PATH = DATA_DIR / "phase_difference_hamiltonian_cross_check.csv"
CROSS_CHECK_ARTIFACT_PATH = ARTIFACTS_DIR / "free_dispersive_hamiltonian_cross_check.json"
CROSS_CHECK_FIGURE_STEM = FIGURES_DIR / "phase_difference_hamiltonian_cross_check"

OMEGA_Q_RAD_S = TWO_PI * OMEGA_Q_HZ
OMEGA_C_RAD_S = TWO_PI * OMEGA_C_HZ
ALPHA_RAD_S = TWO_PI * ALPHA_HZ
CHI_RAD_S = TWO_PI * CHI_HZ
KERR_RAD_S = TWO_PI * KERR_HZ
MODEL_FRAME = FrameSpec(omega_q_frame=OMEGA_Q_RAD_S, omega_c_frame=OMEGA_C_RAD_S)


def build_model(cavity_dim: int) -> DispersiveTransmonCavityModel:
    return DispersiveTransmonCavityModel(
        omega_q=OMEGA_Q_RAD_S,
        omega_c=OMEGA_C_RAD_S,
        alpha=ALPHA_RAD_S,
        chi=CHI_RAD_S,
        kerr=KERR_RAD_S,
        n_cav=int(cavity_dim),
        n_tr=QUBIT_DIM,
    )


def wrap_phase(phase_values: np.ndarray) -> np.ndarray:
    return (np.asarray(phase_values, dtype=float) + math.pi) % TWO_PI - math.pi


def relative_phase_from_unitary(*, unitary, excited_zero, excited_one) -> float:
    amplitude_zero = complex(excited_zero.overlap(unitary * excited_zero))
    amplitude_one = complex(excited_one.overlap(unitary * excited_one))
    return float(np.angle(amplitude_one * np.conjugate(amplitude_zero)))


def relative_phase_helper_scan(*, model: DispersiveTransmonCavityModel, times_s: np.ndarray, cavity_dim: int) -> np.ndarray:
    excited_zero = model.basis_state(1, 0)
    excited_one = model.basis_state(1, 1)
    raw_relative_phases: list[float] = []
    for time_s in times_s:
        unitary = dispersive_phase(
            chi=CHI_RAD_S,
            time=float(time_s),
            cavity_dim=int(cavity_dim),
            qubit_dim=QUBIT_DIM,
            convention="n_e",
        )
        raw_relative_phases.append(
            relative_phase_from_unitary(unitary=unitary, excited_zero=excited_zero, excited_one=excited_one)
        )
    return np.unwrap(np.asarray(raw_relative_phases, dtype=float))


def relative_phase_hamiltonian_scan(*, model: DispersiveTransmonCavityModel, times_s: np.ndarray) -> np.ndarray:
    excited_zero = model.basis_state(1, 0)
    excited_one = model.basis_state(1, 1)
    static_hamiltonian = model.static_hamiltonian(frame=MODEL_FRAME)
    raw_relative_phases: list[float] = []
    for time_s in times_s:
        unitary = (-1.0j * float(time_s) * static_hamiltonian).expm()
        raw_relative_phases.append(
            relative_phase_from_unitary(unitary=unitary, excited_zero=excited_zero, excited_one=excited_one)
        )
    return np.unwrap(np.asarray(raw_relative_phases, dtype=float))


def save_phase_samples(
    *,
    times_ns: np.ndarray,
    analytic_phase_rad: np.ndarray,
    numerical_phase_dim2_rad: np.ndarray,
    numerical_phase_dim3_rad: np.ndarray,
) -> None:
    with DATA_PATH.open("w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "time_ns",
                "analytic_phase_rad",
                "numerical_phase_dim2_rad",
                "numerical_phase_dim3_rad",
            ]
        )
        for row in zip(
            times_ns,
            analytic_phase_rad,
            numerical_phase_dim2_rad,
            numerical_phase_dim3_rad,
            strict=True,
        ):
            writer.writerow([f"{float(value):.12g}" for value in row])


def save_cross_check_samples(
    *,
    times_ns: np.ndarray,
    analytic_phase_rad: np.ndarray,
    helper_phase_rad: np.ndarray,
    hamiltonian_phase_rad: np.ndarray,
) -> None:
    helper_error = wrap_phase(helper_phase_rad - analytic_phase_rad)
    hamiltonian_error = wrap_phase(hamiltonian_phase_rad - analytic_phase_rad)
    helper_vs_hamiltonian = wrap_phase(hamiltonian_phase_rad - helper_phase_rad)
    with CROSS_CHECK_DATA_PATH.open("w", newline="", encoding="ascii") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "time_ns",
                "analytic_phase_rad",
                "helper_phase_rad",
                "hamiltonian_phase_rad",
                "helper_wrapped_error_rad",
                "hamiltonian_wrapped_error_rad",
                "helper_vs_hamiltonian_wrapped_difference_rad",
            ]
        )
        for row in zip(
            times_ns,
            analytic_phase_rad,
            helper_phase_rad,
            hamiltonian_phase_rad,
            helper_error,
            hamiltonian_error,
            helper_vs_hamiltonian,
            strict=True,
        ):
            writer.writerow([f"{float(value):.12g}" for value in row])


def save_summary(payload: dict[str, object]) -> None:
    with ARTIFACT_PATH.open("w", encoding="ascii") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")


def save_figure(
    *,
    times_ns: np.ndarray,
    analytic_phase_rad: np.ndarray,
    numerical_phase_dim2_rad: np.ndarray,
    numerical_phase_dim3_rad: np.ndarray,
    t_pi_ns: float,
) -> None:
    plt.style.use(str(STYLE_PATH))
    figure, axis = plt.subplots(figsize=(6.4, 4.2))
    axis.plot(times_ns, analytic_phase_rad, label="Closed-form phase law", linewidth=2.2)
    axis.plot(times_ns, numerical_phase_dim2_rad, label="Framework helper, cavity cutoff 2", linewidth=1.8)
    axis.plot(
        times_ns,
        numerical_phase_dim3_rad,
        label="Framework helper, cavity cutoff 3",
        linewidth=1.5,
        linestyle="--",
    )
    axis.axhline(math.pi, color="0.35", linestyle=":", linewidth=1.0, label="Pi target")
    axis.axvline(t_pi_ns, color="0.5", linestyle="-.", linewidth=1.0, label="Closed-form pi crossing")
    axis.set_xlabel("Idle time (ns)")
    axis.set_ylabel("Relative excited-branch phase (rad)")
    axis.set_xlim(float(times_ns[0]), float(times_ns[-1]))
    axis.legend(frameon=False)
    figure.savefig(f"{FIGURE_STEM}.png", dpi=300, bbox_inches="tight")
    figure.savefig(f"{FIGURE_STEM}.pdf", bbox_inches="tight")
    plt.close(figure)


def save_cross_check_figure(
    *,
    times_ns: np.ndarray,
    analytic_phase_rad: np.ndarray,
    helper_phase_rad: np.ndarray,
    hamiltonian_phase_rad: np.ndarray,
) -> None:
    plt.style.use(str(STYLE_PATH))
    figure, (phase_axis, error_axis) = plt.subplots(
        2,
        1,
        figsize=(6.4, 6.0),
        sharex=True,
        gridspec_kw={"height_ratios": [3.0, 1.3]},
    )
    phase_axis.plot(times_ns, analytic_phase_rad, label="Closed-form phase law", linewidth=2.2)
    phase_axis.plot(times_ns, helper_phase_rad, label="Framework helper", linewidth=1.8)
    phase_axis.plot(
        times_ns,
        hamiltonian_phase_rad,
        label="Static-Hamiltonian evolution",
        linewidth=1.8,
        linestyle="--",
    )
    phase_axis.axhline(math.pi, color="0.35", linestyle=":", linewidth=1.0)
    phase_axis.set_ylabel("Relative excited-branch phase (rad)")
    phase_axis.legend(frameon=False)

    helper_error = np.maximum(np.abs(wrap_phase(helper_phase_rad - analytic_phase_rad)), 1.0e-18)
    hamiltonian_error = np.maximum(np.abs(wrap_phase(hamiltonian_phase_rad - analytic_phase_rad)), 1.0e-18)
    error_axis.plot(times_ns, helper_error, label="Helper residual", linewidth=1.6)
    error_axis.plot(times_ns, hamiltonian_error, label="Hamiltonian residual", linewidth=1.6, linestyle="--")
    error_axis.set_yscale("log")
    error_axis.set_xlabel("Idle time (ns)")
    error_axis.set_ylabel("Absolute wrapped residual (rad)")
    error_axis.legend(frameon=False)

    figure.savefig(f"{CROSS_CHECK_FIGURE_STEM}.png", dpi=300, bbox_inches="tight")
    figure.savefig(f"{CROSS_CHECK_FIGURE_STEM}.pdf", bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    start_time = time.perf_counter()
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    ARTIFACTS_DIR.mkdir(parents=True, exist_ok=True)

    t_pi_s = math.pi / abs(CHI_RAD_S)
    times_s = np.linspace(0.0, 2.0 * t_pi_s, NUM_SAMPLES)
    times_ns = times_s * 1.0e9
    analytic_phase_rad = -CHI_RAD_S * times_s

    model_dim2 = build_model(CAVITY_DIMS[0])
    model_dim3 = build_model(CAVITY_DIMS[1])
    hamiltonian_model = build_model(HAMILTONIAN_CROSS_CHECK_DIM)
    numerical_phase_dim2_rad = relative_phase_helper_scan(
        model=model_dim2,
        times_s=times_s,
        cavity_dim=CAVITY_DIMS[0],
    )
    numerical_phase_dim3_rad = relative_phase_helper_scan(
        model=model_dim3,
        times_s=times_s,
        cavity_dim=CAVITY_DIMS[1],
    )
    hamiltonian_phase_rad = relative_phase_hamiltonian_scan(model=hamiltonian_model, times_s=times_s)

    max_wrapped_error_dim2 = float(np.max(np.abs(wrap_phase(numerical_phase_dim2_rad - analytic_phase_rad))))
    max_wrapped_error_dim3 = float(np.max(np.abs(wrap_phase(numerical_phase_dim3_rad - analytic_phase_rad))))
    max_wrapped_difference = float(
        np.max(np.abs(wrap_phase(numerical_phase_dim3_rad - numerical_phase_dim2_rad)))
    )
    max_wrapped_hamiltonian_error = float(np.max(np.abs(wrap_phase(hamiltonian_phase_rad - analytic_phase_rad))))
    max_wrapped_helper_vs_hamiltonian = float(
        np.max(np.abs(wrap_phase(hamiltonian_phase_rad - numerical_phase_dim3_rad)))
    )
    manifold_shift_rad_s = float(
        model_dim3.manifold_transition_frequency(1, None) - model_dim3.manifold_transition_frequency(0, None)
    )
    hamiltonian_shift_rad_s = float(
        hamiltonian_model.basis_energy(1, 1, MODEL_FRAME)
        - hamiltonian_model.basis_energy(1, 0, MODEL_FRAME)
        - hamiltonian_model.basis_energy(0, 1, MODEL_FRAME)
        + hamiltonian_model.basis_energy(0, 0, MODEL_FRAME)
    )

    save_phase_samples(
        times_ns=times_ns,
        analytic_phase_rad=analytic_phase_rad,
        numerical_phase_dim2_rad=numerical_phase_dim2_rad,
        numerical_phase_dim3_rad=numerical_phase_dim3_rad,
    )
    save_cross_check_samples(
        times_ns=times_ns,
        analytic_phase_rad=analytic_phase_rad,
        helper_phase_rad=numerical_phase_dim3_rad,
        hamiltonian_phase_rad=hamiltonian_phase_rad,
    )
    save_figure(
        times_ns=times_ns,
        analytic_phase_rad=analytic_phase_rad,
        numerical_phase_dim2_rad=numerical_phase_dim2_rad,
        numerical_phase_dim3_rad=numerical_phase_dim3_rad,
        t_pi_ns=t_pi_s * 1.0e9,
    )
    save_cross_check_figure(
        times_ns=times_ns,
        analytic_phase_rad=analytic_phase_rad,
        helper_phase_rad=numerical_phase_dim3_rad,
        hamiltonian_phase_rad=hamiltonian_phase_rad,
    )

    summary = {
        "study_name": STUDY_NAME,
        "date_created": _timestamp_utc(),
        "description": "Minimal model-level verification that free dispersive idle evolution reaches a pi relative phase between the n=0 and n=1 manifolds at t_pi = pi / |chi|.",
        "parameters": {
            "omega_q_hz": OMEGA_Q_HZ,
            "omega_c_hz": OMEGA_C_HZ,
            "alpha_hz": ALPHA_HZ,
            "chi_hz": CHI_HZ,
            "kerr_hz": KERR_HZ,
            "qubit_dim": QUBIT_DIM,
            "cavity_dims": list(CAVITY_DIMS),
            "hamiltonian_cross_check_cavity_dim": HAMILTONIAN_CROSS_CHECK_DIM,
            "num_samples": NUM_SAMPLES,
            "dispersive_phase_convention": "n_e",
        },
        "load_instructions": "Use json.loads(Path('artifacts/free_dispersive_pi_probe_summary.json').read_text(encoding='ascii')) for metadata and numpy/genfromtxt on data/phase_difference_samples.csv for the sampled phases.",
        "analytic": {
            "phase_model": "delta_phi(t) = -chi * t",
            "t_pi_s": t_pi_s,
            "t_pi_ns": t_pi_s * 1.0e9,
        },
        "numerical": {
            "max_wrapped_error_dim2_rad": max_wrapped_error_dim2,
            "max_wrapped_difference_dim2_vs_dim3_rad": max_wrapped_difference,
            "max_wrapped_hamiltonian_error_rad": max_wrapped_hamiltonian_error,
            "max_wrapped_helper_vs_hamiltonian_rad": max_wrapped_helper_vs_hamiltonian,
            "manifold_shift_rad_s": manifold_shift_rad_s,
            "manifold_shift_hz": manifold_shift_rad_s / TWO_PI,
            "hamiltonian_shift_rad_s": hamiltonian_shift_rad_s,
            "hamiltonian_shift_hz": hamiltonian_shift_rad_s / TWO_PI,
            "sample_count": NUM_SAMPLES,
        },
        "files": {
            "data_csv": str(DATA_PATH.relative_to(STUDY_ROOT)).replace("\\", "/"),
            "figure_png": str(FIGURE_STEM.with_suffix(".png").relative_to(STUDY_ROOT)).replace("\\", "/"),
            "figure_pdf": str(FIGURE_STEM.with_suffix(".pdf").relative_to(STUDY_ROOT)).replace("\\", "/"),
            "cross_check_data_csv": str(CROSS_CHECK_DATA_PATH.relative_to(STUDY_ROOT)).replace("\\", "/"),
            "cross_check_artifact_json": str(CROSS_CHECK_ARTIFACT_PATH.relative_to(STUDY_ROOT)).replace("\\", "/"),
            "cross_check_figure_png": str(CROSS_CHECK_FIGURE_STEM.with_suffix(".png").relative_to(STUDY_ROOT)).replace("\\", "/"),
            "cross_check_figure_pdf": str(CROSS_CHECK_FIGURE_STEM.with_suffix(".pdf").relative_to(STUDY_ROOT)).replace("\\", "/"),
        },
    }
    save_summary(summary)

    cross_check_summary = {
        "study_name": STUDY_NAME,
        "date_created": _timestamp_utc(),
        "description": "Independent static-Hamiltonian evolution cross-check for the model-level free-dispersive phase probe.",
        "parameters": {
            "chi_hz": CHI_HZ,
            "kerr_hz": KERR_HZ,
            "qubit_dim": QUBIT_DIM,
            "helper_cavity_dim": CAVITY_DIMS[1],
            "hamiltonian_cavity_dim": HAMILTONIAN_CROSS_CHECK_DIM,
            "num_samples": NUM_SAMPLES,
        },
        "load_instructions": "Use json.loads(Path('artifacts/free_dispersive_hamiltonian_cross_check.json').read_text(encoding='ascii')) for metadata and numpy/genfromtxt on data/phase_difference_hamiltonian_cross_check.csv for the sampled comparison.",
        "comparison": {
            "t_pi_ns": t_pi_s * 1.0e9,
            "helper_max_wrapped_error_rad": max_wrapped_error_dim3,
            "hamiltonian_max_wrapped_error_rad": max_wrapped_hamiltonian_error,
            "helper_vs_hamiltonian_max_wrapped_difference_rad": max_wrapped_helper_vs_hamiltonian,
            "hamiltonian_shift_hz": hamiltonian_shift_rad_s / TWO_PI,
        },
        "files": {
            "data_csv": str(CROSS_CHECK_DATA_PATH.relative_to(STUDY_ROOT)).replace("\\", "/"),
            "figure_png": str(CROSS_CHECK_FIGURE_STEM.with_suffix(".png").relative_to(STUDY_ROOT)).replace("\\", "/"),
            "figure_pdf": str(CROSS_CHECK_FIGURE_STEM.with_suffix(".pdf").relative_to(STUDY_ROOT)).replace("\\", "/"),
        },
    }
    with CROSS_CHECK_ARTIFACT_PATH.open("w", encoding="ascii") as handle:
        json.dump(cross_check_summary, handle, indent=2)
        handle.write("\n")

    elapsed_s = time.perf_counter() - start_time
    print(f"Analytic t_pi: {t_pi_s * 1.0e9:.6f} ns")
    print(f"Max wrapped analytic-vs-helper error (n_cav=2): {max_wrapped_error_dim2:.3e} rad")
    print(f"Max wrapped dim-2 vs dim-3 difference: {max_wrapped_difference:.3e} rad")
    print(f"Max wrapped analytic-vs-Hamiltonian error: {max_wrapped_hamiltonian_error:.3e} rad")
    print(f"Max wrapped helper-vs-Hamiltonian difference: {max_wrapped_helper_vs_hamiltonian:.3e} rad")
    print(f"Model manifold shift: {manifold_shift_rad_s / TWO_PI:.6f} Hz")
    print(f"Hamiltonian branch shift: {hamiltonian_shift_rad_s / TWO_PI:.6f} Hz")
    print(f"Runtime: {elapsed_s:.3f} s")
    print(f"Saved data: {DATA_PATH}")
    print(f"Saved cross-check data: {CROSS_CHECK_DATA_PATH}")
    print(f"Saved artifact: {ARTIFACT_PATH}")
    print(f"Saved cross-check artifact: {CROSS_CHECK_ARTIFACT_PATH}")
    print(f"Saved figure: {FIGURE_STEM}.png and {FIGURE_STEM}.pdf")
    print(f"Saved cross-check figure: {CROSS_CHECK_FIGURE_STEM}.png and {CROSS_CHECK_FIGURE_STEM}.pdf")


if __name__ == "__main__":
    main()