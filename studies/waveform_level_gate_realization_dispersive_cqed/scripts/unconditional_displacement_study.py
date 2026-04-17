"""Focused study of unconditional cavity displacement in dispersive cQED.

This script extends the existing waveform-level study with the protocol families
requested in the unconditional-displacement prompt:

1. Naive single-tone displacement
2. Fast single-tone pulses with multiple envelope families
3. Two-tone branch-compensated displacement
4. Echoed displacement
5. A bounded hardware-aware optimal-control benchmark
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Any, Callable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import qutip as qt
from scipy.linalg import norm

import common
from common import (
    ARTIFACTS_DIR,
    CHI,
    CHI_PRIME,
    DEFAULT_DT,
    FIGURES_DIR,
    KERR,
    N_CAV,
    N_TR,
    SIGMA_FRACTION,
    TOL_BRIGHT,
    annihilation_expectation,
    apply_plot_style,
    build_frame,
    build_model,
    cavity_branch_transition_frequency,
    cavity_state_from_joint,
    compile_and_prepare,
    displacement_envelope,
    displacement_fidelity,
    displacement_op,
    entanglement_entropy,
    envelope_area,
    make_gaussian_qubit_pulse,
    make_shaped_displacement_pulse,
    save_json,
    simulate_state,
    state_trace_distance,
)
from cqed_sim import (
    BoundaryWindowHardwareMap,
    FirstOrderLowPassHardwareMap,
    GrapeConfig,
    GrapeSolver,
    HardwareModel,
    HeldSampleParameterization,
    ModelControlChannelSpec,
    PiecewiseConstantTimeGrid,
    SmoothIQRadiusLimitHardwareMap,
    build_control_problem_from_model,
    multi_state_transfer_objective,
)
from cqed_sim.pulses.hardware import apply_first_order_lowpass


apply_plot_style()

ALPHA_TARGETS = np.array([0.5, 1.0, 1.5, 2.0], dtype=float)
SINGLE_TONE_DURATIONS_NS = np.array([5.0, 10.0, 20.0, 40.0, 80.0, 160.0], dtype=float)
ECHO_TOTAL_DURATIONS_NS = np.array([60.0, 80.0, 120.0, 160.0], dtype=float)
OPTIMAL_DURATIONS_NS = np.array([40.0, 80.0], dtype=float)
PULSE_FAMILIES = ("square", "gaussian", "cosine")
FILTER_BWS_MHZ = (None, 80.0, 40.0)
MODEL_VARIANTS = ("minimal", "higher_order", "full")
CHI_SCALES = np.array([0.5, 1.0, 1.5, 2.0], dtype=float)
TWO_TONE_CHI_ERRORS = np.array([-0.10, -0.05, 0.0, 0.05, 0.10], dtype=float)
PI_PULSE_ERROR_SCALES = np.array([0.95, 1.0, 1.05], dtype=float)
LOGICAL_CAVITY_LEVELS = tuple(range(4))
WIGNER_XMAX = 4.5
WIGNER_POINTS = 121
DRAG_PI = 0.5e-9
OPTIMAL_N_CAV = 12

FULL_STATE_TEST_SET = (
    ("g", "vacuum"),
    ("e", "vacuum"),
    ("plus_x", "vacuum"),
    ("plus_y", "vacuum"),
    ("g", "fock1"),
    ("e", "fock1"),
    ("g", "fock2"),
    ("e", "fock2"),
    ("g", "fock3"),
    ("e", "fock3"),
    ("g", "coherent"),
    ("e", "coherent"),
    ("plus_x", "coherent"),
    ("plus_y", "coherent"),
)

OPTIMIZATION_STATE_SET = (
    ("g", "vacuum"),
    ("e", "vacuum"),
    ("plus_x", "vacuum"),
    ("plus_y", "vacuum"),
    ("g", "fock1"),
    ("e", "fock1"),
    ("g", "coherent_opt"),
    ("e", "coherent_opt"),
)


def variant_model(variant: str, *, chi_scale: float = 1.0, n_cav: int = N_CAV, n_tr: int = N_TR):
    if variant == "minimal":
        return build_model(n_cav=n_cav, n_tr=n_tr, chi=CHI * chi_scale, chi_prime=None, kerr=0.0)
    if variant == "higher_order":
        return build_model(n_cav=n_cav, n_tr=n_tr, chi=CHI * chi_scale, chi_prime=CHI_PRIME, kerr=0.0)
    if variant == "full":
        return build_model(n_cav=n_cav, n_tr=n_tr, chi=CHI * chi_scale, chi_prime=CHI_PRIME, kerr=KERR)
    raise ValueError(f"Unknown model variant '{variant}'.")


def qubit_state(model, label: str) -> qt.Qobj:
    if label == "g":
        return qt.basis(model.n_tr, 0)
    if label == "e":
        return qt.basis(model.n_tr, 1)
    if label == "plus_x":
        return (qt.basis(model.n_tr, 0) + qt.basis(model.n_tr, 1)).unit()
    if label == "plus_y":
        return (qt.basis(model.n_tr, 0) + 1j * qt.basis(model.n_tr, 1)).unit()
    raise ValueError(f"Unsupported qubit state '{label}'.")


def cavity_state(model, label: str) -> qt.Qobj:
    if label == "vacuum":
        return qt.basis(model.n_cav, 0)
    if label == "fock1":
        return qt.basis(model.n_cav, 1)
    if label == "fock2":
        return qt.basis(model.n_cav, 2)
    if label == "fock3":
        return qt.basis(model.n_cav, 3)
    if label == "coherent":
        return qt.coherent(model.n_cav, 0.75)
    if label == "coherent_opt":
        return qt.coherent(model.n_cav, 0.5)
    raise ValueError(f"Unsupported cavity state '{label}'.")


def joint_state(model, qubit_label: str, cavity_label: str) -> qt.Qobj:
    return qt.tensor(qubit_state(model, qubit_label), cavity_state(model, cavity_label))


def ideal_displaced_state(model, qubit_label: str, cavity_label: str, alpha_target: complex) -> qt.Qobj:
    disp = displacement_op(model.n_cav, complex(alpha_target))
    cav = disp * cavity_state(model, cavity_label)
    return qt.tensor(qubit_state(model, qubit_label), cav)


def branch_overlap(final_state: qt.Qobj, model) -> float:
    rho_g, pop_g = cavity_state_from_joint(final_state, model.n_tr, model.n_cav, 0)
    rho_e, pop_e = cavity_state_from_joint(final_state, model.n_tr, model.n_cav, 1)
    if pop_g <= 1e-12 or pop_e <= 1e-12:
        return 0.0
    return float(qt.metrics.fidelity(rho_g, rho_e) ** 2)


def cavity_logical_population(final_state: qt.Qobj, model, *, max_level: int = 3) -> float:
    rho = final_state if final_state.isoper else final_state.proj()
    rho_c = rho.ptrace(1)
    diag = np.real(np.diag(np.asarray(rho_c.full(), dtype=np.complex128)))
    return float(np.sum(diag[: max_level + 1]))


def logical_probe_states(model) -> list[tuple[str, qt.Qobj]]:
    states: list[tuple[str, qt.Qobj]] = []
    for q_label in ("g", "e"):
        for level in LOGICAL_CAVITY_LEVELS:
            cav_label = "vacuum" if level == 0 else f"fock{level}"
            states.append((f"{q_label}_n{level}", joint_state(model, q_label, cav_label)))
    for q_label in ("g", "e", "plus_x", "plus_y"):
        pieces = [joint_state(model, q_label, "vacuum")]
        for level in LOGICAL_CAVITY_LEVELS[1:]:
            pieces.append(joint_state(model, q_label, f"fock{level}"))
        states.append((f"{q_label}_uniform", sum(pieces).unit()))
    return states


def power_bandwidth_mhz(samples: np.ndarray, dt_s: float, *, fraction: float = 0.95) -> float:
    arr = np.asarray(samples, dtype=np.complex128).reshape(-1)
    if arr.size <= 1:
        return 0.0
    spec = np.fft.fftshift(np.fft.fft(arr))
    freqs = np.fft.fftshift(np.fft.fftfreq(arr.size, d=dt_s))
    power = np.abs(spec) ** 2
    order = np.argsort(np.abs(freqs))
    cumulative = np.cumsum(power[order])
    total = cumulative[-1] if cumulative.size else 1.0
    idx = int(np.searchsorted(cumulative, fraction * total))
    idx = min(idx, len(order) - 1)
    return float(abs(freqs[order[idx]]) / 1.0e6)


def build_filtered_displacement_pulse(
    *,
    alpha: complex,
    duration_s: float,
    family: str,
    carrier: float = 0.0,
    filter_bw_mhz: float | None = None,
    t0: float = 0.0,
    label: str | None = None,
) -> tuple[Any, dict[str, Any]]:
    env = displacement_envelope(family, sigma_fraction=SIGMA_FRACTION, rise_fraction=0.15)
    area = envelope_area(env)
    amp_complex = 1j * complex(alpha) / (duration_s * max(area, 1.0e-15))
    n_steps = max(2, int(round(duration_s / DEFAULT_DT)))
    t_rel = np.arange(n_steps, dtype=float) / float(n_steps)
    baseband = amp_complex * np.asarray(env(t_rel), dtype=np.complex128)
    filtered = apply_first_order_lowpass(
        baseband,
        DEFAULT_DT,
        None if filter_bw_mhz is None else float(filter_bw_mhz) * 1.0e6,
    )
    pulse = common.Pulse(
        channel="storage",
        t0=float(t0),
        duration=float(duration_s),
        envelope=np.asarray(filtered, dtype=np.complex128),
        carrier=float(carrier),
        phase=0.0,
        amp=1.0,
        sample_rate=1.0 / DEFAULT_DT,
        label=label or f"{family}_filtered_disp",
    )
    meta = {
        "command_baseband": baseband,
        "physical_baseband": filtered,
        "filter_bw_mhz": None if filter_bw_mhz is None else float(filter_bw_mhz),
        "bandwidth_mhz": power_bandwidth_mhz(filtered, DEFAULT_DT),
        "distortion_l2_fraction": float(norm(filtered - baseband) / max(norm(baseband), 1.0e-15)),
    }
    return pulse, meta


def single_tone_pulse(
    *,
    model,
    frame,
    alpha: complex,
    duration_s: float,
    family: str,
    filter_bw_mhz: float | None,
    carrier: float = 0.0,
    t0: float = 0.0,
    label: str | None = None,
) -> tuple[list[Any], dict[str, Any]]:
    if filter_bw_mhz is None:
        pulse = make_shaped_displacement_pulse(
            model,
            frame,
            alpha=alpha,
            duration_s=duration_s,
            family=family,
            carrier=carrier,
            t0=t0,
            label=label,
        )
        env = displacement_envelope(family, sigma_fraction=SIGMA_FRACTION, rise_fraction=0.15)
        area = envelope_area(env)
        amp_complex = 1j * complex(alpha) / (duration_s * max(area, 1.0e-15))
        n_steps = max(2, int(round(duration_s / DEFAULT_DT)))
        t_rel = np.arange(n_steps, dtype=float) / float(n_steps)
        baseband = amp_complex * np.asarray(env(t_rel), dtype=np.complex128)
        meta = {
            "command_baseband": baseband,
            "physical_baseband": baseband,
            "filter_bw_mhz": None,
            "bandwidth_mhz": power_bandwidth_mhz(baseband, DEFAULT_DT),
            "distortion_l2_fraction": 0.0,
        }
        return [pulse], meta
    pulse, meta = build_filtered_displacement_pulse(
        alpha=alpha,
        duration_s=duration_s,
        family=family,
        carrier=carrier,
        filter_bw_mhz=filter_bw_mhz,
        t0=t0,
        label=label,
    )
    return [pulse], meta


def effective_branch_alpha(final_state: qt.Qobj) -> complex:
    return annihilation_expectation(final_state)


def branch_vacuum_metrics(model, session, *, alpha_target: complex) -> dict[str, Any]:
    psi_g = simulate_state(session, model.basis_state(0, 0))
    psi_e = simulate_state(session, model.basis_state(1, 0))
    psi_px = simulate_state(session, joint_state(model, "plus_x", "vacuum"))
    psi_py = simulate_state(session, joint_state(model, "plus_y", "vacuum"))

    target_g = ideal_displaced_state(model, "g", "vacuum", alpha_target)
    target_e = ideal_displaced_state(model, "e", "vacuum", alpha_target)
    target_px = ideal_displaced_state(model, "plus_x", "vacuum", alpha_target)
    target_py = ideal_displaced_state(model, "plus_y", "vacuum", alpha_target)

    alpha_g = effective_branch_alpha(psi_g)
    alpha_e = effective_branch_alpha(psi_e)
    return {
        "alpha_g": alpha_g,
        "alpha_e": alpha_e,
        "delta_alpha": float(abs(alpha_g - alpha_e)),
        "g_fidelity": displacement_fidelity(psi_g, target_g),
        "e_fidelity": displacement_fidelity(psi_e, target_e),
        "plus_x_fidelity": displacement_fidelity(psi_px, target_px),
        "plus_y_fidelity": displacement_fidelity(psi_py, target_py),
        "plus_x_entanglement_bits": entanglement_entropy(psi_px, model.n_tr, model.n_cav),
        "plus_y_entanglement_bits": entanglement_entropy(psi_py, model.n_tr, model.n_cav),
        "plus_x_branch_overlap": branch_overlap(psi_px, model),
        "plus_y_branch_overlap": branch_overlap(psi_py, model),
    }


def protocol_metrics_from_session(
    model,
    session,
    *,
    alpha_target: complex,
    state_pairs: tuple[tuple[str, str], ...],
) -> dict[str, Any]:
    results: dict[str, Any] = {"state_fidelities": {}}
    state_fidelity_values: list[float] = []
    for qubit_label, cavity_label in state_pairs:
        psi0 = joint_state(model, qubit_label, cavity_label)
        psi_f = simulate_state(session, psi0)
        psi_t = ideal_displaced_state(model, qubit_label, cavity_label, alpha_target)
        fidelity = displacement_fidelity(psi_f, psi_t)
        state_fidelity_values.append(float(fidelity))
        entry: dict[str, Any] = {"fidelity": fidelity}
        if qubit_label in {"plus_x", "plus_y"}:
            entry["entanglement_entropy_bits"] = entanglement_entropy(psi_f, model.n_tr, model.n_cav)
            entry["branch_overlap"] = branch_overlap(psi_f, model)
        rho_c = psi_f.ptrace(1)
        ideal_c = psi_t.ptrace(1)
        entry["cavity_trace_distance"] = state_trace_distance(rho_c, ideal_c)
        results["state_fidelities"][f"{qubit_label}|{cavity_label}"] = entry
    results["state_test_mean_fidelity"] = float(np.mean(state_fidelity_values))
    results["state_test_min_fidelity"] = float(np.min(state_fidelity_values))

    logical_fids: list[float] = []
    logical_pops: list[float] = []
    for _label, psi0 in logical_probe_states(model):
        psi_f = simulate_state(session, psi0)
        rho_q = np.asarray(psi0.ptrace(0).full(), dtype=np.complex128)
        evals_q, evecs_q = np.linalg.eigh(rho_q)
        q_state = qt.Qobj(evecs_q[:, np.argmax(evals_q)], dims=[[model.n_tr], [1]])
        rho_c = np.asarray(psi0.ptrace(1).full(), dtype=np.complex128)
        evals_c, evecs_c = np.linalg.eigh(rho_c)
        c_state = qt.Qobj(evecs_c[:, np.argmax(evals_c)], dims=[[model.n_cav], [1]])
        psi_t = qt.tensor(q_state, displacement_op(model.n_cav, alpha_target) * c_state)
        logical_fids.append(displacement_fidelity(psi_f, psi_t))
        logical_pops.append(cavity_logical_population(psi_f, model, max_level=max(LOGICAL_CAVITY_LEVELS)))
    results["logical_average_fidelity"] = float(np.mean(logical_fids))
    results["logical_min_fidelity"] = float(np.min(logical_fids))
    results["logical_retained_population"] = float(np.mean(logical_pops))
    return results


def calibrate_two_tone(
    *,
    model,
    frame,
    alpha_target: complex,
    duration_s: float,
    family: str = "gaussian",
    chi_error_fraction: float = 0.0,
) -> dict[str, Any]:
    det_g = cavity_branch_transition_frequency(model, frame, qubit_level=0)
    det_e = cavity_branch_transition_frequency(model, frame, qubit_level=1) * (1.0 + float(chi_error_fraction))
    carriers = (-det_g, -det_e)
    responses = np.zeros((2, 2), dtype=np.complex128)
    for tone_index, carrier in enumerate(carriers):
        pulses, _meta = single_tone_pulse(
            model=model,
            frame=frame,
            alpha=1.0 + 0.0j,
            duration_s=duration_s,
            family=family,
            filter_bw_mhz=None,
            carrier=carrier,
            label=f"two_tone_basis_{tone_index}",
        )
        session = compile_and_prepare(model, frame, pulses)
        responses[0, tone_index] = effective_branch_alpha(simulate_state(session, model.basis_state(0, 0)))
        responses[1, tone_index] = effective_branch_alpha(simulate_state(session, model.basis_state(1, 0)))
    target_vec = np.array([complex(alpha_target), complex(alpha_target)], dtype=np.complex128)
    weights = np.linalg.solve(responses, target_vec)
    pulses: list[Any] = []
    for tone_index, (carrier, weight) in enumerate(zip(carriers, weights, strict=True)):
        tone_pulses, _meta = single_tone_pulse(
            model=model,
            frame=frame,
            alpha=weight,
            duration_s=duration_s,
            family=family,
            filter_bw_mhz=None,
            carrier=carrier,
            label=f"two_tone_{tone_index}",
        )
        pulses.extend(tone_pulses)
    return {
        "pulses": pulses,
        "meta": {
            "weights": weights,
            "responses": responses,
            "carrier_detunings": np.array(carriers, dtype=float),
        },
    }


def build_echo_sequence(
    *,
    model,
    frame,
    alpha_target: complex,
    total_duration_s: float,
    pi_duration_s: float,
    drag: float,
    pi_amp_scale: float = 1.0,
) -> tuple[list[Any], dict[str, Any]]:
    displacement_duration = total_duration_s - 2.0 * pi_duration_s
    if displacement_duration <= 0.0:
        raise ValueError("Echo total duration must exceed 2*pi_duration.")
    half_disp = 0.5 * displacement_duration
    d1, _ = single_tone_pulse(
        model=model,
        frame=frame,
        alpha=0.5 * alpha_target,
        duration_s=half_disp,
        family="square",
        filter_bw_mhz=None,
        carrier=0.0,
        t0=0.0,
        label="echo_D1",
    )
    d2, _ = single_tone_pulse(
        model=model,
        frame=frame,
        alpha=0.5 * alpha_target,
        duration_s=half_disp,
        family="square",
        filter_bw_mhz=None,
        carrier=0.0,
        t0=half_disp + pi_duration_s,
        label="echo_D2",
    )
    pi1 = make_gaussian_qubit_pulse(
        model,
        frame,
        theta=np.pi,
        phase=0.0,
        duration_s=pi_duration_s,
        manifold_level=0,
        drag=drag,
        t0=half_disp,
        label="echo_X1",
    )
    pi2 = make_gaussian_qubit_pulse(
        model,
        frame,
        theta=np.pi,
        phase=0.0,
        duration_s=pi_duration_s,
        manifold_level=0,
        drag=drag,
        t0=2.0 * half_disp + pi_duration_s,
        label="echo_X2",
    )
    pi1 = common.Pulse(
        pi1.channel, pi1.t0, pi1.duration, pi1.envelope, pi1.carrier, pi1.phase,
        pi1.amp * float(pi_amp_scale), pi1.drag, pi1.sample_rate, pi1.label,
    )
    pi2 = common.Pulse(
        pi2.channel, pi2.t0, pi2.duration, pi2.envelope, pi2.carrier, pi2.phase,
        pi2.amp * float(pi_amp_scale), pi2.drag, pi2.sample_rate, pi2.label,
    )
    return [d1[0], pi1, d2[0], pi2], {
        "total_duration_ns": total_duration_s * 1.0e9,
        "pi_duration_ns": pi_duration_s * 1.0e9,
        "pi_amp_scale": float(pi_amp_scale),
    }


def solve_optimal_control_case(*, duration_s: float, alpha_target: complex) -> dict[str, Any]:
    model = variant_model("full", n_cav=OPTIMAL_N_CAV)
    frame = build_frame(model)
    initial_states = [joint_state(model, q_label, c_label) for q_label, c_label in OPTIMIZATION_STATE_SET]
    target_states = [ideal_displaced_state(model, q_label, c_label, alpha_target) for q_label, c_label in OPTIMIZATION_STATE_SET]
    time_grid = PiecewiseConstantTimeGrid.uniform(steps=8, dt_s=float(duration_s / 8.0))
    sample_period_s = max(duration_s / 4.0, 10.0e-9)
    held_steps = max(1, int(np.ceil(duration_s / sample_period_s - 1.0e-15)))
    q_init = float(np.real(1j * alpha_target / duration_s))
    initial_schedule = np.zeros((2, held_steps), dtype=float)
    initial_schedule[1, :] = q_init
    problem = build_control_problem_from_model(
        model,
        frame=frame,
        time_grid=time_grid,
        channel_specs=(
            ModelControlChannelSpec(
                name="storage",
                target="storage",
                quadratures=("I", "Q"),
                amplitude_bounds=(-1.8e8, 1.8e8),
                export_channel="storage",
            ),
        ),
        objectives=(multi_state_transfer_objective(initial_states, target_states, name="unconditional_displacement"),),
        parameterization_cls=HeldSampleParameterization,
        parameterization_kwargs={"sample_period_s": sample_period_s},
        hardware_model=HardwareModel(
            maps=(
                FirstOrderLowPassHardwareMap(cutoff_hz=60.0e6, export_channels=("storage",)),
                SmoothIQRadiusLimitHardwareMap(amplitude_max=1.8e8, export_channels=("storage",)),
                BoundaryWindowHardwareMap(ramp_slices=1, export_channels=("storage",)),
            )
        ),
    )
    result = GrapeSolver(
        GrapeConfig(maxiter=60, seed=11, random_scale=0.05, report_command_reference=True)
    ).solve(problem, initial_schedule=initial_schedule)
    pulses, drive_ops, pulse_meta = result.to_pulses(waveform="physical")
    session = compile_and_prepare(model, frame, pulses, drive_ops=drive_ops)
    return {
        "duration_ns": duration_s * 1.0e9,
        "success": bool(result.success),
        "objective_value": float(result.objective_value),
        "metrics": dict(result.metrics),
        "hardware_metrics": dict(result.hardware_metrics),
        "command_values": np.asarray(result.command_values, dtype=float),
        "physical_values": np.asarray(result.physical_values, dtype=float),
        "schedule_values": np.asarray(result.schedule.values, dtype=float),
        "vacuum_metrics": branch_vacuum_metrics(model, session, alpha_target=alpha_target),
        "full_metrics": protocol_metrics_from_session(
            model, session, alpha_target=alpha_target, state_pairs=FULL_STATE_TEST_SET
        ),
        "pulse_export_meta": pulse_meta,
    }


@dataclass
class ProtocolSnapshot:
    label: str
    family: str
    duration_ns: float
    alpha_target: float
    build: Callable[[Any, Any], tuple[list[Any], dict[str, Any]]]
    complexity: str


def evaluate_snapshot(snapshot: ProtocolSnapshot) -> dict[str, Any]:
    model = variant_model("full")
    frame = build_frame(model)
    pulses, meta = snapshot.build(model, frame)
    drive_ops = {"storage": "cavity", "qubit": "qubit"} if any(p.channel == "qubit" for p in pulses) else {"storage": "cavity"}
    session = compile_and_prepare(model, frame, pulses, drive_ops=drive_ops)
    return {
        "label": snapshot.label,
        "family": snapshot.family,
        "duration_ns": snapshot.duration_ns,
        "alpha_target": snapshot.alpha_target,
        "complexity": snapshot.complexity,
        "vacuum_metrics": branch_vacuum_metrics(model, session, alpha_target=snapshot.alpha_target),
        "full_metrics": protocol_metrics_from_session(
            model, session, alpha_target=snapshot.alpha_target, state_pairs=FULL_STATE_TEST_SET
        ),
        "pulse_meta": meta,
    }


def run_single_tone_sweeps() -> dict[str, Any]:
    print("=" * 72)
    print("Single-tone unconditional-displacement sweeps")
    print("=" * 72)
    t_start = time.time()
    model = variant_model("full")
    frame = build_frame(model)
    baseline: dict[str, Any] = {}
    for family in PULSE_FAMILIES:
        baseline[family] = {}
        for filter_bw in FILTER_BWS_MHZ:
            key = "nofilter" if filter_bw is None else f"lp_{int(filter_bw)}mhz"
            rows = []
            for duration_ns in SINGLE_TONE_DURATIONS_NS:
                duration_s = duration_ns * 1.0e-9
                for alpha_target in ALPHA_TARGETS:
                    pulses, pulse_meta = single_tone_pulse(
                        model=model,
                        frame=frame,
                        alpha=alpha_target,
                        duration_s=duration_s,
                        family=family,
                        filter_bw_mhz=filter_bw,
                    )
                    session = compile_and_prepare(model, frame, pulses)
                    metrics = branch_vacuum_metrics(model, session, alpha_target=alpha_target)
                    rows.append(
                        {
                            "duration_ns": float(duration_ns),
                            "alpha_target": float(alpha_target),
                            **metrics,
                            "bandwidth_mhz": pulse_meta["bandwidth_mhz"],
                            "distortion_l2_fraction": pulse_meta["distortion_l2_fraction"],
                        }
                    )
                print(f"  {family:8s} {key:10s} T={duration_ns:6.1f} ns complete")
            baseline[family][key] = rows

    hierarchy_rows = []
    for variant in MODEL_VARIANTS:
        v_model = variant_model(variant)
        v_frame = build_frame(v_model)
        for duration_ns in SINGLE_TONE_DURATIONS_NS:
            for alpha_target in ALPHA_TARGETS:
                pulses, _ = single_tone_pulse(
                    model=v_model,
                    frame=v_frame,
                    alpha=alpha_target,
                    duration_s=duration_ns * 1.0e-9,
                    family="square",
                    filter_bw_mhz=None,
                )
                session = compile_and_prepare(v_model, v_frame, pulses)
                hierarchy_rows.append(
                    {
                        "variant": variant,
                        "duration_ns": float(duration_ns),
                        "alpha_target": float(alpha_target),
                        **branch_vacuum_metrics(v_model, session, alpha_target=alpha_target),
                    }
                )

    chi_rows = []
    for chi_scale in CHI_SCALES:
        c_model = variant_model("full", chi_scale=chi_scale)
        c_frame = build_frame(c_model)
        for duration_ns in SINGLE_TONE_DURATIONS_NS:
            pulses, _ = single_tone_pulse(
                model=c_model,
                frame=c_frame,
                alpha=1.0,
                duration_s=duration_ns * 1.0e-9,
                family="square",
                filter_bw_mhz=None,
            )
            session = compile_and_prepare(c_model, c_frame, pulses)
            chi_rows.append(
                {
                    "chi_scale": float(chi_scale),
                    "duration_ns": float(duration_ns),
                    **branch_vacuum_metrics(c_model, session, alpha_target=1.0),
                }
            )

    payload = {
        "baseline": baseline,
        "model_hierarchy": hierarchy_rows,
        "chi_scaling": chi_rows,
        "wall_time_s": time.time() - t_start,
    }
    save_json(ARTIFACTS_DIR / "unconditional_single_tone_summary.json", payload)
    return payload


def run_two_tone_sweep() -> dict[str, Any]:
    print("=" * 72)
    print("Two-tone branch-compensated displacement sweep")
    print("=" * 72)
    t_start = time.time()
    model = variant_model("full")
    frame = build_frame(model)
    rows = []
    robustness_rows = []
    for duration_ns in SINGLE_TONE_DURATIONS_NS[2:]:
        for alpha_target in ALPHA_TARGETS:
            calibration = calibrate_two_tone(
                model=model,
                frame=frame,
                alpha_target=alpha_target,
                duration_s=duration_ns * 1.0e-9,
            )
            session = compile_and_prepare(model, frame, calibration["pulses"])
            rows.append(
                {
                    "duration_ns": float(duration_ns),
                    "alpha_target": float(alpha_target),
                    **branch_vacuum_metrics(model, session, alpha_target=alpha_target),
                    **calibration["meta"],
                }
            )
        print(f"  calibrated two-tone duration {duration_ns:6.1f} ns")
    for error_frac in TWO_TONE_CHI_ERRORS:
        calibration = calibrate_two_tone(
            model=model,
            frame=frame,
            alpha_target=1.0,
            duration_s=80.0e-9,
            chi_error_fraction=error_frac,
        )
        session = compile_and_prepare(model, frame, calibration["pulses"])
        robustness_rows.append(
            {
                "chi_error_fraction": float(error_frac),
                **branch_vacuum_metrics(model, session, alpha_target=1.0),
            }
        )
    payload = {
        "calibrated_sweep": rows,
        "chi_robustness": robustness_rows,
        "wall_time_s": time.time() - t_start,
    }
    save_json(ARTIFACTS_DIR / "unconditional_two_tone_summary.json", payload)
    return payload


def run_echo_sweep() -> dict[str, Any]:
    print("=" * 72)
    print("Echoed displacement sweep")
    print("=" * 72)
    t_start = time.time()
    model = variant_model("full")
    frame = build_frame(model)
    rows = []
    sensitivity_rows = []
    for duration_ns in ECHO_TOTAL_DURATIONS_NS:
        for alpha_target in ALPHA_TARGETS:
            pulses, meta = build_echo_sequence(
                model=model,
                frame=frame,
                alpha_target=alpha_target,
                total_duration_s=duration_ns * 1.0e-9,
                pi_duration_s=20.0e-9,
                drag=DRAG_PI,
            )
            session = compile_and_prepare(
                model, frame, pulses, drive_ops={"storage": "cavity", "qubit": "qubit"}
            )
            rows.append(
                {
                    "duration_ns": float(duration_ns),
                    "alpha_target": float(alpha_target),
                    **branch_vacuum_metrics(model, session, alpha_target=alpha_target),
                    **meta,
                }
            )
        print(f"  echo total duration {duration_ns:6.1f} ns complete")
    for scale in PI_PULSE_ERROR_SCALES:
        pulses, meta = build_echo_sequence(
            model=model,
            frame=frame,
            alpha_target=1.0,
            total_duration_s=120.0e-9,
            pi_duration_s=20.0e-9,
            drag=DRAG_PI,
            pi_amp_scale=scale,
        )
        session = compile_and_prepare(
            model, frame, pulses, drive_ops={"storage": "cavity", "qubit": "qubit"}
        )
        sensitivity_rows.append(
            {
                "pi_amp_scale": float(scale),
                **branch_vacuum_metrics(model, session, alpha_target=1.0),
                **meta,
            }
        )
    payload = {
        "echo_sweep": rows,
        "pi_error_sensitivity": sensitivity_rows,
        "wall_time_s": time.time() - t_start,
    }
    save_json(ARTIFACTS_DIR / "unconditional_echo_summary.json", payload)
    return payload


def run_optimal_control_benchmark() -> dict[str, Any]:
    print("=" * 72)
    print("Hardware-aware optimal-control benchmark")
    print("=" * 72)
    t_start = time.time()
    cases = []
    for duration_ns in OPTIMAL_DURATIONS_NS:
        print(f"  solving optimal-control case at {duration_ns:.1f} ns")
        cases.append(solve_optimal_control_case(duration_s=duration_ns * 1.0e-9, alpha_target=1.0))
    payload = {"cases": cases, "wall_time_s": time.time() - t_start}
    save_json(ARTIFACTS_DIR / "unconditional_optimal_control_summary.json", payload)
    return payload


def build_selected_snapshots(optimal_data: dict[str, Any]) -> list[ProtocolSnapshot]:
    best_optimal = max(
        optimal_data["cases"],
        key=lambda item: item["full_metrics"]["state_test_mean_fidelity"],
    )

    def optimal_builder(_model, _frame):
        physical = np.asarray(best_optimal["physical_values"], dtype=float)
        n_steps = physical.shape[1]
        dt_s = float(best_optimal["duration_ns"] * 1.0e-9 / n_steps)
        boundaries = PiecewiseConstantTimeGrid.uniform(steps=n_steps, dt_s=dt_s).boundaries_s()
        pulses = []
        for step in range(n_steps):
            coeff = physical[0, step] + 1j * physical[1, step]
            if abs(coeff) <= 1.0e-14:
                continue
            pulses.append(
                common.Pulse(
                    "storage",
                    float(boundaries[step]),
                    dt_s,
                    common.square_envelope,
                    0.0,
                    float(np.angle(coeff)),
                    float(abs(coeff)),
                    0.0,
                    None,
                    f"oc_{step}",
                )
            )
        return pulses, {"source_duration_ns": float(best_optimal["duration_ns"])}

    return [
        ProtocolSnapshot(
            label="Naive square",
            family="single_square",
            duration_ns=80.0,
            alpha_target=1.0,
            build=lambda model, frame: single_tone_pulse(
                model=model,
                frame=frame,
                alpha=1.0,
                duration_s=80.0e-9,
                family="square",
                filter_bw_mhz=None,
            ),
            complexity="Low",
        ),
        ProtocolSnapshot(
            label="Fast Gaussian",
            family="single_gaussian",
            duration_ns=20.0,
            alpha_target=1.0,
            build=lambda model, frame: single_tone_pulse(
                model=model,
                frame=frame,
                alpha=1.0,
                duration_s=20.0e-9,
                family="gaussian",
                filter_bw_mhz=None,
            ),
            complexity="Low",
        ),
        ProtocolSnapshot(
            label="Two-tone compensated",
            family="two_tone",
            duration_ns=80.0,
            alpha_target=1.0,
            build=lambda model, frame: (lambda cal: (cal["pulses"], cal["meta"]))(
                calibrate_two_tone(model=model, frame=frame, alpha_target=1.0, duration_s=80.0e-9)
            ),
            complexity="Medium",
        ),
        ProtocolSnapshot(
            label="Echoed displacement",
            family="echo",
            duration_ns=120.0,
            alpha_target=1.0,
            build=lambda model, frame: build_echo_sequence(
                model=model,
                frame=frame,
                alpha_target=1.0,
                total_duration_s=120.0e-9,
                pi_duration_s=20.0e-9,
                drag=DRAG_PI,
            ),
            complexity="Medium",
        ),
        ProtocolSnapshot(
            label="Optimal control",
            family="optimal",
            duration_ns=float(best_optimal["duration_ns"]),
            alpha_target=1.0,
            build=optimal_builder,
            complexity="High",
        ),
    ]


def run_selected_protocol_comparison(optimal_data: dict[str, Any]) -> dict[str, Any]:
    print("=" * 72)
    print("Selected protocol comparison")
    print("=" * 72)
    payload = {
        "protocols": [evaluate_snapshot(snapshot) for snapshot in build_selected_snapshots(optimal_data)]
    }
    save_json(ARTIFACTS_DIR / "unconditional_protocol_comparison.json", payload)
    return payload


def build_wigner_data(protocol_payload: dict[str, Any]) -> dict[str, Any]:
    xvec = np.linspace(-WIGNER_XMAX, WIGNER_XMAX, WIGNER_POINTS)
    model = variant_model("full")
    frame = build_frame(model)
    selected = {
        "naive_fail": next(
            item for item in protocol_payload["protocols"] if item["family"] == "single_square"
        ),
        "two_tone_success": next(
            item for item in protocol_payload["protocols"] if item["family"] == "two_tone"
        ),
    }
    cases = {}
    for key, summary in selected.items():
        if summary["family"] == "single_square":
            pulses, _ = single_tone_pulse(
                model=model,
                frame=frame,
                alpha=1.0,
                duration_s=80.0e-9,
                family="square",
                filter_bw_mhz=None,
            )
            drive_ops = {"storage": "cavity"}
        else:
            cal = calibrate_two_tone(
                model=model, frame=frame, alpha_target=1.0, duration_s=80.0e-9
            )
            pulses = cal["pulses"]
            drive_ops = {"storage": "cavity"}
        session = compile_and_prepare(model, frame, pulses, drive_ops=drive_ops)
        final_state = simulate_state(session, joint_state(model, "plus_x", "vacuum"))
        rho_c = final_state.ptrace(1)
        ideal = ideal_displaced_state(model, "plus_x", "vacuum", 1.0).ptrace(1)
        cases[key] = {
            "xvec": xvec,
            "actual_wigner": qt.wigner(rho_c, xvec, xvec),
            "ideal_wigner": qt.wigner(ideal, xvec, xvec),
            "trace_distance": state_trace_distance(rho_c, ideal),
        }
    save_json(ARTIFACTS_DIR / "unconditional_wigner_cases.json", cases)
    return cases


def rows_to_grid(
    rows: list[dict[str, Any]],
    *,
    metric: str,
    durations: np.ndarray,
    alphas: np.ndarray,
) -> np.ndarray:
    grid = np.zeros((len(durations), len(alphas)), dtype=float)
    for row in rows:
        i = int(np.where(np.isclose(durations, row["duration_ns"]))[0][0])
        j = int(np.where(np.isclose(alphas, row["alpha_target"]))[0][0])
        grid[i, j] = float(row[metric])
    return grid


def generate_figures(
    single_tone_data: dict[str, Any],
    two_tone_data: dict[str, Any],
    echo_data: dict[str, Any],
    protocol_payload: dict[str, Any],
    wigner_data: dict[str, Any],
) -> None:
    print("=" * 72)
    print("Generating unconditional-displacement figures")
    print("=" * 72)

    square_rows = single_tone_data["baseline"]["square"]["nofilter"]
    two_tone_rows = two_tone_data["calibrated_sweep"]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2))
    for ax, rows, title, durations in (
        (axes[0], square_rows, "Naive single tone", SINGLE_TONE_DURATIONS_NS),
        (axes[1], two_tone_rows, "Two-tone compensated", SINGLE_TONE_DURATIONS_NS[2:]),
    ):
        im = ax.pcolormesh(
            ALPHA_TARGETS,
            durations,
            rows_to_grid(rows, metric="delta_alpha", durations=durations, alphas=ALPHA_TARGETS),
            shading="nearest",
            cmap="viridis",
        )
        ax.set_xlabel(r"Target displacement $|\alpha|$")
        ax.set_ylabel("Duration (ns)")
        ax.set_title(title)
        plt.colorbar(im, ax=ax, label=r"$\delta \alpha$")
    fig.suptitle("Branch-resolved displacement mismatch", y=1.02)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "unconditional_branch_mismatch.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "unconditional_branch_mismatch.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    for family, color in zip(PULSE_FAMILIES, TOL_BRIGHT, strict=False):
        rows = [
            row
            for row in single_tone_data["baseline"][family]["nofilter"]
            if abs(row["alpha_target"] - 1.0) < 1e-12
        ]
        xs = np.array([row["duration_ns"] for row in rows], dtype=float)
        ys = np.array([row["plus_x_entanglement_bits"] for row in rows], dtype=float)
        order = np.argsort(xs)
        ax.plot(xs[order], ys[order], "o-", color=color, label=f"{family.title()} single tone")
    echo_rows = [
        row for row in echo_data["echo_sweep"] if abs(row["alpha_target"] - 1.0) < 1e-12
    ]
    ax.plot(
        [row["duration_ns"] for row in echo_rows],
        [row["plus_x_entanglement_bits"] for row in echo_rows],
        "s--",
        color=TOL_BRIGHT[4],
        label="Echoed",
    )
    tt_rows = [row for row in two_tone_rows if abs(row["alpha_target"] - 1.0) < 1e-12]
    ax.plot(
        [row["duration_ns"] for row in tt_rows],
        [row["plus_x_entanglement_bits"] for row in tt_rows],
        "d-",
        color=TOL_BRIGHT[2],
        label="Two-tone",
    )
    ax.set_xlabel("Duration (ns)")
    ax.set_ylabel("Entanglement entropy (bits)")
    ax.set_title(r"Residual entanglement for $|{+}x\rangle \otimes |0\rangle$")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(
        FIGURES_DIR / "unconditional_superposition_entanglement.png",
        dpi=300,
        bbox_inches="tight",
    )
    fig.savefig(
        FIGURES_DIR / "unconditional_superposition_entanglement.pdf",
        bbox_inches="tight",
    )
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.0))
    for family, color in zip(PULSE_FAMILIES, TOL_BRIGHT, strict=False):
        rows_nf = [
            row
            for row in single_tone_data["baseline"][family]["nofilter"]
            if abs(row["alpha_target"] - 1.0) < 1e-12
        ]
        rows_lp = [
            row
            for row in single_tone_data["baseline"][family]["lp_40mhz"]
            if abs(row["alpha_target"] - 1.0) < 1e-12
        ]
        xs = np.array([row["duration_ns"] for row in rows_nf], dtype=float)
        order = np.argsort(xs)
        axes[0].plot(
            xs[order],
            np.array([row["delta_alpha"] for row in rows_nf])[order],
            "o-",
            color=color,
            label=f"{family.title()} no filter",
        )
        axes[0].plot(
            xs[order],
            np.array([row["delta_alpha"] for row in rows_lp])[order],
            "--",
            color=color,
            alpha=0.6,
            label=f"{family.title()} + 40 MHz LPF",
        )
        axes[1].plot(
            xs[order],
            np.array([row["bandwidth_mhz"] for row in rows_nf])[order],
            "o-",
            color=color,
            label=f"{family.title()}",
        )
    axes[0].set_xlabel("Duration (ns)")
    axes[0].set_ylabel(r"$\delta \alpha$")
    axes[0].set_title("Conditionality vs duration")
    axes[1].set_xlabel("Duration (ns)")
    axes[1].set_ylabel("95% bandwidth (MHz)")
    axes[1].set_title("Command bandwidth")
    axes[0].legend(fontsize=7, ncol=2)
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "unconditional_filter_tradeoff.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "unconditional_filter_tradeoff.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    for chi_scale, color in zip(CHI_SCALES, TOL_BRIGHT[: len(CHI_SCALES)], strict=False):
        rows = [
            row
            for row in single_tone_data["chi_scaling"]
            if abs(row["chi_scale"] - chi_scale) < 1e-12
        ]
        xs = np.array([row["duration_ns"] for row in rows], dtype=float)
        ys = np.array([row["delta_alpha"] for row in rows], dtype=float)
        order = np.argsort(xs)
        ax.plot(xs[order], ys[order], "o-", color=color, label=fr"$|\chi|$ scale = {chi_scale:.1f}")
    ax.set_xlabel("Duration (ns)")
    ax.set_ylabel(r"$\delta \alpha$")
    ax.set_title(r"Naive displacement error scaling with $|\chi|$")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "unconditional_chi_scaling.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "unconditional_chi_scaling.pdf", bbox_inches="tight")
    plt.close(fig)

    labels = [item["label"] for item in protocol_payload["protocols"]]
    mean_fid = [item["full_metrics"]["state_test_mean_fidelity"] for item in protocol_payload["protocols"]]
    delta_alpha = [item["vacuum_metrics"]["delta_alpha"] for item in protocol_payload["protocols"]]
    ent = [item["vacuum_metrics"]["plus_x_entanglement_bits"] for item in protocol_payload["protocols"]]
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 4.3))
    x = np.arange(len(labels))
    axes[0].bar(x, mean_fid, color=TOL_BRIGHT[: len(labels)])
    axes[0].set_title("State-test mean fidelity")
    axes[0].set_ylabel("Average fidelity")
    axes[1].bar(x, delta_alpha, color=TOL_BRIGHT[: len(labels)])
    axes[1].set_title(r"Branch mismatch $\delta \alpha$")
    axes[2].bar(x, ent, color=TOL_BRIGHT[: len(labels)])
    axes[2].set_title(r"Residual entanglement ($|{+}x\rangle$)")
    for ax in axes:
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=25, ha="right")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "unconditional_protocol_summary.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "unconditional_protocol_summary.pdf", bbox_inches="tight")
    plt.close(fig)

    if wigner_data:
        fig, axes = plt.subplots(2, 2, figsize=(9.0, 7.0))
        for col, (key, title) in enumerate((("naive_fail", "Naive"), ("two_tone_success", "Two-tone"))):
            payload = wigner_data[key]
            for row, field in enumerate(("actual_wigner", "ideal_wigner")):
                ax = axes[row, col]
                im = ax.imshow(
                    np.asarray(payload[field], dtype=float),
                    extent=[-WIGNER_XMAX, WIGNER_XMAX, -WIGNER_XMAX, WIGNER_XMAX],
                    origin="lower",
                    cmap="RdBu_r",
                )
                ax.set_xlabel(r"$x$")
                ax.set_ylabel(r"$p$")
                ax.set_title(f"{title} {'actual' if row == 0 else 'ideal'}")
                plt.colorbar(im, ax=ax, fraction=0.046)
        fig.tight_layout()
        fig.savefig(FIGURES_DIR / "unconditional_wigner_comparison.png", dpi=300, bbox_inches="tight")
        fig.savefig(FIGURES_DIR / "unconditional_wigner_comparison.pdf", bbox_inches="tight")
        plt.close(fig)


def compile_summary(
    single_tone_data: dict[str, Any],
    two_tone_data: dict[str, Any],
    echo_data: dict[str, Any],
    optimal_data: dict[str, Any],
    protocol_payload: dict[str, Any],
) -> dict[str, Any]:
    square_rows = [
        row for row in single_tone_data["baseline"]["square"]["nofilter"]
        if abs(row["alpha_target"] - 1.0) < 1e-12
    ]
    fast_gaussian_rows = [
        row for row in single_tone_data["baseline"]["gaussian"]["nofilter"]
        if abs(row["alpha_target"] - 1.0) < 1e-12
    ]
    best_fast = max(fast_gaussian_rows, key=lambda item: item["plus_x_fidelity"])
    best_two_tone = max(
        [row for row in two_tone_data["calibrated_sweep"] if abs(row["alpha_target"] - 1.0) < 1e-12],
        key=lambda item: item["plus_x_fidelity"],
    )
    best_echo = max(
        [row for row in echo_data["echo_sweep"] if abs(row["alpha_target"] - 1.0) < 1e-12],
        key=lambda item: item["plus_x_fidelity"],
    )
    best_optimal = max(
        optimal_data["cases"],
        key=lambda item: item["full_metrics"]["state_test_mean_fidelity"],
    )
    best_protocol = max(
        protocol_payload["protocols"],
        key=lambda item: item["full_metrics"]["state_test_mean_fidelity"],
    )
    return {
        "baseline_square_80ns": next(
            row for row in square_rows if abs(row["duration_ns"] - 80.0) < 1e-12
        ),
        "best_fast_gaussian": best_fast,
        "best_two_tone": best_two_tone,
        "best_echo": best_echo,
        "best_optimal": best_optimal,
        "best_protocol_snapshot": best_protocol,
        "chi_guideline_ns": float(1.0 / abs(CHI) * 1.0e9),
    }


def main() -> None:
    started = time.time()
    single_tone_data = run_single_tone_sweeps()
    two_tone_data = run_two_tone_sweep()
    echo_data = run_echo_sweep()
    optimal_data = run_optimal_control_benchmark()
    protocol_payload = run_selected_protocol_comparison(optimal_data)
    wigner_data = build_wigner_data(protocol_payload)
    generate_figures(single_tone_data, two_tone_data, echo_data, protocol_payload, wigner_data)
    summary = compile_summary(
        single_tone_data, two_tone_data, echo_data, optimal_data, protocol_payload
    )
    summary["wall_time_s"] = time.time() - started
    save_json(ARTIFACTS_DIR / "unconditional_displacement_summary.json", summary)
    print("=" * 72)
    print("Unconditional displacement study complete")
    print("=" * 72)
    print(f"Total wall time: {summary['wall_time_s']:.1f} s")
    print(f"Best protocol: {summary['best_protocol_snapshot']['label']}")


if __name__ == "__main__":
    main()
