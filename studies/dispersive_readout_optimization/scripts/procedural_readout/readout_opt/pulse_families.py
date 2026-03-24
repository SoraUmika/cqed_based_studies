"""Pulse-family parameterizations for the procedural readout study."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from .bounds import propagate_piecewise_constant, solve_two_segment_nulling


@dataclass(frozen=True)
class PulseDesign:
    """A concrete pulse waveform plus its drive detuning choice."""

    family: str
    params: np.ndarray
    waveform: np.ndarray
    dt: float
    duration: float
    delta_g: float
    metadata: dict[str, float | str | bool | list[float]]


@dataclass(frozen=True)
class PulseFamily:
    """Registry entry for a pulse family."""

    name: str
    dimension: int
    bounds: tuple[tuple[float, float], ...]
    builder: Callable[[np.ndarray, float, float, float, float], PulseDesign]


def midpoint_offset_to_delta_g(offset: float, chi: float) -> float:
    """Map an offset parameter onto the drive detuning from the ground-state cavity."""
    return float(-0.5 * chi + offset * abs(chi))


def duration_weights_to_counts(
    weights: np.ndarray,
    *,
    total_count: int,
    min_count: int,
) -> np.ndarray:
    """Convert positive weights to integer segment sample counts."""
    raw = np.maximum(np.asarray(weights, dtype=float), 1.0e-9)
    n_seg = len(raw)
    if total_count < n_seg * min_count:
        raise ValueError("total_count is too small for the requested minimum segment count.")
    free = total_count - n_seg * min_count
    counts = min_count + np.floor(free * raw / np.sum(raw)).astype(int)
    while int(np.sum(counts)) < total_count:
        idx = int(np.argmax(raw / np.maximum(counts, 1)))
        counts[idx] += 1
    while int(np.sum(counts)) > total_count:
        idx = int(np.argmax(counts))
        if counts[idx] > min_count:
            counts[idx] -= 1
        else:
            break
    return counts


def clip_waveform_amplitude(waveform: np.ndarray, amp_max: float) -> np.ndarray:
    """Pointwise complex-amplitude clipping."""
    wf = np.asarray(waveform, dtype=np.complex128).copy()
    scale = np.maximum(1.0, np.abs(wf) / float(amp_max))
    return wf / scale


def smooth_square_envelope(sample_count: int, rise_fraction: float) -> np.ndarray:
    """Return a unit-amplitude smooth square with cosine edges."""
    n = int(sample_count)
    rise_fraction = float(np.clip(rise_fraction, 0.02, 0.48))
    t_rel = np.linspace(0.0, 1.0, n, endpoint=False)
    env = np.ones(n, dtype=float)
    left = t_rel < rise_fraction
    right = t_rel > 1.0 - rise_fraction
    env[left] = 0.5 * (1.0 - np.cos(np.pi * t_rel[left] / rise_fraction))
    tail = (1.0 - t_rel[right]) / rise_fraction
    env[right] = 0.5 * (1.0 - np.cos(np.pi * np.clip(tail, 0.0, 1.0)))
    return env.astype(np.complex128)


def build_square(params: np.ndarray, duration: float, dt: float, amp_max: float, chi: float) -> PulseDesign:
    n_steps = int(round(duration / dt))
    amp = float(np.clip(params[0], 0.0, 1.0)) * amp_max
    phase = float(params[1])
    delta_g = midpoint_offset_to_delta_g(float(params[2]), chi)
    waveform = np.full(n_steps, amp * np.exp(1j * phase), dtype=np.complex128)
    return PulseDesign(
        family="square",
        params=np.asarray(params, dtype=float),
        waveform=waveform,
        dt=float(dt),
        duration=float(duration),
        delta_g=delta_g,
        metadata={"phase": phase, "amp": amp},
    )


def build_smooth_square(params: np.ndarray, duration: float, dt: float, amp_max: float, chi: float) -> PulseDesign:
    n_steps = int(round(duration / dt))
    amp = float(np.clip(params[0], 0.0, 1.0)) * amp_max
    phase = float(params[1])
    rise_fraction = float(np.clip(params[2], 0.04, 0.45))
    delta_g = midpoint_offset_to_delta_g(float(params[3]), chi)
    waveform = amp * np.exp(1j * phase) * smooth_square_envelope(n_steps, rise_fraction)
    return PulseDesign(
        family="smooth_square",
        params=np.asarray(params, dtype=float),
        waveform=waveform,
        dt=float(dt),
        duration=float(duration),
        delta_g=delta_g,
        metadata={"phase": phase, "amp": amp, "rise_fraction": rise_fraction},
    )


def build_ring_hold(params: np.ndarray, duration: float, dt: float, amp_max: float, chi: float) -> PulseDesign:
    n_steps = int(round(duration / dt))
    counts = duration_weights_to_counts(params[3:6], total_count=n_steps, min_count=max(1, int(0.06 * n_steps)))
    amps = amp_max * np.clip(np.asarray(params[:3], dtype=float), 0.0, 1.0)
    waveform = np.concatenate(
        [np.full(count, amp, dtype=np.complex128) for amp, count in zip(amps, counts, strict=True)]
    )
    waveform = clip_waveform_amplitude(waveform, amp_max)
    delta_g = midpoint_offset_to_delta_g(float(params[6]), chi)
    return PulseDesign(
        family="ring_hold",
        params=np.asarray(params, dtype=float),
        waveform=waveform,
        dt=float(dt),
        duration=float(duration),
        delta_g=delta_g,
        metadata={
            "segment_counts": counts.tolist(),
            "amps": amps.tolist(),
        },
    )


def build_procedural_segments(
    params: np.ndarray,
    duration: float,
    dt: float,
    amp_max: float,
    chi: float,
) -> PulseDesign:
    n_steps = int(round(duration / dt))
    counts = duration_weights_to_counts(params[6:9], total_count=n_steps, min_count=max(1, int(0.06 * n_steps)))
    segments = amp_max * np.array(
        [
            complex(params[0], params[1]),
            complex(params[2], params[3]),
            complex(params[4], params[5]),
        ],
        dtype=np.complex128,
    )
    waveform = np.concatenate(
        [np.full(count, amp, dtype=np.complex128) for amp, count in zip(segments, counts, strict=True)]
    )
    waveform = clip_waveform_amplitude(waveform, amp_max)
    delta_g = midpoint_offset_to_delta_g(float(params[9]), chi)
    return PulseDesign(
        family="procedural_segments",
        params=np.asarray(params, dtype=float),
        waveform=waveform,
        dt=float(dt),
        duration=float(duration),
        delta_g=delta_g,
        metadata={
            "segment_counts": counts.tolist(),
            "segments_real_imag": [float(v) for z in segments for v in (z.real, z.imag)],
        },
    )


def build_nulling_tail(params: np.ndarray, duration: float, dt: float, amp_max: float, chi: float) -> PulseDesign:
    kappa_value = getattr(build_nulling_tail, "_kappa", None)
    if kappa_value is None:
        raise RuntimeError("Nulling-tail builder requires set_nulling_tail_kappa(kappa) before use.")
    n_steps = int(round(duration / dt))
    counts = duration_weights_to_counts(params[4:8], total_count=n_steps, min_count=max(1, int(0.05 * n_steps)))
    delta_g = midpoint_offset_to_delta_g(float(params[8]), chi)
    durations = counts * float(dt)

    seg1 = amp_max * complex(params[0], params[1])
    seg2 = amp_max * complex(params[2], params[3])
    alpha_g = propagate_piecewise_constant(
        np.array([seg1, seg2]),
        durations[:2],
        kappa=float(kappa_value),
        delta=delta_g,
    )
    alpha_e = propagate_piecewise_constant(
        np.array([seg1, seg2]),
        durations[:2],
        kappa=float(kappa_value),
        delta=delta_g + chi,
    )
    tail_1, tail_2 = solve_two_segment_nulling(
        alpha_g[-1],
        alpha_e[-1],
        tau_3=durations[2],
        tau_4=durations[3],
        kappa=float(kappa_value),
        chi=chi,
        delta_g=delta_g,
    )
    waveform = np.concatenate(
        [
            np.full(counts[0], seg1, dtype=np.complex128),
            np.full(counts[1], seg2, dtype=np.complex128),
            np.full(counts[2], tail_1, dtype=np.complex128),
            np.full(counts[3], tail_2, dtype=np.complex128),
        ]
    )
    clipped = np.any(np.abs(waveform) > amp_max + 1.0e-12)
    waveform = clip_waveform_amplitude(waveform, amp_max)
    return PulseDesign(
        family="nulling_tail",
        params=np.asarray(params, dtype=float),
        waveform=waveform,
        dt=float(dt),
        duration=float(duration),
        delta_g=delta_g,
        metadata={
            "segment_counts": counts.tolist(),
            "tail_clipped": bool(clipped),
            "tail_1_real": float(np.real(tail_1)),
            "tail_1_imag": float(np.imag(tail_1)),
            "tail_2_real": float(np.real(tail_2)),
            "tail_2_imag": float(np.imag(tail_2)),
        },
    )


def build_fourier_basis(params: np.ndarray, duration: float, dt: float, amp_max: float, chi: float) -> PulseDesign:
    n_steps = int(round(duration / dt))
    t_rel = np.linspace(0.0, 1.0, n_steps, endpoint=False)
    c0 = complex(params[0], params[1])
    c1 = complex(params[2], params[3])
    c2 = complex(params[4], params[5])
    waveform = c0 + c1 * np.sin(np.pi * t_rel) + c2 * np.sin(2.0 * np.pi * t_rel)
    norm = max(np.max(np.abs(waveform)), 1.0e-12)
    waveform = amp_max * waveform / norm
    delta_g = midpoint_offset_to_delta_g(float(params[6]), chi)
    return PulseDesign(
        family="fourier_basis",
        params=np.asarray(params, dtype=float),
        waveform=clip_waveform_amplitude(waveform, amp_max),
        dt=float(dt),
        duration=float(duration),
        delta_g=delta_g,
        metadata={"norm": float(norm)},
    )


def build_piecewise_reference(
    params: np.ndarray,
    duration: float,
    dt: float,
    amp_max: float,
    chi: float,
    *,
    n_segments: int = 8,
) -> PulseDesign:
    n_steps = int(round(duration / dt))
    weights = np.ones(n_segments, dtype=float)
    counts = duration_weights_to_counts(weights, total_count=n_steps, min_count=max(1, int(0.04 * n_steps)))
    complex_segments = []
    for idx in range(n_segments):
        re = params[2 * idx]
        im = params[2 * idx + 1]
        complex_segments.append(amp_max * complex(re, im))
    waveform = np.concatenate(
        [np.full(count, amp, dtype=np.complex128) for amp, count in zip(complex_segments, counts, strict=True)]
    )
    waveform = clip_waveform_amplitude(waveform, amp_max)
    delta_g = midpoint_offset_to_delta_g(float(params[-1]), chi)
    return PulseDesign(
        family="piecewise_reference",
        params=np.asarray(params, dtype=float),
        waveform=waveform,
        dt=float(dt),
        duration=float(duration),
        delta_g=delta_g,
        metadata={"segment_counts": counts.tolist(), "n_segments": float(n_segments)},
    )


def set_nulling_tail_kappa(kappa: float) -> None:
    """Inject the cavity linewidth used by the analytical nulling builder."""
    build_nulling_tail._kappa = float(kappa)  # type: ignore[attr-defined]


FAMILIES: dict[str, PulseFamily] = {
    "square": PulseFamily(
        name="square",
        dimension=3,
        bounds=((0.0, 1.0), (-np.pi, np.pi), (-1.0, 1.0)),
        builder=build_square,
    ),
    "smooth_square": PulseFamily(
        name="smooth_square",
        dimension=4,
        bounds=((0.0, 1.0), (-np.pi, np.pi), (0.04, 0.40), (-1.0, 1.0)),
        builder=build_smooth_square,
    ),
    "ring_hold": PulseFamily(
        name="ring_hold",
        dimension=7,
        bounds=((0.0, 1.0), (0.0, 1.0), (0.0, 1.0), (0.1, 1.0), (0.1, 1.0), (0.1, 1.0), (-1.0, 1.0)),
        builder=build_ring_hold,
    ),
    "procedural_segments": PulseFamily(
        name="procedural_segments",
        dimension=10,
        bounds=(
            (-1.0, 1.0), (-1.0, 1.0),
            (-1.0, 1.0), (-1.0, 1.0),
            (-1.0, 1.0), (-1.0, 1.0),
            (0.1, 1.0), (0.1, 1.0), (0.1, 1.0),
            (-1.0, 1.0),
        ),
        builder=build_procedural_segments,
    ),
    "nulling_tail": PulseFamily(
        name="nulling_tail",
        dimension=9,
        bounds=(
            (-1.0, 1.0), (-1.0, 1.0),
            (-1.0, 1.0), (-1.0, 1.0),
            (0.1, 1.0), (0.1, 1.0), (0.1, 1.0), (0.1, 1.0),
            (-1.0, 1.0),
        ),
        builder=build_nulling_tail,
    ),
    "fourier_basis": PulseFamily(
        name="fourier_basis",
        dimension=7,
        bounds=(
            (-1.0, 1.0), (-1.0, 1.0),
            (-1.0, 1.0), (-1.0, 1.0),
            (-1.0, 1.0), (-1.0, 1.0),
            (-1.0, 1.0),
        ),
        builder=build_fourier_basis,
    ),
}


def get_family(name: str) -> PulseFamily:
    if name == "piecewise_reference":
        raise ValueError("piecewise_reference is configured directly through build_piecewise_reference().")
    return FAMILIES[name]
