"""Optimization wrappers for the procedural readout pulse families."""

from __future__ import annotations

from dataclasses import dataclass
from functools import partial
from typing import Callable

import numpy as np
from scipy.optimize import minimize

from .config import ReadoutStudyConfig
from .pulse_families import (
    PulseDesign,
    build_piecewise_reference,
    get_family,
    set_nulling_tail_kappa,
)
from .simulate import ReplayEvaluation, evaluate_full_design, evaluate_linear_design


@dataclass(frozen=True)
class OptimizationOutcome:
    """Best design and evaluation found for a given family/duration/regime."""

    family: str
    regime: str
    objective_name: str
    duration: float
    x_best: np.ndarray
    design: PulseDesign
    evaluation: ReplayEvaluation
    history: list[dict[str, float]]


def objective_value(evaluation: ReplayEvaluation, objective_name: str) -> float:
    metrics = evaluation.metrics
    if objective_name == "balanced":
        return float(metrics.score_balanced)
    if objective_name == "info":
        return float(metrics.score_info)
    if objective_name == "emptying":
        return float(metrics.score_emptying)
    raise ValueError(f"Unsupported objective '{objective_name}'.")


def _builder_for_family(
    family: str,
    cfg: ReadoutStudyConfig,
    *,
    n_reference_segments: int = 8,
) -> tuple[Callable[[np.ndarray, float], PulseDesign], tuple[tuple[float, float], ...]]:
    set_nulling_tail_kappa(cfg.kappa)
    if family == "piecewise_reference":
        bounds = tuple([(-1.0, 1.0)] * (2 * n_reference_segments) + [(-1.0, 1.0)])
        return (
            partial(
                build_piecewise_reference,
                dt=cfg.dt,
                amp_max=cfg.amp_max,
                chi=cfg.chi,
                n_segments=n_reference_segments,
            ),
            bounds,
        )
    family_spec = get_family(family)
    return (
        partial(family_spec.builder, dt=cfg.dt, amp_max=cfg.amp_max, chi=cfg.chi),
        family_spec.bounds,
    )


def _evaluator_for_regime(
    regime: str,
    cfg: ReadoutStudyConfig,
) -> Callable[[PulseDesign], ReplayEvaluation]:
    if regime == "linear":
        return partial(evaluate_linear_design, cfg=cfg)
    if regime == "full":
        return partial(evaluate_full_design, cfg=cfg)
    raise ValueError(f"Unsupported regime '{regime}'.")


def _clip_to_bounds(x: np.ndarray, bounds: tuple[tuple[float, float], ...]) -> np.ndarray:
    arr = np.asarray(x, dtype=float).copy()
    for idx, (low, high) in enumerate(bounds):
        arr[idx] = np.clip(arr[idx], low, high)
    return arr


def _random_points(
    bounds: tuple[tuple[float, float], ...],
    *,
    rng: np.random.Generator,
    count: int,
) -> np.ndarray:
    lows = np.array([b[0] for b in bounds], dtype=float)
    highs = np.array([b[1] for b in bounds], dtype=float)
    return lows + (highs - lows) * rng.random((count, len(bounds)))


def optimize_family(
    *,
    family: str,
    regime: str,
    objective_name: str,
    duration: float,
    cfg: ReadoutStudyConfig,
    warm_start: np.ndarray | None = None,
    seed: int | None = None,
    n_random: int = 18,
    n_local: int = 3,
    maxiter: int = 48,
) -> OptimizationOutcome:
    """Optimize a family for the requested regime and scalarized objective."""
    builder, bounds = _builder_for_family(family, cfg)
    evaluator = _evaluator_for_regime(regime, cfg)
    rng = np.random.default_rng(cfg.seed if seed is None else seed)
    history: list[dict[str, float]] = []
    cache: dict[tuple[float, ...], tuple[PulseDesign, ReplayEvaluation]] = {}

    def evaluate_params(x: np.ndarray) -> tuple[PulseDesign, ReplayEvaluation]:
        x_clipped = _clip_to_bounds(x, bounds)
        key = tuple(np.round(x_clipped, decimals=9))
        cached = cache.get(key)
        if cached is not None:
            return cached
        design = builder(x_clipped, duration)
        evaluation = evaluator(design)
        cache[key] = (design, evaluation)
        history.append(
            {
                "score_balanced": float(evaluation.metrics.score_balanced),
                "score_info": float(evaluation.metrics.score_info),
                "score_emptying": float(evaluation.metrics.score_emptying),
            }
        )
        return design, evaluation

    def objective_fn(x: np.ndarray) -> float:
        _, evaluation = evaluate_params(x)
        return -objective_value(evaluation, objective_name)

    samples = _random_points(bounds, rng=rng, count=n_random)
    if warm_start is not None:
        samples = np.vstack([_clip_to_bounds(np.asarray(warm_start, dtype=float), bounds), samples])

    scored_samples: list[tuple[float, np.ndarray, PulseDesign, ReplayEvaluation]] = []
    for sample in samples:
        design, evaluation = evaluate_params(sample)
        score = objective_value(evaluation, objective_name)
        scored_samples.append((score, np.asarray(sample, dtype=float), design, evaluation))
    scored_samples.sort(key=lambda item: item[0], reverse=True)

    best_score, best_x, best_design, best_evaluation = scored_samples[0]
    for _, candidate_x, _, _ in scored_samples[: max(1, n_local)]:
        result = minimize(
            objective_fn,
            x0=candidate_x,
            method="Powell",
            bounds=bounds,
            options={"maxiter": int(maxiter), "xtol": 1.0e-3, "ftol": 1.0e-3},
        )
        design, evaluation = evaluate_params(result.x)
        score = objective_value(evaluation, objective_name)
        if score > best_score:
            best_score = score
            best_x = _clip_to_bounds(result.x, bounds)
            best_design = design
            best_evaluation = evaluation

    return OptimizationOutcome(
        family=family,
        regime=regime,
        objective_name=objective_name,
        duration=float(duration),
        x_best=np.asarray(best_x, dtype=float),
        design=best_design,
        evaluation=best_evaluation,
        history=history,
    )
