"""Study-local readout optimization helpers."""

from .bounds import (
    LinearResponse,
    assignment_fidelity_from_snr2,
    classification_probability,
    matched_filter_snr2,
    solve_linear_response,
    solve_two_segment_nulling,
    t1_limited_assignment_bound,
)
from .config import DEFAULT_CONFIG, ReadoutStudyConfig
from .metrics import ReadoutMetrics, build_metrics, repeated_consistency_from_signals
from .pulse_families import (
    FAMILIES,
    PulseDesign,
    build_piecewise_reference,
    get_family,
    set_nulling_tail_kappa,
)

__all__ = [
    "DEFAULT_CONFIG",
    "FAMILIES",
    "LinearResponse",
    "PulseDesign",
    "ReadoutMetrics",
    "ReadoutStudyConfig",
    "assignment_fidelity_from_snr2",
    "build_metrics",
    "build_piecewise_reference",
    "classification_probability",
    "get_family",
    "matched_filter_snr2",
    "repeated_consistency_from_signals",
    "set_nulling_tail_kappa",
    "solve_linear_response",
    "solve_two_segment_nulling",
    "t1_limited_assignment_bound",
]
