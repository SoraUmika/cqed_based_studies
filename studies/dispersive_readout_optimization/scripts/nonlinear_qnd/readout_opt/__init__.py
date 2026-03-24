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
from .config import DEFAULT_CONFIG, EffectiveMixingConfig, HardwareProfile, ReadoutStudyConfig
from .metrics import ReadoutMetrics, build_metrics, repeated_consistency_from_signals
from .pulse_families import (
    FAMILIES,
    PulseDesign,
    build_piecewise_reference,
    get_family,
    set_nulling_tail_kappa,
)
from .simulate import (
    evaluate_full_design,
    evaluate_hardware_design,
    evaluate_linear_design,
    evaluate_multilevel_design,
    evaluate_nonlinear_design,
    evaluate_rich_design,
    transport_analysis,
)

__all__ = [
    "DEFAULT_CONFIG",
    "EffectiveMixingConfig",
    "FAMILIES",
    "HardwareProfile",
    "LinearResponse",
    "PulseDesign",
    "ReadoutMetrics",
    "ReadoutStudyConfig",
    "assignment_fidelity_from_snr2",
    "build_metrics",
    "build_piecewise_reference",
    "classification_probability",
    "evaluate_full_design",
    "evaluate_hardware_design",
    "evaluate_linear_design",
    "evaluate_multilevel_design",
    "evaluate_nonlinear_design",
    "evaluate_rich_design",
    "get_family",
    "matched_filter_snr2",
    "repeated_consistency_from_signals",
    "set_nulling_tail_kappa",
    "solve_linear_response",
    "solve_two_segment_nulling",
    "t1_limited_assignment_bound",
    "transport_analysis",
]
