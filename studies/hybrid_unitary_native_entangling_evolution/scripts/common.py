"""Shared definitions for the native-entangling hybrid-unitary study."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import os
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence
import json
import re

import numpy as np

WORKSPACE_ROOT = Path(__file__).resolve().parents[3]
STUDY_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = STUDY_ROOT / "data"
FIG_DIR = STUDY_ROOT / "figures"
REPORT_DIR = STUDY_ROOT / "report"
ARTIFACT_DIR = STUDY_ROOT / "artifacts"
LEGACY_STUDY_ROOT = WORKSPACE_ROOT / "studies" / "hybrid_qubit_cavity_control"
STYLE_PATH = (
    WORKSPACE_ROOT
    / ".github"
    / "skills"
    / "publication-figures"
    / "assets"
    / "cqed_style.mplstyle"
)

LOGICAL_BASIS_LABELS = ["|g,0>", "|g,1>", "|e,0>", "|e,1>"]


def resolve_sim_root() -> Path:
    """Locate the sibling cqed_sim checkout used by this workspace."""
    raw_env = os.environ.get("CQED_SIM_ROOT")
    candidates: list[Path] = []
    if raw_env:
        candidates.append(Path(raw_env))
    candidates.append(WORKSPACE_ROOT.parent / "cQED_simulation")
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError(
        "Could not locate cQED_simulation. Set CQED_SIM_ROOT or place the checkout beside the workspace."
    )


SIM_ROOT = resolve_sim_root()


def ensure_sim_root_on_path() -> None:
    """Ensure the local cqed_sim checkout is importable."""
    if str(SIM_ROOT) not in sys.path:
        sys.path.insert(0, str(SIM_ROOT))

OMEGA_Q = 2.0 * np.pi * 6.150e9
OMEGA_C = 2.0 * np.pi * 5.241e9
ALPHA = 2.0 * np.pi * (-255e6)
CHI = 2.0 * np.pi * (-2.84e6)
CHIP = 2.0 * np.pi * (-21e3)
KERR = 2.0 * np.pi * (-28e3)
N_TR = 2
NATIVE_WAIT_DURATION_NS = 256.0

ENTANGLING_GATE_TYPES = {
    "SQR",
    "ConditionalPhaseSQR",
    "FreeEvolveCondPhase",
    "ConditionalDisplacement",
    "JaynesCummingsExchange",
    "BlueSidebandExchange",
}

LOCAL_GATE_TYPES = {
    "QubitRotation",
    "Displacement",
    "SNAP",
    "CavityBlockPhase",
    "PrimitiveGate",
}


@dataclass(frozen=True)
class CostWeights:
    """Default weights for the entangling-aware ranking metric."""

    infidelity: float = 2.0
    entangling_gate_count: float = 0.60
    entangling_time: float = 0.40
    depth: float = 0.05
    implementation_complexity: float = 0.05
    leakage: float = 0.50


def logical_subspace_indices(n_cav: int) -> tuple[int, int, int, int]:
    """Return the logical {|g,0>, |g,1>, |e,0>, |e,1>} indices."""
    n_cav = int(n_cav)
    return (0, 1, n_cav, n_cav + 1)


def target_unitary_matrix() -> np.ndarray:
    """Return the 4x4 logical target unitary used in the hybrid studies."""
    scale = 1.0 / np.sqrt(2.0)
    return np.array(
        [
            [scale, 0.0, scale, 0.0],
            [scale, 0.0, -scale, 0.0],
            [0.0, scale, 0.0, scale],
            [0.0, -scale, 0.0, scale],
        ],
        dtype=np.complex128,
    )


def build_model(*, n_cav: int = 8, n_tr: int = N_TR):
    """Build the dispersive transmon-cavity model used across the hybrid studies."""
    from cqed_sim.core import DispersiveTransmonCavityModel

    return DispersiveTransmonCavityModel(
        omega_c=OMEGA_C,
        omega_q=OMEGA_Q,
        alpha=ALPHA,
        chi=CHI,
        chi_higher=(CHIP,),
        kerr=KERR,
        n_cav=int(n_cav),
        n_tr=int(n_tr),
    )


def build_frame(model: Any):
    """Return the natural rotating frame for the hybrid-study device point."""
    from cqed_sim.core import FrameSpec

    return FrameSpec(omega_c_frame=float(model.omega_c), omega_q_frame=float(model.omega_q))


def embed_logical_state(logical_vector: Sequence[complex], *, n_cav: int = 8):
    """Embed a 4D logical state into the full qubit-cavity Hilbert space."""
    import qutip as qt

    vector = np.asarray(logical_vector, dtype=np.complex128).reshape(-1)
    if vector.shape != (4,):
        raise ValueError(f"logical_vector must have shape (4,), got {vector.shape}")
    full = np.zeros(2 * int(n_cav), dtype=np.complex128)
    full[0] = vector[0]
    full[1] = vector[1]
    full[int(n_cav)] = vector[2]
    full[int(n_cav) + 1] = vector[3]
    return qt.Qobj(full, dims=[[2, int(n_cav)], [1, 1]])


def average_gate_fidelity_from_process(f_process: float, dim: int) -> float:
    """Convert process fidelity into average gate fidelity."""
    return float((dim * float(f_process) + 1.0) / (dim + 1.0))


def dump_json(path: Path, payload: Mapping[str, Any] | Sequence[Any]) -> None:
    """Write a JSON payload with stable formatting."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=False), encoding="utf-8")


def load_json(path: Path) -> Any:
    """Load JSON from disk."""
    return json.loads(path.read_text(encoding="utf-8"))


def apply_publication_style() -> None:
    """Apply the shared study plotting style when available."""
    import matplotlib.pyplot as plt

    if STYLE_PATH.exists():
        plt.style.use(str(STYLE_PATH))


def normalize_label(raw_label: str) -> str:
    """Map legacy study labels onto stable short identifiers."""
    label = raw_label.strip()
    lowered = label.lower()
    canonical = {
        "a_ent": "A_ent",
        "a_local": "A_local",
        "b_ent": "B_ent",
        "b_local": "B_local",
        "l1d": "L1d",
        "l2d": "L2d",
        "l2c": "L2c",
    }
    for prefix, normalized in canonical.items():
        if lowered.startswith(prefix):
            return normalized
    token = re.split(r"[\s(:-]", label, maxsplit=1)[0]
    return token.strip() or label


def is_entangling_gate(gate_type: str) -> bool:
    """Return whether a gate type should count as an expensive entangler."""
    return gate_type in ENTANGLING_GATE_TYPES


def count_gate_types(sequence: Sequence[Mapping[str, Any]]) -> dict[str, int]:
    """Count gate types in a serialized sequence."""
    counts: dict[str, int] = {}
    for gate in sequence:
        gate_type = str(gate.get("type", "UNKNOWN"))
        counts[gate_type] = counts.get(gate_type, 0) + 1
    return counts


def entangling_time_from_sequence(sequence: Sequence[Mapping[str, Any]]) -> float:
    """Compute entangling time in ns from a serialized sequence."""
    total = 0.0
    for gate in sequence:
        if is_entangling_gate(str(gate.get("type", ""))):
            total += 1.0e9 * float(gate.get("duration", 0.0))
    return float(total)


def local_gate_count_from_sequence(sequence: Sequence[Mapping[str, Any]]) -> int:
    """Count the cheap local gates in a serialized sequence."""
    return sum(1 for gate in sequence if str(gate.get("type", "")) in LOCAL_GATE_TYPES)


def infer_family(sequence: Sequence[Mapping[str, Any]] | None, fallback: str | None = None) -> str:
    """Infer a short family label from a sequence when the source does not provide one."""
    if fallback:
        return fallback
    if not sequence:
        return "unknown"
    gate_types = {str(gate.get("type", "")) for gate in sequence}
    if "FreeEvolveCondPhase" in gate_types:
        return "native"
    if "ConditionalPhaseSQR" in gate_types:
        return "CPSQR"
    if "SNAP" in gate_types:
        return "SNAP"
    if "SQR" in gate_types:
        return "SQR"
    return "mixed"


def implementation_complexity_score(
    *,
    sequence: Sequence[Mapping[str, Any]] | None = None,
    active_tones: int | float = 0,
    gate_types: Sequence[str] | None = None,
) -> float:
    """Build a simple calibration-complexity score for cross-candidate ranking."""
    unique_types = {str(gate_type) for gate_type in (gate_types or [])}
    selective_instances = 0
    if sequence:
        unique_types.update(str(gate.get("type", "")) for gate in sequence)
        selective_instances = sum(
            1
            for gate in sequence
            if str(gate.get("type", "")) in {"SQR", "ConditionalPhaseSQR", "SNAP", "FreeEvolveCondPhase"}
        )
    else:
        selective_instances = sum(
            1
            for gate_type in unique_types
            if gate_type in {"SQR", "ConditionalPhaseSQR", "SNAP", "FreeEvolveCondPhase"}
        )
    return float(active_tones) + 0.50 * len(unique_types) + 0.25 * selective_instances


def best_available_fidelity(candidate: Mapping[str, Any]) -> tuple[float, str]:
    """Return the best available comparison fidelity and its provenance."""
    ordered_keys = (
        ("pulse_fidelity", "pulse"),
        ("compiled_estimated_fidelity", "compiled-estimate"),
        ("process_fidelity", "ideal-restart"),
        ("ideal_fidelity", "ideal-frontier"),
    )
    for key, level in ordered_keys:
        value = candidate.get(key)
        if value is not None:
            return float(value), level
    return 0.0, "missing"


def preferred_leakage(candidate: Mapping[str, Any]) -> float:
    """Return the most physical leakage value available for a candidate."""
    for key in ("pulse_leakage_average", "leakage_average", "ideal_leakage"):
        value = candidate.get(key)
        if value is not None:
            return float(value)
    return 0.0


def candidate_weighted_cost(candidate: Mapping[str, Any], weights: CostWeights | None = None) -> float:
    """Compute the entangling-aware ranking metric for a candidate."""
    weights = weights or CostWeights()
    fidelity, _ = best_available_fidelity(candidate)
    entangling_count = float(candidate.get("entangling_gate_count", 0.0) or 0.0)
    entangling_time_ns = float(candidate.get("entangling_time_ns", 0.0) or 0.0)
    depth = float(candidate.get("depth", 0.0) or 0.0)
    complexity = float(candidate.get("implementation_complexity", 0.0) or 0.0)
    leakage = preferred_leakage(candidate)
    return float(
        weights.infidelity * (1.0 - fidelity)
        + weights.entangling_gate_count * entangling_count
        + weights.entangling_time * (entangling_time_ns / NATIVE_WAIT_DURATION_NS)
        + weights.depth * depth
        + weights.implementation_complexity * complexity
        + weights.leakage * leakage
    )


def as_plain_dict(weights: CostWeights) -> dict[str, float]:
    """Expose dataclass weights in a JSON-friendly form."""
    return asdict(weights)