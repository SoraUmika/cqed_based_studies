from __future__ import annotations

import json
from pathlib import Path


STUDY_NAME = "realistic_universal_hybrid_control_dispersive_cqed"
PROJECT_ROOT = Path(__file__).resolve().parents[3]
STUDY_ROOT = PROJECT_ROOT / "studies" / STUDY_NAME
DATA_DIR = STUDY_ROOT / "data"
FIGURES_DIR = STUDY_ROOT / "figures"
ARTIFACTS_DIR = STUDY_ROOT / "artifacts"
REPORT_DIR = STUDY_ROOT / "report"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def save_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def chi_scaled_duration(duration_ns: float, chi_over_2pi_mhz: float) -> float:
    """Return |chi| T / (2 pi) using chi/2pi in MHz and T in ns."""

    return abs(chi_over_2pi_mhz) * duration_ns * 1e-3
