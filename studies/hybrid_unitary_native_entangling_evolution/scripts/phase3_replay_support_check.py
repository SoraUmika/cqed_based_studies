"""Phase 3.1: verify replay support for the top native-heavy candidates.

This script does not run the full physical replay. It checks whether the
installed cqed_sim waveform bridge can translate the shortlisted Phase 2
candidate sequences into waveform-backed primitives, which is the prerequisite
for pulse-level replay through the simulator stack.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

import runtime_compat  # noqa: F401

from common import DATA_DIR, LEGACY_STUDY_ROOT, dump_json

sys.stdout.reconfigure(encoding="utf-8")

from common import ensure_sim_root_on_path

ensure_sim_root_on_path()

from cqed_sim.unitary_synthesis.waveform_bridge import waveform_sequence_from_gates  # noqa: E402

import phase2_native_block_search as phase2  # noqa: E402

INPUT_JSON = DATA_DIR / "phase2_native_block_search.json"
OUTPUT_JSON = DATA_DIR / "phase3_replay_support_check.json"

SHORTLIST = [
    "N2_exact_hc_to_exact_hc",
    "N2_exact_hc_to_A_local",
    "N2_A_local_to_A_local",
]


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def check_candidate(label: str) -> dict[str, Any]:
    rows = load_json(INPUT_JSON)["candidates"]
    row = next(candidate for candidate in rows if candidate["label"] == label)
    sequence = phase2.build_candidate_sequence(
        waits=int(row["waits"]),
        inner_kind=str(row["inner_local_kind"]),
        outer_kind=str(row["outer_local_kind"]),
    )
    gate_types = [gate.__class__.__name__ for gate in sequence.gates]
    try:
        converted = waveform_sequence_from_gates(sequence)
        converted_types = [gate.metadata.get("waveform_family", "unknown") for gate in converted.gates]
        return {
            "label": label,
            "supported": True,
            "gate_types": gate_types,
            "converted_gate_count": len(converted.gates),
            "converted_waveform_families": converted_types,
            "message": "waveform bridge succeeded",
        }
    except Exception as exc:
        return {
            "label": label,
            "supported": False,
            "gate_types": gate_types,
            "converted_gate_count": 0,
            "converted_waveform_families": [],
            "message": f"{type(exc).__name__}: {exc}",
        }


def main() -> None:
    rows = [check_candidate(label) for label in SHORTLIST]
    payload = {
        "metadata": {
            "description": "Replay-support check for shortlisted native-heavy candidates.",
            "legacy_study_root": str(LEGACY_STUDY_ROOT),
        },
        "results": rows,
    }
    dump_json(OUTPUT_JSON, payload)
    print("Phase 3 replay-support check")
    print("=" * 80)
    for row in rows:
        print(f"{row['label']}: supported={row['supported']} :: {row['message']}")
    print(f"\nWrote {OUTPUT_JSON}")


if __name__ == "__main__":
    main()