"""Run the echoed-family duration refinement missing from the main study pass."""

from __future__ import annotations

import json
from pathlib import Path

from resume_study import execute_stage_with_resume
from run_study import duration_case_requests


def main() -> None:
    rows = execute_stage_with_resume(duration_case_requests(), ("echoed_independent",))
    output_path = Path(__file__).resolve().parents[1] / "data" / "duration_echoed_results.json"
    output_path.write_text(json.dumps({"rows": rows}, indent=2), encoding="utf-8")
    print(f"Saved echoed duration rows to {output_path}")


if __name__ == "__main__":
    main()
