from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd


def write_validation_outputs(
    df: pd.DataFrame,
    *,
    output_dir: Path,
    source_root: Path,
) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    baseline = df[df["study"] == "ideal_sqr_direct_vs_echoed_multitone"].copy()
    native = df[df["study"] == "native_rich_multitone_sqr_cpsqr_feasibility"].copy()

    sanity = {
        "checks": [
            {
                "name": "baseline_best_direct_matches_expected_scale",
                "value": float(baseline[baseline["construction"] == "direct_multitone"]["strict_process_fidelity"].max()),
                "expected_comment": "Baseline direct strict process should remain around 0.65 process fidelity / 0.72 average gate fidelity.",
            },
            {
                "name": "native_rich_reaches_strict_gt_0p99",
                "value": float(native["strict_process_fidelity"].max()),
                "expected_comment": "The native-rich extension should exceed 0.99 on at least one strict ideal-SQR case.",
            },
            {
                "name": "native_rich_reaches_cpsqr_gt_0p999",
                "value": float(native["cpsqr_process_fidelity"].max()),
                "expected_comment": "The native-rich extension should show near-unit CPSQR on echoed cases.",
            },
            {
                "name": "row_count_by_study",
                "value": {str(key): int(val) for key, val in df.groupby("study").size().to_dict().items()},
                "expected_comment": "All four source studies must load into the definitive aggregation.",
            },
        ]
    }
    convergence = {
        "sources": [
            str((source_root / "ideal_sqr_direct_vs_echoed_multitone" / "data" / "validation_summary.json").resolve()),
            str((source_root / "multitone_sqr_arbitrary_fock_conditional_rotations" / "data" / "validation_summary.json").resolve()),
            str((source_root / "multitone_sqr_arbitrary_fock_conditional_rotations" / "data" / "echo_comparison_validation.json").resolve()),
            str((source_root / "ideal_sqr_direct_vs_echoed_multitone" / "data" / "analytic_summary.json").resolve()),
        ],
        "note": "This definitive-study pass aggregates source-study convergence/sanity artifacts. No new cross-study convergence rerun was performed in this iteration.",
    }
    (output_dir / "sanity_checks.json").write_text(json.dumps(sanity, indent=2), encoding="utf-8")
    (output_dir / "convergence_checks.json").write_text(json.dumps(convergence, indent=2), encoding="utf-8")
    return {"sanity": sanity, "convergence": convergence}
