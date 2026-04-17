from __future__ import annotations

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import common as c
import run_closed_system_followup as followup


JOB_DIR = c.DATA_DIR / "closed_system_job_specs"


def _write_job(stem: str, job: dict) -> None:
    c.save_json(JOB_DIR / f"{stem}.json", job)


def _job(case: dict, *, maxiter: int, seed: int, tag: str) -> dict:
    return {
        "case": c.json_ready(case),
        "n_cav": 12,
        "search_phase": tag,
        "seed": int(seed),
        "init_guess": "heuristic",
        "maxiter": int(maxiter),
        "multistart": 1,
        "record_updates": {"neighborhood_tag": tag},
    }


def main() -> None:
    JOB_DIR.mkdir(parents=True, exist_ok=True)
    sqr = followup._load_best_record("drsqr")
    cps = followup._load_best_record("drcpsqr")

    sqr_cases = {
        "drsqr_exact_s17": _job(followup._make_case(sqr, blocks=5, max_tones=4, levels=(0, 1, 2, 3), order_tokens=("R", "D", "SQR")), maxiter=24, seed=17, tag="exact"),
        "drsqr_exact_s42": _job(followup._make_case(sqr, blocks=5, max_tones=4, levels=(0, 1, 2, 3), order_tokens=("R", "D", "SQR")), maxiter=24, seed=42, tag="exact"),
        "drsqr_exact_s73": _job(followup._make_case(sqr, blocks=5, max_tones=4, levels=(0, 1, 2, 3), order_tokens=("R", "D", "SQR")), maxiter=24, seed=73, tag="exact"),
        "drsqr_blk4_a4_rds": _job(followup._make_case(sqr, blocks=4, max_tones=4, levels=(0, 1, 2, 3), order_tokens=("R", "D", "SQR")), maxiter=24, seed=17, tag="block_neighbor"),
        "drsqr_blk5_a3_rds_012": _job(followup._make_case(sqr, blocks=5, max_tones=3, levels=(0, 1, 2), order_tokens=("R", "D", "SQR")), maxiter=24, seed=17, tag="level_neighbor"),
        "drsqr_blk5_a3_rds_013": _job(followup._make_case(sqr, blocks=5, max_tones=3, levels=(0, 1, 3), order_tokens=("R", "D", "SQR")), maxiter=24, seed=17, tag="level_neighbor"),
        "drsqr_blk5_a5_rds_01234": _job(followup._make_case(sqr, blocks=5, max_tones=5, levels=(0, 1, 2, 3, 4), order_tokens=("R", "D", "SQR")), maxiter=24, seed=17, tag="level_neighbor"),
        "drsqr_blk5_a4_drs": _job(followup._make_case(sqr, blocks=5, max_tones=4, levels=(0, 1, 2, 3), order_tokens=("D", "R", "SQR")), maxiter=24, seed=17, tag="ordering"),
        "drsqr_blk5_a4_dsr": _job(followup._make_case(sqr, blocks=5, max_tones=4, levels=(0, 1, 2, 3), order_tokens=("D", "SQR", "R")), maxiter=24, seed=17, tag="ordering"),
    }

    cps_cases = {
        "drcpsqr_exact_s17": _job(followup._make_case(cps, blocks=5, max_tones=2, levels=(1, 2), order_tokens=("D", "R", "CPSQR")), maxiter=24, seed=17, tag="exact"),
        "drcpsqr_exact_s42": _job(followup._make_case(cps, blocks=5, max_tones=2, levels=(1, 2), order_tokens=("D", "R", "CPSQR")), maxiter=24, seed=42, tag="exact"),
        "drcpsqr_exact_s73": _job(followup._make_case(cps, blocks=5, max_tones=2, levels=(1, 2), order_tokens=("D", "R", "CPSQR")), maxiter=24, seed=73, tag="exact"),
        "drcpsqr_blk4_a2_drcp": _job(followup._make_case(cps, blocks=4, max_tones=2, levels=(1, 2), order_tokens=("D", "R", "CPSQR")), maxiter=24, seed=17, tag="block_neighbor"),
        "drcpsqr_blk5_a1_drcp_1": _job(followup._make_case(cps, blocks=5, max_tones=1, levels=(1,), order_tokens=("D", "R", "CPSQR")), maxiter=24, seed=17, tag="level_neighbor"),
        "drcpsqr_blk5_a1_drcp_2": _job(followup._make_case(cps, blocks=5, max_tones=1, levels=(2,), order_tokens=("D", "R", "CPSQR")), maxiter=24, seed=17, tag="level_neighbor"),
        "drcpsqr_blk5_a3_drcp_012": _job(followup._make_case(cps, blocks=5, max_tones=3, levels=(0, 1, 2), order_tokens=("D", "R", "CPSQR")), maxiter=24, seed=17, tag="level_neighbor"),
        "drcpsqr_blk5_a3_drcp_123": _job(followup._make_case(cps, blocks=5, max_tones=3, levels=(1, 2, 3), order_tokens=("D", "R", "CPSQR")), maxiter=24, seed=17, tag="level_neighbor"),
        "drcpsqr_blk5_a2_rdcp": _job(followup._make_case(cps, blocks=5, max_tones=2, levels=(1, 2), order_tokens=("R", "D", "CPSQR")), maxiter=24, seed=17, tag="ordering"),
        "drcpsqr_blk5_a2_dcpr": _job(followup._make_case(cps, blocks=5, max_tones=2, levels=(1, 2), order_tokens=("D", "CPSQR", "R")), maxiter=24, seed=17, tag="ordering"),
    }

    manifest = {**sqr_cases, **cps_cases}
    for stem, job in manifest.items():
        _write_job(stem, job)

    c.save_json(JOB_DIR / "manifest.json", {"jobs": sorted(manifest)})
    print(f"wrote {len(manifest)} job specs to {JOB_DIR}")


if __name__ == "__main__":
    main()