from __future__ import annotations

import common as c
import run_structured_extension_analysis as sea


def main() -> None:
    record = sea._load_best_record("drcpsqr")
    n_cav = 12
    logical_maxiter = 3
    augmented_maxiter = 2

    print("[manual-basis] drcpsqr logical stable start", flush=True)
    logical_warm_start = sea._warm_start_from_record(record, n_cav=n_cav)
    logical_sequence = sea._sequence_from_record(record, n_cav=n_cav)
    logical_fit = c.fit_sequence(
        logical_sequence,
        n_cav=n_cav,
        seed=17,
        init_guess="heuristic",
        multistart=1,
        maxiter=logical_maxiter,
        use_fast_path=False,
        warm_start=logical_warm_start,
    )
    logical_fit["maxiter"] = logical_maxiter
    logical_record = sea._record_from_fit(
        record,
        logical_fit,
        n_cav=n_cav,
        search_phase="logical_refit",
        seed=17,
        subspace_label="logical_refit_3x2mix",
    )
    logical_wigner = sea._wigner_metrics_for_record(logical_record, n_cav=n_cav)
    print("[manual-basis] drcpsqr logical stable complete", flush=True)

    spectator_levels = sea._augmented_spectator_levels(logical_record, n_cav=n_cav)
    augmented_subspace, augmented_target = sea._augmented_target_for_spectators(spectator_levels, n_cav=n_cav)
    print(f"[manual-basis] drcpsqr augmented fast start spectators={spectator_levels}", flush=True)
    augmented_sequence = sea._sequence_from_record(record, n_cav=n_cav)
    augmented_fit = c.fit_sequence(
        augmented_sequence,
        n_cav=n_cav,
        seed=17,
        init_guess="heuristic",
        multistart=1,
        maxiter=augmented_maxiter,
        use_fast_path=True,
        warm_start=logical_warm_start,
        target_unitary=augmented_target,
        subspace=augmented_subspace,
    )
    augmented_fit["maxiter"] = augmented_maxiter
    augmented_record = sea._record_from_fit(
        record,
        augmented_fit,
        n_cav=n_cav,
        search_phase="augmented_refit",
        seed=17,
        subspace_label="augmented_refit_3x2mix",
    )
    augmented_wigner = sea._wigner_metrics_for_record(augmented_record, n_cav=n_cav)
    augmented_target_eval = c.evaluate_sequence(
        sea._sequence_from_record(augmented_record, n_cav=n_cav),
        n_cav=n_cav,
        target_unitary=augmented_target,
        subspace=augmented_subspace,
    )
    print("[manual-basis] drcpsqr augmented fast complete", flush=True)

    payload = {
        "base_case_id": str(record["case_id"]),
        "spectator_levels": [int(level) for level in spectator_levels],
        "logical_refit": {
            "record": logical_record,
            "wigner": logical_wigner,
        },
        "augmented_refit": {
            "record": augmented_record,
            "wigner": augmented_wigner,
            "augmented_target_fidelity": float(augmented_target_eval["fidelity"]),
        },
        "interpretation": (
            "supports"
            if augmented_wigner["mean_wigner_rms"] < logical_wigner["mean_wigner_rms"]
            and augmented_record["physical"]["by_n_cav"][str(n_cav)]["fidelity"] >= logical_record["physical"]["by_n_cav"][str(n_cav)]["fidelity"]
            else "does_not_support"
        ),
        "figure_stem": "drcpsqr_wigner_basis_extension_3x2mix",
    }
    sea._plot_wigner_triptych(
        title="Best D + R + CPSQR: logical(3) stable vs augmented(2) fast",
        baseline_record=logical_record,
        augmented_record=augmented_record,
        stem="drcpsqr_wigner_basis_extension_3x2mix",
        n_cav=n_cav,
    )
    c.save_json(c.DATA_DIR / "drcpsqr_basis_extension_detail_3x2mix.json", payload)
    print("__DONE__", flush=True)


if __name__ == "__main__":
    main()